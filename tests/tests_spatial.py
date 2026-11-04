"""Independent limits for wall fluxes, injected columns and photon sheets."""
from copy import deepcopy
import csv
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
from mixture_kinetics import ROOT,RetryStep
from particle_kinetics import EV,PhaseGrid
from slab_transport import solve_transfer
from spatial_sources import deposition
from spatial_editor import validate
from tests_mixture import config,build
from tests_mixture import synthetic_bank
from tests_chemistry import configuration,build as chemical_build
from mixture_chemistry import records,rate_fit
from mixture_kinetics import MixtureSlab


class InternalPhotons(unittest.TestCase):
    def solve(self,opacity=0.,scatter=0.,backend='numpy'):
        return solve_transfer([.25,.5,.75,1.],np.full((4,2),opacity),np.zeros((4,2)),np.zeros(2),
            scatter=np.full(4,scatter),angular_points=8,backend=backend,workers=2,
            internal_beams=[dict(position_cm=.4,direction='right',mu=.5,spectrum=np.array([2.,3.]))])

    def test_transparent_partial_source_cell_and_face_flux(self):
        out=self.solve();np.testing.assert_allclose(out['mean_4piJ'],np.array([0.,.4,1.,1.])[:,None]*[2.,3.],atol=1e-14)
        np.testing.assert_allclose(out['outgoing_right_flux'],[1.,1.5]);np.testing.assert_allclose(out['outgoing_left_flux'],0.)
        self.assertLess(out['maximum_relative_balance_residual'],1e-14)

    def test_absorption_uses_distance_from_actual_sheet(self):
        out=self.solve(opacity=2.)
        np.testing.assert_allclose(out['outgoing_right_flux'],.5*np.array([2.,3.])*np.exp(-2*.6/.5),rtol=1e-14)
        np.testing.assert_allclose(out['true_absorbed_flux'],.5*np.array([2.,3.])*(1-np.exp(-2*.6/.5)),rtol=1e-14)
        self.assertLess(out['maximum_relative_balance_residual'],2e-14)

    def test_scattering_conserves_internal_source_power_and_backends_match(self):
        a=self.solve(scatter=2.);b=self.solve(scatter=2.,backend='numba')
        self.assertGreater(a['outgoing_left_flux'][0],0.)
        np.testing.assert_allclose(a['outgoing_left_flux']+a['outgoing_right_flux'],[1.,1.5],rtol=2e-13)
        np.testing.assert_allclose(a['mean_4piJ'],b['mean_4piJ'],rtol=3e-13)
        self.assertLess(b['maximum_relative_balance_residual'],3e-13)

    def test_left_right_symmetry_and_sheet_at_boundary(self):
        args=([.5,1.],np.ones((2,1)),np.zeros((2,1)),np.zeros(1))
        left=solve_transfer(*args,internal_beams=[dict(position_cm=.3,direction='left',spectrum=np.ones(1))])
        right=solve_transfer(*args,internal_beams=[dict(position_cm=.7,direction='right',spectrum=np.ones(1))])
        np.testing.assert_allclose(left['mean_4piJ'],right['mean_4piJ'][::-1],rtol=1e-14)
        edge=solve_transfer(*args,internal_beams=[dict(position_cm=0.,direction='right',spectrum=np.ones(1))])
        beam=solve_transfer(*args[:3],np.ones(1))
        np.testing.assert_allclose(edge['mean_4piJ'],beam['mean_4piJ'],rtol=1e-14)


class ParticleSources(unittest.TestCase):
    def test_bgk_entropy_handles_narrow_injected_beam_without_changing_moments(self):
        g=PhaseGrid(dict(energy_min_eV=1e-6,energy_max_eV=50.,energy_points=16,angular_points=4),1.008*1.66053906892e-24)
        f=g.ion_initial(dict(Ti_K=650.,uz_cm_s=3e6),1e8)+g.ion_initial(dict(Ti_K=650.),1e7)
        equilibrium=g.maximum_entropy(f)
        self.assertGreaterEqual(equilibrium.min(),0.)
        for feature in [np.ones(g.shape),np.broadcast_to(g.e[:,None],g.shape),g.vz]:
            a=float((f*feature).sum());b=float((equilibrium*feature).sum());self.assertLess(abs(a-b)/max(abs(a),1.),3e-12)
    def test_sheet_column_independent_of_grid_and_finite_width(self):
        for cells in [1,4,11]:
            for width in [0.,.2]:
                weights=deposition(.5,width,1.,cells)
                self.assertAlmostEqual(float(weights.sum()/cells),1.,places=14)
        with self.assertRaises(ValueError):deposition(0.,.1,1.,4)

    def example(self,cells=4):
        cfg=config();cfg.update(cells=cells,thickness_cm=1.,spatial_sources=[
            dict(id='ion',mode='ions',species='H+',position_cm=.2,start_s=.0017,stop_s=.0063,rate_cm2_s=1e8,distribution=dict(Ti_K=500.)),
            dict(id='e',mode='electrons',position_cm=.7,start_s=.002,stop_s=.007,rate_cm2_s=2e8,distribution=dict(Te_eV=.1,hot_Te_eV=.1))])
        return cfg

    def test_pulse_knots_and_exact_external_columns_in_all_balances(self):
        model=build(self.example());model.solve();final=model.snapshots[-1]
        np.testing.assert_allclose(final['source_column'],[1e8*(.0063-.0017),2e8*(.007-.002)],rtol=2e-14)
        self.assertIn(.0017,[s['time_s'] for s in model.snapshots]);self.assertIn(.0063,[s['time_s'] for s in model.snapshots])
        self.assertLess(abs(final['charge'])/model.initial_nuclei.sum(),1e-12)
        self.assertLess(abs(final['nuclei']).max()/model.initial_nuclei.sum(),1e-12)
        self.assertLess(abs(final['energy'])/model.initial_energy,1e-12)
        self.assertLessEqual(model.steps,18)  # output/source knots must not collapse the proposed step

    def test_rejected_step_does_not_duplicate_injected_particles(self):
        model=build(self.example());real=model._local;count=[0]
        def local(z,dt,psi):
            count[0]+=1
            if count[0]==1:raise RetryStep('test retry after injection')
            return real(z,dt,psi)
        with patch.object(model,'_local',side_effect=local):model.solve()
        self.assertGreater(model.rejected,0)
        np.testing.assert_allclose(model.sources.column,[4.6e5,1e6],rtol=3e-14)

    def test_source_seed_selects_initially_empty_excited_level(self):
        cfg=self.example(1);cfg['level_selection']['mode']='auto';cfg['spatial_sources']=[dict(id='upper',mode='particles',species='Hup',position_cm=.5,rate_cm2_s=1e4,distribution=dict(Ti_K=300.))]
        m=build(cfg,upper=True);self.assertTrue(m.selection.active[m.index['Hup']])

    def test_ion_source_rejects_neutral_and_editor_validates_geometry(self):
        cfg=self.example();cfg['spatial_sources'][0]['species']='H'
        with self.assertRaisesRegex(ValueError,'положительный'):build(cfg)
        cfg=self.example();validate(cfg);cfg['spatial_sources'][0]['position_cm']=2.
        with self.assertRaises(ValueError):validate(cfg)


class SurfacePhysics(unittest.TestCase):
    def model(self,quench=1.,neutralize=1.,sticking=0.,electron_absorption=0.):
        cfg=configuration([]);cfg.update(cells=3,thickness_cm=1.)
        cfg['initial_cm3']={'He_metastable':1e7,'He+':1e6,'H':1e7}
        cfg['particle_boundary']=dict(kind='walls',walls={side:dict(quenching_probability=quench,neutralization_probability=neutralize,sticking_probability=sticking,electron_absorption_probability=electron_absorption) for side in ['left','right']})
        return chemical_build(cfg)

    def test_quenching_at_two_faces_matches_incident_velocity_flux(self):
        m=self.model();j=m.index['He_metastable'];g=m.gh['He'];dt=.05*m.dz/max(abs(m.ge.vz).max(),max(abs(x.vz).max() for x in m.gh.values()))
        incident=float((m.P[0,j]*np.maximum(-g.vz,0)).sum())*dt
        old=m.n[:,j].sum()*m.dz;m._transport(dt)
        self.assertAlmostEqual(m.walls.counts[0,j,1]/incident,1.,places=13)
        self.assertAlmostEqual(m.walls.counts[1,j,1]/incident,1.,places=13)
        self.assertAlmostEqual((old-m.n[:,j].sum()*m.dz)/(2*incident),1.,places=10)
        snap=m.snapshot(dt)
        self.assertLess(abs(snap['nuclei']).max()/m.initial_nuclei.sum(),1e-13)
        self.assertLess(abs(snap['charge'])/m.initial_nuclei.sum(),1e-13)
        self.assertLess(abs(snap['energy'])/m.initial_energy,1e-13)
        self.assertTrue(np.all(m.walls.heat>0))

    def test_reflecting_walls_preserve_species_energy_and_charge(self):
        m=self.model(quench=0.,neutralize=0.);old=m.n.sum(axis=0)*m.dz;energy=m.total_energy();charge=m.total_charge()
        m._transport(2e-8)
        np.testing.assert_allclose(m.n.sum(axis=0)*m.dz,old,rtol=1e-13)
        self.assertAlmostEqual(m.total_energy()/energy,1.,places=13)
        self.assertAlmostEqual(m.total_charge()/charge if charge else m.total_charge(),1. if charge else 0.,places=6)

    def test_absorbing_wall_losses_close_charge_nuclei_energy(self):
        m=self.model(sticking=1.,electron_absorption=1.);m._transport(2e-8);snap=m.snapshot(2e-8)
        self.assertGreater(m.species_out.sum(),0.);self.assertGreater(m.electron_out,0.)
        self.assertLess(abs(snap['nuclei']).max()/m.initial_nuclei.sum(),1e-13)
        self.assertLess(abs(snap['charge'])/m.initial_nuclei.sum(),1e-13)
        self.assertLess(abs(snap['energy'])/m.initial_energy,1e-13)


class PublishedChemistry(unittest.TestCase):
    def test_misspelled_alternative_dataset_is_rejected(self):
        cfg=configuration([]);cfg['chemistry']['datasets']={'He_H_chemiionization':'GM52'}
        with self.assertRaisesRegex(ValueError,'Неизвестный набор'):chemical_build(cfg)
    def test_parallel_product_operator_matches_single_worker(self):
        cfg=configuration(['HE25:15:PI','HE25:14:AI']);cfg['initial_cm3']={'He':1e9,'He_metastable':1e8}
        a=chemical_build(cfg);cfg['performance']['workers']=2;b=chemical_build(cfg)
        x=a.chemistry.evolve(0,1e-4);y=b.chemistry.evolve(0,1e-4)
        np.testing.assert_allclose(x,y,rtol=2e-13,atol=1e-40)
        np.testing.assert_allclose(a.P,b.P,rtol=2e-13,atol=1e-40);np.testing.assert_allclose(a.F,b.F,rtol=2e-13,atol=1e-40)

    def test_digitized_cross_section_integrates_to_independent_primary_rate_curve(self):
        from scipy.integrate import simpson
        from particle_kinetics import KB
        rows=records(ROOT/'data/reactions/He--Li/cross_sections.csv');fits=records(ROOT/'data/reactions/He--Li/rate_coefficients.csv')
        for ident in ['BR15:3->5','BR15:4->6','BR15:5->6']:
            tab=np.array([(float(r['energy_eV']),float(r['sigma_cm2'])) for r in rows if r['id']==ident]);T=10000.;theta=KB*T/EV
            e=np.unique(np.r_[tab[:,0],np.geomspace(max(tab[0,0],1e-8),tab[-1,0],10001)]);sig=np.interp(e,tab[:,0],tab[:,1])
            mu=4.002602*6.94/(4.002602+6.94)*1.66053906892e-24
            value=np.sqrt(8*EV/(np.pi*mu))*simpson(e*sig*np.exp(-e/theta),x=e)/theta**1.5
            primary=rate_fit([r for r in fits if r['id']==ident],T)
            self.assertLess(abs(value/primary-1),.05)
    def test_updated_penning_total_and_branching(self):
        cfg=configuration(['GM25:He3S_H:PI','GM25:He3S_H:AI']);cfg['chemistry']['datasets']={'He_H_chemiionization':'GM25'}
        m=chemical_build(cfg);rates=m.chemistry.event_rates(0);a,b=[rates[r['channel']] for r in m.chemistry.reactions]
        self.assertAlmostEqual(a/b,9.,places=12)
        selected=records(ROOT/'data/reactions/H--He/rate_coefficients.csv');rows=[r for r in selected if r['id']=='GM25:He3S_H:PI']
        self.assertAlmostEqual(rate_fit(rows,500.)/(.9*1.02e-9),1.,places=13)

    def test_radiative_CT_emits_energy_instead_of_heating_heavy_products(self):
        cfg=configuration(['GM25:Heplus_H:RCT']);cfg['initial_cm3']={'He+':1e8,'H':1e9};m=chemical_build(cfg)
        before=m.total_energy();rates=m.chemistry.event_rates(0);emission=m.chemistry.photon_emission(0)
        expected=sum(rates[r['channel']]*r['photon_eV']*EV for r in m.chemistry.reactions)
        self.assertAlmostEqual(float(emission@m.w)/expected,1.,places=13)
        events=m.chemistry.evolve(0,1e-6);self.assertGreater(events.sum(),0.)
        self.assertLess(abs(m.total_energy()-before-m.photon_work)/before,1e-12)
        self.assertEqual(int(m.network.matrix[-1,m.chemistry.reactions[0]['channel']]),1)

    def test_inverse_nonradiative_CT_detailed_balance(self):
        path=ROOT/'data/reactions/He--Li/cross_sections.csv';rows=records(path)
        a=np.array([(float(r['energy_eV']),float(r['sigma_cm2'])) for r in rows if r['id']=='BR15:2->3'])
        b=np.array([(float(r['energy_eV']),float(r['sigma_cm2'])) for r in rows if r['id']=='BR15:3->2'])
        np.testing.assert_allclose(a[1:,0]*a[1:,1]*4,b[:,0]*b[:,1]*3,rtol=3e-14,atol=1e-50)

    def test_tabulated_rate_does_not_extrapolate_and_reproduces_knots(self):
        rows=[r for r in records(ROOT/'data/reactions/H--He/rate_coefficients.csv') if r['id']=='GM25:He1S_H:PI']
        for r in rows:self.assertAlmostEqual(rate_fit(rows,float(r['reference_T_K']))/float(r['alpha_cm3_s']),1.,places=13)
        with self.assertRaises(ValueError):rate_fit(rows,100.)


class SpatialOutput(unittest.TestCase):
    def test_line_columns_use_cell_width_and_rank_each_position_independently(self):
        cfg=config();cfg.update(cells=2,thickness_cm=1.)
        cfg['initial_cm3']['Hup']=1e4;cfg['diagnostics']={'phase_history_csv':False}
        cfg['processes'].update(spontaneous=True,photo_bound_bound=True)
        bank,photo=synthetic_bank(upper=True);bank.lines=[dict(id='SYNTHETIC_LINE',lower='H0',upper='Hup',A_s=1e3,wavelength_nm=1239.841984,source='SYNTHETIC TEST ONLY',spectral_kind='line',pump=True)]
        with patch('mixture_kinetics.MixtureData',return_value=bank),patch('mixture_kinetics.PhotoBank',return_value=photo):m=MixtureSlab(cfg,ROOT)
        m.snapshots=[m.snapshot(0.)];m.elapsed=0.
        with tempfile.TemporaryDirectory() as folder:
            m.export_full(folder);rows=records(Path(folder)/'lines_by_position.csv')
            self.assertEqual(len(rows),2);self.assertTrue(all(int(r['emission_rank'])==1 and int(r['absorption_rank'])==1 for r in rows))
            line=m.lines[0];up,down,sigma,ratio=line['optical'];full_tau=(m.n[0,line['i']]*ratio-m.n[0,line['j']])*sigma*m.length
            self.assertAlmostEqual(float(rows[0]['maximum_cell_tau'])/float(full_tau.max()),.5,places=13)
            self.assertEqual([float(r['x_cm']) for r in rows],[.25,.75])
            self.assertTrue((Path(folder)/'spatial_radiation_history.npz').exists())


if __name__=='__main__':unittest.main()
