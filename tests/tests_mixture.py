"""Analytic limits for the coupled mixture. Synthetic tables are labelled.

These small tests run without network; no synthetic coefficient is shipped
as physical atomic data or used by the real H/He/Li examples.
"""
from copy import deepcopy
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import numpy as np
import mpmath as mp
from threadpoolctl import threadpool_limits
from atomic_database import HC_EV_CM
from mixture_photo import hydrogen_sigma,BOHR_CM2,HARTREE_EV
from level_selection import LevelSelection,frozen_evolution,positive_metzler_exponential
from mixture_kinetics import MixtureSlab,ROOT,EV,conservative_profile
from particle_kinetics import PhaseGrid,ME
import kinetic_kernels as K


def synthetic_bank(upper=False,photo_sigma=1e-18):
    states=[]
    for element,mass,g,ip in [('H',1.008,2,2.),('He',4.002602,1,3.)]:
        for q,energy,weight in [(0,0.,g),(1,ip,1 if element=='H' else 2)]:
            states.append(dict(id=f'{element}{q}',label=f'SYNTHETIC {element} {q} 1 {weight}S',element=element,species=element+str(q),charge=q,n=1,multiplicity=weight,L=0,g=weight,energy_cm1=energy/HC_EV_CM,
                excitation_eV=0.,mass_u=mass,nuclei={element:1},is_ground=True,configuration='SYNTHETIC'))
    if upper:states.append(dict(states[0],id='Hup',label='SYNTHETIC H 2 2P',n=2,L=1,g=6,energy_cm1=1./HC_EV_CM,excitation_eV=1.,is_ground=False))
    index={s['id']:j for j,s in enumerate(states)};cross={};source={}
    for a,b,gap in [('H0','H1',2.),('He0','He1',3.)]:
        cross[(a,b,'ionization')]=np.array([[gap,0.],[gap+.001,2e-17],[100.,2e-17]])
        source[(a,b,'ionization')]='SYNTHETIC test table, not physical data'
    if upper:
        cross[('H0','Hup','excitation')]=np.array([[1.,0.],[1.001,2e-16],[100.,2e-16]])
        cross[('Hup','H0','deexcitation')]=np.array([[1e-6,2e-16],[100.,2e-16]])
        source[('H0','Hup','excitation')]=source[('Hup','H0','deexcitation')]='SYNTHETIC test'
    aliases={'H':'H0','H+':'H1','He':'He0','He+':'He1','Hup':'Hup'}
    bank=SimpleNamespace(states=states,index=index,elements=['H','He'],lines=[],cross=cross,cross_source=source,elastic={},
        parent_index={0:1,2:3},manifests={},by={s['id']:s for s in states})
    if upper:bank.parent_index[4]=1
    bank.resolve=lambda name:index[aliases.get(name,name)]
    photo=SimpleNamespace(tables={name:np.array([[edge,photo_sigma],[100.,photo_sigma]]) for name,edge in [('H0',2.),('He0',3.)]},
        sources={'H0':'SYNTHETIC PI','He0':'SYNTHETIC PI'},gaps=[],receipts={})
    return bank,photo


def config():
    return dict(schema='MixtureKineticSlab-1.8',initial_cm3={'H':1e8,'He':1e8,'H+':1e7,'He+':1e7},cells=1,thickness_cm=.01,
        level_selection=dict(mode='all',relative_tolerance=1e-6,missing_data='report'),
        electron_initial=dict(Te_eV=.7,hot_fraction=0.,hot_Te_eV=.5,angular_bias=.2),heavy_initial={'H':{'Ti_K':300.},'He':{'Ti_K':300.}},
        phase_grid=dict(energy_min_eV=1e-6,energy_max_eV=10.,energy_points=32,angular_points=4),
        heavy_phase_grid=dict(energy_min_eV=1e-7,energy_max_eV=1.,energy_points=16,angular_points=4),
        processes={key:False for key in ['electron_bound_bound','electron_ionization','spontaneous','two_photon','photo_bound_bound','photo_bound_free','radiative_recombination','stimulated_recombination']},
        radiation=dict(continuum_points=80,line_points=25,electron_scattering=False),
        collisions=dict(ee_bgk=False,ii_bgk=False,ei_binary=False,neutral_elastic=False),performance=dict(workers=2,radiation_backend='numpy'),
        solver=dict(t_end_s=.01,max_step_s=.001,min_step_s=1e-12,samples=3,max_steps=10000))


def build(cfg,upper=False,photo_sigma=1e-18):
    bank,photo=synthetic_bank(upper,photo_sigma)
    with patch('mixture_kinetics.MixtureData',return_value=bank),patch('mixture_kinetics.PhotoBank',return_value=photo):
        return MixtureSlab(cfg,ROOT)


class PhotoPhysics(unittest.TestCase):
    def test_ground_coulomb_integral_matches_independent_stobbe_formula(self):
        sigma0=512*np.pi**2/(3*np.e**4)/137.035999084*BOHR_CM2
        for excess in [1e-5,.1,1.,10.,50.]:
            eps=np.sqrt(2*excess/HARTREE_EV)
            exact=sigma0*np.exp(4-4*np.arctan(eps)/eps)/(1+eps*eps)**4/(-np.expm1(-2*np.pi/eps))
            self.assertAlmostEqual(hydrogen_sigma(1,0,excess)/exact,1.,places=8)

    def test_excited_2s_matches_direct_radial_quadrature(self):
        from scipy.special import roots_laguerre
        with mp.workdps(25):
            excess=.1;k=mp.sqrt(2*excess)
            nodes,weights=roots_laguerre(160)
            integral=mp.mpf(0)
            for node,weight in zip(nodes,weights):
                if weight<1e-32:continue
                r=mp.mpf(float(2*node))
                integral+=float(weight)*mp.sqrt(2/(mp.pi*k))*mp.coulombf(1,-1/k,k*r)*r*r*(2-r)/mp.sqrt(2)
            sigma=4*mp.pi**2/137.035999084*(excess+mp.mpf(1)/8)/3*abs(integral)**2*BOHR_CM2
        self.assertAlmostEqual(hydrogen_sigma(2,0,.1*HARTREE_EV)/float(sigma),1.,places=8)

    def test_hydrogenic_charge_scaling_and_positive_high_l(self):
        self.assertAlmostEqual(hydrogen_sigma(3,2,4.,2)*4/hydrogen_sigma(3,2,1.,1),1.,places=12)
        self.assertGreater(hydrogen_sigma(5,4,2.),0.)


class BasisPhysics(unittest.TestCase):
    def test_positive_uniformization_matches_stiff_analytic_cascade(self):
        from scipy.linalg import expm
        matrix=np.zeros((6,6));matrix[0,0]=-1e10;matrix[1,0]=1e10;matrix[1,1]=-5e4;matrix[2,1]=5e4;matrix[3:,:3]=np.eye(3)
        result=positive_metzler_exponential(matrix,2e-8)
        self.assertGreaterEqual(result.min(),0.)
        np.testing.assert_allclose(result,expm(matrix*2e-8),rtol=3e-13,atol=1e-26)

    def test_conservative_photon_projection_preserves_count_and_mean(self):
        energy=np.geomspace(.01,20.,201);weights=np.gradient(energy)
        shape=np.exp(-((energy-7.1)/2.)**2)
        projected=conservative_profile(shape,weights,energy,7.)
        self.assertGreaterEqual(projected.min(),0.)
        self.assertAlmostEqual(float(projected@weights),1.,places=14)
        self.assertAlmostEqual(float((projected*weights)@energy),7.,places=12)

    def states(self):return [dict(id=str(i),label=str(i),element='H',charge=0,n=i+1,excitation_eV=float(i)) for i in range(3)]

    def test_sequential_pumping_reaches_energy_above_single_photon(self):
        channels=[dict(i=0,j=1,kind='photoexcitation'),dict(i=1,j=2,kind='photoexcitation')]
        sel=LevelSelection(self.states(),channels,np.array([[1.,0.,0.]]),np.array([[10.,10.]]),1.)
        self.assertTrue(sel.active[2])

    def test_short_lived_upper_is_kept_by_event_visits(self):
        channels=[dict(i=0,j=1,kind='photoexcitation'),dict(i=1,j=0,kind='spontaneous')]
        sel=LevelSelection(self.states(),channels,np.array([[1.,0.,0.]]),np.array([[1e3,1e9]]),.01,dict(relative_tolerance=1e-4))
        self.assertTrue(sel.active[1]);self.assertGreater(sel.visit[1],1.)

    def test_trace_metastable_is_an_independent_selection_seed(self):
        channels=[dict(i=1,j=2,kind='excitation')]
        sel=LevelSelection(self.states(),channels,np.array([[1e12,1.,0.]]),np.array([[1.]]),1.)
        self.assertTrue(sel.active[2])

    def test_later_leakage_expands_without_population_renormalization(self):
        channels=[dict(i=0,j=1,kind='excitation')]
        sel=LevelSelection(self.states(),channels,np.array([[1.,0.,0.]]),np.array([[0.]]),1.)
        leak=sel.leakage(np.array([1.,0.,0.]),np.array([1.]),.1)
        mask=sel.candidates(leak);self.assertTrue(mask[1]);sel.expand(mask)
        self.assertTrue(sel.active[1]);self.assertEqual(sel.restarts,1)

    def test_matrix_residence_has_analytic_two_level_event_budget(self):
        channels=[dict(i=0,j=1,kind='spontaneous')]
        initial=np.array([[1.,3.],[0.,0.]])
        final,residence=frozen_evolution(initial,channels,np.array([20.]),.2)
        np.testing.assert_allclose(final[0],initial[0]*np.exp(-4),rtol=1e-14)
        np.testing.assert_allclose(20*residence[0],initial[0]-final[0],rtol=1e-14)
        np.testing.assert_allclose(final.sum(axis=0),initial.sum(axis=0),rtol=1e-14)


class CoupledPhysics(unittest.TestCase):
    def test_two_elements_really_ionize_and_share_dynamic_electrons(self):
        c=config();c['processes']['electron_ionization']=True
        m=build(c);initial=m.n.copy();energy=m.total_energy()
        with threadpool_limits(limits=1):m.solve()
        self.assertGreater(m.n[0,1],initial[0,1]);self.assertGreater(m.n[0,3],initial[0,3])
        self.assertGreater(m.F.sum(),m.initial_F.sum())
        self.assertLess(abs(m.total_energy()/energy-1),1e-9)
        self.assertLess(abs(m.snapshots[-1]['charge']),1e-6)
        np.testing.assert_allclose(m.n[0]-initial[0],m.matrix@m.integral[0],rtol=1e-7,atol=1e-6)

    def test_photoionization_without_seed_electrons_matches_exponential(self):
        c=config();c['initial_cm3']={'H':1e8,'He':1e8};c['processes']['photo_bound_free']=True
        m=build(c);m.source=np.where((m.energy>=3.5)&(m.energy<=4.),1e-8,0.)
        rates=m.rates(0,m.field(0.)['mean_4piJ'][0])[0]
        k=next(k for k,r in enumerate(m.channels) if r['id']=='PI:H0');gamma=rates[k]
        m.solve();self.assertGreater(m.F.sum(),0.)
        self.assertAlmostEqual(m.n[0,0]/(1e8*np.exp(-gamma*m.end)),1.,places=9)
        self.assertLess(abs(m.snapshots[-1]['energy'])/m.initial_energy,1e-8)

    def test_full_distribution_rr_reduces_ions_and_conserves_matter_plus_photons(self):
        c=config();c['initial_cm3']={'H+':1e8,'He+':1e8};c['processes']['radiative_recombination']=True
        m=build(c,photo_sigma=1e-12);initial=m.n.copy();m.solve()
        self.assertLess(m.n[0,1],initial[0,1]);self.assertLess(m.n[0,3],initial[0,3]);self.assertLess(m.F.sum(),m.initial_F.sum())
        self.assertLess(abs(m.snapshots[-1]['energy'])/m.initial_energy,1e-7)
        np.testing.assert_allclose(m.n[0]-initial[0],m.matrix@m.integral[0],rtol=1e-7,atol=1e-5)

    def test_hot_electron_tail_changes_selected_upper_level(self):
        c=config();c['level_selection']['mode']='auto';c['processes']['electron_bound_bound']=True;c['electron_initial'].update(Te_eV=.02,hot_Te_eV=1.,hot_fraction=0.)
        cold=build(c,upper=True);c['electron_initial']['hot_fraction']=.1;hot=build(c,upper=True)
        self.assertFalse(cold.selection.active[4]);self.assertTrue(hot.selection.active[4])

    def test_nonzero_reservoir_level_is_included_before_transport(self):
        c=config();c['level_selection']['mode']='auto';c['particle_boundary']={'kind':'reservoir','left_heavy':{'Hup':{'density_cm3':1e5,'Ti_K':300.}}}
        m=build(c,upper=True);self.assertTrue(m.selection.active[4]);self.assertAlmostEqual(m.selection.reference['H']/1e5,1.,places=14)

    def test_missing_excitation_to_omitted_accessible_level_is_reported(self):
        c=config();c['level_selection']['mode']='auto';c['processes']['electron_bound_bound']=True
        bank,photo=synthetic_bank(upper=True);bank.cross.pop(('H0','Hup','excitation'))
        with patch('mixture_kinetics.MixtureData',return_value=bank),patch('mixture_kinetics.PhotoBank',return_value=photo):m=MixtureSlab(c,ROOT)
        self.assertFalse(m.selection.active[4]);m._check_data(0,m.field(0.)['mean_4piJ'][0])
        self.assertIn('missing electron channel H0 -> Hup',m.data_issues)

    def test_missing_rr_to_omitted_candidate_is_not_hidden_by_selection(self):
        c=config();c['level_selection']['mode']='auto';c['processes']['radiative_recombination']=True
        m=build(c,upper=True);self.assertFalse(m.selection.active[4]);m._check_data(0,m.field(0.)['mean_4piJ'][0])
        self.assertIn('missing level-resolved RR: Hup',m.data_issues)

    def test_one_and_multiple_workers_give_same_trajectory(self):
        c=config();c['processes']['electron_ionization']=True;c['performance']['workers']=1
        a=build(c);a.solve();c['performance']['workers']=2;b=build(c);b.solve()
        np.testing.assert_allclose(a.F,b.F,rtol=2e-14,atol=1e-40);np.testing.assert_allclose(a.P,b.P,rtol=2e-14,atol=1e-40)

    def test_automatic_basis_restart_matches_full_basis_without_lost_history(self):
        c=config();c['processes']['electron_bound_bound']=True
        full=build(c,upper=True);full.solve()
        adaptive=build(c,upper=True);adaptive.selection.active[4]=False;adaptive.solve()
        self.assertGreaterEqual(adaptive.selection.restarts,1)
        np.testing.assert_allclose(full.F,adaptive.F,rtol=2e-14,atol=1e-40)
        np.testing.assert_allclose(full.P,adaptive.P,rtol=2e-14,atol=1e-40)
        np.testing.assert_allclose(full.integral,adaptive.integral,rtol=2e-14,atol=1e-40)

    def test_omitted_photo_parent_does_not_absorb_without_a_kinetic_event(self):
        c=config();c['processes']['photo_bound_free']=True;m=build(c)
        m.source=np.where((m.energy>3.5)&(m.energy<4.),1e-8,0.)
        m.selection.active[[1,3]]=False
        field=m.field(0.)
        np.testing.assert_array_equal(field['opacity'],0.)
        self.assertEqual(m.snapshot(0.)['reaction_photon_residual'],0.)
        m.selection.active[:]=True;self.assertGreater(m.field(0.)['opacity'].max(),0.)

    def test_all_populations_and_histories_are_exported_with_data_status(self):
        c=config();c['processes']['electron_ionization']=True;c['diagnostics']={'reaction_levels':['H','He+']}
        m=build(c);m.solve()
        with tempfile.TemporaryDirectory() as d:
            q=m.export(d);self.assertEqual(q['charge_status'],'passed');self.assertEqual(q['energy_status'],'passed')
            from result_store import ResultData
            saved=ResultData(d);self.assertEqual(saved['heavy_cm3'].shape[2],4)
            self.assertIn('SYNTHETIC H',saved.states[0]['label']);self.assertIn('e⁻',(Path(d)/'реакции_по_времени.csv').read_text(encoding='utf-8-sig'))

    def test_ion_mass_and_charge_are_actual_kernel_parameters(self):
        ez=np.array([1e7]);ep=np.array([1e7]);iz=np.array([2e5]);ip=np.array([3e5]);az=np.array([0.,np.pi])
        for mass in [1.008*1.66053906892e-24,6.94*1.66053906892e-24]:
            rel,cx,cy,cz=K.pair_geometry(ez,ep,iz,ip,az,mass)
            mu=ME*mass/(ME+mass)
            expected=mu*((ep[0]-ip[0])**2+(ez[0]-iz[0])**2)/(2*EV)
            self.assertAlmostEqual(rel[0,0,0]/expected,1.,places=14)
            one=K.pair_rates(rel,10.,2.,mass,1.);two=K.pair_rates(rel,10.,2.,mass,2.)
            np.testing.assert_allclose(two/one,4.,rtol=1e-14)


if __name__=='__main__':unittest.main()
