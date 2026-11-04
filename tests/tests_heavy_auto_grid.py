"""Numerical grid coverage; synthetic atomic coefficients only."""
from copy import deepcopy
from pathlib import Path
import tempfile
import unittest
import numpy as np
from particle_kinetics import KB,EV,GridOverflow
from phase_grid_setup import resolve_heavy_grid,AMU
from tests_mixture import config,build,synthetic_bank,ROOT


class HeavyAutoGrid(unittest.TestCase):
    def setUp(self):self.cfg=config();self.data,_=synthetic_bank();self.cfg['heavy_phase_grid']['energy_max_eV']='auto'

    def resolve(self):return resolve_heavy_grid(self.cfg,self.data,ROOT)

    def test_hot_gas_auto_grid_preserves_temperature_and_inputs(self):
        self.cfg['heavy_initial']={name:dict(Ti_K=300000.,match_moments=True) for name in ['H','He']};self.cfg['phase_grid']['energy_max_eV']='auto';before=deepcopy(self.cfg)
        m=build(self.cfg);self.assertGreater(next(iter(m.gh.values())).e[-1],8*KB*300000/EV);self.assertEqual(self.cfg,before)
        for element,ident in [('H','H0'),('He','He0')]:
            moments=m.gh[element].moments(m.P[0,m.index[ident]]);self.assertAlmostEqual(moments['T_moment_eV']*EV/KB/300000.,1.,places=10)

    def test_manual_insufficient_grid_is_rejected_with_actual_bounds(self):
        self.cfg['heavy_phase_grid']['energy_max_eV']=50.;self.cfg['heavy_initial']={name:dict(Ti_K=300000.) for name in ['H','He']}
        with self.assertRaisesRegex(GridOverflow,'Emax=50.*206.816.*300000.*auto'):build(self.cfg)

    def test_hot_source_is_included_and_disabled_source_is_ignored(self):
        self.cfg['spatial_sources']=[dict(id='hot',mode='ions',enabled=True,species='H+',distribution=dict(Ti_K=1e6)),dict(id='off',mode='ions',enabled=False,species='He+',distribution=dict(Ti_K=1e9))]
        maximum=self.resolve()['energy_max_eV'];self.assertAlmostEqual(maximum,1.1*8*KB*1e6/EV)

    def test_beam_bound_uses_actual_species_mass_and_thermal_speed(self):
        self.cfg['spatial_sources']=[dict(mode='ions',species='He+',distribution=dict(beams=[dict(fraction=1.,Ti_K=300.,uz_cm_s=1e7)]))]
        mass=next(s['mass_u'] for s in self.data.states if s['element']=='He')*AMU;bulk=mass*1e14/(2*EV);thermal=300*KB/EV
        self.assertAlmostEqual(self.resolve()['energy_max_eV'],1.1*(np.sqrt(bulk)+np.sqrt(8*thermal))**2)

    def test_energy_file_bounds_expand_the_automatic_grid(self):
        with tempfile.TemporaryDirectory() as folder:
            file=Path(folder)/'iedf.csv';file.write_text('E;P\n0.00000001;0\n150;1\n300;0\n')
            self.cfg['heavy_initial']['H0']=dict(kind='energy_table',energy_distribution_file=str(file));grid=self.resolve();self.assertLessEqual(grid['energy_min_eV'],1e-8);self.assertGreaterEqual(grid['energy_max_eV'],300.)

    def test_phase_file_domain_is_not_extended_or_extrapolated(self):
        with tempfile.TemporaryDirectory() as folder:
            file=Path(folder)/'phase.csv';file.write_text('energy_eV;mu;probability_per_eV_mu\n0.000001;-1;0.5000005000005\n0.000001;1;0.5000005000005\n1;-1;0.5000005000005\n1;1;0.5000005000005\n')
            self.cfg['heavy_initial']['H0']=dict(kind='phase_table',distribution_file=str(file));self.assertEqual(self.resolve()['energy_max_eV'],1.)
            self.cfg['heavy_initial']['He']={'Ti_K':300000.}
            with self.assertRaisesRegex(GridOverflow,'Фазовые файлы'):self.resolve()

    def test_reservoir_temperature_is_included(self):
        self.cfg['particle_boundary']=dict(kind='reservoir',left_heavy={'H+':dict(density_cm3=1.,Ti_K=900000.)})
        self.assertGreaterEqual(self.resolve()['energy_max_eV'],8*KB*900000/EV)

    def test_manual_valid_grid_is_kept_exactly(self):
        self.cfg['heavy_phase_grid']['energy_max_eV']=300.;self.assertEqual(self.resolve(),self.cfg['heavy_phase_grid'])

    def test_auto_keeps_the_previous_normal_energy_range_for_reaction_products(self):
        self.assertEqual(self.resolve()['energy_max_eV'],50.)

    def test_one_ulp_progress_advance_does_not_blow_up_eta(self):
        from progress import Progress,parse_event
        clock=[0.];events=[];progress=Progress(mode='json',sink=events.append,now=lambda:clock[0]);progress.start('test')
        for fraction in [.1,.2,.25]:clock[0]+=1.;progress.update(fraction,1.,force=True)
        clock[0]+=.1;progress.update(np.nextafter(.25,1.),1.,force=True)
        eta=parse_event(events[-1])['eta_s'];self.assertIsNotNone(eta);self.assertLess(eta,30.)
        clock[0]+=1.;progress.update(.3,1.,force=True);self.assertLess(parse_event(events[-1])['eta_s'],30.)


class HeavyGridEditor(unittest.TestCase):
    def test_auto_heavy_bound_roundtrips_through_the_actual_editor(self):
        import json,tkinter as tk
        from experiment_editor import ExperimentEditor
        try:root=tk.Tk();root.withdraw()
        except tk.TclError as exc:raise unittest.SkipTest(str(exc))
        try:
            with tempfile.TemporaryDirectory() as folder:
                file=Path(folder)/'config.json';file.write_bytes((ROOT/'configs/gas_probe.json').read_bytes());editor=ExperimentEditor(root,file)
                try:
                    editor.vars['temperature'].set('300000');editor.distributions_editor.grid_vars['h_max'].set('auto');saved=editor.save()
                    self.assertEqual(saved['heavy_phase_grid']['energy_max_eV'],'auto');self.assertEqual(saved['gas']['temperature_K'],300000.);self.assertEqual(json.loads(file.read_text())['heavy_phase_grid']['energy_max_eV'],'auto')
                finally:editor.window.destroy()
        finally:root.destroy()
