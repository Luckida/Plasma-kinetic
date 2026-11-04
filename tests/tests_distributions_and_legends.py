"""Regressions for visible inventories, labelled plots and particle input.

All fixture coefficients here are synthetic, never scientific source data.
"""
import json
from pathlib import Path
import tempfile
import tkinter as tk
import unittest
import numpy as np
from energy_input import read_energy_table,energy_moments,energy_phase
from particle_kinetics import PhaseGrid,GridOverflow,KB,EV
from result_plots import PlotOptions,build_plots
from reaction_catalog import ReactionCatalog
from reaction_fixture import FixtureData
from tests_mixture import config,build,ROOT


class DistributionProjection(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);self.file=Path(self.temp.name)/'energy.txt';self.file.write_text('E_eV P\n0.000001 0\n0.35 2\n1.2 0.6\n2.7 0\n');self.grid=PhaseGrid(dict(energy_min_eV=1e-7,energy_max_eV=10.,energy_points=32,angular_points=4))

    def test_arbitrary_energy_shape_preserves_number_energy_and_isotropy(self):
        x,y=read_energy_table(self.file);_,mean=energy_moments(x,y);phase=energy_phase(self.grid,self.file,17.);moments=self.grid.moments(phase)
        self.assertAlmostEqual(phase.sum(),17.,places=12);self.assertAlmostEqual(moments['energy_eV_cm3']/17.,mean,places=12);self.assertLess(abs(moments['uz_cm_s']),1e-6)
        np.testing.assert_allclose(phase,phase[:,::-1],rtol=1e-14,atol=1e-20)

    def test_eepf_conversion_is_explicit_and_normalised(self):
        x,pdf=read_energy_table(self.file);_,eepf=read_energy_table(self.file,'eepf');self.assertNotAlmostEqual(energy_moments(x,pdf)[1],energy_moments(x,eepf)[1]);self.assertAlmostEqual(np.trapezoid(eepf,x),1.)

    def test_uncovered_energies_are_not_silently_clipped(self):
        self.file.write_text('1 0\n20 1\n40 0\n')
        with self.assertRaises(GridOverflow):energy_phase(self.grid,self.file,1.)

    def test_maxwellian_default_has_requested_temperature_on_coarse_grid(self):
        grid=PhaseGrid(dict(energy_min_eV=1e-6,energy_max_eV=50.,energy_points=16,angular_points=4),4*1.66053906892e-24)
        phase=grid.ion_initial(dict(kind='maxwellian',Ti_K=300.,match_moments=True),11.)
        self.assertAlmostEqual(grid.moments(phase)['T_moment_eV']*EV/KB,300.,places=8)

    def test_electron_maxwellian_matches_input_mean_without_hot_tail(self):
        phase=self.grid.initial(dict(kind='maxwellian',Te_eV=.1,match_moments=True),11.)
        self.assertAlmostEqual(self.grid.moments(phase)['energy_eV_cm3']/11.,.15,places=12)

    def test_initial_electron_file_extends_auto_energy_grid_and_is_used(self):
        self.file.write_text('20 0\n25 1\n40 0\n');cfg=config();cfg['phase_grid']['energy_max_eV']='auto';cfg['electron_initial']=dict(kind='energy_table',energy_distribution_file=str(self.file));m=build(cfg)
        self.assertGreaterEqual(m.ge.e[-1],40.);self.assertGreater(m.ge.moments(m.F[0])['energy_eV_cm3']/m.F[0].sum(),20.)

    def test_each_level_can_override_the_common_heavy_distribution(self):
        cfg=config();cfg['heavy_phase_grid']['energy_max_eV']=5.;cfg['heavy_initial']['H0']=dict(kind='energy_table',energy_distribution_file=str(self.file));cfg['heavy_initial']['H1']=dict(Ti_K=600.)
        # Exact bank IDs are accepted, including IDs without the production ':' convention.
        m=build(cfg);mean=m.gh['H'].moments(m.P[0,m.index['H0']])['energy_eV_cm3']/m.P[0,m.index['H0']].sum();self.assertAlmostEqual(mean,energy_moments(*read_energy_table(self.file))[1],places=12)


class PlotLegendRegression(unittest.TestCase):
    def test_long_reaction_legends_fit_inside_png_at_large_font(self):
        data=FixtureData();ids=[r['id'] for r in data.reactions]
        for row in data.reactions:row['reaction']+=' · длинное название состояния '*3
        pages=build_plots(data,PlotOptions(metrics=['rates'],reactions=ids,fontsize=22,curves_per_figure=8))
        for page in pages:
            fig=page.figure;fig.canvas.draw();legend=fig.legends[0];box=legend.get_window_extent();canvas=fig.bbox
            self.assertEqual(len(legend.get_texts()),len(page.curves));self.assertGreaterEqual(box.x0,0);self.assertGreaterEqual(box.y0,0);self.assertLessEqual(box.x1,canvas.x1);self.assertLessEqual(box.y1,canvas.y1)
            self.assertLess(box.y1,fig.axes[0].get_window_extent().y0)

    def test_zero_curve_on_log_axis_still_has_an_explanation(self):
        data=FixtureData(all_zero=True);page=build_plots(data,PlotOptions(metrics=['rates'],reactions=['PI_H'],y_scale='log'))[0]
        texts=[text.get_text() for text in page.figure.legends[0].get_texts()];caption=' '.join(' '.join(texts).split());self.assertIn('PI_H',caption);self.assertIn('нет положительных',caption)

    def test_all_selected_reactions_get_a_description_not_just_first(self):
        data=FixtureData();catalog=ReactionCatalog(data);text=catalog.describe_many(range(len(data.reactions)))
        for reaction in data.reactions:self.assertIn('ID: '+reaction['id']+'\n',text)
        self.assertEqual(text.count('Источник данных:'),len(data.reactions))


class DistributionWindows(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:cls.root=tk.Tk();cls.root.withdraw()
        except tk.TclError as exc:raise unittest.SkipTest(str(exc))

    @classmethod
    def tearDownClass(cls):cls.root.destroy()

    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);self.folder=Path(self.temp.name)

    def browser(self,data):
        from result_browser import ResultBrowser
        data.save(self.folder);browser=ResultBrowser(self.root,self.folder);self.addCleanup(browser.close);self.root.update();return browser

    def editor(self):
        from experiment_editor import ExperimentEditor
        file=self.folder/'config.json';file.write_bytes((ROOT/'configs/gas_probe.json').read_bytes());editor=ExperimentEditor(self.root,file);self.addCleanup(editor.window.destroy);self.root.update();return editor

    def test_zero_population_levels_are_visible_and_selectable_by_default(self):
        data=FixtureData();data.arrays['n_cm3'][:,:,1:]=0.;browser=self.browser(data)
        self.assertFalse(browser.nonzero.get());self.assertEqual(len(browser.level_keys),len(data.states));browser.select_all_components();self.assertEqual(len(browser.level_list.curselection()),len(data.states))

    def test_both_reagent_menus_keep_all_levels_available(self):
        data=FixtureData();browser=self.browser(data);panel=browser.reaction_panel
        for key in ('partner_a','partner_b'):
            choices=panel.combos[key]['values']
            for state in data.states:self.assertTrue(any(panel.choice_labels[label]=='level:'+state['id'] for label in choices))

    def test_select_all_describes_every_reaction_and_reference_opens(self):
        data=FixtureData();browser=self.browser(data);panel=browser.reaction_panel;panel.select_all();text=panel.detail.get('1.0','end')
        for reaction in data.reactions:self.assertIn('ID: '+reaction['id']+'\n',text)
        window=panel.show_catalogue()
        try:self.root.update();self.assertIn('все доступные',window.title())
        finally:window.destroy()

    def test_electron_file_survives_default_mode_save_and_reopen(self):
        e=self.editor();file=self.folder/'eedf.txt';file.write_text('E P\n0.000001 0\n1 1\n4 0\n');form=e.distributions_editor.electron;form.kind.set('Файл P(E) / EEPF · изотропный');form.file.set(str(file));saved=e.save()
        self.assertEqual(saved['electron_initial']['kind'],'energy_table');self.assertEqual(saved['electron_initial']['energy_distribution_file'],str(file));self.assertEqual(json.loads(e.path.read_text())['electron_initial']['kind'],'energy_table')

    def test_ion_file_and_per_species_override_survive_composition_translation(self):
        e=self.editor();file=self.folder/'iedf.txt';file.write_text('E P\n0.000001 0\n0.1 1\n1 0\n');editor=e.distributions_editor;editor.target.set('He+');editor.switch();editor.use_default.set(False);editor.heavy.kind.set('Файл P(E) / EEPF · изотропный');editor.heavy.file.set(str(file));saved=e.save()
        self.assertEqual(saved['heavy_initial']['He II']['kind'],'energy_table');self.assertEqual(saved['heavy_initial']['He']['kind'],'maxwellian')

    def test_source_defaults_expose_and_save_particle_distribution(self):
        e=self.editor();source=e.sources_editor;source.vars['density'].set('1e7');source.vars['velocity'].set('1e6');saved=e.save();distribution=saved['spatial_sources'][0]['distribution']
        self.assertEqual(distribution['kind'],'maxwellian');self.assertTrue(distribution['match_moments']);self.assertEqual(distribution['Ti_K'],300.)

    def test_save_button_stays_inside_window_on_every_tab(self):
        e=self.editor();e.window.geometry('1000x760')
        for tab in e.book.tabs():
            e.book.select(tab);self.root.update();button=e.save_button
            self.assertTrue(button.winfo_viewable());self.assertLessEqual(button.winfo_rooty()+button.winfo_height(),e.window.winfo_rooty()+e.window.winfo_height())
