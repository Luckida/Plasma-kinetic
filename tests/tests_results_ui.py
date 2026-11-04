"""Real source forms and embedded result charts, exercised with a display."""
import json
from pathlib import Path
import tempfile
import tkinter as tk
import unittest
from PIL import Image
from experiment_editor import ExperimentEditor
from result_browser import ResultBrowser
from mixture_kinetics import ROOT
from tests_mixture import config,build


class NewSourceControls(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:cls.root=tk.Tk();cls.root.withdraw()
        except tk.TclError as exc:raise unittest.SkipTest(str(exc))

    @classmethod
    def tearDownClass(cls):cls.root.destroy()

    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);self.path=Path(self.temp.name)/'experiment.json';self.path.write_bytes((ROOT/'configs/gas_probe.json').read_bytes());self.editor=ExperimentEditor(self.root,self.path);self.addCleanup(self.editor.window.destroy);self.root.update()

    def test_density_and_velocity_save_without_a_flux_field(self):
        e=self.editor;s=e.sources_editor;s.vars['density'].set('12,5');s.vars['velocity'].set('-2e5');s.vars['id'].set('Пучок протонов');saved=e.save();row=saved['spatial_sources'][0]
        self.assertEqual(row['density_cm3'],12.5);self.assertEqual(row['velocity_cm_s'],-2e5);self.assertEqual(row['rate_cm2_s'],2.5e6);self.assertNotIn('rate',s.vars)
        self.assertEqual(len(e.save()['spatial_sources']),1)

    def test_source_with_a_new_element_extends_the_inventory(self):
        e=self.editor;s=e.sources_editor;s.vars['species'].set('Li+');s.vars['density'].set('1e8');s.vars['velocity'].set('1e6');saved=e.save();self.assertIn('Li II',saved['atomic_data']['species'])

    def test_normalised_light_uses_probe_and_displays_energy(self):
        e=self.editor;e.vars['pump_enabled'].set(True);e.vars['pump_use_probe'].set(True);e.vars['pump_scale'].set('.002');e.pump_energy(quiet=True);saved=e.save();r=saved['radiation']
        self.assertEqual(r['normalization'],'integral');self.assertEqual(r['spectrum_mode'],'probe');self.assertEqual(r['multiplier_W_cm2'],.002);self.assertIn('эВ',e.pump_info.get());self.assertIn('Дж/см²',e.pump_info.get())

    def test_disabled_internal_light_accepts_an_unavailable_file(self):
        e=self.editor;s=e.sources_editor;s.vars['mode'].set('Свет');s.enabled.set(False);s.use_probe.set(False);s.vars['file'].set('missing.txt');saved=e.save();self.assertFalse(saved['spatial_sources'][0]['enabled'])

    def test_independent_electron_and_ion_sources_can_both_be_saved(self):
        e=self.editor;s=e.sources_editor;s.vars['id'].set('Ионы');s.apply();s.new();s.vars['id'].set('Электроны');s.vars['mode'].set('Электроны');s.vars['temperature'].set('.1');s.vars['velocity'].set('1e7');s.apply();saved=e.save();self.assertEqual([row['mode'] for row in saved['spatial_sources']],['ions','electrons'])


class ResultWindow(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:cls.root=tk.Tk();cls.root.withdraw()
        except tk.TclError as exc:raise unittest.SkipTest(str(exc))

    @classmethod
    def tearDownClass(cls):cls.root.destroy()

    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);self.folder=Path(self.temp.name);c=config();c['processes']['electron_ionization']=True;c['probe']=dict(mode='flat',minimum_nm=400.,maximum_nm=800.,points=31)
        c['spatial_sources']=[dict(id='Пучок e',mode='electrons',density_cm3=10.,velocity_cm_s=1e7,position_cm=.005,start_s=.002,stop_s=.006,distribution=dict(Te_eV=.1))]
        model=build(c);model.solve();model.export(self.folder);self.browser=ResultBrowser(self.root,self.folder);self.addCleanup(self.browser.close);self.root.update()

    def choose(self,*keys):
        for key,var in self.browser.metric_vars.items():var.set(key in keys)

    def test_only_checked_graphs_and_selected_components_are_built(self):
        b=self.browser;self.choose('components');b.component_list.selection_clear(0,'end');b.component_list.selection_set(0);pages=b.build_sync();self.root.update()
        self.assertEqual(len(pages),1);self.assertEqual(len(pages[0].curves),1);self.assertEqual(pages[0].curves[0][0],b.component_keys[0]);self.assertEqual(len(b.preview.tabs()),1)

    def test_quick_component_choice_and_always_visible_log_controls(self):
        b=self.browser;b.controls.select(b.reactions_tab);self.root.update();b.quick_buttons['components'].invoke();b.vars['y_scale'].set('Логарифмический');self.root.update()
        self.assertEqual(b.options().metrics,['components']);self.assertTrue(b.scale_widgets['y_scale'].winfo_ismapped());self.assertEqual(b.options().y_scale,'log')
        pages=b.build_sync();self.assertTrue(pages);self.assertEqual(pages[0].figure.axes[0].get_yscale(),'log')

    def test_log_scales_start_marker_and_export_keep_large_text(self):
        b=self.browser;self.choose('components');b.vars['x_scale'].set('Логарифмический');b.vars['y_scale'].set('Логарифмический');b.vars['fontsize'].set('18');b.component_list.selection_clear(0,'end');b.component_list.selection_set(0)
        pages=b.build_sync();self.root.update();ax=pages[0].figure.axes[0];self.assertEqual(ax.get_xscale(),'log');self.assertEqual(ax.get_yscale(),'log');self.assertEqual(ax.xaxis.label.get_fontsize(),18)
        self.assertTrue(any(line.get_gid()=='pump_start' for line in ax.lines));b.window.geometry('1050x720');self.root.update();files=b.save_to(self.folder/'chosen')
        with Image.open(files[0]) as image:self.assertEqual(image.width,2200)

    def test_reaction_details_and_spatial_axis_are_selectable(self):
        b=self.browser;self.choose('rates');b.vars['axis'].set('По расстоянию');first=b.reaction_list.get_children()[0];b.reaction_list.selection_set(first);b.show_reaction();pages=b.build_sync();self.root.update()
        self.assertTrue(b.reaction_panel.detail.get('1.0','end').strip());self.assertEqual(pages[0].figure.axes[0].get_xlabel(),'Расстояние, см');self.assertEqual(len(pages[0].curves),1)

    def test_preserving_source_times_does_not_trigger_unselected_figures(self):
        b=self.browser;self.choose('mass','charge');pages=b.build_sync();self.root.update();self.assertEqual(len(pages),2);self.assertTrue(all(page.key.startswith(('mass_','charge_')) for page in pages))
