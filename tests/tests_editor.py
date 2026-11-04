"""Actual Tk widgets: edit/save composition, import and disabled-source flow."""
import json
from pathlib import Path
import tempfile
import tkinter as tk
import unittest
from experiment_editor import ExperimentEditor,MODES
from mixture_kinetics import ROOT


class GasSpectrumEditor(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:cls.root=tk.Tk();cls.root.withdraw()
        except tk.TclError as exc:raise unittest.SkipTest('Для GUI-тестов нужен дисплей: '+str(exc))

    @classmethod
    def tearDownClass(cls):cls.root.destroy()

    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);self.path=Path(self.temp.name)/'experiment.json'
        self.path.write_bytes((ROOT/'configs/gas_probe.json').read_bytes());self.editor=ExperimentEditor(self.root,self.path);self.addCleanup(self.editor.window.destroy);self.root.update()

    def test_widgets_edit_and_save_concentrations_without_json(self):
        e=self.editor;e.vars['species'].set('Li');e.vars['amount'].set('2,5e8');e.apply_row();cfg=e.save()
        self.assertEqual(cfg['initial_cm3']['Li'],2.5e8);self.assertIn('Li I',cfg['atomic_data']['species']);self.assertEqual(json.loads(self.path.read_text())['initial_cm3']['Li'],2.5e8)

    def test_both_angular_grids_are_visible_and_saved(self):
        e=self.editor;e.distributions_editor.grid_vars['e_angles'].set('8');e.distributions_editor.grid_vars['h_angles'].set('6');saved=e.save()
        self.assertEqual(saved['phase_grid']['angular_points'],8);self.assertEqual(saved['heavy_phase_grid']['angular_points'],6)

    def test_invalid_hidden_legacy_angle_cannot_be_saved(self):
        e=self.editor;before=self.path.read_bytes();e.distributions_editor.grid_vars['e_angles'].set('1')
        with self.assertRaisesRegex(ValueError,'angular_points=1'):e.save()
        self.assertEqual(self.path.read_bytes(),before)

    def test_pending_row_edit_is_saved_and_deleted_row_stays_deleted(self):
        e=self.editor;item=e.tree.get_children()[0];e.tree.selection_set(item);e.select_row();name=e.vars['species'].get();e.vars['amount'].set('7e8');cfg=e.save();self.assertEqual(cfg['initial_cm3'][name],7e8)
        e.tree.selection_set(item);e.remove_row();cfg=e.save();self.assertNotIn(name,cfg['initial_cm3'])

    def test_fraction_mode_converts_densities_and_normalizes(self):
        e=self.editor;before=sum(json.loads(self.path.read_text())['initial_cm3'].values());e.vars['gas_mode'].set('Доли, % + полная концентрация');e.change_mode();e.normalize();cfg=e.save();self.assertEqual(cfg['gas']['mode'],'fractions');self.assertAlmostEqual(sum(cfg['initial_cm3'].values())/before,1.,places=10)

    def test_plain_spectrum_preview_and_file_mode_save(self):
        p=Path(self.temp.name)/'measured.dat';p.write_text('1084 2\n1082 1\n');e=self.editor;e.vars['probe_file'].set(str(p));e.vars['probe_mode'].set('file');s=e.preview('probe');self.assertEqual(s.values,(1.,2.));cfg=e.save();self.assertEqual(cfg['probe']['spectrum_file'],str(p))

    def test_disabled_bad_light_files_do_not_prevent_saving(self):
        e=self.editor;e.vars['probe_mode'].set('off');e.vars['pump_enabled'].set(False);e.vars['probe_file'].set('missing.txt');e.vars['pump_file'].set('missing.csv');cfg=e.save();self.assertFalse(cfg['radiation']['enabled']);self.assertEqual(cfg['probe']['mode'],'off')

    def test_invalid_spectrum_error_keeps_original_file_intact(self):
        e=self.editor;before=self.path.read_bytes();e.vars['probe_mode'].set('file');e.vars['probe_file'].set('missing.txt')
        with self.assertRaisesRegex(ValueError,'не найден'):e.save()
        self.assertEqual(self.path.read_bytes(),before)

    def test_ordinary_mode_clears_demo_injection_bias_and_estimates(self):
        cfg=json.loads(self.path.read_text());cfg['spatial_sources']=[dict(mode='electrons')];cfg['electron_initial'].update(angular_bias=.3,uz_cm_s=1e6);cfg['chemistry']['allow_estimated_rates']=True;cfg['collisions']['ee_bgk']=True;self.path.write_text(json.dumps(cfg));saved=self.editor.save()
        self.assertEqual(saved['spatial_sources'],[]);self.assertNotIn('angular_bias',saved['electron_initial']);self.assertFalse(saved['chemistry']['allow_estimated_rates']);self.assertFalse(saved['collisions']['ee_bgk'])

    def test_window_layout_keeps_save_button_visible(self):
        e=self.editor;self.root.update_idletasks();self.assertGreaterEqual(e.window.winfo_width(),850)
        buttons=[]
        def visit(widget):
            if widget.winfo_class()=='TButton' and widget.cget('text')=='Сохранить':buttons.append(widget)
            for child in widget.winfo_children():visit(child)
        visit(e.window);self.assertEqual(len(buttons),1);b=buttons[0];self.assertLess(b.winfo_rooty()+b.winfo_height(),e.window.winfo_rooty()+e.window.winfo_height()+1)

    def test_inactive_modes_disable_unnecessary_controls(self):
        e=self.editor;self.assertEqual(str(e.widgets['pump_file']['state']),'disabled');self.assertEqual(str(e.widgets['total']['state']),'disabled');self.assertEqual(str(e.widgets['ne']['state']),'disabled')
        e.vars['probe_mode'].set('file');self.assertEqual(str(e.widgets['probe_min']['state']),'disabled');self.assertEqual(str(e.widgets['probe_file']['state']),'normal')
        e.vars['pump_enabled'].set(True);e.vars['pump_use_probe'].set(False);self.assertEqual(str(e.widgets['pump_file']['state']),'normal')
