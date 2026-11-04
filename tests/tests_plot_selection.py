"""Real result-option methods without a display, plus numeric plot checks."""
import tempfile
from types import SimpleNamespace
import unittest
import warnings
from unittest.mock import Mock
import tkinter as tk
import numpy as np

from result_browser import ResultBrowser
from result_store import ResultData
from result_plots import build_plots
from tests_mixture import config,build


class PlotSelection(unittest.TestCase):
    def setUp(self):
        self.interp=tk.Tcl();self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        m=build(config());m.solve();m.export(self.temp.name);self.data=ResultData(self.temp.name)
        b=object.__new__(ResultBrowser);b.data=self.data;b.vars={}
        defaults={'axis':'По времени','cell':'Среднее по слою','time':'последний','x_scale':'Линейный','y_scale':'Линейный','fontsize':'16','curves':'5','time_min':'','time_max':'','wave_min':'','wave_max':'','probe_wave':'1083.3'}
        b.vars={k:tk.StringVar(self.interp,value=v) for k,v in defaults.items()};b.cell_values=['Среднее по слою'];b.time_values=['последний']
        b.metric_vars={k:tk.BooleanVar(self.interp,value=k=='components') for k in ['components','levels','rates','top5']}
        b.component_keys=list(self.data.meta['components']);b.level_keys=[s['id'] for s in self.data.states]
        b.component_list=Mock();b.component_list.curselection.return_value=[0];b.level_list=Mock();b.level_list.curselection.return_value=[0]
        b.controls=Mock();b.reactions_tab='reaction';b.components_tab='component';b.graphs_tab='graphs';b.controls.select.return_value='reaction'
        b.metric_widgets={k:Mock() for k in b.metric_vars};b.status=tk.StringVar(self.interp)
        for key,widget in b.metric_widgets.items():widget.cget.return_value=key
        b.reaction_panel=Mock();b.reaction_panel.vars={'mode':tk.StringVar(self.interp,value='Top-5 по времени')}
        b.reaction_panel.request.return_value=SimpleNamespace(selected_ids=[],matching_ids=[],target='all',direction='source')
        b.mark=tk.BooleanVar(self.interp,value=True);b.fractions=tk.BooleanVar(self.interp,value=False);self.browser=b

    def test_browsing_reaction_tab_does_not_override_checked_concentrations(self):
        b=self.browser;self.assertEqual(b.options().metrics,['components']);b.reaction_panel.request.assert_not_called()

    def test_quick_choice_selects_populations_and_retains_scales(self):
        b=self.browser;b.vars['y_scale'].set('Логарифмический');b.select_plot('levels')
        self.assertEqual(b.options().metrics,['levels']);self.assertEqual(b.options().y_scale,'log');b.controls.select.assert_called_with('component')

    def test_concentrations_are_the_saved_populations_with_log_axes(self):
        b=self.browser;b.vars['x_scale'].set('Логарифмический');b.vars['y_scale'].set('Логарифмический')
        with warnings.catch_warnings():
            warnings.simplefilter('error',UserWarning);pages=build_plots(self.data,b.options())
        self.assertTrue(pages)
        self.assertEqual(pages[0].figure.axes[0].get_xscale(),'log');self.assertEqual(pages[0].figure.axes[0].get_yscale(),'log')
        key=b.component_keys[0];expected=self.data['n_cm3'][:,:,self.data.meta['components'][key]].sum(axis=2).mean(axis=1)
        np.testing.assert_allclose(pages[0].curves[0][1],expected)

    def test_multiple_checked_quantities_survive_tab_switches(self):
        b=self.browser;b.metric_vars['levels'].set(True);self.assertEqual(b.options().metrics,['components','levels'])

    def test_concentrations_and_rates_can_be_requested_together(self):
        b=self.browser;b.metric_vars['rates'].set(True);self.assertEqual(b.options().metrics,['components','rates']);b.reaction_panel.request.assert_called_once_with('rates')
