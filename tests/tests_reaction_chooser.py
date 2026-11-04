"""User flows missing in 2.6: live suggestions, identities and plot routing."""
from pathlib import Path
import tempfile
import tkinter as tk
import unittest
import numpy as np
from reaction_catalog import ReactionCatalog
from result_browser import ResultBrowser
from result_plots import PlotOptions,build_plots,top_history
from reaction_fixture import FixtureData


class ReactionCatalogueTests(unittest.TestCase):
    def setUp(self):self.data=FixtureData();self.catalog=ReactionCatalog(self.data)

    def ids(self,indices):return {self.data.reactions[k]['id'] for k in indices}

    def test_aliases_accept_charge_and_metastable_notation(self):
        for query in ['He(2^3S)','He 2³S','HeI 2 3S']:
            self.assertIn('PI_H',self.ids(self.catalog.matching(search=query)))
        self.assertIn('CX',self.ids(self.catalog.matching(search='H+')));self.assertNotIn('PUMP',self.ids(self.catalog.matching(search='H+')))

    def test_zero_rate_channels_are_available_until_explicitly_filtered(self):
        self.assertIn('PI_Li',self.ids(self.catalog.matching()));self.assertNotIn('PI_Li',self.ids(self.catalog.matching(nonzero=True)))

    def test_two_reagents_mean_two_input_particles_not_a_product(self):
        self.assertEqual(self.ids(self.catalog.matching(partner_a='level:HeS',partner_b='level:H0')),{'PI_H'})
        self.assertNotIn('A_S',self.ids(self.catalog.matching(partner_a='level:HeS',partner_b='photon')))

    def test_equal_reagent_choices_require_stoichiometric_multiplicity_two(self):
        self.assertEqual(self.ids(self.catalog.matching(partner_a='level:HeS',partner_b='level:HeS')),{'SELF'})

    def test_compatible_partner_choices_are_built_from_the_network(self):
        choices={c.key for c in self.catalog.choices_for_partners('level:HeS')}
        self.assertIn('level:H0',choices);self.assertIn('level:Li0',choices);self.assertNotIn('photon',choices)

    def test_component_participation_shows_internal_excitation_and_resonant_exchange(self):
        involved=self.ids(self.catalog.matching(target='component:He I',direction='participation'));net=self.ids(self.catalog.matching(target='component:He I',direction='both'))
        self.assertIn('EX',involved);self.assertNotIn('EX',net)
        self.assertIn('CX',self.ids(self.catalog.matching(target='component:H I',direction='participation')));self.assertNotIn('CX',self.ids(self.catalog.matching(target='component:H I',direction='source')))

    def test_suggestions_follow_cell_time_and_interval(self):
        all_rows=self.catalog.matching();early=self.catalog.ranked(all_rows,context={'cell':0,'time_max':0.});late=self.catalog.ranked(all_rows,context={'cell':0,'time_min':1.})
        self.assertEqual(self.data.reactions[early[0]]['id'],'EX');self.assertEqual(self.data.reactions[late[0]]['id'],'DX')
        right=self.catalog.ranked(all_rows,context={'cell':1});self.assertEqual(self.data.reactions[right[0]]['id'],'PUMP')

    def test_filtered_top_five_cannot_include_other_partners(self):
        options=PlotOptions(metrics=['top5'],reaction_filter=['PI_H'],target='level:HeS',direction='sink');history=top_history(self.data,options)
        self.assertEqual([self.data.reactions[k]['id'] for k,_ in history],['PI_H'])

    def test_both_population_directions_produce_separate_labelled_plots(self):
        pages=build_plots(self.data,PlotOptions(metrics=['top5'],target='level:HeS',direction='both'));titles=' '.join(p.title for p in pages)
        self.assertIn('заселение',titles);self.assertIn('расселение',titles)

    def test_zero_rates_are_not_misreported_as_positive_leaders(self):
        data=FixtureData(all_zero=True);pages=build_plots(data,PlotOptions(metrics=['top5']));self.assertEqual(pages[0].curves,[])
        self.assertTrue(any('Нет ненулевых реакций' in text.get_text() for text in pages[0].figure.axes[0].texts))

    def test_a_time_range_without_samples_has_a_clear_error(self):
        with self.assertRaisesRegex(ValueError,'нет сохранённых'):self.catalog.context_rates(time_min=2.)

    def test_zero_reaction_explains_a_saved_temperature_domain_problem(self):
        self.data.meta['quality']['data_issues']=['PI_Li: no unique k(T) interval (диапазона) for the occupied distributions; see chemical_domain_audit.csv']
        details=self.catalog.describe(self.catalog.by_id['PI_Li']);self.assertIn('нет применимого опубликованного диапазона',details);self.assertIn('Скорость равна нулю',details)

    def test_physical_code_locations_are_real_ast_definitions(self):
        from process_map import process_locations
        for kind in ['excitation','photoionization','penning','charge_exchange']:
            for location in process_locations(kind):
                source=Path(location['path']).read_text().splitlines();self.assertIn('def '+location['function'].split('.')[-1],source[location['line']-1]);self.assertTrue(location['code'])


class ReactionWindowTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:cls.root=tk.Tk();cls.root.withdraw()
        except tk.TclError as exc:raise unittest.SkipTest(str(exc))

    @classmethod
    def tearDownClass(cls):cls.root.destroy()

    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);self.data=FixtureData();self.data.save(self.temp.name);self.browser=ResultBrowser(self.root,self.temp.name);self.addCleanup(self.browser.close);self.root.update();self.panel=self.browser.reaction_panel

    def label(self,key,target=False):
        mapping=self.panel.target_labels if target else self.panel.choice_labels
        return next(label for label,value in mapping.items() if value==key)

    def test_browsing_reaction_tab_preserves_the_checked_concentrations(self):
        b=self.browser;b.controls.select(b.reactions_tab);self.root.update();self.assertEqual(b.options().metrics,['components']);pages=b.build_sync();self.assertTrue(pages);self.assertTrue(all(p.key.startswith('components') for p in pages));self.assertEqual(b.build_button['text'],'Построить выбранное')

    def test_partner_pair_automatically_selects_matching_rows(self):
        p=self.panel;p.vars['partner_a'].set(self.label('level:HeS'));p.refresh();p.vars['partner_b'].set(self.label('level:H0'));p.refresh();p.vars['mode'].set('Скорости выбранных реакций')
        self.browser.select_plot('rates');request=p.request();self.assertEqual(request.selected_ids,['PI_H']);self.assertEqual(request.matching_ids,['PI_H']);self.assertIn('Реагенты:',p.detail.get('1.0','end'));self.assertIn('Источник данных:',p.detail.get('1.0','end'))
        self.assertEqual(len(self.browser.build_sync()[0].curves),1)

    def test_searching_beyond_global_top_five_still_gets_an_automatic_selection(self):
        p=self.panel;p.vars['search'].set('PI_Li');p.refresh();request=p.request('rates');self.assertEqual(request.selected_ids,['PI_Li']);self.assertIn('0',str(p.tree.item(str(p.catalog.by_id['PI_Li']),'values')))

    def test_selecting_a_component_proposes_participation_and_levels_propose_both(self):
        p=self.panel;p.vars['target'].set(self.label('component:He I',True));p.refresh();self.assertEqual(p.direction,'participation');self.assertIn('EX',p.matching_ids)
        p.vars['target'].set(self.label('level:HeS',True));p.refresh();self.assertEqual(p.direction,'both')

    def test_filtering_preserves_reaction_identity_instead_of_reusing_list_positions(self):
        p=self.panel;k=p.catalog.by_id['PI_H'];p.tree.selection_set(str(k));p.vars['search'].set('He(2^3S)');p.refresh();self.assertEqual(p.selected_ids,['PI_H'])

    def test_an_empty_network_reports_absence_and_disables_its_plot_button(self):
        folder=Path(self.temp.name)/'empty';folder.mkdir();FixtureData(empty=True).save(folder);b=ResultBrowser(self.root,folder)
        try:self.root.update();self.assertEqual(str(b.reaction_panel.build_button['state']),'disabled');self.assertIn('нет реакций',b.reaction_panel.detail.get('1.0','end'))
        finally:b.close()

    def test_component_tab_keeps_its_original_plot_selection(self):
        b=self.browser;b.controls.select(b.reactions_tab);self.root.update();b.controls.select(b.controls.tabs()[0]);self.root.update();self.assertEqual(b.options().metrics,['components'])

    def test_selected_process_opens_its_actual_code_in_the_application(self):
        window=self.panel.show_code()
        try:self.root.update();self.assertIsNotNone(window);self.assertIn('исходный код',window.title())
        finally:
            if window is not None:window.destroy()

    def test_real_build_button_completes_with_the_selected_reaction(self):
        import time
        p=self.panel;p.vars['search'].set('PI_H');p.vars['mode'].set('Скорости выбранных реакций');p.refresh();self.browser.build_reactions()
        end=time.monotonic()+8
        while self.browser.busy and time.monotonic()<end:self.root.update();time.sleep(.02)
        self.assertFalse(self.browser.busy);self.assertTrue(self.browser.pages);self.assertTrue(all(page.key.startswith('rates') for page in self.browser.pages))
