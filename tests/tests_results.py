"""Compact outputs, exact beam fluence, calibrated light and selected plots."""
import csv
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from light_input import light_spectrum,particle_inputs,integral
from result_store import ResultData,FILE_STORE,top_indices
from result_plots import PlotOptions,build_plots,save_pages,top_history
from tests_mixture import config,build


class InputSources(unittest.TestCase):
    def test_density_and_signed_velocity_define_column_rate_and_distribution(self):
        for velocity in [-2e5,2e5]:
            c=config();c['spatial_sources']=[dict(id='beam',mode='ions',species='H+',density_cm3=12.5,velocity_cm_s=velocity,position_cm=.005,start_s=.002,stop_s=.006,distribution=dict(Ti_K=300.))]
            m=build(c);row=m.sources.items[0];self.assertEqual(row['rate_cm2_s'],12.5*abs(velocity))
            moments=m.gh['H'].moments(row['phase']);self.assertAlmostEqual(moments['uz_cm_s']/velocity,1.,places=10)
            m.solve();self.assertAlmostEqual(m.sources.column[0]/(12.5*abs(velocity)*.004),1.,places=12)
            self.assertLess(abs(m.snapshots[-1]['charge'])/m.initial_nuclei.sum(),1e-12)

    def test_conflicting_old_flux_and_new_density_cannot_silently_override(self):
        with self.assertRaisesRegex(ValueError,'Не согласованы'):particle_inputs(dict(mode='ions',density_cm3=2.,velocity_cm_s=3.,rate_cm2_s=20.))

    def test_zero_speed_adds_no_particles(self):
        c=config();c['spatial_sources']=[dict(id='still',mode='electrons',density_cm3=1e10,velocity_cm_s=0.,position_cm=.005,distribution=dict(Te_eV=.1))];m=build(c);m.sources.apply(0.,.001);self.assertEqual(m.sources.column[0],0.)

    def test_normalised_spectral_shape_has_exact_energy_and_photon_receipt(self):
        c=dict(normalization='integral',spectrum_mode='flat',minimum_nm=400.,maximum_nm=800.,points=501,multiplier_W_cm2=2e-4,mu=.4)
        data=light_spectrum(c,Path('.'));power=integral(data.wavelength_A,data.values)*.4*1e-7
        self.assertAlmostEqual(power/2e-4,1.,places=12);meta=data.metadata
        self.assertLess(meta['minimum_photon_eV'],meta['mean_photon_eV']);self.assertLess(meta['mean_photon_eV'],meta['maximum_photon_eV'])
        self.assertAlmostEqual(meta['photon_flux_cm2_s']*meta['mean_photon_eV']*1.602176634e-19/(2e-4),1.,places=12)

    def test_relative_file_amplitude_cancels_during_normalisation(self):
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)/'shape.txt';cfg=dict(normalization='integral',spectrum_mode='file',spectrum_file=str(p),multiplier_W_cm2=.002,spectrum_options={'x_unit':'nm','y_kind':'relative'})
            p.write_text('400 1\n500 2\n800 1\n');a=light_spectrum(cfg,folder)
            p.write_text('400 10\n500 20\n800 10\n');b=light_spectrum(cfg,folder);np.testing.assert_allclose(a.values,b.values,rtol=1e-14)

    def test_normalised_light_rejects_negative_shape_or_zero_integral(self):
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)/'shape.txt';cfg=dict(normalization='integral',spectrum_mode='file',spectrum_file=str(p),multiplier_W_cm2=1.)
            for text in ['400 -1\n800 2','400 0\n800 0']:
                p.write_text(text)
                with self.assertRaises(ValueError):light_spectrum(cfg,folder)

    def test_boundary_and_internal_light_powers_survive_the_photon_grid(self):
        light=dict(normalization='integral',spectrum_mode='flat',minimum_nm=400.,maximum_nm=800.,points=31,multiplier_W_cm2=3e-6)
        c=config();c['radiation'].update(light,enabled=True,beam_mu=.6,start_s=.002,stop_s=.006);m=build(c)
        self.assertAlmostEqual(float(m.source@m.w)*.6*1e-7/3e-6,1.,places=13)
        self.assertEqual(float(m.field(0.)['outgoing_right_flux']@m.w),0.)
        self.assertAlmostEqual(float(m.field(.003)['outgoing_right_flux']@m.w)*1e-7/3e-6,1.,places=12)
        for direction in ['right','both']:
            c=config();c['spatial_sources']=[dict(light,id='light',mode='radiation',position_cm=.005,width_cm=0.,mu=.4,direction=direction)]
            m=build(c);out=m.field(0.);power=float((out['outgoing_right_flux']+out['outgoing_left_flux'])@m.w)*1e-7;self.assertAlmostEqual(power/3e-6,1.,places=12)


class CompactResults(unittest.TestCase):
    def result(self,folder):
        c=config();c['processes']['electron_ionization']=True;c['probe']=dict(mode='flat',minimum_nm=400.,maximum_nm=800.,points=31)
        c['spatial_sources']=[dict(id='electron pulse',mode='electrons',density_cm3=10.,velocity_cm_s=1e7,position_cm=.005,start_s=.002,stop_s=.006,distribution=dict(Te_eV=.1))]
        m=build(c);m.solve();quality=m.export(folder);return m,ResultData(folder),quality

    def test_default_export_has_only_four_tables_one_store_and_a_readme(self):
        with tempfile.TemporaryDirectory() as folder:
            m,data,q=self.result(folder);files={p.name for p in Path(folder).iterdir()}
            self.assertEqual(files,{'сохранение.csv','кинетика_уровней.csv','реакции_по_времени.csv','топ5_по_времени.csv',FILE_STORE,'О_результатах.txt'})
            self.assertEqual(data['heavy_cm3'].shape[2],len(m.states));self.assertEqual(data['reaction_rates'].shape[2],len(m.channels));self.assertEqual(q['mass_status'],'passed');self.assertEqual(q['kinetic_status'],'passed')

    def test_all_species_reagents_and_products_remain_explicit(self):
        with tempfile.TemporaryDirectory() as folder:
            m,data,q=self.result(folder)
            with (Path(folder)/'кинетика_уровней.csv').open(encoding='utf-8-sig') as stream:rows=list(csv.DictReader(stream,delimiter=';'))
            ids={r['ID уровня'] for r in rows};self.assertEqual(ids,{s['id'] for s in m.states}|{'e'})
            reaction=next(r for r in data.reactions if r['kind']=='ionization');self.assertEqual(reaction['reactants']['e'],1);self.assertEqual(reaction['products']['e'],2)

    def test_reexport_removes_named_old_diagnostics_and_keeps_user_files(self):
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder);(path/'quality.json').write_text('{}');(path/'collision_partners.csv').write_text('old');(path/'my_notes.txt').write_text('keep');self.result(folder)
            self.assertFalse((path/'quality.json').exists());self.assertFalse((path/'collision_partners.csv').exists());self.assertEqual((path/'my_notes.txt').read_text(),'keep')

    def test_electron_kinetics_closes_locally_with_transport_and_injection(self):
        with tempfile.TemporaryDirectory() as folder:
            c=config();c.update(cells=2,thickness_cm=1.);c['processes']['electron_ionization']=True;c['spatial_sources']=[dict(id='e',mode='electrons',density_cm3=10.,velocity_cm_s=1e7,position_cm=.2,start_s=.002,stop_s=.006,distribution=dict(Te_eV=.1))]
            m=build(c);m.solve();quality=m.export(folder);data=ResultData(folder)
            self.assertLess(abs(data['kinetic_residual_e_cm3']).max()/m.initial_F.sum(),1e-11);self.assertEqual(quality['kinetic_status'],'passed');self.assertGreater(abs(m.snapshots[-1]['electron_transport']).max(),0.)
            self.assertEqual(len(build_plots(data,PlotOptions(metrics=['kinetic_e','top5'],target='electron'))),2)

    def test_top_five_are_ranked_from_all_rates_and_have_full_denominator(self):
        with tempfile.TemporaryDirectory() as folder:
            m,data,q=self.result(folder)
            with (Path(folder)/'топ5_по_времени.csv').open(encoding='utf-8-sig') as stream:rows=list(csv.DictReader(stream,delimiter=';'))
            initial=[r for r in rows if float(r['Время, с'])==0. and r['Цель']=='Все реакции' and r['Ячейка']=='1'];expected=top_indices(data['reaction_rates'][0,0]);self.assertEqual([r['ID реакции'] for r in initial],[data.reactions[k]['id'] for k in expected])
            for row,k in zip(initial,expected):self.assertAlmostEqual(float(row['Доля от всех процессов цели']),data['reaction_rates'][0,0,k]/data['reaction_rates'][0,0].sum())

    def test_selected_plots_have_large_text_start_markers_and_log_scale(self):
        with tempfile.TemporaryDirectory() as folder:
            m,data,q=self.result(folder);o=PlotOptions(metrics=['components','mass','top5'],y_scale='log',fontsize=18,curves_per_figure=3);pages=build_plots(data,o)
            self.assertTrue(pages)
            for page in pages:
                ax=page.figure.axes[0];self.assertEqual(ax.get_yscale(),'log');self.assertGreaterEqual(ax.xaxis.label.get_fontsize(),18)
                starts=[line for line in ax.lines if line.get_gid()=='pump_start'];self.assertTrue(any(abs(float(line.get_xdata()[0])-.002)<1e-15 for line in starts));self.assertLessEqual(len(page.curves),3)
            before=tuple(pages[0].figure.get_size_inches());pages[0].figure.set_size_inches(5,4,forward=False)
            saved=save_pages(pages,Path(folder)/'chosen',True);self.assertEqual(len(saved),len(pages));self.assertTrue(saved[0].with_suffix('.csv').is_file())
            again=save_pages(pages,Path(folder)/'chosen');self.assertNotEqual(saved[0],again[0]);self.assertTrue(saved[0].is_file())

    def test_log_time_reports_a_zero_start_without_putting_zero_on_log_axis(self):
        with tempfile.TemporaryDirectory() as folder:
            _,data,_=self.result(folder);data.meta['sources'][0]['start_s']=0.;page=build_plots(data,PlotOptions(metrics=['components'],x_scale='log'))[0];ax=page.figure.axes[0]
            self.assertTrue(any('0 с' in t.get_text() for t in ax.texts));self.assertGreater(ax.get_xlim()[0],0.)

    def test_top_history_has_gaps_instead_of_drawing_nonleaders(self):
        from types import SimpleNamespace
        rates=np.array([[[9,8,7,6,5,1]],[[1,8,7,6,5,9]],[[9,8,7,6,5,1]]],dtype=float)
        class Data:
            reactions=[dict(id=str(i)) for i in range(6)]
            def __getitem__(self,key):return rates
        curves=dict(top_history(Data(),PlotOptions()));self.assertTrue(np.isnan(curves[0][1]));self.assertTrue(np.isnan(curves[5][0]));self.assertEqual(curves[5][1],9.)

    def test_spatial_and_spectral_selection_use_the_chosen_snapshot(self):
        with tempfile.TemporaryDirectory() as folder:
            _,data,_=self.result(folder);o=PlotOptions(metrics=['levels','eedf','heavy_energy','velocity','absorption','transmission','emission','probe_distance','heavy_distribution'],levels=['H0'],axis='space',time_index=1,wavelength_min=450.,wavelength_max=700.);pages=build_plots(data,o)
            self.assertEqual(len(pages),9);spectral=next(p for p in pages if p.key.startswith('absorption'));self.assertGreaterEqual(spectral.x.min(),450.);self.assertLessEqual(spectral.x.max(),700.)

    def test_clearing_selection_does_not_silently_add_other_components(self):
        with tempfile.TemporaryDirectory() as folder:
            _,data,_=self.result(folder)
            with self.assertRaisesRegex(ValueError,'компонент'):build_plots(data,PlotOptions(metrics=['components'],components=[]))
