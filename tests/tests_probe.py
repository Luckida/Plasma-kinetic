"""Analytic weak-probe limits and pump integration; no physical fit added."""
import csv
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
from scipy.special import voigt_profile
from mixture_kinetics import MixtureSlab,ROOT
from probe_spectrum import H,C,EV,KB,AMU
from tests_mixture import config,build,synthetic_bank


def with_line(cfg):
    bank,photo=synthetic_bank(upper=True)
    bank.lines=[dict(id='SYNTHETIC:test',lower='H0',upper='Hup',A_s=1e7,pump=True,source='SYNTHETIC analytic test')]
    with patch('mixture_kinetics.MixtureData',return_value=bank),patch('mixture_kinetics.PhotoBank',return_value=photo):return MixtureSlab(cfg,ROOT)


class WeakProbe(unittest.TestCase):
    def cfg(self):
        c=config();c['probe']=dict(mode='flat',minimum_nm=1220,maximum_nm=1260,points=101);return c

    def test_beer_lambert_line_center_with_statistical_weights(self):
        m=with_line(self.cfg());tau,tr,_=m.probe.evaluate(m.n);line=m.probe.lines[0];center=np.argmin(abs(m.probe.wavelength-line['center']))
        width=line['nu0']/C*np.sqrt(KB*m.Tn/(m.states[line['j']]['mass_u']*AMU));cross=1e7*(C/line['nu0'])**2/(8*np.pi)*voigt_profile(0.,width,1e7/(4*np.pi))
        expected=m.initial[0,line['i']]*3*cross*m.length
        self.assertAlmostEqual(tau[center]/expected,1.,places=9);self.assertAlmostEqual(tr[center],np.exp(-expected),places=12)

    def test_transparent_gas_and_population_inversion(self):
        m=with_line(self.cfg());n=np.zeros_like(m.n);tau,tr,_=m.probe.evaluate(n);np.testing.assert_array_equal(tr,np.ones_like(tr))
        n[:,m.index['Hup']]=1e3;tau,tr,_=m.probe.evaluate(n);self.assertGreater(tr.max(),1.);self.assertLess(tau.min(),0.)

    def test_line_wings_are_not_cut_when_center_is_outside_probe_band(self):
        c=self.cfg();c['probe'].update(minimum_nm=1220,maximum_nm=1230);m=with_line(c);tau,tr,_=m.probe.evaluate(m.n)
        self.assertTrue(np.all(tau>0));self.assertFalse(m.probe.lines[0]['in_band'])
        nu=C/(m.probe.wavelength[0]*1e-8);line=m.probe.lines[0];cross=line['A']*(C/line['nu0'])**2/(8*np.pi)*voigt_profile(nu-line['nu0'],line['sigma'],line['gamma'])
        self.assertAlmostEqual(tau[0]/(m.initial[0,line['i']]*line['ratio']*cross*m.length),1.,places=12)

    def test_position_and_oblique_path_are_integrated(self):
        c=self.cfg();c['cells']=2;c['probe']['mu']=.5;m=with_line(c);tau,tr,centers=m.probe.evaluate(m.n)
        np.testing.assert_allclose(centers[0],tau/4);np.testing.assert_allclose(centers[1],tau*3/4)
        c['probe']['mu']=1.;normal=with_line(c);np.testing.assert_allclose(tau,2*normal.probe.evaluate(normal.n)[0])

    def test_probe_changes_no_kinetics_even_during_radiative_decay(self):
        c=self.cfg();c['initial_cm3']['Hup']=1e5;c['processes']['spontaneous']=True;c['solver'].update(t_end_s=2e-8,max_step_s=1e-9,min_step_s=1e-16)
        a=with_line(c);c['probe']={'mode':'off'};b=with_line(c);a.solve();b.solve()
        np.testing.assert_array_equal(a.P,b.P);np.testing.assert_array_equal(a.F,b.F);self.assertLess(a.n[0,a.index['Hup']],1e5)

    def test_export_contains_wavelength_time_and_position(self):
        m=with_line(self.cfg());m.solve()
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d);m.probe.export(folder)
            with (folder/'probe_transmission.csv').open(encoding='utf-8-sig') as f:rows=list(csv.DictReader(f,delimiter=';'))
            self.assertGreater(len(rows),101);self.assertTrue(all(0<=float(r['transmission'])<=1 for r in rows));self.assertTrue((folder/'probe_by_position_final.csv').is_file());self.assertTrue((folder/'probe_lines.csv').is_file())

    def test_disabled_bad_pump_file_and_history_do_not_block(self):
        c=config();c['phase_grid']['energy_max_eV']='auto';c['radiation'].update(enabled=False,scale=1.,spectrum_file='missing.csv',history_file='missing_history.csv');m=build(c);self.assertEqual(m.scale,0.);np.testing.assert_array_equal(m.source,np.zeros_like(m.source))

    def test_zero_scale_internal_light_does_not_read_bad_file(self):
        c=config();c['phase_grid']['energy_max_eV']='auto';c['spatial_sources']=[dict(id='off',mode='radiation',scale=0.,spectrum_file='missing.txt')];m=build(c);self.assertEqual(m.sources.items,[])

    def test_noncanonical_driving_spectrum_and_exact_pump_switches(self):
        with tempfile.TemporaryDirectory() as d:
            file=Path(d)/'light.csv';file.write_text('wavelength_nm, detector\n600, 1e-6\n500, 2e-6\n500, 2e-6\n')
            c=config();c['radiation'].update(spectrum_file=str(file),spectrum_options={'y_kind':'energy'},start_s=.002,stop_s=.006);m=build(c)
            self.assertGreater(m.source.max(),0.);np.testing.assert_array_equal(m.field(0.)['mean_4piJ'],0.);self.assertGreater(m.field(.003)['mean_4piJ'].max(),0.);np.testing.assert_array_equal(m.field(.006)['mean_4piJ'],0.)
            m.solve();times=[s['time_s'] for s in m.snapshots];self.assertIn(.002,times);self.assertIn(.006,times)

    def test_file_probe_accepts_plain_relative_counts(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'probe.txt';p.write_text('1220 2\n1260 4\n');c=config();c['probe']=dict(mode='file',spectrum_file=str(p));m=with_line(c)
            self.assertEqual(m.probe.incident[0],2.);self.assertEqual(m.probe.incident[-1],4.);self.assertEqual(m.probe.metadata['input_y_kind'],'relative')
