"""Offline analytic and adversarial checks for names, sources and progress."""
import csv
from copy import deepcopy
import hashlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import zipfile
import numpy as np
from atomic_core import Species, selected_species, CrossSection, Network
from atomic_database import SourceCache, parse_adf04_ls3, ccc_zip_tables, prepare_species, lithium_tables, BOHR_CM2
from progress import Progress, parse_event, PREFIX


SYNTHETIC_ADF = b''' H+ 0 1 1 109679.(1S)
 1 1S1 (2)0( 0.5) 0.
 2 2P1 (2)1( 2.5) 82303.
 -1
 1.00 3 1.00+03 2.00+03
 2 1 1.00+08 2.00-01 4.00-01 9.99+00
 -1
'''


class AtomicDataChecks(unittest.TestCase):
    def test_names_and_roman_stages_are_charge_consistent(self):
        for name, key in [('H','H_0'),('Н','H_0'),('водород','H_0'),('Li+','Li_1'),('Li II','Li_1'),('Li+2','Li_2'),('Li2+','Li_2'),('Li++','Li_2'),('He II','He_1')]:
            self.assertEqual(Species.parse(name).key,key)
        for name in ['He IV','Li IIII','H-','Li 1','H/../../Li']:
            with self.assertRaises(ValueError):Species.parse(name)

    def test_duplicate_species_aliases_are_rejected(self):
        with self.assertRaises(ValueError):selected_species({'atomic_data':{'species':['Li+','Li II']}})

    def test_adf_temperature_tail_is_not_a_ninth_temperature(self):
        levels,rates,lines,meta=parse_adf04_ls3(SYNTHETIC_ADF,Species('H'),1.008,'SYNTHETIC',3)
        self.assertEqual(len(levels),2);self.assertEqual(len(rates),2)
        self.assertEqual([r['value'] for r in rates],[.2,.4])
        self.assertEqual(levels[1]['g'],6)
        self.assertIn('2P',levels[1]['label']);self.assertEqual(lines[0]['A_s'],1e8)
        self.assertEqual(meta['parents'][0]['threshold_cm1'],109679.)

    def test_wrong_species_and_ic_data_are_rejected(self):
        with self.assertRaises(ValueError):parse_adf04_ls3(SYNTHETIC_ADF,Species('Li'),6.94,'SYNTHETIC',3)
        with self.assertRaises(ValueError):parse_adf04_ls3(SYNTHETIC_ADF.replace(b'( 2.5)',b'( 1.5)'),Species('H'),1.008,'SYNTHETIC',3)

    def test_uppercase_adas_element_header(self):
        raw=b'HE+ 0 2 1 198310.7(2S)\n 1 1S2 (1)0(0.0) 0.\n -1\n 1.0 3 1000 2000\n -1\n'
        levels,*_=parse_adf04_ls3(raw,Species('He'),4.002602,'SYNTHETIC',3)
        self.assertEqual(levels[0]['species'],'He I');self.assertEqual(levels[0]['g'],1)

    def test_hydrogenic_lowercase_configuration_annotation_is_not_principal_n(self):
        raw=b'HE+ 1 2 2 438908.(1S)\n 1 1s(2S) (2)0(0.5) 0.\n 2 5p(2P) (2)1(2.5) 420000.\n -1\n 1.0 3 1000 2000\n -1\n'
        levels,*_=parse_adf04_ls3(raw,Species('He',1),4.002602,'SYNTHETIC',9)
        self.assertEqual([r['n'] for r in levels],[1,5])

    def test_cross_section_rate_is_direct_sigma_v_integral(self):
        table=CrossSection([1.,10.],[2e-16,2e-16],1.,'SYNTHETIC')
        e=np.array([.5,2.,5.]);n=np.array([3e6,4e6,8e6])
        expected=np.sum(n[1:]*2e-16*np.sqrt(2*e[1:]*1.602176634e-12/9.1093837139e-28))
        self.assertAlmostEqual(table.rate(e,n)/expected,1.,places=14)

    def test_missing_energy_interval_is_not_filled_with_zero(self):
        table=CrossSection([2.,10.],[1e-16,1e-16],1.,'SYNTHETIC')
        for e in [1.5,11.]:
            with self.assertRaises(ValueError):table.at([e])
        self.assertEqual(table.at([.5])[0],0.)
        self.assertEqual(table.rate([.5,100.],[1.,0.]),0.)

    def test_cross_element_penning_network_conserves_each_element_and_charge(self):
        states=[dict(id='He*',label='He I 2 3S',nuclei={'He':1},charge=0),
                dict(id='He',label='He I 1 1S',nuclei={'He':1},charge=0),
                dict(id='Li',label='Li I 2 2S',nuclei={'Li':1},charge=0),
                dict(id='Li+',label='Li II 1S2 1S',nuclei={'Li':1},charge=1),
                dict(id='e',label='e',nuclei={},charge=-1)]
        reactions=[dict(reactants={'He*':1,'Li':1},products={'He':1,'Li+':1,'e':1})]
        net=Network(states,reactions);dn=net.derivative([3.])
        np.testing.assert_array_equal(net.nuclei@dn,0.);self.assertEqual(net.charge@dn,0.)
        self.assertIn('3S',net.labels[0]);self.assertIn('Li II',net.labels[0])
        reactions[0]['products'].pop('e')
        with self.assertRaises(ValueError):Network(states,reactions)

    def test_fractional_nuclear_composition_and_charge_are_rejected(self):
        with self.assertRaises(ValueError):Network([dict(id='x',label='x',nuclei={'H':.5},charge=0)],[])
        with self.assertRaises(ValueError):Network([dict(id='x',label='x',nuclei={'H':1},charge=.5)],[])

    def test_ccc_reads_header_direction_not_filename_and_converts_a0(self):
        levels,*_=parse_adf04_ls3(SYNTHETIC_ADF,Species('H'),1.008,'SYNTHETIC',3)
        stream=io.BytesIO()
        with zipfile.ZipFile(stream,'w') as z:
            z.writestr('1S.2P','# 2P -> 1S; Units: a0^2\n0.01 2.\n1. 3.\n')
        rows,catalogue=ccc_zip_tables(stream.getvalue(),levels,'hydrogenic','SYNTHETIC')
        self.assertEqual(catalogue[0]['kind'],'electron_deexcitation')
        self.assertEqual(catalogue[0]['initial'],levels[1]['id'])
        self.assertAlmostEqual(rows[0]['sigma_cm2']/BOHR_CM2,2.)

    def test_ccc_unknown_units_are_rejected(self):
        levels,*_=parse_adf04_ls3(SYNTHETIC_ADF,Species('H'),1.008,'SYNTHETIC',3)
        stream=io.BytesIO()
        with zipfile.ZipFile(stream,'w') as z:z.writestr('a','# 1S -> 2P; Units: mystery\n1 1\n2 2\n')
        with self.assertRaises(ValueError):ccc_zip_tables(stream.getvalue(),levels,'hydrogenic','SYNTHETIC')

    def test_ccc_cm2_is_not_converted_a_second_time(self):
        levels,*_=parse_adf04_ls3(SYNTHETIC_ADF,Species('H'),1.008,'SYNTHETIC',3)
        stream=io.BytesIO()
        with zipfile.ZipFile(stream,'w') as z:z.writestr('a','# 1S -> 2P; Units: cm^2\n1 1e-16\n2 2e-16\n')
        rows,_=ccc_zip_tables(stream.getvalue(),levels,'hydrogenic','SYNTHETIC')
        self.assertEqual(rows[0]['sigma_cm2'],1e-16)

    def test_cache_reuses_verified_source_and_detects_tampering(self):
        raw=b'SYNTHETIC bytes'
        info=dict(url='https://example.invalid/data',sha256=hashlib.sha256(raw).hexdigest())
        with tempfile.TemporaryDirectory() as d,patch('atomic_database.urlopen',return_value=io.BytesIO(raw)) as fetch:
            cache=SourceCache(d);self.assertEqual(cache.get('a.dat',info),raw)
            self.assertEqual(SourceCache(d,offline=True).get('a.dat',info),raw)
            self.assertEqual(fetch.call_count,1)
            path=next(p for p in Path(d).iterdir() if p.name.endswith('a.dat'));path.write_bytes(b'changed')
            with self.assertRaises(ValueError):cache.get('a.dat',info)

    def test_changed_pinned_source_is_not_written(self):
        with tempfile.TemporaryDirectory() as d,patch('atomic_database.urlopen',return_value=io.BytesIO(b'new')):
            cache=SourceCache(d)
            with self.assertRaises(ValueError):cache.get('a.dat',dict(url='https://example.invalid/data',sha256='0'*64))
            self.assertFalse(list(Path(d).iterdir()))

    def make_registry(self,root):
        info=dict(url='https://example.invalid/synthetic.dat',sha256=hashlib.sha256(SYNTHETIC_ADF).hexdigest(),hash_mode='bytes',format='adf04_ls3')
        (root/'atomic_sources.json').write_text(json.dumps(dict(schema='AtomicSources-1',sources={'test':info},species={'H_0':dict(mass_u=1.008,adf04='test')})))
        return dict(atomic_data=dict(species=['H'],discover=False))

    def test_pack_reuse_offline_and_manual_edits_preserved(self):
        with tempfile.TemporaryDirectory() as d,patch('atomic_database.urlopen',return_value=io.BytesIO(SYNTHETIC_ADF)) as fetch:
            root=Path(d);config=self.make_registry(root)
            rows=prepare_species(config,root);self.assertEqual(rows[0]['levels'],2)
            config['atomic_data']['offline']=True
            prepare_species(config,root);self.assertEqual(fetch.call_count,1)
            pack=root/'atomic_data/H_0_n3'
            report=json.loads((pack/'manifest.json').read_text(encoding='utf-8-sig'));self.assertFalse(report['kinetic_ready'])
            table=pack/'levels.csv';edited=table.read_bytes()+b'\n';table.write_bytes(edited)
            with self.assertRaisesRegex(ValueError,'ручные правки'):prepare_species(config,root)
            self.assertEqual(table.read_bytes(),edited)

    def test_unavailable_source_produces_explicit_missing_report(self):
        with tempfile.TemporaryDirectory() as d,patch('atomic_database.urlopen',side_effect=OSError('offline test')):
            root=Path(d);config=self.make_registry(root);rows=prepare_species(config,root)
            self.assertEqual(rows[0]['levels'],0)
            m=json.loads((root/'atomic_data/H_0_n3/manifest.json').read_text(encoding='utf-8-sig'))
            self.assertTrue(m['download_error']);self.assertIn('offline test',m['gaps'][0])

    def test_one_missing_lithium_channel_does_not_discard_good_tables(self):
        from types import SimpleNamespace
        form=b'<input name="low" value="2s_dS"><input name="low" value="3s_dS"><input name="upp" value="2s_dS"><input name="upp" value="3s_dS">'
        good=b'Cross section in units of cm<SUP>2</SUP><TR><TD ALIGN=RIGHT>1.0<TD><TD><TD ALIGN=RIGHT>1e-16\n<TR><TD ALIGN=RIGHT>2.0<TD><TD><TD ALIGN=RIGHT>2e-16'
        receipts={}
        def fetch(name,info):
            if name.lower()=='ccc_li_form.html':return form
            if '3S_3S' in name:raise OSError('SYNTHETIC missing channel')
            receipts[name]={'sha256_bytes':'synthetic'};return good
        levels=[dict(id='Li_a',species='Li I',label='Li I 2 2S',n=2,L=0,energy_cm1=0.),dict(id='Li_b',species='Li I',label='Li I 3 2S',n=3,L=0,energy_cm1=100.)]
        rows,catalogue,gaps=lithium_tables(SimpleNamespace(get=fetch,receipts=receipts),levels)
        self.assertEqual(len(catalogue),2);self.assertEqual(len(rows),4)
        self.assertTrue(any('SYNTHETIC missing channel' in gap for gap in gaps))

    def test_negative_ccc_channel_is_rejected_without_losing_valid_channels(self):
        from types import SimpleNamespace
        form=b'<input name="low" value="4s_dS"><input name="low" value="4p_dP"><input name="upp" value="4s_dS"><input name="upp" value="4p_dP">'
        good=b'Cross section in units of cm<SUP>2</SUP><TR><TD ALIGN=RIGHT>1<TD><TD><TD ALIGN=RIGHT>1e-16\n<TR><TD ALIGN=RIGHT>2<TD><TD><TD ALIGN=RIGHT>2e-16'
        receipts={};urls=[]
        def fetch(name,info):
            urls.append(info['url'])
            if name=='ccc_Na_form.html':return form
            receipts[name]={'sha256_bytes':'SYNTHETIC'}
            return good.replace(b'2e-16',b'-2e-16') if '4S_4P' in name else good
        levels=[dict(id='a',label='Na 4s',n=4,L=0,energy_cm1=0.),dict(id='b',label='Na 4p',n=4,L=1,energy_cm1=100.)]
        from atomic_database import alkali_tables
        rows,catalogue,gaps=alkali_tables(SimpleNamespace(get=fetch,receipts=receipts),levels,'Na')
        self.assertEqual(len(catalogue),2);self.assertTrue(all(float(r['sigma_cm2'])>0 for r in rows))
        self.assertTrue(any('nonphysical' in gap for gap in gaps))
        self.assertTrue(all('Electron' in u or 'ELECTRON' in u for u in urls))


class ProgressChecks(unittest.TestCase):
    def test_eta_waits_for_measured_work_and_tracks_remaining_fraction(self):
        now=[0.];lines=[];p=Progress('json',lines.append,interval_s=.05,now=lambda:now[0])
        p.start('Интегрирование','test');p.update(0,10)
        for j in range(1,4):
            now[0]=float(j);p.update(j,10,force=True)
            event=parse_event(lines[-1])
            if j<3:self.assertIsNone(event['eta_s'])
        self.assertAlmostEqual(event['fraction'],.3);self.assertAlmostEqual(event['eta_s'],7.)
        p.done();self.assertEqual(parse_event(lines[-1])['status'],'done')

    def test_unmeasurable_stage_has_no_fake_eta(self):
        lines=[];p=Progress('json',lines.append);p.start('JIT','compiler')
        self.assertIsNone(parse_event(lines[-1])['fraction']);self.assertIsNone(parse_event(lines[-1])['eta_s'])

    def test_invalid_progress_lines_and_reversal(self):
        self.assertIsNone(parse_event('ordinary text'))
        self.assertIsNone(parse_event(PREFIX+'{bad'))
        self.assertIsNone(parse_event(PREFIX+json.dumps(dict(protocol='kinetic-progress-1',stage='x',fraction=2))))
        p=Progress();p.start('x');p.update(.5,1)
        with self.assertRaises(ValueError):p.update(.4,1)

    def test_launcher_routes_data_preparation_and_old_schemas(self):
        from launcher import configuration,settings
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'data.json';path.write_text(json.dumps(dict(schema='AtomicData-1.0',atomic_data=dict(prepared_directory='tables'))))
            self.assertEqual(configuration(path)[1],'prepare_species.py')
            self.assertEqual(settings(path)[1],Path(d)/'tables')


if __name__=='__main__':unittest.main()
