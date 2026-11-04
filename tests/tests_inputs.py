"""Input regressions: real user layouts, units and gas inventory changes."""
import json
from pathlib import Path
import tempfile
import unittest
import zipfile
import numpy as np
from experiment_config import gas_initial,prepare_configuration,KB_SI
from spectrum_input import load_spectrum,HC_EV_A,EV_ERG
from mixture_kinetics import ROOT
from launcher import PRESETS


class SpectrumImport(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name)

    def read(self,text,suffix='.txt',options=None,purpose='probe',encoding='utf-8'):
        p=self.root/('spectrum'+suffix);p.write_text(text,encoding=encoding);return load_spectrum(p,options,purpose)

    def test_plain_whitespace_no_header_defaults_nm(self):
        s=self.read('1082 1\n1084 2\n');self.assertEqual(s.wavelength_A,(10820.,10840.));self.assertTrue(s.metadata['warnings'])

    def test_comma_csv_spaces_and_arbitrary_header(self):
        s=self.read('wavelength_nm, detector\n1082.0, 1.0\n1084.0, 2.0\n','.csv');self.assertEqual(s.values,(1.,2.));self.assertEqual(s.metadata['input_x_unit'],'nm')

    def test_russian_decimal_comma_semicolon_and_bom(self):
        s=self.read('Длина волны, нм;Интенсивность\n1082,5;1,25\n1083,5;2,5\n','.csv',encoding='utf-8-sig');self.assertEqual(s.wavelength_A,(10825.,10835.));self.assertEqual(s.values,(1.25,2.5))

    def test_decimal_comma_whitespace_and_fortran_exponents(self):
        s=self.read('# nm counts\n1082,5 1D0\n1083,5 2D0\n');self.assertEqual(s.values,(1.,2.))

    def test_descending_duplicate_coordinates_are_sorted_averaged(self):
        s=self.read('1084 2\n1082 1\n1084 4\n');self.assertEqual(s.values,(1.,3.));self.assertEqual(len(s.metadata['warnings']),3)

    def test_canonical_legacy_format_keeps_cgs_per_angstrom(self):
        s=self.read('wavelength_A;four_pi_J_lambda_erg_cm2_s_A\n10820;3\n10840;4\n','.csv',purpose='driving');self.assertEqual(s.wavelength_A,(10820.,10840.));self.assertEqual(s.values,(3.,4.))

    def test_selected_columns_and_cp1251(self):
        s=self.read('номер;длина, нм;сигнал\n1;1082;3\n2;1084;4\n','.csv',options={'x_column':1,'y_column':2},encoding='cp1251');self.assertEqual(s.values,(3.,4.))

    def test_utf16_tab_table(self):
        s=self.read('wavelength [nm]\tI\n1082\t1\n1084\t2\n','.tsv',encoding='utf-16');self.assertEqual(s.wavelength_A,(10820.,10840.))

    def test_physical_nm_jacobian_and_watt_conversion(self):
        s=self.read('nm power\n100 2\n200 4\n',options={'y_kind':'energy','energy_unit':'W_m2'},purpose='driving');np.testing.assert_allclose(s.values,[200.,400.])

    def test_energy_and_frequency_jacobians(self):
        for unit,x,expected in [('eV',[1.,2.],lambda lam:HC_EV_A/lam**2),('Hz',[1e14,2e14],lambda lam:2.99792458e18/lam**2)]:
            s=self.read(f'{x[0]} 1\n{x[1]} 1\n',options={'x_unit':unit,'y_kind':'energy'},purpose='driving')
            np.testing.assert_allclose(s.values,[expected(lam) for lam in s.wavelength_A],rtol=1e-14)

    def test_photon_conversion_and_normal_flux_angle(self):
        s=self.read('100 1\n200 1\n',options={'x_unit':'nm','y_kind':'photons','angular_kind':'normal_flux','mu':.5},purpose='driving')
        np.testing.assert_allclose(s.values,[HC_EV_A/1000*EV_ERG/10/.5,HC_EV_A/2000*EV_ERG/10/.5])

    def test_relative_counts_cannot_silently_drive_kinetics(self):
        with self.assertRaisesRegex(ValueError,'абсолютная'):self.read('1 1\n2 2',options={'y_kind':'relative'},purpose='driving')

    def test_negative_absolute_rejected_relative_counts_preserved(self):
        s=self.read('1 -1\n2 2');self.assertEqual(s.values[0],-1.)
        with self.assertRaisesRegex(ValueError,'отрицательная'):self.read('1 -1\n2 2',purpose='driving')

    def test_explicit_negative_clipping(self):
        s=self.read('1 -1\n2 2',options={'negative_policy':'clip'});self.assertEqual(s.values[0],0.)

    def test_invalid_values_have_explanatory_errors(self):
        for text in ['1 nan\n2 2','0 1\n2 2','1 bad\n2 2','1 1','1 1\n1 2']:
            with self.subTest(text=text),self.assertRaises(ValueError):self.read(text)

    def test_xlsx_first_sheet_inline_strings_and_numbers(self):
        p=self.root/'probe.xlsx'
        with zipfile.ZipFile(p,'w') as z:
            z.writestr('xl/workbook.xml','<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets><sheet name="probe" sheetId="1" r:id="rId1"/></sheets></workbook>')
            z.writestr('xl/_rels/workbook.xml.rels','<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Target="worksheets/sheet1.xml"/></Relationships>')
            z.writestr('xl/worksheets/sheet1.xml','<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><sheetData><row r="1"><c r="A1" t="inlineStr"><is><t>wavelength_nm</t></is></c><c r="B1" t="inlineStr"><is><t>counts</t></is></c></row><row r="2"><c r="A2"><v>1082</v></c><c r="B2"><v>1</v></c></row><row r="3"><c r="A3"><v>1084</v></c><c r="B3"><v>2</v></c></row></sheetData></worksheet>')
        s=load_spectrum(p);self.assertEqual(s.wavelength_A,(10820.,10840.));self.assertEqual(s.values,(1.,2.))


class GasInput(unittest.TestCase):
    def gas(self,**kw):return dict(mode='densities',temperature_K=300.,components=[dict(species='He',density_cm3=9e9),dict(species='H+',density_cm3=1e9)],**kw)

    def test_density_inventory_includes_ions(self):
        n,stages,T=gas_initial(self.gas(),ROOT/'configs');self.assertEqual(n,{'He':9e9,'H+':1e9});self.assertEqual([s.label for s in stages],['H II','He I']);self.assertEqual(T,300.)

    def test_fractions_and_pressure_fill_use_particle_numbers(self):
        parts=[dict(species='He',fraction=.9),dict(species='H',fraction=.1)]
        for mode,kwargs,total in [('fractions',{'total_density_cm3':1e10},1e10),('pressure',{'pressure_Pa':1.},1/(KB_SI*300)*1e-6)]:
            n,_,_=gas_initial(dict(mode=mode,temperature_K=300,components=parts,**kwargs),ROOT/'configs');self.assertAlmostEqual(n['He']/total,.9);self.assertAlmostEqual(sum(n.values())/total,1.)

    def test_metastable_alias_and_molecule_inventory(self):
        n,stages,T=gas_initial(dict(components=[dict(species='HeI_2_3S',density_cm3=1e6),dict(species='H2',density_cm3=1e8)]),ROOT/'configs');self.assertEqual([s.element for s in stages],['H','He'])
        cfg=prepare_configuration(dict(gas=dict(temperature_K=650.,components=[dict(species='H2',density_cm3=1e8)])),ROOT/'configs');self.assertEqual(cfg['heavy_initial']['H2']['Ti_K'],650.)

    def test_invalid_compositions_fail_before_calculation(self):
        cases=[{'components':[]},{'components':[dict(species='He',density_cm3=-1)]},{'components':[dict(species='He',density_cm3=0)]},{'components':[dict(species='bad!',density_cm3=1)]},{'components':[dict(species='He',density_cm3=1),dict(species='He',density_cm3=2)]},dict(mode='fractions',total_density_cm3=1e9,components=[dict(species='He',fraction=.9)])]
        for gas in cases:
            with self.subTest(gas=gas),self.assertRaises(ValueError):gas_initial(gas,ROOT/'configs')

    def test_editing_mixture_replaces_stale_states_and_drift(self):
        old={'gas':self.gas(),'initial_cm3':{'Li':1e10},'atomic_data':{'species':['Li']},'heavy_initial':{'Li':{'Ti_K':700,'uz_cm_s':1e6}},'level_selection':{'forced_levels':['Li_0:L002','He_0:L002']}}
        new=prepare_configuration(old,ROOT/'configs');self.assertNotIn('Li',new['initial_cm3']);self.assertEqual(new['heavy_initial'],{'H':{'Ti_K':300.},'He':{'Ti_K':300.}});self.assertEqual(new['level_selection']['forced_levels'],['He_0:L002']);self.assertIn('Li',old['initial_cm3'])

    def test_default_launch_has_no_forced_sources_drift_or_estimates(self):
        self.assertEqual(PRESETS[0],'configs/gas_probe.json');c=json.loads((ROOT/PRESETS[0]).read_text());c=prepare_configuration(c,ROOT/'configs')
        self.assertEqual(c['spatial_sources'],[]);self.assertEqual(c['electron_initial'],{'Te_eV':.8,'kind':'maxwellian','match_moments':True});self.assertEqual(c['radiation']['scale'],0.);self.assertFalse(c['chemistry']['allow_estimated_rates']);self.assertTrue(all('uz_cm_s' not in v for v in c['heavy_initial'].values()))
