"""Conservation, analytic kinetics and import edge cases, entirely offline."""
from copy import deepcopy
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import numpy as np
from atomic_core import Network
from atomic_database import HC_EV_CM,write_csv,write_json,load_catalogue
from atomic_database import native_tables,LEVEL_FIELDS,LINE_FIELDS,XS_FIELDS
from mixture_chemistry import append_components,rate_fit
from mixture_kinetics import MixtureSlab,ROOT
from data_import import parse_umist,save_network,convert_lxcat


def bank():
    states=[]
    for element,mass,ip in [('H',1.008,13.598434),('He',4.002602,24.587389)]:
        for q in [0,1]:
            states.append(dict(id=f'{element}_{q}:ground',label=f'{element} {q} ground',element=element,species=f'{element} {q}',charge=q,
                mass_u=mass,n=1,multiplicity=2 if element=='H' else 1,L=0,g=2 if element=='H' else 1,nuclei={element:1},configuration='test ground',
                energy_cm1=q*ip/HC_EV_CM,excitation_eV=0.,is_ground=True))
    states.append(dict(states[2],id='He_metastable',label='He 2^3S',n=2,multiplicity=3,g=3,energy_cm1=19.819615/HC_EV_CM,excitation_eV=19.819615,is_ground=False))
    obj=SimpleNamespace(states=states,elements=['H','He'],stage_ground={f'{e}_{q}':f'{e}_{q}:ground' for e in ['H','He'] for q in [0,1]},
                        lines=[],cross={},cross_source={},elastic={},parent_index={},manifests={})
    append_components(obj,dict(schema='MixtureKineticSlab-2.1'),ROOT)
    obj.index={s['id']:i for i,s in enumerate(obj.states)};obj.by={s['id']:s for s in obj.states}
    aliases={'H':'H_0:ground','H+':'H_1:ground','He':'He_0:ground','He+':'He_1:ground'}
    obj.resolve=lambda n:obj.index[aliases.get(n,n)]
    return obj


def configuration(ids):
    all_ids=[r['id'] for f in (ROOT/'data/reactions').glob('*/reactions.json') for r in json.loads(f.read_text(encoding='utf-8-sig'))['reactions']]
    return dict(schema='MixtureKineticSlab-2.1',initial_cm3={'H':1e9,'He':1e8,'He_metastable':1e5},cells=1,thickness_cm=.01,
        level_selection=dict(mode='all',relative_tolerance=1e-6,missing_data='report'),electron_initial=dict(Te_eV=.01,hot_Te_eV=.01),
        heavy_initial={s['element']:dict(Ti_K=650.) for s in bank().states},
        phase_grid=dict(energy_min_eV=1e-8,energy_max_eV=50.,energy_points=32,angular_points=4),
        heavy_phase_grid=dict(energy_min_eV=1e-8,energy_max_eV=50.,energy_points=16,angular_points=4),
        processes={k:False for k in ['electron_bound_bound','electron_ionization','spontaneous','two_photon','photo_bound_bound','photo_bound_free','radiative_recombination','stimulated_recombination']},
        collisions=dict(ee_bgk=False,ii_bgk=False,ei_binary=False,neutral_elastic=False),radiation=dict(continuum_points=30,electron_scattering=False),
        chemistry=dict(enabled=True,rate_policy='maxwellian_moments',disabled_ids=[n for n in all_ids if n not in ids],product_angular_points=2,product_azimuth_points=2),
        performance=dict(workers=1,radiation_backend='numpy'),solver=dict(t_end_s=1e-3,max_step_s=1e-5,min_step_s=1e-12,samples=3))


def build(cfg):
    photo=SimpleNamespace(tables={},sources={},gaps=[],receipts={})
    with patch('mixture_kinetics.MixtureData',return_value=bank()),patch('mixture_kinetics.PhotoBank',return_value=photo):return MixtureSlab(cfg,ROOT)


class ChemicalPhysics(unittest.TestCase):
    def test_penning_and_association_share_total_and_different_products(self):
        ids=['RD82:He3S_H:PI','RD82:He3S_H:AI'];model=build(configuration(ids))
        rates=model.chemistry.event_rates(0)
        self.assertAlmostEqual(sum(rates.values())/(5e-10*1e9*1e5),1.,places=13)
        model.solve();final=model.snapshots[-1];lost=1e5-final['n'][0,model.index['He_metastable']]
        emitted=final['F'].sum();hydrogen_ion=final['n'][0,model.data.resolve('H+')];molecular=final['n'][0,model.data.resolve('HeH+')]
        self.assertGreater(molecular,0.);self.assertGreater(hydrogen_ion,0.)
        self.assertAlmostEqual(emitted/lost,1.,places=11)
        self.assertAlmostEqual((hydrogen_ion+molecular)/lost,1.,places=11)
        a,b,k,t=1e5,1e9,5e-10,model.end;difference=b-a
        exact=a*difference*np.exp(-k*difference*t)/(b-a*np.exp(-k*difference*t))
        self.assertLess(abs(final['n'][0,model.index['He_metastable']]/exact-1),.002)
        self.assertLess(abs(final['charge'])/model.initial_nuclei.sum(),1e-12)
        self.assertLess(abs(final['nuclei']).max()/model.initial_nuclei.sum(),1e-12)
        self.assertLess(abs(final['energy'])/abs(model.initial_energy),2e-7)
        with tempfile.TemporaryDirectory() as folder:
            quality=model.export(folder)
            self.assertEqual(quality['kinetic_status'],'passed');self.assertTrue((Path(folder)/'реакции_по_времени.csv').exists())

    def test_step_refinement_converges_to_bimolecular_solution(self):
        ids=['RD82:He3S_H:PI','RD82:He3S_H:AI'];errors=[]
        for step in [2e-5,1e-5]:
            cfg=configuration(ids);cfg['solver'].update(t_end_s=2e-4,max_step_s=step)
            model=build(cfg);model.solve();a,b=1e5,1e9;d=b-a;e=np.exp(-5e-10*d*model.end);exact=a*d*e/(b-a*e)
            errors.append(abs(model.n[0,model.index['He_metastable']]-exact))
        self.assertLess(errors[1],.6*errors[0])

    def test_resonant_exchange_changes_velocities_without_changing_populations(self):
        cfg=configuration(['LS05:Hplus_H:CX']);cfg['initial_cm3']={'H':1e7,'H+':1e7};cfg['heavy_initial']['H 1']=dict(Ti_K=300.,uz_cm_s=3e6)
        cfg['heavy_phase_grid']['energy_max_eV']=100.;model=build(cfg)
        i=model.data.resolve('H');j=model.data.resolve('H+');before=model.n.copy();energy=model.total_energy()
        initial_neutral=model.gh['H'].moments(model.P[0,i])['uz_cm_s'];events=model.chemistry.evolve(0,1e-6)
        self.assertIsNotNone(events);np.testing.assert_allclose(model.n,before,rtol=1e-14)
        self.assertGreater(model.gh['H'].moments(model.P[0,i])['uz_cm_s'],initial_neutral)
        self.assertAlmostEqual(model.total_energy()/energy,1.,places=13)

    def test_same_reactant_stoichiometry_consumes_two_particles(self):
        model=build(configuration([]))
        states=[dict(s,charge=int(s['charge'])) for s in model.states]+[dict(id='e',label='electron',nuclei={},charge=-1)]
        network=Network(states,[dict(reactants={'He_metastable':2},products={'He_0:ground':1,'He_1:ground':1,'e':1})])
        self.assertEqual(network.matrix[model.index['He_metastable'],0],-2)
        np.testing.assert_array_equal(network.nuclei@network.matrix,0.)

    def test_actual_self_penning_has_two_losses_and_one_electron_per_event(self):
        cfg=configuration(['HE25:15:PI','HE25:14:AI']);cfg['initial_cm3']={'He':1e9,'He_metastable':1e8}
        model=build(cfg);before=model.n.copy();energy=model.total_energy();nuclei=model.totals().copy()
        counts=model.chemistry.evolve(0,1e-4);total=counts.sum();lost=before[0,model.index['He_metastable']]-model.n[0,model.index['He_metastable']]
        self.assertGreater(total,0.)
        self.assertAlmostEqual(lost/(2*total),1.,places=10)
        self.assertAlmostEqual(model.F.sum()/total,1.,places=10)
        np.testing.assert_allclose(model.totals(),nuclei,rtol=1e-14)
        self.assertAlmostEqual(model.total_energy()/energy,1.,places=10)

    def test_outside_kT_domain_is_an_explicit_gap(self):
        cfg=configuration(['RD82:He3S_H:PI']);cfg['heavy_initial']={s['element']:dict(Ti_K=20.) for s in bank().states}
        model=build(cfg);self.assertTrue(any('диапазона' in s for s in model.chemistry.gaps))
        cfg['chemistry']['coverage_policy']='error';model=build(cfg)
        with self.assertRaisesRegex(ValueError,'Недостаточные'):model.chemistry.evolve(0,1e-9)

    def test_reaction_can_be_disabled_by_id(self):
        model=build(configuration([]));self.assertFalse(model.chemistry.reactions)

    def test_charge_or_nuclei_violation_is_rejected(self):
        with self.assertRaises(ValueError):Network([dict(id='H',label='H',charge=0,nuclei={'H':1}),dict(id='He',label='He',charge=0,nuclei={'He':1})],[dict(reactants={'H':2},products={'He':1})])


class ImportPhysics(unittest.TestCase):
    def test_native_atomic_source_has_a_real_common_energy_reference(self):
        import csv,io
        from atomic_core import Species
        low=dict(id='H_0:L001',species='H I',configuration='1s',n=1,multiplicity=2,L=0,g=2,energy_cm1=0.,energy_reference='ground_of_this_ion',mass_u=1.008,source_index=1,label='TEST H 1s',source='SYNTHETIC')
        high=dict(low,id='H_0:L002',n=2,L=1,g=6,energy_cm1=82303.,source_index=2)
        sigma=[dict(id='ion',kind='electron_ionization',initial=low['id'],final='parent',energy_eV=e,sigma_cm2=s,source='SYNTHETIC') for e,s in [(13.6,0.),(20.,1e-17)]]
        def encoded(rows,fields):
            f=io.StringIO();w=csv.DictWriter(f,fieldnames=fields,delimiter=';');w.writeheader();w.writerows(rows);return f.getvalue().encode()
        raw={'levels':encoded([low,high],LEVEL_FIELDS),'xs':encoded(sigma,XS_FIELDS)}
        cache=SimpleNamespace(get=lambda name,info:raw[info['fixture']])
        entry=dict(mass_u=1.008,ionization_eV=13.6,native_csv={'levels':'levels','electron_cross_sections':'xs'})
        registry=dict(sources={k:dict(fixture=k) for k in raw})
        levels,_,_,cross,_,meta=native_tables(cache,entry,registry,Species('H'),9)
        self.assertEqual(len(levels),2);self.assertEqual(len(cross),2);self.assertAlmostEqual(meta['parents'][0]['threshold_cm1']*HC_EV_CM,13.6)
        raw['levels']=encoded([dict(low,energy_reference='absolute'),high],LEVEL_FIELDS)
        with self.assertRaises(ValueError):native_tables(cache,entry,registry,Species('H'),9)

    def test_all_gui_presets_point_to_existing_current_configs(self):
        from launcher import PRESETS,configuration
        for name in PRESETS:
            cfg,engine=configuration(ROOT/name)
            self.assertEqual(cfg['schema'],'MixtureKineticSlab-2.1');self.assertEqual(engine,'mixture_kinetics.py')
            base=(ROOT/name).parent
            self.assertTrue((base/cfg['atomic_data']['element_directory']).is_dir())
            self.assertTrue((base/cfg['chemistry']['directory']).is_dir())

    def test_umist_multiple_ranges_and_colons_inside_quotes(self):
        text='74:AD:H-:H:H2:e-:::2:4.82e-09:0.02:4.3:10:100:M:A:"https://doi.org/test":"note:with:colons":4.32e-09:-0.39:39.4:101:3000:M:A:"https://doi.org/test":"second":\n'
        groups,report=parse_umist(text.encode(),['H']);self.assertFalse(report['rejected']);rows=groups['H']['rates']
        self.assertEqual(len(rows),2);self.assertIn('note:with:colons',rows[0]['notes'])
        self.assertEqual(groups['H']['reactions'][0]['products'],{'H2':1,'e':1})

    def test_photorates_are_not_imported_as_two_body_rates(self):
        text='1:PH:H2:PHOTON:H:H:::1:5.7e-11:0:4.18:10:41000:C:B:"doi":"note":\n'
        groups,_=parse_umist(text.encode(),['H']);self.assertEqual(len(groups['H']['reactions']),0);self.assertEqual(len(groups['H']['unsupported']),1)

    def test_network_import_stays_in_candidate_state(self):
        text='1:CE:H:H2+:H2:H+:::1:6.4e-10:0:0:10:41000:M:A:"doi":"note":\n'
        with tempfile.TemporaryDirectory() as folder:
            result=save_network(text.encode(),['H'],folder)
            document=json.loads((Path(folder)/'H/reactions.json').read_text(encoding='utf-8-sig'))
            self.assertIs(document['reviewed'],False);self.assertEqual(result['candidate_reactions'],1)

    def test_lxcat_requires_explicit_units_and_level_mapping(self):
        raw=b'EXCITATION\nH -> H(2p)\n10.2\n----------------\n10.2 0\n20 1e-20\n----------------\n'
        mapping=dict(projectile='electron',energy_scale_to_eV=1.,area_scale_to_cm2=1e4,channels=[dict(block_index=1,id='test',kind='electron_excitation',initial='H0',final='H2p')])
        rows,missing=convert_lxcat(raw,mapping,'SYNTHETIC');self.assertFalse(missing);self.assertAlmostEqual(rows[1]['sigma_cm2'],1e-16)
        mapping['channels'][0]['kind']='electron_ionization'
        with self.assertRaises(ValueError):convert_lxcat(raw,mapping,'SYNTHETIC')

    def test_registry_is_separate_per_element(self):
        registry=load_catalogue(ROOT,dict(element_directory='data/elements'))
        self.assertTrue({'H_0','He_0','Li_0'}<=set(registry['species']))
        self.assertTrue({'H','He','Li','Na'}<=set(registry['ground_stages']))


if __name__=='__main__':unittest.main()
