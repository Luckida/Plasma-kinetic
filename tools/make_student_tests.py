"""Student physics configurations; run from any directory."""
from copy import deepcopy
from pathlib import Path
import json,math

WORK=Path(__file__).resolve().parent.parent
MODEL=WORK
FLAGS=['electron_bound_bound','electron_ionization','spontaneous','two_photon',
       'photo_bound_bound','photo_bound_free','radiative_recombination',
       'stimulated_recombination','penning','associative_ionization',
       'charge_exchange','dissociative_recombination','mutual_neutralization',
       'heavy_chemistry','heavy_quenching']
BASE={
 'schema':'MixtureKineticSlab-2.1',
 'atomic_data':{'species':['H'],'nmax':9,'element_directory':'../data/elements',
  'prepared_directory':'../data/elements','cache_directory':'../cache/atomic_sources','offline':True},
 'gas':{'mode':'densities','temperature_K':300.,'components':[{'species':'H','density_cm3':1e10}],
  'distributions':{'default':{'kind':'maxwellian','Ti_K':300.,'match_moments':True}}},
 'electron_initial':{'kind':'maxwellian','Te_eV':1.,'density_cm3':0.,'match_moments':True},
 'phase_grid':{'energy_min_eV':1e-4,'energy_max_eV':200.,'energy_points':128,'angular_points':8},
 'heavy_phase_grid':{'energy_min_eV':1e-6,'energy_max_eV':50.,'energy_points':64,'angular_points':8},
 'level_selection':{'mode':'all','relative_tolerance':1e-6,'forced_levels':[],'missing_data':'report'},
 'cells':1,'thickness_cm':1.,'particle_boundary':{'kind':'periodic'},
 'electric_field':{'mode':'prescribed','Ez_V_cm':0.},
 'radiation':{'enabled':False,'scale':0.,'Tn_K':300.,'continuum_points':80,'line_points':21,'angular_points':4,'electron_scattering':False},
 'processes':{**{key:False for key in FLAGS},'disabled_ids':[]},
 'collisions':{'ee_bgk':False,'ii_bgk':False,'ei_binary':False,'neutral_elastic':False,
  'coulomb_log':10.,'maximum_phase_loss':.2,'relative_azimuth_points':2,'scattering_azimuth_points':2,'neutral_low_energy':'error'},
 'chemistry':{'enabled':False,'rate_policy':'cross_sections_only','directory':'../data/reactions',
  'components_directory':'../data/components','allow_estimated_rates':False,'disabled_ids':[]},
 'spatial_sources':[],'probe':{'mode':'off'},
 'performance':{'workers':1,'radiation_backend':'numpy','collision_cache_mb':128},
 'solver':{'t_end_s':1e-6,'max_step_s':1e-8,'min_step_s':1e-18,'samples':51,'max_steps':1000000},
 'diagnostics':{'output_mode':'compact','plots_on_finish':False,'phase_history_csv':False,'reaction_levels':[]},
 'progress':{'mode':'text','interval_s':.5},'output_directory':'../outputs/study_base'
}

def merge(base,changes):
 out=deepcopy(base)
 for key,value in changes.items():
  out[key]=merge(out[key],value) if isinstance(value,dict) and isinstance(out.get(key),dict) else deepcopy(value)
 return out

def gas(**values):return [{'species':key,'density_cm3':value} for key,value in values.items()]

CASES=[]
def add(name,title,changes,keep=(),group=''):
 cfg=merge(BASE,changes);cfg['output_directory']='../outputs/'+name
 ei=cfg['electron_initial']
 if ei.get('kind')=='maxwellian':
  ei.pop('energy_distribution_file',None);ei.pop('energy_form',None)
 elif ei.get('kind')=='energy_table':
  ei.pop('Te_eV',None);ei.pop('match_moments',None)
 CASES.append({'name':name,'title':title,'config':cfg,'changes':changes,'keep':list(keep),'group':group})
 return cfg

A=627000000.;TAU=1/A
add('study_01_decay','Спонтанный распад H 2p',{
 'gas':{'components':gas(**{'H_0:L003':1e8})},'processes':{'spontaneous':True},
 'solver':{'t_end_s':5*TAU,'max_step_s':TAU/100,'min_step_s':1e-20,'samples':51}},
 ['A:A:H_0:L003:H_0:L001'],'T01')
add('study_02_energy','Возбуждение и девозбуждение H',{
 'gas':{'components':gas(**{'H':1e10,'H+':1e8})},
 'electron_initial':{'Te_eV':20.,'density_cm3':1e8},'processes':{'electron_bound_bound':True},
 'solver':{'t_end_s':1e-3,'max_step_s':1e-5,'min_step_s':1e-16,'samples':51}},
 ['CE:H_0:L001->H_0:L003','CE:H_0:L003->H_0:L001'],'T02')
add('study_02_reverse','Обратный контроль H 2p',{
 'gas':{'components':gas(**{'H_0:L003':1e10,'H+':1e8})},
 'electron_initial':{'Te_eV':20.,'density_cm3':1e8},'processes':{'electron_bound_bound':True},
 'solver':{'t_end_s':1e-3,'max_step_s':1e-5,'min_step_s':1e-16,'samples':51}},
 ['CE:H_0:L001->H_0:L003','CE:H_0:L003->H_0:L001'],'T02')

ION_CHANGES={'gas':{'components':gas(**{'H':1e11,'H+':1e8})},
 'electron_initial':{'kind':'energy_table','energy_distribution_file':'study_e_above.csv','energy_form':'probability','density_cm3':1e8},
 'processes':{'electron_ionization':True},'solver':{'t_end_s':1e-4,'max_step_s':1e-6,'min_step_s':1e-16,'samples':51}}
add('study_03_below','Электроны ниже порога',merge(ION_CHANGES,{'electron_initial':{'energy_distribution_file':'study_e_below.csv'}}),['ION:H_0:L001'],'T03')
ion=add('study_03_above','Электроны выше порога',ION_CHANGES,['ION:H_0:L001'],'T03')
maxw=merge(ION_CHANGES,{'electron_initial':{'kind':'maxwellian','Te_eV':14.,'match_moments':True}})
add('study_03_maxwell','Та же средняя энергия при максвелловской форме',maxw,['ION:H_0:L001'],'T03')
for label,Te,ne,parts in [('cold',1.,0.,{'H':1e11}),('hot',30.,0.,{'H':1e11}),('control',30.,1e8,{'H':1e11,'H+':1e8})]:
 add('study_04_'+label,'Нулевое число электронов '+label,{
  'gas':{'components':gas(**parts)},'electron_initial':{'kind':'maxwellian','Te_eV':Te,'density_cm3':ne},
  'phase_grid':{'energy_max_eV':300.},'processes':{'electron_ionization':True},
  'solver':{'t_end_s':1e-4,'max_step_s':1e-6,'min_step_s':1e-16,'samples':51}},['ION:H_0:L001'],'T04')
for label,n,L in [('base',1e11,1.),('density',2e11,1.),('length',1e11,2.)]:
 add('study_05_'+label,'Поглощение '+label,{
  'atomic_data':{'species':['He']},'gas':{'components':gas(**{'He_0:L002':n})},'thickness_cm':L,
  'probe':{'mode':'flat','minimum_nm':1083.28,'maximum_nm':1083.35,'points':1001,'mu':1.},
  'solver':{'t_end_s':1e-9,'max_step_s':1e-10,'min_step_s':1e-20,'samples':2}},group='T05')

for index,dt in enumerate([1e-6,5e-7,2.5e-7],1):
 add('study_06_dt'+str(index),'Сходимость по времени '+str(index),merge(ION_CHANGES,{'solver':{'max_step_s':dt}}),['ION:H_0:L001'],'T06')
for points in [64,128,256]:
 add('study_06_e'+str(points),'Сходимость по энергии '+str(points),merge(ION_CHANGES,{'phase_grid':{'energy_points':points},'solver':{'max_step_s':2.5e-7}}),['ION:H_0:L001'],'T06')
for cells in [16,32,64]:
 parts=[]
 for z in range(cells):
  x=(z+.5)/cells;sinc=math.sin(math.pi/cells)/(math.pi/cells)
  parts.append({'initial_cm3':{'H':1e10*(1+.1*sinc*math.cos(2*math.pi*x))}})
 add('study_06_z'+str(cells),'Свободный перенос '+str(cells)+' ячеек',{
  'cells':cells,'gas':{'temperature_K':100.,'components':gas(H=1e10),
   'distributions':{'default':{'Ti_K':100.,'uz_cm_s':1e5}}},
  'cell_overrides':parts,'heavy_phase_grid':{'energy_min_eV':1e-6,'energy_max_eV':.5,'energy_points':128,'angular_points':16},
  'radiation':{'Tn_K':100.},'solver':{'t_end_s':2e-6,'max_step_s':2e-8,'min_step_s':1e-18,'samples':5}},group='T06')
for label,bound in [('manual50',50.),('manual300',300.),('auto','auto')]:
 add('study_07_'+label,'Горячий газ '+label,{
  'atomic_data':{'species':['H','He']},'gas':{'temperature_K':300000.,'components':gas(H=1e9,He=9e9),
   'distributions':{'default':{'Ti_K':300000.}}},'electron_initial':{'Te_eV':30.,'density_cm3':0.},
  'phase_grid':{'energy_max_eV':300.},'heavy_phase_grid':{'energy_max_eV':bound},'radiation':{'Tn_K':300000.},
  'solver':{'t_end_s':2e-7,'max_step_s':2e-9,'samples':5}},group='T07')
for T in [250.,300.,1000.,4500.,5000.]:
 add('study_08_T'+str(int(T)),'Диапазон химического коэффициента '+str(int(T))+' K',{
  'atomic_data':{'species':['H','He']},'gas':{'temperature_K':T,'components':gas(**{'H':1e8,'He_0:L002':1e8}),
   'distributions':{'default':{'Ti_K':T}}},'radiation':{'Tn_K':T},
  'chemistry':{'enabled':True,'rate_policy':'maxwellian_moments','datasets':{'He_H_chemiionization':'RD82'}},
  'processes':{'penning':True,'associative_ionization':True},
  'solver':{'t_end_s':1e-3,'max_step_s':1e-5,'min_step_s':1e-16,'samples':51}},
  ['RD82:He3S_H:PI','RD82:He3S_H:AI'],'T08')

def write_cases(isolate=True):
 import sys
 sys.path.insert(0,str(MODEL/'src'))
 from mixture_kinetics import MixtureSlab
 folder=MODEL/'configs'
 (folder/'study_e_below.csv').write_text('energy_eV;P_per_eV\n5;0\n6;1\n7;0\n', encoding='utf-8')
 (folder/'study_e_above.csv').write_text('energy_eV;P_per_eV\n20;0\n21;1\n22;0\n', encoding='utf-8')
 summary=[]
 for case in CASES:
  cfg=deepcopy(case['config']);keep=set(case['keep'])
  if isolate and keep:
   m=MixtureSlab(cfg,folder);available={r['id'] for r in m.channels}
   if not keep<=available:raise ValueError(case['name']+': missing '+str(keep-available))
   cfg['processes']['disabled_ids']=sorted(available-keep)
  path=folder/(case['name']+'.json');path.write_text(json.dumps(cfg,ensure_ascii=False,indent=2)+'\n', encoding='utf-8')
  summary.append({'name':case['name'],'kept':case['keep'],'disabled':len(cfg['processes']['disabled_ids'])})
 (WORK/'validation'/'case_inventory.json').write_text(json.dumps(summary,indent=2), encoding='utf-8')
 return summary

if __name__=='__main__':
 print(json.dumps(write_cases(),indent=2))
