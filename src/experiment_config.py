"""Editable gas composition -> the existing kinetic initial conditions.

Only user input is translated here; no collision rates or new physical
operators are introduced. Gas fractions count heavy particles, not nuclei.
"""
from copy import deepcopy
import json
import math
from pathlib import Path
import re
from atomic_core import Species

KB_SI=1.380649e-23


def finite(value,label):
    try:x=float(str(value).strip().replace(',','.'))
    except (TypeError,ValueError):raise ValueError(label+': нужно число, например 1e10.') from None
    if not math.isfinite(x):raise ValueError(label+': нужно конечное число.')
    return x


def carrier_species(name,root,components_directory='../data/components'):
    """Resolve the element/stage inventory without loading numerical packages."""
    if not isinstance(name,str) or not name.strip():raise ValueError('Укажите компонент газа: He, H, He+, HeI_2_3S или точный ID.')
    name=name.strip()
    match=re.fullmatch(r'([A-Z][a-z]?)_(\d+):.+',name)
    if match:return [Species(match[1],int(match[2]))]
    if re.fullmatch(r'HeI_\d+_\d[SPDF]',name):return [Species('He')]
    try:return [Species.parse(name)]
    except ValueError:pass
    for p in sorted((Path(root)/components_directory).glob('*.json')):
        for component in json.loads(p.read_text(encoding='utf-8-sig')).get('components',[]):
            if component['id']==name:return [Species(e) for e in component['nuclei']]
    raise ValueError('Неизвестный компонент «'+name+'». Используйте He/H/Li, He+, Li II, HeI_2_3S или ID из банка данных.')


def gas_initial(gas,root,chemistry=None):
    allowed={'mode','temperature_K','total_density_cm3','pressure_Pa','components','distributions'}
    if set(gas)-allowed:raise ValueError('Неизвестные поля состава газа: '+', '.join(sorted(set(gas)-allowed)))
    temperature=finite(gas.get('temperature_K',300.),'Температура газа')
    if temperature<=0:raise ValueError('Температура газа должна быть >0 K.')
    mode=gas.get('mode','densities')
    if mode not in ('densities','fractions','pressure'):raise ValueError('Состав: densities, fractions или pressure.')
    items=gas.get('components',[])
    if not isinstance(items,list) or not items:raise ValueError('Добавьте хотя бы один компонент газа с ненулевой концентрацией.')
    initial={};inventory=[];fractions=0.
    for row in items:
        key='density_cm3' if mode=='densities' else 'fraction'
        if set(row)-{'species',key}:raise ValueError('Строка состава: species и '+key+'.')
        name=row.get('species','').strip();value=finite(row.get(key,''),'Количество '+name)
        if value<0:raise ValueError('Концентрация/доля '+name+' не может быть отрицательной.')
        if value==0:continue
        species=carrier_species(name,root,(chemistry or {}).get('components_directory','../data/components'))
        if name in initial:raise ValueError('Компонент «'+name+'» повторён. Измените его значение в существующей строке.')
        initial[name]=value;fractions+=value;inventory.extend(species)
    if not initial:raise ValueError('Все концентрации равны нулю. Укажите газ в таблице состава.')
    if mode!='densities':
        if abs(fractions-1.)>1e-6:raise ValueError(f'Сумма численных долей должна быть 1 (100%), сейчас {100*fractions:.8g}%. В редакторе есть «Нормировать до 100%».')
        total=finite(gas.get('total_density_cm3',''),'Полная концентрация') if mode=='fractions' else finite(gas.get('pressure_Pa',''),'Давление заполнения')/(KB_SI*temperature)*1e-6
        if total<=0:raise ValueError('Полная концентрация/давление заполнения должны быть >0.')
        initial={name:value*total for name,value in initial.items()}
    return initial,sorted(set(inventory),key=lambda s:(s.element,s.charge)),temperature


def prepare_configuration(config,root):
    out=deepcopy(config)
    from light_input import particle_inputs
    out['spatial_sources']=[particle_inputs(item) for item in out.get('spatial_sources',[])]
    if 'gas' not in out:return out
    initial,inventory,T=gas_initial(out['gas'],root,out.get('chemistry'))
    for source in out.get('spatial_sources',[]):
        if source.get('enabled',True) and source.get('mode') in ('ions','particles'):
            inventory.extend(carrier_species(source.get('species',''),root,out.get('chemistry',{}).get('components_directory','../data/components')))
    inventory=sorted(set(inventory),key=lambda s:(s.element,s.charge))
    out['initial_cm3']=initial
    out.setdefault('atomic_data',{})['species']=[s.label for s in inventory]
    # A common gas temperature is an explicit input; no drift is introduced.
    distributions=out['gas'].get('distributions',{})
    if not isinstance(distributions,dict):raise ValueError('gas.distributions: объект компонент -> распределение.')
    common=deepcopy(distributions.get('default',{'Ti_K':T}))
    if common.get('kind','maxwellian')=='maxwellian':common['Ti_K']=T
    out['heavy_initial']={s.element:deepcopy(common) for s in inventory}
    for name in initial:
        try:Species.parse(name)
        except ValueError:
            if ':' not in name and not name.startswith('HeI_'):out['heavy_initial'][name]=deepcopy(common)
    for name,options in distributions.items():
        if name=='default':continue
        if not isinstance(options,dict):raise ValueError('Распределение '+name+': нужен объект параметров.')
        carrier_species(name,root,out.get('chemistry',{}).get('components_directory','../data/components'))
        try:key=Species.parse(name).label
        except ValueError:key=name
        out['heavy_initial'][key]=deepcopy(options)
    out.setdefault('radiation',{})['Tn_K']=T
    elements={s.element for s in inventory}
    options=out.setdefault('level_selection',{})
    options['forced_levels']=[name for name in options.get('forced_levels',[]) if any(name.startswith(e+'_') or name.startswith(e+'I_') for e in elements)]
    out.setdefault('diagnostics',{})['reaction_levels']=list(initial)
    return out


def editor_composition(config):
    """Return density rows for old JSON and a gas object for the new editor."""
    if 'gas' in config:return deepcopy(config['gas'])
    temperatures=[v.get('Ti_K',300.) for v in config.get('heavy_initial',{}).values()]
    return dict(mode='densities',temperature_K=temperatures[0] if temperatures else 300.,
                components=[dict(species=name,density_cm3=value) for name,value in config.get('initial_cm3',{}).items()])
