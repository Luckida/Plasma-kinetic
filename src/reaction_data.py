"""Common equations and provenance for every atomic and chemical channel.

This module exposes supplied data and the already implemented rate laws.
It does not estimate missing coefficients or enable an alternative dataset.
"""
import json
from pathlib import Path
import numpy as np

LAWS={
    'excitation':'Свёртка σ(E)·v с текущей EEDF; умножение на населённость исходного уровня.',
    'deexcitation':'Обратное электронное сечение по микроскопическому детальному балансу; свёртка с EEDF.',
    'ionization':'Электронное σ(E)·v, свёрнутое с текущей EEDF; один дополнительный электрон.',
    'spontaneous':'R = A·n верхнего уровня; A в с⁻¹.',
    'two_photon':'R = A·n верхнего уровня; на событие испускаются два фотона.',
    'photoexcitation':'Сечение линии из A и статистических весов; интеграл по текущему полю излучения.',
    'stimulated_emission':'Сечение вынужденного перехода из A; интеграл по текущему полю излучения.',
    'photoionization':'Интеграл σ(ν)·4πJν/(hν) по частотам; порог задан разностью энергий.',
    'radiative_recombination':'Обратное фотоионизационное сечение по соотношению Милна; свёртка с EEDF.',
    'stimulated_recombination':'Рекомбинация по Милну с учётом текущего поля излучения.'}


def equation(states,r):
    if r.get('chemistry'):
        products=dict(r['products'])
        if r.get('radiative'):products['photon']=1
        return dict(reactants=dict(r['reactants']),products=products)
    a=states[r['i']]['id'];b=states[r['j']]['id'];kind=r['kind'];left={a:1};right={b:1}
    if kind in ('excitation','deexcitation'):left['e']=right['e']=1
    elif kind=='ionization':left['e']=1;right['e']=2
    elif kind in ('radiative_recombination','stimulated_recombination'):
        left['e']=1;right['photon']=1
        if kind=='stimulated_recombination':left['photon']=1;right['photon']=2
    elif kind in ('spontaneous','two_photon'):right['photon']=2 if kind=='two_photon' else 1
    elif kind in ('photoionization','photoexcitation'):left['photon']=1;right.update({'e':1} if kind=='photoionization' else {})
    elif kind=='stimulated_emission':left['photon']=1;right['photon']=2
    return dict(reactants=left,products=right)


def kinetics_summary(r):
    data=r.get('kinetics',{});lines=['Расчёт скорости: '+data.get('law',LAWS.get(r['kind'],'Свёртка сечений по относительным скоростям либо опубликованный k(T); выбранный способ сохраняется в новых результатах.'))]
    if 'energy_interval_eV' in data:lines.append('Диапазон исходной σ(E), эВ: '+' … '.join(f'{x:.8g}' for x in data['energy_interval_eV'])+'; единицы сечения: см².')
    if 'A_s' in data:lines.append(f"A = {data['A_s']:.9g} с⁻¹.")
    if data.get('rate_id'):lines.append('Набор k(T): '+data['rate_id']+'; единицы коэффициента: см³/с.')
    for fit in data.get('fits',[]):
        lines.append('k(T): '+', '.join(f'{key}={value}' for key,value in fit.items() if key not in ('source','id')))
        if fit.get('source'):lines.append('Источник коэффициента: '+fit['source'])
    if data.get('coefficient_convention'):lines.append('Конвенция коэффициента: '+data['coefficient_convention'])
    if data.get('estimated'):lines.append('ОЦЕНОЧНЫЙ набор; не является измеренным полным состояниеразрешённым сечением.')
    if data.get('model_note'):lines.append('Допущение / ограничение: '+data['model_note'])
    if data.get('file'):lines.append('Таблица / описание: '+data['file'])
    return '\n'.join(lines)


def model_reaction_data(model):
    """Selected network plus inactive data-backed candidates for its elements."""
    m=model;rows={};chem={r['id']:r for r in m.chemistry.reactions}
    def emit(r,kinetics=None):
        row=dict(id=r['id'],kind=r['kind'],source=r.get('source',''),**equation(m.states,r))
        row['kinetics']=kinetics or dict(law=LAWS.get(r['kind'],'Опубликованный k(T) или интеграл σ(Erel)·vrel по распределениям участников.'))
        if r['id'] not in rows:rows[r['id']]=row
        return rows[r['id']]
    for r in m.channels:
        data=dict(law=LAWS.get(r['kind'],'Опубликованный k(T) или интеграл σ(Erel)·vrel по распределениям участников.'))
        if 'A' in r:data['A_s']=float(r['A'])
        if r['id'] in chem:
            raw=chem[r['id']];data.update({key:raw[key] for key in ('rate_id','cross_section_id','fits','model_note','estimated','coefficient_convention') if key in raw})
            data['law']='Интеграл σ(Erel)·vrel по обоим фазовым распределениям.' if 'cross' in raw else 'R = k(Trel)·nA·nB с принятой конвенцией одинаковых частиц; Trel вычисляется из моментов. Это тепловое замыкание, а не восстановленное σ(E).'
            if 'cross' in raw:data['energy_interval_eV']=[float(raw['cross'][0,0]),float(raw['cross'][-1,0])]
        emit(r,data)
    for (a,b,kind),table in m.data.cross.items():
        if kind not in ('excitation','deexcitation','ionization'):continue
        ident='ION:'+a if kind=='ionization' else 'CE:'+a+'->'+b
        data=dict(law=LAWS[kind],energy_interval_eV=[float(table[0,0]),float(table[-1,0])])
        row=emit(dict(id=ident,kind=kind,i=m.index[a],j=m.index[b],source=m.data.cross_source[(a,b,kind)]),data)
        row['kinetics'].update(data)
    for line in m.data.lines:
        two=line.get('spectral_kind')=='two_photon';kinds=[('A:','two_photon' if two else 'spontaneous')]
        if not two and line.get('pump',True):kinds.extend([('PUMP:','photoexcitation'),('STIM:','stimulated_emission')])
        for prefix,kind in kinds:
            i,j=(m.index[line['lower']],m.index[line['upper']]) if kind=='photoexcitation' else (m.index[line['upper']],m.index[line['lower']])
            data=dict(law=LAWS[kind],A_s=float(line['A_s']))
            emit(dict(id=prefix+line['id'],kind=kind,i=i,j=j,source=line['source']),data)['kinetics'].update(data)
    for i,j in m.data.parent_index.items():
        ident=m.states[i]['id']
        if ident not in m.photo.tables:continue
        table=m.photo.tables[ident]
        for prefix,kind in [('PI:','photoionization'),('RR:','radiative_recombination'),('SR:','stimulated_recombination')]:
            a,b=(i,j) if kind=='photoionization' else (j,i)
            emit(dict(id=prefix+ident,kind=kind,i=a,j=b,source=m.photo.sources[ident]),dict(law=LAWS[kind],energy_interval_eV=[float(table[0,0]),float(table[-1,0])]))
    folder=m.root/m.cfg.get('chemistry',{}).get('directory','data/reactions')
    from mixture_chemistry import records
    for file in sorted(folder.glob('*/reactions.json')):
        doc=json.loads(file.read_text(encoding='utf-8-sig'))
        if doc.get('reviewed') is not True or not set(doc.get('elements',[]))<=set(m.data.elements):continue
        fits={};cross={}
        if (file.parent/'rate_coefficients.csv').exists():
            for fit in records(file.parent/'rate_coefficients.csv'):fits.setdefault(fit['id'],[]).append(fit)
        if (file.parent/'cross_sections.csv').exists():
            for point in records(file.parent/'cross_sections.csv'):cross.setdefault(point['id'],[]).append(float(point['energy_eV']))
        for raw in doc.get('reactions',[]):
            try:sides={side:{m.chemistry.resolve(name,doc.get('selectors',{})):n for name,n in raw[side].items()} for side in ('reactants','products')}
            except (ValueError,KeyError):sides={side:dict(raw[side]) for side in ('reactants','products')}
            data={key:raw[key] for key in ('rate_id','cross_section_id','model_note','estimated','coefficient_convention') if key in raw}
            data['file']=str(file.relative_to(m.root)) if file.is_relative_to(m.root) else file.name
            if raw.get('rate_id') in fits:data['fits']=fits[raw['rate_id']]
            if raw.get('cross_section_id') in cross:data['energy_interval_eV']=[min(cross[raw['cross_section_id']]),max(cross[raw['cross_section_id']])]
            if raw['id'] in rows:rows[raw['id']]['kinetics'].update(data)
            else:emit(dict(raw,**sides,chemistry=True),data)
    active={r['id'] for r in m.channels};disabled=set(m.cfg.get('processes',{}).get('disabled_ids',[]))|set(m.chemistry.cfg.get('disabled_ids',[]))
    for ident,row in rows.items():
        row['calculated']=ident in active
        row['status']='В расчётной сети' if ident in active else 'Отключена пользователем' if ident in disabled else 'Оценочный набор отключён' if row['kinetics'].get('estimated') and not m.chemistry.cfg.get('allow_estimated_rates',False) else 'Есть в банке; в эту сеть не включена (флаги, альтернативный набор или ограничения данных)'
    return list(rows.values())
