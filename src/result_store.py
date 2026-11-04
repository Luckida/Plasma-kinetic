"""One plotting store and four readable tables; no per-cell diagnostic tree."""
from project_version import VERSION
from pathlib import Path
import csv
import hashlib
import json
import os
import numpy as np
from light_input import light_enabled,light_spectrum

AMU=1.66053906892e-24
FILE_STORE='данные_для_графиков.npz'


def state_caption(state):
    if state.get('molecular'):return state['id']
    element=state['element'];charge=int(state['charge'])
    if element not in ('H','He','Li','Na'):return state['label']
    from atomic_core import SYMBOLS
    if charge==SYMBOLS.index(element)+1:return element+('⁺' if charge==1 else str(charge).translate(str.maketrans('0123456789','⁰¹²³⁴⁵⁶⁷⁸⁹'))+'⁺')
    n=int(state.get('n',0));L=int(state.get('L',0));mult=int(state.get('multiplicity',1))
    if n<=0 or not 0<=L<len('SPDFGHIKLM'):return state['label']
    return f"{state['species']} {n}"+str(mult).translate(str.maketrans('0123456789','⁰¹²³⁴⁵⁶⁷⁸⁹'))+'SPDFGHIKLM'[L]


def write_rows(path,fields,rows):
    with Path(path).open('w',encoding='utf-8-sig',newline='') as stream:
        writer=csv.DictWriter(stream,fields,delimiter=';');writer.writeheader();writer.writerows(rows)


def pretty_side(side,labels):
    return ' + '.join((str(count)+' ' if count!=1 else '')+labels.get(name,name) for name,count in side.items()) or '∅'


def top_indices(values,n=5):
    values=np.asarray(values);return [int(k) for k in np.argsort(-values,kind='stable')[:n] if values[k]>0]


def remove_old_generated_files(folder):
    """Remove only named legacy outputs, preserving unrelated user files."""
    names='populations_by_cell particle_moments conservation component_budgets reaction_top5 collision_diagnostics electron_energy_history electron_phase_history heavy_phase_history collision_partners states level_selection selected_energy_limits upper_data_boundary excitation_energy_diagnostics reaction_catalogue lines_final spectrum_final data_gaps chemical_reactions chemical_event_history chemical_domain_audit chemical_known_omissions pair_coverage spatial_source_history wall_history spatial_profiles spatial_kinetics lines_by_position strongest_lines_by_position spectra_by_position_final probe_transmission probe_history probe_by_position_final probe_lines'.split()
    names=[name+'.csv' for name in names]+['quality.json','input_config.json','resolved_phase_grid.json','atomic_bank_provenance.json','run_metadata.json','chemical_data_provenance.json','probe_input_metadata.json','mixture_history.npz','spatial_radiation_history.npz','summary.txt']
    names+=['populations.png','electron_eedf.png','ion_velocity.png','spatial_profiles.png','spatial_lines.png','spatial_chemistry.png','probe_transmission.png']
    for name in names:
        path=folder/name
        if path.is_file():path.unlink()
    for cell in folder.glob('cell_[0-9][0-9]'):
        dynamic=cell/'reaction_dynamics'
        if dynamic.is_dir():
            for path in dynamic.iterdir():
                if path.is_file() and (path.name=='files.csv' or path.name.endswith(('_rates.csv','_fractions.csv'))):path.unlink()
            if not any(dynamic.iterdir()):dynamic.rmdir()
        if cell.is_dir() and not any(cell.iterdir()):cell.rmdir()


def source_metadata(model):
    rows=[];root=model.root;probe=model.cfg.get('probe');cfg=model.cfg.get('radiation',{})
    candidates=[dict(cfg,id='Свет на границе',mode='radiation',position_cm=0. if cfg.get('beam_side','left')=='left' else model.length)] if light_enabled(cfg) else []
    candidates.extend(item for item in model.sources.items)
    for item in candidates:
        mode=item.get('mode');row={key:item.get(key,'') for key in ('id','mode','position_cm','density_cm3','velocity_cm_s','species')}
        row.update(start_s=float(item.get('start_s',0.)),stop_s=float(item.get('stop_s',model.end)))
        if mode=='radiation':
            data=light_spectrum(item,root,probe);row.update(data.metadata)
            if item.get('normalization')!='integral':
                x=np.asarray(data.wavelength_A);y=np.asarray(data.values);mu=float(item.get('mu',item.get('beam_mu',1.)));scale=float(item.get('scale',1.))
                from spectrum_input import HC_EV_A,EV_ERG
                power=float(np.trapezoid(y,x))*mu*scale*1e-7;photons=float(np.trapezoid(y/(HC_EV_A/x*EV_ERG),x))*mu*scale
                row.update(normal_power_W_cm2=power,photon_flux_cm2_s=photons,minimum_photon_eV=float(HC_EV_A/x[-1]),maximum_photon_eV=float(HC_EV_A/x[0]),mean_photon_eV=power*1e7/photons/EV_ERG if photons else 0.)
            row['pulse_energy_J_cm2']=row['normal_power_W_cm2']*(row['stop_s']-row['start_s'])
        else:row.update(rate_cm2_s=float(item.get('rate_cm2_s',0.)),mean_particle_energy_eV=float(item.get('energy_eV',0.)))
        # No phase arrays or absolute local paths in the user-facing source list.
        rows.append({k:v for k,v in row.items() if k not in ('file','headers')})
    return rows


def export_results(model,folder):
    m=model;folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
    m.progress.start('Результаты','Сохранение четырёх таблиц и данных для графиков')
    snaps=m.snapshots;times=np.asarray([s['time_s'] for s in snaps]);n=np.asarray([s['n'] for s in snaps]);events=np.asarray([s['events'] for s in snaps]);integrals=np.asarray([s['integral'] for s in snaps])
    state_fields=('id','label','configuration','element','species','charge','mass_u','n','L','multiplicity','g','energy_cm1','molecular','nuclei')
    states=[{k:s.get(k,False if k=='molecular' else {s['element']:1} if k=='nuclei' else '') for k in state_fields} for s in m.states]
    captions=[state_caption(s) for s in states]
    for s,caption in zip(states,captions):s['display_label']=caption if captions.count(caption)==1 else caption+' · '+s['id']
    labels={s['id']:s['display_label'] for s in states};labels.update(e='e⁻',photon='hν')
    from reaction_data import model_reaction_data
    reaction_bank=model_reaction_data(m);bank_by_id={r['id']:r for r in reaction_bank}
    for row in reaction_bank:
        row['left']=pretty_side(row['reactants'],labels);row['right']=pretty_side(row['products'],labels);row['reaction']=row['left']+' → '+row['right']
    reactions=[]
    for r in m.channels:
        left=pretty_side(r['reactants'],labels);right=pretty_side(r['products'],labels)
        reactions.append(dict(id=r['id'],kind=r['kind'],reactants=r['reactants'],products=r['products'],left=left,right=right,reaction=left+' → '+right,source=r.get('source',''),kinetics=bank_by_id[r['id']]['kinetics']))
    residual=np.asarray([s['n']-m.initial-np.einsum('sk,zk->zs',m.matrix,s['integral'])-s['transport']-s['source_heavy'] for s in snaps])
    e_residual=np.asarray([s['F'].sum(axis=(1,2))-m.initial_F.sum(axis=(1,2))-s['integral']@m.electron_change-s['electron_transport']-s['source_electrons'] for s in snaps])
    norm=max(float(m.initial_nuclei.sum()),1.);mass_weights=np.array([next(s['mass_u'] for s in m.states if s['element']==e)*AMU for e in m.data.elements])
    massnorm=max(float(m.initial_nuclei@mass_weights),1e-300);mass=np.array([s['nuclei']@mass_weights for s in snaps]);charge=np.array([s['charge'] for s in snaps]);energy=np.array([s['energy'] for s in snaps]);enorm=max(abs(m.initial_energy),1e-100)
    kinetic=np.maximum(np.max(abs(residual),axis=(1,2)),np.max(abs(e_residual),axis=1))/max(float(m.initial.sum()),float(m.initial_F.sum()),1.)
    nuclei=max(float(abs(s['nuclei']).max()) for s in snaps)/norm
    m.data_issues.update(m.chemistry.gaps);m.model_notes.update(m.chemistry.notes)
    values=dict(charge=float(max(abs(charge))/norm),mass=float(max(abs(mass))/massnorm),nuclei=nuclei,energy=float(max(abs(energy))/enorm),kinetic=float(max(kinetic)),reaction_photon=max(s['reaction_photon_residual'] for s in snaps))
    quality=dict(version=VERSION,atomic_data_status='incomplete' if m.data_issues else 'available_channels_checked',data_issues=sorted(m.data_issues),declared_data_models=sorted(m.model_notes),accepted_steps=m.steps,rejected_steps=m.rejected,elapsed_s=m.elapsed,
        spatial_particle_sources=sum(r['mode']!='radiation' for r in m.sources.items),spatial_photon_sources=sum(r['mode']=='radiation' for r in m.sources.items))
    for key,value in values.items():quality[key+'_status']='passed' if value<(1e-6 if key=='energy' else 1e-8 if key=='reaction_photon' else 1e-9) else 'failed';quality[key+'_max_relative_residual']=value
    conservation=np.column_stack([mass/massnorm,charge/norm,kinetic,energy/enorm])
    write_rows(folder/'сохранение.csv',['Время, с','Невязка массы / начальная масса','Невязка заряда / число исходных ядер','Максимальная невязка кинетики / исходная концентрация','Невязка энергии / начальная энергия'],
        (dict(zip(['Время, с','Невязка массы / начальная масса','Невязка заряда / число исходных ядер','Максимальная невязка кинетики / исходная концентрация','Невязка энергии / начальная энергия'],[float(t),*map(float,row)])) for t,row in zip(times,conservation)))
    pop_fields=['Время, с','Ячейка','Расстояние, см','Компонент','ID уровня','Уровень / реагент','Ядерный состав','Заряд','Концентрация, см^-3','Невязка кинетики, см^-3']
    def populations():
        for a,t in enumerate(times):
            for z in range(m.cells):
                for j,s in enumerate(states):yield dict(zip(pop_fields,[float(t),z+1,(z+.5)*m.dz,s['species'],s['id'],s['display_label'],json.dumps(s['nuclei'],ensure_ascii=False),s['charge'],float(n[a,z,j]),float(residual[a,z,j])]))
                yield dict(zip(pop_fields,[float(t),z+1,(z+.5)*m.dz,'e⁻','e','Электроны','{}',-1,float(snaps[a]['F'][z].sum()),float(e_residual[a,z])]))
    write_rows(folder/'кинетика_уровней.csv',pop_fields,populations())
    rate_fields=['Время, с','Ячейка','Расстояние, см','ID реакции','Реакция','Реагенты','Продукты','Скорость событий, см^-3 с^-1','Интеграл событий, см^-3']
    def reaction_rows():
        for a,t in enumerate(times):
            for z in range(m.cells):
                for k,r in enumerate(reactions):yield dict(zip(rate_fields,[float(t),z+1,(z+.5)*m.dz,r['id'],r['reaction'],r['left'],r['right'],float(events[a,z,k]),float(integrals[a,z,k])]))
    write_rows(folder/'реакции_по_времени.csv',rate_fields,reaction_rows())
    groups={name:[j for j,s in enumerate(states) if s['species']==name] for name in sorted({s['species'] for s in states})}
    targets=[(s['id'],s['display_label'],m.matrix[j]) for j,s in enumerate(states)]+[('компонент:'+name,name,m.matrix[ids].sum(axis=0)) for name,ids in groups.items()]+[('e','e⁻',m.electron_change)]
    top_fields=['Время, с','Ячейка','Цель','Направление','Место','ID реакции','Реакция','Реагенты','Продукты','Скорость для цели, см^-3 с^-1','Доля от всех процессов цели','Покрытие top-5']
    def rankings():
        for a,t in enumerate(times):
            for z in range(m.cells):
                scopes=[('Все реакции','Все реакции',None)]+targets
                for ident,label,stoich in scopes:
                    directions=[('события',np.ones(m.nr))] if stoich is None else [('заселение',np.maximum(stoich,0)),('расселение',np.maximum(-stoich,0))]
                    for direction,weights in directions:
                        speeds=events[a,z]*weights;selected=top_indices(speeds);total=float(speeds.sum());coverage=float(speeds[selected].sum()/total) if total else 0.
                        for rank,k in enumerate(selected,1):
                            r=reactions[k];yield dict(zip(top_fields,[float(t),z+1,label,direction,rank,r['id'],r['reaction'],r['left'],r['right'],float(speeds[k]),float(speeds[k]/total),coverage]))
    write_rows(folder/'топ5_по_времени.csv',top_fields,rankings())
    sources=source_metadata(m);gh=next(iter(m.gh.values()));probe=m.probe
    metadata=dict(version=VERSION,states=states,reactions=reactions,reaction_bank=reaction_bank,components=groups,quality=quality,sources=sources,config=m.cfg,
        chemical_omissions=m.chemistry.omissions,chemical_domain_audit=m.chemistry.domain_audit,atomic_provenance=m.data.manifests,photo_sources=m.photo.sources,
        limitations=['Finite atomic/reaction database; see quality.data_issues.','Scalar LS populations and fixed-temperature Voigt profiles.','Weak probe: bound-bound transmission; no continuum absorption or added spontaneous light.','Particle beam density n and signed mean normal speed u prescribe column injection n*abs(u).','Passed balances do not establish grid convergence or complete physical data.'],
        source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in Path(__file__).parent.glob('*.py')})
    arrays=dict(metadata_json=np.asarray(json.dumps(metadata,ensure_ascii=False)),time_s=times,position_cm=(np.arange(m.cells)+.5)*m.dz,n_cm3=n,electrons_cm3=np.asarray([s['F'] for s in snaps]),heavy_cm3=np.asarray([s['P'] for s in snaps]),
        energy_eV=m.ge.e,energy_weights_eV=m.ge.we,electron_vz_cm_s=m.ge.vz,heavy_energy_eV=gh.e,heavy_mu=gh.mu,
        reaction_rates=events,reaction_integrals=integrals,stoichiometry=m.matrix,electron_stoichiometry=m.electron_change,kinetic_residual_cm3=residual,kinetic_residual_e_cm3=e_residual,conservation_relative=conservation,
        photon_nu_Hz=m.nu,emission=np.asarray([s['field']['emission'] for s in snaps]),opacity=np.asarray([s['field']['opacity'] for s in snaps]),
        outgoing_left=np.asarray([s['field']['outgoing_left_flux'] for s in snaps]),outgoing_right=np.asarray([s['field']['outgoing_right_flux'] for s in snaps]))
    if probe.enabled:
        optical=[probe.evaluate(s['n']) for s in snaps]
        arrays.update(probe_wavelength_nm=probe.wavelength/10,probe_incident=probe.incident,probe_transmission=np.asarray([r[1] for r in optical]),probe_tau=np.asarray([r[0] for r in optical]),probe_tau_position=np.asarray([r[2] for r in optical]))
        metadata['probe']=probe.metadata;arrays['metadata_json']=np.asarray(json.dumps(metadata,ensure_ascii=False))
    temporary=folder/(FILE_STORE+'.tmp')
    with temporary.open('wb') as stream:np.savez_compressed(stream,**arrays)
    temporary.replace(folder/FILE_STORE)
    statuses='; '.join(f'{label}: '+('пройдено' if quality[key+'_status']=='passed' else 'НАРУШЕНО') for key,label in [('mass','масса'),('charge','заряд'),('kinetic','кинетика всех уровней и электронов'),('energy','энергия')])
    text=f'Результаты Plasma kinetics {VERSION}\n\n'+statuses+'\n\nсохранение.csv — относительные невязки балансов.\nкинетика_уровней.csv — все компоненты, реагенты и уровни, включая электроны.\nреакции_по_времени.csv — полные реагенты, продукты, скорости и интегралы.\nтоп5_по_времени.csv — общие лидеры и заселение/расселение каждого уровня и компонента; доля от полного процесса.\nданные_для_графиков.npz — исходные массивы для окна «Результаты»; вручную открывать не требуется.\n\nГрафики создаются только по выбранным галочкам в окне «Результаты». Все описания реакций и данные об их происхождении доступны во вкладке «Реакции» через «Справочник».\n'
    if not np.any(arrays['electrons_cm3'][0]):text+='\nНачальная концентрация электронов равна нулю. Температура EEDF задаёт форму, а не число частиц. Задайте nₑ, исходные ионы или электронный источник, если они нужны в эксперименте.\n'
    if sources:
        text+='\nИсточники:\n'
        for source in sources:
            mode_label={'ions':'ионы','electrons':'электроны','radiation':'свет','particles':'нейтральные частицы'}[source['mode']]
            text+=f"{source['id']}: {mode_label}, начало {source['start_s']:.6g} с, конец {source['stop_s']:.6g} с"
            if source['mode']=='radiation':text+=f", мощность {source['normal_power_W_cm2']:.6g} Вт/см², фотоны {source['minimum_photon_eV']:.6g}…{source['maximum_photon_eV']:.6g} эВ, средняя {source['mean_photon_eV']:.6g} эВ"
            else:text+=f", n={source.get('density_cm3','старый формат')} см⁻³, u={source.get('velocity_cm_s','старый формат')} см/с"
            text+='\n'
    if m.data_issues:text+='\nОграничения физических данных:\n'+'\n'.join(sorted(m.data_issues))+'\n'
    (folder/'О_результатах.txt').write_text(text,encoding='utf-8');remove_old_generated_files(folder);m.progress.done(str(folder));return quality


class ResultData:
    def __init__(self,folder):
        self.folder=Path(folder)
        with np.load(self.folder/FILE_STORE,allow_pickle=False) as store:self.arrays={key:store[key] for key in store.files}
        self.meta=json.loads(str(self.arrays.pop('metadata_json')));self.states=self.meta['states'];self.reactions=self.meta['reactions']

    def __getitem__(self,key):return self.arrays[key]
