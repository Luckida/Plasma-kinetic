"""One result-backed catalogue for targets, partners and automatic suggestions.

No coefficients are calculated here. Identities and stoichiometry come from
the saved network, and numerical suggestions come from its saved event rates.
"""
from dataclasses import dataclass
import re
import unicodedata
import numpy as np

KIND_LABELS={'excitation':'Возбуждение электроном','deexcitation':'Девозбуждение электроном','ionization':'Ионизация электроном',
    'spontaneous':'Спонтанное излучение','two_photon':'Двухфотонный распад','photoexcitation':'Фотовозбуждение',
    'photoionization':'Фотоионизация','stimulated_emission':'Вынужденное излучение','radiative_recombination':'Радиационная рекомбинация',
    'stimulated_recombination':'Вынужденная рекомбинация','penning':'Пеннинговская ионизация','associative_ionization':'Ассоциативная ионизация',
    'charge_exchange':'Перезарядка','dissociative_recombination':'Диссоциативная рекомбинация','mutual_neutralization':'Взаимная нейтрализация',
    'heavy_chemistry':'Реакция тяжёлых частиц','heavy_quenching':'Столкновительное тушение'}


def readable_issue(issue):
    replacements={': explicitly labelled estimate disabled; chemistry.allow_estimated_rates=true required':': оценочный коэффициент отключён.',': no unique k(T) interval (диапазона) for the occupied distributions; see chemical_domain_audit.csv':': нет применимого опубликованного диапазона k(T) для текущего распределения.',': sigma(E) does not cover occupied relative energies':': σ(E) не покрывает занятые относительные энергии.'}
    for old,new in replacements.items():issue=issue.replace(old,new)
    return issue


def folded(text):
    text=unicodedata.normalize('NFKC',str(text)).casefold().replace('−','-').replace('ν','v')
    return re.sub(r'[^a-zа-яё0-9+*]', '',text)


@dataclass(frozen=True)
class Choice:
    key:str
    label:str
    members:frozenset


@dataclass
class ReactionRequest:
    metric:str
    target:str
    direction:str
    matching_ids:list
    selected_ids:list


class ReactionCatalog:
    def __init__(self,data):
        self.data=data;self.reactions=data.reactions;self.states={s['id']:s for s in data.states};self.by_id={r['id']:i for i,r in enumerate(self.reactions)}
        self.left=[r.get('reactants',{}) for r in self.reactions];self.right=[r.get('products',{}) for r in self.reactions]
        self.participants=[set(a)|set(b) for a,b in zip(self.left,self.right)]
        used=set().union(*self.participants) if self.participants else set()
        self.choices=[Choice('all','Любой реагент',frozenset(used)),Choice('electron','Электроны · e⁻',frozenset({'e'})),Choice('photon','Фотоны · hν',frozenset({'photon'}))]
        self.targets=[Choice('all','Вся смесь · скорости событий',frozenset(used)),Choice('electron','Электроны · e⁻',frozenset({'e'}))]
        for name,indices in data.meta['components'].items():
            members=frozenset(data.states[i]['id'] for i in indices);choice=Choice('component:'+name,'Компонент: '+name,members)
            self.choices.append(choice);self.targets.append(choice)
        elements=sorted({e for state in data.states for e in state.get('nuclei',{state['element']:1})})
        for element in elements:
            members=frozenset(s['id'] for s in data.states if element in s.get('nuclei',{s['element']:1}))
            self.choices.append(Choice('element:'+element,'Элемент: '+element+' · все стадии',members))
        ground={name:min(float(s.get('energy_cm1',0.)) for s in data.states if s['species']==name) for name in data.meta['components']}
        for state in data.states:
            choice=Choice('level:'+state['id'],'Уровень: '+state.get('display_label',state['label']),frozenset({state['id']}))
            self.choices.append(choice);self.targets.append(choice)
        self.choice_map={c.key:c for c in self.choices};self.target_map={c.key:c for c in self.targets}
        self.aliases={};self.search_text=[]
        def alias(text,ident):self.aliases.setdefault(folded(text),set()).add(ident)
        for state in data.states:
            ident=state['id'];element=state['element'];charge=int(state['charge']);caption=state.get('display_label',state['label'])
            for text in (ident,state['label'],caption,state['species']):alias(text,ident)
            if not state.get('molecular'):
                stage=element+('+'*charge);alias(stage,ident)
                if charge==0:
                    n=int(state.get('n',0));mult=int(state.get('multiplicity',1));L=int(state.get('L',0));term='SPDFGHIKLM'[L] if 0<=L<10 else ''
                    alias(f'{element}({n}^{mult}{term})',ident)
                    if float(state.get('energy_cm1',0.))>ground.get(state['species'],0.)+1e-9:alias(element+'*',ident)
                    if element=='He' and n==2 and L==0:
                        alias('метастабильный гелий',ident);alias('метастабильный '+caption,ident)
            for element in state.get('nuclei',{element:1}):alias('элемент '+element,ident)
        for text in ('e','e-','e⁻','электрон','электроны','electron','electrons'):alias(text,'e')
        for text in ('hv','hν','фотон','фотоны','свет','излучение','photon'):alias(text,'photon')
        for r in self.reactions:self.search_text.append(folded(r.get('reaction','')+' '+r['id']+' '+KIND_LABELS.get(r['kind'],r['kind'])+' '+str(r.get('source',''))))
        self.all_peaks=data['reaction_rates'].max(axis=(0,1)) if self.reactions else np.zeros(0)

    def weights(self,target,direction):
        if target=='all':return np.ones(len(self.reactions))
        if direction=='participation':return np.array([bool(side&self.target_map[target].members) for side in self.participants],dtype=float)
        if target=='electron':stoich=self.data['electron_stoichiometry']
        elif target.startswith('level:'):
            j=next(i for i,s in enumerate(self.data.states) if s['id']==target[6:]);stoich=self.data['stoichiometry'][j]
        else:stoich=self.data['stoichiometry'][self.data.meta['components'][target[10:]]].sum(axis=0)
        return abs(stoich) if direction=='both' else np.maximum(stoich if direction=='source' else -stoich,0)

    def context_rates(self,cell=-1,axis='time',time_index=-1,time_min=None,time_max=None,**_):
        rates=self.data['reaction_rates']
        if axis=='space':return rates[time_index%len(rates)]
        rates=rates.mean(axis=1) if cell<0 else rates[:,cell]
        times=self.data['time_s'];mask=np.ones(len(times),dtype=bool)
        if time_min is not None:mask&=times>=time_min
        if time_max is not None:mask&=times<=time_max
        if not mask.any():raise ValueError('В этом интервале нет сохранённых моментов. Измените границы времени.')
        return rates[mask]

    def choices_for_partners(self,first='all',target='all',direction='participation'):
        eligible=self.matching(target=target,direction=direction)
        used=set().union(*(set(self.left[k]) for k in eligible)) if eligible else set()
        compatible=set()
        for k in eligible:
            for name,count in self.left[k].items():
                if first=='all':compatible.add(name);continue
                for a,amount in self.left[k].items():
                    if a in self.choice_map[first].members and (a!=name or amount>=2):compatible.add(name)
        return [c for c in self.choices if c.key=='all' or bool(c.members&(compatible if first!='all' else used))]

    def pair_matches(self,left,a,b):
        if a=='all' and b=='all':return True
        if a=='all':return bool(set(left)&self.choice_map[b].members)
        if b=='all':return bool(set(left)&self.choice_map[a].members)
        return any(x in self.choice_map[a].members and y in self.choice_map[b].members and (x!=y or left[x]>=2) for x in left for y in left)

    def matching(self,target='all',direction='participation',partner_a='all',partner_b='all',kind='all',search='',nonzero=False,context=None):
        weights=self.weights(target,direction);query=folded(search);alias_members=self.aliases.get(query);terms=[folded(term) for term in search.split() if folded(term)]
        peaks=self.context_rates(**(context or {})).max(axis=0) if nonzero and self.reactions else self.all_peaks
        selected=[]
        for k,r in enumerate(self.reactions):
            if weights[k]<=0 or (kind!='all' and r['kind']!=kind):continue
            if not self.pair_matches(self.left[k],partner_a,partner_b):continue
            if nonzero and peaks[k]<=0:continue
            if query:
                if alias_members is not None:
                    if not (self.participants[k]&alias_members):continue
                elif not all(term in self.search_text[k] for term in terms):continue
            selected.append(k)
        return selected

    def ranked(self,indices,target='all',direction='participation',context=None):
        if not indices:return []
        peaks=self.context_rates(**(context or {})).max(axis=0);weights=self.weights(target,direction)
        return sorted(indices,key=lambda k:(-float(peaks[k]*weights[k]),self.reactions[k]['id']))

    def describe(self,k,context=None,peak=None):
        r=self.reactions[k];peak=float(self.context_rates(**(context or {}))[:,k].max()) if peak is None else float(peak)
        left=r.get('left') or r.get('reaction','').split('→')[0];right=r.get('right') or r.get('reaction','').split('→')[-1]
        status='Ненулевая скорость' if peak>0 else 'Скорость равна нулю в выбранной области'
        details=(r['reaction']+'\nРеагенты: '+left+'\nПродукты: '+right+'\nПроцесс: '+KIND_LABELS.get(r['kind'],r['kind'])+
            f'\nМаксимальная скорость событий: {peak:.6g} см⁻³ с⁻¹. '+status+'.\nИсточник данных: '+str(r.get('source') or 'Не указан в сохранённом результате.')+'\nID: '+r['id'])
        from reaction_data import kinetics_summary
        details+='\n'+kinetics_summary(r)
        issues=[readable_issue(issue) for issue in self.data.meta.get('quality',{}).get('data_issues',[]) if issue.startswith(r['id']+':')]
        if issues:details+='\nУстановленные ограничения:\n'+'\n'.join(issues)
        audit=self.data.meta.get('chemical_domain_audit',{}).get(r['id'])
        if audit:details+=f"\nВстречавшаяся температура коэффициента: {audit['encountered_min_K']:.6g}…{audit['encountered_max_K']:.6g} K. Диапазоны данных: {audit['supported_intervals_K']} K."
        return details

    def describe_many(self,indices,context=None):
        peaks=self.context_rates(**(context or {})).max(axis=0)
        return '\n\n'+'\n\n'.join(self.describe(k,context,peaks[k]) for k in indices)

    def absence_message(self):
        text='В сохранённой сети нет реакций для этих условий выбора. Попробуйте «Любой реагент» или «Участие».'
        disabled=[name for name,value in self.data.meta.get('config',{}).get('processes',{}).items() if value is False]
        if disabled:text+=' Выключенные процессы расчёта: '+', '.join(disabled)+'.'
        gaps=self.data.meta.get('quality',{}).get('data_issues',[])
        if gaps:
            translated=[readable_issue(issue) for issue in gaps]
            text+='\n\nОграничения данных этого расчёта:\n'+'\n'.join(translated)
        omissions=self.data.meta.get('chemical_omissions',[])
        if omissions:text+='\n\nКаналы, не включённые в эту сеть:\n'+'\n'.join(str(row.get('id',''))+': '+str(row.get('reason','')) for row in omissions)
        return text
