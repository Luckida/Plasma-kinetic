"""Join actual atomic packs into one global energy/charge/state catalogue.

No helium substitution and no conversion of Maxwellian Upsilon into sigma.
All stored levels are candidates; level_selection chooses the working basis.
"""
from copy import deepcopy
import csv
import hashlib
import json
from pathlib import Path
import re
import numpy as np
from atomic_core import Species,selected_species,SYMBOLS
from atomic_database import prepare_species,HC_EV_CM,load_catalogue,pack_directory

ROOT=Path(__file__).resolve().parent.parent


def read_csv(path):
    with Path(path).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter=';'))


class MixtureData:
    def __init__(self,config,root=ROOT,progress=None):
        self.root=Path(root);options=deepcopy(config.get('atomic_data',{}));requested=selected_species(config)
        local_tables=options.pop('local_tables',False)
        if not isinstance(local_tables,bool):raise ValueError('atomic_data.local_tables должен быть bool')
        registry=load_catalogue(self.root,options)
        self.descriptor_hashes=registry.get('descriptor_files_sha256',{})
        grounds=registry.get('ground_stages') or json.loads((self.root/config.get('ion_stages_file','ion_stages.json')).read_text(encoding='utf-8-sig'))
        self.elements=sorted({s.element for s in requested});imported=[]
        for element in self.elements:
            if element not in grounds:raise ValueError('Нужны ground-stage данные '+element+' в ion_stages_file')
            maximum=max(s.charge for s in requested if s.element==element)
            # Known lower stages are needed for recombination of an entered ion.
            for q in range(min(maximum+2,len(grounds[element]))):
                s=Species(element,q)
                if s.key in registry['species']:imported.append(s)
                elif q==0:raise ValueError('Нет проверенного атомного импортёра '+s.label)
        self.imported=imported
        request=deepcopy(config);request['atomic_data']=options
        request['atomic_data'].update(species=[s.label for s in imported],nmax=options.get('nmax',9),
            prepared_directory=options.get('prepared_directory','atomic_runtime_bank'),discover=False)
        if not local_tables:prepare_species(request,self.root,progress)
        folder=self.root/request['atomic_data']['prepared_directory'];cap=request['atomic_data']['nmax']
        packs={s.key:pack_directory(folder,s,cap,options) for s in imported}
        self.states=[];self.lines=[];self.cross={};self.cross_source={};self.elastic={};self.gaps=[];self.manifests={};self.parent_index={}
        stage_rows={};stage_edges={};masses={}
        for s in imported:
            pack=packs[s.key];m=json.loads((pack/'manifest.json').read_text(encoding='utf-8-sig'));self.manifests[s.key]=m
            if local_tables:
                m['local_tables']=True;m['local_modified_files']={name:hashlib.sha256((pack/name).read_bytes()).hexdigest() if (pack/name).is_file() else 'missing' for name,digest in m.get('files_sha256',{}).items()
                    if not (pack/name).is_file() or hashlib.sha256((pack/name).read_bytes()).hexdigest()!=digest}
            if not (pack/'levels.csv').exists():raise ValueError('Не подготовлены уровни '+s.label)
            rows=read_csv(pack/'levels.csv')
            if not rows:raise ValueError('Нет уровней '+s.label+'; '+str(m['gaps']))
            for row in rows:
                for key in ('n','multiplicity','L','g','source_index'):row[key]=int(row[key])
                for key in ('energy_cm1','mass_u'):row[key]=float(row[key])
                if not np.all(np.isfinite([row['energy_cm1'],row['mass_u'],row['g']])) or row['energy_cm1']<0 or row['mass_u']<=0 or row['g']<=0:raise ValueError('Неверные уровни/массы: '+s.label)
                row['element']=s.element;row['charge']=s.charge;row['nuclei']={s.element:1}
                row['is_ground']=row['energy_cm1']==min(float(r['energy_cm1']) for r in rows)
                row['excitation_eV']=row['energy_cm1']*HC_EV_CM
            if s.element in masses and abs(rows[0]['mass_u']/masses[s.element]-1)>1e-9:raise ValueError('Несогласованные массы стадий '+s.element+'; нужен один изотопный/ядерный массовый стандарт')
            stage_rows[s.key]=rows;masses[s.element]=rows[0]['mass_u']
            stage_edges[s.key]=m['adf04']['parents'][0]['threshold_cm1']
        self.stage_ground={}
        for element in self.elements:
            offset=0.
            for ground in grounds[element]:
                q=ground['charge'];s=Species(element,q);rows=stage_rows.get(s.key)
                if rows:
                    for row in rows:row['energy_cm1']+=offset
                    self.states.extend(rows);ident=min(rows,key=lambda r:r['excitation_eV'])['id']
                else:
                    ident=s.key+':ground';row=dict(ground,id=ident,element=element,species=s.label,
                        label=s.label+' '+ground['configuration'],nuclei={element:1},mass_u=masses[element],
                        energy_cm1=offset,excitation_eV=0.,is_ground=True,source=ground.get('source','Element descriptor '+element+'/atom.json; electronic ground'))
                    self.states.append(row)
                self.stage_ground[s.key]=ident
                if 'ionization_eV' in ground:offset+=stage_edges.get(s.key,ground['ionization_eV']/HC_EV_CM)
        from mixture_chemistry import append_components
        append_components(self,config,self.root)
        self.index={s['id']:j for j,s in enumerate(self.states)}
        if len(self.index)!=len(self.states):raise ValueError('Повтор глобального уровня')
        self.by={s['id']:s for s in self.states}
        for j,state in enumerate(self.states):
            key=f'{state["element"]}_{state["charge"]+1}'
            if key in self.stage_ground:self.parent_index[j]=self.index[self.stage_ground[key]]
        for s in imported:
            pack=packs[s.key]
            for row in read_csv(pack/'lines.csv'):
                row['A_s']=float(row['A_s']);row['wavelength_nm']=float(row['wavelength_nm']);row['spectral_kind']='line';row['pump']=True
                self.lines.append(row)
            grouped={}
            for row in read_csv(pack/'electron_cross_sections.csv'):
                grouped.setdefault(row['id'],[]).append(row)
            for ident,rows in grouped.items():
                row=rows[0];a=row['initial'];kind=row['kind'];b=row['final']
                table=np.array([(float(r['energy_eV']),float(r['sigma_cm2'])) for r in rows])
                if len(table)<2 or not np.all(np.isfinite(table)) or np.any(table<0) or np.any(np.diff(table[:,0])<=0):raise ValueError('Неверная таблица sigma(E): '+ident)
                if kind in ('electron_elastic','electron_momentum_transfer'):
                    old=self.elastic.get(a)
                    if old is None or kind=='electron_momentum_transfer':self.elastic[a]=(table,row['source'],kind)
                    continue
                if kind=='electron_ionization':b=self.states[self.parent_index[self.index[a]]]['id'];kind='ionization'
                else:kind=kind.replace('electron_','')
                gap=(self.by[b]['energy_cm1']-self.by[a]['energy_cm1'])*HC_EV_CM
                # Only a supplied zero threshold point identifies the threshold.
                if kind=='excitation' and gap>0 and table[0,1]==0:table[:,0]+=gap-table[0,0]
                self.cross[(a,b,kind)]=table;self.cross_source[(a,b,kind)]=row['source']
        # Known inverse excitation via detailed balance, not inversion of Upsilon.
        for (a,b,kind),tab in list(self.cross.items()):
            if kind!='excitation':continue
            gap=(self.by[b]['energy_cm1']-self.by[a]['energy_cm1'])*HC_EV_CM
            if gap<=0 or (b,a,'deexcitation') in self.cross:continue
            ee=tab[:,0]-gap;mask=ee>0
            if mask.sum()<2:continue
            inverse=np.column_stack((ee[mask],self.by[a]['g']/self.by[b]['g']*tab[mask,0]/ee[mask]*tab[mask,1]))
            self.cross[(b,a,'deexcitation')]=inverse;self.cross_source[(b,a,'deexcitation')]=self.cross_source[(a,b,kind)]+'; microscopic detailed balance'
        self._weak_decays()

    def _weak_decays(self):
        lookup={(s['element'],s['charge'],s.get('n'),s.get('multiplicity'),s.get('L')):s['id'] for s in self.states}
        path=ROOT/'data/radiation/supplementary_transitions.csv'
        for row in read_csv(path):
            def translate(name):
                m=re.fullmatch(r'HeI_(\d+)_(\d+)([SPDF])',name)
                return lookup.get(('He',0,int(m[1]),int(m[2]),'SPDF'.index(m[3]))) if m else None
            lower,upper=translate(row['lower']),translate(row['upper'])
            if lower and upper:
                self.lines.append(dict(row,lower=lower,upper=upper,A_s=float(row['A_s']),pump=False,
                    wavelength_nm=1e7/(self.by[upper]['energy_cm1']-self.by[lower]['energy_cm1'])))
        for element in self.elements:
            Z=SYMBOLS.index(element)+1;q=Z-1
            low=lookup.get((element,q,1,2,0));up=lookup.get((element,q,2,2,0))
            if low and up:self.lines.append(dict(id=f'2PH_{element}_{q}',lower=low,upper=up,A_s=8.2293*Z**6,
                wavelength_nm=1e7/(self.by[up]['energy_cm1']-self.by[low]['energy_cm1']),spectral_kind='two_photon',pump=False,
                source='Chluba & Sunyaev 2008, A&A 480, 629; NR 2s-1s Z^6 scaling'))

    def resolve(self,name):
        if name in self.index:return self.index[name]
        if name=='HeII':name='He II'
        old=re.fullmatch(r'HeI_(\d+)_(\d+)([SPDF])',name)
        if old:
            matches=[j for j,s in enumerate(self.states) if s['element']=='He' and s['charge']==0 and s.get('n')==int(old[1]) and s.get('multiplicity')==int(old[2]) and s.get('L')=='SPDF'.index(old[3])]
            if matches:return matches[0]
        species=Species.parse(name)
        if species.key not in self.stage_ground:raise ValueError('Не включён вид '+species.label)
        return self.index[self.stage_ground[species.key]]
