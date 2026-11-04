"""Data-driven reactions between particles, including Penning and charge exchange.

Stoichiometry, rate fits and energy cross sections are separate data files.
No unknown coefficient is silently replaced by zero or a gas-kinetic guess.
"""
from collections import Counter
from itertools import combinations_with_replacement
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from numba.typed import List as TypedList
from atomic_database import HC_EV_CM,write_csv,write_json
from atomic_core import Network
from particle_kinetics import ME,EV,KB,GridOverflow
import chemistry_kernels as CK
from kinetic_kernels import kernel_threads


AMU=1.66053906892e-24


def records(path):
    with Path(path).open(encoding='utf-8-sig',newline='') as stream:
        return list(csv.DictReader(stream,delimiter=';'))


def append_components(bank,config,root):
    """Add declared ground-state carriers; never pretend to resolve v,J levels."""
    cfg=config.get('chemistry',{})
    bank.component_files={}
    if not cfg.get('enabled',config.get('schema')=='MixtureKineticSlab-2.1'):return
    folder=Path(root)/cfg.get('components_directory','data/components')
    for file in sorted(folder.glob('*.json')):
        document=json.loads(file.read_text(encoding='utf-8'))
        if document.get('schema')!='ChemicalComponents-1':raise ValueError('Неверная схема '+str(file))
        for item in document.get('components',[]):
            if not set(item['nuclei'])<=set(bank.elements):continue
            if item['id'] in {s['id'] for s in bank.states}:raise ValueError('Повтор компоненты '+item['id'])
            anchors=item['energy_anchors'];energy=float(item.get('energy_offset_eV',0.))
            for stage,weight in anchors.items():
                if stage not in bank.stage_ground:raise ValueError('Нет энергетического якоря '+stage)
                energy+=weight*next(s['energy_cm1'] for s in bank.states if s['id']==bank.stage_ground[stage])*HC_EV_CM
            mass=sum(next(s['mass_u'] for s in bank.states if s['element']==element)*number for element,number in item['nuclei'].items())
            if not np.isfinite(energy) or not np.isfinite(mass) or mass<=0:raise ValueError('Неверная энергия/масса компоненты '+item['id'])
            row=dict(item,element=item['id'],species=item['id'],mass_u=mass,energy_cm1=energy/HC_EV_CM,
                     excitation_eV=0.,is_ground=True,molecular=True,n=0,multiplicity=1,L=0,g=item.get('g',1),configuration='effective ground carrier')
            bank.states.append(row)
            bank.component_files[file.name]=hashlib.sha256(file.read_bytes()).hexdigest()


def electron_release_bound(bank,config,root):
    """Include known chemi-ionization in the automatic electron grid budget."""
    cfg=config.get('chemistry',{})
    if not cfg.get('enabled',config.get('schema')=='MixtureKineticSlab-2.1'):return 0.
    bound=0.
    for file in sorted((Path(root)/cfg.get('directory','data/reactions')).glob('*/reactions.json')):
        document=json.loads(file.read_text(encoding='utf-8-sig'))
        if document.get('reviewed') is not True or not set(document.get('elements',[]))<=set(bank.elements):continue
        for row in document.get('reactions',[]):
            if not row['products'].get('e') or not config.get('processes',{}).get(row['kind'],True) or row['id'] in cfg.get('disabled_ids',[]) or row['id'] in config.get('processes',{}).get('disabled_ids',[]):continue
            def energy(name):
                if name in ('e','e-'):return 0.
                selector=document.get('selectors',{}).get(name)
                if selector:
                    found=[s for s in bank.states if all(s.get(k)==v for k,v in selector.items())]
                    if len(found)!=1:raise ValueError('Нет уровня')
                    return found[0]['energy_cm1']*HC_EV_CM
                return bank.states[bank.resolve(name)]['energy_cm1']*HC_EV_CM
            try:q=sum(energy(n)*v for n,v in row['reactants'].items())-sum(energy(n)*v for n,v in row['products'].items())
            except (ValueError,KeyError):continue
            bound=max(bound,q)
    return bound


def rate_fit(rows,temperature,gas_temperature=None):
    """Evaluate only inside a published interval. k(T), not sigma(E)."""
    if not np.isfinite(temperature) or temperature<=0:raise ValueError('T реакции должна быть >0')
    if {r['model'] for r in rows}=={'tabulated_loglog'}:
        points=sorted((float(r['reference_T_K']),float(r['alpha_cm3_s'])) for r in rows)
        nodes=np.asarray(points)
        if len(points)<2 or np.any(np.diff(nodes[:,0])<=0) or np.any(nodes<=0):raise ValueError('Неверная таблица k(T)')
        if not nodes[0,0]<=temperature<=nodes[-1,0]:raise ValueError('Табличный k(T) вне опубликованного интервала')
        return float(np.exp(np.interp(np.log(temperature),np.log(nodes[:,0]),np.log(nodes[:,1]))))
    valid=[r for r in rows if float(r['Tmin_K'])<=temperature<=float(r['Tmax_K'])]
    if len(valid)!=1:raise ValueError('Нет единственного температурного диапазона k(T) при T='+str(temperature))
    row=valid[0];alpha=float(row['alpha_cm3_s']);beta=float(row['beta']);gamma=float(row['gamma_K'])
    if row['model']=='arrhenius':value=alpha*(temperature/float(row.get('reference_T_K') or 300.))**beta*np.exp(-gamma/temperature)
    elif row['model']=='power_exponential_temperature':
        if gamma<=0:raise ValueError('Для exp(-T/gamma) gamma_K>0')
        value=alpha*(temperature/float(row.get('reference_T_K') or 1.))**beta*np.exp(-temperature/gamma)
    elif row['model']=='two_temperature_power':
        if gas_temperature is None or not float(row['Tgmin_K'])<=gas_temperature<=float(row['Tgmax_K']):raise ValueError('Вне заданного газового диапазона двухтемпературного коэффициента')
        value=alpha*(temperature/gas_temperature)**beta
    elif row['model'] in ('roberge_penning_branch','roberge_associative_branch'):
        # Explicit interpolation model through the two branching fractions
        # quoted by Roberge & Dalgarno: 0.3 at 300 K, 0.02 at 4500 K.
        branch=90./temperature
        value=alpha*(branch if row['model']=='roberge_associative_branch' else 1-branch)
    else:raise ValueError('Неизвестный вид k(T): '+row['model'])
    if not np.isfinite(value) or value<0:raise ValueError('Неверный коэффициент реакции')
    return float(value)


class Chemistry:
    def __init__(self,model):
        self.model=model;self.cfg=model.cfg.get('chemistry',{});self.reactions=[];self.gaps=set();self.notes=set();self.files={};self.geometries={};self.participants=set();self.omissions=[];self.domain_audit={}
        allowed={'enabled','directory','components_directory','rate_policy','disabled_ids','relative_azimuth_points','product_angular_points','product_azimuth_points','coverage_policy','datasets','allow_estimated_rates','coefficient_preference'}
        if set(self.cfg)-allowed:raise ValueError('Неизвестные chemistry: '+str(set(self.cfg)-allowed))
        if not isinstance(self.cfg.get('enabled',True),bool):raise ValueError('chemistry.enabled должен быть bool')
        if not isinstance(self.cfg.get('allow_estimated_rates',False),bool):raise ValueError('chemistry.allow_estimated_rates должен быть bool')
        choices=self.cfg.get('datasets',{})
        if not isinstance(choices,dict) or any(not isinstance(k,str) or not isinstance(v,str) for k,v in choices.items()):raise ValueError('chemistry.datasets: объект группа -> имя набора')
        if not isinstance(self.cfg.get('disabled_ids',[]),list) or any(not isinstance(v,str) for v in self.cfg.get('disabled_ids',[])):raise ValueError('chemistry.disabled_ids: список строк')
        if not self.cfg.get('enabled',model.cfg.get('schema')=='MixtureKineticSlab-2.1'):return
        self.policy=self.cfg.get('rate_policy','cross_sections_only')
        if self.policy not in ('cross_sections_only','maxwellian_moments'):raise ValueError('chemistry.rate_policy: cross_sections_only или maxwellian_moments')
        if self.cfg.get('coefficient_preference','cross_sections') not in ('cross_sections','rates'):raise ValueError('coefficient_preference: cross_sections или rates')
        if self.cfg.get('coefficient_preference')=='rates' and self.policy!='maxwellian_moments':raise ValueError('Тепловые rates требуют maxwellian_moments')
        if self.cfg.get('coverage_policy','report') not in ('report','error'):raise ValueError('chemistry.coverage_policy: report или error')
        naz=self.cfg.get('relative_azimuth_points',2);nmu=self.cfg.get('product_angular_points',2);nphi=self.cfg.get('product_azimuth_points',2)
        if any(isinstance(n,bool) or not isinstance(n,int) or n<2 or n>32 for n in [naz,nmu,nphi]):raise ValueError('Квадратуры chemistry: целые 2..32')
        self.azimuths=2*np.pi*(np.arange(naz)+.5)/naz
        self.smu,self.sw=np.polynomial.legendre.leggauss(nmu);self.sw*=.5
        self.sphi=2*np.pi*(np.arange(nphi)+.5)/nphi
        folder=model.root/self.cfg.get('directory','data/reactions')
        if not folder.is_dir():raise ValueError('Нет папки химических данных '+str(folder))
        signatures=set();identifiers=set();registered_datasets={}
        for file in sorted(folder.glob('*/reactions.json')):
            document=json.loads(file.read_text(encoding='utf-8'))
            if document.get('schema')!='ParticleReactions-1' or document.get('reviewed') is not True:raise ValueError('Нужен проверенный ParticleReactions-1: '+str(file))
            for raw in document.get('reactions',[]):
                if raw.get('alternative_group'):registered_datasets.setdefault(raw['alternative_group'],set()).add(raw['alternative'])
            if not set(document.get('elements',[]))<=set(model.data.elements):continue
            if (file.parent/'unsupported.csv').is_file():self.omissions.extend(records(file.parent/'unsupported.csv'))
            for part in file.parent.iterdir():
                if part.is_file():self.files[str(part.relative_to(model.root)) if part.is_relative_to(model.root) else str(part)]=hashlib.sha256(part.read_bytes()).hexdigest()
            fits={};cross={}
            if (file.parent/'rate_coefficients.csv').is_file():
                for row in records(file.parent/'rate_coefficients.csv'):
                    fits.setdefault(row['id'],[]).append(row)
                    for key in ('alpha_cm3_s','beta','gamma_K','Tmin_K','Tmax_K'):
                        if not np.isfinite(float(row[key])):raise ValueError('Неконечный k(T): '+row['id'])
                    if float(row['alpha_cm3_s'])<0 or not 0<float(row['Tmin_K'])<=float(row['Tmax_K']):raise ValueError('Неверный диапазон k(T)')
            if (file.parent/'cross_sections.csv').is_file():
                for row in records(file.parent/'cross_sections.csv'):cross.setdefault(row['id'],[]).append(row)
            for raw in document.get('reactions',[]):
                if raw.get('alternative_group'):
                    selected=self.cfg.get('datasets',{}).get(raw['alternative_group'],document.get('default_datasets',{}).get(raw['alternative_group']))
                    if selected!=raw.get('alternative'):continue
                if raw.get('estimated') and not self.cfg.get('allow_estimated_rates',False):
                    self.gaps.add(raw['id']+': explicitly labelled estimate disabled; chemistry.allow_estimated_rates=true required');continue
                if any(isinstance(v,bool) or not isinstance(v,int) or v<=0 for side in ('reactants','products') for v in raw[side].values()):raise ValueError('Нужны целые положительные стехиометрические числа '+raw.get('id',''))
                ident=raw['id']
                if not model.enabled(raw['kind']):continue
                if ident in self.cfg.get('disabled_ids',[]) or ident in model.process.get('disabled_ids',[]):continue
                if ident in identifiers:raise ValueError('Повтор id реакции '+ident)
                identifiers.add(ident)
                try:
                    left=Counter();right=Counter()
                    for name,value in raw['reactants'].items():left[self.resolve(name,document.get('selectors',{}))]+=value
                    for name,value in raw['products'].items():right[self.resolve(name,document.get('selectors',{}))]+=value
                    left=dict(left);right=dict(right)
                except ValueError as exc:
                    self.gaps.add(ident+': '+str(exc));continue
                if any(isinstance(v,bool) or not isinstance(v,int) or v<=0 for v in [*left.values(),*right.values()]):raise ValueError('Неверная стехиометрия '+ident)
                if sum(left.values())!=2 or sum(right.values()) not in (2,3):raise ValueError('Реакционный оператор поддерживает 2 -> 2 или 2 -> 3: '+ident)
                if all(n=='e' for n in left) or all(n=='e' for n in right):raise ValueError('Нужны тяжёлые продукты и реагенты')
                signature=(tuple(sorted(left.items())),tuple(sorted(right.items())),bool(raw.get('radiative')))
                if signature in signatures:raise ValueError('Повтор одной стехиометрии: '+ident+'; выберите один источник')
                signatures.add(signature)
                reaction=dict(raw,reactants=left,products=right)
                network_states=[dict(s,charge=int(s['charge'])) for s in model.states]+[dict(id='e',label='e⁻',nuclei={},charge=-1)]
                Network(network_states,[reaction])
                if raw.get('model_note'):self.notes.add(ident+': '+raw['model_note'])
                if raw.get('cross_section_id') and not (raw.get('rate_id') and self.cfg.get('coefficient_preference')=='rates'):
                    points=cross.get(raw['cross_section_id'],[])
                    if not points:raise ValueError('Нет sigma(E) '+ident)
                    if {r.get('energy_frame') for r in points}!={'centre_of_mass'} or {r.get('units') for r in points}!={'cm2'}:raise ValueError('sigma(E): нужны centre_of_mass и cm2')
                    tab=np.array([(float(r['energy_eV']),float(r['sigma_cm2'])) for r in points])
                    if len(tab)<2 or np.any(np.diff(tab[:,0])<=0) or np.any(tab<0) or not np.all(np.isfinite(tab)):raise ValueError('Неверная sigma(E) '+ident)
                    reaction['cross']=tab
                else:
                    if raw.get('rate_id') not in fits:raise ValueError('Нет k(T) '+ident)
                    reaction['fits']=fits[raw['rate_id']]
                    if self.policy=='cross_sections_only':
                        self.gaps.add(ident+': только k(T); нужен явно выбранный maxwellian_moments или sigma(E)');continue
                    self.notes.add('k(T) at relative moment temperature, not an EEDF/velocity-resolved cross section: '+ident)
                    if raw.get('model_note'):self.notes.add(ident+': '+raw['model_note'])
                reaction['a'],reaction['b']=[name for name,count in left.items() for _ in range(count)]
                reaction['outputs']=[name for name,count in right.items() if name!='e' for _ in range(count)]+['e']*right.get('e',0)
                if right.get('e',0)>1:raise ValueError('Не более одного нового электрона в chemistry; используйте атомный impact-оператор')
                reaction['q']=sum(self.internal(n)*v for n,v in left.items())-sum(self.internal(n)*v for n,v in right.items())
                reaction['kinematic_q']=reaction['q']
                if raw.get('radiative'):
                    if reaction['q']<=0 or len(reaction['outputs'])!=2 or 'e' in left or 'e' in right:raise ValueError('Radiative CT: нужны экзотермические 2 heavy -> 2 heavy')
                    reaction['photon_eV']=reaction['q'];reaction['kinematic_q']=0.
                    if not model.energy[0]<=reaction['photon_eV']<=model.energy[-1]:raise ValueError('Фотонная сетка не покрывает химическое излучение '+ident)
                    self.notes.add(ident+': photon carries the internal energy defect; relative kinetic energy conserved, photon recoil/continuum shape unresolved')
                if raw.get('thermal_inverse'):
                    if 'cross' in reaction or reaction['q']>=0:raise ValueError('thermal_inverse: нужен эндотермический табличный k(T)')
                    def statistical(side):
                        return np.prod([model.states[model.index[n]]['g']**v for n,v in side.items() if n!='e'])
                    reaction['inverse_weight']=float(statistical(right)/statistical(left))
                if raw.get('kinematics')=='resonant_velocity_exchange':
                    names=[reaction['a'],reaction['b'],*reaction['outputs']]
                    if len(reaction['outputs'])!=2 or 'e' in names or abs(reaction['q'])>1e-9 or max(self.grid(n).mass for n in names)/min(self.grid(n).mass for n in names)-1>1e-12:
                        raise ValueError('Резонансная перезарядка требует одинаковых ядерных масс и Q=0')
                convention=raw.get('coefficient_convention')
                if convention not in ('event_rate','unordered_cross_section'):raise ValueError('Нужно coefficient_convention для '+ident)
                if 'cross' in reaction and convention!='unordered_cross_section':raise ValueError('Для sigma нужна unordered_cross_section')
                reaction['symmetry']=.5 if convention=='unordered_cross_section' and reaction['a']==reaction['b'] else 1.
                heavy_left=[model.index[n] for n in left if n!='e'];heavy_right=[model.index[n] for n in right if n!='e']
                index=model.add(ident,heavy_left[0],heavy_right[0],raw['kind'],raw['source'],reactants=left,products=right,chemistry=True,radiative=bool(raw.get('radiative')))
                if index is not None:reaction['channel']=index;self.reactions.append(reaction)
        self.participants={model.index[n] for r in self.reactions for side in ('reactants','products') for n in r[side] if n!='e'}
        for group,value in choices.items():
            if group not in registered_datasets or value not in registered_datasets[group]:raise ValueError('Неизвестный набор chemistry.datasets: '+group+'='+value)
        self.notes.add('Isotropic outgoing Jacobi angles; 3-product Penning electron gets min(Q, available) relative energy; no measured differential electron spectrum.')

    def resolve(self,name,selectors):
        if name in ('e','e-'):return 'e'
        if name in selectors:
            selector=selectors[name];found=[s['id'] for s in self.model.states if all(s.get(k)==v for k,v in selector.items())]
            if len(found)!=1:raise ValueError('Нет единственного уровня '+name)
            return found[0]
        return self.model.states[self.model.data.resolve(name)]['id']

    def grid(self,name):return self.model.ge if name=='e' else self.model.gh[self.model.states[self.model.index[name]]['element']]
    def population(self,name,z):return self.model.F[z] if name=='e' else self.model.P[z,self.model.index[name]]
    def internal(self,name):return 0. if name=='e' else float(self.model.internal[self.model.index[name]])

    def kernel(self,r,z):
        a,b=r['a'],r['b'];ga,gb=self.grid(a),self.grid(b);fa,fb=self.population(a,z),self.population(b,z)
        if not fa.sum() or not fb.sum():return None
        key=(id(ga),id(gb))
        if key not in self.geometries:self.geometries[key]=CK.geometry(ga.vz.ravel(),ga.vperp.ravel(),gb.vz.ravel(),gb.vperp.ravel(),self.azimuths,ga.mass,gb.mass)
        relative=self.geometries[key];ma,mb=ga.moments(fa),gb.moments(fb);reduced=ga.mass*gb.mass/(ga.mass+gb.mass)
        T=reduced*(ma['T_moment_eV']*EV/ga.mass+mb['T_moment_eV']*EV/gb.mass)/KB+reduced*(ma['uz_cm_s']-mb['uz_cm_s'])**2/(3*KB)
        weights=np.broadcast_to(fa.ravel()[:,None,None]*fb.ravel()[None,:,None],relative.shape)/len(self.azimuths)
        # Exact kinematic feasibility includes the small free-electron mass.
        mout=sum(self.grid(n).mass for n in r['outputs'])
        za,pa=ga.vz.ravel()[:,None,None],ga.vperp.ravel()[:,None,None]
        zb,pb=gb.vz.ravel()[None,:,None],gb.vperp.ravel()[None,:,None]
        p2=(ga.mass*pa)**2+(gb.mass*pb)**2+2*ga.mass*gb.mass*pa*pb*np.cos(self.azimuths)[None,None,:]+(ga.mass*za+gb.mass*zb)**2
        kinetic=.5*(ga.mass*(za*za+pa*pa)+gb.mass*(zb*zb+pb*pb))/EV
        allowed=kinetic+r['kinematic_q']-.5*p2/mout/EV>=-1e-10
        if 'cross' in r:
            table=r['cross'];threshold=max(-r['q'],0.)
            outside=((relative<table[0,0])|(relative>table[-1,0]))&(relative>=threshold)
            fraction=float((weights*outside).sum()/max(weights.sum(),1e-300))
            tolerance=self.model.selection.tolerance if hasattr(self.model,'selection') else 1e-6
            if fraction>tolerance:
                self.gaps.add(r['id']+': sigma(E) does not cover occupied relative energies')
            value=np.interp(relative,table[:,0],table[:,1],left=0.,right=0.)*np.sqrt(2*relative*EV/reduced)*allowed
        else:
            gasT=(mb if a=='e' else ma)['T_moment_eV']*EV/KB
            fitT=(ma if a=='e' else mb)['T_moment_eV']*EV/KB if 'e' in (a,b) and any(row['model']=='two_temperature_power' for row in r['fits']) else T
            audit=self.domain_audit.setdefault(r['id'],dict(id=r['id'],encountered_min_K=fitT,encountered_max_K=fitT,maximum_pair_density_cm6=0.,out_of_domain=False,
                supported_intervals_K=' | '.join(str(row['Tmin_K'])+'..'+str(row['Tmax_K']) for row in r['fits'])))
            audit['encountered_min_K']=min(audit['encountered_min_K'],fitT);audit['encountered_max_K']=max(audit['encountered_max_K'],fitT)
            audit['maximum_pair_density_cm6']=max(audit['maximum_pair_density_cm6'],float(fa.sum()*fb.sum()))
            try:k=rate_fit(r['fits'],fitT,gasT)
            except ValueError as exc:
                audit['out_of_domain']=True
                self.gaps.add(r['id']+': no unique k(T) interval (диапазона) for the occupied distributions; see chemical_domain_audit.csv');return None
            if r.get('thermal_inverse'):k*=r['inverse_weight']*np.exp(r['q']*EV/(KB*fitT))
            supported=float((weights*allowed).sum())
            if supported==0:
                self.gaps.add(r['id']+': no energetically accessible pair on the phase grid');return None
            value=allowed.astype(float)*(k*float(weights.sum())/supported)
        return value

    def event_rates(self,z):
        result={}
        for r in self.reactions:
            kernel=self.kernel(r,z)
            if kernel is not None:
                fa,fb=self.population(r['a'],z),self.population(r['b'],z)
                result[r['channel']]=float(np.einsum('a,b,abc->',fa.ravel(),fb.ravel(),kernel)/len(self.azimuths)*r['symmetry'])
        return result

    def evolve(self,z,dt):
        model=self.model;lossP=np.zeros_like(model.P[z]);gainP=np.zeros_like(lossP);lossF=np.zeros_like(model.F[z]);gainF=np.zeros_like(lossF);events=np.zeros(model.nr)
        for r in self.reactions:
            kernel=self.kernel(r,z)
            if kernel is None:continue
            a,b=r['a'],r['b'];ga,gb=self.grid(a),self.grid(b);grids=[self.grid(n) for n in r['outputs']]
            if r.get('kinematics')=='resonant_velocity_exchange':
                fa,fb=self.population(a,z).ravel(),self.population(b,z).ravel()
                pair_events=fa[:,None,None]*fb[None,:,None]*kernel*dt*r['symmetry']/len(self.azimuths)
                la=pair_events.sum(axis=(1,2));lb=pair_events.sum(axis=(0,2));count=la.sum()
                lossP[model.index[a]]+=la.reshape(ga.shape);lossP[model.index[b]]+=lb.reshape(gb.shape)
                gainP[model.index[r['outputs'][0]]]+=la.reshape(ga.shape);gainP[model.index[r['outputs'][1]]]+=lb.reshape(gb.shape)
                events[r['channel']]=count
                continue
            with kernel_threads(model.workers):
                la,lb,products,count,forbidden,overflow=CK.products(self.population(a,z).ravel(),self.population(b,z).ravel(),ga.mass,gb.mass,
                    ga.vz.ravel(),ga.vperp.ravel(),gb.vz.ravel(),gb.vperp.ravel(),self.azimuths,kernel,dt,r['symmetry'],r['kinematic_q'],np.array([g.mass for g in grids]),
                    TypedList([g.e for g in grids]),TypedList([g.mu for g in grids]),self.smu,self.sw,self.sphi,max(r['q'],0.) if r.get('kind')=='penning' else -1.)
            if forbidden>max(count*1e-12,1e-100):raise ValueError('Энергетически запрещённые события '+r['id'])
            if overflow>max(count*1e-12,1e-100):raise GridOverflow('chemistry '+r['id']+': расширьте энергетическую сетку продуктов')
            events[r['channel']]=count
            for name,values,g in ((a,la,ga),(b,lb,gb)):
                if name=='e':lossF+=values.reshape(g.shape)
                else:lossP[model.index[name]]+=values.reshape(g.shape)
            for j,(name,g) in enumerate(zip(r['outputs'],grids)):
                values=products[j,:len(g.e)*len(g.mu)].reshape(g.shape)
                if name=='e':gainF+=values
                else:gainP[model.index[name]]+=values
        if np.any(lossP>model.loss_limit*model.P[z]) or np.any(lossF>model.loss_limit*model.F[z]):
            return None
        if self.gaps and self.cfg.get('coverage_policy','report')=='error':raise ValueError('Недостаточные данные chemistry: '+sorted(self.gaps)[0])
        momentum=float(((gainF-lossF)*model.ge.vz).sum())*model.ge.mass
        for j,state in enumerate(model.states):
            grid=model.gh[state['element']]
            momentum+=float(((gainP[j]-lossP[j])*grid.vz).sum())*grid.mass
        model.chemical_momentum_error+=momentum*model.dz
        model.P[z]+=gainP-lossP;model.F[z]+=gainF-lossF
        model.photon_work-=sum(events[r['channel']]*r.get('photon_eV',0.) for r in self.reactions)*EV*model.dz
        model.data_issues.update(self.gaps);model.model_notes.update(self.notes)
        return events

    def photon_emission(self,z):
        """Photon-count/energy-conserving two-node internal-defect spectrum."""
        m=self.model;emission=np.zeros_like(m.energy)
        for r in self.reactions:
            if not r.get('photon_eV'):continue
            kernel=self.kernel(r,z)
            if kernel is None:continue
            fa,fb=self.population(r['a'],z),self.population(r['b'],z)
            rate=float(np.einsum('a,b,abc->',fa.ravel(),fb.ravel(),kernel)/len(self.azimuths)*r['symmetry'])
            energy=r['photon_eV'];i=min(max(int(np.searchsorted(m.energy,energy))-1,0),len(m.energy)-2)
            frac=(energy-m.energy[i])/(m.energy[i+1]-m.energy[i])
            for k,weight in [(i,1-frac),(i+1,frac)]:emission[k]+=rate*weight*m.energy[k]*EV/m.w[k]
        return emission

    def export(self,output):
        rows=[]
        for r in self.reactions:
            rows.append(dict(id=r['id'],kind=r['kind'],Q_eV=r['q'],coefficient='sigma(E)' if 'cross' in r else 'k(T); relative Maxwellian moment closure',
                convention=r['coefficient_convention'],radiative=bool(r.get('radiative')),photon_energy_eV=r.get('photon_eV',0.),estimated=bool(r.get('estimated')),source=r['source']))
        write_csv(Path(output)/'chemical_reactions.csv',rows,['id','kind','Q_eV','coefficient','convention','radiative','photon_energy_eV','estimated','source'])
        history=[]
        for snapshot in self.model.snapshots:
            for z in range(self.model.cells):
                for reaction in self.reactions:
                    k=reaction['channel']
                    history.append(dict(time_s=snapshot['time_s'],cell=z+1,id=reaction['id'],events_cm3_s=float(snapshot['events'][z,k]),integral_events_cm3=float(snapshot['integral'][z,k])))
        write_csv(Path(output)/'chemical_event_history.csv',history,['time_s','cell','id','events_cm3_s','integral_events_cm3'])
        coverage=[]
        for a,b in combinations_with_replacement(self.model.data.elements,2):
            active=[]
            for r in self.reactions:
                elements={e for n in r['reactants'] if n!='e' for e in self.model.states[self.model.index[n]]['nuclei']}
                if {a,b}==elements:active.append(r['id'])
            coverage.append(dict(element_A=a,element_B=b,implemented_channels=len(active),ids=' | '.join(active),status='partial' if active else 'missing',
                                 note='Inventory only: all other state-resolved, molecular and three-body reactions remain unproven/missing.'))
        write_csv(Path(output)/'pair_coverage.csv',coverage,['element_A','element_B','implemented_channels','ids','status','note'])
        write_csv(Path(output)/'chemical_known_omissions.csv',self.omissions,['id','reason','source'])
        write_csv(Path(output)/'chemical_domain_audit.csv',list(self.domain_audit.values()),['id','encountered_min_K','encountered_max_K','maximum_pair_density_cm6','out_of_domain','supported_intervals_K'])
        write_json(Path(output)/'chemical_data_provenance.json',dict(files_sha256=self.files,missing=sorted(self.gaps),known_omissions=self.omissions,models=sorted(self.notes)))
