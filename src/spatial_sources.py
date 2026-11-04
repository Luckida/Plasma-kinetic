"""Independent sheet/finite-width particle sources and internal photon sheets.

Particle strength is a column rate [cm^-2 s^-1]. A point is deposited into its
containing finite-volume cell, so its volume rate is strength/dz, not strength.
No neutralizing electrons are added unless an electron source is declared.
"""
from copy import deepcopy
import numpy as np
from atomic_database import write_csv
from mixture_data import read_csv
from spectrum_input import driving_rows
from light_input import light_rows,light_enabled
from particle_kinetics import ME, EV, KB

HPLANCK=6.62607015e-27
CLIGHT=2.99792458e10


def deposition(position,width,length,cells):
    if not np.isfinite(position+width) or not 0<=position<=length or width<0:
        raise ValueError('Источник: 0<=position_cm<=L, width_cm>=0')
    edges=np.linspace(0,length,cells+1);dz=length/cells
    if width==0:
        weights=np.zeros(cells);weights[min(int(position/dz),cells-1)]=1/dz
    else:
        lo,hi=position-width/2,position+width/2
        if lo<0 or hi>length:raise ValueError('Источник конечной ширины должен целиком лежать в слое')
        weights=np.maximum(0,np.minimum(edges[1:],hi)-np.maximum(edges[:-1],lo))/(width*dz)
    return weights


def electron_distribution(grid,options,root):
    opts=deepcopy(options)
    if 'uz_cm_s' in opts or 'beams' in opts:
        temperature=opts.pop('Te_eV',1.)*EV/KB;opts['Ti_K']=temperature
        if 'beams' in opts:
            for beam in opts['beams']:beam['Ti_K']=beam.pop('Te_eV',temperature*KB/EV)*EV/KB
        return grid.ion_initial(opts,1.,root)
    return grid.initial(opts,1.,root)


class SpatialSources:
    def __init__(self,model):
        self.model=model;self.items=[];self.times=[];self.participants=set()
        raw=model.cfg.get('spatial_sources',[])
        if not isinstance(raw,list):raise ValueError('spatial_sources: список источников')
        identifiers=set()
        for item in raw:
            allowed={'id','mode','enabled','position_cm','width_cm','start_s','stop_s','rate_cm2_s',
                     'species','distribution','spectrum_file','spectrum_options','scale','direction','mu','density_cm3','velocity_cm_s',
                     'normalization','multiplier_W_cm2','spectrum_mode','minimum_nm','maximum_nm','points'}
            if set(item)-allowed:raise ValueError('Неизвестные параметры источника: '+str(set(item)-allowed))
            if not isinstance(item.get('enabled',True),bool):raise ValueError('source.enabled: bool')
            if not item.get('enabled',True):continue
            if item.get('mode')=='radiation' and not light_enabled(item):continue
            row=deepcopy(item);ident=row.get('id','')
            if not isinstance(ident,str) or not ident or ident in identifiers:raise ValueError('Нужен уникальный id источника')
            identifiers.add(ident)
            mode=row.get('mode');position=float(row.get('position_cm',0.));width=float(row.get('width_cm',0.))
            if mode not in ('ions','electrons','radiation','particles'):raise ValueError('source.mode: ions, electrons, radiation или particles')
            row['weights']=deposition(position,width,model.length,model.cells)
            row['position_cm']=position
            start=float(row.get('start_s',0.));stop=float(row.get('stop_s',model.end))
            if not np.isfinite(start+stop) or not 0<=start<stop<=model.end:raise ValueError('Источник: 0<=start_s<stop_s<=t_end_s')
            row.update(start_s=start,stop_s=stop);self.times.extend([start,stop])
            if mode=='radiation':
                if width!=0:raise ValueError('Излучение: источник-лист имеет width_cm=0')
                if any(k in row for k in ('rate_cm2_s','species','distribution')):raise ValueError('Излучение: задайте spectrum_file и scale')
                mu=float(row.get('mu',1.));scale=float(row.get('scale',1.));direction=row.get('direction','both')
                if row.get('normalization')=='integral':scale=1.
                if not np.isfinite(mu+scale) or not 0<mu<=1 or scale<0 or direction not in ('left','right','both'):raise ValueError('Источник излучения: неверные mu/scale/direction')
                data=light_rows(row,model.root,model.cfg.get('probe'))
                lam=np.array([float(r['wavelength_A']) for r in data]);flux=np.array([float(r['four_pi_J_lambda_erg_cm2_s_A']) for r in data])
                if len(lam)<2 or np.any(np.diff(lam)<=0) or np.any(lam<=0) or np.any(flux<0) or not np.all(np.isfinite(lam+flux)):raise ValueError('Неверный спектр внутреннего источника')
                energy=HPLANCK*CLIGHT/(lam*1e-8)/EV;occupied=energy[flux>0]
                if len(occupied) and (occupied.min()<model.energy[0] or occupied.max()>model.energy[-1]):raise ValueError('Сетка фотонов не покрывает внутренний источник')
                wavelength=CLIGHT/model.nu*1e8
                row['spectrum']=np.interp(wavelength,lam,flux,left=0.,right=0.)*(wavelength*1e-8)**2/CLIGHT*1e8*scale
                if row.get('normalization')=='integral':
                    target=float(row['multiplier_W_cm2'])*1e7/mu;actual=float(row['spectrum']@model.w)
                    if actual<=0:raise ValueError('Сетка фотонов не разрешает форму источника.')
                    row['spectrum']*=target/actual
                row.update(mu=mu,direction=direction)
            else:
                if any(k in row for k in ('spectrum_file','scale','direction','mu')):raise ValueError('Источник частиц: задайте rate_cm2_s и distribution')
                strength=float(row.get('rate_cm2_s',0.))
                if not np.isfinite(strength) or strength<0:raise ValueError('rate_cm2_s>=0')
                row['rate_cm2_s']=strength
                distribution=deepcopy(row.get('distribution',{}))
                if 'density_cm3' in row:distribution['match_moments']=True
                if mode=='electrons':
                    if row.get('species','e') not in ('e','e-'):raise ValueError('Режим electrons создаёт только электроны')
                    row['index']=None;row['phase']=electron_distribution(model.ge,distribution,model.root)
                    row['energy_eV']=float((row['phase']*model.ge.e[:,None]).sum())
                else:
                    j=model.data.resolve(row['species']);state=model.states[j]
                    if mode=='ions' and state['charge']<=0:raise ValueError('Режим ions: нужен положительный ион')
                    row['index']=j;row['phase']=model.gh[state['element']].ion_initial(distribution,1.,model.root)
                    row['energy_eV']=float((row['phase']*model.gh[state['element']].e[:,None]).sum()+model.internal[j])
                    self.participants.add(j)
                from energy_input import distribution_kind
                if 'density_cm3' in row and distribution_kind(distribution)=='phase_table':
                    grid=model.ge if mode=='electrons' else model.gh[state['element']]
                    actual=grid.moments(row['phase'])['uz_cm_s'];requested=float(row['velocity_cm_s'])
                    if abs(actual-requested)>1e-3*max(abs(requested),float(grid.v[-1])*1e-9):raise ValueError(f'{ident}: средняя скорость фазового файла {actual:.6g} см/с не совпадает с заданной {requested:.6g} см/с. Исправьте файл или скорость источника.')
            self.items.append(row)
        self.reset()

    def reset(self):
        m=self.model
        self.heavy=np.zeros((m.cells,m.ns));self.electrons=np.zeros(m.cells)
        self.column=np.zeros(len(self.items));self.energy_erg_cm2=0.

    def checkpoint(self):return self.heavy.copy(),self.electrons.copy(),self.column.copy(),self.energy_erg_cm2
    def restore(self,state):self.heavy,self.electrons,self.column,self.energy_erg_cm2=state

    def active(self,row,t):return row['start_s']<=t<row['stop_s']

    def apply(self,t,dt):
        for k,row in enumerate(self.items):
            if row['mode']=='radiation' or not self.active(row,t):continue
            column=row['rate_cm2_s']*dt;density=column*row['weights'];phase=density[:,None,None]*row['phase'][None,:,:]
            if row['index'] is None:
                self.model.F+=phase;self.electrons+=density
            else:
                self.model.P[:,row['index']]+=phase;self.heavy[:,row['index']]+=density
            self.column[k]+=column;self.energy_erg_cm2+=column*row['energy_eV']*EV

    def photons(self,t):
        beams=[]
        for row in self.items:
            if row['mode']!='radiation' or not self.active(row,t):continue
            directions=('left','right') if row['direction']=='both' else (row['direction'],)
            for direction in directions:
                beams.append(dict(position_cm=row['position_cm'],direction=direction,mu=row['mu'],
                                  spectrum=row['spectrum']/len(directions)))
        return beams

    def nuclei(self):return self.heavy.sum(axis=0)@self.model.nuclear_counts*self.model.dz
    def charge(self):return float((self.heavy@self.model.charges-self.electrons).sum()*self.model.dz)

    def selection_seeds(self):
        seeds=[]
        for row in self.items:
            if row['mode'] not in ('ions','particles'):continue
            seed=np.zeros(self.model.ns);seed[row['index']]=row['rate_cm2_s']*(row['stop_s']-row['start_s'])/self.model.length
            seeds.append(seed)
        return seeds

    def export(self,folder):
        rows=[]
        for snapshot in self.model.snapshots:
            for k,item in enumerate(self.items):
                rows.append(dict(time_s=snapshot['time_s'],id=item['id'],mode=item['mode'],position_cm=item['position_cm'],
                    start_s=item['start_s'],stop_s=item['stop_s'],active=self.active(item,snapshot['time_s']),
                    injected_particles_cm2=float(snapshot['source_column'][k]) if item['mode']!='radiation' else '',
                    photon_power_erg_cm2_s=float(item['mu']*item['spectrum']@self.model.w) if item['mode']=='radiation' and self.active(item,snapshot['time_s']) else 0.))
        write_csv(folder/'spatial_source_history.csv',rows,['time_s','id','mode','position_cm','start_s','stop_s','active','injected_particles_cm2','photon_power_erg_cm2_s'])
