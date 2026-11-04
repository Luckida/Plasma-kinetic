"""Coupled kinetic populations, electrons and heavy species in a thin slab.

Every included heavy level has its own F_s[z,E,mu]. All elements share the
same evolving electrons and radiation. No hydrodynamic pressure closure.
"""
from project_version import VERSION
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import sys
import time
from types import SimpleNamespace
import numpy as np
from scipy.special import voigt_profile
from scipy.optimize import brentq
from threadpoolctl import threadpool_limits
from numba.typed import List as TypedList
from atomic_core import Network,SYMBOLS
from atomic_database import HC_EV_CM,write_csv,write_json
from mixture_data import MixtureData,read_csv,ROOT
from mixture_photo import PhotoBank
from mixture_chemistry import Chemistry,electron_release_bound
from spatial_sources import SpatialSources
from experiment_config import prepare_configuration
from spectrum_input import driving_rows
from light_input import light_rows,light_enabled
from wall_boundaries import WallBoundaries
from level_selection import LevelSelection,frozen_evolution
from particle_kinetics import PhaseGrid,GridOverflow,ME,EV,KB,QE,stream_periodic,stream_open
from collision_cache import CollisionMaps
from collision_evolution import evolve_frozen
from slab_transport import solve_transfer
import kinetic_kernels as K
import radiation as R
import diagnostics as D
from progress import Progress,progress_from_config
from grid_validation import validate_configuration_grids

HPLANCK=6.62607015e-27
CLIGHT=2.99792458e10
AMU=1.66053906892e-24
LINE_WING_WIDTHS=300.


def isolated_line_profile(energy,center,width,gamma):
    """Voigt density per Hz in the resolved isolated-line band.

    energy/center are in eV; width (sigma_D) and gamma are in Hz.

    The existing line quadrature spans +/-300 max(sigma_D, gamma). Extending
    its Lorentz tails across the *continuum* grid creates unmodelled far-wing
    photons and photoionization. This is a numerical isolated-line window,
    not a measured physical cutoff or a model of far-wing redistribution.
    """
    half=LINE_WING_WIDTHS*max(width,gamma)*HPLANCK/EV
    profile=voigt_profile((energy-center)*EV/HPLANCK,width,gamma)
    return np.where(abs(energy-center)<=half,profile,0.)


def conservative_profile(profile,weights,energy,target):
    """Positive spectral quadrature with unit photon count and exact mean E.

    A minimal exponential tilt corrects finite-grid quadrature asymmetry.
    This projects a photon source, never an evolved atomic population.
    """
    mass=np.asarray(profile)*weights;good=mass>0
    if not good.any() or not energy[good].min()<=target<=energy[good].max():raise GridOverflow('Спектральная сетка не покрывает центр профиля')
    x=energy[good]-target;log=np.log(mass[good])
    def probabilities(beta):
        value=log+beta*x;value=np.exp(value-value.max());return value/value.sum()
    def defect(beta):return float(probabilities(beta)@x)
    if abs(defect(0.))<=1e-14*max(abs(target),1.):beta=0.
    else:
        bound=1/max(float(np.ptp(energy[good])),1e-100);left=-bound;right=bound
        for _ in range(100):
            if defect(left)<=0 and defect(right)>=0:break
            left*=2;right*=2
        else:raise GridOverflow('Не удалось согласовать фотоны и энергию профиля')
        beta=brentq(defect,left,right,xtol=1e-14)
    out=np.zeros_like(profile);out[good]=probabilities(beta)/weights[good];return out


class RetryStep(RuntimeError):pass
class ExpandBasis(RuntimeError):
    def __init__(self,mask):self.mask=mask


class MixtureSlab:
    def __init__(self,config,root=ROOT,progress=None):
        validate_configuration_grids(config)
        config=prepare_configuration(config,root)
        self.cfg=deepcopy(config);self.root=Path(root);self.progress=progress or Progress()
        allowed={'schema','atomic_data','ion_stages_file','level_selection','initial_cm3','electron_initial','heavy_initial',
            'phase_grid','heavy_phase_grid','cells','thickness_cm','cell_overrides','solver','radiation','processes',
            'collisions','electric_field','particle_boundary','performance','diagnostics','output_directory','progress','chemistry','spatial_sources','gas','probe'}
        if config.get('schema') not in ('MixtureKineticSlab-1.8','MixtureKineticSlab-2.1') or set(config)-allowed:raise ValueError('Нужна схема MixtureKineticSlab-2.1; неизвестные поля '+str(set(config)-allowed))
        self.end=float(config.get('solver',{}).get('t_end_s',1e-7))
        if not np.isfinite(self.end) or self.end<=0:raise ValueError('t_end_s>0')
        self.cells=config.get('cells',1);self.length=float(config.get('thickness_cm',.01))
        if not isinstance(self.cells,int) or self.cells<1 or self.length<=0 or not np.isfinite(self.length):raise ValueError('Неверная геометрия слоя')
        self.dz=self.length/self.cells;self.edges=np.arange(1,self.cells+1)*self.dz
        perf=config.get('performance',{});self.workers=perf.get('workers',4)
        if set(perf)-{'workers','radiation_backend','collision_cache_mb'} or not isinstance(self.workers,int) or self.workers<1:raise ValueError('Неверные performance')
        self.data=MixtureData(config,self.root,self.progress)
        self.progress.start('Подготовка смеси','Энергетические сетки и распределения частиц')
        from phase_grid_setup import resolve_heavy_grid
        heavy_phase=resolve_heavy_grid(config,self.data,self.root)
        phase=deepcopy(config.get('phase_grid',{}));automatic_electron_grid=phase.get('energy_max_eV','auto')=='auto'
        if phase.get('energy_max_eV','auto')=='auto':
            opts=config.get('electron_initial',{});T=float(opts.get('Te_eV',1.))
            if float(opts.get('hot_fraction',0.))>0:T=max(T,float(opts.get('hot_Te_eV',5.)))
            speed=abs(float(opts.get('uz_cm_s',0.)))+np.sqrt(16*T*EV/ME)
            speed+=QE*1e7*abs(float(config.get('electric_field',{}).get('Ez_V_cm',0.)))/ME*self.end
            maximum=ME*speed*speed/(2*EV)
            for source in config.get('spatial_sources',[]):
                if not source.get('enabled',True) or source.get('mode')!='electrons':continue
                eopts=source.get('distribution',{})
                for beam in eopts.get('beams',[eopts]):
                    Tsrc=float(beam.get('Te_eV',1.))
                    if float(beam.get('hot_fraction',0.))>0:Tsrc=max(Tsrc,float(beam.get('hot_Te_eV',Tsrc)))
                    vsrc=abs(float(beam.get('uz_cm_s',0.)))+np.sqrt(16*Tsrc*EV/ME)
                    maximum=max(maximum,ME*vsrc*vsrc/(2*EV))
            from energy_input import energy_support,distribution_kind
            inputs=[config.get('electron_initial',{})]+[source.get('distribution',{}) for source in config.get('spatial_sources',[]) if source.get('enabled',True) and source.get('mode')=='electrons']
            inputs.extend(local.get('electron_initial',config.get('electron_initial',{})) for local in config.get('cell_overrides',[]))
            for options in inputs:
                if distribution_kind(options)=='energy_table':
                    low,high=energy_support(options,self.root);maximum=max(maximum,high)
                    if low>0:phase['energy_min_eV']=min(phase.get('energy_min_eV',1e-4),low)
            radiation=config.get('radiation',{})
            if light_enabled(radiation):
                rows=light_rows(radiation,self.root,config.get('probe'));energies=[HPLANCK*CLIGHT/(float(r['wavelength_A'])*1e-8)/EV for r in rows if float(r['four_pi_J_lambda_erg_cm2_s_A'])>0]
                if energies:maximum=max(maximum,max(energies)+8*T)
            for source in config.get('spatial_sources',[]):
                if source.get('mode')=='radiation' and light_enabled(source):
                    rows=light_rows(source,self.root,config.get('probe'))
                    energies=[HPLANCK*CLIGHT/(float(r['wavelength_A'])*1e-8)/EV for r in rows if float(r['four_pi_J_lambda_erg_cm2_s_A'])>0]
                    if energies:maximum=max(maximum,max(energies)+8*T)
            phase['energy_max_eV']=1.1*maximum
        self.ge=PhaseGrid(phase)
        if automatic_electron_grid:
            release=electron_release_bound(self.data,config,self.root)
            if release>0:
                phase['energy_max_eV']=max(phase['energy_max_eV'],1.1*(release+2*heavy_phase.get('energy_max_eV',50.)))
                self.ge=PhaseGrid(phase)
        self.states=self.data.states;self.index=self.data.index;self.ns=len(self.states)
        self.gh={element:PhaseGrid(heavy_phase,
            next(s['mass_u'] for s in self.states if s['element']==element)*AMU, label='heavy_phase_grid') for element in sorted({s['element'] for s in self.states})}
        self.initial=np.zeros((self.cells,self.ns));self.F=np.zeros((self.cells,*self.ge.shape))
        hg=next(iter(self.gh.values()));self.P=np.zeros((self.cells,self.ns,*hg.shape))
        level_distributions={}
        for name,options in config.get('heavy_initial',{}).items():
            if name in self.index or ':' in name or name.startswith('HeI_'):level_distributions[self.data.resolve(name)]=options
        overrides=config.get('cell_overrides',[])
        if overrides and len(overrides)!=self.cells:raise ValueError('cell_overrides: по одному объекту на ячейку')
        for z in range(self.cells):
            local=overrides[z] if overrides else {}
            for name,value in local.get('initial_cm3',config.get('initial_cm3',{})).items():
                if not np.isfinite(value) or value<0:raise ValueError('Неверная начальная населённость '+name)
                self.initial[z,self.data.resolve(name)]+=value
            for j,s in enumerate(self.states):
                opts=level_distributions.get(j,config.get('heavy_initial',{}).get(s['species'],config.get('heavy_initial',{}).get(s['element'],{'Ti_K':300.})))
                self.P[z,j]=self.gh[s['element']].ion_initial(opts,self.initial[z,j],self.root)
            eopts=local.get('electron_initial',config.get('electron_initial',{}))
            ne=eopts.get('density_cm3',float(self.initial[z]@np.array([s['charge'] for s in self.states])))
            if 'uz_cm_s' in eopts or 'beams' in eopts:
                opts=dict(eopts);opts['Ti_K']=opts.pop('Te_eV',1.)*EV/KB
                if 'beams' in opts:opts['beams']=[dict(b,Ti_K=b.pop('Te_eV',opts['Ti_K']*KB/EV)*EV/KB) for b in deepcopy(opts['beams'])]
                self.F[z]=self.ge.ion_initial(opts,ne,self.root)
            else:self.F[z]=self.ge.initial(eopts,ne,self.root)
        if self.initial.sum()<=0:raise ValueError('Нужна непустая смесь')
        self.charges=np.array([s['charge'] for s in self.states]);self.internal=np.array([s['energy_cm1']*HC_EV_CM for s in self.states])
        self.nuclear_counts=np.array([[s.get('nuclei',{s['element']:1}).get(e,0) for e in self.data.elements] for s in self.states])
        self.process=config.get('processes',{})
        process_keys={'electron_bound_bound','electron_ionization','spontaneous','two_photon','photo_bound_bound','photo_bound_free','radiative_recombination','stimulated_recombination','disabled_ids',
            'penning','associative_ionization','charge_exchange','dissociative_recombination','mutual_neutralization','heavy_chemistry','heavy_quenching'}
        if set(self.process)-process_keys:raise ValueError('Неизвестные processes '+str(set(self.process)-process_keys))
        if any(not isinstance(v,bool) for k,v in self.process.items() if k!='disabled_ids'):raise ValueError('Флаги процессов должны быть bool')
        if not isinstance(self.process.get('disabled_ids',[]),list) or any(not isinstance(v,str) for v in self.process.get('disabled_ids',[])):raise ValueError('processes.disabled_ids: список строк')
        self.progress.start('Подготовка смеси','Спектр и атомные операторы')
        self._spectrum()
        from probe_spectrum import Probe
        self.probe=Probe(self)
        if any(self.enabled(key) for key in ('photo_bound_free','radiative_recombination','stimulated_recombination')):
            self.photo=PhotoBank(self.root,self.states,self.data.parent_index,float(self.energy[-1]),config.get('atomic_data',{}).get('offline',False),self.progress,config.get('radiation',{}).get('photoionization_file'))
        else:self.photo=SimpleNamespace(tables={},sources={},gaps=[],receipts={})
        self.channels=[];self.impact=[];self.lines=[];self.cont=[];self.gaps=[]
        self._atomic_channels()
        self.sources=SpatialSources(self)
        self.chemistry=Chemistry(self)
        self.nr=len(self.channels);self._network()
        self._pairs()
        self.coll=config.get('collisions',{})
        if set(self.coll)-{'ee_bgk','ii_bgk','ei_binary','neutral_elastic','neutral_low_energy','coulomb_log','maximum_phase_loss','relative_azimuth_points','scattering_azimuth_points','uniformization_tail_probability'}:raise ValueError('Неизвестные collisions')
        if self.coll.get('neutral_low_energy','error') not in ('error','constant'):raise ValueError('neutral_low_energy: error или constant (явно выбранная модель продолжения)')
        if any(not isinstance(self.coll.get(key,False),bool) for key in ('ee_bgk','ii_bgk','ei_binary','neutral_elastic')):raise ValueError('Столкновительные флаги должны быть bool')
        self.loss_limit=float(self.coll.get('maximum_phase_loss',.2))
        if not 0<self.loss_limit<=.5:raise ValueError('maximum_phase_loss: (0,0.5]')
        self.policy=config.get('level_selection',{}).get('missing_data','error')
        if self.policy not in ('error','report'):raise ValueError('missing_data: error или report')
        self.bound=config.get('particle_boundary',{'kind':'periodic'})
        if set(self.bound)-{'kind','left_electrons','right_electrons','left_heavy','right_heavy','walls'}:raise ValueError('Неизвестные particle_boundary')
        if self.bound.get('kind','periodic') not in ('periodic','vacuum','reservoir','walls'):raise ValueError('particle_boundary.kind: periodic, vacuum, reservoir или walls')
        self.walls=WallBoundaries(self)
        self.boundary_P=self.P.copy();self.boundary_F=self.F.copy()
        self.incoming_F=[];self.incoming_P=[]
        for side,z in [('left',0),('right',self.cells-1)]:
            eopts=self.bound.get(side+'_electrons')
            self.incoming_F.append(self.ge.initial(eopts,float(eopts.get('density_cm3',0.)),self.root) if eopts else self.F[z].copy())
            heavy=self.P[z].copy()
            for name,opts in self.bound.get(side+'_heavy',{}).items():
                j=self.data.resolve(name);heavy[j]=self.gh[self.states[j]['element']].ion_initial(opts,float(opts.get('density_cm3',0.)),self.root)
            self.incoming_P.append(heavy)
        self._reset()
        fields=[self.field(0.,all_levels=True)]
        rates=np.array([self.rates(z,fields[0]['mean_4piJ'][z])[0] for z in range(self.cells)])
        selection_initial=self.initial
        if self.bound.get('kind')=='reservoir':
            selection_initial=np.vstack([self.initial,*[p.sum(axis=(1,2))[None,:] for p in self.incoming_P]])
        seeds=self.sources.selection_seeds()
        if seeds:selection_initial=np.vstack([selection_initial,*seeds])
        with threadpool_limits(limits=1):
            selection_rates=rates.copy()
            selection_rates[:,[r['channel'] for r in self.chemistry.reactions]]=0.
            self.selection=LevelSelection(self.states,self.channels,selection_initial,selection_rates,self.end,config.get('level_selection'))
            self.selection.active[list(self.chemistry.participants)]=True
            self.selection.active[list(self.sources.participants)]=True
            if self.walls.enabled:self.selection.active[list(set(self.walls.ground.values()))]=True
            self.selection.close_radiative_cascades()
        self.initial_F=self.F.copy();self.initial_P=self.P.copy()
        self.progress.done(f'Смесь: {len(self.data.elements)} элементов, {self.ns} доступных / {self.selection.active.sum()} выбранных состояний')

    @property
    def n(self):return self.P.sum(axis=(2,3))

    def enabled(self,key):return self.process.get(key,True)

    def add(self,ident,i,j,kind,source='',**extra):
        if ident in self.process.get('disabled_ids',[]):return None
        reaction=dict(id=ident,i=i,j=j,kind=kind,source=source,**extra)
        self.channels.append(reaction);return len(self.channels)-1

    def _spectrum(self):
        cfg=self.cfg.get('radiation',{})
        allowed={'spectrum_file','scale','minimum_eV','maximum_eV','continuum_points','line_points','Tn_K','angular_points','backend','beam_mu','beam_side','left_diffuse_scale','right_diffuse_scale','history_file','photoionization_file','electron_scattering','enabled','spectrum_options','start_s','stop_s','normalization','multiplier_W_cm2','spectrum_mode','minimum_nm','maximum_nm','points'}
        if set(cfg)-allowed:raise ValueError('Неизвестные radiation '+str(set(cfg)-allowed))
        self.Tn=float(cfg.get('Tn_K',300.));self.scale=float(cfg.get('scale',1.))
        if cfg.get('normalization')=='integral':self.scale=1. if light_enabled(cfg) else 0.
        if not cfg.get('enabled',True):self.scale=0.
        if not np.isfinite(self.Tn+self.scale) or self.Tn<=0 or self.scale<0:raise ValueError('Температура излучателя >0 K, масштаб накачки >=0; оба числа конечны.')
        if light_enabled(cfg) and self.scale>0:
            rows=light_rows(cfg,self.root,self.cfg.get('probe'));lam=np.array([float(r['wavelength_A']) for r in rows]);fl=np.array([float(r['four_pi_J_lambda_erg_cm2_s_A']) for r in rows])
        else:lam=np.array([1.,1e9]);fl=np.zeros(2)
        if len(lam)<2 or np.any(np.diff(lam)<=0) or np.any(lam<=0) or np.any(fl<0) or not np.all(np.isfinite(lam+fl)):raise ValueError('Неверный спектр накачки')
        energies=HPLANCK*CLIGHT/(lam*1e-8)/EV
        source_energies=[];source_edges=[]
        if np.any(fl>0):source_edges.extend([energies.min()*(1-1e-10),energies.max()*(1+1e-10)])
        for item in self.cfg.get('spatial_sources',[]):
            if item.get('mode')=='radiation' and light_enabled(item):
                rows=light_rows(item,self.root,self.cfg.get('probe'))
                local=[HPLANCK*CLIGHT/(float(r['wavelength_A'])*1e-8)/EV for r in rows]
                source_energies.extend(local);source_edges.extend([min(local)*(1-1e-10),max(local)*(1+1e-10)])
        high=max([1.]+[(self.states[j]['energy_cm1']-self.states[i]['energy_cm1'])*HC_EV_CM for i,j in self.data.parent_index.items()])+self.ge.e[-1]+next(iter(self.gh.values())).e[-1]+1.
        nonzero=np.r_[energies[fl>0],source_energies];high=max(high,float(nonzero.max())*1.001 if len(nonzero) else 1.)
        line_energies=[self.internal[self.index[r['upper']]]-self.internal[self.index[r['lower']]] for r in self.data.lines]
        positive=[e for e in line_energies if e>0]
        high=max(high,max(positive,default=0.)*1.001)
        high=float(cfg.get('maximum_eV',high));low=float(cfg.get('minimum_eV',min(1e-4,.5*min(positive,default=1.))))
        if not 0<low<high or (len(nonzero) and (nonzero.min()<low or nonzero.max()>high)):raise ValueError('Полоса излучения не покрывает спектр')
        parts=[np.geomspace(low,high,int(cfg.get('continuum_points',400))),energies[(energies>=low)&(energies<=high)],np.asarray(source_energies),np.asarray(source_edges)]
        width=[]
        for row in self.data.lines:
            i,j=self.index[row['lower']],self.index[row['upper']];e=self.internal[j]-self.internal[i];nu=e*EV/HPLANCK
            if not low<e<high:continue
            sigma=nu/CLIGHT*np.sqrt(KB*self.Tn/(self.states[j]['mass_u']*AMU))
            gamma=float(row['A_s'])/(4*np.pi)
            if row.get('spectral_kind')=='two_photon':
                parts.append(np.linspace(low,e-low,100));continue
            u=np.r_[-np.geomspace(5,LINE_WING_WIDTHS,12)[::-1],np.linspace(-5,5,int(cfg.get('line_points',41))),np.geomspace(5,LINE_WING_WIDTHS,12)]
            parts.append((nu+u*max(sigma,gamma))*HPLANCK/EV)
            # Zero guard nodes prevent trapezoidal cells across large continuum
            # gaps from assigning their weight to the last nonzero wing node.
            half=LINE_WING_WIDTHS*max(sigma,gamma)*HPLANCK/EV
            parts.append(e+np.array([-1.,1.])*half*(1+1e-6))
        for i,j in self.data.parent_index.items():
            edge=self.internal[j]-self.internal[i];parts.extend([np.array([np.nextafter(edge,0.),edge]),edge+self.ge.e])
        raw=np.unique(np.concatenate(parts));self.energy=raw[(raw>=low)&(raw<=high)];self.nu=self.energy*EV/HPLANCK;self.w=R.weights(self.nu)
        wavelength=CLIGHT/self.nu*1e8
        self.source=np.interp(wavelength,lam,fl,left=0.,right=0.)*(wavelength*1e-8)**2/CLIGHT*1e8*self.scale
        if cfg.get('normalization')=='integral' and self.scale>0:
            target=float(cfg['multiplier_W_cm2'])*1e7/float(cfg.get('beam_mu',1.));actual=float(self.source@self.w)
            if actual<=0:raise ValueError('Сетка фотонов не разрешает форму накачки.')
            self.source*=target/actual
        self.knots=np.array([0.]);self.scales=np.array([1.])
        if cfg.get('history_file') and self.scale>0:
            rows=read_csv(self.root/cfg['history_file']);self.knots=np.array([float(r['time_s']) for r in rows]);self.scales=np.array([float(r['scale']) for r in rows])
            if len(rows)<2 or np.any(np.diff(self.knots)<=0) or self.knots[0]!=0 or np.any(self.scales<0):raise ValueError('Неверная история накачки')
        self.pump_start=float(cfg.get('start_s',0.));self.pump_stop=float(cfg.get('stop_s',self.end))
        if self.scale>0:
            if not np.isfinite(self.pump_start+self.pump_stop) or not 0<=self.pump_start<self.pump_stop<=self.end:raise ValueError('Накачка: 0<=начало<конец<=время расчёта')
            knots=np.unique(np.r_[self.knots,self.pump_start,self.pump_stop])
            self.scales=np.interp(knots,self.knots,self.scales);self.knots=knots

    def _atomic_channels(self):
        self.sv=[];self.xs_unknown=[]
        for (a,b,kind),table in self.data.cross.items():
            if kind not in ('excitation','deexcitation','ionization'):continue
            if not self.enabled('electron_ionization' if kind=='ionization' else 'electron_bound_bound'):continue
            i,j=self.index[a],self.index[b];ident=f'CE:{a}->{b}' if kind!='ionization' else f'ION:{a}'
            k=self.add(ident,i,j,kind,self.data.cross_source[(a,b,kind)])
            if k is None:continue
            gap=self.internal[j]-self.internal[i]
            sig=np.interp(self.ge.e,table[:,0],table[:,1],left=0.,right=0.)
            sig[self.ge.e<max(gap,0.)]=0.
            unknown=(self.ge.e>table[-1,0])|((self.ge.e<table[0,0])&(self.ge.e>=max(gap,0.)))
            if kind=='deexcitation' and (b,a,'excitation') in self.data.cross:
                upward=self.data.cross[(b,a,'excitation')];energy=self.ge.e-gap
                sig=self.states[j]['g']/self.states[i]['g']*energy/self.ge.e*np.interp(energy,upward[:,0],upward[:,1],left=0.,right=0.)
                unknown=(energy<upward[0,0])|(energy>upward[-1,0]);self.channels[k]['source']+='; inverse at runtime nodes from excitation via detailed balance'
            self.sv.append(sig*self.ge.v);self.xs_unknown.append(unknown);self.impact.append(k)
        self.sv=np.array(self.sv).reshape(len(self.impact),len(self.ge.e));self.xs_unknown=np.array(self.xs_unknown).reshape(len(self.impact),len(self.ge.e))
        self.absorb=np.zeros((self.ns,len(self.nu)));self.emit=np.zeros_like(self.absorb)
        for row in self.data.lines:
            i,j=self.index[row['lower']],self.index[row['upper']];two=row.get('spectral_kind')=='two_photon'
            gap=self.internal[j]-self.internal[i];A=float(row['A_s']);nu0=gap*EV/HPLANCK
            if two:
                if not self.enabled('two_photon'):continue
                y=self.energy/gap
                if self.states[j]['element']=='He' and self.states[j]['charge']==0:
                    shape=read_csv(ROOT/'data/radiation/two_photon_shape.csv');x=np.array([float(r['y']) for r in shape]);v=np.array([float(r['shape']) for r in shape]);profile=np.interp(np.minimum(y,1-y),x,v,left=0.,right=0.)
                else:
                    Z=SYMBOLS.index(self.states[j]['element'])+1
                    shape=[r for r in read_csv(ROOT/'data/radiation/hydrogenic_two_photon.csv') if int(r['Z'])==Z]
                    profile=np.interp(np.minimum(y,1-y),[float(r['y']) for r in shape],[float(r['shape']) for r in shape],left=0.,right=0.)
            else:
                sigma=nu0/CLIGHT*np.sqrt(KB*self.Tn/(self.states[j]['mass_u']*AMU));profile=isolated_line_profile(self.energy,gap,sigma,A/(4*np.pi))
            if profile@self.w<=0:continue
            try:profile=conservative_profile(profile,self.w,self.energy,gap/2 if two else gap)
            except GridOverflow as exc:raise GridOverflow(f'{row["id"]}: E_transition={gap:.12g} eV; '+str(exc)) from exc
            sp=self.add('A:'+row['id'],j,i,'two_photon' if two else 'spontaneous',row['source'],A=A) if self.enabled('two_photon' if two else 'spontaneous') else None
            if sp is not None:self.emit[j]+=A*HPLANCK*self.nu*profile*(2 if two else 1)
            optical=None
            if not two and row.get('pump',True) and self.enabled('photo_bound_bound'):
                sigma=A*(CLIGHT/nu0)**2/(8*np.pi)*profile;ratio=self.states[j]['g']/self.states[i]['g']
                up=self.add('PUMP:'+row['id'],i,j,'photoexcitation',row['source']);down=self.add('STIM:'+row['id'],j,i,'stimulated_emission',row['source'])
                if up is not None:self.absorb[i]+=ratio*sigma
                if down is not None:self.absorb[j]-=sigma
                optical=(up,down,sigma,ratio)
            self.lines.append(dict(row,i=i,j=j,profile=profile,sp=sp,optical=optical,two=two))
        for i,j in self.data.parent_index.items():
            ident=self.states[i]['id'];threshold=self.internal[j]-self.internal[i]
            if ident not in self.photo.tables:continue
            table=self.photo.tables[ident];sigma=np.interp(self.energy,table[:,0],table[:,1],left=0.,right=0.);sigma[self.energy<threshold]=0.
            pi=self.add('PI:'+ident,i,j,'photoionization',self.photo.sources[ident]) if self.enabled('photo_bound_free') else None
            rr=self.add('RR:'+ident,j,i,'radiative_recombination',self.photo.sources[ident]+'; Milne') if self.enabled('radiative_recombination') else None
            sr=self.add('SR:'+ident,j,i,'stimulated_recombination',self.photo.sources[ident]+'; Milne') if self.enabled('stimulated_recombination') else None
            if pi is None and rr is None and sr is None:continue
            if pi is not None:self.absorb[i]+=sigma
            self.cont.append(dict(i=i,j=j,threshold=threshold,table=table,sigma=sigma,pi=pi,rr=rr,sr=sr))

    def _network(self):
        states=[dict(s,charge=int(s['charge'])) for s in self.states]+[dict(id='e',label='e⁻',nuclei={},charge=-1),dict(id='photon',label='hν',nuclei={},charge=0)]
        from reaction_data import equation
        reactions=[equation(self.states,r) for r in self.channels]
        self.network=Network(states,reactions);self.matrix=self.network.matrix[:self.ns];self.electron_change=self.network.matrix[self.ns]
        for r,label,equation in zip(self.channels,self.network.labels,reactions):
            r['reaction']=label;r['reactants']=dict(equation['reactants']);r['products']=dict(equation['products'])
        self.linear=np.array([k for k,r in enumerate(self.channels) if not r.get('chemistry') and r['kind'] not in ('radiative_recombination','stimulated_recombination')],dtype=int)

    def _pairs(self):
        self.pairs={};self.maps={};self.rr_groups={}
        cfg=self.cfg.get('collisions',{});naz=int(cfg.get('relative_azimuth_points',2));nphi=int(cfg.get('scattering_azimuth_points',2))
        if naz<2 or nphi<2:raise ValueError('Азимутальные квадратуры >=2')
        az=2*np.pi*(np.arange(naz)+.5)/naz;self.scatter_phi=2*np.pi*(np.arange(nphi)+.5)/nphi
        for element,g in self.gh.items():
            self.progress.module('Парные скорости и RR: '+element,force=True)
            rows=[r for r in self.cont if self.states[r['i']]['element']==element]
            self.rr_groups[element]=rows
            if not rows and not cfg.get('ei_binary',False) and not cfg.get('neutral_elastic',False):continue
            with K.kernel_threads(self.workers):
                rel,cx,cy,cz=K.pair_geometry(self.ge.vz.ravel(),self.ge.vperp.ravel(),g.vz.ravel(),g.vperp.ravel(),az,g.mass)
            pair=SimpleNamespace(rel=rel,cx=cx,cy=cy,cz=cz,ee=np.repeat(self.ge.e,len(self.ge.mu)),ie=np.repeat(g.e,len(g.mu)))
            self.pairs[element]=pair
            if rows:
                with K.kernel_threads(self.workers):
                    pair.rate,pair.index,pair.mix=K.build_rr(rel,np.array([r['threshold'] for r in rows]),np.array([self.states[r['i']]['g']/self.states[r['j']]['g'] for r in rows]),
                        TypedList([np.ascontiguousarray(r['table'][:,0]) for r in rows]),TypedList([np.ascontiguousarray(r['table'][:,1]) for r in rows]),self.energy,HPLANCK,CLIGHT,g.mass,1.)
                for k,row in enumerate(rows):
                    outside=row['threshold']+rel.ravel()>row['table'][-1,0];pair.rate[k,outside]=0.;pair.index[k,outside]=-1
                    row['pair_row']=k
            self.rr_groups[element]=rows
            if cfg.get('ei_binary',False) or cfg.get('neutral_elastic',False):
                self.maps[element]=CollisionMaps(pair,self.ge,g,self.ge.mu,self.ge.wmu,self.scatter_phi,
                    self.cfg.get('performance',{}).get('collision_cache_mb',256)/len(self.gh),self.workers,self.root,self.progress)

    def _reset(self):
        self.integral=np.zeros((self.cells,len(self.channels)));self.transport=np.zeros((self.cells,self.ns));self.electron_out=0.;self.species_out=np.zeros(self.ns)
        self.electron_transport=np.zeros(self.cells)
        self.energy_out=0.;self.field_work=0.;self.photon_work=0.;self.profile_heat=0.;self.projection_error=0.
        self.chemical_momentum_error=0.
        self.steps=0;self.rejected=0;self.data_issues=set();self.model_notes={'Kinetic bound-bound Voigt profiles use an isolated-line window of +/-300 max(sigma_D,gamma); far-wing continuum/redistribution is not modelled.'};self.leak=np.zeros(self.ns);self.snapshots=[]
        self.sources.reset();self.walls.reset()
        self.initial_energy=self.total_energy();self.initial_charge=self.total_charge();self.initial_nuclei=self.totals()

    def totals(self):return self.n.sum(axis=0)@self.nuclear_counts*self.dz
    def total_charge(self):return float((self.n@self.charges-self.F.sum(axis=(1,2))).sum()*self.dz)
    def total_energy(self):
        return float((self.F*self.ge.e[None,:,None]).sum()+ (self.n@self.internal).sum()+sum((self.P[:,j]*self.gh[s['element']].e[None,:,None]).sum() for j,s in enumerate(self.states)))*EV*self.dz

    def field_E(self):
        cfg=self.cfg.get('electric_field',{});mode=cfg.get('mode','prescribed')
        if mode=='prescribed':return np.full(self.cells,float(cfg.get('Ez_V_cm',0.)))
        if mode!='poisson' or self.bound.get('kind','periodic')=='periodic':raise ValueError('Poisson: нужны открытые границы vacuum и потенциалы')
        rho=(self.n@self.charges-self.F.sum(axis=(1,2)))*QE*1e6;dz=self.dz*.01
        A=2*np.eye(self.cells)-np.eye(self.cells,k=1)-np.eye(self.cells,k=-1);A[0,0]+=1;A[-1,-1]+=1
        b=rho*dz*dz/8.8541878128e-12;b[0]+=2*cfg.get('left_potential_V',0.);b[-1]+=2*cfg.get('right_potential_V',0.)
        potential=np.linalg.solve(A,b);faces=np.r_[-2*(potential[0]-cfg.get('left_potential_V',0.))/dz,-np.diff(potential)/dz,-2*(cfg.get('right_potential_V',0.)-potential[-1])/dz]
        return .5*(faces[:-1]+faces[1:])/100

    def rr_values(self,z,psi):
        result={};f=self.F[z];ne=f.sum();occupation=psi*CLIGHT**2/(8*np.pi*HPLANCK*self.nu**3)
        for element,rows in self.rr_groups.items():
            if not rows:continue
            pair=self.pairs[element];g=self.gh[element]
            for parent in sorted({r['j'] for r in rows}):
                ions=self.P[z,parent];ni=ions.sum()
                active=np.array([r['j']==parent for r in rows])
                if ne and ni and np.any(pair.rel==0):raise GridOverflow('Точная нулевая относительная скорость в RR-квадратуре; измените фазовые сетки')
                with K.kernel_threads(1 if len(rows)<32 else self.workers):
                    spec,el,il,pe,cm,missing=K.rr_emission(f.ravel(),ions.ravel(),pair.rate,pair.index,pair.mix,pair.rel,pair.ee,pair.ie,self.energy,active)
                    es,iss,ps,cms=K.rr_stimulated(f.ravel(),ions.ravel(),pair.rate,pair.index,pair.mix,pair.rel,pair.ee,pair.ie,self.energy,occupation,active)
                for k,row in enumerate(rows):
                    if row['j']!=parent:continue
                    result[id(row)]=(spec[k],el[k],il[k],pe[k],cm[k],es[k],iss[k],ps[k],cms[k],missing[k])
        return result

    def field(self,t,all_levels=False):
        self.progress.module('Общий перенос излучения всех элементов')
        opacity=[];emission=[]
        active=np.ones(self.ns,dtype=bool) if all_levels or not hasattr(self,'selection') else self.selection.active
        cache_key=active.tobytes()
        cached=getattr(self,'_optical_cache',None)
        if cached is not None and cached[0]==cache_key:absorb,emit=cached[1:]
        else:
            absorb,emit=self.optical_matrices(active)
            self._optical_cache=(cache_key,absorb,emit)
        for z in range(self.cells):
            n=self.n[z]*active;ka=n@absorb;eps=n@emit
            eps+=self.chemistry.photon_emission(z)
            rr=self.rr_values(z,np.zeros_like(self.nu));ne=self.F[z].sum()
            for row in self.cont:
                if not active[row['i']] or not active[row['j']]:continue
                values=rr[id(row)];spec=values[0];ion=n[row['j']]
                if row['rr'] is not None:eps+=ne*ion*HPLANCK*self.nu*spec/self.w
                if row['sr'] is not None:ka-=ne*ion*spec/self.w*CLIGHT**2/(8*np.pi*self.nu**2)
            opacity.append(ka);emission.append(eps)
        cfg=self.cfg.get('radiation',{});scale=np.interp(t,self.knots,self.scales)
        if not self.pump_start<=t<self.pump_stop:scale=0.
        scatter=self.F.sum(axis=(1,2))*6.6524587051e-25 if cfg.get('electron_scattering',True) else None
        result=solve_transfer(self.edges,np.array(opacity),np.array(emission),self.source*scale,scatter=scatter,
            angular_points=cfg.get('angular_points',4),beam_mu=cfg.get('beam_mu',1.),beam_side=cfg.get('beam_side','left'),
            left_4piI=self.source*cfg.get('left_diffuse_scale',0.)*scale,right_4piI=self.source*cfg.get('right_diffuse_scale',0.)*scale,
            backend=self.cfg.get('performance',{}).get('radiation_backend','numpy'),workers=self.workers,internal_beams=self.sources.photons(t))
        result['opacity']=np.array(opacity);result['emission']=np.array(emission)
        return result

    def optical_matrices(self,active):
        absorb=np.zeros_like(self.absorb);emit=np.zeros_like(self.emit)
        for line in self.lines:
            i,j=line['i'],line['j']
            if not active[i] or not active[j]:continue
            if line['sp'] is not None:emit[j]+=float(line['A_s'])*HPLANCK*self.nu*line['profile']*(2 if line['two'] else 1)
            if line['optical']:
                up,down,sigma,ratio=line['optical']
                if up is not None:absorb[i]+=ratio*sigma
                if down is not None:absorb[j]-=sigma
        for row in self.cont:
            if active[row['i']] and active[row['j']] and row['pi'] is not None:absorb[row['i']]+=row['sigma']
        return absorb,emit

    def rates(self,z,psi):
        rates=np.zeros(len(self.channels));energies=np.zeros_like(rates);counts=self.F[z].sum(axis=1)
        with K.kernel_threads(self.workers):rates[self.impact]=K.impact_rates(self.sv,counts)
        for line in self.lines:
            if line['sp'] is not None:
                rates[line['sp']]=line['A_s'];energies[line['sp']]=-(2 if line['two'] else 1)*float((line['profile']*self.w)@self.energy)
            if line['optical']:
                up,down,sigma,ratio=line['optical'];raw=sigma*psi/(HPLANCK*self.nu)*self.w
                if up is not None:rates[up]=ratio*raw.sum();energies[up]=float(raw@self.energy/raw.sum()) if raw.sum() else 0.
                if down is not None:rates[down]=raw.sum();energies[down]=-float(raw@self.energy/raw.sum()) if raw.sum() else 0.
        rr=self.rr_values(z,psi);ne=counts.sum();n=self.n[z]
        for row in self.cont:
            raw=row['sigma']*psi/(HPLANCK*self.nu)*self.w;pi=row['pi'];values=rr[id(row)]
            if pi is not None:rates[pi]=raw.sum();energies[pi]=float(raw@self.energy/raw.sum()) if raw.sum() else 0.
            if row['rr'] is not None:rates[row['rr']]=ne*values[1].sum();energies[row['rr']]=-values[3]/max(values[1].sum(),1e-300)
            if row['sr'] is not None:rates[row['sr']]=ne*values[5].sum();energies[row['sr']]=-values[7]/max(values[5].sum(),1e-300)
        for index,event_rate in self.chemistry.event_rates(z).items():
            density=n[self.channels[index]['i']]
            rates[index]=event_rate/density if density else 0.
        return rates,energies,rr

    def _check_data(self,z,psi):
        n=self.n[z];reference={e:max(self.initial[:,[s['element']==e for s in self.states]].sum(),1.) for e in self.data.elements}
        relevant=[i for i,s in enumerate(self.states) if self.selection.active[i] and n[i]>self.selection.tolerance*self.selection.reference[s['element']]]
        for k,index in enumerate(self.impact):
            r=self.channels[index]
            if r['i'] in relevant and (self.F[z].sum(axis=1)[self.xs_unknown[k]].sum()/max(self.F[z].sum(),1.))>self.selection.tolerance:
                self.data_issues.add('sigma(E) does not cover occupied electron energies: '+r['id'])
        for i in relevant:
            state=self.states[i];a=state['id']
            if state.get('molecular'):
                self.model_notes.add('Effective ground carrier, no resolved rovibrational kinetics/opacity: '+a)
                continue
            same_stage=[s for s in self.states if s['species']==state['species']]
            Z=SYMBOLS.index(state['element'])+1
            if len(same_stage)==1 and Z-state['charge']==1:
                edge=.75*13.605693122994*Z*Z
                tail=self.F[z].sum(axis=1)[self.ge.e>=edge].sum()/max(self.F[z].sum(),1.)
                if tail>self.selection.tolerance and self.enabled('electron_bound_bound'):self.data_issues.add('excited H-like ion states accessible but only ground is loaded: '+state['species'])
            for j,target in enumerate(self.states):
                if state['species']!=target['species'] or i==j or not self.enabled('electron_bound_bound'):continue
                kind='excitation' if self.internal[j]>=self.internal[i] else 'deexcitation'
                threshold=max(self.internal[j]-self.internal[i],0.)
                tail=self.F[z].sum(axis=1)[self.ge.e>=threshold].sum()/max(self.F[z].sum(),1.)
                if tail>self.selection.tolerance and (a,target['id'],kind) not in self.data.cross:self.data_issues.add('missing electron channel '+a+' -> '+target['id'])
            parent=self.data.parent_index.get(i)
            if parent is not None:
                gap=self.internal[parent]-self.internal[i];tail=self.F[z].sum(axis=1)[self.ge.e>=gap].sum()/max(self.F[z].sum(),1.)
                if self.enabled('electron_ionization') and tail>self.selection.tolerance and (a,self.states[parent]['id'],'ionization') not in self.data.cross:self.data_issues.add('missing ionization sigma(E): '+a)
                if self.enabled('photo_bound_free') and np.any(psi[self.energy>=gap]>0) and a not in self.photo.tables:self.data_issues.add('missing photoionization: '+a)
        if self.enabled('radiative_recombination'):
            for i,j in self.data.parent_index.items():
                if self.selection.active[j] and n[j]>self.selection.tolerance*self.selection.reference[self.states[i]['element']] and self.states[i]['id'] not in self.photo.tables:
                    self.data_issues.add('missing level-resolved RR: '+self.states[i]['id'])
        for row in self.cont:
            i,j=row['i'],row['j']
            if not self.selection.active[i] or not self.selection.active[j]:continue
            if i in relevant and row['pi'] is not None:
                photons=psi/(HPLANCK*self.nu)*self.w;ionizing=self.energy>=row['threshold']
                outside=self.energy>row['table'][-1,0]*(1+1e-12)
                if float(photons[outside].sum())>self.selection.tolerance*max(float(photons[ionizing].sum()),1e-300):self.data_issues.add('PI spectrum outside cross-section domain: '+self.states[i]['id'])
            if (row['rr'] is not None or row['sr'] is not None) and n[j]>self.selection.tolerance*self.selection.reference[self.states[i]['element']]:
                pair=self.pairs[self.states[i]['element']];weights=self.F[z].ravel()[:,None,None]*self.P[z,j].ravel()[None,:,None]
                outside=row['threshold']+pair.rel>row['table'][-1,0]
                if float((weights*outside).sum())>self.selection.tolerance*max(float(weights.sum())*pair.rel.shape[2],1.):self.data_issues.add('RR relative energies outside PI table: '+self.states[i]['id'])
        if self.data_issues and self.policy=='error':raise ValueError('Атомные данные не замкнуты: '+sorted(self.data_issues)[0]+'; см. missing_data=report только для явно неполной исследовательской модели')

    def upper_boundary_report(self):
        rows=[]
        for species in sorted({s['species'] for s in self.states}):
            indices=[j for j,s in enumerate(self.states) if s['species']==species]
            if len(indices)<=1:continue
            highest=max(self.states[j].get('n',0) for j in indices);top={j for j in indices if self.states[j].get('n',0)==highest}
            events=sum(float(self.integral[:,k].sum()) for k,r in enumerate(self.channels) if r['j'] in top and r['i'] not in top)
            element=self.states[indices[0]]['element'];fraction=events/self.selection.reference[element]
            rows.append(dict(species=species,highest_available_n=highest,inflow_to_highest_shell_cm3=events,fraction_of_smallest_seed=fraction,exposed=fraction>self.selection.tolerance))
        return rows

    def _binary(self,z,dt):
        if self.coll.get('ee_bgk',True):self.F[z],_=self.ge.ee_bgk(self.F[z],dt,self.coll.get('coulomb_log',10.))
        for j,s in enumerate(self.states):
            ions=self.P[z,j];ni=ions.sum();ne=self.F[z].sum()
            if not ni or not ne:continue
            charged=s['charge']!=0
            if not self.coll.get('ei_binary' if charged else 'neutral_elastic',False):continue
            element=s['element'];pair=self.pairs[element];g=self.gh[element]
            if charged:
                moment=self.ge.moments(self.F[z]);screen=743*np.sqrt(max(moment['T_moment_eV'],1e-20)/ne)
                with K.kernel_threads(self.workers):sv=K.pair_rates(pair.rel,screen,self.coll.get('coulomb_log',10.),g.mass,s['charge'])
            else:
                if s['id'] not in self.data.elastic:
                    self.data_issues.add('missing neutral elastic sigma(E): '+s['id']);continue
                table,_,_=self.data.elastic[s['id']];reduced=ME*g.mass/(ME+g.mass)
                weights=self.F[z].ravel()[:,None,None]*ions.ravel()[None,:,None]
                low=pair.rel<table[0,0];outside=pair.rel>table[-1,0]
                continuation=self.coll.get('neutral_low_energy','error')=='constant'
                if not continuation:outside|=low
                elif np.any(low):self.model_notes.add('constant low-energy elastic continuation below '+str(table[0,0])+' eV: '+s['id'])
                if ni>self.selection.tolerance*self.selection.reference[element] and float((weights*outside).sum())>self.selection.tolerance*max(float(weights.sum())*pair.rel.shape[2],1.):self.data_issues.add('neutral elastic sigma does not cover pair energies: '+s['id'])
                sv=np.interp(pair.rel,table[:,0],table[:,1],left=table[0,1] if continuation else 0.,right=0.)*np.sqrt(2*pair.rel*EV/reduced)
            cache=self.maps.get(element)
            if cache is not None and cache.maps:gain,ig,il,frequency=cache.generator_rates(ions.ravel(),sv)
            else:
                with K.kernel_threads(self.workers):gain,ig,il,frequency=K.binary_generator_rates(ions.ravel(),pair.rel,pair.cx,pair.cy,pair.cz,self.ge.e,self.ge.mu,g.e,g.mu,sv,self.ge.mu,self.ge.wmu,self.scatter_phi,g.mass)
            evolved,residence,_,_,_=evolve_frozen(gain,frequency,self.F[z].ravel(),dt,self.coll.get('uniformization_tail_probability',1e-14))
            loss=il.T@residence;addition=ig.T@residence
            if np.any(loss-addition>self.loss_limit*ions.ravel()):raise RetryStep('Большой фазовый сток бинарных столкновений')
            self.F[z]=evolved.reshape(self.ge.shape);self.P[z,j]=(ions.ravel()+addition-loss).reshape(g.shape)
        if self.coll.get('ii_bgk',False):
            for j,s in enumerate(self.states):
                if not s['charge'] or not self.P[z,j].sum():continue
                g=self.gh[s['element']];moment=g.moments(self.P[z,j]);T=moment['T_moment_eV']
                if T<=0:raise ValueError('ii BGK требует ненулевого разброса скоростей')
                nu=4.80e-8*s['charge']**4/np.sqrt(g.mass/1.67262192595e-24)*moment['ne_cm3']*self.coll.get('coulomb_log',10.)/T**1.5
                w=-np.expm1(-nu*dt);self.P[z,j]=(1-w)*self.P[z,j]+w*g.maximum_entropy(self.P[z,j])

    def _local(self,z,dt,psi):
        self.progress.module('Совместные реакции атомов и ионов')
        self._check_data(z,psi);rates,phenergy,_=self.rates(z,psi);mask=self.selection.channel_mask()
        leaked=self.selection.leakage(self.n[z],rates,dt);self.leak+=leaked
        expand=self.selection.candidates(self.leak)
        if expand.any():raise ExpandBasis(expand)
        coefficients=rates[self.linear]*mask[self.linear]
        active=np.flatnonzero(self.selection.active);local={j:k for k,j in enumerate(active)}
        used=[(k,c) for k,c in zip(self.linear,coefficients) if c>0]
        linear=[dict(self.channels[k],i=local[self.channels[k]['i']],j=local[self.channels[k]['j']]) for k,c in used]
        evolved,integrated=frozen_evolution(self.P[z,active].reshape(len(active),-1),linear,np.array([c for k,c in used]),dt)
        new=self.P[z].reshape(self.ns,-1).copy();new[active]=evolved
        residence=np.zeros_like(new);residence[active]=integrated
        events=np.zeros(self.nr)
        for k,c in zip(self.linear,coefficients):events[k]=c*residence[self.channels[k]['i']].sum()
        before=self.F[z].copy()
        with K.kernel_threads(self.workers):
            loss,gain,err=K.impact_events(before,events[self.impact],self.sv,
                np.array([self.internal[self.channels[k]['j']]-self.internal[self.channels[k]['i']] for k in self.impact]),
                np.array([2. if self.channels[k]['kind']=='ionization' else 1. for k in self.impact]),self.ge.e,self.ge.wmu)
        lost=loss.sum(axis=0);added=gain.sum(axis=0)
        for row in self.cont:
            if row['pi'] is None or not events[row['pi']]:continue
            raw=row['sigma']*psi/(HPLANCK*self.nu)*self.w;total=raw.sum()
            if total:
                deposited,error=self.ge.deposit(events[row['pi']]*raw/total,np.maximum(self.energy-row['threshold'],0.),allow_tail=True)
                added+=deposited;self.projection_error+=error*EV*self.dz
        if np.any(lost>self.loss_limit*before):raise RetryStep('Большой неупругий фазовый сток электронов')
        self.F[z]=before-lost+added;self.P[z]=new.reshape(self.P[z].shape)
        self.projection_error+=float(err.sum())*EV*self.dz
        for k,r in enumerate(self.channels):
            if not events[k]:continue
            if r['kind'] in ('spontaneous','two_photon','stimulated_emission','photoexcitation'):
                gap=self.internal[r['j']]-self.internal[r['i']]
                self.profile_heat+=events[k]*(phenergy[k]-gap)*EV*self.dz
        self.photon_work+=float(events@phenergy)*EV*self.dz
        # RR is selective in BOTH phase distributions; one frozen step with
        # the same event counts in populations, photons and charged particles.
        rr=self.rr_values(z,psi);ne=self.F[z].sum();occupation=psi*CLIGHT**2/(8*np.pi*HPLANCK*self.nu**3)
        eloss=np.zeros_like(self.F[z]);hloss=np.zeros_like(self.P[z]);hgain=np.zeros_like(self.P[z])
        for element,rows in self.rr_groups.items():
            if not rows:continue
            pair=self.pairs[element];g=self.gh[element]
            for parent in sorted({r['j'] for r in rows}):
                ni=self.P[z,parent].sum()
                with K.kernel_threads(self.workers):
                    prod,cm=K.rr_heavy_products(self.F[z].ravel(),self.P[z,parent].ravel(),pair.rate,pair.index,pair.mix,pair.rel,pair.cx,pair.cy,pair.cz,g.e,g.mu,g.mass,occupation,False)
                    prods,cms=K.rr_heavy_products(self.F[z].ravel(),self.P[z,parent].ravel(),pair.rate,pair.index,pair.mix,pair.rel,pair.cx,pair.cy,pair.cz,g.e,g.mu,g.mass,occupation,True)
                for row in rows:
                    if row['j']!=parent or not self.selection.active[row['i']] or not self.selection.active[parent]:continue
                    values=rr[id(row)];krow=row['pair_row'];factor=ne*ni*dt
                    for index,ee,ii,pe,product in [(row['rr'],values[1],values[2],values[3],prod[krow]),(row['sr'],values[5],values[6],values[7],prods[krow])]:
                        if index is None:continue
                        event=factor*ee.sum();events[index]+=event
                        eloss+=factor*ee.reshape(self.ge.shape);hloss[parent]+=factor*ii.reshape(g.shape);hgain[row['i']]+=factor*product.reshape(g.shape)
                        self.photon_work-=factor*pe*EV*self.dz
        if np.any(eloss>self.loss_limit*self.F[z]) or np.any(hloss>self.loss_limit*self.P[z]):raise RetryStep('Большой RR фазовый сток')
        self.F[z]-=eloss;self.P[z]+=hgain-hloss
        chemical_events=self.chemistry.evolve(z,dt)
        if chemical_events is None:raise RetryStep('Большой фазовый сток межкомпонентной химии')
        events+=chemical_events
        if np.any(self.F[z]<0) or np.any(self.P[z]<0):raise RetryStep('Отрицательный фазовый узел')
        self.integral[z]+=events

    def _transport(self,dt):
        self.progress.module('Ускорение и кинетический перенос всех видов')
        before_electrons=self.F.sum(axis=(1,2)).copy()
        E=self.field_E()
        for z in range(self.cells):
            for array,g,charge in [(self.F,self.ge,-1)]+[(self.P[:,j],self.gh[s['element']],s['charge']) for j,s in enumerate(self.states) if s['charge']]:
                before=g.moments(array[z]);dv=charge*QE*1e7*E[z]/g.mass*dt
                array[z],err,_=g.kick(array[z],dv);self.field_work+=g.mass*before['ne_cm3']*(before['uz_cm_s']*dv+.5*dv*dv)*self.dz;self.projection_error+=err*EV*self.dz
        if self.walls.enabled:
            self.walls.apply(dt);self.electron_transport+=self.F.sum(axis=(1,2))-before_electrons;return
        periodic=self.bound.get('kind','periodic')=='periodic'
        if periodic:self.F=stream_periodic(self.F,self.ge.vz,dt,self.dz)
        else:
            left,right=self.incoming_F if self.bound.get('kind')=='reservoir' else (np.zeros(self.ge.shape),np.zeros(self.ge.shape))
            self.F,out=stream_open(self.F,self.ge.vz,dt,self.dz,left,right);self.electron_out+=out.sum();self.energy_out+=float((out*self.ge.e[:,None]).sum())*EV
        for j,s in enumerate(self.states):
            g=self.gh[s['element']];old=self.P[:,j].sum(axis=(1,2))
            if periodic:self.P[:,j]=stream_periodic(self.P[:,j],g.vz,dt,self.dz)
            else:
                left,right=(self.incoming_P[0][j],self.incoming_P[1][j]) if self.bound.get('kind')=='reservoir' else (np.zeros(g.shape),np.zeros(g.shape))
                self.P[:,j],out=stream_open(self.P[:,j],g.vz,dt,self.dz,left,right);self.species_out[j]+=out.sum();self.energy_out+=float((out*(g.e[:,None]+self.internal[j])).sum())*EV
            self.transport[:,j]+=self.P[:,j].sum(axis=(1,2))-old
        self.electron_transport+=self.F.sum(axis=(1,2))-before_electrons

    def snapshot(self,t):
        field=self.field(t);events=[];photon_residual=0.
        for z in range(self.cells):
            rates,_,_=self.rates(z,field['mean_4piJ'][z]);events.append(rates*self.n[z,[r['i'] for r in self.channels]]*self.selection.channel_mask())
            photons=self.network.matrix[-1];spectral=float(((field['opacity'][z]*field['mean_4piJ'][z]-field['emission'][z])/(HPLANCK*self.nu))@self.w)
            predicted=-float(photons@events[-1]);norm=max(float(abs(photons)@events[-1]),1e-100)
            photon_residual=max(photon_residual,abs(spectral-predicted)/norm)
        return dict(time_s=float(t),F=self.F.copy(),P=self.P.copy(),n=self.n.copy(),events=np.array(events),integral=self.integral.copy(),transport=self.transport.copy(),electron_transport=self.electron_transport.copy(),
            charge=self.total_charge()-self.initial_charge+float(self.species_out@self.charges-self.electron_out)-self.sources.charge(),
            nuclei=self.totals()-self.initial_nuclei+self.species_out@self.nuclear_counts-self.sources.nuclei(),
            energy=self.total_energy()-self.initial_energy+self.energy_out+self.profile_heat-self.field_work-self.photon_work-self.sources.energy_erg_cm2,
            source_heavy=self.sources.heavy.copy(),source_electrons=self.sources.electrons.copy(),source_column=self.sources.column.copy(),source_energy=self.sources.energy_erg_cm2,
            wall_counts=self.walls.counts.copy(),wall_heat=self.walls.heat.copy(),wall_charge=self.walls.charge.copy(),
            field=field,field_work=self.field_work,photon_work=self.photon_work,profile_heat=self.profile_heat,projection_error=self.projection_error,reaction_photon_residual=photon_residual,
            chemical_momentum_error=self.chemical_momentum_error)

    def solve(self):
        opt=self.cfg.get('solver',{});maximum=float(opt.get('max_step_s',self.end/20));minimum=float(opt.get('min_step_s',maximum*1e-8));samples=int(opt.get('samples',11))
        if not 0<minimum<=maximum<=self.end or samples<2:raise ValueError('Неверные solver')
        targets=np.unique(np.r_[np.linspace(0,self.end,samples),opt.get('extra_times_s',[]),self.knots[(self.knots>=0)&(self.knots<=self.end)],self.sources.times])
        if targets[0]<0 or targets[-1]>self.end:raise ValueError('extra_times_s вне интервала')
        start=time.monotonic()
        with threadpool_limits(limits=1),K.kernel_threads(self.workers):
            while True:
                self.F=self.initial_F.copy();self.P=self.initial_P.copy();self._reset();t=0.;step=maximum
                self.progress.start('Интегрирование смеси','Общие электроны и отдельные тяжёлые распределения',f'активных состояний {self.selection.active.sum()}/{self.ns}; повторов базиса {self.selection.restarts}')
                self.progress.update(0,self.end);self.snapshots=[self.snapshot(0.)]
                try:
                    for target in targets[1:]:
                        while t<target-1e-14*self.end:
                            dt=min(step,target-t)
                            if self.cfg.get('electric_field',{}).get('mode')=='poisson':dt=min(dt,.15/(5.64146e4*np.sqrt(max(float(self.F.sum(axis=(1,2)).max()),1.))))
                            if dt<minimum:raise RuntimeError('Шаг ниже min_step_s')
                            saved=(self.F.copy(),self.P.copy(),self.integral.copy(),self.transport.copy(),self.electron_transport.copy(),self.species_out.copy(),self.leak.copy(),
                                {key:getattr(self,key) for key in ('electron_out','energy_out','field_work','photon_work','profile_heat','projection_error','chemical_momentum_error')})
                            source_saved=self.sources.checkpoint();wall_saved=self.walls.checkpoint()
                            try:
                                self.sources.apply(t,dt/2)
                                self._transport(dt)
                                for z in range(self.cells):self._binary(z,dt)
                                field=self.field(t)
                                for z in range(self.cells):self._local(z,dt,field['mean_4piJ'][z])
                                self.sources.apply(t,dt/2)
                            except RetryStep:
                                self.F,self.P,self.integral,self.transport,self.electron_transport,self.species_out,self.leak,scalars=saved
                                for key,value in scalars.items():setattr(self,key,value)
                                self.sources.restore(source_saved);self.walls.restore(wall_saved)
                                self.rejected+=1;step=dt/2;continue
                            t+=dt;self.steps+=1;step=min(maximum,max(step,dt)*1.25)
                            if self.steps>opt.get('max_steps',100000):raise RuntimeError('Превышено max_steps')
                            self.progress.update(t,self.end,f't={t:.5g}/{self.end:.5g} с; шагов {self.steps}, повторов {self.rejected}')
                        t=float(target)  # source switches occur on exact event times
                        self.snapshots.append(self.snapshot(target));self.progress.update(target,self.end,force=True)
                    for row in self.upper_boundary_report():
                        if row['exposed']:self.data_issues.add('highest available shell appreciably visited: '+row['species']+' n='+str(row['highest_available_n'])+'; need a larger atomic bank')
                    break
                except ExpandBasis as exc:
                    self.selection.expand(exc.mask);self.progress.module('Расширение базиса и повтор от исходного состояния',force=True)
            self.selection.omitted_integrals=self.leak.copy()
        self.elapsed=time.monotonic()-start;self.progress.done('Интегрирование завершено')
        return self.snapshots

    def export(self,output):
        if self.cfg.get('diagnostics',{}).get('output_mode','compact')!='full':
            from result_store import export_results
            return export_results(self,output)
        return self.export_full(output)

    def export_full(self,output):
        output=Path(output);output.mkdir(parents=True,exist_ok=True)
        self.progress.start('Диагностики смеси','Населённости, распределения и независимые балансы')
        population=[];moments=[];conservation=[];components=[];top=[];collisions=[];ephase=[];hphase=[];eedf=[];partner_collisions=[]
        phase_csv=self.cfg.get('diagnostics',{}).get('phase_history_csv',True)
        requested=self.cfg.get('diagnostics',{}).get('reaction_levels',list(self.cfg.get('initial_cm3',{})))
        targets=list(dict.fromkeys(self.data.resolve(name) for name in requested))
        for number,snap in enumerate(self.snapshots):
            t=snap['time_s']
            conservation.append(dict(time_s=t,charge_residual_e_cm2=snap['charge'],energy_residual_erg_cm2=snap['energy'],
                mass_residual_g_cm2=float(sum(value*next(s['mass_u'] for s in self.states if s['element']==element)*AMU for value,element in zip(snap['nuclei'],self.data.elements))),
                **{e+'_nuclei_residual_cm2':float(v) for e,v in zip(self.data.elements,snap['nuclei'])},
                electric_work_erg_cm2=snap['field_work'],photon_work_erg_cm2=snap['photon_work'],line_profile_reservoir_gain_erg_cm2=snap['profile_heat'],
                energy_residual_without_profile_reservoir_erg_cm2=snap['energy']-snap['profile_heat'],projection_error_erg_cm2=snap['projection_error'],
                chemical_projection_momentum_error_g_cm_s=snap['chemical_momentum_error']))
            for z in range(self.cells):
                population.append(dict(time_s=t,cell=z+1,**{s['id']:float(snap['n'][z,j]) for j,s in enumerate(self.states)}))
                em=self.ge.moments(snap['F'][z]);moments.append(dict(time_s=t,cell=z+1,id='e',element='',charge=-1,**em))
                enodes=snap['F'][z].sum(axis=1)
                for k,e in enumerate(self.ge.e):
                    eedf.append(dict(time_s=t,cell=z+1,energy_eV=float(e),node_density_cm3=float(enodes[k]),probability_per_eV=float(enodes[k]/max(em['ne_cm3'],1e-300)/self.ge.we[k])))
                    if phase_csv:
                        for a,mu in enumerate(self.ge.mu):ephase.append(dict(time_s=t,cell=z+1,energy_eV=float(e),mu=float(mu),speed_cm_s=float(self.ge.v[k]),vz_cm_s=float(self.ge.vz[k,a]),node_density_cm3=float(snap['F'][z,k,a])))
                for j,s in enumerate(self.states):
                    hm=self.gh[s['element']].moments(snap['P'][z,j]);moments.append(dict(time_s=t,cell=z+1,id=s['id'],element=s['element'],charge=s['charge'],**hm))
                    g=self.gh[s['element']]
                    if phase_csv:
                        for k,e in enumerate(g.e):
                            for a,mu in enumerate(g.mu):hphase.append(dict(time_s=t,cell=z+1,id=s['id'],energy_eV=float(e),mu=float(mu),speed_cm_s=float(g.v[k]),vz_cm_s=float(g.vz[k,a]),node_density_cm3=float(snap['P'][z,j,k,a])))
                    if hm['ne_cm3'] and em['ne_cm3'] and s['element'] in self.pairs:
                        pair=self.pairs[s['element']];sv=None;source=''
                        if s['charge']:
                            screen=743*np.sqrt(max(em['T_moment_eV'],1e-20)/em['ne_cm3'])
                            with K.kernel_threads(self.workers):sv=K.pair_rates(pair.rel,screen,self.coll.get('coulomb_log',10.),g.mass,s['charge'])
                            source='screened isotropic transport-equivalent Coulomb'
                        elif s['id'] in self.data.elastic:
                            table,source,kind=self.data.elastic[s['id']];left=table[0,1] if self.coll.get('neutral_low_energy','error')=='constant' else 0.
                            sv=np.interp(pair.rel,table[:,0],table[:,1],left=left,right=0.)*np.sqrt(2*pair.rel*EV/(ME*g.mass/(ME+g.mass)));source=kind+'; '+source
                        if sv is not None:
                            rate=float((snap['F'][z].ravel()[:,None,None]*snap['P'][z,j].ravel()[None,:,None]*sv).sum()/pair.rel.shape[2]);nu=rate/em['ne_cm3'];meanv=float(self.ge.v@enodes/em['ne_cm3'])
                            partner_collisions.append(dict(time_s=t,cell=z+1,partner=s['id'],partner_label=s['label'],events_cm3_s=rate,electron_frequency_s=nu,heavy_frequency_s=rate/hm['ne_cm3'],
                                electron_mean_free_path_cm=meanv/nu if nu else '',electron_Knudsen=meanv/nu/self.length if nu else '',
                                enabled=self.coll.get('ei_binary' if s['charge'] else 'neutral_elastic',False),stationary_inelastic_target_speed_ratio=abs(hm['uz_cm_s'])/max(meanv,1.),source=source))
                    defect=snap['n'][z,j]-self.initial[z,j]-self.matrix[j]@snap['integral'][z]-snap['transport'][z,j]-snap['source_heavy'][z,j]
                    components.append(dict(time_s=t,cell=z+1,id=s['id'],label=s['label'],source_cm3_s=float(np.maximum(self.matrix[j],0)@snap['events'][z]),
                        sink_cm3_s=float(np.maximum(-self.matrix[j],0)@snap['events'][z]),kinetic_residual_cm3=float(defect),transport_integral_cm3=float(snap['transport'][z,j]),external_source_integral_cm3=float(snap['source_heavy'][z,j])))
                for j in targets:
                    for direction,stoich in [('source',np.maximum(self.matrix[j],0)),('sink',np.maximum(-self.matrix[j],0))]:
                        values=stoich*snap['events'][z];total=float(values.sum());indices=np.argsort(-values)[:5];indices=indices[values[indices]>0]
                        coverage=float(values[indices].sum()/total) if total else 0.
                        for rank,k in enumerate(indices,1):
                            top.append(dict(time_s=t,cell=z+1,target=self.states[j]['id'],target_label=self.states[j]['label'],direction=direction,rank=rank,id=self.channels[k]['id'],reaction=self.channels[k]['reaction'],
                                events_cm3_s=float(values[k]),fraction_of_total=float(values[k]/total),total_cm3_s=total,top5_coverage=coverage,omitted_cm3_s=total-float(values[indices].sum()),
                                population_cm3=float(snap['n'][z,j]),net_cm3_s=float(self.matrix[j]@snap['events'][z])))
                inelastic=sum(snap['events'][z,k] for k,r in enumerate(self.channels) if r['kind'] in ('excitation','deexcitation','ionization'))
                collisions.append(dict(time_s=t,cell=z+1,electron_inelastic_frequency_s=inelastic/max(em['ne_cm3'],1.),electron_mean_speed_cm_s=float(self.ge.v@snap['F'][z].sum(axis=1)/max(em['ne_cm3'],1.)),
                    charged_states_with_binary_operator=sum(bool(s['charge']) for s in self.states) if self.coll.get('ei_binary',False) else 0,
                    ion_ion_model='per-state BGK; no cross-element heavy-heavy friction' if self.coll.get('ii_bgk',False) else 'off',
                    electron_ion_model='full pair distributions; isotropic transport-equivalent Coulomb cross section' if self.coll.get('ei_binary',False) else 'off'))
            self.progress.update(number+1,len(self.snapshots)+3)
        for name,rows,fields in [('populations_by_cell',population,None),('particle_moments',moments,None),('conservation',conservation,None),('component_budgets',components,None),
            ('reaction_top5',top,['time_s','cell','target','target_label','direction','rank','id','reaction','events_cm3_s','fraction_of_total','total_cm3_s','top5_coverage','omitted_cm3_s','population_cm3','net_cm3_s']),('collision_diagnostics',collisions,None)]:
            write_csv(output/(name+'.csv'),rows,fields or list(rows[0]))
        write_csv(output/'electron_energy_history.csv',eedf,list(eedf[0]))
        if phase_csv:
            write_csv(output/'electron_phase_history.csv',ephase,list(ephase[0]));write_csv(output/'heavy_phase_history.csv',hphase,list(hphase[0]))
        write_csv(output/'collision_partners.csv',partner_collisions,list(partner_collisions[0]) if partner_collisions else ['time_s','cell','partner','events_cm3_s','electron_frequency_s','heavy_frequency_s','enabled','source'])
        states=[{k:s.get(k,'') for k in ('id','label','element','species','charge','n','multiplicity','L','g','energy_cm1','excitation_eV','mass_u','configuration')} for s in self.states]
        for j,row in enumerate(states):
            row['active']=bool(self.selection.active[j]);row['molecular_carrier']=bool(self.states[j].get('molecular',False))
            row['nuclei']=json.dumps(self.states[j].get('nuclei',{self.states[j]['element']:1}),sort_keys=True)
        write_csv(output/'states.csv',states,list(states[0]));write_csv(output/'level_selection.csv',self.selection.report(),list(self.selection.report()[0]))
        bounds=[]
        for species in sorted({s['species'] for s in self.states}):
            included=[s for j,s in enumerate(self.states) if s['species']==species and self.selection.active[j]]
            available=[s for s in self.states if s['species']==species]
            bounds.append(dict(species=species,active_levels=len(included),available_levels=len(available),selected_max_n=max((s.get('n',0) for s in included),default=0),
                selected_excitation_limit_eV=max((s['excitation_eV'] for s in included),default=0.),data_max_n=max(s.get('n',0) for s in available)))
        write_csv(output/'selected_energy_limits.csv',bounds,list(bounds[0]))
        upper=self.upper_boundary_report();write_csv(output/'upper_data_boundary.csv',upper,['species','highest_available_n','inflow_to_highest_shell_cm3','fraction_of_smallest_seed','exposed'])
        energy_diagnostics=[]
        for number,snap in enumerate([self.snapshots[0],self.snapshots[-1]]):
            for z in range(self.cells):
                counts=snap['F'][z].sum(axis=1);total=counts.sum();mom=self.ge.moments(snap['F'][z]);cumulative=np.cumsum(counts)/total if total else np.zeros_like(counts)
                index=min(int(np.searchsorted(cumulative,1-self.selection.tolerance)),len(self.ge.e)-1)
                energy_diagnostics.append(dict(snapshot='initial' if number==0 else 'final',cell=z+1,electron_energy_quantile_eV=float(self.ge.e[index]) if total else 0.,
                    retained_probability=1-self.selection.tolerance,electron_mean_energy_eV=float(counts@self.ge.e/total) if total else 0.,electron_drift_energy_eV=ME*mom['uz_cm_s']**2/(2*EV),
                    source_max_nonzero_photon_eV=float(self.energy[self.source>0].max()) if np.any(self.source>0) else 0.))
        write_csv(output/'excitation_energy_diagnostics.csv',energy_diagnostics,list(energy_diagnostics[0]))
        catalogue=[{k:r[k] for k in ('id','reaction','kind','source')} for r in self.channels]
        write_csv(output/'reaction_catalogue.csv',catalogue,['id','reaction','kind','source'])
        for z in range(self.cells):
            folder=output/f'cell_{z+1:02d}'/'reaction_dynamics';folder.mkdir(parents=True,exist_ok=True);manifest=[]
            for k,history in enumerate(D.reaction_history([row for row in top if row['cell']==z+1])):
                stem=f'{k+1:02d}_'+history['direction'];columns=list(history['columns'].values())
                write_csv(folder/(stem+'_rates.csv'),history['rates'],['time_s','top5_coverage','total_cm3_s']+columns)
                write_csv(folder/(stem+'_fractions.csv'),history['fractions'],['time_s','top5_coverage']+columns)
                manifest.append(dict(target=history['target'],direction=history['direction'],rates_file=stem+'_rates.csv',fractions_file=stem+'_fractions.csv'))
            write_csv(folder/'files.csv',manifest,['target','direction','rates_file','fractions_file'])
        self.progress.module('Сжатие полных функций распределения',force=True)
        self.progress.update(len(self.snapshots)+1,len(self.snapshots)+3)
        gh=next(iter(self.gh.values()))
        np.savez_compressed(output/'mixture_history.npz',time_s=[s['time_s'] for s in self.snapshots],energy_eV=self.ge.e,mu=self.ge.mu,energy_weights_eV=self.ge.we,angular_weights=self.ge.wmu,electrons_cm3=[s['F'] for s in self.snapshots],
            heavy_energy_eV=gh.e,heavy_mu=gh.mu,heavy_phase_cm3=[s['P'] for s in self.snapshots],state_ids=np.array([s['id'] for s in self.states]),
            heavy_energy_weights_eV=gh.we,heavy_angular_weights=gh.wmu,populations_cm3=[s['n'] for s in self.snapshots],reaction_integrals_cm3=[s['integral'] for s in self.snapshots],
            position_cm=(np.arange(self.cells)+.5)*self.dz,external_source_integrals_cm3=[s['source_heavy'] for s in self.snapshots])
        final=self.snapshots[-1];self.F=final['F'].copy();self.P=final['P'].copy();field=final['field'];ranked=[]
        for z in range(self.cells):
            n=final['n'][z]
            for line in self.lines:
                i,j=line['i'],line['j'];row=dict(cell=z+1,transition=self.states[j]['label']+' → '+self.states[i]['label'],wavelength_nm=float(line['wavelength_nm']),spectral_kind=line['spectral_kind'],
                    emission_thin_erg_cm2_s_sr=float(n[j]*float(line['A_s'])*(self.internal[j]-self.internal[i])*EV*self.dz/(4*np.pi)) if line['sp'] is not None else 0.,isolated_equivalent_width_nm=0.,maximum_line_tau=0.)
                if line['optical']:
                    up,down,sigma,ratio=line['optical'];tau=((n[i]*ratio if up is not None else 0.)-(n[j] if down is not None else 0.))*sigma*self.dz
                    row['maximum_line_tau']=float(tau.max());wavelength=CLIGHT/self.nu*1e7
                    value=1-np.exp(-np.clip(tau,-50,700));row['isolated_equivalent_width_nm']=float(-np.sum((value[1:]+value[:-1])*np.diff(wavelength)/2))
                ranked.append(row)
        ranked.sort(key=lambda r:-r['emission_thin_erg_cm2_s_sr'])
        for k,r in enumerate(ranked,1):r['emission_rank']=k
        absorption=sorted(range(len(ranked)),key=lambda k:-ranked[k]['isolated_equivalent_width_nm'])
        for k,index in enumerate(absorption,1):ranked[index]['absorption_rank']=k
        write_csv(output/'lines_final.csv',ranked,list(ranked[0]) if ranked else ['transition','wavelength_nm','spectral_kind','emission_thin_erg_cm2_s_sr','isolated_equivalent_width_nm','maximum_line_tau','emission_rank','absorption_rank'])
        write_csv(output/'spectrum_final.csv',[dict(wavelength_A=float(CLIGHT/nu*1e8),outgoing_left_erg_cm2_s_Hz=float(field['outgoing_left_flux'][k]),outgoing_right_erg_cm2_s_Hz=float(field['outgoing_right_flux'][k]),
            direct_transmitted_erg_cm2_s_Hz=float(field['transmitted_direct_beam_flux'][k]),normal_optical_depth=float(field['normal_optical_depth'][k])) for k,nu in enumerate(self.nu)],
            ['wavelength_A','outgoing_left_erg_cm2_s_Hz','outgoing_right_erg_cm2_s_Hz','direct_transmitted_erg_cm2_s_Hz','normal_optical_depth'])
        norm=max(float(self.initial_nuclei.sum()),1.);energy_norm=max(abs(self.initial_energy),1e-100)
        charge=max(abs(s['charge']) for s in self.snapshots)/norm;nuclei=max(float(abs(s['nuclei']).max()) for s in self.snapshots)/norm
        energy=max(abs(s['energy']) for s in self.snapshots)/energy_norm;kinetic=max(abs(r['kinetic_residual_cm3']) for r in components)/max(float(self.initial.sum()),1.)
        self.data_issues.update(self.chemistry.gaps);self.model_notes.update(self.chemistry.notes)
        quality=dict(version=VERSION,model='coupled electrons, atomic/ionic levels and declared ground molecular carriers F_s(z,E,mu)',
            chemical_network_status='partial_registered_network; see pair_coverage.csv and unsupported.csv in source directories',
            chemical_channels=len(self.chemistry.reactions),chemical_rate_policy=self.cfg.get('chemistry',{}).get('rate_policy','cross_sections_only'),
            chemical_projection_momentum_error_g_cm_s=self.chemical_momentum_error,
            charge_status='passed' if charge<1e-9 else 'failed',charge_max_relative_residual=charge,nuclei_status='passed' if nuclei<1e-9 else 'failed',nuclei_max_relative_residual=nuclei,
            mass_status='passed' if nuclei<1e-9 else 'failed',mass_scope='same elemental nuclear-mass approximation for all ion stages',
            energy_status='passed' if energy<1e-6 else 'failed',energy_max_relative_residual=energy,energy_scope='electron + heavy kinetic + atomic internal; field work, photons, source injection, wall/outflow energy and declared line-profile reservoir',
            kinetic_status='passed' if kinetic<1e-9 else 'failed',kinetic_max_relative_residual=kinetic,
            energy_raw_max_relative_residual=max(abs(s['energy']-s['profile_heat']) for s in self.snapshots)/energy_norm,
            line_profile_reservoir_max_fraction=max(abs(s['profile_heat']) for s in self.snapshots)/energy_norm,
            reaction_photon_status='passed' if max(s['reaction_photon_residual'] for s in self.snapshots)<1e-8 else 'failed',
            reaction_photon_max_relative_residual=max(s['reaction_photon_residual'] for s in self.snapshots),
            atomic_data_status='incomplete' if self.data_issues else 'available_channels_checked',data_issues=sorted(self.data_issues),declared_data_models=sorted(self.model_notes),
            known_chemical_omissions_count=len(self.chemistry.omissions),wall_model='specular reflection, user-supplied probabilities, atomic ground conversion' if self.walls.enabled else 'off',
            spatial_particle_sources=sum(r['mode']!='radiation' for r in self.sources.items),spatial_photon_sources=sum(r['mode']=='radiation' for r in self.sources.items),
            active_states=int(self.selection.active.sum()),candidate_states=self.ns,basis_restarts=self.selection.restarts,
            basis_leakage_max_fraction=max((r['omitted_fraction_of_smallest_seed'] for r in self.selection.report()),default=0.),
            atomic_bank_completeness_status='finite_data_bank; omitted higher physical levels and unknown channels are not proven negligible',
            field_energy_status='not_checked' if self.cfg.get('electric_field',{}).get('mode')=='poisson' else 'prescribed field work accounted',
            accuracy_status='not_established; refine dt,E,mu,z,photon quadratures and atomic bank',actual_workers=min(self.workers,K.config.NUMBA_NUM_THREADS),
            elapsed_s=self.elapsed,accepted_steps=self.steps,rejected_steps=self.rejected,
            limitations=['Scalar LS, no density matrix/alignment/Zeeman physics.',
                'Inelastic electron impact uses stationary-target integral sigma(E), isotropic output and equal ionization sharing; heavy recoil/momentum in these channels omitted.',
                'RR uses both full velocity distributions; photon recoil and electron rest-mass changes of atomic mass omitted.',
                'ee/ii BGK and isotropic binary transport-equivalent e-ion operator; no full multicomponent Landau operator.',
                'Fixed-temperature CRD line profiles; moving-heavy Doppler shifts in transfer not implemented.',
                'Heavy chemistry uses supplied sigma(E), or explicitly selected k(T) moment closure. Isotropic products and Penning electron energy partition are approximations.',
                'Li-He ion-atom data are digitized primary curves (not CDS tables), with finite energy range and printed floor. He*-Li neutral Penning cross sections are explicitly labelled capture estimates.',
                'He*-Li+ radiative charge transfer/association and LiH molecular chemistry remain missing; consult chemical_known_omissions.csv.',
                'Walls use specified material probabilities and specular return; molecular wall dissociation and a self-consistent sheath/surface charging dynamics are not implemented.',
                'Ground molecular carriers have no resolved v,J spectroscopy; three-body/dielectronic recombination and nonthermal free-free remain absent.',
                'Level selection controls leakage into known omitted states; an incomplete database is not completed by the selector.'])
        self.chemistry.export(output)
        self.sources.export(output);self.walls.export(output)
        from spatial_diagnostics import export_spatial
        export_spatial(self,output)
        self.probe.export(output)
        write_json(output/'quality.json',quality);write_json(output/'input_config.json',self.cfg)
        write_json(output/'resolved_phase_grid.json',dict(energy_min_eV=float(self.ge.e[0]),energy_max_eV=float(self.ge.e[-1]),energy_points=len(self.ge.e),angular_points=len(self.ge.mu)))
        write_csv(output/'data_gaps.csv',[dict(issue=i) for i in sorted(self.data_issues)],['issue'])
        write_json(output/'atomic_bank_provenance.json',dict(imported=self.data.manifests,photo_sources=self.photo.sources,photo_source_receipts=self.photo.receipts,
            element_descriptors_sha256=getattr(self.data,'descriptor_hashes',{}),component_files_sha256=getattr(self.data,'component_files',{}),
            photo_missing_candidates=self.photo.gaps,bank_directory=str(self.root/self.cfg.get('atomic_data',{}).get('prepared_directory','atomic_runtime_bank'))))
        hashes={'src/'+path.name:hashlib.sha256(path.read_bytes()).hexdigest() for path in Path(__file__).resolve().parent.glob('*.py')}
        radiation_hashes={path.name:hashlib.sha256(path.read_bytes()).hexdigest() for path in (ROOT/'data/radiation').iterdir() if path.is_file()}
        write_json(output/'run_metadata.json',dict(version=VERSION,files_sha256=hashes,radiation_files_sha256=radiation_hashes,
            input_config_sha256=hashlib.sha256(json.dumps(self.cfg,sort_keys=True).encode()).hexdigest(),
            actual_species=sorted({s['species'] for s in self.states}),shared_electrons=True,shared_radiation=True))
        (output/'summary.txt').write_text(f"MixtureKineticSlab 2.1: {', '.join(self.data.elements)}\nЗаряд: {quality['charge_status']}; ядра: {quality['nuclei_status']}; энергия: {quality['energy_status']}; кинетика: {quality['kinetic_status']}.\n"+
            f"Атомные данные: {quality['atomic_data_status']}. Активных состояний {quality['active_states']}/{self.ns}; повторов базиса {self.selection.restarts}.\n"+
            'Читайте data_gaps.csv, level_selection.csv и docs/MODEL_RU.md. Балансы не доказывают полноту атомных данных/сходимость.\n',encoding='utf-8')
        self.progress.done(str(output));return quality


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--config',type=Path,default=ROOT/'configs/gas_probe.json');parser.add_argument('--output',type=Path);parser.add_argument('--workers',type=int);parser.add_argument('--progress',choices=('text','json','off'));parser.add_argument('--check',action='store_true');args=parser.parse_args(argv)
    cfg=json.loads(args.config.read_text(encoding='utf-8-sig'))
    if args.workers is not None:cfg.setdefault('performance',{})['workers']=args.workers
    progress=progress_from_config(cfg,args.progress)
    try:
        model=MixtureSlab(cfg,args.config.resolve().parent,progress)
        if args.check:progress.done('Смесь, распределения и выбранный базис построены');return 0
        model.solve();quality=model.export(args.output or args.config.resolve().parent/cfg.get('output_directory','output_mixture'))
        print('MixtureKineticSlab:',model.steps,'шагов;',round(model.elapsed,3),'с;',quality['atomic_data_status'],flush=True)
        return 2 if any(v=='failed' for k,v in quality.items() if k.endswith('_status')) or quality['atomic_data_status']=='incomplete' else 0
    except Exception as exc:progress.fail(str(exc));raise


if __name__=='__main__':
    try:sys.exit(main())
    except (ValueError,OSError,RuntimeError,KeyError) as exc:sys.exit('Ошибка: '+str(exc))
