"""Passive, weak-probe transmission through the computed populations.

The probe never contributes to rates, phase grids or population evolution.
Profiles use the same fixed-temperature Voigt/Einstein model as the solver.
Only resolved bound-bound lines are included; model wavelengths are vacuum.
"""
import numpy as np
from scipy.special import voigt_profile
from spectrum_input import load_spectrum
from atomic_database import write_csv,write_json,HC_EV_CM
from particle_kinetics import EV,KB

H=6.62607015e-27
C=2.99792458e10
AMU=1.66053906892e-24


class Probe:
    def __init__(self,model):
        self.model=model;cfg=model.cfg.get('probe',{});self.mode=cfg.get('mode','off')
        allowed={'mode','spectrum_file','spectrum_options','minimum_nm','maximum_nm','points','mu','path_cm'}
        if set(cfg)-allowed:raise ValueError('Неизвестные параметры пробного света: '+str(set(cfg)-allowed))
        if self.mode not in ('off','flat','file'):raise ValueError('Пробный свет: off, flat или file.')
        self.enabled=self.mode!='off';self.metadata={}
        if not self.enabled:return
        self.mu=float(cfg.get('mu',1.));self.path=float(cfg.get('path_cm',model.length))
        if not np.isfinite(self.mu+self.path) or not 0<self.mu<=1 or self.path<=0:raise ValueError('Пробный свет: 0<mu<=1, длина пути >0 см.')
        if self.mode=='file':
            data=load_spectrum(model.root/cfg.get('spectrum_file',''),cfg.get('spectrum_options'))
            original=np.asarray(data.wavelength_A);values=np.asarray(data.values);self.metadata=data.metadata
        else:
            low=float(cfg.get('minimum_nm',1082.));high=float(cfg.get('maximum_nm',1084.));points=int(cfg.get('points',501))
            if not np.isfinite(low+high) or not 0<low<high or not 2<=points<=200000:raise ValueError('Пробный свет: 0<λ_min<λ_max, число точек 2…200000.')
            original=np.linspace(low*10,high*10,points);values=np.ones(points)
            self.metadata=dict(input_x_unit='nm',input_y_kind='relative',points=points,warnings=[],purpose='probe')
        self.lines=[];parts=[original]
        for row in model.data.lines:
            if row.get('spectral_kind')=='two_photon' or not row.get('pump',True):continue
            i,j=model.index[row['lower']],model.index[row['upper']]
            nu0=(model.states[j]['energy_cm1']-model.states[i]['energy_cm1'])*HC_EV_CM*EV/H
            A=float(row['A_s'])
            if nu0<=0 or A<=0:continue
            sigma=nu0/C*np.sqrt(KB*model.Tn/(model.states[j]['mass_u']*AMU));gamma=A/(4*np.pi)
            center=C/nu0*1e8
            in_band=original[0]<=center<=original[-1]
            if in_band:
                offsets=np.r_[-np.geomspace(8,200,24)[::-1],np.linspace(-8,8,161),np.geomspace(8,200,24)]
                nu=nu0+max(sigma,gamma)*offsets;local=C/nu[nu>0]*1e8
                parts.append(local[(local>=original[0])&(local<=original[-1])])
            self.lines.append(dict(id=row['id'],i=i,j=j,nu0=nu0,sigma=sigma,gamma=gamma,A=A,center=center,in_band=in_band,source=row.get('source','')))
        self.wavelength=np.unique(np.concatenate(parts));self.incident=np.interp(self.wavelength,original,values)
        nu=C/(self.wavelength*1e-8)
        for line in self.lines:
            line['cross']=line['A']*(C/line['nu0'])**2/(8*np.pi)*voigt_profile(nu-line['nu0'],line['sigma'],line['gamma'])
            line['ratio']=model.states[line['j']]['g']/model.states[line['i']]['g']
        self.metadata.update(mode=self.mode,path_cm=self.path,mu=self.mu,resolved_points=len(self.wavelength),
            wavelength_basis='vacuum, imported LS energies; not fine-structure experimental wavelengths',
            scope='passive weak probe; bound-bound absorption minus stimulated emission; no continuum absorption or added spontaneous light',
            population_feedback=False)

    def evaluate(self,populations):
        """Return cumulative center/exit optical depths and exact beam ratios."""
        kappa=np.zeros((self.model.cells,len(self.wavelength)))
        for line in self.lines:
            net=populations[:,line['i']]*line['ratio']-populations[:,line['j']]
            kappa+=net[:,None]*line['cross'][None,:]
        cell_tau=kappa*(self.path/self.model.cells/self.mu)
        tau=cell_tau.sum(axis=0)
        if np.any(tau<-700):raise ValueError('Пробный свет: усиление exp(-τ) превышает диапазон чисел; приближение слабого зонда неприменимо.')
        return tau,np.exp(-tau),np.cumsum(cell_tau,axis=0)-.5*cell_tau

    def export(self,folder):
        if not self.enabled:return
        m=self.model;history=[];position=[];ranking=[]
        columns=['time_s','wavelength_nm','incident','transmitted','transmission','absorption_fraction','optical_depth']
        final=[]
        for snap in m.snapshots:
            tau,ratio,centers=self.evaluate(snap['n']);rows=[]
            for k,lam in enumerate(self.wavelength):
                rows.append(dict(time_s=snap['time_s'],wavelength_nm=float(lam/10),incident=float(self.incident[k]),transmitted=float(self.incident[k]*ratio[k]),transmission=float(ratio[k]),absorption_fraction=float(1-ratio[k]),optical_depth=float(tau[k])))
            history.extend(rows);final=rows
        write_csv(folder/'probe_transmission.csv',final,columns)
        write_csv(folder/'probe_history.csv',history,columns)
        for z in range(m.cells):
            for k,lam in enumerate(self.wavelength):
                position.append(dict(x_cm=(z+.5)*m.dz,wavelength_nm=float(lam/10),optical_depth=float(centers[z,k]),transmission=float(np.exp(-centers[z,k]))))
        write_csv(folder/'probe_by_position_final.csv',position,['x_cm','wavelength_nm','optical_depth','transmission'])
        n=m.snapshots[-1]['n']
        for line in self.lines:
            if not line['in_band']:continue
            net=(n[:,line['i']]*line['ratio']-n[:,line['j']]).mean();depth=net*self.path/self.mu*line['cross']
            ranking.append(dict(id=line['id'],wavelength_nm=line['center']/10,peak_optical_depth=float(depth.max()),equivalent_width_nm=float(np.trapezoid(-np.expm1(-depth),self.wavelength/10)),source=line['source']))
        ranking.sort(key=lambda row:-abs(row['equivalent_width_nm']))
        write_csv(folder/'probe_lines.csv',ranking,['id','wavelength_nm','peak_optical_depth','equivalent_width_nm','source'])
        write_json(folder/'probe_input_metadata.json',self.metadata)


def plot_probe(folder):
    import csv
    import matplotlib.pyplot as plt
    path=folder/'probe_transmission.csv'
    if not path.exists():return
    with path.open(encoding='utf-8-sig',newline='') as stream:rows=list(csv.DictReader(stream,delimiter=';'))
    if not rows:return
    x=np.array([float(r['wavelength_nm']) for r in rows]);fig,axes=plt.subplots(2,1,figsize=(10,6.5),sharex=True)
    for key,label in [('incident','Вход'),('transmitted','Выход')]:axes[0].plot(x,[float(r[key]) for r in rows],label=label,lw=1.2)
    axes[0].set(ylabel='Интенсивность в выбранных единицах',title='Пробный свет · слабый зонд, без влияния на кинетику');axes[0].legend()
    axes[1].plot(x,[100*float(r['absorption_fraction']) for r in rows]);axes[1].set(xlabel='Длина волны в вакууме, нм',ylabel='Поглощение, %')
    for ax in axes:ax.grid(alpha=.25)
    fig.tight_layout();fig.savefig(folder/'probe_transmission.png',dpi=160);plt.close(fig)
