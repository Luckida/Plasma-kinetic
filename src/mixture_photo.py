"""Photoionization data for the generic mixture solver.

H/hydrogenic one-electron ions: nonrelativistic Coulomb dipole integrals,
not a Kramers substitution. He I: HS98 tables. Other ionic ground states:
Verner et al. 1996. Unknown excited states remain explicitly missing.
"""
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import re
import numpy as np
import mpmath as mp
from scipy import constants as C
from atomic_core import SYMBOLS
from atomic_database import SourceCache, HC_EV_CM
SOURCE_DOCUMENT=json.loads((Path(__file__).resolve().parent.parent/'data/radiation/sources.json').read_text(encoding='utf-8-sig'))
SOURCE_SHA256={r['file']:r['sha256'] for r in SOURCE_DOCUMENT['sources']}

EV=1.602176634e-12
HARTREE_EV=C.physical_constants['Hartree energy in eV'][0]
BOHR_CM2=(C.physical_constants['Bohr radius'][0]*100)**2
VERNER=SOURCE_DOCUMENT['verner']


@lru_cache(maxsize=32768)
def hydrogen_sigma(n,l,excess_eV,Z=1):
    """LS-averaged sigma in cm2; continuum states normalized per Hartree.

    R_bound uses the normalized Laguerre polynomial. The radial integral of
    the regular Coulomb wave (DLMF 33.2.4-5) is evaluated as a finite Laguerre
    sum of 2F1 functions. H threshold -> 6.3043e-18 cm2. No density matrices.
    """
    if not isinstance(n,int) or not isinstance(l,int) or not 0<=l<n or Z<1 or excess_eV<0:
        raise ValueError('Неверные квантовые числа/энергия hydrogenic PI')
    with mp.workdps(35):
        energy=mp.mpf(max(excess_eV,1e-8))/(HARTREE_EV*Z*Z)
        k=mp.sqrt(2*energy);eta=-1/k;lam=mp.mpf(1)/n+1j*k
        bound=mp.mpf(2)/n**2*mp.sqrt(mp.factorial(n-l-1)/mp.factorial(n+l))
        strength=mp.mpf(0)
        for lp,weight in [(l-1,l),(l+1,l+1)]:
            if lp<0:continue
            norm=bound*mp.sqrt(2/(mp.pi*k))*mp.coulombc(lp,eta)*k**(lp+1)*(mp.mpf(2)/n)**l
            value=mp.mpc(0)
            for j in range(n-l):
                p=l+lp+j+4
                coefficient=(-1)**j*mp.binomial(n+l,n-l-1-j)/mp.factorial(j)*(mp.mpf(2)/n)**j
                value+=coefficient*mp.gamma(p)/lam**p*mp.hyp2f1(lp+1-1j*eta,p,2*lp+2,2j*k/lam)
            strength+=weight*abs(norm*value)**2
        sigma=4*mp.pi**2*C.alpha*(energy+mp.mpf(1)/(2*n*n))*strength/(3*(2*l+1))*BOHR_CM2/Z**2
        return float(sigma)


def verner_sigma(parameters,energy):
    eth,emax,e0,s0,ya,p,yw,y0,y1=parameters
    energy=np.asarray(energy,dtype=float)
    if np.any(energy>emax):raise ValueError('Выход за область Verner PI fit')
    x=energy/e0-y0;y=np.sqrt(x*x+y1*y1)
    value=s0*((x-1)**2+yw*yw)*y**(.5*p-5.5)*(1+np.sqrt(y/ya))**(-p)*1e-18
    return np.where(energy>=eth,value,0.)


class PhotoBank:
    def __init__(self,root,states,parent_index,maximum_eV,offline=False,progress=None,overrides=None):
        self.root=Path(root);self.states=states;self.parents=parent_index
        self.tables={};self.sources={};self.gaps=[]
        bundled=Path(__file__).resolve().parent.parent/'data/radiation/photoionization_sources'
        cache=SourceCache(self.root/'.atomic_photo_sources',offline,progress,bundled_directory=bundled)
        generated=self.root/'.atomic_photo';generated.mkdir(exist_ok=True)
        he_tables={};verner=None
        if overrides:
            import csv
            with (self.root/overrides).open(encoding='utf-8-sig',newline='') as f:
                records=list(csv.DictReader(f,delimiter=';'))
            for row in records:
                self.tables.setdefault(row['level'],[]).append((float(row['photon_eV']),float(row['sigma_cm2'])))
            for ident,table in self.tables.items():
                self.tables[ident]=np.array(sorted(table));self.sources[ident]='user CSV: '+overrides
        for j,state in enumerate(states):
            if j not in parent_index:continue
            ident=state['id'];parent=states[parent_index[j]]
            threshold=(parent['energy_cm1']-state['energy_cm1'])*HC_EV_CM
            if threshold<=0:raise ValueError('Неположительный порог фотоионизации')
            if ident in self.tables:
                tab=self.tables[ident]
                if len(tab)<2 or np.any(np.diff(tab[:,0])<=0) or np.any(tab[:,1]<0) or tab[0,0]>threshold*(1+1e-6):raise ValueError('Неверная пользовательская PI-таблица '+ident)
                continue
            Z=SYMBOLS.index(state['element'])+1;electrons=Z-state['charge']
            if electrons==1:
                n,l=int(state['n']),int(state['L'])
                excess=np.r_[0.,np.geomspace(1e-6,max(maximum_eV-threshold,1e-5),160)]
                code_hash=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
                key=hashlib.sha256((str((n,l,Z,threshold,maximum_eV))+'Coulomb-1.0'+code_hash).encode()).hexdigest()[:24]
                file=generated/(key+'.npz')
                if file.exists():
                    with np.load(file) as data:tab=data['table']
                else:
                    if progress:progress.module('Кулоновские PI-сечения H-like',state['label'],force=True)
                    tab=np.column_stack((threshold+excess,[hydrogen_sigma(n,l,float(e),Z) for e in excess]))
                    np.savez_compressed(file,table=tab)
                self.tables[ident]=tab;self.sources[ident]='NR Coulomb dipole radial integral; DLMF 33.2; threshold aligned preserving excess energy'
            elif state['element']=='He' and state['charge']==0:
                suffix=str(state['multiplicity'])+'spdfghi'[state['L']]+('e' if state['L']%2==0 else 'o')
                name='p0202.'+suffix
                if name not in he_tables:
                    info=dict(url='https://gitlab.nublado.org/cloudy/cloudy/-/raw/master/data/sh98_he1/pi/'+name,
                              reference='Hummer & Storey 1998; Cloudy data, license in radiation_data')
                    if name in SOURCE_SHA256:info['sha256']=SOURCE_SHA256[name]
                    text=cache.get(name,info).decode();lines=text.splitlines();pos=0;found={}
                    while pos<len(lines):
                        if lines[pos].split()==['0','0','0']:break
                        mult,angular,_,number=map(int,lines[pos+3].split())
                        nn=angular+number+int(angular==0 and mult==3)
                        count=int(lines[pos+4].split()[1]);tab=np.array([list(map(float,r.split())) for r in lines[pos+6:pos+6+count]])
                        found[nn]=tab;pos+=6+count
                    he_tables[name]=found
                if state['n'] not in he_tables[name]:
                    self.gaps.append(ident+': HS98 level absent');continue
                tab=he_tables[name][state['n']].copy();ryd=C.physical_constants['Rydberg constant times hc in eV'][0]
                tab[:,0]=threshold+(tab[:,0]-tab[0,0])*ryd;tab[:,1]*=1e-18
                if maximum_eV>tab[-1,0]:
                    tail=np.geomspace(np.nextafter(tab[-1,0],np.inf),maximum_eV,80)
                    tab=np.vstack((tab,np.column_stack((tail,tab[-1,1]*(tab[-1,0]/tail)**3))))
                self.tables[ident]=tab;self.sources[ident]='HS98; '+name+'; threshold aligned; matched E^-3 high-energy tail as in legacy He radiation'
            elif state.get('is_ground'):
                if verner is None:
                    raw=cache.get('verner_photo.dat',VERNER)
                    verner={(int(r[0]),int(r[1])):np.array(r[2:],dtype=float) for r in (line.split() for line in raw.decode().splitlines()) if r}
                parameters=verner.get((Z,electrons))
                if parameters is not None:
                    high=min(maximum_eV,parameters[1]-parameters[0]+threshold)
                    energy=np.geomspace(threshold,high,220)
                    self.tables[ident]=np.column_stack((energy,verner_sigma(parameters,energy-threshold+parameters[0])))
                    self.sources[ident]='Verner et al. 1996 ground-state fit; threshold aligned'
                else:self.gaps.append(ident+': Verner ground absent')
            else:self.gaps.append(ident+': no excited-state PI/RR cross section')
        self.receipts=cache.receipts

    def sigma(self,ident,energy):
        table=self.tables[ident];energy=np.asarray(energy)
        if np.any(energy>table[-1,0]*(1+1e-12)):
            raise ValueError('PI/RR энергия вне исходной таблицы '+ident)
        return np.interp(energy,table[:,0],table[:,1],left=0.,right=table[-1,1])
