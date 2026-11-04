"""Positive energy/angle particle operators, cgs and eV; no atomic data here.

F[E,mu] is number per cm3 in a quadrature node, NOT a phase-space PDF.
Axisymmetry about z retains v_perp, v_z and two propagation directions.
"""
import numpy as np
from scipy.linalg import solve_banded
from scipy.special import logsumexp
from numpy.polynomial.legendre import Legendre
from pathlib import Path
import csv
from radiation import weights
from grid_validation import validate_phase_grid

ME = 9.1093837139e-28
MI = 4.002602*1.66053906892e-24
EV = 1.602176634e-12
KB = 1.380649e-16
QE = 1.602176634e-19  # E is V/cm, hence force = q*E*1e7 dyn


class GridOverflow(ValueError):
    pass


class PhaseGrid:
    def __init__(self, config, mass_g=ME, label='phase_grid'):
        self.mass=mass_g
        validate_phase_grid(config, label)
        lo, hi = config.get('energy_min_eV',1e-4), config.get('energy_max_eV',50.)
        ne, na = config.get('energy_points',100), config.get('angular_points',8)
        self.e = np.geomspace(lo, hi, ne)
        # Gauss-Lobatto includes beams mu=+-1 with correct small endpoint weights.
        poly=Legendre.basis(na-1)
        self.mu=np.r_[-1.,poly.deriv().roots(),1.]
        self.we=weights(self.e);self.wmu=1/(na*(na-1)*poly(self.mu)**2)
        self.v = np.sqrt(2*self.e*EV/self.mass)
        self.vz = self.v[:,None]*self.mu[None,:]
        self.vperp = self.v[:,None]*np.sqrt(1-self.mu[None,:]**2)
        self.dos = np.sqrt(self.e)*self.we
        self.base = self.dos[:,None]*self.wmu
        self.shape = (ne,na)

    def initial(self, config, density, root=Path('.')):
        allowed = {'Te_eV','hot_fraction','hot_Te_eV','angular_bias','density_cm3','distribution_file','kind','energy_distribution_file','energy_form','match_moments'}
        if set(config)-allowed: raise ValueError('Неизвестные параметры electron_initial')
        from energy_input import distribution_kind,energy_phase
        kind=distribution_kind(config)
        if kind=='phase_table':return self.from_csv(Path(root)/config['distribution_file'],density)
        if kind=='energy_table':return energy_phase(self,Path(root)/config['energy_distribution_file'],density,config.get('energy_form','probability'))
        T, hot, Th = config.get('Te_eV',1.), config.get('hot_fraction',0.), config.get('hot_Te_eV',5.)
        bias = config.get('angular_bias',0.)
        if not np.all(np.isfinite([T,hot,Th,bias,density])) or min(T,Th)<=0 or not 0<=hot<=1 or abs(bias)>1 or density<0:
            raise ValueError('Неверное начальное распределение электронов')
        def component(temp):
            raw=self.dos*np.exp(-self.e/temp)
            # Finite domain normalization is explicit; tail reported separately.
            if self.e[-1]/temp < 8:
                raise GridOverflow('Emax < 8 kTe: начальный хвост не разрешён; увеличьте сетку или уменьшите Te')
            return raw/raw.sum()
        if config.get('match_moments',False) and bias==0:
            if self.e[-1]<8*max(T,Th if hot else T):raise GridOverflow('Emax < 8 kTe: увеличьте сетку.')
            result=(1-hot)*self.equilibrium_moments(density,1.5*T,0.)
            if hot:result+=hot*self.equilibrium_moments(density,1.5*Th,0.)
            return result
        fe=(1-hot)*component(T)+hot*component(Th)
        return density*fe[:,None]*self.wmu*(1+bias*self.mu)

    def ion_initial(self, config, density, root=Path('.')):
        allowed={'uz_cm_s','Ti_K','beams','distribution_file','density_cm3','match_moments','kind','energy_distribution_file','energy_form'}
        if set(config)-allowed:raise ValueError('Неизвестные ion_initial: '+str(set(config)-allowed))
        from energy_input import distribution_kind,energy_phase
        kind=distribution_kind(config)
        if kind=='phase_table':return self.from_csv(Path(root)/config['distribution_file'],density)
        if kind=='energy_table':
            if config.get('uz_cm_s',0.)!=0:raise ValueError('Файл P(E) изотропен и имеет u_z=0. Для направленного пучка выберите максвелловское распределение или фазовый CSV E, μ.')
            return energy_phase(self,Path(root)/config['energy_distribution_file'],density,config.get('energy_form','probability'))
        beams=config.get('beams',[dict(fraction=1.,uz_cm_s=config.get('uz_cm_s',0.),Ti_K=config.get('Ti_K',300.))])
        if not beams or abs(sum(b.get('fraction',0.) for b in beams)-1)>1e-10:raise ValueError('Сумма ion_initial.beams.fraction должна быть 1')
        out=np.zeros(self.shape)
        for beam in beams:
            if set(beam)-{'fraction','uz_cm_s','Ti_K'}:raise ValueError('Неизвестные параметры ионного пучка')
            w=float(beam['fraction']);u=float(beam.get('uz_cm_s',0.));T=float(beam.get('Ti_K',300.))*KB/EV
            if not np.all(np.isfinite([w,u,T,density])) or w<0 or T<=0 or density<0:raise ValueError('Неверный ионный пучок')
            bulk=self.mass*u*u/(2*EV)
            if self.e[-1]<bulk+8*T:
                motion=f', u_z={u:.7g} см/с' if u else ''
                raise GridOverflow(f'Сетка тяжёлых частиц: Emax={self.e[-1]:.7g} эВ, требуется не менее {bulk+8*T:.7g} эВ при Ti={T*EV/KB:.7g} K (kTi={T:.7g} эВ){motion}. В «EEDF / ионы → Тяжёлые частицы → Emax» задайте auto или увеличьте численную границу. Температура не изменена.')
            logp=np.log(self.base)-(self.e[:,None]+bulk-self.mass*u*self.vz/EV)/T
            p=np.exp(logp-logsumexp(logp))
            if config.get('match_moments',False):p=self.equilibrium_moments(1.,bulk+1.5*T,u)
            out+=density*w*p
        return out

    def from_csv(self,path,density):
        with Path(path).open(encoding='utf-8-sig',newline='') as f:rows=list(csv.DictReader(f,delimiter=';'))
        if not rows or set(rows[0])!={'energy_eV','mu','probability_per_eV_mu'}:raise ValueError('Фазовый CSV: energy_eV;mu;probability_per_eV_mu')
        ee=np.unique([float(r['energy_eV']) for r in rows]);mm=np.unique([float(r['mu']) for r in rows])
        if len(ee)<2 or len(mm)<2 or ee[0]<0 or mm[0]<-1 or mm[-1]>1 or len(rows)!=len(ee)*len(mm):raise ValueError('Нужна прямоугольная сетка energy,mu без пропусков/дубликатов')
        raw=np.zeros((len(ee),len(mm)));seen=set()
        for row in rows:
            key=(float(row['energy_eV']),float(row['mu']))
            if key in seen:raise ValueError('Дубликат узла фазового CSV')
            seen.add(key);raw[np.searchsorted(ee,key[0]),np.searchsorted(mm,key[1])]=float(row['probability_per_eV_mu'])
        if np.any(raw<0) or not np.all(np.isfinite(raw)):raise ValueError('Фазовый CSV требует конечной неотрицательной PDF')
        integral=float(weights(ee)@raw@weights(mm))
        if abs(integral-1)>1e-3:raise ValueError('Интеграл входного фазового CSV должен быть 1')
        if ee[0]>self.e[0] or ee[-1]<self.e[-1] or mm[0]>-1 or mm[-1]<1:raise ValueError('Фазовый CSV должен покрывать сетку расчёта')
        tab=np.array([np.interp(self.mu,mm,row) for row in raw])
        p=np.array([np.interp(self.e,ee,tab[:,j]) for j in range(len(self.mu))]).T
        p=p*self.we[:,None]*(2*self.wmu)
        if abs(p.sum()-1)>.02:raise ValueError('Сетка расчёта недостаточно разрешает входную PDF')
        # Initial-condition quadrature normalization only. Never used on a solution.
        return density*p/p.sum()

    def moments(self, f):
        ne=float(f.sum()); en=float((f*self.e[:,None]).sum()); pz=float((f*self.vz).sum())
        uz=pz/ne if ne else 0.
        thermal=en-self.mass*ne*uz**2/(2*EV)
        ppar=self.mass*float((f*(self.vz-uz)**2).sum())
        pperp=self.mass/2*float((f*self.vperp**2).sum())
        return dict(ne_cm3=ne,energy_eV_cm3=en,uz_cm_s=uz,P_parallel_erg_cm3=ppar,P_perp_erg_cm3=pperp,
                    T_moment_eV=2*thermal/(3*ne) if ne else 0.,momentum_g_cm2_s=self.mass*pz,
                    forward_fraction=float(f[:,self.mu>0].sum()/ne) if ne else 0.,
                    tail_fraction=float(f[-2:].sum()/ne) if ne else 0.)

    def deposit(self, amount, energies, mus=None, allow_tail=False):
        """Linear number/energy conserving projection within [Emin,Emax].

        Tiny sub-Emin energies go to Emin; the returned energy error is a
        diagnostic, not a reservoir used to make an energy audit pass.
        Above-Emax significant contributions stop the calculation.
        """
        a,e=np.broadcast_arrays(np.asarray(amount,dtype=float),np.asarray(energies,dtype=float))
        if np.any(a < 0) or np.any(~np.isfinite(a)) or np.any(~np.isfinite(e)):
            raise ValueError('Неверный источник фазового распределения')
        if not allow_tail and np.any((e>self.e[-1]) & (a>max(float(a.sum())*1e-13,1e-100))):
            raise GridOverflow('Продукты реакции/ускорение вышли за Emax; нужна более широкая сетка/He++')
        clipped=np.clip(e,self.e[0],self.e[-1]); err=float((a*(clipped-e)).sum())
        j=np.clip(np.searchsorted(self.e,clipped,side='right')-1,0,len(self.e)-2)
        q=(clipped-self.e[j])/(self.e[j+1]-self.e[j])
        out=np.zeros(self.shape)
        if mus is None:
            v=np.bincount(j.ravel(),weights=(a*(1-q)).ravel(),minlength=len(self.e))
            v+=np.bincount((j+1).ravel(),weights=(a*q).ravel(),minlength=len(self.e))
            out[:]=v[:,None]*self.wmu
        else:
            mu=np.broadcast_to(np.asarray(mus),a.shape)
            vz=np.sqrt(2*np.maximum(e,0.)*EV/self.mass)*mu
            vlo=self.v[j];vhi=self.v[j+1]
            mlo=np.clip(vz/vlo,-1,1)
            mhi=np.divide(vz-(1-q)*vlo*mlo,q*vhi,out=np.zeros_like(vz),where=q>0)
            mhi=np.clip(mhi,-1,1)
            for dj,ve,mm in ((0,1-q,mlo),(1,q,mhi)):
                k=np.clip(np.searchsorted(self.mu,mm,side='right')-1,0,len(self.mu)-2)
                rr=np.clip((mm-self.mu[k])/(self.mu[k+1]-self.mu[k]),0,1)
                for dk,vm in ((0,1-rr),(1,rr)):
                    np.add.at(out,(j.ravel()+dj,k.ravel()+dk),(a*ve*vm).ravel())
        return out,err

    def kick(self, f, dv):
        if dv==0: return f.copy(),0.,0.
        z=self.vz+dv; en=self.mass*(z*z+self.vperp*self.vperp)/(2*EV)
        speed=np.sqrt(z*z+self.vperp*self.vperp)
        mu=np.divide(z,speed,out=np.zeros_like(z),where=speed>0)
        out,err=self.deposit(f,en,mu)
        expected=self.mass*float((f*(self.vz+dv)).sum())
        defect=self.moments(out)['momentum_g_cm2_s']-expected
        return out,err,defect

    def maximum_entropy(self, f):
        """Discrete drifting equilibrium with exactly the same N,E,pz (BGK)."""
        n=float(f.sum())
        if n==0:return f.copy()
        return self.equilibrium_moments(n,float((f*self.e[:,None]).sum()/n),float((f*self.vz).sum()/n))

    def equilibrium_moments(self,n,mean_energy_eV,uz_cm_s):
        """Existing entropy quadrature with prescribed number, energy and speed."""
        if n==0:return np.zeros(self.shape)
        escale=max(float(mean_energy_eV),self.e[0])
        vscale=np.sqrt(2*escale*EV/self.mass)
        x=np.broadcast_to(self.e[:,None]/escale,self.shape).ravel()
        y=(self.vz/vscale).ravel();target=np.array([mean_energy_eV/escale,uz_cm_s/vscale])
        feat=np.array([x,y]); logbase=np.log(self.base.ravel())
        T=max(2*(target[0]-target[1]**2)/3,1e-8)
        theta=np.array([-1/T,2*target[1]/T])
        def evaluate(p):
            logw=logbase+p@feat; prob=np.exp(logw-logsumexp(logw))
            mean=feat@prob; d=feat-mean[:,None]
            return prob,mean-target,(d*prob)@d.T
        # A narrow injected beam can make the continuous-Maxwellian guess
        # nearly singular on a coarse phase grid. Use a second, interior start
        # and the convex entropy dual, rather than an unconstrained moment fix.
        def dual(p):return float(logsumexp(logbase+p@feat)-p@target)
        for start in (theta,np.zeros(2)):
            theta=start.copy()
            for _ in range(180):
                prob,res,cov=evaluate(theta)
                if np.max(abs(res))<2e-12:return (n*prob).reshape(self.shape)
                step=np.linalg.lstsq(cov,res,rcond=1e-14)[0]
                if not np.all(np.isfinite(step)) or np.dot(res,step)<=0:break
                magnitude=float(np.linalg.norm(step))
                if magnitude>100:step*=100/magnitude
                loss=dual(theta);slope=float(res@step)
                for power in range(35):
                    factor=2.**(-power);trial=theta-factor*step
                    trialres=evaluate(trial)[1]
                    if dual(trial)<loss-1e-4*factor*slope or np.linalg.norm(trialres)<np.linalg.norm(res)*(1-1e-4*factor):theta=trial;break
                else:break
        raise ValueError('Не удалось согласовать концентрацию, энергию и скорость на сетке; увеличьте energy_points/angular_points или диапазон энергии.')

    def ee_bgk(self, f, dt, coulomb_log):
        m=self.moments(f)
        if m['ne_cm3']==0:return f.copy(),0.
        T=m['T_moment_eV']
        if T<=0:raise ValueError('Неположительная дисперсия электронных скоростей')
        # NRL thermal collision-time scale, used as an explicitly approximate
        # BGK relaxation rate; this is NOT a Landau collision operator.
        nu=2.91e-6*m['ne_cm3']*coulomb_log/T**1.5
        a=-np.expm1(-nu*dt)
        return (1-a)*f+a*self.maximum_entropy(f),nu

    def neutral_elastic(self, f, dt, rate, Tn):
        """Isotropic integral-XS closure + leading-mass-ratio energy diffusion.

        Scharfetter-Gummel detailed balance: stationary pdf sqrt(E)exp(-E/Tn).
        rate=N_He sigma_total v, valid as sigma_mt only under stated isotropy.
        Energy lost by electrons is recorded in the fixed neutral reservoir.
        """
        iso=f.sum(axis=1)[:,None]*self.wmu
        out=iso+(f-iso)*np.exp(-rate[:,None]*dt)
        temp=KB*Tn/EV; ef=(self.e[:-1]+self.e[1:])/2
        D=2*ME/MI*(rate[:-1]+rate[1:])/2*temp*ef
        de=np.diff(self.e); logM=.5*np.log(self.e)-self.e/temp
        x=-(np.diff(logM))
        def bernoulli(v):
            ans=np.empty_like(v);small=abs(v)<1e-4;pos=v>50;neg=v<-50
            ans[small]=1-v[small]/2+v[small]**2/12
            ans[pos]=v[pos]*np.exp(-v[pos]);ans[neg]=-v[neg]
            mid=~(small|pos|neg);ans[mid]=v[mid]/np.expm1(v[mid])
            return ans
        up=D/de/self.we[:-1]*bernoulli(x)
        down=D/de/self.we[1:]*bernoulli(-x)
        ab=np.zeros((3,len(self.e)))
        ab[1]=1+dt*(np.r_[up,0.]+np.r_[0.,down])
        ab[0,1:]=-dt*down;ab[2,:-1]=-dt*up
        return solve_banded((1,1),ab,out)


def stream_periodic(f, vz, dt, dz):
    """Conservative translation with minmod linear reconstruction, any CFL.

    It is positive and conservative. Limited reconstruction reduces spreading
    of ballistic packets; cell/timestep refinement is still needed.
    """
    s=vz*dt/dz; integer=np.floor(s).astype(np.int64); a=s-integer
    cells=len(f); z=np.arange(cells)[:,None,None]
    j=np.arange(f.shape[1])[None,:,None]; k=np.arange(f.shape[2])[None,None,:]
    left=f-np.roll(f,1,axis=0);right=np.roll(f,-1,axis=0)-f
    slope=np.where(left*right>0,np.sign(left)*np.minimum(abs(left),abs(right)),0.)
    current=(z-integer)%cells;previous=(z-integer-1)%cells
    return (1-a)*f[current,j,k]+a*f[previous,j,k]+.5*a*(1-a)*(slope[previous,j,k]-slope[current,j,k])


def stream_open(f, vz, dt, dz, left, right):
    """Finite-volume characteristic integrals with fixed incoming reservoirs.

    Returns net *outward* column number per phase node (cm^-2), including
    incoming boundary particles; never restores quasineutrality algebraically.
    """
    cells=len(f); out=np.empty_like(f)
    for j in range(f.shape[1]):
        for k in range(f.shape[2]):
            speed=vz[j,k]; shift=speed*dt/dz
            vals=f[:,j,k]; prefix=np.r_[0.,np.cumsum(vals)]
            def primitive(x):
                if x<0:return x*left[j,k]
                if x>=cells:return prefix[-1]+(x-cells)*right[j,k]
                i=int(np.floor(x));return prefix[i]+(x-i)*vals[i]
            for z in range(cells):
                out[z,j,k]=primitive(z+1-shift)-primitive(z-shift)
    return out,dz*(f.sum(axis=0)-out.sum(axis=0))

