"""CPU kernels for kinetic e/He+ distributions. No Euler equations.

Parallel axes have private outputs, followed by a deterministic sum. No
fastmath, random sampling, shared histogram writes or tolerance changes.
"""
from contextlib import contextmanager
import math
import numpy as np
from numba import njit,prange,get_num_threads,set_num_threads,config
from particle_kinetics import ME,MI,EV


@contextmanager
def kernel_threads(workers):
    before=get_num_threads();wanted=min(workers,config.NUMBA_NUM_THREADS)
    if before==wanted:
        yield
    else:
        set_num_threads(wanted)
        try:yield
        finally:set_num_threads(before)


@njit(cache=True)
def bracket(x,value):
    i=np.searchsorted(x,value,side='right')-1
    i=max(0,min(len(x)-2,i))
    return i,min(1.,max(0.,(value-x[i])/(x[i+1]-x[i])))


@njit(cache=True)
def put_phase(out,amount,energy,vz,ee,mu,mass):
    """Positive number/E projection, exact pz inside the nodal convex hull.

    In the small polar cone outside that hull, pz has a reported projection
    defect. No correction of populations is used to hide this defect.
    """
    k,a=bracket(ee,energy)
    vlo=math.sqrt(2*ee[k]*EV/mass);vhi=math.sqrt(2*ee[k+1]*EV/mass)
    mlo=min(1.,max(-1.,vz/vlo))
    mhi=min(1.,max(-1.,(vz-(1-a)*vlo*mlo)/(a*vhi))) if a>0 else 0.
    j,b=bracket(mu,mlo);l,c=bracket(mu,mhi);na=len(mu)
    out[k*na+j]+=amount*(1-a)*(1-b);out[k*na+j+1]+=amount*(1-a)*b
    out[(k+1)*na+l]+=amount*a*(1-c);out[(k+1)*na+l+1]+=amount*a*c


@njit(cache=True)
def phase_projection(energy,vz,ee,mu,mass):
    k,a=bracket(ee,energy);na=len(mu)
    vlo=math.sqrt(2*ee[k]*EV/mass);vhi=math.sqrt(2*ee[k+1]*EV/mass)
    lo=min(1.,max(-1.,vz/vlo))
    hi=min(1.,max(-1.,(vz-(1-a)*vlo*lo)/(a*vhi))) if a>0 else 0.
    j,b=bracket(mu,lo);l,c=bracket(mu,hi)
    return (k*na+j,k*na+j+1,(k+1)*na+l,(k+1)*na+l+1),((1-a)*(1-b),(1-a)*b,a*(1-c),a*c)


@njit(parallel=True,cache=True)
def scatter_map(rel,cx,cy,cz,p,ee,emu,ie,imu,scatter_mu,scatter_weight,scatter_phi,mass_i=MI):
    """Geometry only: no density, distribution, screening or rate is cached."""
    ne,ni,_=rel.shape;columns=len(scatter_mu)*len(scatter_phi)*4
    erows=np.empty((ne*ni,columns),dtype=np.int32);irows=np.empty_like(erows)
    ew=np.empty((ne*ni,columns));iw=np.empty_like(ew)
    reduced=ME*mass_i/(ME+mass_i);mass=ME+mass_i
    for a in prange(ne):
        for b in range(ni):
            pair=a*ni+b;speed=math.sqrt(2*rel[a,b,p]*EV/reduced);slot=0
            for t in range(len(scatter_mu)):
                z=speed*scatter_mu[t];perp=speed*math.sqrt(1-scatter_mu[t]**2)
                for phi in scatter_phi:
                    x=perp*math.cos(phi);y=perp*math.sin(phi)
                    ex=cx[a,b,p]+mass_i/mass*x;ey=cy[a,b,p]+mass_i/mass*y;ez=cz[a,b,p]+mass_i/mass*z
                    ix=cx[a,b,p]-ME/mass*x;iy=cy[a,b,p]-ME/mass*y;iz=cz[a,b,p]-ME/mass*z
                    ep,wp=phase_projection(ME*(ex*ex+ey*ey+ez*ez)/(2*EV),ez,ee,emu,ME)
                    ip,vi=phase_projection(mass_i*(ix*ix+iy*iy+iz*iz)/(2*EV),iz,ie,imu,mass_i)
                    weight=scatter_weight[t]/len(scatter_phi)
                    for k in range(4):
                        erows[pair,slot]=a*ne+ep[k];irows[pair,slot]=a*ni+ip[k]
                        ew[pair,slot]=wp[k]*weight;iw[pair,slot]=vi[k]*weight;slot+=1
    return erows,ew,irows,iw


@njit(parallel=True,cache=True)
def pair_rates(rel,screen,logcap,mass_i=MI,ion_charge=1.):
    out=np.empty_like(rel)
    for a in prange(rel.shape[0]):
        for b in range(rel.shape[1]):
            for p in range(rel.shape[2]):out[a,b,p]=momentum_transfer_rate(rel[a,b,p],screen,logcap,mass_i,ion_charge)
    return out


@njit(parallel=True,cache=True)
def sparse_product(data,indices,indptr,vector):
    out=np.empty(len(indptr)-1)
    for row in prange(len(out)):
        value=0.
        for k in range(indptr[row],indptr[row+1]):value+=data[k]*vector[indices[k]]
        out[row]=value
    return out


@njit(parallel=True,cache=True)
def impact_rates(matrix,energy_counts):
    out=np.zeros(len(matrix))
    for r in prange(len(matrix)):
        value=0.
        for k in range(matrix.shape[1]):value+=matrix[r,k]*energy_counts[k]
        out[r]=value
    return out


@njit(parallel=True,cache=True)
def impact_events(f,events,sv,gap,multiplicity,ee,angular_weights):
    """Each reaction has an independent loss/gain array at the same frozen f.

    events come from the coupled heavy-state matrix exponential, not from
    independent reaction trajectories. Isotropic outgoing electrons, equal
    ionization sharing; true sigma(E), no Maxwellian rate substitution.
    """
    nr=len(events);ne,na=f.shape
    losses=np.zeros((nr,ne,na));gains=np.zeros_like(losses);errors=np.zeros(nr)
    for r in prange(nr):
        if events[r]==0:continue
        total=0.
        for k in range(ne):
            for a in range(na):total+=sv[r,k]*f[k,a]
        if total==0:continue
        for k in range(ne):
            count=0.
            for a in range(na):
                number=events[r]*sv[r,k]*f[k,a]/total
                losses[r,k,a]=number;count+=number
            e=(ee[k]-gap[r])/multiplicity[r]
            j,q=bracket(ee,e)
            represented=(1-q)*ee[j]+q*ee[j+1]
            errors[r]+=count*multiplicity[r]*(represented-e)
            for a in range(na):
                gains[r,j,a]+=count*multiplicity[r]*(1-q)*angular_weights[a]
                gains[r,j+1,a]+=count*multiplicity[r]*q*angular_weights[a]
    return losses,gains,errors


@njit(parallel=True,cache=True)
def pair_geometry(ez,ep,iz,ip,azimuths,mass_i=MI):
    shape=(len(ez),len(iz),len(azimuths));rel=np.empty(shape)
    cx=np.empty(shape);cy=np.empty(shape);cz=np.empty(shape)
    reduced=ME*mass_i/(ME+mass_i)
    for a in prange(len(ez)):
        for b in range(len(iz)):
            for p in range(len(azimuths)):
                ix=ip[b]*math.cos(azimuths[p]);iy=ip[b]*math.sin(azimuths[p])
                gx=ep[a]-ix;gy=-iy;gz=ez[a]-iz[b]
                rel[a,b,p]=reduced*(gx*gx+gy*gy+gz*gz)/(2*EV)
                cx[a,b,p]=(ME*ep[a]+mass_i*ix)/(ME+mass_i)
                cy[a,b,p]=mass_i*iy/(ME+mass_i);cz[a,b,p]=(ME*ez[a]+mass_i*iz[b])/(ME+mass_i)
    return rel,cx,cy,cz


@njit(parallel=True,cache=True)
def build_rr(rel,threshold,g,table_energy,table_sigma,photons,h,c,mass_i=MI,parent_g=2.):
    nr=len(threshold);shape=(nr,rel.size)
    rate=np.zeros(shape);index=np.full(shape,-1,dtype=np.int32);mix=np.zeros(shape)
    flat=rel.ravel();reduced=ME*mass_i/(ME+mass_i)
    for r in prange(nr):
        for k in range(len(flat)):
            e=flat[k];ph=threshold[r]+e
            if e<=0:continue
            x=table_energy[np.int64(r)];s=table_sigma[np.int64(r)]
            sig=np.interp(ph,x,s)
            if ph>x[-1]:sig=s[-1]*(x[-1]/ph)**3
            speed=math.sqrt(2*e*EV/reduced)
            # Milne, with reduced mass and relative speed.
            rate[r,k]=g[r]/parent_g*(ph*EV/(reduced*c*speed))**2*sig*speed
            if photons[0]<=ph<=photons[-1]:
                j,a=bracket(photons,ph);index[r,k]=j;mix[r,k]=a
    return rate,index,mix


@njit(parallel=True,cache=True)
def rr_emission(fe,fi,rate,index,mix,rel,ee,ie,photons,active):
    nr=rate.shape[0];ne=len(fe);ni=len(fi);naz=rel.shape[2]
    spectrum=np.zeros((nr,len(photons)));eloss=np.zeros((nr,ne));iloss=np.zeros((nr,ni))
    photon_energy=np.zeros(nr);cm_energy=np.zeros(nr);missing=np.zeros(nr)
    density=fe.sum()*fi.sum();threshold=np.zeros(nr)
    if density==0:return spectrum,eloss,iloss,photon_energy,cm_energy,missing
    for r in prange(nr):
        if not active[r]:continue
        whole=0.;outside=0.
        for a in range(ne):
            if fe[a]==0:continue
            for b in range(ni):
                if fi[b]==0:continue
                for p in range(naz):
                    k=(a*ni+b)*naz+p
                    w=fe[a]*fi[b]/density/naz*rate[r,k];whole+=w
                    j=index[r,k]
                    if j<0:outside+=w;continue
                    t=mix[r,k]
                    spectrum[r,j]+=w*(1-t);spectrum[r,j+1]+=w*t
                    eloss[r,a]+=w;iloss[r,b]+=w
                    ph=photons[j]*(1-t)+photons[j+1]*t
                    photon_energy[r]+=w*ph
                    cm_energy[r]+=w*(ee[a]+ie[b]-rel[a,b,p])
        missing[r]=outside/max(whole,1e-300)
    return spectrum,eloss,iloss,photon_energy,cm_energy,missing


@njit(parallel=True,cache=True)
def rr_stimulated(fe,fi,rate,index,mix,rel,ee,ie,photons,occupation,active):
    nr=rate.shape[0];ne=len(fe);ni=len(fi);naz=rel.shape[2]
    eloss=np.zeros((nr,ne));iloss=np.zeros((nr,ni));photon_energy=np.zeros(nr);cm_energy=np.zeros(nr)
    density=fe.sum()*fi.sum()
    if density==0:return eloss,iloss,photon_energy,cm_energy
    for r in prange(nr):
        if not active[r]:continue
        for a in range(ne):
            if fe[a]==0:continue
            for b in range(ni):
                if fi[b]==0:continue
                for p in range(naz):
                    k=(a*ni+b)*naz+p;j=index[r,k]
                    if j<0:continue
                    t=mix[r,k];occ=(1-t)*occupation[j]+t*occupation[j+1]
                    w=fe[a]*fi[b]/density/naz*rate[r,k]*occ
                    eloss[r,a]+=w;iloss[r,b]+=w
                    photon_energy[r]+=w*(photons[j]*(1-t)+photons[j+1]*t)
                    cm_energy[r]+=w*(ee[a]+ie[b]-rel[a,b,p])
    return eloss,iloss,photon_energy,cm_energy


@njit(parallel=True,cache=True)
def rr_heavy_products(fe,fi,rate,index,mix,rel,cx,cy,cz,heavy_e,heavy_mu,mass_i,occupation,stimulated):
    """Recombined heavy-state phase source from the SAME pair event weights.

    Energy is the pair centre-of-mass energy. Nuclear-mass approximation and
    photon recoil remain declared; projection defects are measured outside.
    """
    nr=rate.shape[0];out=np.zeros((nr,len(heavy_e)*len(heavy_mu)));cm=np.zeros(nr)
    density=fe.sum()*fi.sum();naz=rel.shape[2]
    if density==0:return out,cm
    for r in prange(nr):
        for a in range(len(fe)):
            if fe[a]==0:continue
            for b in range(len(fi)):
                if fi[b]==0:continue
                for p in range(naz):
                    k=(a*len(fi)+b)*naz+p;j=index[r,k]
                    if j<0:continue
                    w=fe[a]*fi[b]/density/naz*rate[r,k]
                    if stimulated:w*=occupation[j]*(1-mix[r,k])+occupation[j+1]*mix[r,k]
                    energy=(ME+mass_i)*(cx[a,b,p]**2+cy[a,b,p]**2+cz[a,b,p]**2)/(2*EV)
                    vz=(ME+mass_i)/mass_i*cz[a,b,p]
                    put_phase(out[r],w,energy,vz,heavy_e,heavy_mu,mass_i);cm[r]+=w*energy
    return out,cm


@njit(cache=True)
def momentum_transfer_rate(relative_energy,screen_cm,log_cap,mass_i=MI,ion_charge=1.):
    """g*sigma_mt, classical Rutherford with b<=screen_cm, capped log.

    sigma_mt=2*pi*b0^2*log(1+(screen/b0)^2), b0=e^2/(mu*g^2).
    Effective ISOTROPIC binary scattering uses this transport XS. Higher
    angular moments are a model approximation, not a Landau operator.
    """
    if relative_energy<=0:return 0.
    reduced=ME*mass_i/(ME+mass_i);speed=math.sqrt(2*relative_energy*EV/reduced)
    b0=abs(ion_charge)*2.30707755234e-19/(2*relative_energy*EV)
    logterm=min(.5*math.log1p((screen_cm/b0)**2),log_cap)
    return 4*math.pi*b0*b0*logterm*speed


@njit(parallel=True,cache=True)
def binary_ei_generator(fi,rel,cx,cy,cz,ee,emu,ie,imu,screen,logcap,scatter_mu,scatter_weight,scatter_phi,mass_i=MI,ion_charge=1.):
    """Frozen-ion electron Markov generator and ion event responses.

    Output rows correspond to the initial electron phase node. Repeated
    electron collisions are integrated by positive uniformization outside
    this kernel; ion updates use the SAME electron residence integrals.
    """
    ne=rel.shape[0];ni=len(fi);naz=rel.shape[2]
    egain=np.zeros((ne,ne));igain=np.zeros((ne,ni));iloss=np.zeros((ne,ni));frequency=np.zeros(ne)
    reduced=ME*mass_i/(ME+mass_i);totalmass=ME+mass_i
    for a in prange(ne):
        for b in range(ni):
            if fi[b]==0:continue
            for p in range(naz):
                pair=fi[b]*momentum_transfer_rate(rel[a,b,p],screen,logcap,mass_i,ion_charge)/naz
                if pair==0:continue
                frequency[a]+=pair;iloss[a,b]+=pair
                speed=math.sqrt(2*rel[a,b,p]*EV/reduced)
                for t in range(len(scatter_mu)):
                    z=speed*scatter_mu[t];perp=speed*math.sqrt(1-scatter_mu[t]**2)
                    for phi in scatter_phi:
                        x=perp*math.cos(phi);y=perp*math.sin(phi)
                        ex=cx[a,b,p]+mass_i/totalmass*x;ey=cy[a,b,p]+mass_i/totalmass*y;ez=cz[a,b,p]+mass_i/totalmass*z
                        ix=cx[a,b,p]-ME/totalmass*x;iy=cy[a,b,p]-ME/totalmass*y;iz=cz[a,b,p]-ME/totalmass*z
                        en=ME*(ex*ex+ey*ey+ez*ez)/(2*EV);inrg=mass_i*(ix*ix+iy*iy+iz*iz)/(2*EV)
                        w=pair*scatter_weight[t]/len(scatter_phi)
                        put_phase(egain[a],w,en,ez,ee,emu,ME)
                        put_phase(igain[a],w,inrg,iz,ie,imu,mass_i)
    return egain,igain,iloss,frequency


@njit(parallel=True,cache=True)
def binary_generator_rates(fi,rel,cx,cy,cz,ee,emu,ie,imu,rates,scatter_mu,scatter_weight,scatter_phi,mass_i):
    """Frozen-ion electron Markov generator and ion event responses.

    Output rows correspond to the initial electron phase node. Repeated
    electron collisions are integrated by positive uniformization outside
    this kernel; ion updates use the SAME electron residence integrals.
    """
    ne=rel.shape[0];ni=len(fi);naz=rel.shape[2]
    egain=np.zeros((ne,ne));igain=np.zeros((ne,ni));iloss=np.zeros((ne,ni));frequency=np.zeros(ne)
    reduced=ME*mass_i/(ME+mass_i);totalmass=ME+mass_i
    for a in prange(ne):
        for b in range(ni):
            if fi[b]==0:continue
            for p in range(naz):
                pair=fi[b]*rates[a,b,p]/naz
                if pair==0:continue
                frequency[a]+=pair;iloss[a,b]+=pair
                speed=math.sqrt(2*rel[a,b,p]*EV/reduced)
                for t in range(len(scatter_mu)):
                    z=speed*scatter_mu[t];perp=speed*math.sqrt(1-scatter_mu[t]**2)
                    for phi in scatter_phi:
                        x=perp*math.cos(phi);y=perp*math.sin(phi)
                        ex=cx[a,b,p]+mass_i/totalmass*x;ey=cy[a,b,p]+mass_i/totalmass*y;ez=cz[a,b,p]+mass_i/totalmass*z
                        ix=cx[a,b,p]-ME/totalmass*x;iy=cy[a,b,p]-ME/totalmass*y;iz=cz[a,b,p]-ME/totalmass*z
                        en=ME*(ex*ex+ey*ey+ez*ez)/(2*EV);inrg=mass_i*(ix*ix+iy*iy+iz*iz)/(2*EV)
                        w=pair*scatter_weight[t]/len(scatter_phi)
                        put_phase(egain[a],w,en,ez,ee,emu,ME)
                        put_phase(igain[a],w,inrg,iz,ie,imu,mass_i)
    return egain,igain,iloss,frequency
