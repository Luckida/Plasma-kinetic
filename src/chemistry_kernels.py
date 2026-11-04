"""Velocity-resolved two-body chemistry. cgs masses/velocities; energies eV.

Product angles are an explicitly isotropic closure, not differential data.
Jacobi coordinates conserve event energy and vector momentum before the
positive projection onto each species' own energy-angle grid.
"""
import math
import numpy as np
from numba import njit,prange
from kinetic_kernels import put_phase
from particle_kinetics import EV


@njit(cache=True)
def geometry(za,pa,zb,pb,azimuths,ma,mb):
    relative=np.zeros((len(za),len(zb),len(azimuths)))
    for a in range(len(za)):
        for b in range(len(zb)):
            for c in range(len(azimuths)):
                x=pa[a]-pb[b]*math.cos(azimuths[c])
                y=-pb[b]*math.sin(azimuths[c]);z=za[a]-zb[b]
                relative[a,b,c]=.5*ma*mb/(ma+mb)*(x*x+y*y+z*z)/EV
    return relative


@njit(parallel=True,cache=True)
def products(fa,fb,ma,mb,za,pa,zb,pb,azimuths,kernel,dt,factor,q,
             masses,energies,mus,scatter_mu,scatter_w,scatter_phi,electron_release):
    """Return selective losses, products, events, forbidden and overflow counts.

    Two or three outgoing particles. For 3 products, first two form a Jacobi
    pair. The last particle takes min(release, available) eV of Jacobi energy;
    a negative release requests equal splitting (explicit closure).
    """
    nprod=len(masses);size=0
    for i in range(nprod):size=max(size,len(energies[i])*len(mus[i]))
    all_out=np.zeros((len(fa),nprod,size));la=np.zeros(len(fa));all_lb=np.zeros((len(fa),len(fb)))
    events=np.zeros(len(fa));forbidden=np.zeros(len(fa));overflow=np.zeros(len(fa))
    mass_sum=masses.sum()
    for a in prange(len(fa)):
        out=all_out[a];lb=all_lb[a]
        if fa[a]==0:continue
        for b in range(len(fb)):
            if fb[b]==0:continue
            for c in range(len(azimuths)):
                amount=fa[a]*fb[b]*kernel[a,b,c]*dt*factor/len(azimuths)
                if amount==0:continue
                va=np.empty(3);va[0]=pa[a];va[1]=0.;va[2]=za[a]
                vb=np.empty(3);vb[0]=pb[b]*math.cos(azimuths[c]);vb[1]=pb[b]*math.sin(azimuths[c]);vb[2]=zb[b]
                momentum=ma*va+mb*vb;centre=momentum/mass_sum
                kin=.5*(ma*np.dot(va,va)+mb*np.dot(vb,vb))/EV
                available=kin+q-.5*np.dot(momentum,momentum)/mass_sum/EV
                if available < -1e-10:
                    forbidden[a]+=amount;continue
                available=max(available,0.)
                la[a]+=amount;lb[b]+=amount;events[a]+=amount
                if nprod==2:
                    reduced=masses[0]*masses[1]/mass_sum
                    v=math.sqrt(2*available*EV/reduced)
                    for u in range(len(scatter_mu)):
                        transverse=math.sqrt(max(1-scatter_mu[u]**2,0.))
                        for phi in scatter_phi:
                            direction=np.empty(3);direction[0]=transverse*math.cos(phi);direction[1]=transverse*math.sin(phi);direction[2]=scatter_mu[u]
                            velocities=np.empty((2,3))
                            velocities[0]=centre+masses[1]/mass_sum*v*direction
                            velocities[1]=centre-masses[0]/mass_sum*v*direction
                            w=amount*scatter_w[u]/len(scatter_phi)
                            for j in range(2):
                                energy=.5*masses[j]*np.dot(velocities[j],velocities[j])/EV
                                if energy>energies[j][-1]*(1+1e-12):overflow[a]+=w
                                put_phase(out[j],w,energy,velocities[j,2],energies[j],mus[j],masses[j])
                else:
                    pair_mass=masses[0]+masses[1]
                    outer=min(electron_release,available) if electron_release>=0 else .5*available
                    inner=available-outer
                    vouter=math.sqrt(2*outer*EV/(pair_mass*masses[2]/mass_sum))
                    vinner=math.sqrt(2*inner*EV/(masses[0]*masses[1]/pair_mass))
                    for u in range(len(scatter_mu)):
                        transverse=math.sqrt(max(1-scatter_mu[u]**2,0.))
                        for phi in scatter_phi:
                            direction=np.empty(3);direction[0]=transverse*math.cos(phi);direction[1]=transverse*math.sin(phi);direction[2]=scatter_mu[u]
                            pair_centre=centre-masses[2]/mass_sum*vouter*direction
                            third=centre+pair_mass/mass_sum*vouter*direction
                            for v in range(len(scatter_mu)):
                                tr=math.sqrt(max(1-scatter_mu[v]**2,0.))
                                for beta in scatter_phi:
                                    internal=np.empty(3);internal[0]=tr*math.cos(beta);internal[1]=tr*math.sin(beta);internal[2]=scatter_mu[v]
                                    velocities=np.empty((3,3))
                                    velocities[0]=pair_centre+masses[1]/pair_mass*vinner*internal
                                    velocities[1]=pair_centre-masses[0]/pair_mass*vinner*internal
                                    velocities[2]=third
                                    w=amount*scatter_w[u]*scatter_w[v]/len(scatter_phi)**2
                                    for j in range(3):
                                        energy=.5*masses[j]*np.dot(velocities[j],velocities[j])/EV
                                        if energy>energies[j][-1]*(1+1e-12):overflow[a]+=w
                                        put_phase(out[j],w,energy,velocities[j,2],energies[j],mus[j],masses[j])
    return la,np.sum(all_lb,axis=0),np.sum(all_out,axis=0),events.sum(),forbidden.sum(),overflow.sum()
