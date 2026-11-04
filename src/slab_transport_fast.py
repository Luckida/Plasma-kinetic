"""Compiled, frequency-parallel implementation of the reference slab operator.

No fastmath, changed tolerances, reduced angles or dropped scattering. Each
frequency has private working arrays. A small pivoted solve avoids BLAS and
nested threading inside the parallel loop.
"""
import math
import numpy as np
from numba import njit,prange,get_num_threads,set_num_threads,config
from slab_transport import quadrature


@njit(cache=True)
def scalar_kernels(t):
    if abs(t)<1e-3:
        f=1+t*(-.5+t*(1/6+t*(-1/24+t*(1/120-t/720))))
        g=.5+t*(-1/6+t*(1/24+t*(-1/120+t*(1/720-t/5040))))
    else:
        f=-math.expm1(-t)/t;g=(1-f)/t
    return math.exp(-t),f,g


@njit(cache=True)
def pivot_solve(a,b):
    n=len(b)
    for k in range(n):
        pivot=k
        for i in range(k+1,n):
            if abs(a[i,k])>abs(a[pivot,k]):pivot=i
        if abs(a[pivot,k])<1e-30:return np.full(n,np.nan)
        if pivot!=k:
            for j in range(n):a[k,j],a[pivot,j]=a[pivot,j],a[k,j]
            b[k],b[pivot]=b[pivot],b[k]
        for i in range(k+1,n):
            multiplier=a[i,k]/a[k,k]
            for j in range(k+1,n):a[i,j]-=multiplier*a[k,j]
            b[i]-=multiplier*b[k]
    x=np.empty(n)
    for i in range(n-1,-1,-1):
        value=b[i]
        for j in range(i+1,n):value-=a[i,j]*x[j]
        x[i]=value/a[i,i]
    return x


@njit(parallel=True,cache=True)
def transfer_kernel(widths,ka,eps,beam,left,right,ks,mus,weights,beam_mu,side):
    cells,freq=ka.shape
    means=np.empty((cells,freq));outs=np.empty((2,freq));beam_out=np.empty(freq);all_lam=np.empty((cells,cells,freq))
    for fidx in prange(freq):
        lam=np.zeros((cells,cells));direct=np.zeros(cells);escape=np.zeros((2,cells));bout=np.zeros(2)
        for incoming_side in range(2):
            for angle in range(len(mus)):
                mu=mus[angle];weight=weights[angle];source=np.zeros(cells)
                boundary=left[fidx] if incoming_side==0 else right[fidx]
                for order in range(cells):
                    cell=order if incoming_side==0 else cells-1-order
                    length=widths[cell]/mu
                    exponential,ff,gg=scalar_kernels((ka[cell,fidx]+ks[cell])*length)
                    for j in range(cells):lam[cell,j]+=.5*weight*source[j]*ff
                    lam[cell,cell]+=.5*weight*length*gg
                    direct[cell]+=.5*weight*boundary*ff
                    source*=exponential;source[cell]+=length*ff;boundary*=exponential
                for j in range(cells):escape[1-incoming_side,j]+=.5*weight*mu*source[j]
                bout[1-incoming_side]+=.5*weight*mu*boundary
        bdirect=beam[fidx]
        for order in range(cells):
            cell=order if side==0 else cells-1-order
            exponential,ff,_=scalar_kernels((ka[cell,fidx]+ks[cell])*widths[cell]/beam_mu)
            direct[cell]+=bdirect*ff;bdirect*=exponential
        bout[1-side]+=beam_mu*bdirect;beam_out[fidx]=beam_mu*bdirect
        matrix=np.eye(cells);rhs=direct.copy()
        for i in range(cells):
            for j in range(cells):
                rhs[i]+=lam[i,j]*eps[j,fidx];matrix[i,j]-=lam[i,j]*ks[j]
        mean=pivot_solve(matrix,rhs)
        for i in range(cells):
            means[i,fidx]=mean[i]
            for j in range(cells):all_lam[i,j,fidx]=lam[i,j]
        for output_side in range(2):
            value=bout[output_side]
            for j in range(cells):value+=escape[output_side,j]*(eps[j,fidx]+ks[j]*mean[j])
            outs[output_side,fidx]=value
    return means,outs,beam_out,all_lam


def solve_fast(args,angular_points,beam_mu,beam_side,workers):
    edges,widths,ka,eps,beam,left,right,ks,cells,freq=args
    mus,weights=quadrature(angular_points);total=ka+ks[:,None]
    tau=widths@total
    if np.any(tau/min(float(mus.min()),beam_mu)<-20):raise ValueError('Slab ray gain tau<-20; finite lateral geometry required')
    active=min(workers,config.NUMBA_NUM_THREADS)
    old=get_num_threads()
    try:
        set_num_threads(active)
        mean,outgoing,beam_direct,lam=transfer_kernel(widths,ka,eps,beam,left,right,ks,mus,weights,beam_mu,0 if beam_side=='left' else 1)
    finally:set_num_threads(old)
    if not np.all(np.isfinite(mean)) or np.any(mean<0):raise ValueError('Nonphysical slab radiation field')
    incoming=beam_mu*beam+.5*float(weights@mus)*(left+right)
    emitted=widths@eps;absorbed=widths@(ka*mean);residual=incoming+emitted-outgoing.sum(axis=0)-absorbed
    gross=abs(incoming)+abs(emitted)+abs(outgoing).sum(axis=0)+widths@abs(ka*mean)
    optical_depth=widths@ka
    return dict(mean_4piJ=mean,outgoing_left_flux=outgoing[0],outgoing_right_flux=outgoing[1],incoming_flux=incoming,
        true_emitted_flux=emitted,true_absorbed_flux=absorbed,transmitted_direct_beam_flux=beam_direct,
        spectral_balance_residual=residual,maximum_relative_balance_residual=float(np.max(abs(residual)/np.maximum(gross,1e-300))),
        widths_cm=widths,lambda_operator=lam,normal_optical_depth=optical_depth,minimum_normal_optical_depth=float(optical_depth.min()),
        negative_opacity_present=bool(np.any(ka<0)),angular_model='finite discrete ordinates; check grazing-angle convergence',minimum_diffuse_mu=float(mus.min()),backend='numba',workers=active)
