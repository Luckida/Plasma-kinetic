"""Plane-parallel scalar S_N transfer, with an independently directed pump.

z is Cartesian; coefficients are uniform within each cell. Diffuse emission
and elastic isotropic scattering use Gauss ordinates in both hemispheres.
The collimated incident beam is integrated at its specified mu, independently
of the diffuse angular grid. Spectral powers below are per unit face area.
"""
from functools import lru_cache
import numpy as np
from photon_transport import kernels


def validate_inputs(edges_cm,opacity,emission,beam_4piJ,scatter=None,
                   angular_points=8,beam_mu=1.,beam_side='left',
                   left_4piI=None,right_4piI=None):
    edges=np.asarray(edges_cm,dtype=float)
    if edges.ndim!=1 or not len(edges) or not np.all(np.isfinite(edges)) or edges[0]<=0 or np.any(np.diff(edges)<=0):
        raise ValueError('Slab right edges must increase from z=0')
    widths=np.diff(np.r_[0.,edges]);ka=np.asarray(opacity,dtype=float);eps=np.asarray(emission,dtype=float)
    if ka.ndim!=2 or ka.shape!=eps.shape or len(ka)!=len(edges):raise ValueError('Slab spectra must have shape (cell, frequency)')
    cells,freq=ka.shape;beam=np.asarray(beam_4piJ,dtype=float)
    left=np.zeros(freq) if left_4piI is None else np.asarray(left_4piI,dtype=float)
    right=np.zeros(freq) if right_4piI is None else np.asarray(right_4piI,dtype=float)
    ks=np.zeros(cells) if scatter is None else np.asarray(scatter,dtype=float)
    if ks.shape!=(cells,) or np.any(ks<0) or not np.all(np.isfinite(ks)):raise ValueError('Invalid slab scattering')
    if beam.shape!=(freq,) or left.shape!=(freq,) or right.shape!=(freq,):raise ValueError('Invalid slab boundary shape')
    if any(not np.all(np.isfinite(x)) for x in (ka,eps,beam,left,right)) or any(np.any(x<0) for x in (eps,beam,left,right)):
        raise ValueError('Slab boundaries/emission must be finite and nonnegative')
    if isinstance(angular_points,bool) or not isinstance(angular_points,int) or angular_points<2:raise ValueError('At least 2 diffuse angular points per hemisphere')
    if isinstance(beam_mu,bool) or not np.isfinite(beam_mu) or not 0<beam_mu<=1 or beam_side not in ('left','right'):raise ValueError('beam_mu must be in (0,1], side left/right')
    return edges,widths,ka,eps,beam,left,right,ks,cells,freq


@lru_cache(maxsize=16)
def quadrature(points):
    nodes,weights=np.polynomial.legendre.leggauss(points)
    return (nodes+1)/2,weights/2


def solve_transfer(edges_cm,opacity,emission,beam_4piJ,scatter=None,
                   angular_points=8,beam_mu=1.,beam_side='left',
                   left_4piI=None,right_4piI=None,backend='numpy',workers=1,internal_beams=None):
    args=validate_inputs(edges_cm,opacity,emission,beam_4piJ,scatter,angular_points,beam_mu,beam_side,left_4piI,right_4piI)
    if backend not in ('numpy','numba','auto'):raise ValueError('radiation_backend: numpy, numba or auto')
    if isinstance(workers,bool) or not isinstance(workers,int) or workers<1:raise ValueError('workers must be a positive integer')
    if backend!='numpy':
        try:
            from slab_transport_fast import solve_fast
        except ImportError:
            if backend=='numba':raise ValueError('Install requirements_fast.txt for the numba backend')
        else:
            result=solve_fast(args,angular_points,beam_mu,beam_side,workers)
            return add_internal_beams(result,args,internal_beams,angular_points,backend,workers) if internal_beams else result
    edges,widths,ka,eps,beam,left,right,ks,cells,freq=args
    mus,weights=quadrature(angular_points)
    total=ka+ks[:,None]
    lam=np.zeros((cells,cells,freq));direct=np.zeros((cells,freq))
    escape=np.zeros((2,cells,freq));boundary_out=np.zeros((2,freq))
    # q=4*pi*I; mean 4*pi*J = 1/2 int[-1,1] q dmu.
    for side,incoming in enumerate((left,right)):
        order=range(cells) if side==0 else range(cells-1,-1,-1)
        for mu,weight in zip(mus,weights):
            if np.any(np.sum(total*widths[:,None],axis=0)/mu < -20):raise ValueError('Slab ray gain tau<-20; finite lateral geometry required')
            source=np.zeros((cells,freq));boundary=incoming.copy()
            for cell in order:
                length=widths[cell]/mu;exponential,f,g=kernels(total[cell]*length)
                lam[cell]+=.5*weight*source*f[None,:]
                lam[cell,cell]+=.5*weight*length*g
                direct[cell]+=.5*weight*boundary*f
                source*=exponential[None,:];source[cell]+=length*f
                boundary*=exponential
            escape[1-side]+=.5*weight*mu*source
            boundary_out[1-side]+=.5*weight*mu*boundary
    incoming_diffuse=.5*sum(weight*mu for mu,weight in zip(mus,weights))*(left+right)
    # A beam's scalar integral is Psi; normal energy flux is |mu| Psi.
    if np.any(np.sum(total*widths[:,None],axis=0)/beam_mu < -20):raise ValueError('Slab beam gain tau<-20; finite lateral geometry required')
    beam_direct=beam.copy();order=range(cells) if beam_side=='left' else range(cells-1,-1,-1)
    for cell in order:
        exponential,f,_=kernels(total[cell]*widths[cell]/beam_mu)
        direct[cell]+=beam_direct*f;beam_direct*=exponential
    boundary_out[1 if beam_side=='left' else 0]+=beam_mu*beam_direct
    rhs=direct+np.einsum('ijf,jf->if',lam,eps)
    if np.any(ks):
        matrix=np.eye(cells)[None,:,:]-np.moveaxis(lam,2,0)*ks[None,None,:]
        try:mean=np.linalg.solve(matrix,rhs.T[...,None])[...,0].T
        except np.linalg.LinAlgError as error:raise ValueError('Singular slab scattering operator') from error
    else:mean=rhs
    if not np.all(np.isfinite(mean)) or np.any(mean<0):raise ValueError('Nonphysical slab radiation field')
    effective=eps+ks[:,None]*mean
    outgoing=boundary_out+np.einsum('sif,if->sf',escape,effective)
    incoming=beam_mu*beam+incoming_diffuse
    emitted=widths@eps;absorbed=widths@(ka*mean)
    residual=incoming+emitted-outgoing.sum(axis=0)-absorbed
    gross=abs(incoming)+abs(emitted)+abs(outgoing).sum(axis=0)+widths@abs(ka*mean)
    optical_depth=widths@ka
    result=dict(mean_4piJ=mean,outgoing_left_flux=outgoing[0],outgoing_right_flux=outgoing[1],
                incoming_flux=incoming,true_emitted_flux=emitted,true_absorbed_flux=absorbed,
                transmitted_direct_beam_flux=beam_mu*beam_direct,
                spectral_balance_residual=residual,maximum_relative_balance_residual=float(np.max(abs(residual)/np.maximum(gross,1e-300))),
                widths_cm=widths,lambda_operator=lam,normal_optical_depth=optical_depth,
                minimum_normal_optical_depth=float(optical_depth.min()),negative_opacity_present=bool(np.any(ka<0)),
                angular_model='finite discrete ordinates; check grazing-angle convergence',minimum_diffuse_mu=float(mus.min()),backend='numpy',workers=1)
    return add_internal_beams(result,args,internal_beams,angular_points,backend,workers) if internal_beams else result


def add_internal_beams(result,args,beams,angular_points,backend,workers):
    """Exact sheet injection, including a source inside a finite-volume cell.

    Linearity gives (I-Lambda*ks) delta_J=direct_sheet_J. Scattered sheet
    photons then escape through the same formal solver as plasma photons.
    This avoids treating the sheet as uniform isotropic cell emissivity.
    """
    edges,widths,ka,eps,beam,left,right,ks,cells,freq=args
    lower=np.r_[0.,edges[:-1]];total=ka+ks[:,None]
    direct=np.zeros_like(ka);escaped=np.zeros((2,freq));injected=np.zeros(freq)
    for sheet in beams:
        pos=float(sheet['position_cm']);mu=float(sheet.get('mu',1.));direction=sheet['direction'];spectrum=np.asarray(sheet['spectrum'],dtype=float)
        if not np.isfinite(pos+mu) or not 0<=pos<=edges[-1] or not 0<mu<=1 or direction not in ('left','right'):raise ValueError('Invalid internal photon sheet')
        if spectrum.shape!=(freq,) or np.any(spectrum<0) or not np.all(np.isfinite(spectrum)):raise ValueError('Invalid internal source spectrum')
        injected+=mu*spectrum;current=spectrum.copy();optical=np.zeros(freq)
        order=range(cells) if direction=='right' else range(cells-1,-1,-1)
        for z in order:
            traversed=max(0.,edges[z]-max(lower[z],pos)) if direction=='right' else max(0.,min(edges[z],pos)-lower[z])
            if traversed==0:continue
            tau=total[z]*traversed/mu;optical+=tau
            if np.any(optical < -20):raise ValueError('Internal slab ray gain tau<-20')
            exponential,f,_=kernels(tau)
            direct[z]+=current*(traversed/widths[z])*f;current*=exponential
        escaped[1 if direction=='right' else 0]+=mu*current
    lam=result['lambda_operator']
    if np.any(ks):
        matrix=np.eye(cells)[None,:,:]-np.moveaxis(lam,2,0)*ks[None,None,:]
        delta=np.linalg.solve(matrix,direct.T[...,None])[...,0].T
        scattered=solve_transfer(edges,ka+ks[:,None],ks[:,None]*delta,np.zeros(freq),
                                 angular_points=angular_points,backend=backend,workers=workers)
        escaped[0]+=scattered['outgoing_left_flux'];escaped[1]+=scattered['outgoing_right_flux']
    else:delta=direct
    if np.any(delta<0) or not np.all(np.isfinite(delta)):raise ValueError('Nonphysical internal radiation field')
    result['mean_4piJ']=result['mean_4piJ']+delta
    result['outgoing_left_flux']=result['outgoing_left_flux']+escaped[0]
    result['outgoing_right_flux']=result['outgoing_right_flux']+escaped[1]
    result['internal_injected_flux']=injected
    result['true_absorbed_flux']=widths@(ka*result['mean_4piJ'])
    outgoing=result['outgoing_left_flux']+result['outgoing_right_flux']
    residual=result['incoming_flux']+injected+result['true_emitted_flux']-outgoing-result['true_absorbed_flux']
    gross=abs(result['incoming_flux'])+injected+abs(result['true_emitted_flux'])+abs(outgoing)+widths@abs(ka*result['mean_4piJ'])
    result['spectral_balance_residual']=residual
    result['maximum_relative_balance_residual']=float(np.max(abs(residual)/np.maximum(gross,1e-300)))
    return result
