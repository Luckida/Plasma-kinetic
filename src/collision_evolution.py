"""Positive integration of a frozen-partner kinetic collision generator.

Uniformization is a sum of nonnegative Poisson-weighted Markov iterates.
No population renormalization or clipping. The omitted tail is explicit.
"""
import math
import numpy as np
from scipy.sparse import csr_matrix
from scipy.special import gammainc
from scipy.sparse import eye
from scipy.sparse.linalg import spsolve


def evolve_frozen(gain,frequency,initial,dt,tail_probability=1e-14):
    """Positive collision step without work proportional to nu_max*dt.

    For short optical depths use accurate uniformization. For a stiff step,
    backward Euler gives a positive M-matrix solve; residence=dt*f_new makes
    the SAME discrete event budget valid for the ion response. Both paths
    retain the global splitting's first-order accuracy. No clipping/rescaling.
    """
    maximum=float(frequency.max()) if len(frequency) else 0.
    if maximum*dt<=16:
        f,residence,tail,pieces=uniformize(gain,frequency,initial,dt,tail_probability)
        return f,residence,tail,pieces,'uniformization'
    generator=csr_matrix(gain.T)-csr_matrix((frequency,(np.arange(len(frequency)),np.arange(len(frequency)))),shape=gain.shape)
    matrix=eye(len(initial),format='csc')-dt*generator.tocsc()
    final=spsolve(matrix,initial)
    if np.any(~np.isfinite(final)) or np.any(final<0):raise ValueError('Implicit collision solve is not positive/finite; reduce the step')
    return final,dt*final,0.,1,'implicit_euler'


def uniformize(gain,frequency,initial,dt,tail_probability=1e-14):
    """Return f(t), integral f dt, omitted Poisson probability and substeps.

    gain[a,b]: transition gain in b per electron in a, s^-1. Self transitions
    remain in event residence times even though they cancel in the generator.
    """
    if dt==0:return initial.copy(),np.zeros_like(initial),0.,0
    maximum=float(frequency.max()) if len(frequency) else 0.
    if maximum==0:return initial.copy(),initial*dt,0.,0
    if np.any(gain<0) or np.any(frequency<0):raise ValueError('Negative collision generator')
    P=gain.T/maximum
    P[np.diag_indices_from(P)]+=1-frequency/maximum
    if np.any(P<0):raise ValueError('Nonpositive uniformization operator')
    matrix=csr_matrix(P)
    pieces=max(1,int(np.ceil(maximum*dt/16)))
    if pieces>10000:raise RuntimeError('Collision uniformization too stiff: refine the grid/step')
    h=dt/pieces;x=maximum*h;current=initial.copy();residence=np.zeros_like(initial);tail_bound=0.
    for _ in range(pieces):
        vector=current.copy();final=np.zeros_like(current);integral=np.zeros_like(current)
        poisson=math.exp(-x);k=0
        while True:
            tail=float(gammainc(k+1,x))
            final+=poisson*vector;integral+=(tail/maximum)*vector
            if tail<tail_probability and k>=x:break
            k+=1;poisson*=x/k;vector=matrix@vector
            if k>10000:raise RuntimeError('Poisson tail did not converge')
        current=final;residence+=integral;tail_bound+=tail
    return current,residence,tail_bound,pieces
