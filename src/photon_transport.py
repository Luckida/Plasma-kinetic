"""Stable formal-characteristic attenuation kernels used by the slab transfer solver."""
import numpy as np


def kernels(tau):
    """exp(-t), (1-exp(-t))/t, [1-(1-exp(-t))/t]/t."""
    t=np.asarray(tau,dtype=float)
    if np.any(t < -20):raise ValueError('Gain tau < -20 requires a different transfer treatment')
    f=np.empty_like(t);g=np.empty_like(t)
    small=abs(t)<1e-3;x=t[small]
    # Тот же ряд, Horner без дорогих степеней малых отрицательных чисел.
    f[small]=1+x*(-.5+x*(1/6+x*(-1/24+x*(1/120-x/720))))
    g[small]=.5+x*(-1/6+x*(1/24+x*(-1/120+x*(1/720-x/5040))))
    f[~small]=-np.expm1(-t[~small])/t[~small]
    g[~small]=(1-f[~small])/t[~small]
    return np.exp(-t),f,g
