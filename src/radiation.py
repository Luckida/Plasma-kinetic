"""Spectral trapezoidal quadrature; populations and transfer are in the mixture and slab modules."""
import numpy as np


def weights(x):
    w = np.empty_like(x)
    w[0], w[-1] = (x[1] - x[0]) / 2, (x[-1] - x[-2]) / 2
    w[1:-1] = (x[2:] - x[:-2]) / 2
    return w
