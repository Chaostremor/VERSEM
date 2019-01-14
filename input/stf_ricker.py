import numpy as np

def stf_ricker(t,f):
    """Computes the Ricker Wavelet also known as the Mexican hat function

    .. math::

        s(t) (1 - 2 * a_ricker * (t-t0))**2) * np.exp(-(a_ricker * ((t-t[int(t0)]))) ** 2) 

    """

    # delta 
    dt = t[1] - t[0]

    # Origin shift index
    t0 = 1/f/dt

    # Ricker factor
    a_ricker = 4 * f

    s = (1 - 2 * a_ricker * (t-t[int(t0)])**2) * np.exp(-(a_ricker * ((t-t[int(t0)]))) ** 2) 

    return s