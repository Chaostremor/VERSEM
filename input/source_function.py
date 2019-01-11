"""This script contains functions to export source time function given 
a certain time vector and parameters necessary."""

import numpy as np

def gaussian(t,a,b,c):
    """This function computes a gaussian pulse given the parameters a,b,c 
    and a time vector t

    :param t: 1D ``numpy`` array
    :param a: ``float``
    :param b: ``float``
    :param c: ``float``

    :rtype: 1D ``numpy`` array of the same size as t


    Mathematical formulation:

    .. math::
        
        f(x)=ae^{-{\\frac {(x-b)^{2}}{2c^{2}}}}


    """
    
    # General definition of the Gaussian Pulse
    f = a*np.exp( - ( (t-b)**2) / (2*c**2) )

    return f


def ricker(t,t0,p0):
    """This function computes the Ricker Wavelet given a dominant period p0,
    time vector t and origin t0

    :param t: time vector
    :param t0: origin of symmetry vector
    :param p0: dominant period

    :rtype: vector containing ricker function values

    .. math::

        s(t) = -8/p_0(t-t_0)e^{\\frac{(t-t0)^2}{(4/p_0)^2}}

    """

    # Computing the Gaussian
    s = -8/p0*(t-t0)*np.exp( ((t-t0)**2) / ((4/p0)**2))

    return s


def rickerSEM(dt,pt):
    """Taken from seismolive

    """

    nt = int(2 * pt / dt)
    c = np.zeros(nt)
    t0 = pt / dt
    a_ricker = 4 / pt

    t = np.zeros(nt)

    for it in range(0, nt):
        t[it] = ((it + 1) - t0) * dt
        c[it] = -2 * a_ricker * t[it] * np.exp(-(a_ricker * t[it]) ** 2)

    return t+t[int(t0)],c