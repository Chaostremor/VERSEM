"""This script contains functions to export source time function given 
a certain time vector and parameters necessary."""

import numpy as np


def delta_function(t):
    """ This function takes in the time stepping vector and outputs
    a zero vector and the time step that is 1/100th of all time steps.

    :param t: 1D ``numpy`` array

    :rtype: delta source_time function vector

    """

    # Length of t
    nt = len(t)

    # initialize new array
    s = np.zeros(nt)

    # Find index where to :math:`t_i = t(n_t/100)
    s[int(nt/100)] = 1 

    return s


def gaussian(t,f):
    """This function computes a gaussian pulse given the parameters a,b,c 
    and a time vector t

    :param t: 1D ``numpy`` array
    :param f: ``float`` peak frequency

    :rtype: 1D ``numpy`` array of the same size as t


    Mathematical formulation:

    .. math::
        
        f(x)=ae^{-{\\frac {(x-b)^{2}}{2c^{2}}}}


    """

    # delta 
    dt = t[1] - t[0]
    
    # time shift ind
    t0 = 4.5/f/dt

    # Amplitude
    a = 1

    # Half-duration?
    c = 1/f

    # General definition of the Gaussian Pulse
    f = a*np.exp( - ( (t-t[int(t0)])**2) / (2*c**2) )

    return f


def ricker(t,f):
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


def rickerINT(t,f):
    """Modified function from seismo-live.org. Computes the
    Integral of the Ricker Wavelet
    """

    # delta 
    dt = t[1] - t[0]

    # Origin shift index
    t0 = 1/f / dt

    # Ricker factor
    a_ricker = 4 * f

    s = -2 * a_ricker * (t-t[int(t0)]) * np.exp(-(a_ricker * ((t-t[int(t0)]))) ** 2) 

    return s