import numpy as np

def stf_gaussian(t,f):
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