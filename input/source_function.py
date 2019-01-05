"""This script contains functions to export source time function given 
a certain time vector and parameters necessary."""

import numpy as np
def gaussian(t,a,b,c):
    """.. function ::
    
    This function computes a gaussian pulse given the parameters a,b,c 
    and a time vector t

    :param t: 1D ``numpy`` array
    :param a: ``float``
    :param b: ``float``
    :param c: ``float``
    :rtype: 1D ``numpy`` array of the same size as t


    Mathematical formulation:

    .. math::
        
        f(x)=ae^{-{\frac {(x-b)^{2}}{2c^{2}}}}

    """
    
    # General definition of the Gaussian Pulse
    f = a*np.exp( - ( (t-b)**2) / (2*c**2) )

    return f
