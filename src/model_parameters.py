"""Function set to define Model parameters

Author: Lucas Sawade

"""


def velocity_conversion(rho,vp,vs):
    """.. function:: velocity_conversion(rho,vp,vs)

    Converts the density, vp and vs velocity to the Lame parameters

    :param rho: 1D ``numpy`` array containing density values
    :param vp: 1D ``numpy`` array containing density values
    :param vs: 1D ``numpy`` array containing density values
    :rtype: 2 element tupel() containing mu and lambda ``numpy`` arrays 
            of same size as the parameters.
                                                                    


    Mathematical solution
    .. math::

        \\begin{align}
            \\mu &= \\rho \\cdot V s ^ { 2 }\\\\
            \\lambda &= ( V p ) ^ { 2 } \\cdot \\rho - 2 \\mu
        \\end{align}

    """
    
    # Computing mu
    mu  = rho*vs**2

    # Computing lambda
    lmd = vp**2 *rho - 2*mu

    return (mu,lmd)


    
