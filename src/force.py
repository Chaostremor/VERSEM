"""implementation of the force-term"""
import numpy as np
from scipy import interpolate


def genforce(force_term,force_location,gll_coordinates):
    """.. function :: force(source_time_function,force_term,
                                        force_location,gll_coordinates):

    This function implements the force term of a point source for the 
    seismic wave equation.

    :param force_term: 2 element ``numpy`` array with forces in 2D
    :param force_loc: 2 element ``numpy`` numpy array with force location 
    :param gll_coordinates:

    """

    # number of elements
    nel = gll_coordinates.shape[0]

    ## Find closest GLL point to forceterm
    # find distance to all gll points
    dist = np.sqrt( (gll_coordinates[:,0]-force_location[0])**2 
                        + (gll_coordinates[:,1]-force_location[1])**2  )
    
    # find GLL index use argsort()[:n] for n first indeces
    index = dist.argsort()[:1]

    # Computing 
    Fx = np.zeros(nel)
    Fy = np.zeros(nel)
    
    Fx[index] = force_term[0]
    Fy[index] = force_term[1]


    return Fx,Fy



def F(tp,source_time,Fx):
    """.. function:: F(t)

    """
    
    # Length of the function 
    N = len(source_time)
    
    # Create repeat matrix
    repF = np.repeat(Fx.reshape([len(Fx),1]),N,axis=1)

    # repeat vector
    fp = repF*source_time
    
    return interpolate.interp1d(tp,fp)





if __name__ == "__main__":
    pass
