"""implementation of the force-term"""

import numpy as np
from scipy import interpolate
import src.gll_library as gll
import src.loc2glob as l2g

#######################################################################
###            Constructing Element Force Matrices                  ###
#######################################################################

def element_force_mat(f,J,W):
    """Computes the Force Vector for each element

    :param f:  [ngll_el] ``numpy`` array

    :param J: flattened Jacobian determinant vector [ngll_el]
              ``numpy``  corresponding to nodes numbers

    :param W: weight vector [ngll_el] ``numpy``  corresponding to 
              nodes numbers
    
    :rtype: ``numpy`` [ngll_total] 

    The description of the Jacobian can be found on the theory
    documentation.

    """

    #Retrieves the number of gll points in an element
    ngll_el = len(J)

    #the element mass matrix
    fe = np.zeros([ngll_el,ngll_el])
    for i in range(ngll_el):
        fe[i,i] = f[i]*J[i]*W[i]	#basically a diagonal matrix because of the nature of Shape Functions

    return fe

def global_force_vector(gll_coordinates,gll_connect,f,dN_local,W):
    """ Computes the Global Force Vector

    :param gll_coordinates: ``numpy`` array of size [ngll_total]x[2] 
                            containing the coordinates of all the gll points

    :param gll_connect: ``numpy`` array of size [el_no]x[ngll_el]. Contains 
                        the global indexing of gll nodes.

    :param f: force  vector [ngll_el] ``numpy``

    :param dN_local: Local derivative of shape functions at each gll point 
                     in an element. `numpy`` array of size 
                     [total ngll]x[2]x[total ngll]

    :param W: Weight matrix [ngll_el] ``numpy`` corresponding to nodes numbers
    
    :rtype: ``numpy`` [ngll_total]x[ngll_total] array

    """
    
    #Retrieving the number of elements in the domain and the number of gll points per element
    el_no = len(gll_connect)
    ngll_el = len(gll_connect[0])	#Assuming number of gll points per element in constant

    Fg = np.zeros([len(gll_coordinates),len(gll_coordinates)])
    
    for i in range(el_no):    #el_no is the total number of elements
        gll_coords_el = gll_coordinates[gll_connect[i]]
        J_el = np.zeros(len(gll_coords_el))
        for j in range(len(gll_coords_el)):
            #Jacobian for the nodes in a specific element
    	    J_el[j] = np.linalg.det(gll.Jacobian2D(dN_local[j,:,:],gll_coords_el))
	
        #We are now ready to construct the element matrices
        Fe = force_mat(f[gll_connect[i]],J_el,W)
        #Constructing global mass matrix
        Fg += l2g.local2global(Fe,Fg,gll_connect,[i])

    Fg_col = np.diag(Fg)
    return Fg_col




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
    """dfnalskdnjflakjsdf lkjc fasdlkj

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
