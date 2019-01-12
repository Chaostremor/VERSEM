"""Force term implementation script"""

import numpy as np
from scipy import interpolate
import src.gll_library as gll
import src.loc2glob as l2g

# =================== Main Routine ====================================== #
# Takes in all arguments and just uses function from this script to 
# compute the complete necessary force vector for the time stepping.
#

def F(force_term,force_location,t,source_time_function,gll_coordinates,gll_connect,dN_local,W):
    """Computes the whole force vector interpolator corresponding to a solution 
    diplacement vector [ux1,uy1,ux2,uy2,ux3,.....]. Given a time argument within 
    the range of ``t`` the interpolator returns a force vector.


    :param force_term: 2 element ``numpy`` array including the force 
                       direction in 2D

    :param force_location: 2 element ``numpy`` array including the force
                           location in 2D

    :param t: ``nt`` element ``numpy`` array contaning the time steps of
              the source time function.

    :param source_time_function: ``nt`` element ``numpy`` array contaning 
                                   the function values corresponding to t.

    :param gll_coordinates: [ngll_tot]x[dim] ``numpy`` array containing
                            GLL node coordinates.

    :param dN_local: [dim]x[ngll_el] ``numpy`` array that contains the local
                      derivative matrix

    :param W: [ngll_el] ``numpy`` array containing the integration weights

    
    :rtype: Interpolator of the form F(t). It takes in one argument that is 
            time and spits out a [2*ngll] ``numpy`` array corresponding to
            a solution diplacement vector [ux1,uy1,ux2,uy2,ux3,.....].

    """
    # Number of total nodes
    ngll_total = len(gll_coordinates)
    
    # Creating two vectors of the same length as GLL nodes containing the 
    # force in x and y direction at the indeces corresponding to the closest
    # GLL coordinates
    # Uses the genforce function in the 
    Fx_loc,Fy_loc = genforce(force_term,force_location,gll_coordinates)

    # Computing the force vectors for x and y direction
    Fx = global_force_vector(gll_coordinates,gll_connect,Fx_loc,dN_local,W)
    Fy = global_force_vector(gll_coordinates,gll_connect,Fy_loc,dN_local,W)

    # Assigning the vectors in the two directions such that they accommodate
    # alternating displacements [ux1,uy1,ux2,uy2,ux3,.....]
    F_tot = np.zeros(2*ngll_total)
    F_tot[0:2*ngll_total-1:2] = Fx
    F_tot[1:2*ngll_total:2]   = Fy

    # Length of the function 
    N = len(source_time_function)
    
    # Create repeat matrix
    repF = np.repeat(F_tot.reshape([len(F_tot),1]),N,axis=1)

    # repeat vector
    fp = repF*source_time_function
    
    return interpolate.interp1d(t,fp)

# ====================================================================#




# -------- Constructing Element Force Matrices -----------------------#

def element_force_matrix(f,J,W):
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

# -------- Constructing Global Force Vector ----------------------- #

def global_force_vector(gll_coordinates,gll_connect,f,dN_local,W):
    """Computes the Global Force Vector from a given GLL located 
    force vector.

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
        Fe = element_force_matrix(f[gll_connect[i]],J_el,W)
        #Constructing global mass matrix
        Fg += l2g.local2global(Fe,Fg,gll_connect,[i])

    Fg_col = np.diag(Fg)


    return Fg_col




def genforce(force_term,force_location,gll_coordinates):
    """This function implements the force term of a point source for the 
    seismic wave equation. It creates two vectors of the same length as 
    GLL nodes containing the force in x and y direction at the indeces 
    corresponding to the closest GLL coordinates

    :param force_term: 2 element ``numpy`` array with forces in 2D

    :param force_loc: 2 element ``numpy`` numpy array with force location 

    :param gll_coordinates: [ngll] ``numpy`` array containg the gll
                            coordinates

    :rtype: tuple of two [ngll] ``numpy`` array. One vector for the force
            in x direction and one for y direction

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

if __name__ == "__main__":
    pass
