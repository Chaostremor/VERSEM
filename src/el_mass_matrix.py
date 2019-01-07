import numpy as np
import src.gll_library as gll
import src.loc2glob as l2g
#######################################################################
###            Constructing Element Mass Matrices                  ####
#######################################################################

def el_mass_mat(rho,J,W):
    """.. function:: el_mass_mat(ngll_total,rho,J,W)

    Computes the Elemetal Mass Matrix M_e for each element

    :param rho: flattened density matrix [1]x[total Number of GLL points]
              (``numpy``)

    :param J: flattened Jacobian determinant matrix [1]x[total Number of GLL points]
              (``numpy``)

    :param W: flattened weight matrix [1]x[total Number of GLL points]
              (``numpy``)
    
    :rtype: ``numpy`` [ngll_total]x[ngll_total] array

    The description of the Jacobian can be found on the theory
    documentation.
    """
    #Retrieves the number of gll points in an element
    ngll_el = len(rho)

    #the element mass matrix
    M = np.zeros([ngll_el,ngll_el])
    for i in range(ngll_el):
        M[i,i] = rho[i]*J[i]*W[i]	#basically a diagonal matrix because of the nature of Shape Functions

    return M

def glob_el_mass_mat(gll_coordinates,gll_connect,rho,dN_local,W):
    """.. function:: glob_el_mass_mat(gll_coordinates,gll_connect,rho,dN_local,W)

    Computes the Global Mass Matrix Mg

    :param gll_coordinates: ``numpy`` array of size [ngll_total]x[2] containing the coordinates of all the gll points

    :param gll_connect: ``numpy`` array of size [el_no]x[ngll_el]. Contains the global indexing of gll nodes.

    :param rho: flattened density matrix [1]x[ngll_el]
              (``numpy``)

    :param dN_local: Local derivative of shape functions at each gll point in an element. `numpy`` array of size [total ngll]x[2]x[total ngll]

    :param W: flattened weight matrix [1]x[ngll_el]
              (``numpy``)
    
    :rtype: ``numpy`` [ngll_total]x[ngll_total] array

    """
    
    #Retrieving the number of elements in the domain and the number of gll points per element
    el_no = len(gll_connect)
    ngll_el = len(gll_connect[0])	#Assuming number of gll points per element in constant

    Mg = np.zeros([len(gll_coordinates),len(gll_coordinates)])

    M = np.zeros([2*len(gll_coordinates),2*len(gll_coordinates)])
    
    for i in range(el_no):    #el_no is the total number of elements
        gll_coords_el = gll_coordinates[gll_connect[i]]
        J_el = np.zeros(len(gll_coords_el))
        for j in range(len(gll_coords_el)):
            #Jacobian for the nodes in a specific element
    	    J_el[j] = np.linalg.det(gll.Jacobian2D(dN_local[j,:,:],gll_coords_el))
	
        #We are now ready to construct the element matrices
        Me = el_mass_mat(rho[gll_connect[i]],J_el,W)
        #Constructing global mass matrix
        Mg += l2g.local2global(Me,Mg,gll_connect,[i])

    M[0:len(gll_coordinates),0:len(gll_coordinates)] = Mg
    M[len(gll_coordinates):2*len(gll_coordinates),len(gll_coordinates):2*len(gll_coordinates)] = Mg

    return M
    
