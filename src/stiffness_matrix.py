import numpy as np
import src.gll_library as gll
import src.loc2glob as l2g


#######################################################################
###            Constructing Element Stiffness Matrices             ####
#######################################################################

def element_stiffness_matrix(gll_coords_el,dim,ngll_el,dN_local,comp,W,lmd,mu):
    """Computes the Elemetal Stiffness Matrix in three parts A,B,C for each element.

    :param gll_coords: the global coordinates for the gll points in a particular element.

    :param dim: the dimensionaloty of our system. Currently 2D. So dim=2.

    :param dN_local: Local derivative of shape functions at each gll point in an element. `numpy`` array of size [total ngll]x[2]x[total ngll]

    :param comp: the component of u that we are solving for when calling this function

    :param W: flattened weight matrix [1]x[total Number of GLL points] (``numpy``)

    :param lmd: the flattened :math:`\\lambda` array [1]x[total Number of GLL points] (``numpy``)

    :param mu: the flattened :math:`\\mu` array [1]x[total Number of GLL points] (``numpy``)
    
    :rtype: A and B are ``numpy`` [dim]x[ngll_el]x[ngll_el] array and C is ``numpy`` [ngll_el]x[ngll_el] array

    The description of the Jacobian can be found on the theory
    documentation.
    """
    
    # Initialiying empty matrices
    A = np.zeros([dim,ngll_el,ngll_el])
    B = np.zeros([dim,ngll_el,ngll_el])
    C = np.zeros([ngll_el,ngll_el])

    # Initialize empty matrix for storing the determinant of the 
    # Jacobian at each gll point in an element
    J_el = np.zeros(len(gll_coords_el))	

    # Initialize global Derivative
    global_der = np.zeros(np.shape(dN_local))

    for j in range(len(gll_coords_el)):
        #Jacobian for the nodes in a specific element
        Jacob = gll.Jacobian2D(dN_local[j,:,:],gll_coords_el)
        J_el[j] = np.linalg.det(Jacob)

        # Computing the global derivatives
        global_der[j,:,:] = gll.global_derivative(Jacob,dN_local[j,:,:])
    

    # Sum following the Notation in the documentation strictly
    for l in range(ngll_el):
        for m in range(ngll_el):
            for r in range(dim):
                for k in range(ngll_el):
                    A[r,l,m] += -(global_der[k,comp,l]*lmd[k]*global_der[k,r,m]*(J_el[k])*W[k])
                    B[r,l,m] += -(global_der[k,r,l]*mu[k]*global_der[k,comp,m]*(J_el[k])*W[k])
                    C[l,m] += -(global_der[k,r,l]*mu[k]*global_der[k,r,m]*(J_el[k])*W[k])

    return A,B,C



def global_stiffness_matrix(gll_coordinates,gll_connect,dN_local,W,dim,lmd,mu):
    """Computes the Global Mass Matrix Mg

    :param gll_coordinates: ``numpy`` array of size [ngll_total]x[dim] containing the coordinates of all the gll points

    :param gll_connect: ``numpy`` array of size [el_no]x[ngll_el]. Contains the global indexing of gll nodes.

    :param dN_local: Local derivative of shape functions at each gll point in an element. `numpy`` array of size [total ngll]x[2]x[total ngll]

    :param W: flattened weight matrix [1]x[ngll_el]
              (``numpy``)

    :param comp: depends on which component of u we are solving for

    :param dim: the dimensionality of our system. For now its 2.
    
    :param lmd: the flattened `\\lambda` array [1]x[total Number of GLL points] (``numpy``)

    :param mu: the flattened `\\mu` array [1]x[total Number of GLL points] (``numpy``)

    :rtype: Ag and Bg are ``numpy`` [dim]x[ngll_total]x[ngll_total] array and C is ``numpy`` [ngll_total]x[ngll_total] array

    """

    #Retrieving the number of elements in the domain and the number of gll points per element
    ngll_total = len(gll_coordinates)
    el_no = len(gll_connect)
    ngll_el = len(gll_connect[0])	#Assuming number of gll points per element in constant

    #Initializing the three parts of the global stiffness matrix
    Ag0 = np.zeros([dim,len(gll_coordinates),len(gll_coordinates)])
    Bg0 = np.zeros([dim,len(gll_coordinates),len(gll_coordinates)])
    Ag1 = np.zeros([dim,len(gll_coordinates),len(gll_coordinates)])
    Bg1 = np.zeros([dim,len(gll_coordinates),len(gll_coordinates)])
    
    Cg0 = np.zeros([len(gll_coordinates),len(gll_coordinates)])
    Cg1= np.zeros([len(gll_coordinates),len(gll_coordinates)])

    #Looping over all elements to create the global matrices
    for i in range(el_no):

        #The stiffness matrix => Computed only once.
        #We split up the computation of the elemental stiffness matrix into three parts
        #Part A

        #Retrieving the gll coordinates corresponding to the current element
        gll_coords_el = gll_coordinates[gll_connect[i]]
        mu_el = mu[gll_connect[i]]
        lmd_el = lmd[gll_connect[i]]

        #Obtaining the local matrices
        comp = 0
        A0,B0,C0 = element_stiffness_matrix(gll_coords_el,dim,ngll_el,dN_local,comp,W,lmd_el,mu_el)
        comp = 1
        A1,B1,C1 = element_stiffness_matrix(gll_coords_el,dim,ngll_el,dN_local,comp,W,lmd_el,mu_el)

        
        for j in range(dim):
            # comp 0
            Ag0[j,:,:] += l2g.local2global(A0[j,:,:],Ag0[j,:,:],gll_connect,[i])
            Bg0[j,:,:] += l2g.local2global(B0[j,:,:],Bg0[j,:,:],gll_connect,[i])
            # comp 1
            Ag1[j,:,:] += l2g.local2global(A1[j,:,:],Ag1[j,:,:],gll_connect,[i])
            Bg1[j,:,:] += l2g.local2global(B1[j,:,:],Bg1[j,:,:],gll_connect,[i])
            

        Cg0 += l2g.local2global(C0,Cg0,gll_connect,[i])
        Cg1 += l2g.local2global(C1,Cg1,gll_connect,[i])
        


    ## Reassign locations for the Stiffness components
    # Initialize empty Stiffness matrix
    K = np.zeros([2*ngll_total,2*ngll_total])
    # Redistributing Top-Left
    K[0:2*ngll_total:2,0:2*ngll_total:2] = Ag0[0,:,:] + Bg0[0,:,:] + Cg0
    # Redistributing Top-right
    K[0:2*ngll_total:2,1:2*ngll_total:2] = Ag0[1,:,:] + Bg0[1,:,:]
    # Redistributing Bottom-left
    K[1:2*ngll_total:2,0:2*ngll_total:2] = Ag1[0,:,:] + Bg1[0,:,:] 
    # Redistributing Bottom-right
    K[1:2*ngll_total:2,1:2*ngll_total:2] = Ag1[1,:,:] + Bg1[1,:,:] + Cg1
    
    return K



        

