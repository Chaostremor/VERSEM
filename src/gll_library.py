"""This is a function library for GLL interpolation and quadrature

Author: Lucas Sawade

"""

import numpy as np


#######################################################################
###                 Lagrange Polynomial                            ####
#######################################################################

def lagrange(i, x, xi):
    """The algorithm follows the definition of the Lagrange Polynomial
    strictly.

    :param i: Polynomial number.
    :param x: location of evaluation.
    :param xi: numpy 1D array of collocation points.

    :rtype: value of Lagrange polynomial of order N (len(xi)-1)  and 
            polynomial number i [0, N-1] at location x at given 
            collocation points xi (not necessarily the GLL-points).

    """

    # Necessary for generic dimensional code
    """
    if len(xi)==1:
        return 1.
    """

    fac = 1
    
    for j in range(0, len(xi)):
        if j != i:
            fac = fac * ((x - xi[j]) / (xi[i] - xi[j]))
    return fac




def lagrange2D(i,x,xi,j,y,eta):
    """Computes the lagrange polynomial in two dimensions 
    
    :param i: polynomial number 0-th dimension
    :param x: location of evalutation 0-th dimension
    :param xi: collocation points 0-th dimension

    :param j: polynomial number in 1st dimension
    :param y: location of evalutation in 1st dimension
    :param eta: collocation points in 1st dimension

    :rtype: returns value of Lagrange polynomial of order N (len(xi)-1) 
            and polynomial number i [0, N-1] at location x at given 
            collocation points xi (not necessarily the GLL-points).
            calculates the 2D lagrange polynomial given collocation 
            point sets xi and eta the coordinate x,y, and the polynomial 
            numbers i,j
    
    """

    return lagrange(i,x,xi)*lagrange(j,y,eta)   


def lagrange1st(i,x,xi):
    """Computes the 1st derivative of the lagrange polynomial in 1D.

    :param i: polynomial number
    :param x: location of evalutation
    :param xi: collocation points

    
    :rtype: return the value of the derivative of the lagrange 
            polynomial in 1D at point x, of polynomial number i, with 
            collocation points xi.
        
    The derivative of the Lagrange polynomial is the Lagrange polynomial 
    multiplied with the Legendre Polynomial. Through simplification which
    is described in detail in the documentation we can omit one 
    numerator when k equals j.

    """
    
    sum = 0
    
    # Loop over Legendre summing term
    for k in range(len(xi)):
        # Check whether k is equal to i
        if k != i:

            # Reset fac for every term of the sum
            fac = 1

            # Loop over Lagrange multiplication term
            for j in range(len(xi)):
                if j != i:
                    # When k equals j we have cancelling numerator
                    if j == k:
                        fac = fac *      1       / (xi[i] - xi[j])
                    else:
                        fac = fac * ((x - xi[j]) / (xi[i] - xi[j]))
            
            # Note that we don't need a legendre polynomial term, since it 
            # cancels out. See Documentation for details.
            sum = sum + fac
    
    return sum 
    
    
            
def lagrange1st2D(i,x,xi,j,y,eta,d):
    """Computes the 1st derivative of the 2D Lagrange polynomial.

    :param i: polynomial number in 0-th dimension
    :param x: location of evalutation in 0-th dimension
    :param xi: collocation points in 0-th dimension

    :param j: polynomial number in 1st dimension
    :param y: location of evalutation in 1st dimension
    :param eta: collocation points in 1st dimension
    
    :param d: dimension of differentiation (0 or 2)

    :rtype: return the value of the derivative of the lagrange 
            polynomial in 1D at point x, of polynomial number i, with 
            collocation points xi.

    This function computes the 2D derivative of the lagrange polynomial 
    in dimension d (0 or 1). The collocation points in the first and 
    second dimension are xi and eta, respectively. The location is (x,y)
    and the numbers of the polynomials in the different dimensions are 
    i and j.

    """

    # derivative in xi dimension
    if d==0:
        return lagrange1st(i,x,xi)*lagrange(j,y,eta)
    # derivative in eta dimension
    elif d==1:
        return lagrange(i,x,xi)*lagrange1st(j,y,eta)


def lagrangeDerMat2D(x,xi,y,eta):
    """This function computes the 2D dervative matrix of the lagrange 
    polynomial in the form
        
    ::
    
         _         _     _                                          _
        |  dN/dxi   |   |  dN1/dxi  dN2/dxi  dN3/dxi  ...  dNn/dxi   |
        |           | = |                                            |
        |_ dN/deta _|   |_ dN1/deta dN2/deta dN3/deta ...  dNn/deta _|
    

    :param x: location in 0th dimension
    :param xi: ``numpy`` 1D array of collocation points in 0th dimension

    :param y: location in 1st dimension
    :param eta: ``numpy`` 1D array collocation points in 1st dimension
    
    :rtype: ``numpy`` 2x[len(xi)*len(eta)] array 
    
    For the numbering of the shape functions see documentation.

    """
    
    # Initialize numpy array 
    dN_dxi = np.zeros([2,len(xi)*len(eta)])

    ### Compute derivatives of the k-th dimension
    # Counter to set index in derivative matrix.
    for k in range(2):
        counter = 0
        # Loop over y dimension
        for j in range(len(eta)):
            # Loop over x direction
            for i in range(len(xi)):
                dN_dxi[k,counter] = lagrange1st2D(i,x,xi,j,y,eta,k)
                counter += 1
    
    return dN_dxi


#######################################################################
###                The Jacobian Matrix                              ###
#######################################################################


def Jacobian(dN,x):
    """Computes the Jacobian Matrix [dx/dxi] between the global and the 
    local element coordinates

    :param dN: shape function derivative matrix (``numpy``)
               [dim]x[Number of GLL points]
    :param x: global node coordinates [dim]x[Number of GLL points]
              (``numpy``)
    
    :rtype: ``numpy`` [dim]x[dim] array

    The description of the Jacobian can be found on the theory
    documentation.

    """
    pass

def Jacobian2D(dN,x):
    """Computes the Jacobian Matrix [dx/dxi] between the global and the 
    local element coordinates

    :param dN: shape function derivative matrix (``numpy``)
               [dim]x[total Number of GLL points]
    :param x: global node coordinates [dim]x[total Number of GLL points]
              (``numpy``)
    
    :rtype: ``numpy`` [dim]x[dim] array

    The description of the Jacobian can be found on the theory
    documentation.
    """

    return np.matmul(dN,x)


def dN_local_2D(xi,eta):
    """Computes the flattened Jacobian for each GLL point. Numbering as 
    described in the documentation.

    :param xi: ``numpy`` 1D array of collocation points in xi direction
    :param eta: ``numpy`` 1D array of collocation points in eta direction

    :rtype: ``numpy`` array of size [total ngll]x[2]x[total ngll]

    Needs to be computed once the beginning since it is the same for all
    elements. It's the derivative of each shape function in each 
    dimension computed at each gll point of the element.

    
    """
    # Inititalize the array
    dN_local = np.zeros([len(xi)*len(eta),2,len(xi)*len(eta)])

    # We fill this once and for all as this shall stay the same over all 
    # elements
    count = 0
    for j in range(len(eta)):
        for i in range(len(xi)):
            dN_local[count,:,:] = lagrangeDerMat2D(xi[i],xi,\
                                                       eta[j],eta)
            count += 1

    return dN_local



#######################################################################
###                Global Derivatives                             ###
#######################################################################

def global_derivative(jacob,dNdxi):
    """This function computes the global derivatives from the local shape
    function derivates.

    :param jacob: the [dim]x[dim] ``numpy`` array containing the 
                  Jacobian transformation matrix.

    :param dNdxi: local shapefunction derivative ``numpy`` array of size
                  [dim]x[total Number of GLL points]

    :rtype: global shapefunction derivative ``numpy`` array of size
            [dim]x[total Number of GLL points].

    The global derivative is needed for the GLL quadrature on element
    level.
    
    The algorithm is simple since
    
    .. math::
        
        \\frac{d}{dx} = J^{-1} \\frac{d}{d\\xi}

    """
    # One liner since it's a simple inversion and matrix multiplication
    #                    J^(-1)      * dNdxi 
    return np.matmul(np.linalg.inv(jacob),dNdxi)



#######################################################################
###                Legendre Polynomials                             ###
#######################################################################

def legendre(i,x,xi):
    """Computes the legendre polynomial.
    
    :param i: number of polynomial
    :param x: location of evalutation
    :param xi: ``numpy`` 1D array of collocation points.

    :rtype: returns the value of Legendre Polynomial P_i(x) at location 
            x given collocation points xi (not necessarily GLL points) 
            and polynomial number i. 
    
    Extremely simple algorithm.
    
    """

    sum = 0

    for j in range(len(xi)):
        if j != i:
            sum = sum + 1/(x - xi[j])

    return sum



#######################################################################
###                 GLL - Points and Weights                       ####
#######################################################################

def gll_pw(N):
    """Outputs the Gauss-Lobatto-Legendre integration weights as a
    function of polynomial degree.
    
    :param N: Polynomial degree
    
    :rtype: 1x2 tuple of 1D ``numpy`` 1x(N+1) arrays, the first 
            containing the collocation points and the second containing
            quadrature weights

    Takes in polynomial degree and returns the (N+1) points and weights
    Returns GLL (Gauss Lobato Legendre module with collocation points and
    weights)
    
    """

    # Initialization of integration weights and collocation points
    # [xi, weights] =  gll(N)
    # Values taken from Diploma Thesis Bernhard Schuberth
    if N == 1:
        xi = [-1,1]
        weights = [1,1]
    elif N == 2:
        xi = [-1.0, 0.0, 1.0]
        weights = [0.33333333, 1.33333333, 0.33333333]
    elif N == 3:
        xi = [-1.0, -0.447213595499957, 0.447213595499957, 1.0]
        weights = [0.1666666667, 0.833333333, 0.833333333, 0.1666666666]
    elif N == 4:
        xi = [-1.0, -0.6546536707079772, 0.0, 0.6546536707079772, 1.0]
        weights = [0.1, 0.544444444, 0.711111111, 0.544444444, 0.1]
    elif N == 5:
        xi = [-1.0, -0.7650553239294647, -0.285231516480645, 0.285231516480645,
              0.7650553239294647, 1.0]
        weights = [0.0666666666666667,  0.3784749562978470,
                   0.5548583770354862, 0.5548583770354862, 0.3784749562978470,
                   0.0666666666666667]
    elif N == 6:
        xi = [-1.0, -0.8302238962785670, -0.4688487934707142, 0.0,
              0.4688487934707142, 0.8302238962785670, 1.0]
        weights = [0.0476190476190476, 0.2768260473615659, 0.4317453812098627,
                   0.4876190476190476, 0.4317453812098627, 0.2768260473615659,
                   0.0476190476190476]
    elif N == 7:
        xi = [-1.0, -0.8717401485096066, -0.5917001814331423,
              -0.2092992179024789, 0.2092992179024789, 0.5917001814331423,
              0.8717401485096066, 1.0]
        weights = [0.0357142857142857, 0.2107042271435061, 0.3411226924835044,
                   0.4124587946587038, 0.4124587946587038, 0.3411226924835044,
                   0.2107042271435061, 0.0357142857142857]
    elif N == 8:
        xi = [-1.0, -0.8997579954114602, -0.6771862795107377,
              -0.3631174638261782, 0.0, 0.3631174638261782,
              0.6771862795107377, 0.8997579954114602, 1.0]
        weights = [0.0277777777777778, 0.1654953615608055, 0.2745387125001617,
                   0.3464285109730463, 0.3715192743764172, 0.3464285109730463,
                   0.2745387125001617, 0.1654953615608055, 0.0277777777777778]
    elif N == 9:
        xi = [-1.0, -0.9195339081664589, -0.7387738651055050,
              -0.4779249498104445, -0.1652789576663870, 0.1652789576663870,
              0.4779249498104445, 0.7387738651055050, 0.9195339081664589, 1.0]
        weights = [0.0222222222222222, 0.1333059908510701, 0.2248893420631264,
                   0.2920426836796838, 0.3275397611838976, 0.3275397611838976,
                   0.2920426836796838, 0.2248893420631264, 0.1333059908510701,
                   0.0222222222222222]
    elif N == 10:
        xi = [-1.0, -0.9340014304080592, -0.7844834736631444,
              -0.5652353269962050, -0.2957581355869394, 0.0,
              0.2957581355869394, 0.5652353269962050, 0.7844834736631444,
              0.9340014304080592, 1.0]
        weights = [0.0181818181818182, 0.1096122732669949, 0.1871698817803052,
                   0.2480481042640284, 0.2868791247790080, 0.3002175954556907,
                   0.2868791247790080, 0.2480481042640284, 0.1871698817803052,
                   0.1096122732669949, 0.0181818181818182]
    elif N == 11:
        xi = [-1.0, -0.9448992722228822, -0.8192793216440067,
              -0.6328761530318606, -0.3995309409653489, -0.1365529328549276,
              0.1365529328549276, 0.3995309409653489, 0.6328761530318606,
              0.8192793216440067, 0.9448992722228822, 1.0]
        weights = [0.0151515151515152, 0.0916845174131962, 0.1579747055643701,
                   0.2125084177610211, 0.2512756031992013, 0.2714052409106962,
                   0.2714052409106962, 0.2512756031992013, 0.2125084177610211,
                   0.1579747055643701, 0.0916845174131962, 0.0151515151515152]
    elif N == 12:
        xi = [-1.0, -0.9533098466421639, -0.8463475646518723,
              -0.6861884690817575, -0.4829098210913362, -0.2492869301062400,
              0.0, 0.2492869301062400, 0.4829098210913362,
              0.6861884690817575, 0.8463475646518723, 0.9533098466421639,
              1.0]
        weights = [0.0128205128205128, 0.0778016867468189, 0.1349819266896083,
                   0.1836468652035501, 0.2207677935661101, 0.2440157903066763,
                   0.2519308493334467, 0.2440157903066763, 0.2207677935661101,
                   0.1836468652035501, 0.1349819266896083, 0.0778016867468189,
                   0.0128205128205128]
    else:
        raise NotImplementedError

    return np.array(xi), np.array(weights)

def flattened_weights2D(w_xi,w_eta):
    """This function flattens the weights into a 1D array.

    :param xi: ``numpy`` 1D array of weights corresponding to 
               collocation points in xi direction
    :param eta: ``numpy`` 1D array of weights corresponding to
                collocation points in eta direction

    :rtype: ``numpy`` array of size 1x[len(xi)*len(eta)]
    
    Numbering of the weights follows the numbering scheme in the 
    documentation.
    
    """

    # Initialize empty array of Weights
    W = np.zeros(len(w_xi)*len(w_eta))

    # Flattening the Weights
    counter = 0
    
    for j in range(len(w_eta)):
        for i in range(len(w_xi)):
            W[counter] = w_xi[i]*w_eta[j]
            counter += 1

    return W
    














#######################################################################
###                 Next function                                  ####
#######################################################################
