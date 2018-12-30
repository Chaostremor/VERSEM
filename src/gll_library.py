import numpy as np


#######################################################################
###                 Lagrange Polynomial                            ####
#######################################################################



class LagrangePoly1D(object):                                               
    """Callable lagrange interpolator object.                                      
                                                                        
    Example usage: to construct the lagrange polynomial 
    l_0^1(x) =  and evaluate p(3):                                               
                                                                        
    p = Polynomial(1,[-1 1])                                           
                                                                        
    p(3)                                                                
                                                                        
    """
    pass

    def __init__(self, degree, xi):                                         
        """ """
         
        self._degree = degree
        self._xi = xi
                                                                        
    def _f(self,i):                                                     
        """Calculation of the i-th lagrange polynomial of 
        degree <self._degree>.
        """
        fac = 1
        for j in range(-1, self._degree):
            if j != i:
                fac = fac * ((x - xi[j + 1]) / (xi[i + 1] - xi[j + 1]))
        return fac
        

    def _df(self,x):
        """derivative of the lagrange polynomial"""
        pass
                                                   

    # Instances of classes that have a defined __call__ method are      
    # themselves callable, as if they were functions                    
    def __call__(self, x):                                              
        return self._f(x)          



def lagrange(i, x, xi):
    """
    Function that evaluates Lagrange polynomial of order N (len(xi)-1) 
    and polynomial i [0, N-1] at location x at given collocation points
    xi (not necessarily the GLL-points).
    """
    fac = 1
    
    for j in range(0, len(xi)):
        if j != i:
            fac = fac * ((x - xi[j]) / (xi[i] - xi[j]))
    return fac




def lagrange2D(i,x,xi,j,y,eta):
    """lagrange2D(i,x,xi,j,y,eta)
    
    calculates the 2D lagrange polynomial given collocation point sets
    xi and eta the coordinate x,y, and the polynomial numbers i,j
    
    i, x, xi are the parameters for the polynomial in one direction and 
    j, y, eta the parameters for the other direction in orthogonal space

    """

    return lagrange(i,x,xi)*lagrange(j,y,eta)   




#######################################################################
###                Legendre Polynomials                             ###
#######################################################################

def legendre(i,x,xi):
    """
    Returns the value of Legendre Polynomial P_i(x) at location x given
    collocation points xi (not necessarily GLL points) and polynomial
    number i.
    extremely simple algorithm.
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
    """
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

#######################################################################
###                 Next function                                  ####
#######################################################################
