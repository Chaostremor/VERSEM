"""This is a function library for reshape operation and diagonal Matrix operation

Author: Fan Wu

"""

import numpy as np

def reshape(M, C, K, dim):
    """The function reorganize mass matrix and stiff matrix for 
    second order ODEs which uses first order solver

    :param M: original mass matrix, ``numpy`` [dim]x[dim] array.
    :param C: original absorbing matrix, ``numpy`` [dim]x[dim] array.
    :param K: original stiffness matrix, ``numpy`` [dim]x[dim] array.
    :param dim: dimension of the matrix, ``int``.

    :rtype: first term is the reshaped mass matrix, ``numpy`` [2xdim]x[2xdim] array.
            second term is the reshaped stiffness matrix, ``numpy`` [2xdim]x[2xdim] array.

    """
    if dim == 1:
        if C is None:
            C = 0
        M_out = np.array([[1,0],[0,M]])
        K_out = np.array([[0,-1],[K,C]])
    else:
        if C is None:
            C = np.zeros([dim, dim])
        I = np.eye(dim)
        O = np.zeros([dim,dim])
        M1 = np.concatenate((I,O),axis=1)
        M2 = np.concatenate((O,M),axis=1)
        M_out = np.concatenate((M1,M2))

        K1 = np.concatenate((O,-I),axis=1)
        K2 = np.concatenate((K, C), axis=1)
        K_out = np.concatenate((K1,K2))

    return M_out, K_out

def reshape_f(fn, dim):
    """The function reorganize force vector for second order 
    ODEs which uses first order solver

    :param fn: original force vector, ``numpy`` [dim] array.
    :param dim: dimension of the vector.

    :rtype: reshaped force vector, ``numpy`` [2xdim] array.

    """
    if dim == 1:
        return np.array([0,fn])
    else:
        Of = np.zeros(dim)
        fn_out = np.concatenate((Of,fn))
        return fn_out


def normalize(M, K, dim):
    """The function do the calculation M^(-1)*K with diagonal mass matrix m and stiffness matrix K.

    :param M: mass matrix, ``numpy`` [dim]x[dim] array.
    :param K: stiffness matrix, ``numpy`` [dim]x[dim] array.
    :param dim: dimension of the matrix.

    :rtype: M^(-1)*K, ``numpy`` [dim]x[dim] array.

    """

    kk = np.zeros([dim,dim])
    if dim == 1:
        return K / M
    for i in range(dim):
        kk[i,:] = K[i,:] / M[i,i]
    return kk

def normalize2(M,K,dim):
    """The function do the calculation K*M^(-1) with 
    diagonal mass matrix m and stiffness matrix K.

    :param M: mass matrix, ``numpy`` [dim]x[dim] array.
    :param K: stiffness matrix, ``numpy`` [dim]x[dim] array.
    :param dim: dimension of the matrix.

    :rtype: K*M^(-1), ``numpy`` [dim]x[dim] array.

    """

    kk = np.zeros([dim,dim])
    if dim == 1:
        return K / M
    for i in range(dim):
        kk[:,i] = K[:,i] / M[i,i]
    return kk


def normalize_f(M, fn, dim):
    """The function do the calculation M^(-1)*f with diagonal 
    mass matrix m and force vector fn.

    :param M: mass matrix, ``numpy`` [dim]x[dim] array.
    :param fn: force vector, ``numpy`` [dim] array.
    :param dim: dimension of the matrix.

    :rtype: M^(-1)*fn, ``numpy`` [dim] array.

    """
    p = np.zeros(dim)
    if dim == 1:
        return fn / M
    for i in range(dim):
        p[i] = fn[i] / M[i,i]
    return p


def diagmul_f(M,f,dim):
    """The function do the calculation M*f with diagonal 
    ass matrix m and force vector fn.

    :param M: mass matrix, ``numpy`` [dim]x[dim] array.
    :param f: force vector, ``numpy`` [dim] array.
    :param dim: dimension of the matrix.

    :rtype: M*f, ``numpy`` [dim] array.

    """
    p = np.zeros(dim)
    if dim == 1:
        return M * f
    for i in range(dim):
        p[i] = f[i]*M[i,i]
    return p


