"""This is a function library for reshape operation and diagonal Matrix operation

Author: Fan Wu

"""

import numpy as np

def reshape(M, C, K, dim):
    """.. function:: reshape(M, C, K, dim)

    The function reorganize mass matrix and stiff matrix for second order ODEs which uses first order solver

    :param M: original mass matrix.
    :param C: original absorbing matrix.
    :param K: original stiffness matrix.
    :param dim: dimension of the matrix.

    :rtype: M_out is the reorganized mass matrix and K_out is the reorganized stiffness matrix, since we turn a
            second order ODEs into first order, there is no absorbing matrix.

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
    """.. function:: reshape_f(fn, dim)

    The function reorganize force vector for second order ODEs which uses first order solver

    :param fn: original force vector.
    :param dim: dimension of the vector.

    :rtype: fn_out is the reorganized force vector.

    """
    if dim == 1:
        return np.array([0,fn])
    else:
        Of = np.zeros(dim)
        fn_out = np.concatenate((Of,fn))
        return fn_out


def normalize(M, K, dim):
    """.. function:: normalize(M, K, dim)

    The function do the calculation M^(-1)*K with diagonal mass matrix m and stiffness matrix K.

    :param M: mass matrix.
    :param K: stiffness matrix.
    :param dim: dimension of the matrix.

    :rtype: kk is M^(-1)*K

    """

    kk = np.zeros([dim,dim])
    if dim == 1:
        return K / M
    for i in range(dim):
        kk[i,:] = K[i,:] / M[i,i]
    return kk

def normalize2(M,K,dim):
    """.. function:: normalize2(M, K, dim)

    The function do the calculation K*M^(-1) with diagonal mass matrix m and stiffness matrix K.

    :param M: mass matrix.
    :param K: stiffness matrix.
    :param dim: dimension of the matrix.

    :rtype: kk is K*M^(-1)

    """

    kk = np.zeros([dim,dim])
    if dim == 1:
        return K / M
    for i in range(dim):
        kk[:,i] = K[:,i] / M[i,i]
    return kk


def normalize_f(M, fn, dim):
    """.. function:: normalize_f(M, fn, dim)

    The function do the calculation M^(-1)*f with diagonal mass matrix m and force vector fn.

    :param M: mass matrix.
    :param fn: force vector.
    :param dim: dimension of the matrix.

    :rtype: p is M^(-1)*fn

    """
    p = np.zeros(dim)
    if dim == 1:
        return fn / M
    for i in range(dim):
        p[i] = fn[i] / M[i,i]
    return p


def diagmul_f(M,f,dim):
    """.. function:: diagmul_f(M, f, dim)

    The function do the calculation M*f with diagonal mass matrix m and force vector fn.

    :param M: mass matrix.
    :param f: force vector.
    :param dim: dimension of the matrix.

    :rtype: p is M*f

    """
    p = np.zeros(dim)
    if dim == 1:
        return M * f
    for i in range(dim):
        p[i] = f[i]*M[i,i]
    return p


