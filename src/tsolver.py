import numpy as np

import numba as nb


def euler_explicit(xn, dt, f, t):
    # M is the mass matrix, which is diagnol, K is the stiff matrix
    # f is the force term, dt is the time step, xn is the previous step
    fxn = f(xn,t)
    xn1 = xn + dt * fxn
    return xn1,fxn

def euler_adjust(xn, dt, f, t):
    fxn = f(xn,t)
    kn = xn + dt * fxn
    fkn = f(kn,t+dt)
    xn1 = xn + dt * (fxn + fkn)/2
    return xn1,fxn

def rk4(xn, dt, f, t):
    k1 = f(xn,t)
    k2 = f(xn+dt/2*k1,t+dt/2)
    k3 = f(xn+dt/2*k2,t+dt/2)
    k4 = f(xn+dt*k3,t+dt)
    xn1 = xn + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
    return xn1,k1

def ab2(xn, dt, f, t, cache):
    fxn = f(xn,t)
    fxn_1 = cache
    xn1 = xn + dt/2 * (3*fxn - fxn_1)
    return xn1, fxn

def ab3(xn, dt, f, t, cache):
    fxn = f(xn,t)
    fxn_1 = cache[0,:]
    fxn_2 = cache[1,:]
    xn1 = xn + dt/12 * (23*fxn - 16*fxn_1 + 5*fxn_2)
    return xn1, fxn

def ab4(xn, dt, f, t, cache):
    fxn = f(xn,t)
    fxn_1 = cache[0,:]
    fxn_2 = cache[1,:]
    fxn_3 = cache[2,:]
    xn1 = xn + dt/24 * (55*fxn - 59*fxn_1 + 37*fxn_2 - 9*fxn_3)
    return xn1, fxn

## central difference newmark, dt^2 approximation, estimate then adjust
## we restrict to C = diag here


## xn: only xn
# cache first conlumn x', second column x''

def newmark(M ,K ,f,t ,dt ,un, nstep ,cache, gamma):
    if nstep == 0:
        fn = f(t)
        v0 = cache
        a0 = normalize_f(M,fn,len(fn)) - np.dot(normalize(M,K,len(fn)),un)
        cache = np.concatenate(([v0],[a0]))
        un1, cache = newmark(M ,K ,f,t ,dt ,un, 1 ,cache, gamma)
        return un1, cache
    else:
        fn = f(t+dt)
        vn = cache[0,:]
        an = cache[1,:]
        un1 = un + dt*vn + 0.5*dt*dt*an
        fn1 = fn - np.dot(K, un1)
        an1 = normalize_f(M,fn1,len(fn1))
        vn1 = vn + (1-gamma)*dt*an + gamma*dt*an1
        cache = np.concatenate(([vn1],[an1]))
        return un1, cache

def newmark2(K,f,t,dt,un,nstep,cache):

    if nstep == 0:
        un1 = un + cache*dt
        cache = un
        return un1, cache
    else:
        un = un
        un_1 = cache
        fn = f(t)
        ff = np.dot(K,un)
        un1 = (dt*dt)*(fn-ff) + 2*un - un_1
        cache = un
        return un1, cache





def reshape(M, C, K, dim):
    if dim == 1:
        M_out = np.array([[1,0],[0,M]])
        K_out = np.array([[0,-1],[K,0]])
    else:
        I = np.eye(dim)
        O = np.zeros([dim,dim])
        M1 = np.concatenate((I,O),axis=1)
        M2 = np.concatenate((O,M),axis=1)
        M_out = np.concatenate((M1,M2))

        K1 = np.concatenate((O,-I),axis=1)
        K2 = np.concatenate((K, O), axis=1)
        K_out = np.concatenate((K1,K2))

    return M_out, K_out

def reshape_f(fn,dim):
    if dim == 1:
        return np.array([0,fn])
    else:
        Of = np.zeros(dim)
        fn_out = np.concatenate((Of,fn))
        return fn_out


def normalize(M, K, dim):
    kk = np.zeros([dim,dim])
    if dim == 1:
        return K / M
    for i in range(dim):
        kk[i,:] = K[i,:] / M[i,i]
    return K


def normalize_f(M, fn, dim):
    p = np.zeros(dim)
    if dim == 1:
        return fn / M
    for i in range(dim):
        p[i] = fn[i] / M[i,i]
    return p


def diagmul_M(M,f,dim):
    p = np.zeros(dim)
    if dim == 1:
        return M * f
    for i in range(dim):
        p[i] = f[i]*M[i,i]
    return p

def normalize2(M,K,dim):
    kk = np.zeros([dim,dim])
    if dim == 1:
        return K / M
    for i in range(dim):
        kk[:,i] = K[:,i] / M[i,i]
    return kk













