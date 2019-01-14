"""This is a function library for different explicit time schemes

Author: Fan Wu

"""

import numpy as np

import numba as nb


def euler_explicit(xn, dt, f, t):
    """.. function:: euler_explicit(xn, dt, f, t)

    The function do the time stepping of ODEs with explicit Euler method.

    :param xn: value of the vector at t.
    :param dt: time interval between step.
    :param f: force function f(x,t), given x and t return a vector.
    :param t: current time

    :rtype: xn1 is the value of the vector at time t+dt; fxn is value f(xn,t), the value will be stored is for multi-step method

    """

    fxn = f(xn,t)
    xn1 = xn + dt * fxn
    return xn1,fxn

def  rk2(xn, dt, f, t):
    """.. function:: rk2(xn, dt, f, t)

    The function do the time stepping of ODEs with second order Runge-Kutta.

    :param xn: value of the vector at t.
    :param dt: time interval between step.
    :param f: force function f(x,t), given x and t return a vector.
    :param t: current time

    :rtype: xn1 is the value of the vector at time t+dt; fxn is value f(xn,t), the value will be stored is for multi-step method

    """
    fxn = f(xn,t)
    xkn = xn + dt * fxn
    fkn = f(xkn,t+dt)
    xn1 = xn + dt * (fxn + fkn)/2
    return xn1,fxn

def rk4(xn, dt, f, t):
    """.. function:: rk4(xn, dt, f, t)

    The function do the time stepping of ODEs with fourth order Runge-Kutta.

    :param xn: value of the vector at t.
    :param dt: time interval between step.
    :param f: force function f(x,t), given x and t return a vector.
    :param t: current time

    :rtype: xn1 is the value of the vector at time t+dt; fxn is value f(xn,t), the value will be stored is for multi-step method

    """
    k1 = f(xn,t)
    xk1 = xn+dt/2*k1
    k2 = f(xk1,t+dt/2)
    xk2 = xn+dt/2*k2
    k3 = f(xk2,t+dt/2)
    xk3 = xn+dt*k3
    k4 = f(xk3,t+dt)
    xn1 = xn + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
    return xn1,k1

def ab2(xn, dt, f, t, cache):
    """.. function:: ab2(xn, dt, f, t, cache)

    The function do the time stepping of ODEs with second order Adams-Bashforth method.

    :param xn: value of the vector at t.
    :param dt: time interval between step.
    :param f: force function f(x,t), given x and t return a vector.
    :param t: current time
    :param cache: the value of f(x,t) at previous steps

    :rtype: xn1 is the value of the vector at time t+dt; fxn is value f(xn,t), the value will be stored is for multi-step method

    """
    fxn = f(xn,t)
    fxn_1 = cache
    xn1 = xn + dt/2 * (3*fxn - fxn_1)
    return xn1, fxn

def ab3(xn, dt, f, t, cache):
    """.. function:: ab3(xn, dt, f, t, cache)

    The function do the time stepping of ODEs with third order Adams-Bashforth method.

    :param xn: value of the vector at t.
    :param dt: time interval between step.
    :param f: force function f(x,t), given x and t return a vector.
    :param t: current time
    :param cache: the value of f(x,t) at previous steps

    :rtype: xn1 is the value of the vector at time t+dt;
    fxn is value f(xn,t), the value will be stored is for multi-step method

    """
    fxn = f(xn,t)
    fxn_1 = cache[0,:]
    fxn_2 = cache[1,:]
    xn1 = xn + dt/12 * (23*fxn - 16*fxn_1 + 5*fxn_2)
    return xn1, fxn

def ab4(xn, dt, f, t, cache):
    """.. function:: ab4(xn, dt, f, t, cache)

    The function do the time stepping of ODEs with fourth order Adams-Bashforth method.

    :param xn: value of the vector at t.
    :param dt: time interval between step.
    :param f: force function f(x,t), given x and t return a vector.
    :param t: current time
    :param cache: the value of f(x,t) at previous steps

    :rtype: xn1 is the value of the vector at time t+dt;
    fxn is value f(xn,t), the value will be stored is for multi-step method

    """
    fxn = f(xn,t)
    fxn_1 = cache[0,:]
    fxn_2 = cache[1,:]
    fxn_3 = cache[2,:]
    xn1 = xn + dt/24 * (55*fxn - 59*fxn_1 + 37*fxn_2 - 9*fxn_3)
    return xn1, fxn

def newmark(K, f, t, dt, un, nstep, cache, gamma):
    """.. function:: newmark(K, f, t, dt, un, nstep, cache, gamma)

    The function do the time stepping of ODEs with Newmark method, here we our vector have been normalized
    and the mass matrix is unit. Because we only deal with explicit scheme, we don't have the absorbing matrix.

    :param K: stiffness matrix
    :param f: force term f(t), given t return the force .
    :param t: current time
    :param dt: interval between time
    :param un: current displacement vector
    :param nstep: current number of steps
    :param cache: stored current velocity vector vn and acceleration vector an from previous step
    :param gamma: parameter \gamma of newmark scheme


    :rtype: un1 is the value of the displacement vector at time t+dt;
    cache stored velocity and acceleration vector un1 and an1

    """
    if nstep == 0:
        fn = f(t)
        v0 = cache
        a0 = fn - np.dot(K,un)
        cache = np.concatenate(([v0],[a0]))
        un1, cache = newmark(K, f, t, dt, un, 1, cache, gamma)
        return un1, cache
    else:
        fn = f(t+dt)
        vn = cache[0,:]
        an = cache[1,:]
        un1 = un + dt*vn + 0.5*dt*dt*an
        fn1 = fn - np.dot(K, un1)
        an1 = fn1
        vn1 = vn + (1-gamma)*dt*an + gamma*dt*an1
        cache = np.concatenate(([vn1],[an1]))
        return un1, cache















