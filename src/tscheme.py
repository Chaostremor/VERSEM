"""This is a class that intake inputs of mass matrix, stiffness matrix, solver name and so on to do the time marching

Author: Fan Wu

"""

import numpy as np
from src import tsolver
from src import diagopr
import os

class Tscheme:
    """.. class:: Tscheme

    Class Tscheme solve the 1st order or 2nd order ODEs with different explicit solvers. This scheme only deals with
    constant coeffcient ODEs, the mass matrix has to be diagonal.

    """

    def __init__(self, M, K, f, x0, t,outdir, C=None, gamma=0.5, solver='euler_explicit',solver2='rk4',interval=1):
        """.. function:: __init__(self, M, K, f, x0, t,outdir, C=None, gamma=0.5,
        solver='euler_explicit',solver2='rk4',interval=1)

        Initialize the class

        :param M: mass matrix
        :param C: absorbing matrix, default None
        :param K: stiffness matrix
        :param f: force term f(t), given t return the force vector.
        :param t: time vector for stepping, t[0] should be the start time
        :param outdir: result output directory
        :param gamma: parameter for Newmark scheme, default 0.5
        :param solver: name of the main solver, default euler_explicit
        :param solver2: solver for the first few steps of the multistep method, default rk4
        :param interval: interval for output the result, should be an integer, default 1

        """
        self.solver = solver
        self.solver2 = solver2
        self.M = M  ## no time dependence
        self.C = C
        self.K = K  ## no time dependence
        self.f = f  ## time dependent, should be a function
        self.x0 = x0 ## first column is x, second column is x'
        self.t = t ## time for solving, an array
        self.dim = len(M)
        self.nstep = 0
        self.allstep = np.size(t) - 1
        self.gamma = gamma
        self.outdir = outdir
        self.interval = interval

        ## self.cache (multiple steps) ; self.x(store result at certain time

    def _getorder_(self):
        """.. function:: _getorder_(self)

        get the order of the ODEs

        """
        self.order = round(np.size(self.x0)/len(self.M))

    def _getsolverorder_(self):

        """.. function:: _getsolverorder_(self)

        get the information(multi-step or not; solving 1st order ODEs or second order) of the solvers

        """

        if self.solver == 'euler_explicit' or self.solver == 'rk4' or self.solver == 'rk2':
            self.solverorder = 1
            self.steporder = 1
        elif self.solver == 'ab2':
            self.solverorder = 1
            self.steporder = 2
        elif self.solver == 'ab3':
            self.solverorder = 1
            self.steporder = 3
        elif self.solver == 'ab4':
            self.solverorder = 1
            self.steporder = 4
        elif self.solver == 'newmark':
            self.solverorder = 2
        else:
            raise Exception('There is no solver named "%s."' %self.solver)

    def _getinitial(self):

        """.. function:: _getinitial(self)

        get initial value of the displacement vector

        """

        if self.order == 1:
            self.u0 = self.x0
        elif self.order == 2:
            self.u0 = self.x0[0,:]
        print("nsteps: %d, value: %e" % (self.nstep, np.max(np.abs(self.u0))))

    def _reshape_f(self,t):

        """.. function:: _reshape_f(self,t)

        reshape the force function for second order ODEs if use a first order scheme

        :param t: time t

        """

        ff = self.f0(t)
        fn = diagopr.reshape_f(ff, len(ff))
        return fn

    def _f(self,x,t):

        """.. function:: _f(self,t)

        turn equation x' + Mx == f(t) into x' = f(t) - Mx, _f is the right term with unknowns x and t

        :param x: vector x
        :param t: time t

        """

        fxn = self.f1(t) - np.dot(self.K,x)
        return fxn


    def _preprocess(self):

        """.. function:: _preprocess(self)

        preprocess the equations, it has the following procedure:

        1. get the order of equations and solvers
        2. normalize the vector, turn the mass matrix into unit
        3. if equation order is 2 and solver order is 1, the reorganize the matrix and force term to make it first
        order ODEs
        4. turn the equation x' + Mx == f(t) into x' = f(t) - Mx if it is first order ODEs

        """


        self._getorder_()
        self._getsolverorder_()
        self._getinitial()
        if self.solverorder > self.order:
            raise Exception('Higher order scheme cannot solve lower order equations')

        self.K = diagopr.normalize2(self.M, self.K, self.dim)
        if self.C is not None:
            self.C = diagopr.normalize2(self.M,self.C,self.dim)
        self.M0 = self.M
        self.M = np.eye(self.dim)


        if self.order == 1:
            self.x0 = diagopr.diagmul_f(self.M0, self.x0, self.dim)
        elif self.order == 2:
            self.x0[0, :] = diagopr.diagmul_f(self.M0, self.x0[0, :], self.dim)
            self.x0[1, :] = diagopr.diagmul_f(self.M0, self.x0[1, :], self.dim)

        if self.order == 2 and self.solverorder == 1:
            self.M,self.K = diagopr.reshape(self.M, self.C, self.K, len(self.M))
            self.C = None
            self.dim *=2
            self.x0 = np.reshape(self.x0,[self.dim])
            self.f0 = self.f
            self.f = self._reshape_f
        if self.solverorder == 1:
            self.f1 = self.f
            self.f = self._f

        self.xn = self.x0

    def step(self):

        """.. function:: step(self)

        do one time step using different methods, and return the vector at the new step

        :rtype un is the vector at the new step


        """


        ## initialize at the beginning of a step
        dt = self.t[self.nstep+1] - self.t[self.nstep]
        tn = self.t[self.nstep]

        if self.solverorder == 1:
            ## one step method
            if self.steporder == 1:
                xn = self.xn
                result,_ = eval('tsolver.' + self.solver + '(xn,dt,self.f,tn)')
                self.xn = result

            ## two step method
            elif self.steporder == 2:
                if self.nstep == 0:
                    xn = self.xn
                    result, fxn = eval('tsolver.' + self.solver2 + '(xn,dt,self.f,tn)')
                    self.xn = result
                    self.cache = fxn
                else:
                    xn = self.xn
                    result, fxn = eval('tsolver.' + self.solver + '(xn,dt,self.f,tn,self.cache)')
                    self.xn = result
                    self.cache = fxn

            ## multiple step method
            else:
                if self.nstep < self.steporder - 1:
                    xn = self.xn
                    result, fxn = eval('tsolver.' + self.solver2 + '(xn,dt,self.f,tn)')
                    self.xn = result
                    if self.nstep == 0:
                        self.cache = fxn
                    elif self.nstep == 1:
                        self.cache = np.concatenate(([fxn],[self.cache]))
                    else:
                        self.cache = np.concatenate(([fxn],self.cache))
                else:
                    xn = self.xn
                    result, fxn = eval('tsolver.' + self.solver + '(xn,dt,self.f,tn,self.cache)')
                    self.xn = result
                    self.cache = np.concatenate(([fxn],self.cache))
                    self.cache = np.delete(self.cache, -1, axis=0)

        ## newmark family
        elif self.solverorder == 2:
            if self.C is not None:
                raise Exception('There is no explicit newmark scheme for general absorbing matrix')
            #(M, C, K, f, t, dt, xn, nstep, cache, gamma)
            paralist = '(self.K,self.f,tn,dt,self.xn,self.nstep,self.cache,self.gamma)'
            if self.nstep == 0:
                self.cache = self.xn[1,:]
                self.xn = self.xn[0,:]
                result,cache = eval('tsolver.' + self.solver + paralist)
                self.xn = result
                self.cache = cache
            else:
                result, cache = eval('tsolver.' + self.solver + paralist)
                self.xn = result
                self.cache = cache

        self.nstep += 1

        ## reorganize the result
        if self.order == 1:
            un = result

        elif self.order ==2:
            if self.solverorder == 1:
                un = result[0:int(len(result)/2)]

            elif self.solverorder == 2:
                un = result

        return diagopr.normalize_f(self.M0, un, len(un))


    def process(self):

        """.. function:: process(self)

        solving the ODEs and save the data

        """

        self._preprocess()
        if os.path.exists(self.outdir) == False:
            os.makedirs(self.outdir)
        for i in range(self.allstep):
                if self.nstep == 0:
                    pass
                    # np.save(self.outdir+'u'+str(self.t[0]),self.u0)
                un = self.step()
                self.un = un
                if self.nstep % self.interval == 0:
                    pass
                    # np.save(self.outdir+'u'+str(self.t[self.nstep]),un)
                if self.nstep % 1 == 0:
                    print("nsteps: %d, value: %e" %(self.nstep,np.max(np.abs(un))))








