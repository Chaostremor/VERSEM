import numpy as np
from . import tsolver
import os

class Tscheme:
    def __init__(self, M, K, f, x0, t,outdir, C=None, gamma=0.5, solver='euler_explicit',solver2='rk4',interval=1):
        self.solver = solver
        self.solver2 = solver2
        self.M = M  ## no time dependence
        self.C = C  ## no time dependence
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

        self.K = tsolver.normalize2(self.M,self.K,self.dim)
        self.M0 = self.M
        self.x0[0,:] = tsolver.diagmul_M(self.M0,self.x0[0,:],self.dim)
        self.x0[1,:] = tsolver.diagmul_M(self.M0,self.x0[1,:],self.dim)

        self.M = np.eye(self.dim)


        ## self.cache (multiple steps) ; self.x(store result at certain time)

    def process(self):
        self._preprocess()
        if os.path.exists(self.outdir) == False:
            os.makedirs(self.outdir)
        for i in range(self.allstep):
            if self.order == 1:
                if self.nstep == 0:
                    pass
                un,vn = self.step()
                if self.nstep % self.interval == 0:
                    pass
                # print(un,vn)
            elif self.order == 2:
                if self.nstep == 0:
                    #pass
                    np.save(self.outdir+'u'+str(self.t[0]),self.u0)
                    #np.save(self.outdir + 'v' + str(self.t[0]), self.v0)
                    #np.save(self.outdir + 'a' + str(self.t[0]), self.a0)
                un,vn,an = self.step()
                if self.nstep % self.interval == 0:
                    #pass
                    np.save(self.outdir+'u'+str(self.t[self.nstep]),un)
                    #np.save(self.outdir + 'v' + str(self.t[self.nstep]), vn)
                    #np.save(self.outdir + 'a' + str(self.t[self.nstep]), an)
                if self.nstep % 50 == 49:
                    print('Timestep: %d of %d' % (i,self.allstep), 'Max. Displ.: %f ' % np.max(np.abs(un)))


    def _preprocess(self):
        self._getorder_()
        self._getsolverorder_()
        self._getinitial()


        if self.solverorder > self.order:
            raise Exception('Higher order scheme cannot solve lower order equations')
        if self.order == 2 and self.solverorder == 1:
            self.M,self.K = tsolver.reshape(self.M,self.C,self.K,len(self.M))
            self.C = None
            self.dim *=2
            self.x0 = np.reshape(self.x0,[self.dim])
            self.f0 = self.f
            self.f = self._reshape_f
        if self.solverorder == 1:
            self._MK = tsolver.normalize(self.M,self.K,len(self.K))
            self.f1 = self.f
            self.f = self._f

        self.xn = self.x0

    def _getinitial(self):
        if self.order == 1:
            self.u0 = self.x0
            f0 = self.f(self.t[0])
            self.v0 = tsolver.normalize_f(self.M,f0-np.dot(self.K,self.u0),len(f0))
        elif self.order == 2:
            self.u0 = self.x0[0,:]
            self.v0 = self.x0[1,:]
            f0 = self.f(self.t[0])
            self.a0 = tsolver.normalize_f(self.M,f0-np.dot(self.K,self.u0)-np.dot(self.C,self.v0),len(f0))
            print(self.nstep,self.u0,self.v0,self.a0)


    def step(self):
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
            #(M, C, K, f, t, dt, xn, nstep, cache, gamma)
            paralist = '(self.M,self.K,self.f,tn,dt,self.xn,self.nstep,self.cache,self.gamma)'
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
            vn = self.f(un,tn+dt)

            un = tsolver.normalize_f(self.M0, un, self.dim)
            vn = tsolver.normalize_f(self.M0, vn, self.dim)

            return un,vn

        elif self.order ==2:
            if self.solverorder == 1:
                un = result[0:int(len(result)/2)]
                vn = result[int(len(result)/2):]
                an_0 = self.f(result,tn+dt)
                an = an_0[int(len(result)/2):]

            elif self.solverorder == 2:
                un = result
                vn = cache[0,:]
                an = cache[1,:]

            un = tsolver.normalize_f(self.M0,un,len(un))
            vn = tsolver.normalize_f(self.M0, vn, len(vn))
            an = tsolver.normalize_f(self.M0, an, len(an))
            return un,vn,an

    def _getorder_(self):
        if self.C is None:
            self.order = 1
        else:
            self.order = 2

    def _getsolverorder_(self):
        if self.solver == 'euler_explicit' or self.solver == 'rk4' or self.solver == 'euler_adjust':
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

    def _reshape_f(self,t):
        fn = self.f0(t)
        fn = tsolver.reshape_f(fn,len(fn))
        return fn

    def _f(self,x,t):
        fn = self.f1(t)
        fxn = tsolver.normalize_f(self.M,fn,len(fn)) - np.dot(self._MK,x)
        return fxn








