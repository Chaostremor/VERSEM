import unittest
from .context import src
import numpy as np

class TestTscheme(unittest.TestCase):

    """ testing class Tscheme"""

    def test_init(self):
        """ testing the initialization function"""
        M = np.eye(2)*4
        K = np.array([[1,0.],[0.,1]])
        f = lambda t : np.array([0.3 * np.cos(t),-0.1 * np.cos(t)])
        x0 = np.array([[1, 1],[0,0]])
        t = np.array(list(range(100))) * 0.1
        outdir = "../result/timestep/"

        tstep = src.tscheme.Tscheme(M=M, K=K, x0=x0, t=t, f=f, outdir=outdir)
        np.testing.assert_array_almost_equal(tstep.M,M)
        np.testing.assert_array_almost_equal(tstep.K,K)
        np.testing.assert_array_almost_equal(tstep.f(1),f(1))
        np.testing.assert_array_almost_equal(tstep.x0,x0)
        np.testing.assert_array_almost_equal(tstep.t,t)
        self.assertEqual(tstep.outdir,outdir)

        self.assertEqual(tstep.dim,len(M))
        self.assertEqual(tstep.nstep,0)
        self.assertEqual(tstep.allstep,len(t)-1)

        self.assertEqual(tstep.C,None)
        self.assertEqual(tstep.gamma,0.5)
        self.assertEqual(tstep.solver,'euler_explicit')
        self.assertEqual(tstep.solver2,'rk4')
        self.assertEqual(tstep.interval,1)

        C = np.eye(2)*4
        gamma = 0.2
        solver = 'rk4'
        solver2 = 'euler_explicit'
        interval = 2
        tstep = src.tscheme.Tscheme(M=M, K=K, x0=x0, t=t, f=f, outdir=outdir,C=C,gamma=gamma,solver=solver,
                                    solver2=solver2,interval=interval)
        np.testing.assert_array_almost_equal(tstep.C,C)
        self.assertEqual(tstep.gamma,gamma)
        self.assertEqual(tstep.solver,solver)
        self.assertEqual(tstep.solver2,solver2)
        self.assertEqual(tstep.interval,interval)

    def test_getorder(self):
        """test _getorder_ function, if x0 is 2*n matrix, it should return 2 and
        if x0 is a vector it should return 1 """
        M = np.eye(2)*4
        K = np.array([[1,0.],[0.,1]])
        f = lambda t : np.array([0.3 * np.cos(t),-0.1 * np.cos(t)])
        x0 = np.array([[1, 1],[0,0]])
        t = np.array(list(range(100))) * 0.1
        outdir = "../result/timestep/"

        tstep = src.tscheme.Tscheme(M=M, K=K, x0=x0, t=t, f=f, outdir=outdir)
        tstep._getorder_()
        self.assertEqual(tstep.order,2)

        x0 = np.array([1,2])
        tstep = src.tscheme.Tscheme(M=M, K=K, x0=x0, t=t, f=f, outdir=outdir)
        tstep._getorder_()
        self.assertEqual(tstep.order, 1)

    def test_getsolverorder(self):
        """test _getsolverorder_ function"""
        M = np.eye(2)*4
        K = np.array([[1,0.],[0.,1]])
        f = lambda t : np.array([0.3 * np.cos(t),-0.1 * np.cos(t)])
        x0 = np.array([[1, 1],[0,0]])
        t = np.array(list(range(100))) * 0.1
        outdir = "../result/timestep/"

        solver = 'euler_explicit'
        tstep = src.tscheme.Tscheme(solver=solver,M=M, K=K, x0=x0, t=t, f=f, outdir=outdir)
        tstep._getsolverorder_()
        self.assertEqual(tstep.solverorder,1)
        self.assertEqual(tstep.steporder,1)

        solver = 'rk2'
        tstep = src.tscheme.Tscheme(solver=solver,M=M, K=K, x0=x0, t=t, f=f, outdir=outdir)
        tstep._getsolverorder_()
        self.assertEqual(tstep.solverorder,1)
        self.assertEqual(tstep.steporder,1)

        solver = 'rk4'
        tstep = src.tscheme.Tscheme(solver=solver,M=M, K=K, x0=x0, t=t, f=f, outdir=outdir)
        tstep._getsolverorder_()
        self.assertEqual(tstep.solverorder,1)
        self.assertEqual(tstep.steporder,1)

        solver = 'ab2'
        tstep = src.tscheme.Tscheme(solver=solver,M=M, K=K, x0=x0, t=t, f=f, outdir=outdir)
        tstep._getsolverorder_()
        self.assertEqual(tstep.solverorder,1)
        self.assertEqual(tstep.steporder,2)

        solver = 'ab3'
        tstep = src.tscheme.Tscheme(solver=solver,M=M, K=K, x0=x0, t=t, f=f, outdir=outdir)
        tstep._getsolverorder_()
        self.assertEqual(tstep.solverorder,1)
        self.assertEqual(tstep.steporder,3)

        solver = 'ab4'
        tstep = src.tscheme.Tscheme(solver=solver,M=M, K=K, x0=x0, t=t, f=f, outdir=outdir)
        tstep._getsolverorder_()
        self.assertEqual(tstep.solverorder,1)
        self.assertEqual(tstep.steporder,4)

        solver = 'newmark'
        tstep = src.tscheme.Tscheme(solver=solver,M=M, K=K, x0=x0, t=t, f=f, outdir=outdir)
        tstep._getsolverorder_()
        self.assertEqual(tstep.solverorder,2)

    def test_getinitial(self):
        """test _getinitial_ function"""
        M = np.eye(2)*4
        K = np.array([[1,0.],[0.,1]])
        f = lambda t : np.array([0.3 * np.cos(t),-0.1 * np.cos(t)])
        x0 = np.array([[1, 1],[0,0]])
        t = np.array(list(range(100))) * 0.1
        outdir = "../result/timestep/"

        tstep = src.tscheme.Tscheme(M=M, K=K, x0=x0, t=t, f=f, outdir=outdir)
        tstep._getorder_()
        tstep._getinitial()

        np.testing.assert_array_almost_equal(tstep.u0,x0[0,:])
        ##
        x0 = np.array([1,2])
        tstep = src.tscheme.Tscheme(M=M, K=K, x0=x0, t=t, f=f, outdir=outdir)
        tstep._getorder_()
        tstep._getinitial()

        np.testing.assert_array_almost_equal(tstep.u0, x0)

    def test_reshape_f(self):
        """test _reshape_f function"""
        M = np.eye(2)*4
        K = np.array([[1,0.],[0.,1]])
        f = lambda t : np.array([0.3 * t,-0.1 * t])
        x0 = np.array([[1, 1],[0,0]])
        t = np.array(list(range(100))) * 0.1
        outdir = "../result/timestep/"

        tstep = src.tscheme.Tscheme(M=M, K=K, x0=x0, t=t, f=f, outdir=outdir)
        tstep.f0 = f
        np.testing.assert_array_almost_equal(np.array([0,0,0.3,-0.1]),tstep._reshape_f(1))

    def test_f(self):
        """test _f"""
        M = np.eye(2)
        K = np.array([[1,2.],[-1.,1]])
        f = lambda t : np.array([0.3 * t,-0.1 * t])
        x0 = np.array([1, 1])
        t = np.array(list(range(100))) * 0.1
        outdir = "../result/timestep/"

        tstep = src.tscheme.Tscheme(M=M, K=K, x0=x0, t=t, f=f, outdir=outdir)
        tstep.f1 = f

        fn = f(0) - np.dot(K,x0)
        fn_out = tstep._f(x0,0)
        np.testing.assert_array_almost_equal(fn,fn_out)

    def test_preprocess(self):
        """test _preprocess function
        for 1st order ODEs, test the M, K x0 after normalization and f after moving
        for 2nd order ODEs with 1st order method, test M, K, x0, f after reshape and normalization and moving
        for 2nd order ODEs with 2st order method, test M, K, x0 after normalization
        """
        M = np.array([[1,0],[0,2]])
        K = np.array([[1,2.],[-1.,1]])
        f = lambda t : np.array([0.3 * t,-0.1 * t])
        x0 = np.array([1, 1])
        t = np.array(list(range(100))) * 0.1
        outdir = "../result/timestep/"

        tstep = src.tscheme.Tscheme(M=M, K=K, x0=x0, t=t, f=f, outdir=outdir)
        tstep._preprocess()
        x0_exp = np.array([1, 2])
        fn_out = tstep.f(x0_exp,0)
        fn_exp = f(0) - np.dot(K,x0)
        M_exp = np.eye(2)
        M0_exp = M
        K_exp = np.array([[1,1.],[-1.,0.5]])


        np.testing.assert_array_almost_equal(tstep.M,M_exp)
        np.testing.assert_array_almost_equal(tstep.M0, M0_exp)
        np.testing.assert_array_almost_equal(tstep.K, K_exp)
        np.testing.assert_array_almost_equal(tstep.x0, x0_exp)
        self.assertEqual(tstep.order,1)
        np.testing.assert_array_almost_equal(fn_out, fn_exp)
        self.assertEqual(tstep.C,None)

        M = np.array([[1,0],[0,2]])
        K = np.array([[1,2.],[-1.,1]])
        C = np.array([[4,6.],[-1.,4]])
        f = lambda t : np.array([0.3 * t,-0.1 * t])
        x0 = np.array([1, 1])
        t = np.array(list(range(100))) * 0.1
        outdir = "../result/timestep/"

        tstep = src.tscheme.Tscheme(M=M, C=C, K=K, x0=x0, t=t, f=f, outdir=outdir)
        tstep._preprocess()
        C_exp = np.array([[4,3.],[-1.,2]])
        np.testing.assert_array_almost_equal(tstep.C,C_exp)


        M = np.array([[1,0],[0,2]])
        K = np.array([[1,2.],[-1.,1]])
        f = lambda t : np.array([0.3 * t,-0.1 * t])
        x0 = np.array([[1, 1],[2,3]])
        t = np.array(list(range(100))) * 0.1
        outdir = "../result/timestep/"
        solver = 'newmark'

        tstep = src.tscheme.Tscheme(M=M, K=K, x0=x0, t=t, f=f, outdir=outdir,solver=solver)
        tstep._preprocess()
        x0_exp = np.array([[1, 2],[2,6]])
        M_exp = np.eye(2)
        M0_exp = M
        K_exp = np.array([[1,1.],[-1.,0.5]])

        np.testing.assert_array_almost_equal(tstep.M,M_exp)
        np.testing.assert_array_almost_equal(tstep.M0, M0_exp)
        np.testing.assert_array_almost_equal(tstep.K, K_exp)
        np.testing.assert_array_almost_equal(tstep.x0, x0_exp)
        self.assertEqual(tstep.order,2)

        M = np.array([[1,0],[0,2]])
        K = np.array([[1,2.],[-3.,8]])
        f = lambda t : np.array([0.3 * t,-0.1 * t])
        x0 = np.array([[1, 1],[2,3]])
        t = np.array(list(range(100))) * 0.1
        outdir = "../result/timestep/"
        solver = 'rk4'

        tstep = src.tscheme.Tscheme(M=M, K=K, x0=x0, t=t, f=f, outdir=outdir,solver=solver)
        tstep._preprocess()
        x0_exp = np.array([1, 2,2,6])
        M_exp = np.eye(4)
        M0_exp = M
        K_exp = np.array([[0,0,-1,0],[0,0,0,-1],[1,1,0,0],[-3,4,0,0]])
        fn_exp = np.concatenate(([0,0],f(0))) - np.dot(K_exp,x0_exp)
        fn_out = tstep.f(x0_exp,0)

        np.testing.assert_array_almost_equal(tstep.M,M_exp)
        np.testing.assert_array_almost_equal(tstep.M0, M0_exp)
        np.testing.assert_array_almost_equal(tstep.K, K_exp)
        np.testing.assert_array_almost_equal(tstep.x0, x0_exp)
        self.assertEqual(tstep.order,2)
        np.testing.assert_array_almost_equal(fn_exp, fn_out)

    def test_step(self):
        """test step function
        compare the result the calling step with calling directly from tsolver for different solvers,
        the result should all be the same
        """
        M = np.eye(2)
        K = np.array([[2,3],[4,1]])
        xn = np.array([1,3])
        dt = 0.1
        t = np.array(list(range(100))) * 0.1
        f1 = lambda x,t : np.array([t+1,t*t]) - np.dot(K,x)
        f = lambda t : np.array([t+1,t*t])
        outdir = "../result/timestep/"
        solver = 'euler_explicit'

        tstep = src.tscheme.Tscheme(M=M, K=K, x0=xn, t=t, f=f, outdir=outdir,solver=solver)
        tstep._preprocess()
        xn1_exp, _ = src.tsolver.euler_explicit(xn, dt, f1, t[0])
        xn1_out = tstep.step()

        np.testing.assert_array_almost_equal(xn1_exp, xn1_out)


        M = np.eye(2)
        K = np.array([[2,3],[4,1]])
        xn = np.array([1,3])
        dt = 0.1
        t = np.array(list(range(100))) * 0.1
        f1 = lambda x,t : np.array([t+1,t*t]) - np.dot(K,x)
        f = lambda t : np.array([t+1,t*t])
        outdir = "../result/timestep/"
        solver = 'rk4'

        tstep = src.tscheme.Tscheme(M=M, K=K, x0=xn, t=t, f=f, outdir=outdir,solver=solver)
        tstep._preprocess()
        xn1_out = tstep.step()
        xn1_exp, _ = src.tsolver.rk4(xn, dt, f1, t[0])

        np.testing.assert_array_almost_equal(xn1_exp, xn1_out)


        M = np.eye(2)
        K = np.array([[2,3],[4,1]])
        xn = np.array([1,3])
        dt = 0.1
        t = np.array(list(range(100))) * 0.1
        f1 = lambda x,t : np.array([t+1,t*t]) - np.dot(K,x)
        f = lambda t : np.array([t+1,t*t])
        outdir = "../result/timestep/"
        cache = np.array([1, 1])
        solver = 'ab2'

        tstep = src.tscheme.Tscheme(M=M, K=K, x0=xn, t=t, f=f, outdir=outdir,solver=solver)
        tstep._preprocess()
        tstep.cache = cache
        tstep.nstep = 5
        xn1_out = tstep.step()
        xn1_exp, _ = src.tsolver.ab2(xn, dt, f1, t[5],cache)

        np.testing.assert_array_almost_equal(xn1_exp, xn1_out)


        M = np.eye(2)
        K = np.array([[2,3],[4,1]])
        xn = np.array([1,3])
        dt = 0.1
        t = np.array(list(range(100))) * 0.1
        f1 = lambda x,t : np.array([t+1,t*t]) - np.dot(K,x)
        f = lambda t : np.array([t+1,t*t])
        outdir = "../result/timestep/"
        cache = np.array([[1, 1],[3,1]])
        solver = 'ab3'

        tstep = src.tscheme.Tscheme(M=M, K=K, x0=xn, t=t, f=f, outdir=outdir,solver=solver)
        tstep._preprocess()
        tstep.cache = cache
        tstep.nstep = 5
        xn1_out = tstep.step()
        xn1_exp, _ = src.tsolver.ab3(xn, dt, f1, t[5],cache)

        np.testing.assert_array_almost_equal(xn1_exp, xn1_out)


        M = np.eye(2)
        K = np.array([[2,3],[4,1]])
        xn = np.array([1,3])
        dt = 0.1
        t = np.array(list(range(100))) * 0.1
        f1 = lambda x,t : np.array([t+1,t*t]) - np.dot(K,x)
        f = lambda t : np.array([t+1,t*t])
        outdir = "../result/timestep/"
        cache = np.array([[1, 1],[3,1],[4,3]])
        solver = 'ab4'

        tstep = src.tscheme.Tscheme(M=M, K=K, x0=xn, t=t, f=f, outdir=outdir,solver=solver)
        tstep._preprocess()
        tstep.cache = cache
        tstep.nstep = 5
        xn1_out = tstep.step()
        xn1_exp, _ = src.tsolver.ab4(xn, dt, f1, t[5],cache)

        np.testing.assert_array_almost_equal(xn1_exp, xn1_out)


        M = np.eye(2)
        K = np.array([[2,3],[4,1]])
        xn = np.array([[1,3],[2,1]])
        dt = 0.1
        t = np.array(list(range(100))) * 0.1
        f = lambda t : np.array([t+1,t*t])
        outdir = "../result/timestep/"
        solver = 'newmark'
        un = xn[0,:]
        cache = xn[1,:]

        tstep = src.tscheme.Tscheme(M=M, K=K, x0=xn, t=t, f=f, outdir=outdir,solver=solver)
        tstep._preprocess()
        xn1_exp, _ = src.tsolver.newmark(K,f,t[0],dt,un,0,cache,0.5)
        xn1_out = tstep.step()

        np.testing.assert_array_almost_equal(xn1_exp, xn1_out)


        M = np.eye(2)
        K = np.array([[2,3],[4,1]])
        xn = np.array([[1,3],[2,1]])
        dt = 0.1
        t = np.array(list(range(100))) * 0.1
        f = lambda t : np.array([t+1,t*t])
        outdir = "../result/timestep/"
        solver = 'newmark'
        un = xn[0,:]
        cache = np.array([[2,1],[4,5]])

        tstep = src.tscheme.Tscheme(M=M, K=K, x0=xn, t=t, f=f, outdir=outdir,solver=solver)
        tstep._preprocess()
        tstep.nstep = 5
        tstep.xn = np.array([1,3])
        tstep.cache = cache
        xn1_exp, _ = src.tsolver.newmark(K,f,t[5],dt,un,5,cache,0.5)
        xn1_out = tstep.step()

        np.testing.assert_array_almost_equal(xn1_exp, xn1_out)

    def test_process(self):
        """test process
        solver the same equations with different schemes, the result should be very close for different schemes.
        """
        M = np.array([[2,0],[0,3]])
        K = np.array([[1,0.],[2.,4]])
        x0 = np.array([[1, 1],[0,0]])
        t = np.array(list(range(100))) * 0.01
        outdir = "../result/timestep/"
        def f(t):
            if t<0.1:
                return(np.array([0.3 * np.cos(t),-0.1 * np.cos(t)]))
            else:
                return (np.array([0,0]))

        solver = 'euler_explicit'
        tstep = src.tscheme.Tscheme(M=M, K=K, x0=x0, t=t, f=f, outdir=outdir, solver=solver)
        tstep.process()
        result1 = tstep.un


        M = np.array([[2,0],[0,3]])
        K = np.array([[1,0.],[2.,4]])
        x0 = np.array([[1, 1],[0,0]])
        t = np.array(list(range(100))) * 0.01
        outdir = "../result/timestep/"
        def f(t):
            if t<0.1:
                return(np.array([0.3 * np.cos(t),-0.1 * np.cos(t)]))
            else:
                return (np.array([0,0]))
        solver = 'rk2'
        tstep = src.tscheme.Tscheme(M=M, K=K, x0=x0, t=t, f=f, outdir=outdir, solver=solver)
        tstep.process()
        result2 = tstep.un

        np.testing.assert_array_almost_equal(result1,result2,decimal=2)

        M = np.array([[2,0],[0,3]])
        K = np.array([[1,0.],[2.,4]])
        x0 = np.array([[1, 1],[0,0]])
        t = np.array(list(range(100))) * 0.01
        outdir = "../result/timestep/"
        def f(t):
            if t<0.1:
                return(np.array([0.3 * np.cos(t),-0.1 * np.cos(t)]))
            else:
                return (np.array([0,0]))
        solver = 'rk4'
        tstep = src.tscheme.Tscheme(M=M, K=K, x0=x0, t=t, f=f, outdir=outdir, solver=solver)
        tstep.process()
        result3 = tstep.un

        np.testing.assert_array_almost_equal(result2,result3,decimal=3)

        M = np.array([[2,0],[0,3]])
        K = np.array([[1,0.],[2.,4]])
        x0 = np.array([[1, 1],[0,0]])
        t = np.array(list(range(100))) * 0.01
        outdir = "../result/timestep/"
        def f(t):
            if t<0.1:
                return(np.array([0.3 * np.cos(t),-0.1 * np.cos(t)]))
            else:
                return (np.array([0,0]))
        solver = 'ab2'
        tstep = src.tscheme.Tscheme(M=M, K=K, x0=x0, t=t, f=f, outdir=outdir, solver=solver)
        tstep.process()
        result4 = tstep.un

        np.testing.assert_array_almost_equal(result3,result4,decimal=3)

        M = np.array([[2,0],[0,3]])
        K = np.array([[1,0.],[2.,4]])
        x0 = np.array([[1, 1],[0,0]])
        t = np.array(list(range(100))) * 0.01
        outdir = "../result/timestep/"
        def f(t):
            if t<0.1:
                return(np.array([0.3 * np.cos(t),-0.1 * np.cos(t)]))
            else:
                return (np.array([0,0]))
        solver = 'ab3'
        tstep = src.tscheme.Tscheme(M=M, K=K, x0=x0, t=t, f=f, outdir=outdir, solver=solver)
        tstep.process()
        result5 = tstep.un

        np.testing.assert_array_almost_equal(result3,result5,decimal=3)

        M = np.array([[2,0],[0,3]])
        K = np.array([[1,0.],[2.,4]])
        x0 = np.array([[1, 1],[0,0]])
        t = np.array(list(range(100))) * 0.01
        outdir = "../result/timestep/"
        def f(t):
            if t<0.1:
                return(np.array([0.3 * np.cos(t),-0.1 * np.cos(t)]))
            else:
                return (np.array([0,0]))
        solver = 'ab4'
        tstep = src.tscheme.Tscheme(M=M, K=K, x0=x0, t=t, f=f, outdir=outdir, solver=solver)
        tstep.process()
        result6 = tstep.un

        np.testing.assert_array_almost_equal(result3,result6,decimal=3)


        M = np.array([[2,0],[0,3]])
        K = np.array([[1,0.],[2.,4]])
        x0 = np.array([[1, 1],[0,0]])
        t = np.array(list(range(100))) * 0.01
        outdir = "../result/timestep/"
        def f(t):
            if t<0.1:
                return(np.array([0.3 * np.cos(t),-0.1 * np.cos(t)]))
            else:
                return (np.array([0,0]))
        solver = 'newmark'
        tstep = src.tscheme.Tscheme(M=M, K=K, x0=x0, t=t, f=f, outdir=outdir, solver=solver)
        tstep.process()
        result7 = tstep.un

        np.testing.assert_array_almost_equal(result3,result7,decimal=3)

if __name__ == "__main__":
    unittest.main()