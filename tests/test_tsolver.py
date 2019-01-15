import unittest
import numpy as np
from .context import src

class TestTsolver(unittest.TestCase):

    """test all solvers in tsolver, compare the expected vaule from
    my calculation with the output of solver"""

    def test_euler_explicit(self):
        xn = np.array([1,3])
        K = np.array([[2,3],[4,1]])
        dt = 1
        t = 4
        f = lambda x,t : np.array([t+1,t*t]) - np.dot(K,x)
        xn1_out,_ = src.tsolver.euler_explicit(xn, dt, f, t)
        xn1_exp = np.array([-5,12])
        np.testing.assert_array_almost_equal(xn1_exp,xn1_out)

    def test_rk2(self):
        xn = np.array([1,3])
        K = np.array([[2,3],[4,1]])
        dt = 1
        t = 4
        f = lambda x,t : np.array([t+1,t*t]) - np.dot(K,x)
        xn1_out,_ = src.tsolver.rk2(xn, dt, f, t)
        xn1_exp = np.array([-12,24])
        np.testing.assert_array_almost_equal(xn1_exp,xn1_out)

    def test_rk4(self):
        xn = np.array([1,1])
        K = np.array([[1,0],[1,1]])
        dt = 1
        t = 4
        f = lambda x,t : np.array([4*t,2*t]) - 2*np.dot(K,x)
        xn1_out,_ = src.tsolver.rk4(xn, dt, f, t)
        xn1_exp = np.array([7,-19/3])
        np.testing.assert_array_almost_equal(xn1_exp,xn1_out)

    def test_ab2(self):
        xn = np.array([1,1])
        K = np.array([[1,0],[1,1]])
        cache = np.array([1, 1])
        dt = 1
        t = 4
        f = lambda x,t : np.array([4*t,2*t]) - 2*np.dot(K,x)
        xn1_out,_ = src.tsolver.ab2(xn, dt, f, t, cache)
        xn1_exp = np.array([21.5,6.5])
        np.testing.assert_array_almost_equal(xn1_exp,xn1_out)

    def test_ab3(self):
        xn = np.array([1,1])
        K = np.array([[1,0],[1,1]])
        cache = np.array([[1, 1],[0,0]])
        dt = 1
        t = 4
        f = lambda x,t : np.array([4*t,2*t]) - 2*np.dot(K,x)
        xn1_out,_ = src.tsolver.ab3(xn, dt, f, t, cache)
        xn1_exp = np.array([26.5,22/3])
        np.testing.assert_array_almost_equal(xn1_exp,xn1_out)

    def test_ab4(self):
        xn = np.array([1,1])
        K = np.array([[1,0],[1,1]])
        cache = np.array([[1, 0],[0,-1],[-1,1]])
        dt = 1
        t = 4
        f = lambda x,t : np.array([4*t,2*t]) - 2*np.dot(K,x)
        xn1_out,_ = src.tsolver.ab4(xn, dt, f, t, cache)
        xn1_exp = np.array([31,33/4])
        np.testing.assert_array_almost_equal(xn1_exp,xn1_out)

    def test_newmark(self):
        K = np.array([[1,0],[1,1]])
        f = lambda t: np.array([4*t,2*t])
        t = 1
        dt = 1
        un = np.array([1,1])
        cache = np.array([2,1])
        nstep = 0
        gamma = 0.5
        un1, cache = src.tsolver.newmark(K, f, t, dt, un, nstep, cache, gamma)

        un1_exp = np.array([4.5,2])
        vn1_exp = np.array([5.25,-0.25])
        an1_exp = np.array([3.5,-2.5])

        np.testing.assert_array_almost_equal(un1_exp,un1)
        np.testing.assert_array_almost_equal(vn1_exp,cache[0,:])
        np.testing.assert_array_almost_equal(an1_exp,cache[1,:])



if __name__ == "__main__":
    unittest.main()
