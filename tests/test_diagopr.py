import unittest
from .context import src
import numpy as np

class TestDiagopr(unittest.TestCase):

    """test all the functions in diagopr
    for reshape function, compare the matrix of output with expected ones
    for matrix operation functions, compare the result with numpy functions
    """

    def test_reshape(self):
        M = 1
        C = 2
        K = 3
        M_eff = np.array([[1,0],[0,1]])
        K_eff = np.array([[0,-1],[3,2]])
        M_out,K_out = src.diagopr.reshape(M,C,K,1)
        np.testing.assert_array_almost_equal(M_out,M_eff)
        np.testing.assert_array_almost_equal(K_out,K_eff)

        M = np.eye(2)
        K = np.array([[5,2.],[0.,1]])
        C = np.eye(2)*4
        M_eff = np.eye(4)
        K_eff = np.array([[0,0,-1,0],[0,0,0,-1],[5,2,4,0],[0,1,0,4]])
        M_out,K_out = src.diagopr.reshape(M,C,K,len(M))
        np.testing.assert_array_almost_equal(M_out,M_eff)
        np.testing.assert_array_almost_equal(K_out,K_eff)

    def test_reshape_f(self):
        f = 1
        f_eff = np.array([0,1])
        f_out = src.diagopr.reshape_f(f,1)
        np.testing.assert_array_almost_equal(f_out, f_eff)

        f = np.array([1,2,3])
        f_eff = np.array([0,0,0,1,2,3])
        f_out = src.diagopr.reshape_f(f,len(f))
        np.testing.assert_array_almost_equal(f_out, f_eff)

    def test_normalize(self):
        M = 2
        K = 5
        MK = 2.5
        MK_out = src.diagopr.normalize(M,K,1)
        self.assertEqual(MK,MK_out)

        M = np.array([[2,0],[0,3]])
        K = np.array([[3,5],[6,4]])
        K_out = src.diagopr.normalize(M,K,len(M))
        np.testing.assert_array_almost_equal(K_out,np.dot(np.linalg.inv(M),K))

    def test_normalize2(self):
        M = 2
        K = 5
        MK = 2.5
        MK_out = src.diagopr.normalize2(M, K, 1)
        self.assertEqual(MK, MK_out)

        M = np.array([[2, 0], [0, 3]])
        K = np.array([[3, 5], [6, 4]])
        K_out = src.diagopr.normalize2(M, K, len(M))
        np.testing.assert_array_almost_equal(K_out, np.dot(K,np.linalg.inv(M)))

    def test_normalize_f(self):
        M = 2
        f = 1
        f_out = src.diagopr.normalize_f(M,f,1)
        self.assertEqual(f_out,0.5)

        M = np.array([[2, 0], [0, 3]])
        f = np.array([2,6])
        f_out = src.diagopr.normalize_f(M,f,len(f))
        f_exp = np.array([1,2])
        np.testing.assert_array_almost_equal(f_out,f_exp)

    def test_diagmul_f(self):
        M = 2
        f = 1
        f_out = src.diagopr.diagmul_f(M,f,1)
        self.assertEqual(f_out,2)

        M = np.array([[2, 0], [0, 0.5]])
        f = np.array([2,6])
        f_out = src.diagopr.diagmul_f(M,f,len(f))
        f_exp = np.array([4,3])
        np.testing.assert_array_almost_equal(f_out,f_exp)



if __name__ == "__main__":
    unittest.main()