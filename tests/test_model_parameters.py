""" Testing the Model parameter conversion
"""
from .context import src
import numpy as np
import unittest

class testing_Lame(unittest.TestCase):
    """Testing Lame"""
    def testLame(self):
        """Testing velocity conversion in src/model_parameters.py
        """

        # Create artificial vp,vs and rho
        vp = np.array([3,4])
        vs = np.array([2,3])
        rho = np.array([1,2])

        # Create solution 
        mu_sol = np.array([4,18])
        lmd_sol = np.array([1,-4])
        
        # calculate lambda and mu from vp,vs,rho
        mu,lmd = src.model_parameters.velocity_conversion(rho,vp,vs)
        
        # Checking the solutions
        np.testing.assert_array_equal(mu,mu_sol)
        np.testing.assert_array_equal(lmd,lmd_sol)




if __name__ == "__main__":
    unittest.main()
