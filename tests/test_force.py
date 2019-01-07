"""Test for force term"""
from .context import src
import unittest

import numpy as np



class testForce(unittest.TestCase):
    """some doc string for the testclass
    """

    def test_Force(self):
        """Testing the forceterm function for point source at GLL point
        """
        
        # Make up GLL coordinates
        gll_coordinates = np.array([[0,0],[1,0],[0,1],[1,1],[2,0],[2,1]])
        
        # Force location
        force_location = np.array([1.6, 1.6])

        #Force values
        force_term = np.array([1,2])

        Fx,Fy = src.force.genforce(force_term,force_location,
                                                gll_coordinates)
        
        # Define Solution
        Fx_sol = np.array([0,0,0,0,0,1])
        Fy_sol = np.array([0,0,0,0,0,2])
        
        # Check if solution and computation match
        np.testing.assert_array_almost_equal(Fx,Fx_sol)
        np.testing.assert_array_almost_equal(Fy,Fy_sol)
        
    
    def testF(self):
        """Testing force vector interpolator"""
        
        # Create interpolating values
        tp = np.linspace(0,4,10000)
        fp = tp**2

        # Force vector
        Fx = np.array([0.,2.,0.])
        
        # Compute interpolator
        f = src.force.F(tp,fp,Fx)
        
        # Checking solution
        np.testing.assert_array_almost_equal(f(2),np.array([0,8,0]))


if __name__ == "__main__":
    unittest.main()
