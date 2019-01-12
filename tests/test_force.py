"""Test for force term"""
from .context import src
import unittest

import numpy as np



class testForce(unittest.TestCase):
    """some doc string for the testclass
    """

    def test_element_force_matrix(self):
        """ Testing the element force matrix computation performed by
        element_force_matrix() in src/force.py"""

        # Force term
        f = np.array(range(10))
        # Jacobian
        J = np.array(range(10))
        # Weights
        W = np.array(range(10))
        
        # Computation
        el_force_check = src.force.element_force_matrix(f,J,W)
        
        # Solution
        sol_el_force_mat = np.diag(f**3)

        # Check
        np.testing.assert_array_almost_equal(el_force_check,sol_el_force_mat)

    def test_global_mass_matrix(self):     
        ###### 1 ######
        # Setting the order
        N = 1

        # Getting collocation points
        xi,w_xi  = src.gll_library.gll_pw(N)
        eta,w_eta = src.gll_library.gll_pw(N)
        
        # Creating arbitrary coordinate matrix as shown in description
        gll_coordinates = np.array([[0,0],[1,0],[0,1],[1,1],[2,0],[2,1]])
        gll_connect = np.array([[1,2,3,4],[2,5,4,6]]) - 1

        # Local derivative Matrix
        dN_local = src.gll_library.dN_local_2D(xi,eta)
        # 1D Weight array
        W = src.gll_library.flattened_weights2D(w_xi,w_eta)

        # Ones force array
        f = np.ones(len(gll_coordinates))

        # Global computed force vector
        Fg = src.force.global_force_vector(gll_coordinates,gll_connect,f,dN_local,W)

        # Global force vector solution
        Fg_sol = np.array([0.25, 0.5 , 0.25, 0.5, 0.25, 0.25])
                           
        # TEST
        np.testing.assert_array_almost_equal(Fg,Fg_sol)

    def testGenForce(self):
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
        

        ## GLL SETUP
        # Polynomial order
        N = 1
        # Make up GLL coordinates and connectivity matrix
        gll_coordinates = np.array([[0,0],[1,0],[0,1],[1,1],[2,0],[2,1]])
        gll_connect = np.array([[1,2,3,4],[2,5,4,6]]) - 1

        # Getting collocation points
        xi,w_xi  = src.gll_library.gll_pw(N)
        eta,w_eta = src.gll_library.gll_pw(N)

        # Local derivative Matrix
        dN_local = src.gll_library.dN_local_2D(xi,eta)
        # 1D Weight array
        W = src.gll_library.flattened_weights2D(w_xi,w_eta)

        ## FORCE
        # Force location
        force_location = np.array([1.6, 1.6])

        # Force values
        force_term = np.array([1,2])

        # Create interpolating values
        t = np.linspace(0,4,1000)
        source_time_function = t**2
        
        # Compute interpolator
        f = src.force.F(force_term,force_location,t,source_time_function,gll_coordinates,gll_connect,dN_local,W)
        
        # Solution for f
        f_sol= np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 8])

        # Checking solution
        np.testing.assert_array_almost_equal(f(4),f_sol)


if __name__ == "__main__":
    unittest.main()
