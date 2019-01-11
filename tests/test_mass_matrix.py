from .context import src
import unittest
import numpy as np

class testGlobMassMat(unittest.TestCase):
    def test_el_mass_mat(self):
        a = np.array(range(10))
        b = np.array(range(10))
        c = np.array(range(10))
        el_mass_check = src.mass_matrix.element_mass_matrix(a,b,c)
        sol_el_mass_mat = np.diag(a**3)
        np.testing.assert_array_almost_equal(el_mass_check,sol_el_mass_mat)

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
        el_no = 2
        ngll_el = 4 

        dN_local = src.gll_library.dN_local_2D(xi,eta)
        W = src.gll_library.flattened_weights2D(w_xi,w_eta)
        rho = np.ones(len(gll_coordinates))

        Mglob_mass = src.mass_matrix.global_mass_matrix(gll_coordinates,gll_connect,rho,dN_local,W)


        Mg_sol = np.array([[0.25, 0.  , 0.  , 0.  , 0.  , 0.  ,0.  , 0.  , 0.  , 0.  , 0.  , 0   ],
                           [0.  , 0.25, 0.  , 0.  , 0.  , 0.  ,0.  , 0.  , 0.  , 0.  , 0.  , 0   ],
                           [0.  , 0.  , 0.5 , 0.  , 0.  , 0.  ,0.  , 0.  , 0.  , 0.  , 0.  , 0   ],
                           [0.  , 0.  , 0.  , 0.5 , 0.  , 0.  ,0.  , 0.  , 0.  , 0.  , 0.  , 0   ],
                           [0.  , 0.  , 0.  , 0.  , 0.25, 0.  ,0.  , 0.  , 0.  , 0.  , 0.  , 0   ],
                           [0.  , 0.  , 0.  , 0.  , 0.  , 0.25,0.  , 0.  , 0.  , 0.  , 0.  , 0   ],
                           [0.  , 0.  , 0.  , 0.  , 0.  , 0,   0.5 , 0.  , 0.  , 0.  , 0.  , 0   ],
                           [0.  , 0.  , 0.  , 0.  , 0.  , 0,   0.  , 0.5 , 0.  , 0.  , 0.  , 0.  ],
                           [0.  , 0.  , 0.  , 0.  , 0.  , 0,   0.  , 0.  , 0.25, 0.  , 0.  , 0.  ],
                           [0.  , 0.  , 0.  , 0.  , 0.  , 0,   0.  , 0.  , 0.  , 0.25, 0.  , 0.  ],
                           [0.  , 0.  , 0.  , 0.  , 0.  , 0,   0.  , 0.  , 0.  , 0.  , 0.25, 0.  ],
                           [0.  , 0.  , 0.  , 0.  , 0.  , 0,   0.  , 0.  , 0.  , 0.  , 0.  , 0.25]])



        np.testing.assert_array_almost_equal(Mglob_mass,Mg_sol)

if __name__ == "__main__":
    unittest.main()
