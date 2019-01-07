from .context import src
import unittest
import numpy as np

class testGlobMassMat(unittest.TestCase):
    def test_el_mass_glob(self):     
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
        comp = 0
        dim = 2

        dN_local = src.gll_library.dN_local_2D(xi,eta)

        W = src.gll_library.flattened_weights2D(w_xi,w_eta)

        gll_coords_el = gll_coordinates[gll_connect[0]]

        lmd = np.ones(len(gll_coords_el))
        mu = np.ones(len(gll_coords_el))

        #Mglob_A,Mglob_B,Mglob_C = src.el_stiffness_matrix.glob_el_stiff_mat(el_no,ngll_el,gll_coordinates,gll_connect,dN_local,W,comp,dim)
        A,B,C = src.el_stiffness_matrix.el_stiff(gll_coords_el,dim,ngll_el,dN_local,comp,W,lmd,mu)

        A_sol = np.array([[-0.5 ,  0.5,  0  ,  0   ],
                           [0.5 , -0.5,  0  ,  0   ],
                           [0   ,  0,  -0.5 ,  0.5 ],  
                           [0   ,  0,   0.5 , -0.5]])

        A_sol1 = np.array([[-0.25, -0.25, 0.25,   0.25],  
                            [0.25,  0.25, -0.25, -0.25], 
                           [-0.25, -0.25,  0.25,  0.25],  
                            [0.25,  0.25, -0.25, -0.25]])

            

        B_sol = np.array([[-0.5 ,  0.5,  0  ,  0   ],
                           [0.5 , -0.5,  0  ,  0   ],
                           [0   ,  0,  -0.5 ,  0.5 ],  
                           [0   ,  0,   0.5 , -0.5]])


        C_sol = np.array([[-1.,  0.5,  0.5,   0  ],
                          [0.5, -1.,  0,     0.5],  
                          [0.5,  0,  -1.,    0.5],  
                          [0,    0.5,  0.5, -1. ]])

        np.testing.assert_array_almost_equal(A[0,:,:],A_sol)
        np.testing.assert_array_almost_equal(B[0,:,:],B_sol)
        np.testing.assert_array_almost_equal(C,C_sol)


if __name__ == "__main__":
    unittest.main()
