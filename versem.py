import numpy as np

# Source files
import src
import src.gll_library as gll
import src.mass_matrix as Mass
import src.stiffness_matrix as Stiffness
import src.force as Force

# Input Source Time Function
import input

# Timer
import time
from src.time_decorator import timer


if __name__ == '__main__':
    pass

    # In the final version, these should come from the input file.
    #-------------------- INPUTS ------------------------
    ngll_x = 3
    ngll_y = 1
    ngll_z = 3
    velocity_model = 'input/2Layer1000.npy'
    dim = 2

    # Inputs
    solver = 'newmark'
    output_dir = 'results/timesteps/'

    nt = 5000

    # Obtain force make
    force_term = np.array([1,0])
    force_location = np.array([9*1000,-4*1000])


    # reading the mesh
    ngll_el = ngll_x*ngll_y*ngll_z
    X,Y,Z,connect = src.mesh_spec.readEx('input/RectMesh.e')

    #--------------------------------------------------------------------


    # Start timer 
    start_time = time.time()
    print("Start Meshing...")

    ########### Should be part of the Mesh Object as well ?????? I think ########################
    # Obtaining the gll_coordinates and gll indices                                             #
    gll_coordinates, gll_connect = src.mesh_spec.mesh_interp2D(X,Y,Z,connect,ngll_x,ngll_z)     #
                                                                                                #
    # Note, how I multiplied by a thousand below. The mesh itself is only [0,20]x[0,-10]        #
    # Therefore, to make the values more realistic, I decided to multiply by four.              #
    # I did not have time to create a new mesh sadly. Bad, but quick fix                        #
    gll_coordinates = gll_coordinates*1000                                                      #
    ngll_total = len(gll_coordinates)                                                           #
                                                                                                #
                                                                                                #
    # takes in global flattened density matrix. 1D matrix of length = total                     #
    # number of gll points in the mesh.                                                         #
    [rho,vp,vs] = src.mesh_spec.assignSeismicProperties(velocity_model,gll_coordinates)         #

    # Conversion from rho vp vs to lambda and mu
    mu,lmd = src.model_parameters.velocity_conversion(rho,vp,vs)

    # ---------------- GLL Setup ---------------------------
    xi,w_xi = gll.gll_pw(ngll_x-1)
    #zi,w_zi = gll.gll_pw(ngll_y-1)
    eta,w_eta = gll.gll_pw(ngll_z-1)

    # To store the weights for all the gll points in an element.
    W = gll.flattened_weights2D(w_xi,w_eta)

    # derivative of the shape functions with respect to the local coords.                       #
    # array of [no of gll points in an element]x[dim]x[no of gll points in an element]          #
    dN_local = gll.dN_local_2D(xi,eta)                                                          #
                                                                                                #
    #############################################################################################


    

    # ---------------- Stability condition ---------------------------

    # Minimum dx value (only works for rectangular mesh)
    dxmin = gll_coordinates[1,0] - gll_coordinates[0,0]
    
    # Courant value
    eps = 0.025
    dt = eps*dxmin/2/np.amax(vp)

    # Time vector for the time stepping.
    t = np.arange(0,nt,1)*dt

    # ----------------------------------------------------------------



    # ---------------- Source time function --------------------------
    # Peak frequency
    f0 = 5
    source_time_function = input.source_time_function.gaussian(t,f0)
    # ----------------------------------------------------------------

    

    # Print Time Needed for Meshing 
    print("Finished Meshing.")
    print("Time: %s sec" % (time.time() - start_time))




    # ---------------- MASS MATRIX ------------------------------------
    # Computed only once before the time stepping

    M = timer(Mass.global_mass_matrix,start_time,"computing mass matrix",gll_coordinates,gll_connect,rho,dN_local,W)
    print(np.shape(M))

    # -----------------------------------------------------------------


    # ---------------- STIFFNESS MATRIX -------------------------------
    # Computed only once before the time stepping
    K = timer(Stiffness.global_stiffness_matrix,start_time,
              "computing the stiffness matrix",
              gll_coordinates, gll_connect,dN_local,W,dim,lmd,mu)
    print(np.shape(K))
    
    # ----------------------------------------------------------------


    # ---------------- FORCE VECTOR ------------------------------------
    # Computed only once before the time stepping
    F = timer(Force.F,start_time,"computing the force vector",force_term,
              force_location,t,source_time_function,gll_coordinates,
              gll_connect,dN_local,W)
    print(np.shape(F(0.2)))
    # ----------------------------------------------------------------




    # ---------------- Time Stepping ------------------------------------

    # Create zero displacement vector (Last two time steps)
    u = np.zeros([2,2*len(gll_coordinates)])

    # Create random C
    C = np.zeros([len(F(0.2)),len(F(0.2))])

    # Time Stepping with tscheme 
    tstep = src.tscheme.Tscheme(solver=solver,M=M,K=-K,x0=u,t=t,f=F,C=C,outdir=output_dir)
    timer(tstep.process,start_time,"Time Marching")


    # Save meshdata.
    np.save("results/gll_coordinates.npy",gll_coordinates)
    np.save("results/force_term.npy",force_term)
    np.save("results/force_location.npy",force_location)


    print('Done.')