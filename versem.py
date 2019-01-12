import numpy as np
import src.gll_library as gll
import src.mesh_spec as mesh
import src.mass_matrix as Mass
import src.stiffness_matrix as Stiffness
import src.loc2glob as l2g
import src.model_parameters as mp
import src.force
import input.source_function
import matplotlib.pyplot as plt
plt.ion()
import src.tscheme
import sys
import time
from src.time_decorator import timer


if __name__ == '__main__':
    pass

    #In the final version, these should come from the input file.
    #-----------------------------------------------------------------
    ngll_x = 3
    ngll_y = 1
    ngll_z = 3
    velocity_model = 'input/vel_mod_new.npy'
    dim = 2

    # Inputs
    solver = 'newmark'
    output_dir = 'results/timesteps/'

    nt = 1000


    # Obtain force make
    force_term = np.array([1,0])
    force_location = np.array([10*125,-2*125])


    '''
    # a is amplitude, b is location, c is width
    a,b,c = 1,0.5,0.1
    source_time_function = input.source_function.gaussian(t_source,a,b,c)
    '''

    # reading the mesh
    ngll_el = ngll_x*ngll_y*ngll_z
    X,Y,Z,connect = mesh.readEx('input/RectMesh.e')

    #--------------------------------------------------------------------


    # Start timer 
    start_time = time.time()
    print("Start Meshing...")

    # Obtaining the gll_coordinates and gll indices
    gll_coordinates, gll_connect = mesh.mesh_interp2D(X,Y,Z,connect,ngll_x,ngll_z)
    gll_coordinates = gll_coordinates*125
    ngll_total = len(gll_coordinates)


    # takes in global flattened density matrix. 1D matrix of length = total
    # number of gll points in the mesh.
    [rho_real,vp_real,vs_real] = mesh.assignSeismicProperties(velocity_model,gll_coordinates)

    rho_unit = 1.
    v_unit = 1.

    rho = rho_real/rho_unit
    vp = vp_real/v_unit
    vs = vs_real/v_unit

    gll_coordinates

    dxmin = gll_coordinates[1,0] - gll_coordinates[0,0]
    # Stability condition
    # Stability condition  
    eps = 0.25           # Courant value
    dt = eps*dxmin/2/np.amin(vp_real)

    # Source time function
    ## creating source time function

    # Ricker
    t0 = 0.005
    p0 = 0.3/ngll_el*9
    #t_source = t
    #source_time_function = input.source_function.ricker(t,t0,p0)
    t_source,source_time_function = input.source_function.rickerSEM(dt,p0)

    source_time_function = np.concatenate((source_time_function,np.zeros(nt-len(source_time_function))))



    t = np.arange(0,nt,1)*dt

    plt.plot(t,source_time_function)
    plt.show()
    #sys.exit()

    # Create zero displacement vector 
    u = np.zeros([2,2*len(gll_coordinates)])


    # Conversion from rho vp vs to lambda and mu
    mu,lmd = mp.velocity_conversion(rho,vp,vs)

    xi,w_xi = gll.gll_pw(ngll_x-1)
    #zi,w_zi = gll.gll_pw(ngll_y-1)
    eta,w_eta = gll.gll_pw(ngll_z-1)

    # To store the weights for all the gll points in an element.
    W = gll.flattened_weights2D(w_xi,w_eta)

    # derivative of the shape functions with respect to the local coords.
    # array of [no of gll points in an element]x[dim]x[no of gll points in an element]
    dN_local = gll.dN_local_2D(xi,eta)



    # Print Time Needed for Meshing 
    print("Finished Meshing.")
    print("Time: %s sec" % (time.time() - start_time))



    # ---------------- MASS MATRIX ------------------------------------
    #The mass matrix => Computed only once if rho is constant in time.
    #Mglob_mass = np.zeros([len(gll_coordinates),len(gll_coordinates)])
    # M = el_mass.glob_el_mass_mat(gll_coordinates,gll_connect,rho,dN_local,W)

    M = timer(Mass.global_mass_matrix,start_time,"computing mass matrix",gll_coordinates,gll_connect,rho,dN_local,W)
    print(np.shape(M))



    # ---------------- STIFFNESS MATRIX ------------------------------------
    #The stiffness matrix => Computed only once.
    #We divide into three parts for ease of formulation and coding

    comp = 0
    Kglob_A0,Kglob_B0,Kglob_C0 = timer(Stiffness.global_stiffness_matrix,start_time,
                                        "computing 1st component of stiffness matrix",
                                        gll_coordinates, gll_connect,dN_local,W,comp,dim,lmd,mu)

    comp = 1
    Kglob_A1,Kglob_B1,Kglob_C1 = timer(Stiffness.global_stiffness_matrix,start_time,
                                        "computing 2nd component of stiffness matrix",
                                        gll_coordinates, gll_connect,dN_local,W,comp,dim,lmd,mu)


    ## Reassign locations for the Stiffness components
    # Initialize empty Stiffness matrix
    K = np.zeros([2*ngll_total,2*ngll_total])
    # Redistributing Top-Left
    K[0:2*ngll_total:2,0:2*ngll_total:2] = Kglob_A0[0,:,:] + Kglob_B0[0,:,:] + Kglob_C0
    # Redistributing Top-right
    K[0:2*ngll_total:2,1:2*ngll_total:2] = Kglob_A0[1,:,:] + Kglob_B0[1,:,:]
    # Redistributing Bottom-left
    K[1:2*ngll_total:2,0:2*ngll_total:2] = Kglob_A1[0,:,:] + Kglob_B1[0,:,:] 
    # Redistributing Bottom-right
    K[1:2*ngll_total:2,1:2*ngll_total:2] = Kglob_A1[1,:,:] + Kglob_B1[1,:,:] + Kglob_C1

    print(np.shape(K))

    #sys.exit()


    # ---------------- FORCE VECTOR ------------------------------------
    # Force term
    F = timer(src.force.F,start_time,"computing the force vector",force_term,
              force_location,t,source_time_function,gll_coordinates,
              gll_connect,dN_local,W)


    print(np.shape(F(0.2)))

    # Create random c
    C = np.zeros([len(F(0.2)),len(F(0.2))])




    # ---------------- Time Stepping ------------------------------------

    # Time Stepping with tscheme 
    tstep = src.tscheme.Tscheme(solver=solver,M=M,K=-K,x0=u,t=t,f=F,C=C,outdir=output_dir)
    timer(tstep.process,start_time,"Time Marching")

    np.save("results/gll_coordinates.npy",gll_coordinates)


    print('Done.')






