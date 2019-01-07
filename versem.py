import numpy as np
import src.gll_library as gll
import src.mesh_spec as mesh
import src.el_mass_matrix as el_mass
import src.el_stiffness_matrix as el_stiff
import src.loc2glob as l2g
import src.model_parameters as mp
import src.force
import input.source_function
import matplotlib.pyplot as plt
plt.ion()
import src.tscheme
import sys



#In the final version, these should come from the input file.
#-----------------------------------------------------------------
ngll_x = 5 
ngll_y = 1
ngll_z = 5
nelm_x = 20
nelm_y = 1
nelm_z = 10
velocity_model = 'input/vel_mod.npy'
dim = 2

# Inputs
solver = 'newmark'
output_dir = 'results/timesteps/'

nt = 2000




# Obtain force make
force_term = np.array([1,1])
force_location = np.array([10*125,-5*125])


'''
# a is amplitude, b is location, c is width
a,b,c = 1,0.5,0.1
source_time_function = input.source_function.gaussian(t_source,a,b,c)
'''

# reading the mesh
ngll_el = ngll_x*ngll_y*ngll_z
el_no = nelm_x*nelm_y*nelm_z
X,Y,Z,connect = mesh.readEx('input/RectMesh.e')

#--------------------------------------------------------------------


# Obtaining the gll_coordinates and gll indices
gll_coordinates, gll_connect = mesh.mesh_interp2D(X,Y,Z,connect,ngll_x,ngll_z)
gll_coordinates = gll_coordinates*125
ngll_total = len(gll_coordinates)
sys.exit()
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
dt = eps*dxmin/2/np.amax(vp)

# Source time function
## creating source time function

# Ricker
t0 = 0.005
p0 = 0.1
#t_source = t
#source_time_function = input.source_function.ricker(t,t0,p0)
t_source,source_time_function = input.source_function.rickerSEM(dt,p0)

source_time_function = np.concatenate((source_time_function,np.zeros(nt-len(source_time_function))))


t = np.arange(0,nt,1)*dt

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


print('Computing Mass Matrix ...')
#The mass matrix => Computed only once if rho is constant in time.
#Mglob_mass = np.zeros([len(gll_coordinates),len(gll_coordinates)])
# M = el_mass.glob_el_mass_mat(gll_coordinates,gll_connect,rho,dN_local,W)

M_orig = el_mass.glob_el_mass_mat(gll_coordinates,gll_connect,rho,dN_local,W)

M_tmp = np.diag(M_orig)

M = np.zeros([2*ngll_total,2*ngll_total])
M[0:2*ngll_total-1:2,0:2*ngll_total-1:2] = np.diag(M_tmp[0:ngll_total])
M[1:2*ngll_total:2,1:2*ngll_total:2]   = np.diag(M_tmp[ngll_total:2*ngll_total])

print(np.shape(M))

print('Computing Stiffness Matrix ...')
#The stiffness matrix => Computed only once.
#We divide into three parts for ease of formulation and coding
comp = 0
Kglob_A0,Kglob_B0,Kglob_C0 = el_stiff.glob_el_stiff_mat(gll_coordinates,
                                gll_connect,dN_local,W,comp,dim,lmd,mu)

comp = 1
Kglob_A1,Kglob_B1,Kglob_C1 = el_stiff.glob_el_stiff_mat(gll_coordinates,
                                gll_connect,dN_local,W,comp,dim,lmd,mu)


K = np.zeros([2*ngll_total,2*ngll_total])

"""
# Top-left
K[0:ngll_total,0:ngll_total] = Kglob_A0[0,:,:] + Kglob_B0[0,:,:] + Kglob_C0
# Top-right
K[0:ngll_total,ngll_total:2*ngll_total] = Kglob_A0[1,:,:] + Kglob_B0[1,:,:]
# Bottom-left
K[ngll_total:2*ngll_total,0:ngll_total] = Kglob_A1[0,:,:] + Kglob_B1[0,:,:] 
#Bottom-right
K[ngll_total:2*ngll_total,ngll_total:2*ngll_total] = Kglob_A1[1,:,:] + Kglob_B1[1,:,:] + Kglob_C1

# Top-left 
K[0:2*ngll_total-1:2,0:2*ngll_total-1:2] = Kglob_A0[0,:,:] + Kglob_B0[0,:,:] + Kglob_C0
# Top-right
K[0:2*ngll_total-1:2,1:2*ngll_total:2] = Kglob_A0[1,:,:] + Kglob_B0[1,:,:]
# Bottom-left
K[1:2*ngll_total:2,0:2*ngll_total-1:2] = Kglob_A1[0,:,:] + Kglob_B1[0,:,:] 
#Bottom-right
K[1:2*ngll_total:2,1:2*ngll_total:2] = Kglob_A1[1,:,:] + Kglob_B1[1,:,:] + Kglob_C1
"""
# Top-left 
K[0:2*ngll_total:2,0:2*ngll_total:2] = Kglob_A0[0,:,:] + Kglob_B0[0,:,:] + Kglob_C0
# Top-right
K[0:2*ngll_total:2,1:2*ngll_total:2] = Kglob_A0[1,:,:] + Kglob_B0[1,:,:]
# Bottom-left
K[1:2*ngll_total:2,0:2*ngll_total:2] = Kglob_A1[0,:,:] + Kglob_B1[0,:,:] 
#Bottom-right
K[1:2*ngll_total:2,1:2*ngll_total:2] = Kglob_A1[1,:,:] + Kglob_B1[1,:,:] + Kglob_C1

print(np.shape(K))

print('Computing Force Vector ...')
# Force term
Fx_loc,Fy_loc = src.force.genforce(force_term,force_location,gll_coordinates)

# Computing the force vectors for x and y direction
Fx = src.force.glob_force_mat(gll_coordinates,gll_connect,Fx_loc,dN_local,W)
Fy = src.force.glob_force_mat(gll_coordinates,gll_connect,Fy_loc,dN_local,W)


F_tot = np.zeros(2*ngll_total)
F_tot[0:2*ngll_total-1:2] = Fx
F_tot[1:2*ngll_total:2]   = Fy


# Force vector interpolators
print(np.shape(np.concatenate((Fx,Fy))))
#F = src.force.F(t,source_time_function,np.concatenate((Fx,Fy)))
F = src.force.F(t,source_time_function,F_tot)
print(np.shape(F(0.2)))

# Create random c
C = np.zeros([len(F(0.2)),len(F(0.2))])

print('Time Stepping ...')
# Time Stepping with tscheme 
tstep = src.tscheme.Tscheme(solver=solver,M=M,K=-K,x0=u,t=t,f=F,C=C,outdir=output_dir)
tstep.process()

np.save("results/gll_coordinates.npy",gll_coordinates)

print('Done.')






