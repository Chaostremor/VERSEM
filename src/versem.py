import numpy as np
import gll_library as gll
import mesh_spec as mesh
import el_mass_matrix as el_mass
import el_stiffness_matrix as el_stiff
import loc2glob as l2g

#In the final version, these should come from the input file.
#-----------------------------------------------------------------
ngll_x = 5 
ngll_y = 1
ngll_z = 5
nelm_x = 20
nelm_y = 1
nelm_z = 10
velocity_model = '../input/vel_mod.npy'
dim = 2


ngll_el = ngll_x*ngll_y*ngll_z
el_no = nelm_x*nelm_y*nelm_z
X,Y,Z,connect = mesh.readEx('../input/RectMesh.e')

#--------------------------------------------------------------------
# Obtaining the gll_coordinates and gll indices
gll_coordinates, gll_connect = mesh.mesh_interp2D(X,Y,Z,connect,ngll_x,ngll_z)

# takes in global flattened density matrix. 1D matrix of length = total
# number of gll points in the mesh.
[rho,vp,vs] = mesh.assignSeismicProperties(velocity_model,gll_coordinates)

# Conversion from rho vp vs to lambda and mu
lmd,mu = velocity_conversion(rho,vp,vs)

xi,w_xi = gll.gll_pw(ngll_x-1)
#zi,w_zi = gll.gll_pw(ngll_y-1)
eta,w_eta = gll.gll_pw(ngll_z-1)

# To store the weights for all the gll points in an element.
W = gll.flattened_weights2D(w_xi,w_eta)

# derivative of the shape functions with respect to the local coords.
# array of [no of gll points in an element]x[dim]x[no of gll points in an element]
dN_local = gll.dN_local_2D(xi,eta)

#The mass matrix => Computed only once if rho is constant in time.
#Mglob_mass = np.zeros([len(gll_coordinates),len(gll_coordinates)])
Mglob_mass = el_mass.glob_el_mass_mat(gll_coordinates,gll_connect,rho,dN_local,W)

#The stiffness matrix => Computed only once.
#We divide into three parts for ease of formulation and coding
comp = 0
Mglob_A,Mglob_B,Mglob_C = el_stiff.glob_el_stiff_mat(gll_coordinates,gll_connect,dN_local,W,comp,dim,lmd,mu)








