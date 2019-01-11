"""This is script contains functions used to get from a quadrilateral
mesh to a Spectral-Element mesh. File for mesh object creation?"""

# Necessary for calculation and import of exodus file
import numpy as np
import netCDF4

# Import GLL library to get the lagrange polynomials for interpolation
# of the grid
from . import gll_library as gll



def readUniformVelocity(input_mesh,outfilename):
    """.. function testAssignPorperties()
    
    Saves ``.npy`` array file with  
    ::
    
        Columns:    X Y Z rho vp vs
    
    """
    
    # Get XYZ coordinates of the mesh
    X,Y,Z,connect = readEx(input_mesh)

    # Two layered material for test
    rho1 = 2000  # kg/m^3
    vs1  = 1600  # m/s.
    vp1  = vs1*1.7  # m/s 

    rho2 = 2300  # kg/m^3
    vs2  = 2200  # m/s
    vp2  = vs2*1.7  # m/s 

    # Initialize empty array
    prop = np.zeros([len(X),6])

    # Populate property matrix
    prop[:,0] = X*125
    prop[:,1] = Y*125
    prop[:,2] = Z*125

    index1 = np.where(prop[:,2]>=-5*125)
    print(index1)
    index2 = np.where(prop[:,2]<-5*125)
    print(index2)

    prop[index1,3] = rho1
    prop[index1,4] = vp1
    prop[index1,5] = vs1

    prop[index2,3] = rho2
    prop[index2,4] = vp2
    prop[index2,5] = vs2

    # Save materials to .npy file
    np.save(outfilename,prop)


def assProp():
    """testning the property assignment"""
    
    ngllx = 5
    ngllz = 5
    X,Y,Z,connect = readEx('input/RectMesh.e')
    gll_coordinates, gll_connect = \
                            mesh_interp2D(X,Y,Z,connect,ngllx,ngllz)
    
    v_mod = 'input/vel_mod.npy'

    # assign properties
    prop = assignSeismicProperties(v_mod,gll_coordinates)


def assignSeismicProperties(velocity_model,gll_coordinates):
    """Reads an ``.npy`` file that contains a ``numpy`` array with velocity 
    model of the format
    
    ::
    
        Columns: X, Y, Z, rho, vp ,vs


    Then, it assigns the properties to the interpolated mesh using the 
    nearest neighbour principle.

    :param v_name: velocity model array specified as above.

    :param gll_coordinates: coordinate matrix of the interpolated GLL
                            points using the original matrix.
        
    :rtype: Nx3 ``numpy`` array with the seismic properties for each 
            node in gll_coordinates

    ::
        
        Columns: rho,vp,vs
    

    I think this should be correct.

    """
    
    # Load .npy file
    v = np.load(velocity_model)

    # Number of nodes
    N = len(gll_coordinates[:,0])

    # Initialize propertz matrix
    prop = np.zeros([N,3])

    # Check which coordinate is the closest to interpolated node in for 
    # loop

    for i in range(N):
        # Here for now only 2D
        dist = np.sqrt( (v[:,0]-gll_coordinates[i,0])**2  
                        + (v[:,2]-gll_coordinates[i,1])**2)
        index = np.argmin(dist)
        prop[i,:] = v[index,3::]

    return prop




def readUniformVelocity(input_mesh,outfilename):
    """.. function testAssignPorperties()
    
    Saves ``.npy`` array file with  
    ::
    
        Columns:    X Y Z rho vp vs
    
    """
    
    # Get XYZ coordinates of the mesh
    X,Y,Z,connect = readEx(input_mesh)

    # Uniform material for test
    rho = 2700  # kg/m^3
    vp  = 6500  # m/s 
    vs  = 3750  # m/s

    # Initialize empty array
    prop = np.zeros([len(X),6])

    # Populate property matrix
    prop[:,0] = X
    prop[:,1] = Y
    prop[:,2] = Z
    prop[:,3] = rho
    prop[:,4] = vp
    prop[:,5] = vs
    
    # Save materials to .npy file
    np.save(outfilename,prop)


def assProp():
    """testning the property assignment"""
    
    ngllx = 5
    ngllz = 5
    X,Y,Z,connect = readEx('../input/RectMesh.e')
    gll_coordinates, gll_connect = \
                            mesh_interp2D(X,Y,Z,connect,ngllx,ngllz)
    
    v_mod = '../input/vel_mod.npy'

    # assign properties
    prop = assignSeismicProperties(v_mod,gll_coordinates)


def assignSeismicProperties(velocity_model,gll_coordinates):
    """.. function:: assignSeismicProperties(name)
    
    Reads an ``.npy`` file that contains a ``numpy`` array with velocity 
    model of the format
    
    ::
    
        Columns:    X Y Z rho vp vs

    Then, it assigns the properties to the interpolated mesh using the 
    nearest neighbour principle.

    :param v_name: velocity model array specified as above.
    :param gll_coordinates: coordinate matrix of the interpolated GLL
                            points using the original matrix.
        
    :rtype: 3xN ``numpy`` array with the seismic properties for each 
            node in gll_coordinates
    ::
        
        Columns     rho vp vs
    
    Return type 3XN for easier assignment to [rho,vp,vs]


    """
    
    # Load .npy file
    v = np.load(velocity_model)

    # Number of nodes
    N = len(gll_coordinates[:,0])

    # Initialize propertz matrix
    prop = np.zeros([N,3])

    # Check which coordinate is the closest to interpolated node in for 
    # loop

    for i in range(N):
        # Here for now only 2D
        dist = np.sqrt( (v[:,0]-gll_coordinates[i,0])**2  
                        + (v[:,2]-gll_coordinates[i,1])**2)
        index = np.argmin(dist)
        prop[i,:] = v[index,3::]

    return np.transpose(prop)



def readEx(name):
    """readEx(name)
    Function reads Exodus file and outputs Coordinates, connection 
    matrix, boundary coordinates and nodes.
    """
    
    # Loading the Dataset
    nc = netCDF4.Dataset(name) 
    
    # Getting coordinates
    X = nc.variables['coordx']
    Y = nc.variables['coordy']
    Z = nc.variables['coordz']

    # Getting global coordinate numbering                                   
    connect = np.array(nc.variables['connect1'])-1    
    
    return (np.array(X[:]),np.array(Y[:]),np.array(Z[:]),np.array(connect[:]))



def mesh_interp2D(X,Y,Z,connect,ngllx,nglly):
    """.. function :: mesh_interp(X,Y,Z,connect,ngllx,nglly)
    
    Creates new gll point mesh from initial control point mesh for 
    variable accuracy vs. efficiency. 
    
    :param X: x coordinates in Nx1 ``numpy`` array
    :param Y: y coordinates in Nx1 ``numpy`` array
    :param Z: z coordinates in Nx1 ``numpy`` array
    
    :param connect: 1x[Number of control points] ``numpy`` array that
                    defines the connectivity of each element
    
    :param ngllx: number of new gll points in xi direction
    :param nglly: number of new gll points in eta direction

    Function takes in coordinates of meshgrid and its connectivity 
    matrix. Then, it interpolates the GLL points onto the global grid 
    and defines the numbering for the new found set of points.
    """

    # Number of elements (nel) from the number of rows of connectivity 
    # matrix and number of control points (Nn) from the columns
    nel,Nn = connect.shape
    
    # Number of points per side
    Nside = int(np.round(np.sqrt(Nn)))

    # Local control points
    xi_control ,__ = gll.gll_pw(Nside-1) # -1 since gll_pw takes in polynomial degree


    # Node setup in mesh file
    #
    #      3----4               3--7--4
    #      |    |               |  |  |  ??
    #      |    |               6--9--8
    #      2----1               |  |  |
    #                           2--5--1
    #
    # Local number of the lagrange polynomials for the control points
    # entirely dependent on numbering of the 
    if Nn == 4:
        polynum  = np.array([[1,0],[0,0],[0,1],[1,1]])
    elif Nn == 9:
        polynum  = np.array([[2,0],[0,0],[0,2],[2,2],\
                            [1,0],[0,1],[1,2],[2,1],[1,1]])
         



    # GLL Points
    xi ,__ = gll.gll_pw(ngllx-1)  # points in x direction to be interpolated
    eta ,__ = gll.gll_pw(nglly-1) # points in y direction to be interpolated

    # Total number of refined GLL points
    NN = ngllx*nglly
    
    # Save shape function values in transform array
    gll_shape_matrix = np.zeros([NN,Nn])

    # Loop over local GLL points
    # GLL indexing and coordinates
    local_ind   = np.zeros([NN,1])
    local_coord = np.zeros([NN,2])
    
    counter = 0
    for j in range(nglly):
        for i in range(ngllx):
            # GLL point locations:
            local_coord[counter,:] = np.array([xi[i],eta[j]])
            counter+=1

    for i in range(NN):
        # Calculate each control point's shape function at each GLL
        # point
        for k in range(Nn):
            gll_shape_matrix[i,k] = \
                    gll.lagrange2D( polynum[k,0],local_coord[i,0],xi_control,\
                                    polynum[k,1],local_coord[i,1],xi_control)
   
    
    #######
    # Find global GLL points and numbering for EVERY element:

    # Empty connectivity array
    gll_connect = np.zeros([nel,NN],dtype='int')
    tol = 1e-10
    
    for i in range(nel):
        

        # Calculate global GLL coordinates for 1 element:
        glob_gll = np.array([gll_shape_matrix @ X[connect[i,:]].T, \
                             gll_shape_matrix @ Z[connect[i,:]].T]).T
            
        
        if i == 0: 
            # Create new element coordinates with NN number of nodes
            gll_coordinates   = glob_gll
            glob_gll_test     = glob_gll
            gll_connect[i,:]  = np.arange(0,NN)
        
        else: 
                
            for j in range(NN):
                """Looping over variable to check whether there is 
                overlap between the elements"""

                # Check whether coordinate exists already
                xbool = np.logical_and(\
                        np.abs(gll_coordinates[:,0]-glob_gll[j,0])<tol,\
                        np.abs(gll_coordinates[:,1]-glob_gll[j,1])<tol) 

                temp = np.where(xbool)
                ind = temp[0]
                
                # if it doesnt exist
                if ind.size==0:
                    # append new coordinates to global data set 
                    gll_coordinates = np.append(\
                                gll_coordinates, \
                                np.array([glob_gll[j,:]]),axis=0)
                    
                    # numbering at highest recent number +1
                    gll_connect[i,j] = np.amax(gll_connect) + 1
                    
                # if it already exists
                else:
                    # do NOT append a new coordinate, but use index for 
                    # numbering purpose
                    gll_connect[i,j] = ind
    

    return (gll_coordinates,gll_connect)





