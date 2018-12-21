"""This is script contains functions used to get from a quadrilateral
mesh to a Spectral-Element mesh. File for mesh object creation?"""

# Necessary for calculation and import of exodus file
import numpy as np
import netCDF4

# Plotting
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
import matplotlib

# Import GLL library to get the lagrange polynomials for interpolation
# of the grid
import gll_library as gll


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
    """mesh_interp(X,Y,Z,connect,ngllx,nglly)

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

def plot_elements(X,Y,connect,gll_coordinates):
    """plot_elements()

    This function plots the GLL points as well as control points in 2D

    """
    print(connect.shape)
    ##########
    # Calculate the Centre of elements
    try:
        nel,__ = connect.shape
    except ValueError:
        nel = 1
        connect = np.array([connect]) 
    
    # Finding element centres
    el_num_coor = np.zeros([nel,2])
    for i in range(nel):
        el_num_coor[i,0] = np.mean(X[connect[i,:]])
        el_num_coor[i,1] = np.mean(Y[connect[i,:]])

    
    # Creating polygons for each element
    xy = np.array([X[:], Y[:]]).T
    patches = []
    for coords in xy[connect[:]]:
        quad = Polygon(coords, True)
        patches.append(quad)
    
    ##########
    # Plotting 
    fig,ax = plt.subplots()

    # Plot Polygon Patches
    colors = 100 * np.random.rand(len(patches))
    p = PatchCollection(patches, cmap=matplotlib.cm.coolwarm,edgecolor="k", alpha=0.4)
    p.set_array(np.array(colors))
    ax.add_collection(p)

    # GLL Points
    ax.scatter(gll_coordinates[:,0],\
            gll_coordinates[:,1],50,color="k", marker="x")
    # Control Points
    ax.scatter(X, Y, 100, marker="o",
                          edgecolor="k",
                          color="None",
                          linewidth=2)
    # Plot element number
    for i in range(nel):    
        ax.text(el_num_coor[i,0], el_num_coor[i,1], "{0:d}".format(i+1),size=8,
             ha="center", va="center",
             bbox=dict(boxstyle="round",
                       ec=(1., 0.5, 0.5),
                       fc=(1., 0.8, 0.8),
                       )
             )
    plt.xlabel('x')
    plt.ylabel('y')
    ax.legend(['GLL Points', 'Control Points'])
    plt.title('Numbered elements, Nodes and GLL points')
    

def test_interp():
    """test_interp()
    
    Test the functions of this suite and plots them subsequently.
    """
    
    ngllx = 5
    ngllz = 5


    X,Y,Z,connect = readEx('../input/RectMesh.e')
    #print(X)
    #print(Z)
    #print(connect)
    #nel,__ = connect[0,:].shape 
    #print(nel)
    gll_coordinates, gll_connect = mesh_interp2D(X,Y,Z,connect,ngllx,ngllz)
    
    # Plotting First Element
    plot_elements(X[connect[0,:]],
            Z[connect[0,:]],connect[0,:],
            gll_coordinates[gll_connect[0,:],:]) 
    
    # Plotting All elements
    plot_elements(X,Z,connect,gll_coordinates)

    plt.show()

if __name__ == "__main__":
    test_interp()
    




