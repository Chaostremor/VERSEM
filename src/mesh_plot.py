"""This is script contains functions to plot a mesh

"""

# Necessary for calculation and import of exodus file
import numpy as np

# imports from mesh_spec
from src.mesh_spec import *


# Plotting
from matplotlib.patches import Polygon
from  matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
import matplotlib

# Import GLL library to get the lagrange polynomials for interpolation
# of the grid
import src.gll_library as gll



def plot_elements(X,Y,connect,gll_coordinates,gll_connect,\
                    num_o_el="all",plot_el_num="no",plot_node_num="no",\
                    plot_axis="yes"):
    """This function plots the GLL points as well as control points in 2D

    :param X: 1xN vectors with node coordinates
    
    :param Y: 1xN vectors with node coordinates
    
    :param connect: connectivity matrix depending on the number of GLL 
                    points, ([total number of elements] X [GLL^2])

    :param num_o_el: number of elements to be plotted, default is string 
                     valued 'all', set to number small

    :param plot_el_num: default is "no", change to "yes" if element numbers 
                        should be plotted
    :param plot_node_num: default is "no", change to "yes" if node numbers 
                          should be plotted
    :param plot_axis: default is "yes", if you want to omit axes 
                      set to "no" 

    :rtype: A figure plotting the 'num_o_el' 2D elements


    """
    
    # Change matrices depending on the number of elements to be printed
    if num_o_el=="all":
        # Number of elements, catch shape function error with one element 
        num_o_el,__ = connect.shape
    
    else:
        ## Control point change
        # Connectivity matrix 
        connect = connect[0:num_o_el,:]
        # number of used coordinates
        num_o_coor_n = np.max(connect)+1
        X = X[0:num_o_coor_n]
        Y = Y[0:num_o_coor_n]
        
        ## GLL point change
        # Connectivity matrix
        gll_connect = gll_connect[0:num_o_el,:]
        # number of used coordinates
        num_o_coor_gll = np.max(gll_connect)+1
        gll_coordinates= gll_coordinates[0:num_o_coor_gll,:] 

        
        
    ##########
    # Calculate the Centre of elements
    el_num_coor = np.zeros([num_o_el,2])
    for i in range(num_o_el):
        el_num_coor[i,0] = np.mean(X[connect[i,:]])
        el_num_coor[i,1] = np.mean(Y[connect[i,:]])


    ##########
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
            gll_coordinates[:,1],15,color="k", marker="x")
    # Control Points
    # alpha fill
    ax.scatter(X, Y, 50, marker="o",
                          edgecolor="None",
                          color="k",
                          linewidth=2,alpha=0.3)
    # outline
    ax.scatter(X, Y, 50, marker="o",
                          edgecolor="k",
                          color="None",
                          linewidth=2)
    
    # Plot element number
    if plot_el_num=="yes":
        # Plot element number
        for i in range(num_o_el):    
            ax.text(el_num_coor[i,0], el_num_coor[i,1], "{0:d}".format(i+1),size=8,
                 ha="right", va="top",
                 bbox=dict(boxstyle="round",
                           ec=(1., 0.5, 0.5),
                           fc=(1., 0.8, 0.8),
                           )
                 )

    # Plot element number
    if plot_node_num=="yes":
        # Plot element number
        for i in range(num_o_coor_gll):    
            ax.text(gll_coordinates[i,0], \
                    gll_coordinates[i,1], "{0:d}".format(i+1),size=8,
                    ha="left", va="bottom")

    plt.xlabel('x')
    plt.ylabel('y')
    plt.axis('equal')
    if plot_axis=="no":
        plt.axis("off")
    else:
        ax.legend(['GLL Points', 'Control Points'])
        plt.title('Numbered elements, Nodes and GLL points')
    

def test_interp():
    """test_interp()
    
    Running the example and plots them subsequently.

    """
    
    ngllx = 5
    ngllz = 5


    X,Y,Z,connect = readEx('input/RectMesh.e')
    gll_coordinates, gll_connect = mesh_interp2D(X,Y,Z,connect,ngllx,ngllz)

    # Plotting 1 element
    plot_elements(X,Z,connect,gll_coordinates,gll_connect,num_o_el=1,\
                                plot_el_num="yes",plot_node_num="yes")

    # Plotting 2 elements
    plot_elements(X,Z,connect,gll_coordinates,gll_connect,num_o_el="all",\
                                plot_el_num="yes",plot_node_num="no",\
                                plot_axis = "no")
                           


    # Plotting All elements
    #plot_elements(X,Z,connect,gll_coordinates,gll_connect)


    plt.show()

if __name__ == "__main__":
    test_interp()
    




