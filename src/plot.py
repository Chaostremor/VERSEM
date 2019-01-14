import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import os
import sys


def main(frames,dim='X'):
    """Function plots quickly the results in the directory
    results.

    :param dim: diplacement dimension to be plotted. 'X' or 'Y'
    :param frames: total number of frames to be plotted over the 
                total number of time steps time

    Can be run directly from command line with

    ::

        python -m src.plot 50 "X"

    for 50 frames total and displacement in x direction.

    """

    # loading coordinates
    gll_coordinates = np.load("results/gll_coordinates.npy")

    file_list = []

    # Displacement filenames
    for file in os.listdir("results/timesteps"):
        if file.startswith("u"):
            file_list.append(os.path.join("results/timesteps", file))

    # number of characters to delete
    nc = len("results/timesteps")+2

    file_list.sort()
    N = len(file_list)

    # Empty matrices
    ux = np.zeros([len(gll_coordinates),int(frames)])
    uy = np.zeros([len(gll_coordinates),int(frames)])
    new_file_list = []

    counter = 0
    for i in range(0,N,int(N/frames)):

        # Get new,shorter file list
        new_file_list.append(file_list[i])

        # Load and assign values
        u = np.load(file_list[i])
        ux[:,counter] = u[::2]
        uy[:,counter] = u[1::2]

        counter += 1
            

    # scaling factor
    c = 0.15        
    uymin,uymax  = np.amin(uy)*c,np.amax(uy)*c
    uxmin,uxmax  = np.amin(ux)*c,np.amax(ux)*c

    umin,umax = -np.sqrt(uxmin**2 + uymin**2),np.sqrt(uxmax**2 + uymax**2)

    xmin,xmax = np.amin(gll_coordinates[:,0]),np.amax(gll_coordinates[:,0])
    ymin,ymax = np.amin(gll_coordinates[:,1]),np.amax(gll_coordinates[:,1])

    XX = np.linspace(xmin,xmax,1000)
    YY = np.linspace(ymin,ymax,500)

    x,y = np.meshgrid(XX,YY,indexing='xy')
    
    plt.ion()
    plt.axis('equal')

    
    
    #ke = np.array([])

    

    for i in range(frames):
        plt.clf()
        
        if dim=="X":
            # X
            U = interpolate.griddata(gll_coordinates,ux[:,i],(x,y),fill_value=9999)
            plt.pcolormesh(x,y,U,vmin=uxmin,vmax=uxmax)
            plt.axis('equal')
            plt.axis('off')
            plt.title("Time: %10.8s sec" % new_file_list[i][nc:-4])
            plt.pause(0.01) 
            continue

        if dim=="Y":
            # Y
            U = interpolate.griddata(gll_coordinates,uy[:,i],(x,y),fill_value=9999)
            plt.pcolormesh(x,y,U,vmin=uymin,vmax=uymax)
            plt.axis('equal')
            plt.axis('off')
            plt.title("Time: %10.8s sec" % new_file_list[i][nc:-4])
            plt.pause(0.01) 
            continue

        

if __name__ == "__main__":
    print(sys.argv[1],sys.argv[2])
    main(int(sys.argv[1]),sys.argv[2])
