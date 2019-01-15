""" This contains the functions for ploting the wavefield snapshot
    
    Author: Chao Song
    
"""
 
# Necessary for calculation 
import numpy as np

# Necessary for Plotting
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Necessary for interpolation
from scipy import interpolate

# Necessary for retrieving the result file
import os

# Function for plotting 
def plt_snap(xnpt, dtype='u', stype='xdim', fname=None):
    """.. function:: plt_snap(xnpt, dtype='u', stype='xdim', fname=None)
    
    This function is used to read the selected data and plot the snapshot on
    preferred dimension.
    
    :param xnpt: integer number of points along 1st dimension after 
                 interpolation
    :param dtype: data type to read, string, 'u', 'v' or 'a'
    :param stype: type of which dimension to show, string, 'xdim', 'ydim', or
                  'all'
    :param fname: output file name of the snapshot, string, default is None
    
    :rtype: return the object of plotted animation, ``ani``

    This function can be called with proper parameters, and also can be 
    run directly from command line with default settings, in VERSEM, type:

    ::

        python scripts/pltsnapfunc.py

    for plotting displacement wavefield snapshot animation of the 1st
    dimension with 500 points along the 1st dimension after interpolation.
    The resulting figure is saved as ``u_xdim_xnpt500.mp4``. 
    
    """

    # load global gll points coordinates
    crd = np.load('results/gll_coordinates.npy')
    crdx = crd[:, 0]
    crdy = crd[:, 1]
    
    # get the number of gll points
    ngll = len(crdx)

    # get the range of showing region 
    xmin = np.amin(crdx)
    xmax = np.amax(crdx)
    ymin = np.amin(crdy)
    ymax = np.amax(crdy)
    
    ## get the file list of results
    file_list = []
    
    # if data type is displacement
    if dtype == 'u':        
        for file in os.listdir('results/timesteps'):
            if file.startswith("u"):
                file_list.append(os.path.join("results/timesteps", file))
        file_list.sort()
    
    # if data type is velocity        
    if dtype == 'v':
        for file in os.listdir('results/timesteps'):
            if file.startswith("v"):
                file_list.append(os.path.join("results/timesteps", file))
        file_list.sort()

    # if data type is acceleration    
    if dtype == 'a':
        for file in os.listdir('results/timesteps'):
            if file.startswith("a"):
                file_list.append(os.path.join("results/timesteps", file))
        file_list.sort()

    # get the number of total time steps
    nt = len(file_list)
    
    ## load the data file
    
    # if want to show 1st dimension
    if stype == 'xdim':
        datao = np.zeros([nt, ngll])
        for i in range(nt):
            temp = np.load(file_list[i])
            datao[i, :] = temp[0::2]
    
    # if want to show 2nd dimension            
    if stype == 'ydim':
        datao = np.zeros([nt, ngll])
        for i in range(nt):            
            temp = np.load(file_list[i])
            datao[i, :] = temp[1::2]
            
    # if want to show the scalar norm
    if stype == 'all':
        dataxo = np.zeros([nt, ngll])
        datayo = np.zeros([nt, ngll])
        datao = np.zeros([nt, ngll])
        for i in range(nt):
            temp = np.load(file_list[i])
            dataxo[i, :] = temp[0::2]
            datayo[i, :] = temp[1::2]
            
        # scalar norm is the sum of square
        datao = np.sqrt(np.add(np.square(dataxo),np.square(datayo)))      
    
    # get the ratio of axes, and multiply with float to make it float
    ratio = (1.0* (ymax-ymin)) / (xmax-xmin)
    ynpt = np.int(xnpt * ratio)

    # linearly interpolate the location points to denser
    xnew = np.linspace(xmin, xmax, num=xnpt)
    ynew = np.linspace(ymin, ymax, num=ynpt)

    # initiate the new data array    
    datan = np.zeros([nt, ynpt-1, xnpt-1])

    # make a mesh with new coordinates, 'xy' indexing leads to inversed size
    # of array
    xng, yng = np.meshgrid(xnew, ynew, indexing = 'xy')
    
    # interpolation to the data points with new coordinates, use 0 for points 
    # outside the original region
    for i in range(nt):
        tmp = interpolate.griddata(crd, datao[i, :], (xng, yng), fill_value=0)
    
        # xng and yng are bounds, so tmp should be the value *inside* those 
        # bounds. Therefore, remove the last value from the tmp array.
        tmp = tmp[:-1, :-1]    
        datan[i, :, :] = tmp

    # initiate the figure   
    fig = plt.figure(figsize=(10, 10*ratio))
    ax = fig.add_subplot(111)

    # set the axis information
    plt.xlabel('X (m)')
    plt.ylabel('Y (m)')
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    ax.set_aspect('equal')
    ax.autoscale(enable=True, axis='y', tight=True)

    # set the colormap
    cmap = plt.get_cmap('bwr')

    # set the value range of colormap with percentile of data    
    crange = np.percentile(np.abs(datao), 95)
    
    # plot the first snapshot as an initiation
    im = plt.pcolormesh(xng, yng, datan[0, :, :], cmap=cmap, \
                    vmin=-crange, vmax=crange)
    
    # show the colorbar
    plt.colorbar(im)

    # define the function of animation
    def animate(i):
        """..function:: animate(i)
        
        keep alternating data in called function at each step to make animation
        
        ::param i: index denoting the time step
        
        :rtype: return the refreshed frame, ``im``
        
        """
        
        # set the format of time 
        time_form = "Time = %f (s)"
        
        # only alternate the data in the figure every time step
        im.set_array(datan[i, :, :].ravel())
        
        # set the range of colorbar
        im.set_clim(-crange, crange)
        
        # show the time change as well
        temp1 = file_list[i]
        temp2 = temp1.split('/')
        temp3 = temp2[-1]
        time = float(temp3[1:-4])
        plt.title(time_form % (time))
        
        return im

    # make the animation by keep calling function animate
    # interval:     time interval between frames. ms
    # repeat_delay: time interval between repeat, ms
    # save_count:   number of frames to save in cache temporarily
    ani = animation.FuncAnimation(fig, animate, range(nt), interval=100,\
                              repeat_delay=1000, save_count=nt+1) 

    # set the output figure name and path
    outdir = 'results/'
    if fname is None:
        fname = dtype + '_' + stype + '_' + 'xnpt' + str(xnpt) + '_snap.mp4'
    fpath = outdir+fname

    # Set up formatting for the movie files
    # writer: decipher type
    # fps: number of frames per second
    # dpi: resolution
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=60, metadata=dict(artist='Me'), bitrate=1800)

    # save the animation
    ani.save(fpath, writer=writer, dpi=300)

    # open the switch of showing the figure
    plt.show()
    
    
    return ani


# for default one command:
#    python pltsnapshotfunc.py
if __name__ == "__main__":
    plt_snap(500, dtype='u', stype='xdim')
































