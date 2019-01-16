""" This contains the function for ploting the seismograms at specified
    stations at preferred time period with different data type
"""

# Necessary for calculation 
import numpy as np

# Necessary for Plotting
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

# Necessary for interpolation
from scipy import interpolate

# Necessary for retrieving the result file
import os

def plt_seismo(xnew, ynew, stime, etime, times, dtype='u', fname=None):
    """.. function:: plt_seismo(xnew, ynew, stime, etime, times, dtype='u',
                                fname=None)
    
    This function is used to read the selected data and plot the seismogram at
    each given station point at preferred time period
    
    :param xnew: 1D ``numpy`` array of 1st dimension corrdinates of desired
                 stations
    :param ynew: 1D ``numpy`` array of 2nd dimension corrdinates of desired
                 stations
    :param stime: ``numpy`` float start time of period to zoom in 
    :param etime: ``numpy`` float end time of period to zoom in
    :param times: ``numpy`` float or integer amplification factor along time
                  to interpolate
    :param dtype: data type to read, string, 'u', 'v' or 'a'
    :param fname: output file name of the seismogram, string, default is None

    This function can be called with proper parameters, and also can be 
    run directly from command line with default settings, in VERSEM, type:

    ::

        python scripts/pltseismofunc.py

    for plotting displacement seismograms at stations on the surface from 0
    to 2500 seperated by 100. The zoomed-in time is from 4.32 to 4.38 s.
    The resulting figure is saved as ``u_4.32-4.38s_seismo.pdf``.

    """  
  
    # load global gll points coordinates
    crd = np.load('results/gll_coordinates.npy')
    crdx = crd[:, 0]
    crdy = crd[:, 1]

    # get the number of gll points
    ngll = len(crdx)

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
    
    # load the data file
    dataxo = np.zeros([nt, ngll])
    datayo = np.zeros([nt, ngll])
    for i in range(nt):
        temp = np.load(file_list[i])
        dataxo[i, :] = temp[0::2]
        datayo[i, :] = temp[1::2]
        
    # scalar norm is the sum of square
    datao = np.sqrt(np.add(np.square(dataxo),np.square(datayo)))      
    
    # initiate the new data array
    npt = len(xnew)
    dxn = np.zeros([nt, npt])
    dyn = np.zeros([nt, npt])

    # interpolation to the data points with new coordinates, use 0 for points 
    # outside the original region
    for i in range(0,nt,100):
        # interpolation to 1st, 2nd dimension, respectively
        tmp1 = interpolate.griddata(crd, dataxo[i, :], (xnew, ynew), fill_value=0)
        dxn[i, :] = tmp1
        
        tmp2 = interpolate.griddata(crd, datayo[i, :], (xnew, ynew), fill_value=0)
        dyn[i, :] = tmp2

    # scalar norm is the sum of square        
    dan = np.sqrt(np.add(np.square(dxn), np.square(dyn)))
 
    # get the original time array
    t = np.zeros([nt])
    for i in range(nt):
        temp1 = file_list[i]
        temp2 = temp1.split('/')
        temp3 = temp2[-1]
        t[i] = float(temp3[1:-4])
    
    # set up new time array 
    ntn = nt*times
    tn = np.linspace(np.min(t), np.max(t), num=ntn)
    dxf = np.zeros([ntn, npt])
    dyf = np.zeros([ntn, npt])
    for i in range(npt):
        # interpolate to 1st dimension data along time 
        f = interpolate.interp1d(t, dxn[:, i], kind='cubic', \
                                 bounds_error=None, fill_value=0)
        tmp = f(tn)
        dxf[:, i] = tmp
        
        # interpolate to 1st dimension data along time 
        f = interpolate.interp1d(t, dyn[:, i], kind='cubic', \
                                 bounds_error=None, fill_value=0)
        tmp = f(tn)
        dyf[:, i] = tmp
     

    # initiate figure
    fig = plt.figure(figsize=(15, 9))
    
    # axis settings for ax1
    ax1 = fig.add_subplot(221)
    ax1.set_ylabel('Time (s)')
    ax1.set_xlabel('Trace number')
    tit1 = dtype + ' x'
    ax1.set_title(tit1)
    ax1.set_ylim(0, np.max(tn))
    
    # step is the inteval between traces
    step = 2*np.arange(0,npt,1)
    ax1.set_xlim(step[0]-1, step[-1]+1)
    ax1.invert_yaxis()
    
    # normalized and plot 1st dimension with all time 
    dxfn = np.zeros([ntn, npt]) 
    for i in range(npt):
        dxfn[:, i] = dxf[:, i] / np.max(np.abs(dxf[:, i]))
        plt.plot(dxfn[:, i]+step[i], tn, color='blue', linewidth=2)

    # axis settings for ax2
    ax2 = fig.add_subplot(223)
    ax2.set_ylabel('Time (s)')
    ax2.set_xlabel('Trace number')
    tit2 = 'Zoomed-in ' + dtype + ' x'
    ax2.set_title(tit2)
    ax2.set_ylim(stime, etime)
    ax2.set_xlim(step[0]-1, step[-1]+1)
    ax2.invert_yaxis()
    
    # normalized and plot 1st dimension with specified time 
    dxfn = np.zeros([ntn, npt]) 
    for i in range(npt):
        dxfn[:, i] = dxf[:, i] / np.max(np.abs(dxf[:, i]))
        plt.plot(dxfn[:, i]+step[i], tn, color='blue', linewidth=2)


    # axis settings for ax3
    ax3 = fig.add_subplot(222)
    ax3.set_ylabel('Time (s)')
    ax3.set_xlabel('Trace number')
    tit3 = dtype + ' y'
    ax3.set_title(tit3)
    ax3.set_ylim(0, np.max(tn))
    ax3.set_xlim(step[0]-1, step[-1]+1)  
    ax3.invert_yaxis()

    # normalized and plot 2nd dimension with all time   
    dyfn = np.zeros([ntn, npt])  
    for i in range(npt):
        dyfn[:, i] = dyf[:, i] / np.max(np.abs(dyf[:, i]))
        plt.plot(dyfn[:, i]+step[i], tn, color='blue', linewidth=2)

    # axis settings for ax4        
    ax4 = fig.add_subplot(224)
    ax4.set_ylabel('Time (s)')
    ax4.set_xlabel('Trace number')
    tit4 = 'Zoomed-in ' + dtype + ' y'
    ax4.set_title(tit4)
    ax4.set_ylim(stime, etime)
    ax4.set_xlim(step[0]-1, step[-1]+1)  
    ax4.invert_yaxis()

    # set the output dir
    outdir = 'results/seismograms/'
    if os.path.exists(outdir) == False:
            os.makedirs(outdir)

    # normalized and plot 2nd dimension with specified time   
    dyfn = np.zeros([ntn, npt])  
    for i in range(npt):
        dyfn[:, i] = dyf[:, i] / np.max(np.abs(dyf[:, i]))
        plt.plot(dyfn[:, i]+step[i], tn, color='blue', linewidth=2)
        # set up written data name and whole path
        wdata = np.array([tn, dyf[:, i]])
        wdnm = dtype + '_x' + str(xnew[i]) + '_y' + str(ynew[i]) + '.npy'
        wdpath = outdir + wdnm
        np.save(wdpath, wdata)

    # set up figure name and whole path
    if fname is None:
        fname = dtype + '_' + str(stime) + '-' + str(etime) + 's_seismo.pdf'
    fpath = outdir + fname
    
    # save figure as a pdf file
    # dpi: resolution
    plt.savefig(fpath, dpi=300, orientation='portrait')
    
    # open the switch of showing the figure
    plt.show()
    

if __name__ == "__main__":
    plt_seismo(np.linspace(0, 2500, 26), np.zeros([26]), 4.32, 4.38, 10, \
               dtype='u')
