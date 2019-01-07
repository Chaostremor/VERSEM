import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import os

# loading coordinates
gll_coordinates = np.load("results/gll_coordinates.npy")

file_list = []

# Displacement filenames
for file in os.listdir("results/timesteps"):
    if file.startswith("u"):
        file_list.append(os.path.join("results/timesteps", file))

file_list.sort()
N = len(file_list)

ux = np.zeros([len(gll_coordinates),N])
uy = np.zeros([len(gll_coordinates),N])

for i in range(N):
        u = np.load(file_list[i])
        ux[:,i] = u[::2]
        uy[:,i] = u[1::2]
        
        # ux[:,i] = u[:len(gll_coordinates)]
        # uy[:,i] = u[len(gll_coordinates):]
uymin,uymax  = np.amin(uy)*.5,np.amax(uy)*.5
uxmin,uxmax  = np.amin(ux)*.5,np.amax(ux)*.5

umin,umax = -np.sqrt(uxmin**2 + uymin**2)*0.5,np.sqrt(uxmax**2 + uymax**2)*0.5

xmin,xmax = np.amin(gll_coordinates[:,0]),np.amax(gll_coordinates[:,0])
ymin,ymax = np.amin(gll_coordinates[:,1]),np.amax(gll_coordinates[:,1])

XX = np.linspace(xmin,xmax,1000)
YY = np.linspace(ymin,ymax,500)

x,y = np.meshgrid(XX,YY,indexing='xy')

plt.ion()
plt.axis('equal')

ke = np.array([])

for i in range(100,N,50):
    plt.clf()
    # print(file_list[i])
    U = interpolate.griddata(gll_coordinates,ux[:,i],(x,y),fill_value=9999)
    plt.pcolormesh(x,y,U,vmin=uxmin,vmax=uxmax)
    # U = interpolate.griddata(gll_coordinates,np.sqrt(ux[:,i]**2+uy[:,i]**2),(x,y),fill_value=9999)
    # plt.pcolormesh(x,y,U,vmin=umin,vmax=umax)
    plt.axis('equal')
    plt.axis('off')
    plt.title("%d" % i)
    plt.pause(0.01)

    ke = np.append(ke,np.sum(np.ndarray.flatten(U)))


plt.plot(ke)




