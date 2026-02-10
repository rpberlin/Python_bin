import matplotlib.pyplot as plt
import numpy as np
from pint import UnitRegistry
Q_ = UnitRegistry.Quantity

def calc_laplace(arr, stencil):
    arr_lapl = 0*arr 
    for i in range(1,Nx-1):
        for j in range(1,Ny-1):
            arr_lapl[i,j] = np.sum(arr[i-1:i+2, j-1:j+2] * stencil)
            #arr_lapl[i,j] = arr[i-1,j] + arr[i+1,j] + arr[i,j+1] + arr[i,j-1] - 4*arr[i,j]

    arr_lapl[0,1:-1] = 2*arr_lapl[1,1:-1] - arr_lapl[2,1:-1]
    arr_lapl[-1,  1:-1] = 2*arr_lapl[-2,  1:-1] - arr_lapl[-3,  1:-1] 
    arr_lapl[1:-1, 0 ] = 2*arr_lapl[1:-1, 1 ] - arr_lapl[1:-1, 2 ]    
    arr_lapl[1:-1, -1 ] = 2*arr_lapl[1:-1, -2 ] - arr_lapl[1:-1, -3 ] 

    arr_lapl[0,  0 ]  = 0.5*((2*arr_lapl[1,  0]-arr_lapl[2,  0]) + (2*arr_lapl[0,  1]-arr_lapl[0,  2]))        
    arr_lapl[-1, 0 ]  = 0.5*((2*arr_lapl[-2, 0]-arr_lapl[-3, 0]) + (2*arr_lapl[-1, 1]-arr_lapl[-1, 2]))      
    arr_lapl[0,  -1]  = 0.5*((2*arr_lapl[1, -1]-arr_lapl[2, -1]) + (2*arr_lapl[0, -2]-arr_lapl[0, -3]))      
    arr_lapl[-1, -1]  = 0.5*((2*arr_lapl[-2,-1]-arr_lapl[-3,-1]) + (2*arr_lapl[-1,-2]-arr_lapl[-1,-3]))
    return arr_lapl 

Fo = .1
URF = 0.3
N_steps = 2000
R_seed = Q_(2,'cm')
R_spike = Q_(.3,'cm')
Lx = Q_(10,'cm')
Ly = Q_(10,'cm')
x0 = [-.5*Lx,-.5*Ly]
Nx = 40
Ny = 40
dx = Lx/Nx
dy = Ly/Ny
dt = Fo*dx.m**2

arr = np.zeros((Nx, Ny))
eps = 1e-6
x_centers = x0[0] + (np.arange(Nx) + 0.5) * dx
y_centers = x0[1] + (np.arange(Ny) + 0.5) * dy
centroids = np.zeros((Nx,Ny,2))
Xc, Yc = np.meshgrid(x_centers, y_centers, indexing="ij")
centroids = np.stack((Xc, Yc), axis=-1)

#mask = ((centroids[:, :, 0] < Q_(2,'cm')) & (centroids[:, :, 1] < Q_(1,'cm')) & (centroids[:, :, 1] > Q_(-1,'cm')) & (centroids[:, :, 0] > Q_(-2,'cm')))
#mask = (centroids[:,:,0] < dy)
#mask = (centroids[:,:,0]**2 + centroids[:,:,1]**2 < Q_(4,'cm^2'))
mask = (centroids[:,:,0]**2 + centroids[:,:,1]**2 < R_seed**2)  | ((centroids[:,:,0]-R_seed)**2 + centroids[:,:,1]**2 < R_spike**2 )
arr[mask] = 1
arr0 = arr 

compact_stencil = np.array([[1,4,1],[4,-20,4],[1,4,1]])
#compact_stencil = np.array([[0,1,0],[1,-4,1],[0,1,0]])
#compact_stencil = np.array([[1,2,1],[2,-12,2],[1,2,1]])
print(compact_stencil)
for i in range(0,N_steps):
    arr_lapl = -1*calc_laplace(arr,compact_stencil)
    dphi_dt = (1.0-arr)*arr_lapl
    dphi = dphi_dt*dt
    arr_star = arr + dphi
    arr = arr*(1-URF) + URF*arr_star
    arr[arr > 1] =1
    arr[arr <0] = 0

    if i%20 == 0:
        plt.figure()
        plt.contourf(arr)
        plt.colorbar()
        plt.axis('square')
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.title("Contour Plot of arr")
        plt.show()
