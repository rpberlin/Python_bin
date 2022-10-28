#!/Library/Frameworks/Python.framework/Versions/3.10/bin/python3

import sys
import math
from parseGPX import parseGPX
from plot_PMdata import calcPMdata
from VAMify import VAMify
from compare2GPXroutes import compare2GPXroutes
#from compare2GPXroutes import getQtyAtS
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.ticker import StrMethodFormatter, NullFormatter
import numpy as np
import scipy

def squareUVduct(xmin, xmax, ymin, ymax, Nx, Ny, iesTable, wallReflectivity):
    xNodes = np.linspace(xmin, xmax, num=Nx+1)
    yNodes = np.linspace(ymin, ymax, num=Ny+1)
    fluence = np.zeros((Nx,Ny),float)
    maxSteps = 100
    DX = (xmax-xmin)/Nx
    DY = (ymax-ymin)/Ny


    nRays = len(iesTable)
    p0 = np.empty((0, 2), float)
    p1 = np.empty((0, 2), float)
    for iRay in range(nRays):
        x0 = iesTable[iRay,0]
        y0 = iesTable[iRay,1]
        Ix0 = iesTable[iRay,2]
        Iy0 = iesTable[iRay,3]
        i0,j0 = getijFromxy(xNodes,yNodes,x0,y0)
        faceID,xintcp,yintcp,i1,j1 = getNextFace(x0, y0, i0, j0, Ix0, Iy0, xNodes, yNodes, 0)
        p0 = np.append(p0,[[x0, y0]],axis=0)
        p1 = np.append(p1,[[xintcp, yintcp]],axis=0)
        iStep=0
        Imagsq0 = Ix0*Ix0+Iy0*Iy0
        Iresid2 = 1;
        while Iresid2 > .01:
        #while iStep < maxSteps:
            iStep = iStep+1
            x0 = xintcp
            y0 = yintcp
            i0 = i1
            j0 = j1
            faceID,xintcp,yintcp,i1,j1 = getNextFace(x0, y0, i0, j0, Ix0, Iy0, xNodes, yNodes, iStep)
            if xintcp < xmin or xintcp > xmax:
                print("whoooaaaa")
            if faceID ==0:
                if i1 > Nx-1:
                    i1 = i0
                    Ix0 = -1*Ix0*wallReflectivity
                    Iy0 = Iy0*wallReflectivity
                fluence[i1,j1]=fluence[i1,j1]+DY*abs(Ix0)
            if faceID ==1:
                if j1 > Ny-1:
                    j1 = j0
                    Iy0 = -1*Iy0*wallReflectivity
                    Ix0 = Ix0*wallReflectivity
                fluence[i1,j1]=fluence[i1,j1]+DX*abs(Iy0)
            if faceID ==2:
                if i1 < 0:
                    i1 =i0
                    Ix0 = -1*Ix0*wallReflectivity
                    Iy0 = Iy0*wallReflectivity
                fluence[i1,j1]=fluence[i1,j1]+DY*abs(Ix0)
            if faceID ==3 :
                if j1 < 0:
                    j1 = j0
                    Iy0 = -1*Iy0*wallReflectivity
                    Ix0 = Ix0*wallReflectivity
                fluence[i1,j1]=fluence[i1,j1]+DX*abs(Iy0)
            #p1 = np.append(p1,[[xintcp, yintcp]],axis=0)
            Imagsqi = Ix0*Ix0+Iy0*Iy0
            Iresid2 = Imagsqi/Imagsq0
        print('iRay ',iRay,'iSteps ', iStep)

    #print('p0')
    #print(p0)
    #print('p1')
    #print(p1)
    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    line1 = plt.plot(p0[:,0],p0[:,1],'ro',label='P0')
   # line2 = plt.plot(p1[:,0],p1[:,1],'-+',label='P1')
    line3 = plt.plot(xNodes,0*xNodes+ymin*np.ones(Nx+1,float),'-o',label='xNodes')
    line4 = plt.plot(0*yNodes+xmin*np.ones(Ny+1,float),yNodes,'-o',label='yNodes')
    plt.legend()
    #plt.xlim([xmin,xmax])
    #plt.ylim([ymin,ymax])
    ax1 = fig.add_subplot(122)
    plt.contourf(np.log(fluence.T))

    plt.show()
    print('done')


        #print('i0, j0 =',i0,', ',j0, ' faceID: ', faceID,' i1: ',i1,' j1 ',j1)



def getNextFace(x0,y0,i0,j0,Ix0,Iy0,xNodes,yNodes,iStep):
    dxe = xNodes[i0+1]-x0
    dxw = xNodes[i0]-x0
    dyn = yNodes[j0+1]-y0
    dys = yNodes[j0]-y0
    rayTheta = math.atan2(Ix0,Iy0)

    cross0 = Ix0*dyn - Iy0*dxe   #0 - NorthEast
    cross1 = Ix0*dyn - Iy0*dxw   #1 - NorthWest
    cross2 = Ix0*dys - Iy0*dxw   #2 - SouthWest
    cross3 = Ix0*dys - Iy0*dxe   #3 - SouthEast


    if  cross0 >= 0 and cross3 < 0:
        faceID = 0
        i1 = i0+1
        j1 = j0
        xintcp = x0+dxe
        yintcp = y0+dxe*Iy0/Ix0
    elif cross1 >= 0 and cross0 < 0:
        faceID = 1
        i1 = i0
        j1 = j0+1
        xintcp = x0+dyn*Ix0/Iy0
        yintcp = y0+dyn
    elif cross2 >= 0 and cross1 < 0:
        faceID = 2
        i1 = i0-1
        j1 = j0
        xintcp = x0+dxw
        yintcp = y0+dxw*Iy0/Ix0
    elif cross3 >=0 and cross2 < 0:
        faceID = 3
        i1 = i0
        j1 = j0-1
        xintcp = x0+dys*Ix0/Iy0
        yintcp = y0+dys
    else:
        print('No Face Found')
    #print(rayTheta,cross0,cross1,cross2,cross3, faceID)
    if xintcp < xNodes[0] or xintcp > xNodes[-1] or yintcp < yNodes[0] or yintcp > yNodes[-1]:
        print('OUT OF BOUNDS ','faceID ',faceID,'',' xintcp ', xintcp, ' i0 ',i0, ' i1 ',i1, ' x0 ',x0, ' y0 ',y0, ' dxw ', dxw, 'dxe ',dxe)
    return faceID, xintcp,yintcp, i1,j1


def circleIESdata(I,centroid,radius,nPts):
    iesTable = np.empty((0, 4), float)
    for i in range(nPts):
        theta = i*2*math.pi/nPts+.1
        xloc = centroid[0]+radius*math.cos(theta)
        yloc = centroid[1]+radius*math.sin(theta)
        Ix = I*math.cos(theta)/nPts
        Iy = I*math.sin(theta)/nPts
        iesTable = np.append(iesTable,[[xloc, yloc, Ix, Iy]],axis=0)
        #print('iesTable = ',iesTable)
    return iesTable


def getijFromxy(xNodes,yNodes, x,y):
    i = int((x-xmin)/(xNodes[1]-xNodes[0]))
    j = int((y-ymin)/(yNodes[1]-yNodes[0]))
    return i,j


if __name__ == '__main__':
    n = len(sys.argv)
    print("Total arguments passed:", n)
    for i in range(1, n):
        print(i," ",sys.argv[i])

    if n < 2:
        Io = 10
        nPts = 100
        centroid = [.15,.05]
        radius = .02
        iesData = circleIESdata(Io,centroid,radius,nPts)

        xmin = -.3
        xmax = .3
        ymin = -.1
        ymax = .1
        Nx = 100
        Ny = 30
        sigma = 0.9
        squareUVduct(xmin,xmax,ymin,ymax, Nx, Ny, iesData, sigma)
