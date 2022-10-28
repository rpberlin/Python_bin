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

def squareUVduct(xmin, xmax, ymin, ymax, Nx, Ny, rayData, wallReflectivity):
    nRays = len(rayData)
    xNodes = np.linspace(xmin, xmax, num=Nx+1)
    yNodes = np.linspace(ymin, ymax, num=Ny+1)
    fluence = np.zeros((Nx,Ny),float)
    delta_fluence = np.zeros((Nx,Ny),float)
    newRayData = np.zeros((nRays, 4), float)
    maxSteps = 500
    rayConveLimit = 1e-4
    rayConvLimit2 = rayConveLimit*rayConveLimit
    fluenceConvergenceLimit = 1e-6
    fluenceNorm0=1e6
    fluenceResid=1
    residPlot = np.empty((0, 2), float)
    nIts = 0
    nItsMax = 500

    DX = (xmax-xmin)/Nx
    DY = (ymax-ymin)/Ny

    p0 = np.empty((0, 2), float)
    p1 = np.empty((0, 2), float)
    faceIDcounts = np.zeros(4)
    while fluenceResid > fluenceConvergenceLimit*maxSteps:
        for iRay in range(nRays):
            x1 = rayData[iRay,0]
            y1 = rayData[iRay,1]
            Ix1 = rayData[iRay,2]
            Iy1 = rayData[iRay,3]
            i1,j1 = getijFromxy(xNodes,yNodes,x1,y1)
            faceID,x1,y1,i1,j1,Ix1,Iy1 = getNextFace(x1, y1, i1, j1, Ix1, Iy1, xNodes, yNodes,wallReflectivity)
            if faceID ==0 or faceID == 2:
                delta_fluence[i1,j1]=delta_fluence[i1,j1]+DY*abs(Ix1)
            if faceID ==1 or faceID == 3:
                delta_fluence[i1,j1]=delta_fluence[i1,j1]+DX*abs(Iy1)
            iStep=0
            Imagsq0 = Ix1*Ix1+Iy1*Iy1
            Iresid2 = 1;
            #while Iresid2 > rayConvLimit2:
            while iStep < maxSteps:
                iStep = iStep+1
                faceID,x1,y1,i1,j1, Ix1,Iy1 = getNextFace(x1, y1, i1, j1, Ix1, Iy1, xNodes, yNodes,wallReflectivity)
                if faceID ==0 or faceID == 2:
                    delta_fluence[i1,j1]=delta_fluence[i1,j1]+DY*abs(Ix1)
                if faceID ==1 or faceID == 3:
                    delta_fluence[i1,j1]=delta_fluence[i1,j1]+DX*abs(Iy1)
                #p1 = np.append(p1,[[x1, y1]],axis=0)
                faceIDcounts[faceID] +=1
                Imagsqi = Ix1*Ix1+Iy1*Iy1
                Iresid2 = Imagsqi/Imagsq0
            #print('faceCounts: ',faceIDcounts)


            newRayData[iRay,0]=x1
            newRayData[iRay,1]=y1
            newRayData[iRay,2]=Ix1
            newRayData[iRay,3]=Iy1
        fluenceNorm = scipy.linalg.norm(delta_fluence)
        if nIts == 0:
            fluenceNorm0 = fluenceNorm
        nIts=nIts+1
        fluenceResid = fluenceNorm/fluenceNorm0
        print(nIts,fluenceResid)
        residPlot = np.append(residPlot,[[nIts, fluenceResid]],axis=0)

        fluence = fluence + delta_fluence
        delta_fluence = 0*delta_fluence
        rayData = newRayData



    xgrid = np.zeros((Nx+1,Ny+1),float)
    ygrid = np.zeros((Nx+1,Ny+1),float)
    fluenceAtnodes = np.zeros((Nx+1,Ny+1),float)

    for i in range(Nx+1):
        for j in range(Ny+1):
            xgrid[i,j]=xNodes[i]
            ygrid[i,j]=yNodes[j]
            imin = min(Nx-1,i)
            imax = min(Nx-1,i+1)
            jmin = min(Ny-1,j)
            jmax = min(Ny-1,j+1)
            f0 = fluence[imin,jmin]
            f1 = fluence[imax,jmin]
            f2 = fluence[imin,jmax]
            f3 = fluence[imax,jmax]
            fluenceAtnodes[i,j]=0.25*(f0+f1+f2+f3)




    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    #line1 = plt.plot(p0[:,0],p0[:,1],'ro',label='P0')
    #line2 = plt.plot(p1[:,0],p1[:,1],'-+',label='P1',linewidth=1)
    #line3 = plt.plot(xNodes,0*xNodes+ymin*np.ones(Nx+1,float),'-o',label='xNodes')
    #line4 = plt.plot(0*yNodes+xmin*np.ones(Ny+1,float),yNodes,'-o',label='yNodes')
    plt.semilogy(residPlot[:,0],residPlot[:,1],label='L2 Fluence Residual')
    plt.legend()
    #plt.xlim([xmin,xmax])
    #plt.ylim([ymin,ymax])
    ax1 = fig.add_subplot(122)
    #plt.contourf(np.log(fluence.T))
    plt.contourf(xgrid,ygrid,fluenceAtnodes)


    plt.show()
    print('done')
    return fluence


        #print('i0, j0 =',i0,', ',j0, ' faceID: ', faceID,' i1: ',i1,' j1 ',j1)



def getNextFace(x0,y0,i0,j0,Ix0,Iy0,xNodes,yNodes,wallReflectivity):
    Nx = len(xNodes)-1
    Ny = len(yNodes)-1
    if i0 == Nx or y0 == Ny or i0 == -1 or j0 == -1:
        print('whooaaaa')
    dxe = xNodes[i0+1]-x0
    dxw = xNodes[i0]-x0
    dyn = yNodes[j0+1]-y0
    dys = yNodes[j0]-y0
    rayTheta = math.atan2(Ix0,Iy0)

    cross0 = Ix0*dyn - Iy0*dxe   #0 - NorthEast
    cross1 = Ix0*dyn - Iy0*dxw   #1 - NorthWest
    cross2 = Ix0*dys - Iy0*dxw   #2 - SouthWest
    cross3 = Ix0*dys - Iy0*dxe   #3 - SouthEast

    Ix1=Ix0
    Iy1=Iy0
    if  cross0 >= 0 and cross3 < 0:
        faceID = 0
        i1 = i0+1
        j1 = j0
        xintcp = x0+dxe
        yintcp = y0+dxe*Iy0/Ix0
        if i1 == Nx:
            i1 = i0
            Ix1 = -1*Ix0*wallReflectivity
            Iy1 = Iy0*wallReflectivity
    elif cross1 >= 0 and cross0 < 0:
        faceID = 1
        i1 = i0
        j1 = j0+1
        xintcp = x0+dyn*Ix0/Iy0
        yintcp = y0+dyn
        if j1 == Ny:
            j1 = j0
            Iy1 = -1*Iy0*wallReflectivity
            Ix1 = Ix0*wallReflectivity
    elif cross2 >= 0 and cross1 < 0:
        faceID = 2
        i1 = i0-1
        j1 = j0
        xintcp = x0+dxw
        yintcp = y0+dxw*Iy0/Ix0
        if i1 < 0:
            i1 = 0
            Ix1 = -1*Ix0*wallReflectivity
            Iy1 = Iy0*wallReflectivity
    elif cross3 >=0 and cross2 < 0:
        faceID = 3
        i1 = i0
        j1 = j0-1
        xintcp = x0+dys*Ix0/Iy0
        yintcp = y0+dys
        if j1 < 0:
            j1 = 0
            Iy1 = -1*Iy0*wallReflectivity
            Ix1 = Ix0*wallReflectivity
    else:
        print('No Face Found')
    #print(rayTheta,cross0,cross1,cross2,cross3, faceID)
    if xintcp < xNodes[0] or xintcp > xNodes[-1] or yintcp < yNodes[0] or yintcp > yNodes[-1]:
        print('OUT OF BOUNDS ','faceID ',faceID,'',' xintcp ', xintcp, ' i0 ',i0, ' i1 ',i1, ' x0 ',x0, ' y0 ',y0, ' dxw ', dxw, 'dxe ',dxe)
    return faceID, xintcp,yintcp, i1,j1, Ix1, Iy1


def circleRaymaker(I,centroid,radius,nPts):
    rayTable = np.empty((0, 4), float)
    for i in range(nPts):
        theta = i*2*math.pi/nPts+.2
        xloc = centroid[0]+radius*math.cos(theta)
        yloc = centroid[1]+radius*math.sin(theta)
        Ix = I*math.cos(theta)/nPts
        Iy = I*math.sin(theta)/nPts
        rayTable = np.append(rayTable,[[xloc, yloc, Ix, Iy]],axis=0)
        #print('iesTable = ',iesTable)
    return rayTable


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
        nPts = 10000
        centroid = [0.0,0.0]
        radius = 0.3
        rayData = circleRaymaker(Io,centroid,radius,nPts)
        xmin = -1
        xmax = 1
        ymin = -1
        ymax = 1
        Nx = 100
        Ny = 100
        sigma = .8
        fluence = squareUVduct(xmin,xmax,ymin,ymax, Nx, Ny, rayData, sigma)
        print(fluence)
