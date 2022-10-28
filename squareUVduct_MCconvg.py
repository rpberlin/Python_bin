#!/Library/Frameworks/Python.framework/Versions/3.10/bin/python3

import sys
import math
from parseGPX import parseGPX
from plot_PMdata import calcPMdata
from VAMify import VAMify
from compare2GPXroutes import compare2GPXroutes
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.ticker import StrMethodFormatter, NullFormatter
import numpy as np
import scipy
import random

def squareUVductMC(xmin, xmax, ymin, ymax,I_tot, Nx, Ny, wallReflectivity):

    xNodes = np.linspace(xmin, xmax, num=Nx+1)
    yNodes = np.linspace(ymin, ymax, num=Ny+1)
    fluence = np.zeros((Nx,Ny),float)
    delta_fluence = np.zeros((Nx,Ny),float)
    irrLR = np.zeros((Ny,2),float)
    irrUD  = np.zeros((Nx,2),float)
    delta_irrLR = np.zeros((Ny,2),float)
    delta_irrUD  = np.zeros((Nx,2),float)

    maxSteps = 2000
    rayConveLimit = 1e-4
    rayConvLimit2 = rayConveLimit*rayConveLimit
    fluenceConvergenceLimit = 1e-6
    fluenceNorm0=1e6
    fluenceResid=1
    residPlot = np.empty((0, 4), float)
    nIts = 0
    nItsMax = 10

    DX = (xmax-xmin)/Nx
    DY = (ymax-ymin)/Ny

    p0 = np.empty((0, 2), float)
    p1 = np.empty((0, 2), float)
    while fluenceResid > fluenceConvergenceLimit and nIts <nItsMax:
        rayData = circleRaymakerMC(I_tot,centroid,radius,nPts)
        nRays = len(rayData)
        for iRay in range(nRays):
            x1 = rayData[iRay,0]
            y1 = rayData[iRay,1]
            Ix1 = rayData[iRay,2]
            Iy1 = rayData[iRay,3]
            i1,j1 = getijFromxy(xNodes,yNodes,x1,y1)
            faceID,x1,y1,i1,j1,Ix1,Iy1,dF,dILR,dIUD = getNextFace(x1, y1, i1, j1, Ix1, Iy1, xNodes, yNodes,wallReflectivity)
            delta_fluence[dF[0],dF[1]] += dF[2]
            delta_irrLR[dILR[0],dILR[1]]+=dILR[2]
            delta_irrUD[dIUD[0],dIUD[1]]+=dIUD[2]

            iStep=0
            Imagsq0 = Ix1*Ix1+Iy1*Iy1
            Iresid2 = 1;
            while Iresid2 > rayConvLimit2:
                iStep = iStep+1
                faceID,x1,y1,i1,j1, Ix1,Iy1,dF,dILR,dIUD = getNextFace(x1, y1, i1, j1, Ix1, Iy1, xNodes, yNodes,wallReflectivity)
                delta_fluence[dF[0],dF[1]] += dF[2]
                delta_irrLR[dILR[0],dILR[1]]+=dILR[2]
                delta_irrUD[dIUD[0],dIUD[1]]+=dIUD[2]
                Imagsqi = Ix1*Ix1+Iy1*Iy1
                Iresid2 = Imagsqi/Imagsq0


        nIts=nIts+1
        fluenceStar = (1/(nIts+1))*(delta_fluence-fluence)
        irrLRStar = (1/(nIts+1))*(delta_irrLR-irrLR)
        irrUDStar = (1/(nIts+1))*(delta_irrUD-irrUD)
        irrLR = irrLR + irrLRStar
        irrUD = irrUD + irrUDStar
        fluence = fluence+fluenceStar
        fluenceNorm = scipy.linalg.norm(fluenceStar)
        irrLRNorm = scipy.linalg.norm(irrLRStar)
        irrUDNorm = scipy.linalg.norm(irrUDStar)
        if nIts == 1:
            fluenceNorm0 = fluenceNorm
            irrLRNorm0 = irrLRNorm
            irrUDNorm0 = irrUDNorm

        #fluenceResid = fluenceNorm/fluenceNorm0
        fluenceResid = fluenceNorm/scipy.linalg.norm(fluence)
        irrLRResid = irrLRNorm/irrLRNorm0
        irrUDResid = irrUDNorm/irrUDNorm0
        print(nIts,fluenceResid, irrLRResid,irrUDResid)
        residPlot = np.append(residPlot,[[nIts, fluenceResid, irrLRResid,irrUDResid]],axis=0)

        delta_fluence = 0*delta_fluence
        delta_irrLR = 0*delta_irrLR
        delta_irrUD = 0*delta_irrUD


    xgrid = np.zeros((Nx+1,Ny+1),float)
    ygrid = np.zeros((Nx+1,Ny+1),float)
    fluenceAtnodes = np.zeros((Nx+1,Ny+1),float)
    irrAtNodesLR = np.zeros((Ny+1,2),float)
    irrAtNodesUD = np.zeros((Nx+1,2),float)

    for i in range(Nx+1):
        imin = min(Nx-1,i)
        imax = min(Nx-1,i+1)
        irrAtNodesUD[i,0] = 0.5*(irrUD[imin,0]+irrUD[imax,0])
        irrAtNodesUD[i,1] = 0.5*(irrUD[imin,1]+irrUD[imax,1])
        for j in range(Ny+1):
            xgrid[i,j]=xNodes[i]
            ygrid[i,j]=yNodes[j]
            jmin = min(Ny-1,j)
            jmax = min(Ny-1,j+1)
            fluenceAtnodes[i,j]=0.25*(fluence[imin,jmin]+fluence[imax,jmin]+fluence[imin,jmax]+fluence[imax,jmax])
            irrAtNodesLR[j,0] = 0.5*(irrLR[jmin,0]+irrLR[jmax,0])
            irrAtNodesLR[j,1] = 0.5*(irrLR[jmin,1]+irrLR[jmax,1])


    fig = plt.figure()
    ax1 = fig.add_subplot(221)
    plt.semilogy(residPlot[:,0],residPlot[:,1],label='L2 Fluence Residual')
    plt.semilogy(residPlot[:,0],residPlot[:,2],label='LR Irradiance Residual')
    plt.semilogy(residPlot[:,0],residPlot[:,3],label='UD Irradiance Residual')
    plt.legend()
    ax2 = fig.add_subplot(222)
    plt.plot(xNodes,irrAtNodesUD[:,0],label='bottom')
    plt.plot(xNodes,irrAtNodesUD[:,1],label='top')
    plt.legend()

    ax2 = fig.add_subplot(223)
    plt.plot(yNodes,irrAtNodesLR[:,0],label='left')
    plt.plot(yNodes,irrAtNodesLR[:,1],label='right')
    plt.legend()

    ax3 = fig.add_subplot(224)
    #plt.contourf(np.log(fluence.T))
    plt.contourf(xgrid,ygrid,fluenceAtnodes)

    #print(irrUD,irrLR)
    plt.show()
    print('done')
    return fluence, irrLR,irrUD


def getNextFace(x0,y0,i0,j0,Ix0,Iy0,xNodes,yNodes,wallReflectivity):
    Nx = len(xNodes)-1
    Ny = len(yNodes)-1
    DY = (xNodes[1]-xNodes[0])
    DX = (yNodes[1]-yNodes[0])
    dF = [0,0,0]
    dILR = [0,0,0]
    dIUD = [0,0,0]

    if i0<0 or i0>=Nx or j0>=Ny or j0<0:
        print('Holt up')
    dxe = xNodes[i0+1]-x0
    dxw = xNodes[i0]-x0
    dyn = yNodes[j0+1]-y0
    dys = yNodes[j0]-y0

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
            dILR = [j1,1,Ix0*DY]
            i1 = i0
            Ix1 = -1*Ix0*wallReflectivity
            Iy1 = Iy0*wallReflectivity
        dF = [i1,j1,abs(Ix1)*DY]
    elif cross1 >= 0 and cross0 < 0:
        faceID = 1
        i1 = i0
        j1 = j0+1
        xintcp = x0+dyn*Ix0/Iy0
        yintcp = y0+dyn
        if j1 == Ny:
            dIUD = [i1,1,Iy0*DX]
            j1 = j0
            Iy1 = -1*Iy0*wallReflectivity
            Ix1 = Ix0*wallReflectivity
        dF = [i1,j1,abs(Iy1)*DX]
    elif cross2 >= 0 and cross1 < 0:
        faceID = 2
        i1 = i0-1
        j1 = j0
        xintcp = x0+dxw
        yintcp = y0+dxw*Iy0/Ix0
        if i1 < 0:
            dILR = [j1,0,-1.0*Ix0*DY]
            i1 = 0
            Ix1 = -1*Ix0*wallReflectivity
            Iy1 = Iy0*wallReflectivity
        dF = [i1,j1,abs(Ix1)*DY]
    elif cross3 >=0 and cross2 < 0:
        faceID = 3
        i1 = i0
        j1 = j0-1
        xintcp = x0+dys*Ix0/Iy0
        yintcp = y0+dys
        if j1 < 0:
            dIUD = [i1,0,-1.0*Iy0*DX]
            j1 = 0
            Iy1 = -1*Iy0*wallReflectivity
            Ix1 = Ix0*wallReflectivity
        dF = [i1,j1,abs(Iy1)*DX]
    else:
        print('No Face Found')
    #print(rayTheta,cross0,cross1,cross2,cross3, faceID)
    if xintcp < xNodes[0] or xintcp > xNodes[-1] or yintcp < yNodes[0] or yintcp > yNodes[-1]:
        print('OUT OF BOUNDS ','faceID ',faceID,'',' xintcp ', xintcp, ' i0 ',i0, ' i1 ',i1, ' x0 ',x0, ' y0 ',y0, ' dxw ', dxw, 'dxe ',dxe)
    return faceID, xintcp,yintcp, i1,j1, Ix1, Iy1, dF, dILR, dIUD


def circleRaymakerMC(I,centroid,radius,nPts):
    rayTable = np.empty((0, 4), float) #0 xloc 1:yloc 2:Ix 3:Iy
    dtheta = 2*math.pi/nPts
    for i in range(nPts):
        theta = i*dtheta+dtheta*random.uniform(0, 1)
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

    if n < 2:
        Io = 10
        nPts = 3000
        centroid = [0.5,0.2]
        radius = 0.2
        xmin = -1
        xmax = 1
        ymin = -0.6
        ymax = 0.9
        Nx = 400
        Ny = 200
        sigma = .2
        fluence, irrLR,irrUD = squareUVductMC(xmin,xmax,ymin,ymax,Io, Nx, Ny, sigma)
