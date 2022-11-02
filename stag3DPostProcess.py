#!/Library/Frameworks/Python.framework/Versions/3.10/bin/python3
import numpy as np
import sys
import matplotlib.pyplot as plt

import os


def vProfReader(vProfFilename):
    file1 = open(vProfFilename, 'r')
    Lines = file1.readlines()
    nLines = len(Lines)
    Zlocs = []
    Pprof = []
    Uprof = []
    Vprof = []
    Wprof = []
    Tprof = []
    for line in Lines:
        items = line.split(' ')
        #print('Nitems=',len(items))
        Zlocs.append(items[0])
        Pprof.append(items[1])
        Uprof.append(items[2])
        Vprof.append(items[3])
        Wprof.append(items[4])
        Tprof.append(items[5])

    fig = plt.figure(figsize=(18, 10))
    ax2 = fig.add_subplot(131)
    plt.plot(Zlocs,Pprof)
    ax2.set_xlabel('Zlocation (m)')
    ax2.set_ylabel('Pressure (Pa)')

    ax2 = fig.add_subplot(132)
    line0 = plt.plot(Zlocs,Uprof,label='U')
    line0 = plt.plot(Zlocs,Vprof,label='V')
    line0 = plt.plot(Zlocs,Wprof,label='W')
    ax2.set_xlabel('Zlocation (m)')
    ax2.set_ylabel('Vel (m/s)')
    plt.legend()

    ax2 = fig.add_subplot(133)
    plt.plot(Zlocs,Tprof)
    ax2.set_xlabel('Zlocation (m)')
    ax2.set_ylabel('T (K)')



    plt.show()
    #print('vectorSize: ',nVectorDataPts)
    #print('photoNet: ',photoNet)
    return Zlocs, Pprof, Uprof, Vprof, Wprof, Tprof


def sliceReader(sliceFilename):
    file1 = open(sliceFilename, 'r')
    Lines = file1.readlines()
    nLines = len(Lines)
    NX, NZ = Lines[0].split(' ')
    NX, NZ = int(NX), int(NZ)
    print('NX: ',NX,' NZ: ',NZ)
    meshgridX = np.zeros([NX,NZ])
    meshgridZ = np.zeros([NX,NZ])
    meshgridP = np.zeros([NX,NZ])
    meshgridU = np.zeros([NX,NZ])
    meshgridV = np.zeros([NX,NZ])
    meshgridW = np.zeros([NX,NZ])
    meshgridT = np.zeros([NX,NZ])
    meshgridHx = np.zeros([NX,NZ])
    meshgridHy = np.zeros([NX,NZ])
    meshgridHz = np.zeros([NX,NZ])
    meshgridPprime = np.zeros([NX,NZ])
    meshgridMassbal = np.zeros([NX,NZ])
    idx = 1;
    for i in range(0,NX):
        for k in range(0,NZ):
            line = Lines[idx]
            items = line.split(' ')
            #print('nItems: ',len(items))
            meshgridX[i,k] = items[0]
            meshgridZ[i,k] = items[1]
            meshgridP[i,k] = items[2]
            meshgridU[i,k] = items[3]
            meshgridV[i,k] = items[4]
            meshgridW[i,k] = items[5]
            meshgridT[i,k] = items[6]
            meshgridHx[i,k] = items[7]
            meshgridHy[i,k] = items[8]
            meshgridHz[i,k] = items[9]
            meshgridPprime[i,k] = items[10]
            meshgridMassbal[i,k] = items[11]
            idx +=1


    fig = plt.figure(figsize=(16, 11))
    ax2 = fig.add_subplot(331)
    plt.contourf(meshgridX, meshgridZ, meshgridU)
    #ax2.set_xlabel('Xlocation (m)')
    ax2.set_ylabel('Zlocation (m)')
    plt.title('U-Vel')
    plt.colorbar()
    #plt.cmap()

    ax2 = fig.add_subplot(332)
    ax2.contourf(meshgridX, meshgridZ, meshgridV)
    #ax2.set_xlabel('Xlocation (m)')
    #ax2.set_ylabel('Zlocation (m)')
    plt.title('V-Vel')
    plt.colorbar()

    ax2 = fig.add_subplot(333)
    plt.contourf(meshgridX, meshgridZ, meshgridW)
    #ax2.set_xlabel('Xlocation (m)')
    #ax2.set_ylabel('Zlocation (m)')
    plt.title('W-Vel')
    plt.colorbar()


    ax2 = fig.add_subplot(334)
    plt.contourf(meshgridX, meshgridZ, meshgridHx)
    #ax2.set_xlabel('Xlocation (m)')
    ax2.set_ylabel('Zlocation (m)')
    plt.title('HX')
    plt.colorbar()


    ax2 = fig.add_subplot(335)
    ax2.contourf(meshgridX, meshgridZ, meshgridHy)
    #ax2.set_xlabel('Xlocation (m)')
    #ax2.set_ylabel('Zlocation (m)')
    plt.title('HY')
    plt.colorbar()

    ax2 = fig.add_subplot(336)
    plt.contourf(meshgridX, meshgridZ, meshgridHz)
    #ax2.set_xlabel('Xlocation (m)')
    #ax2.set_ylabel('Zlocation (m)')
    plt.title('HZ')
    plt.colorbar()

    ax2 = fig.add_subplot(337)
    plt.contourf(meshgridX, meshgridZ, meshgridP)
    ax2.set_xlabel('Xlocation (m)')
    ax2.set_ylabel('Zlocation (m)')
    plt.title('Pressure')
    plt.colorbar()
    #plt.cmap()

    ax2 = fig.add_subplot(338)
    ax2.contourf(meshgridX, meshgridZ, meshgridPprime)
    ax2.set_xlabel('Xlocation (m)')
    #ax2.set_ylabel('Zlocation (m)')
    plt.title('Pprime')
    plt.colorbar()

    ax2 = fig.add_subplot(339)
    plt.contourf(meshgridX, meshgridZ, meshgridMassbal)
    ax2.set_xlabel('Xlocation (m)')
    #ax2.set_ylabel('Zlocation (m)')
    plt.title('Massbal')
    plt.colorbar()

    plt.show()
    #print('vectorSize: ',nVectorDataPts)
    #print('photoNet: ',photoNet)
    return


if __name__ == '__main__':
    command1 = 'g++ -g main.cpp'
    command2 = './a.out'
    #print(command1)
    os.system(command1)
    os.system(command2)


    sliceReader('vSlice.txt')
    Zlocs, Pprof, Uprof, Vprof, Wprof, Tprof = vProfReader('VerticalProfile.txt')
