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

def granFondoOptExplicit(gpxFilename, basePower, myMass, bikeMass, CdA):

    sPath, zElev, oldPower, cad, hr = parseGPX(gpxFilename)
    mass = myMass+bikeMass;
    eta = 0.5*1.18*CdA;
    driveEff = 0.95
    dt =1;
    slopeFilterDistance = 30;
    nPts=len(zElev);
    maxDist = sPath[nPts-1]
    HRmin = 73;
    HRmax = 187;
    Pmax = 300;
    NFrontMin = 36;
    NFrontMax = 52;
    NRearMin = 11;
    NRearMax = 30;
    MinCadence = 60; #Minimum Cadence at Minimum Gear Ratio
    MaxCadence = 85;
    MinSpeedKPH = 3.6*(1/60)*2.14*MinCadence*NFrontMin/NRearMax;
    MaxSpeedKPH = 3.6*(1/60)*2.14*MaxCadence*NFrontMax/NRearMin;

    #Filter Profile
    print('Begin Filtering')
    zFiltered = zElev;
    slopesPct = 0.0*zElev
    powerFax = 0.0*zElev
    aeroFax = 0.0*zElev
    i0=0
    for i in range(0,nPts):
        sMin = sPath[i]-0.5*slopeFilterDistance
        sMax = sPath[i]+0.5*slopeFilterDistance
        zPlus, itmp = getQtyAtSfromI0(sMax,sPath,zElev,i0)
        zMinus, i0 = getQtyAtSfromI0(sMin,sPath,zElev,i0)
        zFiltered[i]=0.5*(zMinus+zPlus)
    i0=0
    for i in range(0,nPts):
        sMin = sPath[i]-0.5*slopeFilterDistance
        sMax = sPath[i]+0.5*slopeFilterDistance
        zPlusFilt, itmp = getQtyAtSfromI0(sMax,sPath,zFiltered,i0)
        zMinusFilt, i0 = getQtyAtSfromI0(sMin,sPath,zFiltered,i0)
        slopesPct[i]=100*(zPlusFilt-zMinusFilt)/slopeFilterDistance
        powerFax[i],aeroFax[i] = powerSlopeFactor(slopesPct[i])
    print('Filtering Complete, Max Slope: ',max(slopesPct))

    fig = plt.figure()
    ax1 = fig.add_subplot(611)
    line0 = plt.plot(sPath/1000,zElev,label='Raw')
    line1 = plt.plot(sPath/1000,zFiltered,label='Filtered')
    ax2 = fig.add_subplot(612)
    line0 = plt.plot(sPath/1000,slopesPct,label='Raw')
    ax3 = fig.add_subplot(613)
    line0 = plt.plot(sPath/1000,powerFax,label='Raw')
    line0 = plt.plot(sPath/1000,aeroFax,label='Raw')
    #plt.show()



    iSpeed = 0.5
    iAcc = 0
    iDist = 0
    iStep=0
    iPower = 0
    iHR = HRmin
    i0=0
    HRX2 = 0;
    speeds =[]
    sPos =  []
    zPos =  []
    HRs = []
    powers = []
    VAM = []
    while iDist < maxDist:
        iSlope, i0 = getQtyAtSfromI0(iDist,sPath,slopesPct, i0);
        iZ, i0 = getQtyAtSfromI0(iDist,sPath,zFiltered, i0);
        powerFactor, aeroFactor = powerSlopeFactor(iSlope);
        iPower= basePower*powerFactor;
        iForce  = driveEff*iPower/iSpeed
        netForce  = iForce - mass*9.81*iSlope/100 - aeroFactor*eta*iSpeed*iSpeed
        iAcc = netForce/mass
        HRdot, HRX2dot = dHdt(iHR, HRmin, HRmax, iPower, Pmax,HRX2)
        #rint(iStep,iSlope,iPower,3.6*iSpeed,iAcc)
        HRs = np.append(HRs,iHR)
        sPos = np.append(sPos,iDist)
        speeds = np.append(speeds,iSpeed)
        zPos = np.append(zPos,iZ)
        powers = np.append(powers,iPower)
        VAM = np.append(VAM,max(0,iSpeed*iSlope*3600/100))

        iDist = iDist + dt*iSpeed
        iSpeed = max(iSpeed +iAcc*dt,0.1)
        iHR = iHR + HRdot
        HRX2 = HRX2 + HRX2dot
        iStep = iStep+1
    iStep = iStep-1
    hours = int(iStep/3600)
    minutes = int((iStep-hours*3600)/60)
    seconds = int(iStep-hours*3600-minutes*60)

    ax4 = fig.add_subplot(614)
    line0 = plt.plot(sPos/1000,powers,label='Power(W)')
    ax4 = fig.add_subplot(615)
    line0 = plt.plot(sPos/1000,speeds*3.6,label='KPH')
    ax4 = fig.add_subplot(616)
    line0 = plt.plot(sPos/1000,VAM,label='HR')
    plt.legend()
    plt.show()

    VAMify(sPos, zPos, powers, 0*powers, HRs,gpxFilename)
    print('Elapsed Time: ',iStep,' ',hours,':',minutes,':',seconds)
    compare2GPXroutes(sPos, zPos, powers, 0*powers, HRs ,'Sim',sPath, zElev, oldPower, cad, hr,'Real')


def powerSlopeFactor(slope):
    #lookupTable is Nx3: column0=slope in pct column1=powerFactor column2=PPaeroFactor
    slopeLookupTable = [ [-4 , 0, 1],[-2, .3, 0.9],[0, 1.0, 0.9],[1 , 1.05, 0.9], [2,  1.1, 1.0],[3,  1.15, 1.0],[4 ,  1.1, 1.0],[6 , 1.2, 1.0],[8 , 1.2 ,1.0], [10 ,1.3 ,1.0]]
    nn = len(slopeLookupTable)
    #print(nn,len(slopeLookupTable[0]))
    if slope < slopeLookupTable[0][0]:
        return slopeLookupTable[0][1], slopeLookupTable[0][2]
    if slope > slopeLookupTable[nn-1][0]:
        return slopeLookupTable[nn-1][1], slopeLookupTable[nn-1][2]
    for i in range(0,nn-1):
        if slope >= slopeLookupTable[i][0]:
            if slope < slopeLookupTable[i+1][0]:
                x1 = slopeLookupTable[i][0]
                x2 = slopeLookupTable[i+1][0]
                y1 = slopeLookupTable[i][1]
                y2 = slopeLookupTable[i+1][1]
                z1 = slopeLookupTable[i][2]
                z2 = slopeLookupTable[i+1][2]
                powerFactor = y1+(slope-x1)*(y2-y1)/(x2-x1)
                aeroFactor  = z1+(slope-x1)*(z2-z1)/(x2-x1)
                return powerFactor, aeroFactor
    return 0.1,10

def dHdt(H,HRmin, HRmax, power, Pmax,HRX2):
    a1 = .0113
    a2 = .0072
    a3 = .0041
    a4 = .0049
    a5 = 19.8002
    X1 = (H-HRmin)/HRmax
    u = power/Pmax;
    phi = a4*X1/(1+math.exp(-1*(X1-a5)));

    X1dot = -a1*X1+a2*HRX2+a2*u;
    HRX2dot = -a3*HRX2+phi;
    HRdot = X1dot*HRmax;
    return HRdot, HRX2dot



def getQtyAtSfromI0(dist,sPath,Qty,i0):
    nPoints = len(sPath);
    if dist > sPath[nPoints-1]:
        return Qty[nPoints-1], nPoints-1
    if dist < sPath[0]:
        return Qty[0], 0
    elev = Qty[0];
    for i in range(i0,nPoints-1):
        if dist >= sPath[i]:
            elev = Qty[i]
            if dist < sPath[i+1]:
                x1 = sPath[i]
                x2 = sPath[i+1]
                y1 = Qty[i]
                y2 = Qty[i+1]
                elev = y1+(dist-x1)*(y2-y1)/(x2-x1)
                return elev, i
    return elev, 0




if __name__ == '__main__':
    n = len(sys.argv)
    print("Total arguments passed:", n)
    for i in range(1, n):
        print(i," ",sys.argv[i])

    if n < 2:
        print("Correct Usage granFondoOptExplicit.py filename.gpx")

    basePower = 220
    myMass = 86
    bikeMass = 7.5
    CdA = 0.45
    inputfilename = sys.argv[1]

    if n ==2:
        granFondoOptExplicit(inputfilename, basePower, myMass, bikeMass, CdA)
