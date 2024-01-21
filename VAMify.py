#!/Library/Frameworks/Python.framework/Versions/3.10/bin/python3
import matplotlib.pyplot as plt
from   matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import sys
from parseGPX import *
from plot_PMdata import *
from myRegress import myRegress
from krigingInterpolation import krigingInterpolation
import numpy as np
plt.style.use('bmh')

def VAMify(sPath, elev, power, cad, hr,gpxfilename):


    #print(power)
    plot_PMdata(power,gpxfilename)

    nFilt = 30;
    nFilt5 = 150;
    nFilt10 = 300;
    nPts=len(elev)
    VAM = np.zeros((nPts,1),float)
    VAM5 = np.zeros((nPts,1),float)
    VAM10 = np.zeros((nPts,1),float)
    VAM1000ref =1000*np.ones((nPts,1),float)
    maxVAM =0
    maxVAM5 =0
    maxVAM10 =0
    for i in range(1+nFilt,nPts-nFilt):
        VAM[i]=(3600/(2*nFilt+1))*max(0,(elev[i+nFilt]-elev[i-nFilt]))
        maxVAM = max(maxVAM,VAM[i])

    for i in range(1+nFilt5,nPts-nFilt5):
        VAM5[i]=(3600/(2*nFilt5+1))*max(0,(elev[i+nFilt5]-elev[i-nFilt5]))
        maxVAM5 = max(maxVAM5,VAM5[i])

    for i in range(1+nFilt10,nPts-nFilt10):
        VAM10[i]=(3600/(2*nFilt10+1))*max(0,(elev[i+nFilt10]-elev[i-nFilt10]))
        maxVAM10 = max(maxVAM10,VAM10[i])
    print(maxVAM,maxVAM5,maxVAM10)


    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax1.set_ylabel('Elevation (m)')
    ax1.set_xlabel('Distance (km)')
    line1 = plt.plot(sPath/1000,elev,'b', linewidth=1)
    ax1.xaxis.set_major_locator(MultipleLocator(5))
    #ax1.set_xlim([0,105])

    ax2 = fig.add_subplot(212)
    ax2.set_ylabel('VAM (m/hr)')
    ax2.set_xlabel('Distance (km)')
    line1 = plt.plot(sPath/1000,VAM,'r', label='1-Min',linewidth=1)
    line2 = plt.plot(sPath/1000,VAM5,'g', label='5-Min',linewidth=1)
    line3 = plt.plot(sPath/1000,VAM10, 'k',label='10-Min',linewidth=1)
    line4 = plt.plot(sPath/1000,VAM1000ref,'b',linewidth=1)
    plt.legend()

    VAMbins = [120,180,240,300,360,420,480,540,600,720,840,960]
    VAMtable = np.empty((0, 3), float)
    powerIntegral = np.zeros(nPts)
    for i in range(1,nPts):
        powerIntegral[i]=powerIntegral[i-1]+power[i]

    for bin in VAMbins:
        if bin < nPts:
            maxVAMi=0
            maxPOWi=0
            for i in range(0,nPts):
                iMin = max(0,i-int(bin/2))
                iMax = min(nPts-1,i+int(bin/2))
                iVAM = (elev[iMax]-elev[iMin])*3600/bin
                if iVAM > maxVAMi:
                    maxVAMi = iVAM
                    maxPOWi = (powerIntegral[iMax]-powerIntegral[iMin])/bin
            VAMtable= np.append(VAMtable,[[bin, maxPOWi, maxVAMi]],axis=0)

    powerLocs = [[160],[200],[240],[260],[280],[290],[305],[320],[335]]
    VAMinterp, variogram = krigingInterpolation(VAMtable[:,1],VAMtable[:,2],powerLocs)
    print('XtoRegFunc',VAMtable[:,1])
    print('YtoRegFunc',VAMtable[:,2])
    coeffs, yhat, resid = myRegress(VAMtable[:,2],VAMtable[:,1])
    m1 =coeffs[0]+1e-6
    b1 =coeffs[1]
    #m1=1
    #b1=1
    #yhat = 0*VAMtable[:,1]




    coeffs2, yhat2, resid2 = myRegress(VAMtable[:,1],VAMtable[:,2])
    m2 =coeffs2[0]
    b2 =coeffs2[1]
    p1000vam = np.array( [(1000-b1)/m1 , 1000*m2+b2])
    p1000vamstar = p1000vam.sum()/(1e-6+len(p1000vam))
    vam1000 = np.array([1000, 1000])




    fig2 = plt.figure(2)
    ax2 = fig2.add_subplot(131)
    line1 = plt.plot(VAMtable[:,0],VAMtable[:,1],label='Power Curve')
    ax2a = ax2.twinx()
    line2 = plt.plot(VAMtable[:,0],VAMtable[:,2],'k',label='VAM Curve')
    ax2.set_ylabel('Power (W)')
    ax2.set_xlabel('Duration (s)')
    ax2a.set_ylabel('VAM')
    plt.legend()


    ax2 = fig2.add_subplot(132)
    line1 = plt.scatter(VAMtable[:,1],VAMtable[:,2],label='vam vs power')
    line2 = plt.plot(VAMtable[:,1],yhat,label='fit')
    line3 = plt.scatter(p1000vam,vam1000,label='P for 1000VAM ='+str(int(p1000vamstar)))
    line4 = plt.plot(powerLocs,VAMinterp,label='Kriging')
    plt.legend()
    ax2.set_xlabel('Power (W)')

    ax3 = fig2.add_subplot(133)
    line5 = plt.scatter(variogram[:,0],variogram[:,1],label='Variogram')
    ax3.set_xlabel('lag (W)')
    ax3.set_ylabel('Abs Diff (VAM)')

    plt.show(block=False)
    plt.show()
    #print(type(sPath))










if __name__ == '__main__':
    n = len(sys.argv)
    #print("Total arguments passed:", n)
    #for i in range(1, n):
        #print(i," ",sys.argv[i])

    if n < 2:
        print("Correct Usage VAMify.py filename.gpx")
    inputfilename = sys.argv[1];
    sPath, elev, power, cad, hr = parseGPX(inputfilename)
    VAMify(sPath, elev, power, cad, hr, inputfilename)
