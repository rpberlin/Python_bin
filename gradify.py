#!/Library/Frameworks/Python.framework/Versions/3.10/bin/python3
import matplotlib.pyplot as plt
from   matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import sys
from parseGPX import *
from plot_PMdata import *
from granFondoOptExplicit import getQtyAtSfromI0
from myRegress import myRegress
from krigingInterpolation import krigingInterpolation
import numpy as np
plt.style.use('bmh')

def gradify(sPath, zElev, power, cad, hr,gpxfilename):
    #Filter Profile
    print('Begin Filtering')
    nPts = len(zElev)
    slopeFilterDistance = 50
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
    print('Filtering Complete, Max Slope: ',max(slopesPct))

    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    #line0 = plt.plot(sPath/1000,zElev,label='Raw')
    line1 = plt.plot(sPath/1000,zFiltered,color='b',label='Elevation')
    ax1.set_ylabel('Elevation (m)')

    ax2 = fig.add_subplot(212)
    #line0 = plt.plot(sPath/1000,zElev,label='Raw')
    line2 = plt.plot(sPath/1000,slopesPct,color='r',label='Slope')
    ax1.set_ylabel('Slope (%)')

    bins = np.linspace(-10,15,50)
    counts,bins = np.histogram(slopesPct,50)
    fig2 = plt.figure(2)
    plt.stairs(counts,bins)
    plt.title(f'Mean = {np.mean(slopesPct):.3f}  Stdev: {np.std(slopesPct):.3f} (kmh_rms)')

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
    sPath, elev, power, cad, hr, times = parseGPX(inputfilename)
    gradify(sPath, elev, power, cad, hr, inputfilename)
