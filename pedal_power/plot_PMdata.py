#!/Library/Frameworks/Python.framework/Versions/3.10/bin/python3
import scipy
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import sys
import math
import numpy as np
from parseGPX import parseGPX
plt.style.use('bmh')

def plot_PMdata(power,gpxfilename):
    pmData = scipy.io.loadmat('/Users/ryanblanchard/myApplications/Matlab_bin/powerMinutes.mat')
    pmApple21tx=pmData['PMtxApple']
    pmAOTCtx=pmData['PMtxAOTC']
    pmLKNtx=pmData['PMtxLKN']
    pmJan25tx=pmData['PMtxJan25']
    pmTourmtx=pmData['PMtxTourm']
    pmGoodtx=pmData['PMtxgood']
    pmLLanbtx=pmData['PMtxLlanb']

    pmApple21=pmData['pMApple']
    pmAOTC=pmData['pMAOTC']
    pmLKN=pmData['pMLKN']
    pmJan25=pmData['pMJan25']
    pmTourm=pmData['pMTourm']
    pmGood=pmData['pMgood']
    pmLLanb=pmData['pMLlanb']

    pmDataCurrX, pmDataCurr = calcPMdata(power)


    #print(pmDataCurrX)
    fig = plt.figure()
    ax1 = fig.add_subplot(131)
    ax1.set_ylabel('Average Power (W)')
    ax1.set_xlabel('Duration (min)')
    line1 = plt.semilogx(pmApple21tx[:,0]/60,pmApple21tx[:,1], label='Apple21', linewidth=1)
    line2 = plt.semilogx(pmAOTCtx[:,0]/60,pmAOTCtx[:,1], label='AOTC', linewidth=1)
    line3 = plt.semilogx(pmLKNtx[:,0]/60,pmLKNtx[:,1], label='LKN', linewidth=1)
    line4 = plt.semilogx(pmJan25tx[:,0]/60,pmJan25tx[:,1], label='Jan25', linewidth=1)
    line5 = plt.semilogx(pmTourmtx[:,0]/60,pmTourmtx[:,1], label='Tourm', linewidth=1)
    line6 = plt.semilogx(pmGoodtx[:,0]/60,pmGoodtx[:,1], label='Good', linewidth=1)
    line7 = plt.semilogx(pmLLanbtx[:,0]/60,pmLLanbtx[:,1], label='LLanb', linewidth=1)
    linecurr = plt.semilogx(pmDataCurrX[:,0]/60,pmDataCurrX[:,1], 'k',label=gpxfilename, linewidth=2)
    ax1.xaxis.set_major_formatter(mticker.ScalarFormatter())
    plt.legend()

    ax2 = fig.add_subplot(132)
    ax2.set_ylabel('Normalized Average Power (W)')
    ax2.set_xlabel('Duration (min)')
    line1 = plt.semilogx(pmApple21tx[:,0]/60,pmApple21tx[:,2], label='Apple21', linewidth=1)
    line2 = plt.semilogx(pmAOTCtx[:,0]/60,pmAOTCtx[:,2], label='AOTC', linewidth=1)
    line3 = plt.semilogx(pmLKNtx[:,0]/60,pmLKNtx[:,2], label='LKN', linewidth=1)
    line4 = plt.semilogx(pmJan25tx[:,0]/60,pmJan25tx[:,2], label='Jan25', linewidth=1)
    line5 = plt.semilogx(pmTourmtx[:,0]/60,pmTourmtx[:,2], label='Tourm', linewidth=1)
    line6 = plt.semilogx(pmGoodtx[:,0]/60,pmGoodtx[:,2], label='Good', linewidth=1)
    line6 = plt.semilogx(pmLLanbtx[:,0]/60,pmLLanbtx[:,2], label='LLanb', linewidth=1)
    linecurr = plt.semilogx(pmDataCurrX[:,0]/60,pmDataCurrX[:,2], 'k', label=gpxfilename, linewidth=2)
    ax2.xaxis.set_major_formatter(mticker.ScalarFormatter())

    plt.legend()


    ax3 = fig.add_subplot(133)
    ax3.set_ylabel('Duration (min)')
    ax3.set_xlabel('Power (W)')
    line1 = plt.semilogx(pmApple21[:,0],pmApple21[:,1], label='Apple21', linewidth=1)
    line2 = plt.semilogx(pmAOTC[:,0],pmAOTC[:,1], label='AOTC', linewidth=1)
    line3 = plt.semilogx(pmLKN[:,0],pmLKN[:,1], label='LKN', linewidth=1)
    line4 = plt.semilogx(pmJan25[:,0],pmJan25[:,1], label='Jan25', linewidth=1)
    line5 = plt.semilogx(pmTourm[:,0],pmTourm[:,1], label='Tourm', linewidth=1)
    line6 = plt.semilogx(pmGood[:,0],pmGood[:,1], label='Good', linewidth=1)
    line7 = plt.semilogx(pmLLanb[:,0],pmLLanb[:,1], label='Llanb', linewidth=1)
    linecurr = plt.semilogx(pmDataCurr[:,0],pmDataCurr[:,1], 'k', label=gpxfilename, linewidth=2)
    ax3.xaxis.set_minor_formatter(mticker.ScalarFormatter())
    plt.legend()
    plt.show(block=False)
    plt.show()

def calcPMdata(power):
    #Calculate Linear and Normalized Power Profiles
    powerTimeBins=[30,45,60,90,120,180,240,300,360,420,480,540,600,720,840,960,1080,1200,25*60,30*60,35*60,40*60,45*60,50*60,55*60,60*60,1.5*60*60,2*60*60]
    n = len(power)
    if n==0:
        power = np.zeros(2)
        n=2

    powerProfiles = np.empty((0, 3), float)
    powerIntegrals = np.zeros((n, 2))
    for i in range(1,n):
        powerIntegrals[i,0]=powerIntegrals[i-1,0]+power[i]
        powerIntegrals[i,1]=powerIntegrals[i-1,1]+power[i]**4
        #print(i,power[i],powerIntegrals[i,0],powerIntegrals[i,1])
    idx=0
    for bin in powerTimeBins:
        if n>bin:
            filt = int(bin/2)
            bestPow = 0
            bestPowNP = 0
            for i in range(filt,n-filt):
                imax = i+filt
                imin = i-filt
                sumPower = powerIntegrals[imax,0]-powerIntegrals[imin,0]
                sumPower4 =powerIntegrals[imax,1]-powerIntegrals[imin,1]
                bestPow = max(bestPow, sumPower/bin)
                bestPowNP = max(bestPowNP, math.pow(sumPower4/bin,0.25) )
            powerProfiles= np.append(powerProfiles,[[bin, bestPow, bestPowNP]],axis=0)
            #print(bin,bestPow, bestPowNP)


    #Calculate Total Cumulative Power-Time Profile
    maxP = int(max(power))
    powerMinutes = np.empty((0, 2), float)
    for i in range(200,maxP):
        count = (1/60)*sum(j > i for j in power)
        powerMinutes= np.append(powerMinutes,[[i, count]],axis=0)

    return powerProfiles, powerMinutes



if __name__ == '__main__':
    n = len(sys.argv)
    print("Total arguments passed:", n)
    for i in range(1, n):
        print(i," ",sys.argv[i])

    if n < 2:
        print("Correct Usage plot_PMdata.py filename.gpx")
    inputfilename = sys.argv[1];
    sPath, elev, power, cad, hr, time = parseGPX(inputfilename)
    plot_PMdata(power,inputfilename)
