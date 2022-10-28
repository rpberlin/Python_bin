#!/Library/Frameworks/Python.framework/Versions/3.10/bin/python3

import sys
import math
from parseGPX import parseGPX
from plot_PMdata import calcPMdata
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.ticker import StrMethodFormatter, NullFormatter
import numpy as np
import scipy
plt.style.use('bmh')

def compare2GPXroutes(sPath1, elev1, power1, cad1, hr1 ,gpx1,sPath2, elev2, power2, cad2, hr2,gpx2):
    #sPath1, elev1, power1, cad1, hr1 = parseGPX(gpx1)
    #sPath2, elev2, power2, cad2, hr2 = parseGPX(gpx2)

    speeds1 = np.diff(sPath1)*3.6
    speeds2 = np.diff(sPath2)*3.6
    speeds1 = np.append(speeds1[1], speeds1)
    speeds2 = np.append(speeds2[0], speeds2)

    nPts1 = len(sPath1)
    nPts2 = len(sPath2)
    powerIntegral1 = np.zeros(nPts1)
    powerIntegral2 = np.zeros(nPts2)
    if len(power1) == 0:
        power1=np.zeros(nPts1)
    if len(power2) == 0:
        power2=np.zeros(nPts2)
    #print(nPts1," ",nPts2, " sizes of power1 and power2 ",len(elev1)," ",len(elev2))
    for i in range(1,nPts1):
        a = power1[i]
        powerIntegral1[i]=powerIntegral1[i-1]+a
    for i in range(1,nPts2):
        aa = power2[i]
        powerIntegral2[i]=powerIntegral2[i-1]+aa

    pace1 = 0*speeds1;
    pace2 = 0*speeds2;
    nFilt =30;
    for i in range(0,nPts1):
        i1=max(0,i-nFilt)
        i2=min(nPts1-1,i+nFilt)
        speeds1[i]=3.6*(sPath1[i2]-sPath1[i1])/(i2-i1)
        pace1[i]=60/max(speeds1[i],4)

    for i in range(0,nPts2):
        i1=max(1,i-nFilt);
        i2=min(nPts2-1,i+nFilt);
        speeds2[i]=3.6*(sPath2[i2]-sPath2[i1])/(i2-i1);
        pace2[i]=60/max(speeds2[i],4);


    print("leng(power1)= ",len(power1))
    POWfilt1=power1;
    POWfilt2=power2;
    VAM1=np.zeros(nPts1);
    VAM2=np.zeros(nPts2);
    nFilt=30;
    for i in range(0,nPts1):
        i1=max(0,i-nFilt);
        i2=min(nPts1-1,i+nFilt);
        VAM1[i]=(3600/(i2-i1))*max(0,(elev1[i2]-elev1[i1]));
        POWfilt1[i]=(powerIntegral1[i2]-powerIntegral1[i1])/(i2-i1);

    for i in range(0,nPts2):
        i1=max(0,i-nFilt);
        i2=min(nPts2-1,i+nFilt);
        VAM2[i]=(3600/(i2-i1))*max(0,(elev2[i2]-elev2[i1]));
        POWfilt2[i]=(powerIntegral2[i2]-powerIntegral2[i1])/(i2-i1);


    fig = plt.figure()
    ax1 = fig.add_subplot(611)
#    line1 = plt.plot(sPath1/1000,speeds1,label=gpx1)
#    line2 = plt.plot(sPath2/1000,speeds2,label=gpx2)
#    ax1.set_ylabel('kPh');
    line1 = plt.plot(sPath1/1000,pace1,label=gpx1)
    line2 = plt.plot(sPath2/1000,pace2,label=gpx2)
    ax1.set_ylabel('Minutes/Km');
    plt.legend()


    ax2 = fig.add_subplot(612)
    tDiff=np.zeros(nPts2);
    i0=0
    for i in range(0,nPts2):
        tinterp, i0=getQtyAtSfromI0(sPath2[i],sPath1,range(0,nPts1),i0)
        tDiff[i] = tinterp-i
        #tDiff[i]=getQtyAtS(sPath2[i],sPath1,range(0,nPts1))-i;

    plt.plot(sPath2/1000,np.zeros(nPts2))
    plt.plot(sPath2/1000,tDiff/60,label='Time Gap (min)');
    ax2.set_xlabel('km');
    ax2.set_ylabel('Time Gap (min)');
    plt.legend()



    ax3 = fig.add_subplot(613)
    plt.plot(sPath1/1000,elev1)
    ax3.set_ylabel('Elev (m)');


    ax4 = fig.add_subplot(614)
    plt.plot(sPath1/1000,VAM1)
    plt.plot(sPath2/1000,VAM2)
    ax4.set_ylabel('VAM (m/hr)');

    ax5 = fig.add_subplot(615)
    line1  = plt.plot(sPath1/1000,POWfilt1)
    line2  = plt.plot(sPath2/1000,POWfilt2)
    ax5.set_ylabel('1 Min. Power (W)');


    ax6 = fig.add_subplot(616)
    plt.plot(sPath1/1000,hr1,label=gpx1)
    plt.plot(sPath2/1000,hr2,label=gpx2)
    plt.legend()
    ax6.set_ylabel('HR (bpm)');


    powerProfiles1, powerMinutes1 = calcPMdata(power1)
    powerProfiles2, powerMinutes2 = calcPMdata(power2)

    fig2 = plt.figure(2)
    ax1a = fig2.add_subplot(131)
    plt.semilogx(powerProfiles1[:,0]/60,powerProfiles1[:,1],label=gpx1)
    plt.semilogx(powerProfiles2[:,0]/60,powerProfiles2[:,1],label=gpx2)
    plt.legend()
    ax1a.xaxis.set_major_formatter(mticker.ScalarFormatter())
    ax1a.set_xlabel('Duration (minutes)');
    ax1a.set_ylabel('Power (W)');

    ax2a = fig2.add_subplot(132)
    plt.semilogx(powerProfiles1[:,0]/60,powerProfiles1[:,2],label=gpx1)
    plt.semilogx(powerProfiles2[:,0]/60,powerProfiles2[:,2],label=gpx2)
    plt.legend()
    ax2a.xaxis.set_major_formatter(mticker.ScalarFormatter())
    ax2a.set_xlabel('Duration (minutes)');
    ax2a.set_ylabel('Normalized Power (W)');

    ax3a = fig2.add_subplot(133)
    plt.semilogx(powerMinutes1[:,0],powerMinutes1[:,1],label=gpx1)
    plt.semilogx(powerMinutes2[:,0],powerMinutes2[:,1],label=gpx2)
    plt.xticks(np.arange(200,300,step=20))
    plt.legend()
    #plt.grid(visible=True, which='minor', linestyle='--')
    ax3a.xaxis.set_minor_formatter(mticker.ScalarFormatter())
    #ax3a.xaxis.set_major_formatter(mticker.ScalarFormatter())
    ax3a.xaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
    ax3a.set_xlabel('Power (W)');
    ax3a.set_ylabel('Duration (minutes)');
    plt.show()

    Pbar1=sum(power1)/nPts1
    Pbar2=sum(power2)/nPts2
    WPbar1 = math.pow(sum([num ** 4 for num in power1])/nPts1,0.25)
    WPbar2 = math.pow(sum([num ** 4 for num in power2])/nPts2,0.25)
    Ubar1=3.6*sPath1[nPts1-1]/nPts1
    Ubar2=3.6*sPath2[nPts2-1]/nPts2
    DwgtUbar1 = 3.6*sum([num ** 2 for num in speeds1])/[nPts1-1]
    DwgtUbar2 = 3.6*sum([num ** 2 for num in speeds2])/[nPts2-1]




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

def getQtyAtS(dist,sPath,Qty):
    nPoints = len(sPath);
    if dist > sPath[nPoints-1]:
        return Qty[nPoints-1]
    if dist < sPath[0]:
        return Qty[0]
    elev = Qty[0];
    for i in range(0,nPoints-1):
        if dist >= sPath[i]:
            elev = Qty[i]
            if dist < sPath[i+1]:
                x1 = sPath[i]
                x2 = sPath[i+1]
                y1 = Qty[i]
                y2 = Qty[i+1]
                elev = y1+(dist-x1)*(y2-y1)/(x2-x1)
                return elev
    return elev

if __name__ == '__main__':
    n = len(sys.argv)
    print("Total arguments passed:", n)
    for i in range(1, n):
        print(i," ",sys.argv[i])

    if n < 3:
        print("Correct Usage compare2GPXroutes.py filename1.gpx filename2.gpx")

    if n >= 3:
        inputfilename1 = sys.argv[1];
        inputfilename2 = sys.argv[2];
        sPath1, elev1, power1, cad1, hr1 = parseGPX(inputfilename1)
        sPath2, elev2, power2, cad2, hr2 = parseGPX(inputfilename2)
        compare2GPXroutes(sPath1, elev1, power1, cad1, hr1 ,inputfilename1,sPath2, elev2, power2, cad2, hr2,inputfilename2)
