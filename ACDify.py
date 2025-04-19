#!/Library/Frameworks/Python.framework/Versions/3.10/bin/python3


import mySSACD as SSACD
import sys, pathlib
import matplotlib.pyplot as plt
import parseGPX
import pickle
HRobjDirectory = '/Users/ryanblanchard/Documents/gpxFiles/HRobjs/'

if __name__ == '__main__':
    n = len(sys.argv)
    mass_bike = 8
    mass_rider = 85

    print("Total arguments passed:", n)
    #for i in range(1, n):
        #print(i," ",sys.argv[i])
    if n == 1:
        inputfilename = '/Users/ryanblanchard/Downloads/tmp2.tcx'
        #inputfilename = '/Users/ryanblanchard/Documents/gpxFiles/Ardiden309.gpx'
        #inputfilename = '/Users/ryanblanchard/Documents/gpxFiles/Ardiden265.gpx'
    if n==2:
        inputfilename = sys.argv[1]

    sPath, elev, power, cad, hr, times = parseGPX.parseGPX(inputfilename)

    ACD0 = 0.4
    myACDobj0 = SSACD.mySSACD(ACD0,sPath,elev,power,mass_rider,mass_bike)
    plt.ion()
    fig = plt.figure()
    ax1 = fig.add_subplot(411)
    line1 = plt.plot(elev,color='b',label='Elevation')
    ax1.set_ylabel('Elevation (m)')

    ax2 = fig.add_subplot(412)
    line2 = plt.plot(cad,color='r',label='Cadence')
    ax1.set_ylabel('Cadence (m)')
    
    ax2 = fig.add_subplot(413)
    line2 = plt.plot(power,color='k',label='Power')
    ax1.set_ylabel('Power (W)')

    ax2 = fig.add_subplot(414)
    line2 = plt.plot(3.6*myACDobj0.speedRefData,color='g',label='Power')
    ax1.set_ylabel('Speed (W)')
    plt.show()


    idx1 = input("Start Index (4677)")
    idx2 = input("End Index (4747)")
    plt.ioff()

    #idx1 = 4677
    #idx2 = 4747
    idx1 = int(idx1)
    idx2 = int(idx2)

    sPath = sPath[idx1:idx2+1]
    zElev = elev[idx1:idx2+1]
    power = power[idx1:idx2+1]


    myACDobj1 = SSACD.mySSACD(ACD0,sPath,zElev,power,mass_rider,mass_bike)
    
    fig = plt.figure()
    ax1 = fig.add_subplot(311)
    line1 = plt.plot(zElev,color='b',label='Elevation')
    ax1.set_ylabel('Elevation (m)')

    ax2 = fig.add_subplot(312)
    line2 = plt.plot(myACDobj1.speedRefData,color='r',label='Speed')
    ax1.set_ylabel('Cadence (m)')
    
    ax2 = fig.add_subplot(313)
    line2 = plt.plot(power,color='k',label='Power')
    ax1.set_ylabel('Power (W)')
    plt.show()

    myACDobj1.speed_fit_plot(myACDobj1.speedRefData,0.4)

    myACDobj1.ACDfitter_scipy()
    myACDobj1.plot_fitted_ACD()





    #errHist = mySSHRobj1.HRfitter(hr, power)
    #HRhatBestFit = mySSHRobj1.HRsim(power, hr[0])
    #mySSHRobj1.HRerrPlot(hr,HRhatBestFit)
    #mySSHRobj1.HRerrPlotWithBounds(hr, HRhatBestFit, power)
    #mySSHRobj1.HRerrPlotWith300WBounds(hr, HRhatBestFit, power)
    #mySSHRobj1.HRbeforeAfterPlot(hr,HRhatUnfitted,HRhatBestFit,errHist)
    print(myACDobj1)

    #fpathgpx = pathlib.Path(inputfilename)


    #HRobjFilename = HRobjDirectory+"SSHRobj_"+fpathgpx.stem+".hrobj"
    #print(HRobjFilename)
    #with open(HRobjFilename,'wb') as file_pickle:
    #    pickle.dump(mySSHRobj1,file_pickle)

    #file_pickle.close()