#!/Library/Frameworks/Python.framework/Versions/3.10/bin/python3


import mySSHeart as SSheart
import sys, pathlib
import parseGPX
import pickle
HRobjDirectory = '/Users/ryanblanchard/Documents/gpxFiles/HRobjs/'

if __name__ == '__main__':
    n = len(sys.argv)
    HRmin = 90
    HRmax = 187
    Pmax = 250
    print("Total arguments passed:", n)
    #for i in range(1, n):
        #print(i," ",sys.argv[i])
    if n == 1:
        inputfilename = '/Users/ryanblanchard/Documents/gpxFiles/LKN_Fondo_22.gpx'
        inputfilename = '/Users/ryanblanchard/Documents/gpxFiles/Ardiden309.gpx'
        #inputfilename = '/Users/ryanblanchard/Documents/gpxFiles/Ardiden265.gpx'
    if n==2:
        inputfilename = sys.argv[1];

    sPath, elev, power, cad, hr, time = parseGPX.parseGPX(inputfilename)
    mySSHRobj1 = SSheart.mySSHeart(HRmin, HRmax, Pmax)
    HRhatUnfitted = mySSHRobj1.HRsim(power, hr[0])
    
    errHist = mySSHRobj1.HRfitter(hr, power)
    
    HRhatBestFit = mySSHRobj1.HRsim(power, hr[0])
    #mySSHRobj1.HRerrPlot(hr,HRhatBestFit)
    mySSHRobj1.HRerrPlotWithBounds(hr, HRhatBestFit, power)
    mySSHRobj1.HRbeforeAfterPlot(hr,HRhatUnfitted,HRhatBestFit,errHist)    
    print(mySSHRobj1)
    
    fpathgpx = pathlib.Path(inputfilename)
    
    
    HRobjFilename = HRobjDirectory+"SSHRobj_"+fpathgpx.stem+".hrobj"
    print(HRobjFilename)
    with open(HRobjFilename,'wb') as file_pickle:
        pickle.dump(mySSHRobj1,file_pickle)
    
    file_pickle.close()
    
    
