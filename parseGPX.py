#!/Library/Frameworks/Python.framework/Versions/3.10/bin/python3

import xml.etree.ElementTree as ET
import math
import numpy as np
import sys
#def printRecur(root):
#    """Recursively prints the tree."""
#    print ('tag= ',root.tag, ' attrib= ', root.attrib, ' text= ',root.text )
#
#    for elem in root.getchildren():
#        printRecur(elem)
#    indent -= 4

def parseGPX(gpxFilename):
    power = [] 
    lat = []
    long = []
    cad = []
    elev = []
    hr = []
    nPts = 0
    print('Filename = ',gpxFilename)
    tree = ET.parse(gpxFilename)
    root = tree.getroot()

    for x in root.iter():
        if 'trkpt' in str(x.tag):
            nPts+=1
            trkpt = x
            lat.append(float(trkpt.get('lat').strip()))
            long.append(float(trkpt.get('lon').strip()))
            for y in trkpt.iter():
                if 'ele' in str(y.tag):
                    #print('elev= ', y.text)
                    elev.append(float(y.text.strip()))
                if 'power' in str(y.tag):
                    #print('power= ', y.text)
                    power.append(float(y.text))
                if 'cad' in str(y.tag):
                    #print('cad= ', y.text)
                    cad.append(float(y.text))
                if 'hr' in str(y.tag):
                    #print('hr= ', y.text)
                    hr.append(float(y.text));
    #print('nPts = ',nPts)

    Rearth = 6371000
    sPathSum = 0
    sPath = [0]
    dS = []
    for i in range(1,nPts):
        dx = Rearth*math.cos((3.1415/180)*(lat[i]))*(3.1415/180)*(long[i]-long[i-1]);
        dy = Rearth*(3.1415/180)*(lat[i]-lat[i-1]);
        dS =math.sqrt(dx*dx+dy*dy);
        sPathSum = sPathSum+dS;
        sPath.append(sPathSum)
        #print('dx= ',dx,' dy= ',dy,' sPath=', sPathSum)

#    trk = root[1]
#    trkseg = trk[2]
#    for trkpnt in trkseg
#        latlong = trkpnt.attrib
#        elev =

    #printRecur(root)

#    for child in root:
#        print(child.tag,child.attrib)

    #print('Done')
    #print('Npts: ',nPts,' lenPath ',len(sPath),' lenElev ',len(elev),' nPow ',len(power), 'lenCad ',len(cad), 'lenHr', len(hr))
    sPath = np.array(sPath)
    elev = np.array(elev)
    power = np.array(power)
    cad = np.array(cad)
    hr = np.array(hr)

    for i in range(0,nPts-len(power)):
        power = np.append(power,[0])


    #print('Npts: ',nPts,' lenPath ',len(sPath),' lenElev ',len(elev),' nPow ',len(power), 'lenCad ',len(cad), 'lenHr', len(hr))

    return sPath, elev, power, cad, hr

def parseTCX(tcxFilename):
    power = []
    lat = []
    long = []
    cad = []
    elev = []
    hr = []
    nPts = 0
    file1 = open(tcxFilename, 'r')
    Lines = file1.readlines()
    nLines = len(Lines)
    keywords = ['DistanceMeters','HeartRateInBeatsPerMinute_t','Cadence','Watts']
    delimiters = ['>','<','/','Value']
    for i, line in enumerate(Lines):
        for keyword in keywords:
            if keyword in line:
                 bigLine = " ".join(bigLine.split(delimiter))

 


    #print('nPts = ',nPts)

    Rearth = 6371000
    sPathSum = 0
    sPath = [0]
    dS = []
    for i in range(1,nPts):
        dx = Rearth*math.cos((3.1415/180)*(lat[i]))*(3.1415/180)*(long[i]-long[i-1]);
        dy = Rearth*(3.1415/180)*(lat[i]-lat[i-1]);
        dS =math.sqrt(dx*dx+dy*dy);
        sPathSum = sPathSum+dS;
        sPath.append(sPathSum)
        #print('dx= ',dx,' dy= ',dy,' sPath=', sPathSum)

#    trk = root[1]
#    trkseg = trk[2]
#    for trkpnt in trkseg
#        latlong = trkpnt.attrib
#        elev =

    #printRecur(root)

#    for child in root:
#        print(child.tag,child.attrib)

    #print('Done')
    #print('Npts: ',nPts,' lenPath ',len(sPath),' lenElev ',len(elev),' nPow ',len(power), 'lenCad ',len(cad), 'lenHr', len(hr))
    sPath = np.array(sPath)
    elev = np.array(elev)
    power = np.array(power)
    cad = np.array(cad)
    hr = np.array(hr)

    for i in range(0,nPts-len(power)):
        power = np.append(power,[0])


    #print('Npts: ',nPts,' lenPath ',len(sPath),' lenElev ',len(elev),' nPow ',len(power), 'lenCad ',len(cad), 'lenHr', len(hr))

    return sPath, elev, power, cad, hr





if __name__ == '__main__':
    n = len(sys.argv)
    print("Total arguments passed:", n)
    for i in range(1, n):
        print(i," ",sys.argv[i])

    if n < 2:
        print("Correct Usage parseGPX.py filename.gpx")
    inputfilename = sys.argv[1];
    parseGPX(inputfilename)
    
