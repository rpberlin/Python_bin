#!/Library/Frameworks/Python.framework/Versions/3.10/bin/python3

import xml.etree.ElementTree as ET
import gpxpy
import math
import numpy as np
import sys
import re 
from  datetime import datetime
#def printRecur(root):
#    """Recursively prints the tree."""
#    print ('tag= ',root.tag, ' attrib= ', root.attrib, ' text= ',root.text )
#
#    for elem in root.getchildren():
#        printRecur(elem)
#    indent -= 4

import xml.etree.ElementTree as ET

def parseGPXorTCX(filename):
    sPath = None
    elev  = None
    power = None
    cad   = None
    hr    = None
    date  = None 
    elapsed_times = None
    if filename.lower().endswith('.gpx'):
        sPath, elev, power, cad, hr, elapsed_times, = parseGPX(filename)
        date = extract_gpx_date(filename)
    elif filename.lower().endswith('.tcx'):
        sPath, elev, power, cad, hr, elapsed_times = parse_tcx(filename)
        date = extract_tcx_date(filename)
    return sPath, elev, power, cad, hr, elapsed_times, date 

def extract_tcx_date(tcx_file):
    namespaces = {    'tcx': 'http://www.garmin.com/xmlschemas/TrainingCenterDatabase/v2'}
    tree = ET.parse(tcx_file)
    root = tree.getroot()
    # Extract the date again
    timestamp = root.find('.//tcx:Track/tcx:Trackpoint/tcx:Time', namespaces)
    date_str = timestamp.text.split("T")[0]
    return date_str

def extract_gpx_date(gpx_file):
    try:
        # Parse the GPX file
        tree = ET.parse(gpx_file)
        root = tree.getroot()

        # Namespace for GPX elements
        ns = {'gpx': 'http://www.topografix.com/GPX/1/1'}

        # Find the time element within the metadata
        time_element = root.find(".//gpx:metadata/gpx:time", namespaces=ns)

        if time_element is not None:
            # Extract and return the date
            date_str = time_element.text.split("T")[0]
            return date_str
        else:
            print("No date information found in the GPX file.")
            return None

    except ET.ParseError as e:
        print(f"Error parsing GPX file: {e}")
        return None

def parse_gpx2(file_path):
    # Initialize lists to store data
    latitudes = []
    longitudes = []
    elevations = []
    elapsed_times = []
    power_data = []
    cadence_data = []
    heart_rate_data = []

    # Parse GPX file
    with open(file_path, 'r') as gpx_file:
        gpx = gpxpy.parse(gpx_file)

        for track in gpx.tracks:
            for segment in track.segments:
                for point in segment.points:
                    # Extract data from each track point
                    latitudes.append(point.latitude)
                    longitudes.append(point.longitude)
                    elevations.append(point.elevation)
                    elapsed_times.append(point.time.timestamp())  # Convert time to timestamp
                    power_value = 0
                    cadence_value = 0
                    heart_rate_value = 0
                    for extension in point.extensions:
                        if 'power' in extension.tag:
                            power_value = extension.text
                        if 'cadence' in extension.tag:
                            cadence_value = extension.text
                        if 'heartRate' in extension.tag:
                            heart_rate_value = extension.text

                    

                    power_data.append(power_value)
                    cadence_data.append(cadence_value)
                    heart_rate_data.append(heart_rate_value)    

    return latitudes, longitudes, elevations, power_data, elapsed_times


def get_sPath(lat,long):
    Rearth = 6371000
    sPathSum = 0
    sPath = [0]
    dS = []
    for i in range(1,len(lat)):
        dx = Rearth*math.cos((3.1415/180)*(lat[i]))*(3.1415/180)*(long[i]-long[i-1]);
        dy = Rearth*(3.1415/180)*(lat[i]-lat[i-1]);
        dS =math.sqrt(dx*dx+dy*dy);
        sPathSum = sPathSum+dS;
        sPath.append(sPathSum)
    return sPath 


def parseGPX(gpxFilename):
    power = [] 
    lat = []
    long = []
    cad = []
    elev = []
    hr = []
    elapsed_times = []
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
                    if 'none' in y.text.lower():
                        power.append(0)
                    else:
                        power.append(float(y.text))
                if 'cad' in str(y.tag):
                    #print('cad= ', y.text)
                    cad.append(float(y.text))
                if 'hr' in str(y.tag):
                    #print('hr= ', y.text)
                    hr.append(float(y.text));
                if 'time' in str(y.tag):
                    timeslist = re.split('T|Z',y.text)
                    time_point = datetime.strptime(timeslist[1], '%H:%M:%S')
                    #print('time: ',time_point.timestamp())
                    elapsed_times.append(time_point.timestamp())
    #print('nPts = ',nPts)
    sPath = get_sPath(lat,long)
    sPath = np.array(sPath)
    elev = np.array(elev)
    power = np.array(power)
    cad = np.array(cad)
    hr = np.array(hr)
    elapsed_times = np.array(elapsed_times)
    if len(elapsed_times) == 0:
        elapsed_times = np.arange(len(sPath))
    elapsed_times = elapsed_times - elapsed_times[0]

    for i in range(0,nPts-len(power)):
        power = np.append(power,[0])


    #print('Npts: ',nPts,' lenPath ',len(sPath),' lenElev ',len(elev),' nPow ',len(power), 'lenCad ',len(cad), 'lenHr', len(hr))

    return sPath, elev, power, cad, hr, elapsed_times



# Function to parse the TCX file for latitude, longitude, distance, heart rate, cadence, and watts
def parse_tcx(tcxfilename):
    namespaces = {'tcx': 'http://www.garmin.com/xmlschemas/TrainingCenterDatabase/v2'}
    tree = ET.parse(tcxfilename)
    root = tree.getroot()
    power = []
    sPath = []
    cad = []
    elev = []
    hr = []
    elapsed_times = []
    
    # Loop through the trackpoints
    for track in root.findall('.//tcx:Track/tcx:Trackpoint', namespaces):
       
        for y in track.iter():
            #print(y.tag)
            if 'Distance' in str(y.tag):
                sPath.append(float(y.text))
            if 'ele' in str(y.tag):
                #print('elev= ', y.text)
                elev.append(float(y.text.strip()))
            if 'Watts' in str(y.tag):
                #print('power= ', y.text)
                power.append(float(y.text))
            if 'Cadence' in str(y.tag):
                #print('cad= ', y.text)
                cad.append(float(y.text))
            if 'HeartRateBpm' in str(y.tag):
                #print('hr= ', y[0].text)
                hr.append(float(y[0].text))
            if 'Time' in str(y.tag):
                timeslist = re.split('T|Z',y.text)
                time_point = datetime.strptime(timeslist[1], '%H:%M:%S.%f0')
                #print('time: ',time_point.timestamp())
                elapsed_times.append(time_point.timestamp())
    
    elapsed_times = [elapsed_times[i]-elapsed_times[0] for i in range(0,len(elapsed_times))]
           

    return sPath, elev, power, cad, hr, elapsed_times
        


if __name__ == '__main__':
    n = len(sys.argv)
    #tcxfile = '/Users/ryanblanchard/Documents/gpxFiles/2025/FulGaz_Luz_Ardiden_Under_overs_Short260a.tcx'
   #parse_tcx2(tcxfile)

    print("Total arguments passed:", n)
    for i in range(1, n):
        print(i," ",sys.argv[i])

    if n < 2:
        print("Correct Usage parseGPX.py filename.gpx")
    inputfilename = sys.argv[1];
    parseGPX(inputfilename)
    
