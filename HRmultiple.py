#!/Library/Frameworks/Python.framework/Versions/3.10/bin/python3


import mySSHeart as SSheart
import sys, pathlib, os
import parseGPX
import pickle
import matplotlib.pyplot as plt
from  datetime import datetime

HRobjDirectory = '/Users/ryanblanchard/Documents/gpxFiles/HRobjs/'

def make_HR_comparison_plots( date_list , gpx_names_list, watts_at_140bpm_list, watts_at_185bpm_list, bpm_at_300W_list):
    
    date_objects = [datetime.strptime(date, '%Y-%m-%d') for date in date_list]
    sorted_data = sorted(zip(date_objects,date_list, gpx_names_list, watts_at_140bpm_list, watts_at_185bpm_list, bpm_at_300W_list))
    date_objects,date_list, gpx_names_list, watts_at_140bpm_list, watts_at_185bpm_list,bpm_at_300W_list  = zip(*sorted_data)

    # Create a figure with 3 subplots
    fig, axs = plt.subplots(1, 3, figsize=(16, 8))

    # Plot data1 in the first subplot
    axs[0].plot(date_objects, watts_at_140bpm_list, marker='o', linestyle='',label='watts_at_140bpm_list')
    axs[0].set_ylabel('Watts (-)')
    axs[0].set_ylim([100,200])

    # Plot data2 in the second subplot
    axs[1].plot(date_objects, watts_at_185bpm_list, marker='o', linestyle='',label='watts_at_185bpm_list')
    axs[1].set_ylabel('Watts (-)')
    axs[1].set_ylim([230,400])

    # Plot data3 in the third subplot
    axs[2].plot(date_objects, bpm_at_300W_list, marker='o', linestyle='',label='bpm_at_300W_list')
    axs[2].set_ylabel('Beats Per Min (beats/min)')
    axs[2].set_ylim([165,200])

    # Set x-axis labels using names
    for ax in axs:
        #ax.set_xticks(range(len(gpx_names_list)))
        ax.set_xticks(date_list, rotation=60)
        ax.legend()

    # Adjust layout to prevent overlapping
    plt.tight_layout()

    # Show the plots
    plt.show()
    return 


if __name__ == '__main__':
    gpx_names_list = []
    #gpx_names_list = ['/Users/ryanblanchard/Documents/gpxFiles/2025/FulGaz_2022_TdF_Stage_20_Time_Trial_Lacapelle_Marival_to_Rocamadour_7.gpx','/Users/ryanblanchard/Documents/gpxFiles/2025/FulGaz_Luz_Ardiden_Under_overs_Short260a.tcx']
    
    HRmin = 90
    HRmax = 187
    Pmax = 250

    
    date_list = []
    gpx_labels_list = []
    bpm_at_300W_list  = []
    watts_at_140bpm_list = []
    watts_at_185bpm_list = []

    for i in range(1,len(sys.argv)):
        gpx_names_list.append(sys.argv[i])
    
    n = len(gpx_names_list)

    if n < 2:
        print("Not enough arguments, correct usage: HRmultiple.py track1.gpx track2.gpx ...")
        sys.exit(0) 


    for inputfilename in gpx_names_list:

        if not inputfilename.lower().endswith('.gpx') or inputfilename.lower().endswith('.tcx'):
            print('Sketchy filename')
            #sys.exit(0)

        #date_i = parseGPX.extract_gpx_date(inputfilename)
        #print(inputfilename,' : ',date_i)     
        sPath, elev, power, cad, hr, time, date_i = parseGPX.parseGPXorTCX(inputfilename)
        fpathgpx = pathlib.Path(inputfilename)
        HRobjFilename = HRobjDirectory+"SSHRobj_"+fpathgpx.stem+".hrobj"
        if os.path.exists(HRobjFilename):
            print("HRobjFilename exists already: loading from file: ")
            with open(HRobjFilename, 'rb') as pickle_file:
                mySSHRobj_i = pickle.load(pickle_file)
        else:
            #sPath, elev, power, cad, hr, time, date_i = parseGPX.parseGPXorTCX(inputfilename)
            mySSHRobj_i = SSheart.mySSHeart(HRmin, HRmax, Pmax,hr,power)
            if len(hr) > 1:
                HRhatUnfitted = mySSHRobj_i.HRsim(power, hr[0])
                res  = mySSHRobj_i.HRfitter_scipy()
                xFinal = res.x 
                mySSHRobj_i.setHRclassParamsFromX(xFinal)
                with open(HRobjFilename,'wb') as file_pickle:
                    pickle.dump(mySSHRobj_i,file_pickle)
                    file_pickle.close()
        stem = inputfilename.split('.')[0]
        date_list.append(date_i)
        gpx_labels_list.append(stem)
        watts_at_185bpm_list.append(mySSHRobj_i.getSSPowFor185BPM())
        watts_at_140bpm_list.append(mySSHRobj_i.getSSPowFor140BPM())
        bpm_at_300W_list.append(mySSHRobj_i.getHRfor300W())

    
    print(len(date_list))    
    make_HR_comparison_plots(date_list,gpx_labels_list, watts_at_140bpm_list, watts_at_185bpm_list, bpm_at_300W_list)

    

    


        

   
