#!/Library/Frameworks/Python.framework/Versions/3.10/bin/python3
import matplotlib.pyplot as plt
from   matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import datetime as dt
import sys
import numpy as np
plt.style.use('bmh')

def calcDaylight(degNorth):
    earthTilt = 23.5*np.pi/180
    latBerlin = 52.52*np.pi/180
    latCharlotte = 35.22*np.pi/180
    latMuc =  66.5*np.pi/180
    minutes_per_day = 24*60
    sun_orbit_omega = 2*np.pi/365.25 #rad/day
    tilt_phis = []
    berlin_alphas = []
    clt_alphas = []
    muc_alphas = []
    berlin_night_fracs = []
    clt_night_fracs = []
    muc_night_fracs = []
    berlin_night_minutes = []
    clt_night_minutes = []
    muc_night_minutes = []
    berlin_new_minutes = []
    clt_new_minutes = []
    muc_new_minutes = []
    start_date = dt.date(2024, 12, 21)
    
    days = np.arange(0,183)
    dates = [start_date + dt.timedelta(days=int(day)) for day in days] 


    for day in days:
        tilt_phi = earthTilt*np.cos(sun_orbit_omega*day)
        tilt_phis.append(tilt_phi)
        berlin_alpha = np.pi+2*np.arcsin(np.sin(latBerlin)*np.sin(tilt_phi)/np.cos(latBerlin))
        clt_alpha = np.pi+2*np.arcsin(np.sin(latCharlotte)*np.sin(tilt_phi)/np.cos(latCharlotte))
        muc_alpha = np.pi+2*np.arcsin(np.sin(latMuc)*np.sin(tilt_phi)/np.cos(latMuc))
        berlin_alphas.append(berlin_alpha)
        clt_alphas.append(clt_alpha)
        muc_alphas.append(muc_alpha)
        berlin_night_frac = berlin_alpha/(2*np.pi)
        clt_night_frac = clt_alpha/(2*np.pi)
        muc_night_frac = muc_alpha/(2*np.pi)
        berlin_night_fracs.append(berlin_night_frac)
        clt_night_fracs.append(clt_night_frac)
        muc_night_fracs.append(muc_night_frac)
        berlin_night_minutes.append(minutes_per_day*berlin_night_frac)
        clt_night_minutes.append(minutes_per_day*clt_night_frac)
        muc_night_minutes.append(minutes_per_day*muc_night_frac)


    for i in range(len(days)-1):
        berlin_delta = berlin_night_minutes[i]-berlin_night_minutes[i+1]
        clt_delta = clt_night_minutes[i]-clt_night_minutes[i+1]
        muc_delta = muc_night_minutes[i]-muc_night_minutes[i+1]
        berlin_new_minutes.append(berlin_delta)
        clt_new_minutes.append(clt_delta)
        muc_new_minutes.append(muc_delta)
        print(dates[i], berlin_delta, clt_delta, muc_delta)
        

        

    fig2 = plt.figure(2)
    ax1 = fig2.add_subplot(221)
    line1 = plt.plot(days,np.degrees(tilt_phis))
    ax1.set_xlabel('days')
    ax1.set_ylabel('tilt angle (deg)')
    plt.legend()


    ax2 = fig2.add_subplot(222)
    line1 = plt.plot(days,np.degrees(berlin_alphas),label='Ber')
    line2 = plt.plot(days,np.degrees(clt_alphas),label='CLT')
    line3 = plt.plot(days,np.degrees(muc_alphas),label='MUC')
    ax2.set_xlabel('days')
    ax2.set_ylabel('Degrees of path in shadow)')
    plt.legend()

    ax3 = fig2.add_subplot(223)
    line1 = plt.plot(days,24*np.array(berlin_night_fracs),label='Ber')
    line2 = plt.plot(days,24*np.array(clt_night_fracs),label='CLT')
    line3 = plt.plot(days,24*np.array(muc_night_fracs),label='MUC')
    ax3.set_xlabel('days')
    ax3.set_ylabel('berlin night fraction (-)')
    plt.legend()

    ax4 = fig2.add_subplot(224)
    line1 = plt.plot(dates[0:-1],berlin_new_minutes,label='Berlin')
    line2 = plt.plot(dates[0:-1],clt_new_minutes,label='Charlotte')
    line3 = plt.plot(dates[0:-1],muc_new_minutes,label='Hamburg')
    ax4.set_xlabel('Datum')
    ax4.set_ylabel('Zusaetzliche Minuten Tageslicht Pro Tag')
    plt.legend()


    plt.show(block=False)
    plt.show()
    #print(type(sPath))










if __name__ == '__main__':
    n = len(sys.argv)
    print("Total arguments passed:", n)
    #for i in range(1, n):
        #print(i," ",sys.argv[i])
    #inputfilename = '/Users/ryanblanchard/Downloads/FulGaz_Swellendam_to_Bonnievale3.gpx'

    
    degNorth = 0
    calcDaylight(degNorth)

