#!/Library/Frameworks/Python.framework/Versions/3.10/bin/python3

import sys
import math
from parseGPX import parseGPX
from plot_PMdata import calcPMdata
from VAMify import VAMify
from compare2GPXroutes import compare2GPXroutes
#from compare2GPXroutes import getQtyAtS
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.ticker import StrMethodFormatter, NullFormatter
import numpy as np
import scipy

def float_input(input,default_value=None):
    #print("float: ",float(input)," ", float(input)+1)
    try:
        return float(input)
    except:
        return default_value
    
def interactive_speed_from_power():
    power_in = float_input(input("Enter Rider Power: (250W) "),250)
    myMass = float_input(input("Enter Rider Mass in kg (83.0): "),83)
    bikeMass = float_input(input("Enter Bike Mass: (7.5)"),7.5)
    slope = float_input(input("Enter Slope (.03)"),.03)
    CdA = float_input(input("Enter CdA: (0.45) "),0.45)
    speed_kph = speed_solver(myMass+bikeMass,CdA,slope,power_in)
    print("speed KPH:",speed_kph)

def power_calc(mass_total,cda,slope, speed_kph):
    speed_ms = speed_kph/3.6
    mechanical_efficiency = 0.95
    rho_air = 1.18
    loss_factor = 1/mechanical_efficiency
    return loss_factor * ( mass_total * speed_ms * slope * 9.81 + 0.5 * rho_air * cda * speed_ms *speed_ms * speed_ms)

def power_surplus(mass_total,cda,slope,speed_kph, power_in):
    power_out = power_calc(mass_total,cda,slope, speed_kph)
    return power_in-power_out

def speed_solver(mass_total,cda,slope,power_in):
    eps_kph = .01
    speed_0 = 32.0
    urf = .7 
    delta_speed = 1e6
    while abs(delta_speed) > eps_kph:
        err = power_surplus(mass_total,cda,slope,speed_0, power_in)
        derr = err -power_surplus(mass_total,cda,slope, speed_0 + eps_kph , power_in)
        delta_speed = err*eps_kph/derr
        speed_0 = speed_0+urf*delta_speed
        print(speed_0)
    return speed_0
    
def interactive_time_from_power_and_distance():
    distance = float_input(input("Enter Distance: (160.09 km) "),160.09)
    power_in = float_input(input("Enter Rider Power: (250W) "),250)
    myMass = float_input(input("Enter Rider Mass in kg (83.0): "),83)
    bikeMass = float_input(input("Enter Bike Mass: (7.5)"),7.5)
    slope = float_input(input("Enter Slope (.03)"),.03)
    CdA = float_input(input("Enter CdA: (0.45) "),0.45)
    speed_kph = speed_solver(myMass+bikeMass,CdA,slope,power_in)
    elapsed_time_minutes = 60*distance/speed_kph

    print('Elapsed Time: ',elapsed_time_minutes,' minutes at ',speed_kph,' kph') 
    powers = np.linspace(0.8*power_in,1.2*power_in,21)
    drags  = np.linspace(0.8*CdA,1.2*CdA,21)
    psweep_speeds = np.zeros_like(powers)
    dragsweep_speeds = np.zeros_like(powers)
    for i, power in enumerate(powers):
        psweep_speeds[i] = speed_solver(myMass+bikeMass,CdA,slope,power)

    for i, drag in enumerate(drags):
        dragsweep_speeds[i] = speed_solver(myMass+bikeMass,drag,slope,power_in)

    psweep_times = 60*distance/psweep_speeds
    dragsweep_times = 60*distance/dragsweep_speeds

    fig = plt.figure(figsize=(15, 10))
    ax2 = fig.add_subplot(121)
    plt.plot(powers,psweep_times)
    ax2.set_xlabel('Power (W)')
    ax2.set_ylabel('Elapsed Time (min)')
    ax3 = fig.add_subplot(122)
    plt.plot(drags,dragsweep_times)
    ax3.set_xlabel('Drag Coefficient CdA (-)')
    ax3.set_ylabel('Elapsed Time (min)')
    plt.show()

    return elapsed_time_minutes






def speed_given_power_mass_drag_and_slope():
    basePower = input("Enter Rider Power: (250W) ")
    myMass = input("Enter Rider Mass in (83.0kg): ")
    bikeMass = input("Enter Rider Mass: (7.5kg)")
    CdA = input("Enter Rider Mass: (0.45) ")


    return 0 

    





if __name__ == '__main__':
    # speed_kph = speed_solver(83+7.5,0.5,0.03,250)
    # slopes = np.linspace(0,0.1,21)
    # speeds_low = 0*slopes
    # speeds_nom = 0*slopes
    # speeds_high = 0*slopes
    # power_offset_factor = 1.2
    # for i, slope in enumerate(slopes):
    #     speeds_low[i] = speed_solver(83+7.5,0.5,slope,250/power_offset_factor)
    #     speeds_nom[i] = speed_solver(83+7.5,0.5,slope,250)
    #     speeds_high[i] = speed_solver(83+7.5,0.5,slope,250*power_offset_factor)
    #     print(slope,speed_kph)

    # plt.figure()
    # plt.plot(100*slopes,speeds_low*slopes*1000,label=str(250/power_offset_factor))
    # plt.plot(100*slopes,speeds_nom*slopes*1000,label='250')
    # plt.plot(100*slopes,speeds_high*slopes*1000,label = str(250*power_offset_factor))
    # plt.xlabel('Slope (%)')
    # plt.ylabel('Speed (kph)')
    # plt.legend()
    # plt.show()
    
    print("Menu:")
    print("1- Speed From Power, Weight, Slope, and Drag")
    print("2- Power for  Speed, Weight, Slope, and Drag")
    print("3- Duration for Distance, Power, Weight, Slope, and Drag")
    runMode = input("Select your choice:")
    match runMode:
        case "1":
            interactive_speed_from_power()
        case "2":
            interactive_power_for_speed()
        case "3":
            interactive_time_from_power_and_distance()



    # power_in = float_input(input("Enter Rider Power: (250W) "),250)
    # myMass = float_input(input("Enter Rider Mass in kg (83.0): "),83)
    # bikeMass = float_input(input("Enter Bike Mass: (7.5)"),7.5)
    # slope = float_input(input("Enter Slope (.03)"),.03)
    # CdA = float_input(input("Enter CdA: (0.45) "),0.45)
    # speed_kph = speed_solver(myMass+bikeMass,CdA,slope,power_in)
    # print('Speed: ', speed_kph)
