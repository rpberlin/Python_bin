#!/Library/Frameworks/Python.framework/Versions/3.10/bin/python3

import liquidPiston_valveSweep as lpVS
import matplotlib.pyplot as plt
import numpy as np
plt.style.use('bmh')


if __name__ == '__main__':

    inlet_open_angle_Min = 0
    inlet_open_angle_Max = 50
    nInletPts = 10
    inlet_open_pts = np.linspace(inlet_open_angle_Min,inlet_open_angle_Max,nInletPts)

    outlet_open_angle_Min = 510
    outlet_open_angle_Max = 540
    nOutletPts = 10
    outlet_open_pts = np.linspace(outlet_open_angle_Min,outlet_open_angle_Max,nOutletPts)


    injector_open_angle_Min= 340
    injector_open_angle_Max = 360
    nInjectorAnglePts = 10
    injector_open_pts = np.linspace(injector_open_angle_Min,injector_open_angle_Max,nInjectorAnglePts)
    optimalEtas = np.zeros(nInjectorAnglePts)
    optimalInletOpenAngles = np.zeros(nInjectorAnglePts)
    optimalOutletOpenAngles = np.zeros(nInjectorAnglePts)
    for i, injector_open_angle in enumerate(injector_open_pts):
        etaGrid, inletGrid, outletGrid, bestEta, bestInletOpenAngle, bestOutletOpenAngle = lpVS.runValveSweep(inlet_open_pts, outlet_open_pts,140,injector_open_angle)
        optimalEtas[i] = bestEta
        optimalInletOpenAngles[i] = bestInletOpenAngle
        optimalOutletOpenAngles[i] = bestOutletOpenAngle


    fig = plt.figure(figsize=(15, 10))
    ax2 = fig.add_subplot(131)
    plt.plot(inlet_open_pts,optimalEtas)
    ax2.set_xlabel('Fuel Injector Open Angle (deg)')
    ax2.set_ylabel('Maximum Achievable Thermal Efficiency (-)')
    plt.title('Max Thermal Efficiency')

    ax2 = fig.add_subplot(132)
    plt.plot(inlet_open_pts,optimalInletOpenAngles)
    ax2.set_xlabel('Fuel Injector Open Angle (deg)')
    ax2.set_ylabel('Optimal Inlet Valve Open Angle (deg)')
    plt.title('Optimal Intake Valve Timing')

    ax2 = fig.add_subplot(133)
    plt.plot(inlet_open_pts,optimalOutletOpenAngles)
    ax2.set_xlabel('Fuel Injector Open Angle (deg)')
    ax2.set_ylabel('Optimal Exhaust Valve Open Angle (deg)')
    plt.title('Optimal Exhaust Valve Timing')

    plt.show()
