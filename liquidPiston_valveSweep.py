#!/Library/Frameworks/Python.framework/Versions/3.10/bin/python3

import liquidPiston_try1 as p1
import matplotlib.pyplot as plt
import numpy as np
plt.style.use('bmh')


def runValveSweep(inlet_open_pts=np.linspace(0,40,5), outlet_open_pts=np.linspace(512,532,3),inlet_open_delta=140,injector_open_angle= 355, injector_open_delta=50, rpm=3000):
    etaGrid = np.zeros([len(inlet_open_pts),len(outlet_open_pts)])
    inletGrid = np.zeros([len(inlet_open_pts),len(outlet_open_pts)])
    outletGrid = np.zeros([len(inlet_open_pts),len(outlet_open_pts)])

    bestEta = 0
    bestInletOpenAngle = 0
    bestOutletOpenAngle = 0

    for i, inlet_open_angle in enumerate(inlet_open_pts):
        for j, outlet_open_angle in enumerate(outlet_open_pts):
            eta =0.48
            states, Q , W , MEP, eta , CO_emission ,xticks, ticklabels = p1.runLPtry1_Simulation(rpm, inlet_open_angle, inlet_open_delta,outlet_open_angle,injector_open_angle,injector_open_delta)
            etaGrid[i,j]=eta
            inletGrid[i,j]=inlet_open_angle
            outletGrid[i,j]=outlet_open_angle
            if eta > bestEta:
                bestEta = eta
                bestInletOpenAngle = inlet_open_angle
                bestOutletOpenAngle = bestOutletOpenAngle


    print(x,outlet_open_pts)


    return etaGrid, inletGrid, outletGrid, bestEta, bestInletOpenAngle, bestOutletOpenAngle


def plotValveSweepResults(inlet_open_pts,outlet_open_pts,etaGrid, inletGrid, outletGrid, bestEta, bestInletOpenAngle, bestOutletOpenAngle):
    fig = plt.figure(figsize=(15, 10))
    plt.contourf(inletGrid,outletGrid,etaGrid)
    plt.xlabel('Inlet Valve Open Angle (deg)')
    plt.ylabel('Exhaust Valve Open Angle(deg)')
    plt.title(f'Best Thermal Efficiency = {bestEta:.3f} (-)')
    cbar = plt.colorbar()
    plt.show()
    pass




if __name__ == '__main__':
    #rpm = 3000
    inlet_open_angle_Min = 20
    inlet_open_angle_Max = 50
    nInletPts = 5
    inlet_open_pts = np.linspace(inlet_open_angle_Min,inlet_open_angle_Max,nInletPts)

    outlet_open_angle_Min = 500
    outlet_open_angle_Max = 540
    nOutletPts = 5
    outlet_open_pts = np.linspace(outlet_open_angle_Min,outlet_open_angle_Max,nOutletPts)

    etaGrid, inletGrid, outletGrid, bestEta, bestInletOpenAngle, bestOutletOpenAngle = runValveSweep(inlet_open_pts, outlet_open_pts)
    plotValveSweepResults(inlet_open_pts,outlet_open_pts,etaGrid, inletGrid, outletGrid, bestEta, bestInletOpenAngle, bestOutletOpenAngle)
