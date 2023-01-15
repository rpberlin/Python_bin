#!/Library/Frameworks/Python.framework/Versions/3.10/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import sys
import colebrook_pipeflow as cPipe


gravity = 9.81          #m/s^2      @375 C
rho_NaK = .789          #kg/m^3     @375 C
drho_dT_NaK = -.0002356 #kg/m^3-K   @375 C
cp_NaK =    896         #J/kg-K     @375 C
mu_NaK =    2.35e-4     #Pa-s       @375 C
k_NaK  =    22.4        #W/m-K      @375 C
Pr_NaK = cp_NaK*mu_NaK/k_NaK

deltaH  = 10
deltaL  = 10
totalPipeLength = 2*(deltaH+deltaL)

pipeDiam = 1 #insideDiameter of the pipe
pipeA = 0.25*np.pi*np.power(pipeDiam,2)  #pipe Cross Sectional Area

QdotShutdown = 180e3 #Heat Rejection
QdotNominal = 30e3 #Heat Rejection

ThotNominal = 371 #sodium temperatue in deg C under nominal conditions
ThotShutdown = 377 #sodium temperatue in deg C under shutdown conditions
Tambient = 20 #ambient temperatue in deg C

deltaT = 100                 #Good Starting Point for guessed Temperature Drop
mdot0 = QdotShutdown/(deltaT*cp_NaK)
U0 = mdot0/(pipeA*rho_NaK)


mdotShutdown, Ushutdown, TreturnShutdown, f0 = cPipe.solveBuoyantFlowLoop(QdotShutdown, ThotShutdown, rho_NaK, drho_dT_NaK, cp_NaK, mu_NaK, U0, pipeDiam, totalPipeLength, deltaH)
'''
while abs(Err) > tolerance and steps < maxSteps:
    steps +=1 
    mdot0 = rho_NaK*U0*pipeA                   # NaK Flow Rate kg/s
    deltaT = QdotShutdown/(mdot0*cp_NaK)   # NaK temperature change across heat exchanger
    Treturn = ThotShutdown-deltaT
    delta_rho = -1*deltaT*drho_dT_NaK         # NaK density change responsible for driving natural convection
    deltaP = delta_rho*gravity*deltaH       # Natural Convection driving pressure
    U1,f0 = cPipe.solvePipeVelocity(deltaP,rho_NaK,mu_NaK,totalPipeLength,pipeDiam,0,f0)  #Calculate the updated velocity and 
    U0 = U0*(1-solverURF)+U1*solverURF
    Err = U1-U0
    #print('Step: ',steps,' Err: ',Err,' U0: ',U0,' f0: ',f0,' U1: ',U1,' deltaT: ',deltaT,' Treturn: ',Treturn)
'''
