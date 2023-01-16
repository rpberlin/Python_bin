#!/Library/Frameworks/Python.framework/Versions/3.10/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import sys


def ColebrookEqualsZeroIfCorrect(Re,epsOverD,fd):
    return -2.0*np.sqrt(fd)*np.log10(epsOverD/3.7+2.5/(Re*np.sqrt(fd)))-1.0

def solveImplicitColebrook(Re,epsOverD,f0=0.03):
    fstep = 1.001
    solverURF = .3
    tolerance = 1e-6
    Err = ColebrookEqualsZeroIfCorrect(Re,epsOverD,f0)
    nSteps = 0
    while abs(Err) > tolerance:
        Err = ColebrookEqualsZeroIfCorrect(Re,epsOverD,f0)
        dErr_df = (ColebrookEqualsZeroIfCorrect(Re,epsOverD,f0*fstep)-Err)/(f0*(fstep-1.0))
        f1 = f0 - Err/dErr_df
        f0 = f0*(1-solverURF) + solverURF*f1
        
        nSteps += 1
        if nSteps > 1000:
            return 1.0
        #print('f0: ',f0, ' f1: ',f1,' Err: ',Err,' dEdf: ',dErr_df)

    return f0

def getFrictionFactor(Re,epsOverD=0,f0=0.03):
    ReTransition = 1000

    if Re < ReTransition:
        f_d = 64/Re
    else:
        f_d = solveImplicitColebrook(Re,epsOverD,f0)
    return f_d



def myLogSpace(a,b,n):
    alpha = np.power(b/a,1/(n-1))
    vec = np.zeros(n)
    for i in range(0,n):
        vec[i] = a*np.power(alpha,i)
    return vec


def plotMoodyChart(ReVals,epsVals):
    n = len(ReVals)
    m = len(epsVals)
    fvals_forChart = np.zeros([n,m])
    f0= 0.03
    for i,Re in enumerate(ReVals):
        for j,eps in enumerate(epsVals):
            f0 = getFrictionFactor(Re,eps,f0)
            fvals_forChart[i,j] = f0


    print(fvals_forChart)
    print(np.log10(10))
    for j, eps in enumerate(epsVals):
        plt.loglog(ReVals, fvals_forChart[:,j], '-', label=f'eps/D: {eps:.2f}')
    plt.xlabel('Reynolds Number (-)')
    plt.ylabel('Friction Factor (-)')
    plt.title('Moody Chart')
    plt.legend()

    plt.show()

def calcPressureLoss(Velocity,rho,mu,Diam,Length,eps=0,f0=.03):
    Re = rho*Velocity*Diam/mu
    f0 = getFrictionFactor(Re,eps/Diam,f0)
    return rho*Velocity*Velocity*Length*f0/(2*Diam), f0

def solvePipeVelocity(deltaP,rho,mu,Length,Diam,eps=0,fguess=0.03):
    if deltaP < 0:
        deltaP = -1*deltaP

    U0 = 0.01
    dP0,f0 = calcPressureLoss(U0,rho,mu,Diam,Length,eps,fguess)
    Ustep = 1.001
    solverURF = .5
    tolerance = 1e-4
    Err = (dP0 - deltaP)/deltaP
    nSteps =0
    while abs(Err) > tolerance:
        dPstar, f0 = calcPressureLoss(U0*Ustep,rho,mu,Diam,Length,eps,f0)
        ErrStar = (dPstar - deltaP)/deltaP
        dErr_dU = (ErrStar-Err)/(U0*(Ustep-1.0))
        U1 = U0 - Err/dErr_dU
        U0 = U0*(1-solverURF) + solverURF*U1
        U0 = max(1e-4,U0)
        dP0, f0 = calcPressureLoss(U0,rho,mu,Diam,Length,eps,f0)
        Err = (dP0-deltaP)/deltaP
        nSteps += 1
        if nSteps > 100:
            return 1.0,1.0
        #print(U0,U1,f0, Err)

    return U0, f0

def solveBuoyantFlowLoop(Qdot,Thot, rho, drho_dT, cp, mu, dmu_dT, U0, pipeDiam, totalPipeLength, deltaH, f0=.03):
    #f0 = 0.03
    gravity = 9.81                                 # NaK Bulk Velocity First Guess
    solverURF = .3
    tolerance = 1e-4
    Err = 1
    steps = 0
    maxSteps = 50
    Apipe = np.pi*0.25*np.power(pipeDiam,2)
    while abs(Err) > tolerance and steps < maxSteps:
        steps +=1
        mdot0 = rho*U0*Apipe                  # Flow Rate kg/s
        deltaT = Qdot/(mdot0*cp)   #temperature change across heat exchanger
        Tcold = Thot-deltaT
        delta_rho = -1*deltaT*drho_dT         # NaK density change responsible for driving natural convection
        deltaP = delta_rho*gravity*deltaH       # Natural Convection driving pressure
        rhoCold = rho + drho_dT*(Tcold-300)
        rhoHot =  rho + drho_dT*(Thot-300)
        muCold = mu + dmu_dT*(Tcold-300)
        muHot =  mu + dmu_dT*(Thot-300)
        rhoEff = 0.5*(rhoCold + rhoHot)   # average density in the flow loop
        muEff  = 0.5*(muCold + muHot)      # average viscosity in the flow loop    
        U1,f0 = solvePipeVelocity(deltaP,rhoEff,muEff,totalPipeLength,pipeDiam,0,f0)  #Calculate the updated velocity and
        U0 = U0*(1-solverURF)+U1*solverURF
        Err = U1-U0
    #print('Step: ',steps,' Err: ',Err,' U0: ',U0,' f0: ',f0,' U1: ',U1,' deltaT: ',deltaT,' Tcold: ',Tcold)
    return mdot0, U1, Tcold, rhoEff, muEff, f0

def rhoIDG_1atm_Tcelsius(MW,Tc):
    R_g         = 8315         #J/mol-K
    P_atm       = 101325        #Atmospheric Pressure
    rho = P_atm*MW/(R_g*(Tc+273.15))
    return rho
    

def heatExchangerDesign(Qdot,Thot,Treturn,Tambient,Dpipe):
    cp_air = 1040                       #Specific Heat of Air
    h_conv = 7                          #Convective Heat Transfer Coefficient of Natural Convection in Air
    MW_air = 28.2                       #Molecular Weight of Air
    gravity = 9.81     
    tube_over_pipe_sectionArea = 2          #Ratio of total tube to NaK pipe cross-section areas 
    blockage_ratio = 0.5                #Fraction of Air Cross Section Blocked by Tubes        
    rho_air = rhoIDG_1atm_Tcelsius(MW_air, Tambient)    #Air Density at Ambient Temperature and 1atm
    T_tubeAvg = 0.5*(Thot+Treturn)                      #Average Temparature from NaK Available for Heat Xfer
    deltaTairNaK = T_tubeAvg-Tambient                   #Average Temperature Difference Driving Heat Xfer from NaK to Air
    deltaTairair = 0.25*deltaTairNaK                    #Air temperature rise across heat exchanger limited to .25 of available delta T
    rho_airOut = rhoIDG_1atm_Tcelsius(MW_air, Tambient+deltaTairair)  #air density after heat exchanger
    AtubesOuter = Qdot/(deltaTairNaK*h_conv)            #Total Tube Surface Area required to reject Qdot given delta T and h
    mdotAir = Qdot/(deltaTairair*cp_air)                #Mass flow of air required to limit deltaTairair from exceeding constraint
    Across_pipe = 0.25*np.pi*np.power(Dpipe,2)          #Cross-sectional are of NaK pipe
    Across_Alltubes = tube_over_pipe_sectionArea*Across_pipe
    
    dMin = .005
    dMax = .05
    nDiams = 20
    bestMaxDimension = 1e6
    allDiams = myLogSpace(dMin, dMax, nDiams)   #list of all diameters to be checked
    for diam in allDiams:
        A_per_tube = np.pi*0.25*np.power(diam,2)    #Cross sectional Area per tube 
        N_tubes = Across_Alltubes/A_per_tube        #N tubes
        W = diam*np.sqrt(N_tubes)/blockage_ratio             #Total Width of Heat Exchanger 
        L =  AtubesOuter/(np.pi*diam*N_tubes)       #Tube Length To Give Required Surface Area
        A_open = L*W*blockage_ratio                 #Total Area of Heat Exchanger Open for Air Flow
        U_open = mdotAir/(rho_air*A_open)           #Air Velocity at minimum cross-section (Between Tubes)
        KE_air = 0.5*rho_air*np.power(U_open,2)     #Kinetic Energy of Air at minimum cross section
        H = W+KE_air/(0.5*gravity*(rho_air-rho_airOut))   #Stack Height required to supply necessary KE_air
        maxDimension = max(L,W,H)
        #print('Diam**: ',diam, N_tubes, L, W, H, maxDimension)
        if maxDimension < bestMaxDimension:                           #find diameter that minimizes total heat exchanger volume
            bestD, bestN, bestL, bestW, bestH, bestMaxDimension = diam, N_tubes, L, W, H, maxDimension    
            #print('Diam**: ',diam, N_tubes, L, W, H, maxDimension)
    return bestD, bestN, bestL, bestW, bestH, bestMaxDimension
    
    
if __name__ == '__main__':

    rho = 1000
    grav = 9.8
    h = 1.4
    Dmin = .00001
    Dmax = .1
    DVec = myLogSpace(Dmin,Dmax,20)
    UVec = 0*DVec
    eps = 0
    L = 20
    mu = 1e-3
    f0 = .03
    for i, D in enumerate(DVec):
        U,f0 = solvePipeVelocity(rho*grav*h,rho,mu,L,D,eps/D,f0)
        UVec[i] = U
        print(D,U)
    plt.semilogx(DVec,UVec)
    plt.show()



    if '--plot' in sys.argv:
        ReMin = 100
        ReMax = 1e8
        epsMin = 1e-6
        epsMax = .01
        nReSteps = 100
        nepsSteps = 5
        ReVals = myLogSpace(ReMin,ReMax,nReSteps)
        epsVals = myLogSpace(epsMin,epsMax,nepsSteps)
        print(ReVals,epsVals)

        plotMoodyChart(ReVals,epsVals)
