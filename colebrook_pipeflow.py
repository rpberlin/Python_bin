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
    while abs(Err) > tolerance:
        Err = ColebrookEqualsZeroIfCorrect(Re,epsOverD,f0)
        dErr_df = (ColebrookEqualsZeroIfCorrect(Re,epsOverD,f0*fstep)-Err)/(f0*(fstep-1.0))
        f1 = f0 - Err/dErr_df
        f0 = f0*(1-solverURF) + solverURF*f1
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
    sign =1
    if deltaP < 0:
        sign = -1
        deltaP = sign*deltaP
    U0 = 0.01
    dP0,f0 = calcPressureLoss(U0,rho,mu,Diam,Length,eps,fguess)
    Ustep = 1.001
    solverURF = .5
    tolerance = 1e-4
    Err = (dP0 - deltaP)/deltaP
    while abs(Err) > tolerance:
        dPstar, f0 = calcPressureLoss(U0*Ustep,rho,mu,Diam,Length,eps,f0)
        ErrStar = (dPstar - deltaP)/deltaP
        dErr_dU = (ErrStar-Err)/(U0*(Ustep-1.0))
        U1 = U0 - Err/dErr_dU
        U0 = U0*(1-solverURF) + solverURF*U1
        dP0, f0 = calcPressureLoss(U0,rho,mu,Diam,Length,eps,f0)
        Err = (dP0-deltaP)/deltaP
        #print(U0,U1,f0, Err)

    return sign*U0, f0

def solveBuoyantFlowLoop(Qdot,Thot, rho, drho_dT, cp, mu, U0, pipeDiam, totalPipeLength, deltaH ):
    f0 = 0.03
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
        U1,f0 = solvePipeVelocity(deltaP,rho,mu,totalPipeLength,pipeDiam,0,f0)  #Calculate the updated velocity and
        U0 = U0*(1-solverURF)+U1*solverURF
        Err = U1-U0
        print('Step: ',steps,' Err: ',Err,' U0: ',U0,' f0: ',f0,' U1: ',U1,' deltaT: ',deltaT,' Tcold: ',Tcold)
    return mdot0, U1, Tcold, f0






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
