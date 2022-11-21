#!/Library/Frameworks/Python.framework/Versions/3.10/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import random
plt.style.use('bmh')

def newSwitch(switch, threshold):
    if random.randint(0,100) < threshold:
        return int(not switch)
    else:
        return switch




def  runSim():

    qringmax = 1000
    qbroilmax =1800
    qbottommax=700
    ringSwitch = 1
    broilSwitch = 1
    bottomSwitch = 1
    mShroud = 2
    mShroud2 = 3
    mShell1 = 5
    mShell2 = 5
    mShell4 = 6
    mShell5 = 7
    mRing = 3
    mBroil = 8
    mBottom = 3.75
    mLoad = 5
    Tamb = 300
    aShroud =1
    aShroud2 =1
    aShell1 = 5
    aShell2 = 4
    aShell3 = 3
    aShell4 = 2.5
    aShell5 = 2
    aRing = .1
    aBroil = .4
    aBottom = .3
    aLoad = 2
    nSteps =125

    q1 = np.zeros([nSteps,1])
    q2 = np.zeros([nSteps,1])
    q3 = np.zeros([nSteps,1])
    TShroud = np.zeros([nSteps,1])+Tamb
    TShroud2 = np.zeros([nSteps,1])+Tamb
    TShell1 = np.zeros([nSteps,1])+Tamb
    TShell2 = np.zeros([nSteps,1])+Tamb
    TShell3 = np.zeros([nSteps,1])+Tamb
    TShell4 = np.zeros([nSteps,1])+Tamb
    TShell5 = np.zeros([nSteps,1])+Tamb
    Tring = np.zeros([nSteps,1])+Tamb
    Tbroil = np.zeros([nSteps,1])+Tamb
    TBottom = np.zeros([nSteps,1])+Tamb
    Tload = np.zeros([nSteps,1])+Tamb
    time = np.linspace(0,nSteps-1,nSteps)

    for i in range(1,nSteps):
        Tamb = 500+nSteps+50*np.sin(.01*nSteps)
        ringSwitch = newSwitch(ringSwitch,5)
        broilSwitch = newSwitch(broilSwitch,5)
        bottomSwitch = newSwitch(bottomSwitch,5)
        print(ringSwitch," ",broilSwitch," ",bottomSwitch)
        qringshroud1 = (Tring[i-1]-TShroud[i-1])*aShroud*aRing
        qringshroud2 = (Tring[i-1]-TShroud2[i-1])*aShroud2*aRing
        qringShell1 = (Tring[i-1]-TShell1[i-1])*aShell1*aRing
        qringShell2 = (Tring[i-1]-TShell2[i-1])*aShell2*aRing
        qringShell3 = (Tring[i-1]-TShell3[i-1])*aShell3*aRing

        q1[i]=ringSwitch*qringmax
        q2[i]=broilSwitch*qbroilmax
        q3[i]=bottomSwitch*qbottommax

        Tdotring = (ringSwitch*qringmax-qringshroud1-qringshroud2-qringShell1-qringShell2-qringShell3)/mRing

        qbroilshroud1 = (Tbroil[i-1]-TShroud[i-1])*aShroud*aBroil
        qbroilshroud2 = (Tbroil[i-1]-TShroud[i-1])*aShroud2*aBroil
        qbroilShell4  = (Tbroil[i-1]-TShell4[i-1])*aShell4*aBroil
        qbroilShell5  = (Tbroil[i-1]-TShell5[i-1])*aShell5*aBroil

        Tdotbroil = (broilSwitch*qbroilmax-qbroilshroud1-qbroilshroud2-qbroilShell4-qbroilShell5)/mBroil

        qbottomshroud1 = (TBottom[i-1]-TShroud[i-1])*aShroud*aBottom
        qbottomshroud2 = (TBottom[i-1]-TShroud[i-1])*aShroud2*aBottom
        qbottomShell3  = (TBottom[i-1]-TShell3[i-1])*aShell3*aBottom
        qbottomShell4  = (TBottom[i-1]-TShell4[i-1])*aShell4*aBottom
        qbottomShell5  = (TBottom[i-1]-TShell5[i-1])*aShell5*aBottom

        Tdotbottom = (bottomSwitch*qbottommax-qbottomshroud1-qbottomshroud2-qbottomShell3-qbottomShell4-qbottomShell5)/mBottom

        qShell12 = (TShell1[i-1]-TShell2[i-1])*aShell2*aShell1*.05
        qShell13 = (TShell1[i-1]-TShell3[i-1])*aShell3*aShell1*.05
        qShell14 = (TShell1[i-1]-TShell4[i-1])*aShell4*aShell1*.05
        qShell15 = (TShell1[i-1]-TShell5[i-1])*aShell5*aShell1*.05
        qShell1Amb = (TShell1[i-1]-Tamb)*aShell1*aShell1*.05

        TdotShell1 = (qringShell1 - qShell12-qShell13-qShell14-qShell15-qShell1Amb)/mShell1


        qShell23 = (TShell2[i-1]-TShell3[i-1])*aShell3*aShell2*.1
        qShell24 = (TShell2[i-1]-TShell4[i-1])*aShell4*aShell2*.1
        qShell25 = (TShell2[i-1]-TShell5[i-1])*aShell5*aShell2*.1
        qShell2Amb = (TShell2[i-1]-Tamb)*aShell2*aShell2*.01

        TdotShell2 = (qringShell2 + qShell12 - qShell23-qShell24-qShell25-qShell2Amb)/mShell2

        qShell34 = (TShell3[i-1]-TShell4[i-1])*aShell4*aShell3*.1
        qShell35 = (TShell3[i-1]-TShell5[i-1])*aShell5*aShell3*.1
        qShell3Amb = (TShell3[i-1]-Tamb)*aShell3*aShell3*.01

        TdotShell3 = (qringShell3 + qShell13 + qShell23 - qShell34 - qShell35-qShell3Amb)/mShell4

        qShell45 = (TShell4[i-1]-TShell5[i-1])*aShell5*aShell4*.1
        qShell4Amb = (TShell4[i-1]-Tamb)*aShell4*aShell4*.01

        TdotShell4 = (qShell14 + qShell24 + qShell34 - qShell45-qShell4Amb)/mShell4

        qShell5Amb = (TShell5[i-1]-Tamb)*aShell5*aShell5*.01

        TdotShell5 = (qShell15 + qShell25 + qShell35 + qShell45-qShell5Amb)/mShell5

        Tring[i]=Tring[i-1]+Tdotring
        Tbroil[i]=Tbroil[i-1]+Tdotbroil
        TBottom[i]=TBottom[i-1]+Tdotbottom
        TShell1[i]=TShell1[i-1]+TdotShell1
        TShell2[i]=TShell2[i-1]+TdotShell2
        TShell3[i]=TShell3[i-1]+TdotShell3
        TShell4[i]=TShell4[i-1]+TdotShell4
        TShell5[i]=TShell5[i-1]+TdotShell5


    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    line1 = plt.plot(time,q1,label='Heater1')
    line2 = plt.plot(time,q2,label='Heater2')
    line3 = plt.plot(time,q3,label='Heater3')
    #ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('Heat Input (W)')
    plt.title('Inputs')
    plt.legend()

    ax2 = fig.add_subplot(212)
    line1 = plt.plot(time,Tring,label='Tring')
    line2 = plt.plot(time,Tbroil,label='Tbroil')
    line3 = plt.plot(time,TBottom,label='Tbottom')
    line4 = plt.plot(time,TShell1,label='Tshell1')
    #line5 = plt.plot(time,TShell2,label='Tshell2')
    line5 = plt.plot(time,TShell3,label='Tshell3')
    #ine5 = plt.plot(time,TShell4,label='Tshell4')
    line5 = plt.plot(time,TShell5,label='Tshell5')
    ax2.set_xlabel('Time (s)')
    ax2.set_ylabel('Temperature (K)')
    plt.title('System Response')
    plt.legend()
    ax2.legend(loc='lower right')
    plt.show()




    return X_table, T_flame

def plotTandX(T,X,species_to_track):
    nPhiPts, nSpecies = X.shape
    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    line1 = plt.plot(phiPts,T,label='Temperature')
    ax1.set_xlabel('Equivalence Ratio (-)')
    ax1.set_ylabel('Temperature (K)')
    plt.title('Temperature')
    plt.legend()

    ax2 = fig.add_subplot(122)
    for i, specie in enumerate(species_to_track):
        line2 = plt.semilogy(phiPts,X[:,i],label=specie)
    ax2.set_ylim([1e-6,1])
    ax2.set_xlabel('Equivalence Ratio (-)')
    ax2.set_ylabel('Mole Fraction (-)')
    plt.legend()
    plt.title('Equilibrium Mole Fractions')
    plt.show()





if __name__ == '__main__':
    runSim()
