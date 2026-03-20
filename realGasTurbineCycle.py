#!/Library/Frameworks/Python.framework/Versions/3.10/bin/python3
import numpy as np
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

plt.style.use('bmh')

def getGamma(gas):
    return gas.cp_mass/gas.cv_mass

def tauAndEtafromPolytropicEffandPi(gas,states,eta_poly,Pi0):
    nSteps = 100
    Ta,Pa = gas.TP
    ha = gas.enthalpy_mass
    sa = gas.entropy_mass
    sb_isentropic = sa
    Pb = Pa*Pi0
    gas.SP = sb_isentropic,Pb
    hb_isentropic = gas.enthalpy_mass
    Tb_isentropic = gas.T

    gas.TP = Ta,Pa
    Pi_i = np.power(Pi0,1/nSteps)
    Pi = Pa
    Ti = Ta
    print(Pa)
    for i in range(0,nSteps):
        gamma = getGamma(gas)
        if Pi0 >= 1.0:
            tau_i = np.power(Pi_i,(gamma-1)/(gamma*eta_poly))
        else:
            tau_i = np.power(Pi_i,(gamma-1)*eta_poly/gamma)
        Pb=Pi*Pi_i

        Tb = Ti*tau_i
        gas.TP = Tb,Pb
        h_sens,s_sens = sensible_HS(gas)
        states.append(gas.state,stage='',h_sensible=h_sens,s_sensible=s_sens)
        print(i, Pb, Tb, gamma)
        Pi = Pb
        Ti = Tb

    gas.TP = Tb,Pb
    hb = gas.enthalpy_mass
    Tb = gas.T
    eta = (hb_isentropic-ha)/(hb-ha)
    if Pi0 < 1:
        eta = 1/eta
    tau_isentropic = Tb_isentropic/Ta
    tau = Tb/Ta
    print(f'tau: {tau:.3f} tau_isentropic: {tau_isentropic:.3f} eta: {eta:.3f}')
    return tau,eta,Tb,hb

def combustorSimGivenT3andT4(gas,states,Tt3,Tt4,air,fuel):
    P3 = gas.P
    P4 = P3
    gas.TP = Tt3,P4

    firstPhiGuess = 1.3 #starting point equivalence ratio
    nMaxSteps = 20
    convergenceToleranceK  = .01 #maximum permissible combustion temperature tolerance
    gas.set_equivalence_ratio(phi=firstPhiGuess, fuel=fuel, oxidizer=air)
    gas.equilibrate('HP')
    T0 = gas.T
    if T0 < Tt4:
        return  gas.enthalpy_mass
    Ta = T0
    phia = firstPhiGuess
    while abs(Ta - Tt4) > convergenceToleranceK:
        gas.TP = Tt3,P4
        phib = phia*((Tt4-Tt3)/(Ta-Tt3))
        gas.set_equivalence_ratio(phi=phib, fuel=fuel, oxidizer=air)
        gas.equilibrate('HP')
        Tb = gas.T
        print(phia, phib, Tb, Tt4)
        Ta = Tb
        phia = phib

    gas.TP = Tt3,P4
    h3 = gas.enthalpy_mass
    gas.TP = Tt4,P4
    h4 = gas.enthalpy_mass
    qdotburner = h4-h3
    phiburner = phib
    fburner = gas.mixture_fraction(fuel=fuel,oxidizer=air)

    return fburner, phiburner, qdotburner, h4


def plotTSDiagram(gas, states):

    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    line1 = plt.semilogy(states.v,states.P,'-')
    plt.xlabel('Specific Volume (m^3/kg)')
    plt.ylabel('Pressure (Pa)')

    plt.legend()

    for i, txt in enumerate(states.stage):
        plt.annotate(txt,(states.v[i],states.P[i]))

    ax2 = fig.add_subplot(122)
    ax2.plot(states.s_sensible,states.h_sensible,'-')
    plt.xlabel('Specific Entropy (J/kg-K)')
    plt.ylabel('Sensible Enthalpy (J/kg)')

    for i, txt in enumerate(states.stage):
        plt.annotate(txt,(states.s_sensible[i],states.h_sensible[i]))
    ax2.legend()
    plt.show()

def sensible_HS(gas):
    T,P = gas.TP
    H = gas.enthalpy_mass
    S = gas.entropy_mass
    gas.TP  = 298.15,ct.one_atm
    Href = gas.enthalpy_mass
    Sref = gas.entropy_mass
    H = H-Href
    S = S-Sref
    gas.TP = T,P
    return H,S



if __name__ == '__main__':
    T0 = 300
    Tt4 = 2500
    P0 = ct.one_atm
    OPR = 40        #Overall Pressure Ratio
    eta_poly = .9
    eta_poly_list = [0.85,0.9,0.95]
    ctair  =ct.Solution('air.yaml')
    states_list = [ct.SolutionArray(ctair) for tmp in eta_poly_list]
    eta_thermal_list = np.ones_like(eta_poly_list)
    
    fig, axs = plt.subplots(1, 3, figsize=(18, 6))
    for i, eta_p  in enumerate(eta_poly_list):
        ctair.TP = T0,P0
        states_i = states_list[i]
        tau,eta,Tb,hb = tauAndEtafromPolytropicEffandPi(ctair,states_list[i],eta_p,1.8)
        eta_thermal_list[i]=eta
        axs[0].plot(states_i.T,states_i.P/P0,label=f'eta_p ={eta_p}')
        axs[1].plot(states_i.entropy_mass,states_i.enthalpy_mass,label=f'eta_p ={eta_p}')
    axs[0].set_xlabel('Temperature (K)')
    axs[0].set_ylabel('Pressure Ratio (-)')
    axs[1].set_xlabel('Entropy (J/kg-K)')
    axs[1].set_ylabel('Enthalpy (J/kg)')
    axs[2].plot(eta_poly_list,eta_thermal_list)
    axs[2].set_xlabel('Polytropic Efficiency (-)')
    axs[2].set_ylabel('Compressor Isentropic Efficiency (-)')
    axs[0].legend()
    axs[1].legend()
    plt.show()




    gas = ct.Solution('gri30.yaml')
    fuel = "CH4:1"
    gas.TP = T0, P0
    air = "O2:1.00,N2:3.77"
    h0 = gas.enthalpy_mass
    gas.set_equivalence_ratio(phi=0, fuel=fuel, oxidizer=air)
    states = ct.SolutionArray(gas,extra=['stage','h_sensible','s_sensible'])

    h_sens0,s_sens0 = sensible_HS(gas)
    states.append(gas.state,stage='Stage 0',h_sensible=h_sens0,s_sensible=s_sens0)
    tauC, etaC, Tt3,h3 = tauAndEtafromPolytropicEffandPi(gas,states,eta_poly,OPR)
    h_sens3,s_sens3 = sensible_HS(gas)
    states.append(gas.state,stage='Stage 3',h_sensible=h_sens3,s_sensible=s_sens3)

    fburner, phiburner, qdotburner, h4 = combustorSimGivenT3andT4(gas,states,Tt3,Tt4,air,fuel)
    h_sens4,s_sens4 = sensible_HS(gas)
    states.append(gas.state,stage='Stage 4',h_sensible=h_sens4,s_sensible=s_sens4)

    tauT, etaT, Tt5, h5 = tauAndEtafromPolytropicEffandPi(gas,states,eta_poly,1.0/OPR)
    h_sens5,s_sens5 = sensible_HS(gas)
    states.append(gas.state,stage='Stage 5',h_sensible=h_sens5,s_sensible=s_sens5)

    gas.TP = T0, P0
    gas.set_equivalence_ratio(phi=0, fuel=fuel, oxidizer=air)
    h_sens0,s_sens0 = sensible_HS(gas)
    states.append(gas.state,stage='',h_sensible=h_sens0,s_sensible=s_sens0)
    plotTSDiagram(gas, states)


    wdotCompressor = (h_sens3-h_sens0)
    wdotTurbine = (1+fburner)*(h_sens4-h_sens5)
    wdotNet = wdotTurbine - wdotCompressor
    eta_thermal = wdotNet/qdotburner
    print(f'wdotCompressor: {wdotCompressor:.3f} qdotburner: {qdotburner:.3f} wdotTurbine: {wdotTurbine:.3f} wdotNet: {wdotNet:.3f} eta_thermal: {eta_thermal:.3f}')


