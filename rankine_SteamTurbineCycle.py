#!/Library/Frameworks/Python.framework/Versions/3.10/bin/python3
import numpy as np
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('bmh')

def getGamma(gas):
    return gas.cp_mass/gas.cv_mass


def rankine_steamCycle(mdot,Plow,Phigh,Thigh,compressor_eta, turbine_eta):
    rFluid = ct.Water()
    rFluid.PQ = Plow,1
    states = ct.SolutionArray(rFluid,extra=['stage'])

    #State 1 - before compressor
    P1 = Plow
    rFluid.PQ = P1, 0               #Steam is saturated
    T1 = rFluid.T
    h1 = rFluid.enthalpy_mass
    rho1, P1, Q1 = rFluid.DPQ
    s1 = rFluid.entropy_mass
    states.append(rFluid.state,stage='Stage 1')
    #states.append(rFluid,P=P1,T=T1,Q=Q1,H=h1,rho=rho1,s=s1)


    #Perfect compressor
    P2 = Phigh
    s2_isentropic = s1
    rFluid.SP = s2_isentropic, P2
    h2_isentropic = rFluid.enthalpy_mass
    deltah12_isentropic = h2_isentropic - h1


    #State2 after Real Compressor after compressor
    h2 = h1+deltah12_isentropic/compressor_eta
    rFluid.HP = h2,P2
    s2 = rFluid.entropy_mass
    T2 = rFluid.T
    rho2, P2, Q2 = rFluid.DPQ
    states.append(rFluid.state,stage='Stage 2')
    #states.append(rFluid,P=P2,T=T2,Q=Q2,H=h2,rho=rho2,s=s2)


    #Heat Until Onset of Boiling
    rFluid.PQ = P2, 0
    states.append(rFluid.state,stage='')

    #Heat Until Completion of Boiling
    rFluid.PQ = P2, 1
    states.append(rFluid.state,stage='')

    #state 3 - after superheater
    P3 = Phigh
    T3 = Thigh
    rFluid.TPQ = T3,P3,1
    s3 = rFluid.entropy_mass
    h3 = rFluid.enthalpy_mass
    states.append(rFluid.state,stage='Stage 3')
   # states.append(rFluid,P=P3,T=T3,Q=Q3,H=h3,rho=rho3,s=s3)

   #Perfect turbine
    s4_isentropic = s3
    rFluid.SP = s4_isentropic, Plow
    h4_isentropic = rFluid.enthalpy_mass
    deltah34_isentropic = h3-h4_isentropic

    #State 4 - after turbine
    h4 = h3 - deltah34_isentropic*turbine_eta #adiabatic
    P4 = Plow
    rFluid.HP = h4, P4
    s4 = rFluid.entropy_mass
    T4 = rFluid.T
    rho4, P4, Q4 = rFluid.DPQ
    states.append(rFluid.state,stage='Stage 4')
    #states.append(rFluid,P=P4,T=T4,Q=Q4,H=h4,rho=rho4,s=s4)

    # Complete Loop Back to State 1
    rFluid.PQ = P1, 0
    states.append(rFluid.state,stage='')

    #QdotL = mdot*(h1-h4)

    WdotIn = mdot*(h2-h1)
    QdotIn = mdot*(h3-h2)
    WdotOut = mdot*(h3-h4)
    eta_thermal = (WdotOut-WdotIn)/QdotIn
    eta_carnot = 1-T1/T3
    eta_2nd    = eta_thermal/eta_carnot
    output_str = '{:35s}{:>4.2f} {}'
    print(output_str.format('Work Input:',WdotIn/1e3, 'kW/kg'))
    print(output_str.format('Heat Input:',QdotIn/1e3, 'kW/kg'))
    print(output_str.format('Work Output:',WdotOut/1e3, 'kW/kg'))
    print(output_str.format('Thermal Eff.:',eta_thermal, '(-)'))
    print(output_str.format('Max Possible Eff.:',eta_carnot, '(-)'))
    print(output_str.format('2nd Law Eff.:',eta_2nd, '(-)'))
    print(output_str.format('Stage4 Quality:',Q4, '(-)'))
    print(output_str.format('Temperature Min:',T1, '(K)'))

    print(T1,T2,T3,T4)
    print(eta_thermal,WdotIn,WdotOut,QdotIn)
    temp =1



    return rFluid, states, Q4, T1

def plotHSDiagram(rFluid, states):
    extrapFactor = .1
    nPoints = 50

    states2 = ct.SolutionArray(rFluid)

    s1 = 0.5*states.P.min()
    s2 = 0.98*rFluid.critical_pressure
    sVec = np.linspace(s1,s2,nPoints)


    rFluid.TP = rFluid.critical_temperature, rFluid.critical_pressure
    scr = rFluid.entropy_mass

    for i, P in enumerate(sVec):
        rFluid.PQ = P,0
        states2.append(rFluid.state)

    rFluid.TP = rFluid.critical_temperature, rFluid.critical_pressure
    states2.append(rFluid.state)

    for i, P in enumerate(reversed(sVec)):
        rFluid.PQ = P,1
        states2.append(rFluid.state)


    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    plt.plot(states.s,states.T,label='Rankine Cycle')
    plt.plot(states2.s,states2.T)
    for i, txt in enumerate(states.stage):
        plt.annotate(txt,(states.s[i],states.T[i]))
    plt.xlabel('Specific Entropy (J/kg-K)')
    plt.ylabel('Temperature (K)')

    ax2 = fig.add_subplot(122)
    plt.loglog(states.v,states.P,label='Rankine Cycle')
    plt.loglog(states2.v,states2.P)
    for i, txt in enumerate(states.stage):
        plt.annotate(txt,(states.v[i],states.P[i]))
    plt.xlabel('Specific Volume (m^3/kg)')
    plt.ylabel('Pressure (Pa)')
    plt.legend()


    plt.legend()




    plt.show()















if __name__ == '__main__':
    Plow = 75.0e3 #Pa
    Phigh = 3.0e6  #Pa
    Thigh = 350+273.15 #K
    mdot = 1.0 #kg/s
    compressor_efficiency = 1.0
    rFluid, states, Q4, Tmin = rankine_steamCycle(mdot,Plow,Phigh,Thigh,compressor_efficiency,1.0)
    plotHSDiagram(rFluid,states)
