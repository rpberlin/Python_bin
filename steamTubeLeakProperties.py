#!/Library/Frameworks/Python.framework/Versions/3.10/bin/python3
import numpy as np
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('bmh')

def getGamma(gas):
    return gas.cp_mass/gas.cv_mass


def steamTubeLeaks(PhighSteam, PlowNa,Thigh):
    rFluid = ct.Water()
    nIntermediatePoints = 40
    Pintermediates = np.linspace(PhighSteam,PlowNa,nIntermediatePoints)
    statesQ0 = ct.SolutionArray(rFluid,extra=['stage'])
    statesQ1 = ct.SolutionArray(rFluid,extra=['stage'])
    statesSuperHeat = ct.SolutionArray(rFluid,extra=['stage'])
    statesQ1f = ct.SolutionArray(rFluid,extra=['stage'])

    rFluid.PQ = PhighSteam, 0
    h = rFluid.enthalpy_mass
    statesQ0.append(rFluid.state,stage='Steam Side Q = 0')
    for P in Pintermediates:
        rFluid.HP = h,P
        statesQ0.append(rFluid.state,stage='')
    Qf_from_Q0 = rFluid.Q
    Tf_from_Q0 = rFluid.T
    statesQ0.append(rFluid.state,stage='Sodium Side From Qinit = 0')

    rFluid.PQ = PhighSteam, 1
    h = rFluid.enthalpy_mass
    statesQ1.append(rFluid.state,stage='Steam Side Q = 1')
    for P in Pintermediates:
        rFluid.HP = h,P
        statesQ1.append(rFluid.state,stage='')
    Qf_from_Q1 = rFluid.Q
    Tf_from_Q1 = rFluid.T
    statesQ1.append(rFluid.state,stage='Sodium Side From Qinit = 1')

    rFluid.TPQ = Thigh,PhighSteam,1
    h = rFluid.enthalpy_mass
    statesSuperHeat.append(rFluid.state,stage='Steam Side Superheated')
    for P in Pintermediates:
        rFluid.HP = h,P
        statesSuperHeat.append(rFluid.state,stage='')
    Qf_from_Superheater = rFluid.Q
    Tf_from_Superheater = rFluid.T
    statesSuperHeat.append(rFluid.state,stage='Sodium Side From Superheater')


    rFluid.PQ = PlowNa,1
    h = rFluid.enthalpy_mass
    statesQ1f.append(rFluid.state,stage='Sodium Side From Qfinal = 1')
    Qf_from_Q1f = rFluid.Q
    Tf_from_Q1f = rFluid.T
    for P in reversed(Pintermediates):
        rFluid.HP = h,P
        statesQ1f.append(rFluid.state,stage='')

    statesQ1f.append(rFluid.state,stage='Sodium Side From Qfinal = 1')


    output_str = '{:35s}{:>4.2f} {} {:>4.2f}{}'
    print(f'T+Q from Q0: {Tf_from_Q0:.1f} K\t {Qf_from_Q0:.3f} -')
    print(f'T+Q from Q1: {Tf_from_Q1:.1f} K\t {Qf_from_Q1:.3f} -')
    print(f'T+Q from Superheater: {Tf_from_Superheater:.1f} K\t {Qf_from_Superheater:.3f} -')
    print(f'T+Q from Qfinal = 1: {Tf_from_Q1f:.1f} K\t {Qf_from_Q1f:.3f} -')

    return rFluid, statesQ0, statesQ1, statesSuperHeat, statesQ1f

def plotHSDiagram(rFluid, statesQ0, statesQ1, statesSuperHeat,statesQ1f):
    extrapFactor = .1
    nPoints = 100

    states2 = ct.SolutionArray(rFluid)

    s1 = 0.5*statesQ0.P.min()
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
    plt.plot(statesQ0.s,statesQ0.T,label='Q0 Expansion')
    plt.plot(statesQ1.s,statesQ1.T,label='Q1 Expansion')
    plt.plot(statesSuperHeat.s,statesSuperHeat.T,label='Q0 Expansion')
    plt.plot(statesQ1f.s,statesQ1f.T,label='Expansion to Qf =1')
    plt.plot(states2.s,states2.T)
    #for i, txt in enumerate(states.stage):
    #    plt.annotate(txt,(states.s[i],states.T[i]))
    plt.xlabel('Specific Entropy (J/kg-K)')
    plt.ylabel('Temperature (K)')

    ax2 = fig.add_subplot(122)
    plt.loglog(statesQ0.v,statesQ0.P,label='Q0 Expansion')
    plt.loglog(statesQ1.v,statesQ1.P,label='Q1 Expansion')
    plt.loglog(statesSuperHeat.v,statesSuperHeat.P,label='Superheat Expansion')
    plt.loglog(statesQ1f.v,statesQ1f.P,label='Expansion to Qf =1')
    plt.loglog(states2.v,states2.P)
    #for i, txt in enumerate(states.stage):
    #    plt.annotate(txt,(statesQ0.v[i],statesQ0.P[i]))
    plt.xlabel('Specific Volume (m^3/kg)')
    plt.ylabel('Pressure (Pa)')
    plt.legend()


    plt.legend()




    plt.show()







if __name__ == '__main__':
    Plow = 101.3e3 #Pa
    Phigh = 3.0e6  #Pa
    Thigh = 350+273.15 #K

    rFluid, statesQ0, statesQ1, statesSuperHeat, statesQ1f = steamTubeLeaks(Phigh,Plow,Thigh)
    plotHSDiagram(rFluid,statesQ0, statesQ1, statesSuperHeat,statesQ1f )
