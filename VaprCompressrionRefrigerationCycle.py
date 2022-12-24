#!/Library/Frameworks/Python.framework/Versions/3.10/bin/python3
import numpy as np
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('bmh')


def  R134aVaporCompressionCycle(mdot,Plow,Phigh,compressor_efficiency):
    rFluid = ct.Hfc134a()
    rFluid.PQ = Plow,1
    states = ct.SolutionArray(rFluid,extra=['stage'])
    
    #State 1 - before compressor
    P1 = Plow
    rFluid.PQ = P1, 1
    T1 = rFluid.T
    h1 = rFluid.enthalpy_mass
    rho1, P1, Q1 = rFluid.DPQ
    s1 = rFluid.entropy_mass
    states.append(rFluid.state,stage='Stage 1')
    #states.append(rFluid,P=P1,T=T1,Q=Q1,H=h1,rho=rho1,s=s1)
    
    
    #compressor
    P2 = Phigh 
    s2_isentropic = s1 
    rFluid.SP = s2_isentropic, P2
    h2_isentropic = rFluid.enthalpy_mass
    deltah12_isentropic = h2_isentropic - h1
    
    
    
    #state 2 after compressor
    h2 = h1+deltah12_isentropic/compressor_efficiency
    rFluid.HP = h2,P2
    s2 = rFluid.entropy_mass
    T2 = rFluid.T
    rho2, P2, Q2 = rFluid.DPQ
    states.append(rFluid.state,stage='Stage 2')
    #states.append(rFluid,P=P2,T=T2,Q=Q2,H=h2,rho=rho2,s=s2)
    
    #state 3a start condensing
    rFluid.PQ = P2, 1
    states.append(rFluid.state,stage='Stage 3a')
    
    #state 3 - after condenser
    P3 = Phigh 
    rFluid.PQ = P3,0
    h3 = rFluid.enthalpy_mass
    s3 = rFluid.entropy_mass
    T3 = rFluid.T
    rho3, P3, Q3 = rFluid.DPQ
    states.append(rFluid.state,stage='Stage 3')
   # states.append(rFluid,P=P3,T=T3,Q=Q3,H=h3,rho=rho3,s=s3)
    
    #State 4 - after throttling valve
    h4 = h3 #adiabatic
    P4 = Plow
    rFluid.HP = h4, P4
    s4 = rFluid.entropy_mass
    T4 = rFluid.T
    rho4, P4, Q4 = rFluid.DPQ
    states.append(rFluid.state,stage='Stage 4')
    #states.append(rFluid,P=P4,T=T4,Q=Q4,H=h4,rho=rho4,s=s4)
    
    
    #State4star - turbine instead of throttling valve
    turbine_efficiency = .8
    s4star_isentropic = s3
    P4star = Plow
    rFluid.SP = s4star_isentropic, Plow
    h4star_isentropic = rFluid.enthalpy_mass 
    h4star = h3-turbine_efficiency*(h3-h4star_isentropic)
    WdotOutmax = mdot*(h3-h4star)
    #rFluid.HP = h4star,P4star
    #states.append(rFluid.state)
        
    
    
    # Complete Loop Back to State 1
    rFluid.PQ = P1, 1
    states.append(rFluid.state,stage='')
    

    QdotL = mdot*(h1-h4)
    WdotIn = mdot*(h2-h1)
    QdotH = mdot*(h2-h3)
    COP = QdotL/WdotIn

      
    
    return rFluid, states, QdotL, QdotH, WdotIn, COP, WdotOutmax
        
def plotTSDiagramWithDome(rFluid, states):
    extrapFactor = .1
    nPoints = 50
    
    states2 = ct.SolutionArray(rFluid)
    
    P1 = 0.5*states.P.min()
    P2 = 0.98*rFluid.critical_pressure
    
    PVec = np.linspace(P1,P2,nPoints)

  
    rFluid.TP = rFluid.critical_temperature, rFluid.critical_pressure
    scr = rFluid.entropy_mass
    
    for i, P in enumerate(PVec):
        rFluid.PQ = P,0
        states2.append(rFluid.state)
    
    rFluid.TP = rFluid.critical_temperature, rFluid.critical_pressure
    states2.append(rFluid.state)
        
    for i, P in enumerate(reversed(PVec)):
        rFluid.PQ = P,1
        states2.append(rFluid.state)

       
    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    line1 = plt.plot(states.s,states.T,'-o',label='Vapor Compression Cycle')
    plt.plot(states2.s,states2.T)
    plt.xlabel('Specific Entropy (J/kg-K)')
    plt.ylabel('Temperature (K)')
    plt.legend()
    
    for i, txt in enumerate(states.stage):
        plt.annotate(txt,(states.s[i],states.T[i]))
        
    ax2 = fig.add_subplot(122)
    line2 = plt.plot(states.enthalpy_mass,states.P,'-o',label='Vapor Compression Cycle')
    plt.plot(states2.enthalpy_mass,states2.P)
    plt.xlabel('Enthalpy (J/kg)')
    plt.ylabel('Pressure (Pa)')
    plt.legend()
    
    for i, txt in enumerate(states.stage):
        plt.annotate(txt,(states.enthalpy_mass[i],states.P[i]))
    
    plt.show()
    

        
   
if __name__ == '__main__':
    Plow = 0.14e6 #Pa
    Phigh = 0.8e6  #Pa
    mdot = 0.05 #kg/s 
    compressor_efficiency = 0.8
    rFluid, states,  QdotL, QdotH, WdotIn, COP, WdotOutmax = R134aVaporCompressionCycle(mdot,Plow,Phigh,compressor_efficiency)
    plotTSDiagramWithDome(rFluid,states)
    
