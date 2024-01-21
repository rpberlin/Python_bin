#!/Library/Frameworks/Python.framework/Versions/3.10/bin/python3
import numpy as np
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('bmh')


def  R134aVaporCompressionCycle(mdot,Plow,Phigh,compressor_efficiency):
    rFluid = ct.Hfc134a()
    rFluid.PQ = Plow,1
    states = ct.SolutionArray(rFluid)
    
    #State 1 - before compressor
    P1 = Plow
    rFluid.PQ = P1, 1
    T1 = rFluid.T
    h1 = rFluid.enthalpy_mass
    rho1, P1, Q1 = rFluid.DPQ
    s1 = rFluid.entropy_mass
    states.append(rFluid.state)
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
    states.append(rFluid.state)
    #states.append(rFluid,P=P2,T=T2,Q=Q2,H=h2,rho=rho2,s=s2)
    
    #state 3a start condensing
    rFluid.PQ = P2, 1
    states.append(rFluid.state)
    
    #state 3 - after condenser
    P3 = Phigh 
    rFluid.PQ = P3,0
    h3 = rFluid.enthalpy_mass
    s3 = rFluid.entropy_mass
    T3 = rFluid.T
    rho3, P3, Q3 = rFluid.DPQ
    states.append(rFluid.state)
   # states.append(rFluid,P=P3,T=T3,Q=Q3,H=h3,rho=rho3,s=s3)
    
    #State 4 - after throttling valve
    h4 = h3 #adiabatic
    P4 = Plow
    rFluid.HP = h4, P4
    s4 = rFluid.entropy_mass
    T4 = rFluid.T
    rho4, P4, Q4 = rFluid.DPQ
    states.append(rFluid.state)
    #states.append(rFluid,P=P4,T=T4,Q=Q4,H=h4,rho=rho4,s=s4)
    
    # Complete Loop Back to State 1
    rFluid.PQ = P1, 1
    states.append(rFluid.state)
    

    QdotL = mdot*(h1-h4)
    WdotIn = mdot*(h2-h1)
    QdotH = mdot*(h2-h3)
    COP = QdotL/WdotIn
    print(h1,h2,h3,h4)
    print(QdotL,WdotIn,QdotH,COP)
    temp =1    
    
    plt.figure()
    #plt.plot([s1, s2, s3, s4],[h1,h2,h3,h4])
    plt.plot(states.s,states.T)
    plt.show()
    
    
    return rFluid, states
        
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

       
    plt.figure()
    plt.plot(states.s,states.T,label='Vapor Compression Cycle')
    plt.plot(states2.s,states2.T)
    plt.xlabel('Specific Entropy (J/kg-K)')
    plt.ylabel('Temperature (K)')
    plt.legend()
    plt.show()

        
        
    
    
    
        
    

    
    
    



if __name__ == '__main__':
    Plow = 0.14e6 #Pa
    Phigh = 0.8e6  #Pa
    mdot = 0.05 #kg/s 
    compressor_efficiency = 0.8
    rFluid, states = R134aVaporCompressionCycle(mdot,Plow,Phigh,compressor_efficiency)
    plotHSDiagram(rFluid,states)
    
