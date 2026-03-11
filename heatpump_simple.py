import numpy as np
import CoolProp.CoolProp as cp
import pint
ureg = pint.UnitRegistry()
Q_ = ureg.Quantity
import matplotlib.pyplot as plt
plt.style.use('bmh')

def dict_unit_cleanup(dict):
    for key in dict.keys():
        tmplist = dict[key]
        if isinstance(tmplist[0], Q_):
                    units = tmplist[0].units
                    dict[key] = Q_([tmp.m_as(units) for tmp in tmplist], units)
        else:
            # leave floats, strings, etc. as plain lists or convert to array
            dict[key] = tmplist  # or np.array(tmplist) if you want numpy arrays
    return dict

def P_from_TQ(T,q,fluid):
    return Q_(cp.PropsSI('P','T',T.m_as('K'),'Q',q,fluid), 'Pa')

def get_triple_TP(fluid):
    P_triple = Q_(cp.PropsSI('ptriple',fluid),'Pa')
    T_triple = Q_(cp.PropsSI('Ttriple',fluid),'K')
    return T_triple, P_triple 

def get_critical_TP(fluid):
    P_crit = Q_(cp.PropsSI('Pcrit',fluid),'Pa')
    T_crit = Q_(cp.PropsSI('Tcrit',fluid),'K')
    return T_crit, P_crit 

def V_from_TQ(T,q,fluid):
    return Q_(cp.PropsSI('V', 'T', T.m_as('K'), 'Q', q, 'R32'),'m^3/kg')

def V_from_TS(T,s,fluid):
    return Q_(cp.PropsSI('V', 'T', T.m_as('K'), 'S', s.m_as('J/kg/K'), 'R32'),'m^3/kg')

def H_from_TP(T,P,fluid):
    return Q_(cp.PropsSI('H', 'T', T.m_as('K'), 'P', P.m_as('Pa'), 'R32'),'J/kg*K',fluid)

def H_from_PT(T,P,fluid):
    return Q_(cp.PropsSI('H', 'T', T.m_as('K'), 'P', P.m_as('Pa'), 'R32'),'J/kg*K',fluid)

def S_from_TP(T,P,fluid):
    return Q_(cp.PropsSI('S','T',T.m_as('K'),'P',P.m_as('Pa'),fluid),'J/kg/K')

def H_from_TQ(T,q,fluid):
    return Q_(cp.PropsSI('H','T',T.m_as('K'),'Q',q,fluid), 'J/kg')

def H_from_PQ(p,q,fluid):
    return Q_(cp.PropsSI('H','P',p.m_as('Pa'),'Q',q,fluid), 'J/kg')

def S_from_TQ(T,q,fluid):
    return Q_(cp.PropsSI('S','T',T.m_as('K'),'Q',q,fluid), 'J/kg/K')

def H_from_TS(T,s,fluid):
    return Q_(cp.PropsSI('H','T',T.m_as('K'),'S',s.m_as('J/kg/K'),fluid), 'J/kg')

def P_from_TS(T,s,fluid):
    return Q_(cp.PropsSI('P','T',T.m_as('K'),'S',s.m_as('J/kg/K'),fluid), 'Pa')

def T_from_PQ(P,q,fluid):
    return Q_(cp.PropsSI('T','P',P.m_as('Pa'),'Q',q,fluid), 'K')

def Q_from_PH(P,h,fluid):
    return cp.PropsSI('Q','P',P.m_as('Pa'),'H',h.m_as('J/kg'),fluid)


def  HeatPumpCycle(fluid, Q_dot,T_low,T_high):
    #state 0 compressor inlet
    T1 = T_low
    Q1 = 1.0 
    P1 = P_from_TQ(T1,Q1,fluid)
    s1 = S_from_TQ(T1,Q1,fluid)
    h1 = H_from_TQ(T1,Q1,fluid)
    v1 = V_from_TQ(T1,Q1,fluid)

    T2 = T_high
    s2 = s1
    h2 = H_from_TS(T2,s2,fluid)
    P2 = P_from_TS(T2,s2,fluid)
    Q2 = 1 
    v2 = V_from_TS(T2,s2,fluid)

    P2a = P2
    Q2a = 1
    T2a = T_from_PQ(P2a,Q2a,fluid)
    h2a = H_from_TQ(T2a,Q2a,fluid)
    s2a = S_from_TQ(T2a,Q2a,fluid)
    v2a= V_from_TQ(T2a,Q2a,fluid)

    P3 = P2
    Q3 = 0.0
    T3 = T_from_PQ(P3,Q3,fluid)
    h3 = H_from_TQ(T3,Q3,fluid)
    s3 = S_from_TQ(T3,Q3,fluid)
    v3 = V_from_TQ(T3,Q3,fluid)

    P4 = P1
    h4 = h3
    Q4 = Q_from_PH(P4,h4,fluid)
    T4 = T_from_PQ(P4,Q4,fluid)
    s4 = S_from_TQ(T4,Q4,fluid)
    v4 = V_from_TQ(T4,Q4,fluid)

    m_dot = Q_dot/(h2-h3)
    w_dot = m_dot*(h2-h1)
    cop = Q_dot/w_dot 
    print(f"m_dot: {m_dot.to('kg/s'):.3f} w_dot: {w_dot.to('kW'):.3f} COP: {cop.to(''):.2f}")

    states = {
    'T': [T1, T2, T2a, T3, T4,T1],
    'P': [P1,P2, P2a, P3, P4, P1],
    's': [s1,s2,s2a,s3,s4, s1],
    'Q': [Q1,Q2,Q2a,Q3,Q4,Q1],
    'h': [h1,h2,h2a,h3,h4,h1],
    'v': [v1,v2,v2a,v3,v4,v1],
    'labels': ['1','2','2a','3','4'],
    'cop' : [cop],
    'w_dot' : [w_dot],
    'T_amb' : [T1]
    }
    states = dict_unit_cleanup(states)
    return states, w_dot, cop
        
def plotloopandHSDiagram(fluid, dictlist):
    PFactorMin = .1
    PFactorMax = 10
    nPoints = 300
    
    T_crit, P_crit = get_critical_TP(fluid)
    T_triple, P_triple = get_triple_TP(fluid)
    T_sweep = np.linspace(T_crit,T_triple,nPoints)
    h_sweepl = H_from_TQ(T_sweep,0,fluid)
    h_sweepg = H_from_TQ(T_sweep,1.0,fluid)
    P_sweepl = P_from_TQ(T_sweep,0,fluid)
    P_sweepg = P_from_TQ(T_sweep,1.0,fluid)
    s_liq = S_from_TQ(T_sweep,0,fluid)
    s_vap = S_from_TQ(T_sweep,1.0,fluid)
    v_liq = V_from_TQ(T_sweep,0,fluid)
    v_vap = V_from_TQ(T_sweep,1.0,fluid)


    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))
    axes[1,0].plot(s_liq.to('J/kg/K'),h_sweepl.to('J/kg'),'r-',label='liquid')
    axes[1,0].plot(s_vap.to('J/kg/K'),h_sweepg.to('J/kg'),'b-',label='gas')
    axes[1,1].semilogx(v_liq,P_sweepl.to('atm'),'r-',label='liquid')
    axes[1,1].semilogx(v_vap,P_sweepg.to('atm'),'b-',label='gas')
    cops = []
    T_ambs = []
    for dict in dictlist:
        s_pts = dict['s']
        T_pts = dict['T']
        P_pts = dict['P']
        h_pts = dict['h']
        v_pts = dict['v']
        cops.append(dict['cop'])
        T_ambs.append(dict['T_amb'])

        axes[1, 0].plot(s_pts, h_pts.to('J/kg'),label = f"T_ambient = {dict['T_amb'].to('degC')}")
        axes[1, 0].set_xlabel('Entropy [J/kg/K]')
        axes[1, 0].set_ylabel('Enthalpy [J/kg]')
        axes[1, 0].set_title('H-S Diagram')

        # Top right: P vs T
        axes[1, 1].semilogx(v_pts.to('m^3/kg'), P_pts.to('atm'),label = f"T_ambient = {dict['T_amb'].to('degC')}")
        axes[1, 1].set_xlabel('Specific Volume [m^3/kg]')
        axes[1, 1].set_ylabel('Pressure [atm]')
        axes[1, 1].set_title('P-T Diagram')

    

    axes[0, 1].plot(T_ambs, cops, '-o' )
    axes[0, 1].set_xlabel('T_ambient [°C]')
    axes[0, 1].set_ylabel('COP [-]')
    axes[0, 1].set_title('COP vs T_ambient')

    axes[1,0].set_xlim([1200,2500])
    axes[1,0].set_ylim([300000,600000])
    axes[1,0].legend()
    plt.legend()
    plt.show()
    

    

    # Top left: T vs S


    # Bottom right: placeholder
    axes[0, 0].set_visible(False)

    plt.tight_layout()
    plt.show()

    return


     
    



if __name__ == '__main__':

    fluid = "R32"
    Q_dot = Q_(10,'kW')
    T_high1 = Q_(100,'degC')
    T_low1 = Q_(-10,'degC')
    T_high2 = Q_(100,'degC')
    T_low2 = Q_(0,'degC')
    T_low3 = Q_(10,'degC')
    T_high3 = Q_(100,'degC')
    T_low4 = Q_(20  ,'degC')
    T_high4 = Q_(100,'degC')
    states1, w_dot1, cop1 = HeatPumpCycle(fluid, Q_dot,T_low1,T_high1)
    states2, w_dot2, cop2 = HeatPumpCycle(fluid, Q_dot,T_low2,T_high2)
    states3, w_dot3, cop3 = HeatPumpCycle(fluid, Q_dot,T_low3,T_high3)
    states4, w_dot4, cop4 = HeatPumpCycle(fluid, Q_dot,T_low4,T_high4)
    dictlist = [states1,states2, states3, states4]

    
    plotloopandHSDiagram(fluid,dictlist)
    
    
