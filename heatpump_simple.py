"""
heat_pump_cycle.py
==================
Steady-state vapor-compression heat pump model using CoolProp for fluid
properties and Pint for unit-aware arithmetic.

Cycle overview (ideal, no superheat/subcooling losses):
    1  → 2   Isentropic compression (compressor)
    2  → 2a  Desuperheating (entry to condenser saturation curve)
    2a → 3   Condensation at constant pressure (condenser)
    3  → 4   Isenthalpic expansion (expansion valve)
    4  → 1   Evaporation at constant pressure (evaporator)
"""

import numpy as np
import CoolProp.CoolProp as cp
import pint
import matplotlib.pyplot as plt

ureg = pint.UnitRegistry()
Q_ = ureg.Quantity

plt.style.use('bmh')


# ---------------------------------------------------------------------------
# Utility: unit cleanup for state dictionaries
# ---------------------------------------------------------------------------

def dict_unit_cleanup(d):
    """
    Convert lists of Pint Quantities in a state dict into Quantity arrays.

    For each key whose values are a list of Quantities with the same unit,
    the list is collapsed into a single Quantity array (more efficient for
    downstream plotting).  Non-Quantity entries (floats, strings, …) are
    left as plain lists.
    """
    for key, values in d.items():
        if isinstance(values[0], Q_):
            units = values[0].units
            d[key] = Q_([v.m_as(units) for v in values], units)
        # else: leave as-is (floats, strings, etc.)
    return d


# ---------------------------------------------------------------------------
# CoolProp wrappers — all return Pint Quantities
# ---------------------------------------------------------------------------

def P_from_TQ(T, q, fluid):
    """Saturation pressure at temperature T and vapor quality q."""
    return Q_(cp.PropsSI('P', 'T', T.m_as('K'), 'Q', q, fluid), 'Pa')

def T_from_PQ(P, q, fluid):
    """Saturation temperature at pressure P and vapor quality q."""
    return Q_(cp.PropsSI('T', 'P', P.m_as('Pa'), 'Q', q, fluid), 'K')

def H_from_TQ(T, q, fluid):
    """Specific enthalpy at temperature T and vapor quality q."""
    return Q_(cp.PropsSI('H', 'T', T.m_as('K'), 'Q', q, fluid), 'J/kg')

def H_from_PQ(P, q, fluid):
    """Specific enthalpy at pressure P and vapor quality q."""
    return Q_(cp.PropsSI('H', 'P', P.m_as('Pa'), 'Q', q, fluid), 'J/kg')

def H_from_TS(T, s, fluid):
    """Specific enthalpy at temperature T and entropy s (superheated region)."""
    return Q_(cp.PropsSI('H', 'T', T.m_as('K'), 'S', s.m_as('J/kg/K'), fluid), 'J/kg')

def S_from_TQ(T, q, fluid):
    """Specific entropy at temperature T and vapor quality q."""
    return Q_(cp.PropsSI('S', 'T', T.m_as('K'), 'Q', q, fluid), 'J/kg/K')

def P_from_TS(T, s, fluid):
    """Pressure at temperature T and entropy s (superheated region)."""
    return Q_(cp.PropsSI('P', 'T', T.m_as('K'), 'S', s.m_as('J/kg/K'), fluid), 'Pa')

def V_from_TQ(T, q, fluid):
    """Specific volume at temperature T and vapor quality q."""
    return Q_(cp.PropsSI('V', 'T', T.m_as('K'), 'Q', q, fluid), 'm^3/kg')

def V_from_TS(T, s, fluid):
    """Specific volume at temperature T and entropy s (superheated region)."""
    return Q_(cp.PropsSI('V', 'T', T.m_as('K'), 'S', s.m_as('J/kg/K'), fluid), 'm^3/kg')

def Q_from_PH(P, h, fluid):
    """Vapor quality at pressure P and specific enthalpy h."""
    return cp.PropsSI('Q', 'P', P.m_as('Pa'), 'H', h.m_as('J/kg'), fluid)

def get_triple_TP(fluid):
    """Triple-point temperature and pressure for a given fluid."""
    return (Q_(cp.PropsSI('Ttriple', fluid), 'K'),
            Q_(cp.PropsSI('ptriple', fluid), 'Pa'))

def get_critical_TP(fluid):
    """Critical-point temperature and pressure for a given fluid."""
    return (Q_(cp.PropsSI('Tcrit', fluid), 'K'),
            Q_(cp.PropsSI('Pcrit', fluid), 'Pa'))


# ---------------------------------------------------------------------------
# Core cycle model
# ---------------------------------------------------------------------------

def HeatPumpCycle(fluid, Q_dot, T_low, T_high):
    """
    Compute an ideal vapor-compression heat pump cycle.

    Parameters
    ----------
    fluid  : str   CoolProp fluid name (e.g. 'R32')
    Q_dot  : Quantity [power]  Condenser heat output (heating capacity)
    T_low  : Quantity [temperature]  Evaporator saturation temperature
    T_high : Quantity [temperature]  Condenser saturation temperature

    Returns
    -------
    states : dict   Thermodynamic state arrays for the five cycle points
    w_dot  : Quantity [power]   Compressor power input
    cop    : Quantity [-]       Coefficient of performance (heating COP)

    State points
    ------------
    1   Compressor inlet  — saturated vapour at T_low
    2   Compressor outlet — superheated vapour after isentropic compression
    2a  Start of condensation — saturated vapour at condenser pressure
    3   Condenser outlet  — saturated liquid at T_high
    4   Evaporator inlet  — two-phase mixture after isenthalpic expansion
    """

    # --- State 1: compressor inlet (saturated vapour at evaporator temp) ---
    T1 = T_low
    P1 = P_from_TQ(T1, 1.0, fluid)
    s1 = S_from_TQ(T1, 1.0, fluid)
    h1 = H_from_TQ(T1, 1.0, fluid)
    v1 = V_from_TQ(T1, 1.0, fluid)

    # --- State 2: compressor outlet (isentropic → same entropy as state 1) ---
    T2 = T_high
    s2 = s1                               # isentropic process: s2 = s1
    h2 = H_from_TS(T2, s2, fluid)
    P2 = P_from_TS(T2, s2, fluid)
    v2 = V_from_TS(T2, s2, fluid)

    # --- State 2a: saturated vapour at condenser pressure ---
    # This marks where the superheated vapour meets the saturation curve and
    # condensation begins; used to close the loop on the H-S diagram.
    T2a = T_from_PQ(P2, 1.0, fluid)
    h2a = H_from_TQ(T2a, 1.0, fluid)
    P2a = P2
    s2a = S_from_TQ(T2a, 1.0, fluid)
    v2a = V_from_TQ(T2a, 1.0, fluid)

    # --- State 3: condenser outlet (saturated liquid at condenser pressure) ---
    P3 = P2
    T3 = T_from_PQ(P3, 0.0, fluid)
    h3 = H_from_TQ(T3, 0.0, fluid)
    s3 = S_from_TQ(T3, 0.0, fluid)
    v3 = V_from_TQ(T3, 0.0, fluid)

    # --- State 4: evaporator inlet (isenthalpic expansion: h4 = h3) ---
    P4 = P1
    h4 = h3                               # throttling: h conserved
    Q4 = Q_from_PH(P4, h4, fluid)        # vapour quality after expansion
    T4 = T_from_PQ(P4, Q4, fluid)
    s4 = S_from_TQ(T4, Q4, fluid)
    v4 = V_from_TQ(T4, Q4, fluid)

    # --- Performance metrics ---
    m_dot = Q_dot / (h2a - h3)           # mass flow rate from condenser duty
    w_dot = m_dot * (h2 - h1)            # compressor shaft power
    cop   = Q_dot / w_dot                 # heating COP

    print(f"  m_dot: {m_dot.to('kg/s'):.3f} | "
          f"w_dot: {w_dot.to('kW'):.3f} | "
          f"COP: {cop.to(''):.2f}")

    # --- Pack state points into a labelled dictionary ---
    states = {
        'T':      [T1, T2, T2a, T3, T4, T1],
        'P':      [P1, P2, P2a, P3, P4, P1],
        's':      [s1, s2, s2a, s3, s4, s1],
        'Q':      [1.0, None, 1.0, 0.0, Q4, 1.0],  # vapour quality at each point
        'h':      [h1, h2, h2a, h3, h4, h1],
        'v':      [v1, v2, v2a, v3, v4, v1],
        'labels': ['1', '2', '2a', '3', '4'],
        'cop':    [cop],
        'w_dot':  [w_dot],
        'T_amb':  [T1],
    }
    states = dict_unit_cleanup(states)
    return states, w_dot, cop


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plotloopandHSDiagram(fluid, dictlist):
    """
    Plot four diagnostic diagrams for a list of cycle operating conditions.

    Panels
    ------
    Top-left  : (hidden / reserved for future use)
    Top-right : COP vs. ambient (evaporator) temperature
    Bottom-left: H-S diagram with saturation dome and cycle loops
    Bottom-right: P-v diagram (log scale) with saturation dome and cycle loops
    """
    n_points = 300

    # --- Build saturation dome from triple point to critical point ---
    T_crit,  P_crit  = get_critical_TP(fluid)
    T_triple, P_triple = get_triple_TP(fluid)
    T_sweep = np.linspace(T_crit, T_triple, n_points)

    h_liq = H_from_TQ(T_sweep, 0,   fluid)
    h_vap = H_from_TQ(T_sweep, 1.0, fluid)
    s_liq = S_from_TQ(T_sweep, 0,   fluid)
    s_vap = S_from_TQ(T_sweep, 1.0, fluid)
    P_liq = P_from_TQ(T_sweep, 0,   fluid)
    P_vap = P_from_TQ(T_sweep, 1.0, fluid)
    v_liq = V_from_TQ(T_sweep, 0,   fluid)
    v_vap = V_from_TQ(T_sweep, 1.0, fluid)

    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))

    # Draw saturation dome on H-S and P-v panels
    axes[0, 0].plot(s_liq.to('J/kg/K'), T_sweep.to('degC'), 'r-', label='Sat. liquid')
    axes[0, 0].plot(s_vap.to('J/kg/K'), T_sweep.to('degC'), 'b-', label='Sat. vapour')
    axes[1, 0].plot(s_liq.to('J/kg/K'), h_liq.to('J/kg'), 'r-', label='Sat. liquid')
    axes[1, 0].plot(s_vap.to('J/kg/K'), h_vap.to('J/kg'), 'b-', label='Sat. vapour')
    axes[1, 1].semilogx(v_liq, P_liq.to('atm'), 'r-', label='Sat. liquid')
    axes[1, 1].semilogx(v_vap, P_vap.to('atm'), 'b-', label='Sat. vapour')

    cops   = []
    T_ambs = []

    for state in dictlist:
        T_amb_label = state['T_amb'].to('degC')
        cops.append(state['cop'])
        T_ambs.append(state['T_amb'])

        # T-S diagram: cycle loop for this operating condition
        axes[0, 0].plot(state['s'], state['T'].to('degC'),label=f"T_amb = {T_amb_label:.0f}")
        axes[0, 0].set_xlabel('Entropy [J/kg/K]')
        axes[0, 0].set_ylabel('Temperature [°C]')
        axes[0, 0].set_title('T-S Diagram')

        # H-S diagram: cycle loop for this operating condition
        axes[1, 0].plot(state['s'], state['h'].to('J/kg'),label=f"T_amb = {T_amb_label:.0f}")
        axes[1, 0].set_xlabel('Entropy [J/kg/K]')
        axes[1, 0].set_ylabel('Enthalpy [J/kg]')
        axes[1, 0].set_title('H-S Diagram')

        # P-v diagram: cycle loop for this operating condition
        axes[1, 1].semilogx(state['v'].to('m^3/kg'), state['P'].to('atm'),
                            label=f"T_amb = {T_amb_label:.0f}")
        axes[1, 1].set_xlabel('Specific Volume [m³/kg]')
        axes[1, 1].set_ylabel('Pressure [atm]')
        axes[1, 1].set_title('P-v Diagram')

    # COP vs ambient temperature summary plot
    axes[0, 1].plot(T_ambs, cops, '-o')
    axes[0, 1].set_xlabel('T_ambient [°C]')
    axes[0, 1].set_ylabel('COP [-]')
    axes[0, 1].set_title('COP vs T_ambient')

    # Zoom T-S diagram to the relevant region (wet/superheated zone for R32)
    axes[0, 0].set_xlim([1000, 2500])
    axes[0, 0].set_ylim([-20, 100])
    axes[0, 0].legend(fontsize=8)

    # Zoom H-S diagram to the relevant region (wet/superheated zone for R32)
    axes[1, 0].set_xlim([1200, 2500])
    axes[1, 0].set_ylim([300_000, 600_000])
    axes[1, 0].legend(fontsize=8)
    axes[1, 1].legend(fontsize=8)

    plt.tight_layout()
    plt.show()


# ---------------------------------------------------------------------------
# Entry point — sweep over four evaporator temperatures
# ---------------------------------------------------------------------------

if __name__ == '__main__':

    fluid = "R32"
    Q_dot  = Q_(10, 'kW')      # fixed condenser output (heating capacity)
    T_high = Q_(90, 'degC')   # fixed condenser temperature for all cases

    # Evaporator temperatures representing different ambient/source conditions
    T_lows = [Q_(-10, 'degC'), Q_(0, 'degC'), Q_(10, 'degC'), Q_(20, 'degC')]

    print("Cycle performance summary:")
    dictlist = []
    for T_low in T_lows:
        print(f"\n  T_evap = {T_low}:")
        states, w_dot, cop = HeatPumpCycle(fluid, Q_dot, T_low, T_high)
        dictlist.append(states)

    plotloopandHSDiagram(fluid, dictlist)