#!/Library/Frameworks/Python.framework/Versions/3.10/bin/python3


import cantera as ct
import numpy as np


import matplotlib.pyplot as plt
plt.style.use('bmh')


# reaction mechanism, kinetics type and compositions
def runLPtry1_Simulation(rpm, inlet_open_angle=0,inlet_open_angle_delta=180,outlet_open_angle=522,injector_open_angle=355,injector_open_angle_delta=50):

    f = rpm / 60.  # engine speed [1/s] (3000 rpm)
    V_H = .1e-3  # displaced volume [m**3]
    epsilon = 30.  # compression ratio [-]
    d_piston = V_H**.333  # piston diameter [m]
    V_oT = V_H / (epsilon - 1.)
    A_piston = .25 * np.pi * d_piston ** 2
    stroke = V_H / A_piston

    def piston_speed(t):
        return - stroke / 2 * 2 * np.pi * f * np.sin(crank_angle(t))


    def trapz(phi, t):
        tempSum = 0
        for i in range(1,len(t)):
           tempSum += 0.5*(phi[i-1]+phi[i])*(t[i]-t[i-1])
        return tempSum

    def crank_angle(t):
        """Convert time to crank angle"""
        return np.remainder(2 * np.pi * f * t, 4 * np.pi)

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

    #Specify reaction chemistry and composition of fuel + air
    reaction_mechanism = 'nDodecane_Reitz.yaml'
    phase_name = 'nDodecane_IG'
    comp_air = 'o2:1, n2:3.76'
    comp_fuel = 'c12h26:1'

    # turbocharger temperature, pressure, and composition
    T_inlet = 400.  # K
    p_inlet = 1.3e5  # Pa
    comp_inlet = comp_air

    # outlet pressure
    p_outlet = 1.3e5  # Pa

    # fuel properties (gaseous!)
    T_injector = 300.  # K
    p_injector = 1600e5  # Pa
    comp_injector = comp_fuel

    # ambient properties
    T_ambient = 300.  # K
    p_ambient = 1e5  # Pa
    comp_ambient = comp_air

    # Inlet valve friction coefficient, open and close timings
    inlet_valve_coeff = 1e-6
    inlet_open_rad = 0. / 180. * np.pi
    inlet_close_rad = (inlet_open_angle+inlet_open_angle_delta) / 180. * np.pi

    # Outlet valve friction coefficient, open and close timings
    outlet_valve_coeff = 1.e-7
    outlet_open_rad = outlet_open_angle / 180 * np.pi
    outlet_close_rad = 0. / 180. * np.pi

    # Fuel mass, injector open and close timings
    injector_open_rad = injector_open_angle/ 180. * np.pi
    injector_close_rad = (injector_open_angle + injector_open_angle_delta) / 180. * np.pi
    injector_mass = .5e-5  # kg

    # Simulation time and parameters
    sim_n_revolutions = 8
    delta_T_max = 20.
    rtol = 1.e-12
    atol = 1.e-16


    # load reaction mechanism
    gas = ct.Solution(reaction_mechanism, phase_name)
    #gas = ct.Solution(reaction_mechanism)
    # define initial state and set up reactor
    gas.TPX = T_inlet, p_inlet, comp_inlet
    rho_0 = gas.density
    cyl = ct.IdealGasReactor(gas)
    cyl.volume = V_oT

    # define inlet state
    gas.TPX = T_inlet, p_inlet, comp_inlet
    inlet = ct.Reservoir(gas)

    # inlet valve
    inlet_valve = ct.Valve(inlet, cyl)
    inlet_delta = np.mod(inlet_close_rad - inlet_open_rad, 4 * np.pi)
    inlet_valve.valve_coeff = inlet_valve_coeff
    inlet_valve.set_time_function(
        lambda t: np.mod(crank_angle(t) - inlet_open_rad, 4 * np.pi) < inlet_delta)

    # define injector state (gaseous!)
    gas.TPX = T_injector, p_injector, comp_injector
    injector = ct.Reservoir(gas)

    # injector is modeled as a mass flow controller
    injector_mfc = ct.MassFlowController(injector, cyl)
    injector_delta = np.mod(injector_close_rad - injector_open_rad, 4 * np.pi)
    injector_t_open = (injector_close_rad - injector_open_rad) / 2. / np.pi / f
    injector_mfc.mass_flow_coeff = injector_mass / injector_t_open
    injector_mfc.set_time_function(
        lambda t: np.mod(crank_angle(t) - injector_open_rad, 4 * np.pi) < injector_delta)

    # define outlet pressure (temperature and composition don't matter)
    gas.TPX = T_ambient, p_outlet, comp_ambient
    outlet = ct.Reservoir(gas)

    # outlet valve
    outlet_valve = ct.Valve(cyl, outlet)
    outlet_delta = np.mod(outlet_close_rad - outlet_open_rad, 4 * np.pi)
    outlet_valve.valve_coeff = outlet_valve_coeff
    outlet_valve.set_time_function(
        lambda t: np.mod(crank_angle(t) - outlet_open_rad, 4 * np.pi) < outlet_delta)

    # define ambient pressure (temperature and composition don't matter)
    gas.TPX = T_ambient, p_ambient, comp_ambient
    ambient_air = ct.Reservoir(gas)

    # piston is modeled as a moving wall
    piston = ct.Wall(ambient_air, cyl)
    piston.area = A_piston
    piston.set_velocity(piston_speed)

    # create a reactor network containing the cylinder and limit advance step
    sim = ct.ReactorNet([cyl])
    sim.rtol, sim.atol = rtol, atol
    cyl.set_advance_limit('temperature', delta_T_max)

    # simulate with a maximum resolution of 1 deg crank angle
    dt = 1. / (360 * f)
    t_stop = sim_n_revolutions / f

    # set up output data arrays
    states = ct.SolutionArray(
        cyl.thermo,
        extra=('t', 'ca', 'V', 'm', 'mdot_in', 'mdot_out', 'dWv_dt', 'injector_mdot','hs_sensible','ca_total'),
    )

    #####################################################################
    # Run Simulation
    #####################################################################
    step=0
    while sim.time < t_stop:
        step+=1
        if step%100 == 0:
            print('Step: ',step,' Crank Angle: ',crank_angle(sim.time),' time: ', sim.time)
        # perform time integration
        sim.advance(sim.time + dt)

        # calculate results to be stored
        dWv_dt = - (cyl.thermo.P - ambient_air.thermo.P) * A_piston * \
            piston_speed(sim.time)

        # append output data
        states.append(cyl.thermo.state,
                      t=sim.time, ca=crank_angle(sim.time),
                      V=cyl.volume, m=cyl.mass,
                      mdot_in=inlet_valve.mass_flow_rate,
                      mdot_out=outlet_valve.mass_flow_rate,
                      dWv_dt=dWv_dt,injector_mdot=injector_mfc.mass_flow_rate,hs_sensible=sensible_HS(gas),ca_total = f*360*sim.time)


    ######################################################################
    # Integral Results
    ######################################################################
    # heat release
    Q = trapz(states.heat_release_rate * states.V, states.t)
    output_str = '{:45s}{:>4.1f} {}'
    print(output_str.format('Heat release rate per cylinder (estimate):',
                            Q / states.t[-1] / 1000., 'kW'))

    # expansion power
    W = trapz(states.dWv_dt, states.t)
    print(output_str.format('Expansion power per cylinder (estimate):',
                            W / states.t[-1] / 1000., 'kW'))

    # efficiency
    eta = W / Q
    print(output_str.format('Thermal Efficiency:', eta * 100., '%'))

    #Indicated Mean Effective Pressure
    MEP = W/V_H
    print(output_str.format('Mean Effective Pressure: ', MEP/100000,'bar'))

    # CO emissions
    MW = states.mean_molecular_weight
    CO_emission = trapz(MW * states.mdot_out * states('CO').X[:, 0], states.t)
    CO_emission /= trapz(MW * states.mdot_out, states.t)
    print(output_str.format('CO emission (estimate):', CO_emission * 1.e6, 'ppm'))

    def ca_ticks(t):
        return np.round(crank_angle(t) * 180 / np.pi, decimals=1)-1

    xticks = np.arange(states.t[0],states.t[-1], 0.02)

    ticklabels=ca_ticks(xticks)



    return states, Q , W , MEP, eta , CO_emission ,xticks, ticklabels



def plotStates(states, ticks,ticklabels):
    t = states.t
    #xticks = np.arange(t[0], t[1], 0.02)

    # pressure and temperature
    fig = plt.figure(figsize=(16, 11))
    ax1 = fig.add_subplot(421)
    ax1.set_ylabel('Mass Injection Rate (kg/s)')
    ax1.set_xlabel(r'$\phi$ [deg]')
    ax1.plot(t, states.injector_mdot)
    ax1.set_xticks(ticks)
    ax1.set_xticklabels(ticklabels)

    ax1 = fig.add_subplot(422)
    ax1.plot(t, states.m)
    ax1.set_ylabel('Mass In Cylinder (-)')
    ax1.set_xlabel(r'$\phi$ [deg]')
    ax1.set_xticks(ticks)
    ax1.set_xticklabels(ticklabels)

    ax1 = fig.add_subplot(423)
    ax1.plot(t, states.P / 1.e5)
    ax1.set_ylabel('$Pressure$ [bar]')
    ax1.set_xlabel(r'$\phi$ [deg]')
    ax1.set_xticks(ticks)
    ax1.set_xticklabels(ticklabels)

    ax2 = fig.add_subplot(424)
    ax2.plot(states.t, states.T)
    ax2.set_ylabel('$T$ [K]')
    ax2.set_xlabel(r'$\phi$ [deg]')
    ax1.set_xticks(ticks)
    ax1.set_xticklabels(ticklabels)

    # p-V diagram
    ax = fig.add_subplot(425)
    ax.semilogy(states.V * 1000, states.P / 1.e5)
    ax.set_xlabel('$V$ [l]')
    ax.set_ylabel('$Pressure$ [bar]')

    # T-S diagram
    ax = fig.add_subplot(426)
    ax.plot(1.e-3*states.hs_sensible[:,1], 1.e-3*states.hs_sensible[:,0])
    ax.set_xlabel('$S$ [kJ/kg]')
    ax.set_ylabel('$H$ [kJ/kg-K]')

    # heat of reaction and expansion work
    ax = fig.add_subplot(427)
    ax.plot(t, 1.e-3 * states.heat_release_rate * states.V, label=r'$\dot{Q}$')
    ax.plot(t, 1.e-3 * states.dWv_dt, label=r'$\dot{W}_v$')
    #ax.set_ylim(-2e2, 1e3)
    ax.legend(loc=0)
    ax.set_ylabel('Heat and Work Rates [kW]')
    ax.set_xlabel(r'$\phi$ [deg]')
    ax1.set_xticks(ticks)
    ax1.set_xticklabels(ticklabels)


    # gas composition
    ax = fig.add_subplot(428)
    ax.plot(t, states('o2').X, label='O2')
    ax.plot(t, states('co2').X, label='CO2')
    ax.plot(t, states('co').X, label='CO')
    #ax.plot(t, states('CH4').X * 100, label='Fuel x100')
    ax.legend(loc=0)
    ax.set_ylabel('$X_i$ [-]')
    ax.set_xlabel(r'$\phi$ [deg]')
    ax1.set_xticks(ticks)
    ax1.set_xticklabels(ticklabels)

    fig.tight_layout()
    plt.show()






if __name__ == '__main__':
    rpm = 7000
    inlet_open_angle = 0
    inlet_open_delta = 180
    outlet_open_angle = 522
    injector_open_angle = 355
    injector_open_delta = 50
    states, Q , W , MEP, eta , CO_emission ,xticks, ticklabels = runLPtry1_Simulation(rpm, inlet_open_angle, inlet_open_delta,outlet_open_angle,injector_open_angle,injector_open_delta)
    plotStates(states,xticks,ticklabels)
    print('Simulation Complete')
