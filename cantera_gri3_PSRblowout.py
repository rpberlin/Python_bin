#!/Library/Frameworks/Python.framework/Versions/3.10/bin/python3

import time
import numpy as np
import matplotlib.pyplot as plt
import cantera as ct
plt.style.use('bmh')


def run_PSR_BlowoutSweep(gas, t_resMin, residence_time):
    inlet = ct.Reservoir(gas)

    h0 = gas.enthalpy_mass
    gas.equilibrate('HP')
    combustor = ct.IdealGasReactor(gas)
    combustor.volume = 1.0


    # Create a reservoir for the exhaust
    exhaust = ct.Reservoir(gas)


    # A PressureController has a baseline mass flow rate matching the 'master'
    # MassFlowController, with an additional pressure-dependent term. By explicitly
    # including the upstream mass flow rate, the pressure is kept constant without
    # needing to use a large value for 'K', which can introduce undesired stiffness.


    # the simulation only contains one reactor
    sim = ct.ReactorNet([combustor])

    # Run a loop over decreasing residence times, until the reactor is extinguished,
    # saving the state after each iteration.
    states = ct.SolutionArray(gas, extra=['tres','mdot','mass','deltaH','h_sensible','s_sensible'])

    while residence_time > t_resMin:
        mdot = combustor.mass/residence_time
        inlet_mfc = ct.MassFlowController(inlet, combustor, mdot=mdot)
        outlet_mfc = ct.PressureController(combustor, exhaust, master=inlet_mfc, K=0.01)
        sim.set_initial_time(0.0)  # reset the integrator
        mdot = combustor.mass/residence_time
        sim.advance_to_steady_state()

        #print('tres = {:.2e}; mass = {:.3f} ;''T = {:.1f}'.format(residence_time,combustor.mass, combustor.T))
        h_sens,s_sens = sensible_HS(combustor.thermo)
        states.append(combustor.thermo.state, tres=residence_time, mdot=mdot, mass=combustor.mass,deltaH=combustor.thermo.enthalpy_mass-h0,h_sensible=h_sens,s_sensible=s_sens)
        residence_time *= 0.9  # decrease the residence time for the next iteration


    return states

def plot_PSR_BlowoutResults(states, species_to_track):
    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    #ax1.semilogx(states.tres, states.h_sensible*states.mdot, '.-',color='C0')
    #ax2 = ax1.twinx()
    ax1.semilogx(states.tres, states.T, '.-',color='C1')
    ax1.set_xlabel('residence time [s]')
    #ax1.set_ylabel('heat release rate per Mass Flow [W/m$^3$]', color='C0')
    ax1.set_ylabel('Temperature [K]')
    plt.title('PSR Temperature')
    

    ax3 = fig.add_subplot(122)
    for specie in species_to_track:
        idx = states.species_index(specie)
        ax3.loglog(states.tres, states.X[:,idx], '.-', label=specie)
    ax3.set_ylim([1e-6,1])
    ax3.set_xlabel('Residence Time (-)')
    ax3.set_ylabel('Mole Fraction (-)')
    plt.title('Species Profiles')
    plt.legend()

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
    time0 = time.perf_counter()
    gas = ct.Solution('gri30.yaml')
    t_resMax=10
    t_resMin=1e-6

    equiv_ratio = 0.9
    fuel = 'CH4:1.0 NH3:0.0 H2:.00'
    air = 'O2:1.0, N2:3.76'
    species_to_track = ['O2', 'H2', 'H2O', 'NH3', 'NO', 'NO2', 'CO', 'CO2']
    gas.TP = 298.15, ct.one_atm
    gas.set_equivalence_ratio(equiv_ratio, fuel = fuel, oxidizer = air)
    states = run_PSR_BlowoutSweep(gas,t_resMin, t_resMax)
    plot_PSR_BlowoutResults(states,species_to_track)
    time1 = time.perf_counter()
    print('Elapsed Time: ',time1-time0)
