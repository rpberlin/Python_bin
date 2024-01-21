#!/Library/Frameworks/Python.framework/Versions/3.10/bin/python3
import numpy as np
import cantera as ct
import matplotlib.pyplot as plt
plt.style.use('bmh')


def  gri30_equilibPhiSweep(phiPts, species_to_track):
    carbon = ct.Solution('graphite.yaml')
    g1 = ct.Solution('gri30.yaml')
    fuel = "CH4:.1,H2:0.1"
    g1.TP = 273, 101325
    mix_phases = [(g1, 1.0), (carbon, 0.0)]
    air = "O2:1.00,N2:0.1"


    #g1.set_equivalence_ratio(1.0, fuel="CH4:1", oxidizer="O2:0.233,N2:0.767", basis='mass')
    #g1.X = 'CH4:1, O2:2, N2:7.51'
    nSpecies = len(species_to_track)
    #Z_burnt = g1.mixture_fraction(fuel, air)
    T_flame = np.zeros(nPhiPts)
    X_table = np.zeros([nPhiPts, nSpecies])

    for i, phi in enumerate(phiPts):
        g1.set_equivalence_ratio(phi=phi, fuel=fuel, oxidizer=air)
        mix = ct.Mixture(mix_phases)
        mix.T = 300
        mix.P = 101325
        #g1.equilibrate('HP',solver='gibbs')
        mix.equilibrate('UV')
        T_flame[i] = mix.T
        X_table[i,:] = g1[species_to_track].Y

        #print('T_flame: ',T_flame)
        #print('X: ',g1[species_to_track].Y)
    #print(g1.species_names)
    return X_table, T_flame

def plotTandX(T,X,species_to_track):
    nPhiPts, nSpecies = X.shape
    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    line1 = plt.plot(phiPts,T-273,label='Temperature Rise in CH4/Air')
    ax1.set_xlabel('Equivalence Ratio (-)')
    ax1.set_ylabel('Temperature Delta ($^\circ$C)')
    plt.title('CV Ignition Temperature Rise')
    #plt.legend()

    ax2 = fig.add_subplot(122)
    for i, specie in enumerate(species_to_track):
        line2 = plt.semilogy(phiPts,X[:,i],label=specie)
    ax2.set_ylim([1e-6,1])
    ax2.set_xlabel('Equivalence Ratio (-)')
    ax2.set_ylabel('Mole Fraction (-)')
    plt.legend()
    plt.title('Equilibrium Mole Fractions')
    plt.show()

    fig.tight_layout()





if __name__ == '__main__':
    phi_low = 0.4
    phi_high = 1.8
    nPhiPts = 50
    phiPts = np.linspace(phi_low, phi_high, nPhiPts)
    species_to_track = ['O2', 'H2', 'CO', 'CO2']

    X_table, T_flame = gri30_equilibPhiSweep(phiPts,species_to_track)
    plotTandX(T_flame,X_table,species_to_track)
