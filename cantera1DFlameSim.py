#!/Library/Frameworks/Python.framework/Versions/3.10/bin/python3
import numpy as np
import cantera as ct
import matplotlib.pyplot as plt
import time
plt.style.use('bmh')


def run1D_MethaneFlameSim(Tin, Pop, Phi, species_to_track):

    Tin = 300.0  # unburned gas temperature [K]
    reactants = 'H2:1.1, O2:2, N2:5, H2O:.01, NH3:.04, CH4:.5'  # premixed gas composition
    width = 0.03  # m
    loglevel = 0  # amount of diagnostic output (0 to 8)
    
    # Solution object used to compute mixture properties, set to the state of the
    # upstream fuel-air mixture
    gas = ct.Solution('gri30.yaml')
    #gas = ct.Solution('h2o2.yaml')
    gas.TPX = Tin, Pop, reactants 
    
    f = ct.FreeFlame(gas, width=width)
    f.set_refine_criteria(ratio=3, slope=0.06, curve=0.12)
    f.show_solution()
    
    f.transport_model = 'Mix'
    f.solve(loglevel=loglevel, auto=True)
    
    Zlocs = f.grid
    nZlocs = len(Zlocs)
    nSpecies = len(species_to_track)
    Xtable = np.zeros([nZlocs, nSpecies])
    
    for i, specie in enumerate(species_to_track):
        idx = gas.species_index(specie)
        tmp = f.X[idx,:]
        Xtable[:,i] = tmp
                
    

    
    return f.grid, Xtable, f.T , f.velocity[0]

def plot1DflameProps(Zlocs,X,T,species_to_track,S_lam):
    nZpts, nSpecies = X.shape
    zoomWindow = 100 #Larger window shows more data points
    minMoleFraction = 1e-7 #threshold for plotting mole fractions
    diff = [0]
    minDiff = 1e6
    ZminDiff = 0
    for i in range (1,nZpts):
        iDiff = Zlocs[i]-Zlocs[i-1]
        diff.append(iDiff)
        if iDiff < minDiff:
            minDiff = iDiff
            ZminDiff = 0.5*(Zlocs[i]+Zlocs[i-1])
    diff[0]=2*diff[1]-diff[2]

    zoomThreshold = zoomWindow*np.min(diff)

    
    iMin = 0
    iMax = nZpts-1
    for i in range(1,nZpts):
        if diff[i] < zoomThreshold and diff[i-1] >= zoomThreshold:
            iMin = i
        elif diff[i] >= zoomThreshold and diff[i-1] < zoomThreshold:
            iMax = i
    
    
    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    line1 = plt.plot(Zlocs[iMin:iMax]-ZminDiff,T[iMin:iMax],label='Temperature')
    ax1.set_xlabel('Location (-)')
    ax1.set_ylabel('Temperature (K)')
    plt.title(f'Flame Speed = {S_lam:.3f} m/s')
    plt.legend()
            
    
    ax2 = fig.add_subplot(122)
    for j, specie in enumerate(species_to_track):
        line2 = plt.semilogy(Zlocs[iMin:iMax]-ZminDiff,X[iMin:iMax,j],label=specie)
    ax2.set_ylim([minMoleFraction,1])
    ax2.set_xlabel('Equivalence Ratio (-)')
    ax2.set_ylabel('Mole Fraction (-)')
    plt.legend()
    plt.title('Equilibrium Mole Fractions')
  
    
    plt.show()
    
    
    
if __name__ == '__main__':
    time0 = time.perf_counter()
    Phi = 1.0
    Tinlet = 300
    Pop = ct.one_atm
    species_to_track = ['O2', 'H2', 'H2O', 'NH3', 'NO', 'NO2', 'CO', 'CO2']
    
    Zlocs, X_table, T_flame, S_lam = run1D_MethaneFlameSim(Tinlet, Pop, Phi, species_to_track)
    plot1DflameProps(Zlocs, X_table,T_flame,species_to_track, S_lam)
    time1 = time.perf_counter()
    print('Elapsed Time: ',time1-time0)