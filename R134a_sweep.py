#!/Library/Frameworks/Python.framework/Versions/3.10/bin/python3
import numpy as np
import cantera as ct
import matplotlib.pyplot as plt
plt.style.use('bmh')


def  R134aPropertysweep(T0, P0):
    rFluid = ct.Hfc134a()
    rFluid.TP = T0, P0
    

    h0 = rFluid.enthalpy_mass
    rho0 = rFluid.density_mass
    Pcr = rFluid.critical_pressure
    
    
    



if __name__ == '__main__':
    T0 = 300
    P0 = ct.one_atm
    R134aPropertysweep(T0,P0)

