import numpy as np
import matplotlib.pyplot as plt


def Tt_per_T(M,gam):
    return 1 + 0.5*(gam-1) * M **2

def Pt_per_P(M,gam):
    return Tt_per_T(M,gam) ** (gam/(gam-1))

def rhot_per_rho(M,gam):
    return Tt_per_T(M,gam) ** (gam/(gam-1))

def normshock_M2(M1,gam):
    return np.sqrt((1+0.5*(gam-1)*M1 **2)/(gam*M1**2 - 0.5*(gam-1)))

def normshock_P2_per_P1(M1,gam):
    return 1 + ((2*gam)/(gam+1))*(M1**2 - 1)

def normshock_rho2_per_rho1(M1,gam):
    return ((gam+1) * M1 ** 2)/((gam-1) * M1 **2 +2)

def normshock_T2_per_T1(M1,gam):
    return normshock_P2_per_P1(M1,gam)/normshock_rho2_per_rho1(M1,gam)

def normshock_To2_per_To1(M1,gam):
    return 1

def normshock_Po2_per_Po1(M1,gam):
    term1 = (((gam+1)*M1**2)/((gam-1)*M1**2 +2)) **(gam/(gam-1))
    term2 = ((gam+1)/(2*gam*M1**2 - (gam-1))) **(1/(gam-1))
    return term1*term2

def oblique_shock_theta(M1,beta,gam):
    return np.arctan((2/np.tan(beta))*(((M1*np.sin(beta))**2 -1)/(M1**2*(gam+np.cos(2*beta)+2))))

def A_per_Astar(M,gam):
    return np.sqrt((1/M**2)*((2/(gam+1)*(1+0.5*(gam-1)*M**2)))**((gam+1)/(gam-1)))





if __name__ == '__main__':
    M_sweep = np.linspace(0.1,3,20)
    Pt_sweep  = Pt_per_P(M_sweep,1.4)
    Tt_sweep  = Tt_per_T(M_sweep,1.4)
    M2_sweep  = normshock_M2(M_sweep,1.4)
    P2_sweep  = normshock_P2_per_P1(M_sweep,1.4)
    T2_sweep  = normshock_T2_per_T1(M_sweep,1.4)
    P02_sweep = normshock_Po2_per_Po1(M_sweep,1.4)
    A_per_Astar_sweep = A_per_Astar(M_sweep,1.4)
    #plt.plot(M_sweep,Pt_sweep)
    #plt.plot(M_sweep,Tt_sweep)
    #plt.plot(M_sweep,P02_sweep)
    plt.plot(M_sweep,A_per_Astar_sweep)
    plt.show()