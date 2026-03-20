import numpy as np
import matplotlib.pyplot as plt




def get_tau_given_pi_and_e(pi, gam, e):
    return pi ** ((gam-1)*e/gam)

def get_M_from_Pt_P_ratio(Pt_P_ratio,gam):
    return np.sqrt((Pt_P_ratio ** ((gam-1)/gam) - 1)*(2/(gam-1)))

def aspeed(gam,R,T):
    return np.sqrt(gam*R*T)

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


def MFP_Machterm(M,gam,R):
    return M*np.sqrt(gam/R)*np.power(1+0.5*(gam-1)* M **2,-1*(gam+1)/(2*gam -2))


def get_subsonic_Mach_from_mdot_A_Pt_Tt(mdot, A, Pt, Tt, gam, R):
    print(mdot,A,Pt,Tt, gam, R)
    LHS = mdot*np.sqrt(Tt)/(A*Pt)
    URF = 0.5
    eta = 0.001
    maxsteps = 20
    M0 = 0.75
    f_ref = MFP_Machterm(M0,gam,R)
    for i in range(0,maxsteps):
        f0 = MFP_Machterm(M0,gam,R)
        dfdM = (MFP_Machterm(M0+eta,gam,R) - f0)/eta
        err = LHS - f0
        if abs(err/f_ref) < eta:
            print('SUBSONIC DIFFUSER M CONVERGED')
            return M0 
        M1 = min( M0 + err/dfdM,1.0)
        #print(i, M0, f0, err, dfdM, M1)
        M0 = M0 * URF + (1-URF)*M1
    print('*** INCOMPLETE DIFFUSER MACH CONVERGENCE ***')
    return max(0,M0)

def getfanpi_for_given_Uout_Pout_from_Pt_Tt_in(Uy_target, Py, Ptx, Ttx, gam, R, e_poly):
    pi0 = 1.3
    URF = 0.5
    eta = 0.001
    maxsteps = 20
    
    for i in range(0,maxsteps):
        pi0_plus_eta = pi0+eta
        tau = get_tau_given_pi_and_e(pi0,gam,e_poly)
        tau_plus_eta = get_tau_given_pi_and_e(pi0_plus_eta,gam,e_poly)
        Tty = Ttx*tau
        Tty_plus_eta = Ttx*tau_plus_eta
        Pty = Ptx * pi0
        Pty_plus_eta = Ptx*pi0_plus_eta
        My = get_M_from_Pt_P_ratio(Pty/Py, gam)
        My_plus_eta = get_M_from_Pt_P_ratio(Pty_plus_eta/Py, gam)
        Ty = Tty/Tt_per_T(My,gam)
        Ty_plus_eta = Tty_plus_eta/Tt_per_T(My_plus_eta,gam)
        Uy = My*aspeed(gam,R,Ty)
        Uy_plus_eta = My_plus_eta*aspeed(gam,R,Ty_plus_eta)
        dUy_deta = (Uy_plus_eta - Uy)/eta
        err = Uy_target - Uy 
        if abs(err/Uy_target) < eta:   
            print('FAN PI ONVERGED')
            return pi0, tau, Pty, Tty, Ty 
        pi1 = pi0 + err/dUy_deta
        #print(i, Uy_target.to('m/s'), Uy.to('m/s'), pi0, pi1)
        pi0 = pi0 * URF + (1-URF)*pi1
    print(f'*** INCOMPLETE FAN PI CONVERGENCE *** Err/Utarget = {err/Uy_target:.2g}')
    return pi0, tau, Pty, Tty, Ty
        




if __name__ == '__main__':
    M_sweep = np.linspace(0.1,1.6,20)
    Pt_sweep  = Pt_per_P(M_sweep,1.4)
    Tt_sweep  = Tt_per_T(M_sweep,1.4)
    M2_sweep  = normshock_M2(M_sweep,1.4)
    P2_sweep  = normshock_P2_per_P1(M_sweep,1.4)
    T2_sweep  = normshock_T2_per_T1(M_sweep,1.4)
    P02_sweep = normshock_Po2_per_Po1(M_sweep,1.4)
    A_per_Astar_sweep = A_per_Astar(M_sweep,1.4)
    MFP_machtermsweep = MFP_Machterm(M_sweep,1.4,1)
    #plt.plot(M_sweep,Pt_sweep)
    #plt.plot(M_sweep,Tt_sweep)
    #plt.plot(M_sweep,P02_sweep)
    #plt.plot(M_sweep,A_per_Astar_sweep)
    #plt.show()
    #print(A_per_Astar(0.097,1.4))
    plt.plot(M_sweep,MFP_machtermsweep)
    plt.plot([0.5,0.76,1.0,1.5],[0.511,0.6478,.684,.582],'o')
    plt.show()
    print(MFP_Machterm(0.5,1.4,1))



