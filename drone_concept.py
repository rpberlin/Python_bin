import pint
ureg = pint.UnitRegistry()
import compressible_flow as cf
import numpy as np
import matplotlib.pyplot as plt
Q_ = ureg.Quantity

#station_0 = freestream
#station_1 = nacelle throat
#station_2 = fan/comporessor inlet
#station_13 = fan outlet
#station_19 = fan nozzle outlet


M0 = 0.7
M2_theta = 0.2
C_drag = 0.03
T0 = Q_(300,'K')
P0 = Q_(1,'atm')
rho0 = Q_(1.18,'kg/m^3')
rhoPLA = Q_(1000,'kg/m^3')
grav = Q_(9.81,'m/s^2')
gamma_air = 1.4
MW_air = Q_(28.97,'g/mol')
R_u = Q_(8.314,'J/mol/K')
R_air = R_u/MW_air
cp_air = gamma_air*R_air/(gamma_air-1)


a0 = np.sqrt(gamma_air*R_air*T0)

Pt0 = P0*cf.Pt_per_P(M0,gamma_air)
Tt0 = T0*cf.Tt_per_T(M0,gamma_air)

U0 = a0*M0
nacelle_inlet_D_i = Q_(7.1,'cm')
nacelle_inlet_D_o = Q_(8.9,'cm')
D_o_2 = Q_(10,'cm')
D_i_2 = Q_(6,'cm')
nacelle_OD = Q_(11,'cm')
stagnation_streamline_A = 0.25*np.pi * (nacelle_inlet_D_o **2 - nacelle_inlet_D_i **2)
stagnation_streamline_D = np.sqrt(4*stagnation_streamline_A/np.pi)

m_dot_in = rho0 * U0 * stagnation_streamline_A

nacelle_frontal_A = 0.25 * np.pi *nacelle_OD**2
nacelle_diffuser_half_hangle = Q_(4,'deg')
nacelle_diffuser_distance = Q_(8,'cm')
F_drag = 0.5 * rho0 * C_drag * (U0 ** 2) * nacelle_frontal_A
delta_U = F_drag/m_dot_in
U19 = U0 + delta_U

A2_star = stagnation_streamline_A/cf.A_per_Astar(M0,gamma_air)
A_2 = 0.25 * np.pi * (D_o_2 ** 2 - D_i_2 ** 2)
print(stagnation_streamline_A,A2_star,A_2)
if A_2 < A2_star:
    print('A2 TOO SMALL!')


P19 = P0
Pt2 = Pt0
Tt2 = Tt0
pi_c, tau_c, Pt19, Tt19, T19 =  cf.getfanpi_for_given_Uout_Pout_from_Pt_Tt_in(U19, P19, Pt2, Tt2, gamma_air, R_air, 0.9)


w_dot = m_dot_in*cp_air*Tt2*(tau_c-1) 

r_meanline =  0.5*(D_i_2+D_o_2)

M2 = cf.get_subsonic_Mach_from_mdot_A_Pt_Tt(m_dot_in, A_2, Pt2, Tt0, gamma_air, R_air)
T2 = Tt0/cf.Tt_per_T(M2,gamma_air)
a2 = cf.aspeed(gamma_air,R_air,T2)
U2_axial = M2*a2
U2_theta_SF = 0*U2_axial
U2_theta_RF = U2_theta_SF + M2_theta*a2
beta_in = np.atan(U2_theta_RF/U2_axial)
omega = U2_theta_RF/r_meanline
U3_theta_SF = ((w_dot/(m_dot_in*omega)) - r_meanline*U2_theta_SF)/r_meanline
U3_theta_RF = -1.0*U3_theta_SF + omega*r_meanline
U3_axial = U2_axial
beta_out = np.atan(U3_theta_RF/U3_axial)
swirl_angle_out = np.atan(U3_theta_SF/U3_axial)

rho19 = P19.to('Pa')/(R_air.to('J/kg/K')*T19)
A_19 = m_dot_in/(U19*rho19)
D_i_19 = nacelle_inlet_D_o
D_o_19 = np.sqrt( 4*A_19/np.pi + D_i_19 **2)

acc_centrif_mean = r_meanline * (omega.to('rad/s') **2) 
g_outer = acc_centrif_mean/grav
h_blade = 0.5*(D_o_2-D_i_2)
sigma_root = rhoPLA * h_blade * acc_centrif_mean
sigma_yield = Q_(70,'MPa') 
yield_fraction = sigma_root/sigma_yield


print(f"Pi_c = {pi_c:.3f} M2 = {M2:.3f} W_dot= {w_dot.to('W'):.1f}")
print(f"F_drag = {F_drag.to('N'):.1f} : {(F_drag/grav).to('g')}")
print(f"Omega = {omega.to('rad/s'):.0f} beta_in: {beta_in.to('deg'):.3f} beta_out: {beta_out.to('deg'):.3f} swirl_angle_out: {swirl_angle_out.to('deg'):.3f}")
print(f"U2_theta = {U2_theta_RF.to('m/s'):.0f} {U2_theta_SF.to('m/s'):.3f} U3_theta = {U3_theta_RF.to('m/s'):.0f} {U3_theta_SF.to('m/s'):.3f} ")
print(f"D_i_19 = {D_i_19.to('cm'):.2f} D_o_19 = {D_o_19.to('cm'):.2f}")
print(f"Centrifulag Acc : {acc_centrif_mean.to('m/s^2'):.1f} :: {g_outer.to(''):.2f}")
print(f"Stress/Yield Ratio: {yield_fraction.to(''):.4f}")
print(f"Pt0: {Pt0.to('Pa'):.1f} Tto0: {Tt0.to('K'):.1f}")
print('\n\n\n done')















#plt.plot(M_cruise,F_drag.to('N'))
#plt.plot(M_cruise,m_dot_in.to('kg/s'))
#plt.plot(M_cruise,U_out.to('m/s'))




