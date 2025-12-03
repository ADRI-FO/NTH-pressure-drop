#===IMPORT PACKAGES===
from pyXSteam.XSteam import XSteam
steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS) # m/kg/sec/°C/bar/W
import matplotlib.pyplot as plt
import math as math

from iapws import IAPWS97
from CoolProp.CoolProp import PropsSI
from scipy.optimize import root
from scipy.optimize import fsolve
from scipy.integrate import quad
import numpy as np
# import pandas as pd
# import six
# import matplotlib
# import scipy.integrate as integrate
# from scipy.optimize import minimize, differential_evolution

#===INPUT PARAMETERS===
p_in = 55 # [bar]
T_in = steamTable.tsat_p(p_in) - 40 # [°C]
L_heated = 3.1 # [m]
D_rod = 10.3e-3 # [m]
Pitch = 21.2e-3 # [m]
relative_wall_roughness = 1e-3
heat_flux_pin_avg = 2.4e6 # [W/m²]
g = 9.80665 # [m/s²]

#===First exercise: Calculate friction monophase pressure drop through channel without heating===
#=== For a mass flow rate ranging from 0 to 2.5 kg/s ===
rho_subcooled = steamTable.rho_pt(p_in, T_in) # [kg/m³]
n_mass_flow_rate = 50
mass_flow_rate_list = [i * 2.5 / n_mass_flow_rate for i in range(1, n_mass_flow_rate+1)] # [kg/s]
dp_friction_list_1 = []
#===A : First for a pipe diameter equal to 2e-2 [m] and with McAdams correlation
#=== f = 0,184 * Re_C^(-0,2) if 3e4 < Re_C < 1e6 (turbulent) ===
D_pipe_1 = 2e-2 # [m]
A_pipe = math.pi * (D_pipe_1/2)**2 # [m²]

for mass_flow_rate in mass_flow_rate_list:
    u_1 = mass_flow_rate/(rho_subcooled * A_pipe)
    Re_1 = rho_subcooled * u_1 * D_pipe_1 / steamTable.my_pt(p_in, T_in)
    if Re_1 < 30000:
        f_1 = 64 / Re_1  # Laminar flow
    if Re_1 > 30000:
        f_1 = 0.184 * Re_1**(-0.2)
        # if Re_1 > 1e6:
        #     print('Reynolds number out of range for McAdams correlation: Re_C = ', Re_1)
    dp_friction = f_1 * (L_heated / D_pipe_1) * (rho_subcooled * u_1**2 / 2) / 1e5 # [bar]
    dp_friction_list_1.append(dp_friction)
# #===Plot the friction pressure drop===
# plt.plot(mass_flow_rate_list, dp_friction_list_1, label='D_pipe_1 = 2e-2 m, McAdams correlation')
# plt.xlabel('Mass Flow Rate [kg/s]')
# plt.ylabel('Friction Pressure Drop [bar]')
# plt.title('Friction Pressure Drop vs Mass Flow Rate')
# plt.legend()
# plt.grid()
# plt.show()

#===B : Then for the actual BWR subchannel still with McAdams correlation===
A_flow = (Pitch**2 - math.pi * (D_rod/2)**2) # [m²]
Perimeter_flow = math.pi * D_rod # [m]
Hydraulic_diameter = 4 * A_flow / Perimeter_flow # [m] 
dp_friction_list_2 = []
for mass_flow_rate in mass_flow_rate_list:
    u_2 = mass_flow_rate/(rho_subcooled * A_flow)
    Re_2 = rho_subcooled * u_2 * Hydraulic_diameter / steamTable.my_pt(p_in, T_in)
    if Re_2 < 30000:
        f_2 = 64 / Re_2  # Laminar flow
    if Re_2 > 30000:
        f_2 = 0.184 * Re_2**(-0.2)
        # if Re_2 > 1e6:
        #     print('Reynolds number out of range for McAdams correlation: Re_C = ', Re_2)
    dp_friction = f_2 * (L_heated / Hydraulic_diameter) * (rho_subcooled * u_2**2 / 2) / 1e5 # [bar]
    dp_friction_list_2.append(dp_friction)
# #===Plot the friction pressure drop===
# plt.plot(mass_flow_rate_list, dp_friction_list_2, label='BWR subchannel, McAdams correlation')
# plt.xlabel('Mass Flow Rate [kg/s]')
# plt.ylabel('Friction Pressure Drop [bar]')
# plt.title('Friction Pressure Drop vs Mass Flow Rate')
# plt.legend()
# plt.grid()
# plt.show()

#===C : Finally for the actual BWR subchannel with the Colebrook-White equation===

def colebrook_root(Re, rel_rough):
    def eq(f):
        f = f[0]
        return [
            1.0/math.sqrt(f) + 
            2.0*math.log10(rel_rough/3.7 + 2.51/(Re*math.sqrt(f)))
        ]

    sol = root(eq, [0.02])
    return sol.x[0]

dp_friction_list_3 = []
for mass_flow_rate in mass_flow_rate_list:
    u_C = mass_flow_rate / (rho_subcooled * A_flow)
    Re_C = rho_subcooled * u_C * Hydraulic_diameter / steamTable.my_pt(p_in, T_in)
    f = 0.184 * Re_C**(-0.2)   # Initial guess
    rel_rough = relative_wall_roughness   # already ε/D
    for _ in range(50):
        f_old = f
        f = 1.0 / ( -2.0 * math.log10(rel_rough/3.7 + 2.51/(Re_C*math.sqrt(f_old))) )**2
        if abs(f - f_old) < 1e-8:
            break
    friction_factor_C = f
    # Alternatively, use the root-finding function:
    #friction_factor_C = colebrook_root(Re_C, relative_wall_roughness)
    dp_C = friction_factor_C * (L_heated / Hydraulic_diameter) * (rho_subcooled * u_C**2 / 2) / 1e5
    dp_friction_list_3.append(dp_C)
#===Plot the friction pressure drop===
# plt.plot(mass_flow_rate_list, dp_friction_list_3, label='BWR subchannel, Colebrook-White equation')
# plt.xlabel('Mass Flow Rate [kg/s]')
# plt.ylabel('Friction Pressure Drop [bar]')
# plt.title('Friction Pressure Drop vs Mass Flow Rate')
# plt.legend()
# plt.grid()
# plt.show()


#===Second exercise: Calculate two phase pressure drop through the channel ===
#=== D : First assume the inlet equilibrium quality is 0.15 and the channel is adiabatic (no heating) ===
x_in_D = 0.15
rho_L_D = steamTable.rhoL_p(p_in) # [kg/m³]
rho_G_D = steamTable.rhoV_p(p_in) # [kg/m³]
void_fraction_D = 1/(1 + ((1 - x_in_D)/x_in_D) * (rho_G_D/rho_L_D)) # [-]
rho_m_D = void_fraction_D * rho_G_D + (1 - void_fraction_D) * rho_L_D # [kg/m³]
# print(f'Void fraction at inlet quality of 0.15: {void_fraction_D}')
# print(f'Mixture density at inlet quality of 0.15: {rho_m_D} kg/m³')
dp_acc__D = 0 # No acceleration in adiabatic case without change in quality
friction_factor_L0_list_D = []
dp_grav_D = rho_m_D * 9.81 * L_heated / 1e5 # [bar]
dp_tot_list_D = []
dp_fric_list_D = []
for mass_flow_rate in mass_flow_rate_list:
    G_m_D = mass_flow_rate / A_flow # [kg/m²/s]
    mu_L_D = PropsSI('VISCOSITY', 'P', p_in * 1e5, 'Q', 0, 'Water') # [Pa.s]
    mu_G_D = PropsSI('VISCOSITY', 'P', p_in * 1e5, 'Q', 1, 'Water') # [Pa.s]
    mu_ratio_D = 1/(1 + x_in_D * (mu_L_D/mu_G_D - 1)) # [-]
    u_L0_D = mass_flow_rate/(rho_L_D * A_flow) # [m/s]   # Liquid only velocity at inlet conditions (all the mass flow rate is considered to be liquid)
    Re_L0_D = rho_L_D * u_L0_D * Hydraulic_diameter / mu_L_D
    if Re_L0_D < 30000:
        friction_factor_L0_D = 64 / Re_L0_D  # Laminar flow
    if Re_L0_D > 30000:
        friction_factor_L0_D = 0.184 * Re_L0_D**(-0.2)
        # if Re_L0_D > 1e6:
        #     print('Reynolds number out of range for McAdams correlation: Re_C = ', Re_L0_D)
    friction_factor_L0_list_D.append(friction_factor_L0_D)
    friction_factor_TP_D = friction_factor_L0_D * (mu_ratio_D)**0.2
    dp_fric = (friction_factor_TP_D * (L_heated / Hydraulic_diameter) * (G_m_D**2 / (2 * rho_m_D))) / 1e5 # [bar]
    dp_fric_list_D.append(dp_fric)
    dp_tot = dp_fric + dp_grav_D + dp_acc__D
    dp_tot_list_D.append(dp_tot)
# #===Plot friction factor for liquid only at inlet conditions===
# plt.plot(mass_flow_rate_list, friction_factor_L0_list_D, label='Friction factor for liquid only at inlet conditions')
# plt.xlabel('Mass Flow Rate [kg/s]')
# plt.ylabel('Friction Factor [-]')
# plt.ylim(0, 0.03)
# plt.title('Friction Factor vs Mass Flow Rate')
# plt.legend()
# plt.grid()
# plt.show()
# #===Plot the two phase pressure drop===
# plt.plot(mass_flow_rate_list, dp_tot_list_D, label='Inlet quality = 0.15, Adiabatic channel')
# plt.xlabel('Mass Flow Rate [kg/s]')
# plt.ylabel('Two Phase Pressure Drop [bar]')
# plt.title('Two Phase Pressure Drop vs Mass Flow Rate')
# plt.legend()
# plt.grid()
# plt.show()


def rho_m_plus(rho_L, rho_G, x):
    return 1 / (x / rho_G + (1 - x) / rho_L)


#===E : Then assume the inlet and outlet equilibrium quality are 0 and 0.15 respectively.===
#=== Thus instead of using a heat flux as input it is assumed that the quality evolves linearly from inlet to outlet===
x_in_E = 0.0
x_out_E = 0.15
x_E = lambda z: x_in_E + (x_out_E - x_in_E) * (z / L_heated) # [-], linear evolution of quality along the heated length
rho_L_E = steamTable.rhoL_p(p_in) # [kg/m³]
rho_G_E = steamTable.rhoV_p(p_in) # [kg/m³]

dp_tot_list_E = []
dp_acc_list_E = []
dp_fric_list_E = []
dp_grav_list_E = []
# Not a good idea to descretize the lenght. Better to integrate directly over the lenght.
for mass_flow_rate in mass_flow_rate_list:
    void_fraction_E = lambda z: 1/(1 + ((1 - x_E(z))/x_E(z)) * (rho_G_E/rho_L_E)) if x_E(z) != 0 else 0.0  # [-]
    rho_m_E = lambda z: void_fraction_E(z) * rho_G_E + (1 - void_fraction_E(z)) * rho_L_E  # [kg/m³]
    G_m_E = mass_flow_rate / A_flow # [kg/m²/s]
    mu_L_E = PropsSI('VISCOSITY', 'P', p_in * 1e5, 'Q', 0, 'Water') # [Pa.s]
    mu_G_E = PropsSI('VISCOSITY', 'P', p_in * 1e5, 'Q', 1, 'Water') # [Pa.s]
    mu_ratio_E = lambda z: 1/(1 + x_E(z) * (mu_L_E/mu_G_E - 1)) # [-]
    u_L0_E = mass_flow_rate /(rho_L_E * A_flow ) # [m/s]
    Re_L0_E = rho_L_E * u_L0_E * Hydraulic_diameter / mu_L_E
    if Re_L0_E < 30000:
        friction_factor_L0_E = 64 / Re_L0_E  # Laminar flow
    if Re_L0_E > 30000:
        friction_factor_L0_E = 0.184 * Re_L0_E**(-0.2)
        # if Re_L0_E > 1e6:
        #     print('Reynolds number out of range for McAdams correlation: Re_C = ', Re_L0_E)
    #===Integration over the heated length===
    dp_grav_E_dz = lambda z: 9.81 * (void_fraction_E(z) * rho_G_E + (1 - void_fraction_E(z)) *rho_L_E) # [Pa/m]
    dp_grav_E, _  = quad(dp_grav_E_dz, 0, L_heated)
    dp_grav_list_E.append(dp_grav_E / 1e5) # [bar]

    phi_2_L0_E = lambda z: mu_ratio_E(z)**0.2 * (rho_L_E / rho_m_E(z))
    dp_fric_E_dz = lambda z: (friction_factor_L0_E * (G_m_E**2) / (Hydraulic_diameter * 2 * rho_L_E)) * phi_2_L0_E(z) # [Pa/m]
    dp_fric_E, _ = quad(dp_fric_E_dz, 0, L_heated)
    dp_fric_list_E.append(dp_fric_E / 1e5) # [bar]

    dp_acc_E = G_m_E**2 * (1/rho_m_plus(rho_L_E, rho_G_E, x_out_E) - 1/rho_m_plus(rho_L_E, rho_G_E, x_in_E)) # [Pa]
    dp_acc_list_E.append(dp_acc_E / 1e5) # [bar]

    # print(f'Mass Flow Rate: {mass_flow_rate} kg/s, dp_grav: {dp_grav_E/1e5} bar, dp_acc: {dp_acc_E/1e5} bar, dp_fric: {dp_fric_E/1e5} bar}')

    dp_tot_E = (dp_grav_E + dp_acc_E + dp_fric_E) / 1e5 # [bar]
    dp_tot_list_E.append(dp_tot_E)
        
# #===Plot the two phase pressure drop===
# plt.plot(mass_flow_rate_list, dp_tot_list_E, label='Inlet quality = 0.0, Outlet quality = 0.15')
# plt.xlabel('Mass Flow Rate [kg/s]')
# plt.ylabel('Two Phase Pressure Drop [bar]')
# plt.title('Two Phase Pressure Drop vs Mass Flow Rate')
# plt.legend()
# plt.grid()
# plt.show()


#===Third exercise: Axial heat flux distribution===
#===F: Find axial position where the equilibrium quality is zero for a mass flow rate of 0,1kg/s===
#===This bulk boiling point (zB). Assume cosine axial heat distribution neglecting pressure drop for material properties===

#===Cosinusoïdal axial heat flux distribution===
#=== z = 0 at the middle of the heated length ===
#=== Integral of cos(z) from -pi/2 to pi/2 = 2 ===
# Therefore, maximum heat flux is the average heat flux  divided by 2/pi
mass_flow_rate_F = 0.1 # [kg/s]
heat_flux_pin_max = heat_flux_pin_avg / (2 / math.pi) # [W/m²]
def heat_flux_pin(z):
    z = float(np.atleast_1d(z)[0])
    return heat_flux_pin_max * np.cos(np.pi * z / L_heated)


# #=== Graphical representation of the axial heat flux distribution ===
# n_axial = 100
# z_values = [i * 1/n_axial for i in range(-int((L_heated/2)*n_axial), int((L_heated/2)*n_axial))] # [m]
# heat_flux_values = [heat_flux_pin(z) for z in z_values]
# plt.plot(heat_flux_values, z_values)
# plt.ylabel('Axial Position [m]')
# plt.xlabel('Heat Flux [W/m²]')
# plt.title('Axial Heat Flux Distribution')
# plt.grid()
# plt.show()

# #===Verification of the average heat flux calculation===
# integral_heat_flux = sum(heat_flux_values) * 1/n_axial # [W/m²]
# average_heat_flux = integral_heat_flux / L_heated # [W/m²]
# print(f'Calculated Average Heat Flux: {average_heat_flux} W/m²')


def h_m_z(z, mass_flow_rate):
    #===Calculate inlet properties===
    h_in = steamTable.h_pt(p_in, T_in) # [kJ/kg]
    
    #===Calculate enthalpy at position z===
    q_dz = lambda z: heat_flux_pin(z) * (math.pi * D_rod ) # [W]
    h_m_z = h_in + (quad(q_dz, -L_heated/2, z)[0]/ mass_flow_rate)  / 1000 # [kJ/kg]
    return h_m_z

def x_e_z(z, mass_flow_rate, p):
    #===Calculate dynamic quality at position z===
    h_l_sat = steamTable.hL_p(p) # [kJ/kg]
    h_g_sat = steamTable.hV_p(p) # [kJ/kg]
    x_e = (h_m_z(z, mass_flow_rate) - h_l_sat) / (h_g_sat - h_l_sat) # [-]
    return x_e


#===Find axial position of bulk boiling point (zB)===
def bulk_boiling_point_eq(z, mass_flow_rate):
    #===Check if the temperature at position z is equal to the saturation temperature===
    h_l_sat = steamTable.hL_p(p_in) # [kJ/kg]
    return h_m_z(z, mass_flow_rate) - h_l_sat

zB_solution = fsolve(bulk_boiling_point_eq, 0.01, args=(mass_flow_rate_F,)) # Initial guess at the middle of the heated length
zB = zB_solution[0]
#print(f'Axial position of bulk boiling point from the channel inlet (zB): {zB + (L_heated/2)} m')


#===Find axial position where the flow (dynamic quality) is one. This is were we have saturated vapour (zV)===
def saturated_vapour_point_eq(z, mass_flow_rate, p):
    #===Calculate dynamic quality at position z===
    return x_e_z(z, mass_flow_rate, p) - 1.0

zV_solution = fsolve(saturated_vapour_point_eq, 0.02, args=(mass_flow_rate_F, p_in)) # Initial guess at the middle of the heated length
zV = zV_solution[0]
#print(f'Axial position of saturated vapour point from the channel inlet (zV): {zV + (L_heated/2)} m')


#===Find axial position where the flow (dynamic quality) is zero: bubble detachment (zD) (start of subcool boiling)===
#Using Saha and Zuber correlation or Stanton for bubble detachment (low flow rates), and also do you have to conscider the local or the average heat flux? (check in the book)
def bubble_detachment_point_eq(z, mass_flow_rate):

    G_m = mass_flow_rate_F / A_flow # [kg/m²/s]
    Peclet = G_m * Hydraulic_diameter * steamTable.CpL_p(p_in) * 1e3 / steamTable.tcL_p(p_in) # [-] Or insted using the value at satutation, use the local one but need an iterative process
    # print(f'Peclet number at z = {z} m: {Peclet}')
    if Peclet < 7e4:
        #Use Nusselt correlation
        Nusselt_departure = 455
        T_bulk = steamTable.tsat_p(p_in) - (heat_flux_pin(z) * Hydraulic_diameter) / (Nusselt_departure * steamTable.tcL_p(p_in))
    else:
        #Use Stanton correlation
        Stanton_departure = 6.5e-3
        T_bulk = steamTable.tsat_p(p_in) - (heat_flux_pin(z)) / (Stanton_departure * G_m * steamTable.CpL_p(p_in) * 1e3)
    
    #===Find the z where the enthalpy gives a temperature equal to the bulk temperature calculated above===
    h_bulk = steamTable.h_pt(p_in, T_bulk) # [kJ/kg]
    return h_m_z(z, mass_flow_rate) - h_bulk

zD_solution = fsolve(bubble_detachment_point_eq, -1, args=(mass_flow_rate_F,)) # Initial guess at the middle of the heated length
zD = zD_solution[0]
# print(f'value of zD (from middle of heated length): {zD} m')
# print(f'Axial position of bubble detachment point from the channel inlet (zD): {zD + (L_heated/2)} m')



# #===Ze is when x and xe are equal===
# def equilibrium_point_eq(z):
#     #===Calculate inlet properties===
#     h_in = steamTable.h_pt(p_in, T_in) # [kJ/kg]

#     #===Calculate parameters at position z===
#     q_dz = lambda z: heat_flux_pin(z) * (math.pi * D_rod ) # [W]
#     h_m_z = h_in + (quad(q_dz, -L_heated/2, z)[0]/ mass_flow_rate_F)  / 1000 # [kJ/kg]

#     #===Calculate dynamic quality at position z===
#     h_l_sat = steamTable.hL_p(p_in) # [kJ/kg]
#     h_g_sat = steamTable.hV_p(p_in) # [kJ/kg]
#     x_e = (h_m_z - h_l_sat) / (h_g_sat - h_l_sat) # [-]

#     return x_e - 1.0
# zE_solution = fsolve(equilibrium_point_eq, 0.03) # Initial guess at the middle of the heated length
# zE = zE_solution[0]
# print(f'Axial position of equilibrium point from the channel inlet (zE): {zE + (L_heated/2)} m')


def void_fraction_HEM(x, p):
    rho_L = steamTable.rhoL_p(p)
    rho_G = steamTable.rhoV_p(p)
    return x/(x +  (rho_G/rho_L) * (1 - x)) # [-]

def rho_m(x, p):
    rho_L = steamTable.rhoL_p(p)
    rho_G = steamTable.rhoV_p(p)
    return void_fraction_HEM(x, p) * rho_G + (1 - void_fraction_HEM(x, p)) * rho_L # [kg/m³] 

#=== FINAL EXAM CHALLENGE : Plot the total pressure drop through the subchanneland its individual components in function of the mass flow rate ===
#=== for mass flow rates ranging from 0 to 2.5 kg/s.===
#=== Use McAdams correlation for the two phase pressure drop and Coolebrook for the monophase pressure drop===
#=== Assume constant material properties and use the inlet pressure to calculate those for the subcooled region===
#=== E.g. Four graphs (liquid, two phase, vapour, total) with each three components (acceleration, gravity and friction) ===

#===Before : Write enthalpy and equilibrium quality in terms of height, mass flux and heating===

#===Graphs of h and xe  along z for mass flow rate = 2.5 kg/s===
mass_flow_rate_graphs = 2.5 # [kg/s]
n_axial_G = 100
z_values_G = [i * L_heated / n_axial_G - (L_heated/2) for i in range(n_axial_G +1)] # [m]

h_values_G = [h_m_z(z, mass_flow_rate_graphs) for z in z_values_G]
plt.plot(h_values_G, z_values_G, label='Enthalpy along z')
plt.xlabel('Enthalpy [kJ/kg]')
plt.ylabel('Axial Position z [m]')
plt.title('Enthalpy  vs Axial Position z (Mass Flow Rate = 2.5 kg/s)')
plt.legend()
plt.grid()
plt.show()

x_e_values_G = [x_e_z(z, mass_flow_rate_graphs, p_in) for z in z_values_G]  
plt.plot(x_e_values_G, z_values_G, label='Equilibrium Quality along z')
plt.xlabel(' Equilibrium Quality [-]')
plt.ylabel('Axial Position z [m]')
plt.title('Equilibrium Quality vs Axial Position z (Mass Flow Rate = 2.5 kg/s)')
plt.legend()
plt.grid()
plt.show()

n_mass_flow_rate = 10
mass_flow_rate_list_final = [i * 2.5 / n_mass_flow_rate for i in range(1, n_mass_flow_rate+1)] # [kg/s]



#===First graph: Liquid region (from 0 to zB)===

dp_acc_list_liquid = []
dp_fric_list_liquid = []
dp_grav_list_liquid = []
dp_tot_list_liquid = []

for mass_flow_rate in mass_flow_rate_list_final:
    h_m_mid_liquid = h_m_z(zB/2, mass_flow_rate) # [kJ/kg]
    rho_liquid = steamTable.rho_ph(p_in, h_m_mid_liquid) # [kg/m³]
    mu_liquid = steamTable.my_ph(p_in, h_m_mid_liquid) # [Pa.s]
    u_liquid = mass_flow_rate / (rho_liquid * A_flow)
    Re_liquid = rho_liquid * u_liquid * Hydraulic_diameter / mu_liquid
    # print(f'Mass Flow Rate: {mass_flow_rate} kg/s, Liquid Reynolds number: {Re_liquid}')
    #===Acceleration pressure drop===
    dp_acc_liquid = 0 # No acceleration in subcooled region (in reality there is a small acceleration due to change in density with temperature)
    dp_acc_list_liquid.append(dp_acc_liquid)
    #===Friction pressure drop===
    # Look at the dependance of corelation with Reynolds number
    if Re_liquid < 2100:
        f_liquid = 64 / Re_liquid  # Laminar flow
    elif 2100 <= Re_liquid < 3000:
        # Transitional flow, use linear interpolation between laminar and turbulent
        f_laminar = 64 / Re_liquid
        f_turbulent = colebrook_root(Re_liquid, relative_wall_roughness)
        f_liquid = f_laminar + (f_turbulent - f_laminar) * ((Re_liquid - 2100) / (3000 - 2100))
    else:
        f_liquid = colebrook_root(Re_liquid, relative_wall_roughness)
    dp_friction_liquid = f_liquid * ((zB + (L_heated/2)) / Hydraulic_diameter) * (rho_liquid * u_liquid**2 / 2) / 1e5
    dp_fric_list_liquid.append(dp_friction_liquid)
    #===Gravity pressure drop===
    dp_gravity_liquid = rho_liquid * g * (zB + (L_heated/2)) / 1e5
    dp_grav_list_liquid.append(dp_gravity_liquid)
    #===Total pressure drop===
    dp_total_liquid = dp_acc_liquid + dp_friction_liquid + dp_gravity_liquid
    dp_tot_list_liquid.append(dp_total_liquid)

#comparer la température en à zb avec celle de l'enthalpie calculée et celle de la formule slide 16 cours 14

#===Plot the liquid region pressure drop components===
plt.plot(mass_flow_rate_list_final, dp_acc_list_liquid, label='Acceleration Pressure Drop')
plt.plot(mass_flow_rate_list_final, dp_fric_list_liquid, label='Friction Pressure Drop')
plt.plot(mass_flow_rate_list_final, dp_grav_list_liquid, label='Gravity Pressure Drop')
plt.plot(mass_flow_rate_list_final, dp_tot_list_liquid, label='Total Pressure Drop')
plt.xlabel('Mass Flow Rate [kg/s]')
plt.ylabel('Pressure Drop [bar]')
plt.title('Pressure Drop Components in Liquid Region vs Mass Flow Rate')
plt.legend()
plt.grid()
plt.show()


#===Second graph: Two phase region (from zB to zV)===
dp_acc_list_two_phase = []
dp_fric_list_two_phase = []
dp_grav_list_two_phase = []
dp_tot_list_two_phase = []

#===Levy correlation to link x and xe===
#=== rho_L use for friction should be at saturation===



for mass_flow_rate in mass_flow_rate_list_final:
    G_m_TP = mass_flow_rate / A_flow # [kg/m²/s]

    #===Calculate properties at bulk boiling point===
    h_m_zB = h_m_z(zB, mass_flow_rate) # [kJ/kg]
    x_e_in_TP = (x_e_z(zB, mass_flow_rate, p_in)) # [-]
    #ECRIRE UNE FONCTION AU DESSUS QUI DONNE X EN FONCTION DE X_E

    # rho_L_TP = steamTable.rhoL_p(p_in) # [kg/m³]
    # rho_G_TP = steamTable.rhoV_p(p_in) # [kg/m³]
    # void_fraction_TP = void_fraction(x_in_TP, p_in) # [-]
    # rho_m_TP = rho_m(x_in_TP, p_in) # [kg/m³]
    # print(f'Mass Flow Rate: {mass_flow_rate} kg/s, Mixture density at bulk boiling point: {rho_m_TP} kg/m³')

    # #===Acceleration pressure drop===
    # 
    # dp_acc_TP = G_m_TP**2 * (1/rho_m_plus(rho_L_TP, rho_G_TP, 0.15) - 1/rho_m_plus(rho_L_TP, rho_G_TP, x_in_TP)) / 1e5
    # dp_acc_list_two_phase.append(dp_acc_TP)
    # #===Friction pressure drop (using McAdams correlation)===
    # u_m_TP = mass_flow_rate / (rho_m_TP * A_flow)
    # Re_m_TP = rho_m_TP * u_m_TP * Hydraulic_diameter / steamTable.my_pt(p_in, T_in)
    # if Re_m_TP < 30000:
    #     f_TP = 64 / Re_m_TP  # Laminar flow
    # if Re_m_TP > 30000:
    #     f_TP = 0.184 * Re_m_TP**(-0.2)
    # dp_friction_TP = f_TP * ((zV - zB) / Hydraulic_diameter) * (rho_m_TP * u_m_TP**2 / 2) / 1e5
    # dp_fric_list_two_phase.append(dp_friction_TP)
    # #===Gravity pressure drop===
    # dp_gravity_TP = rho_m_TP * g * (zV - zB) / 1e5
    # dp_grav_list_two_phase.append(dp_gravity_TP)
    # #===Total pressure drop===
    # dp_total_TP = dp_acc_TP + dp_friction_TP + dp_gravity_TP
    # dp_tot_list_two_phase.append(dp_total_TP)
# #===Plot the two phase region pressure drop components===   
# plt.plot(mass_flow_rate_list_final, dp_acc_list_two_phase, label='Acceleration Pressure Drop')
# plt.plot(mass_flow_rate_list_final, dp_fric_list_two_phase, label='Friction Pressure Drop')
# plt.plot(mass_flow_rate_list_final, dp_grav_list_two_phase, label='Gravity Pressure Drop')
# plt.plot(mass_flow_rate_list_final, dp_tot_list_two_phase, label='Total Pressure Drop')
# plt.xlabel('Mass Flow Rate [kg/s]')
# plt.ylabel('Pressure Drop [bar]')
# plt.title('Pressure Drop Components in Two Phase Region vs Mass Flow Rate')
# plt.legend()
# plt.grid()
# plt.show()
