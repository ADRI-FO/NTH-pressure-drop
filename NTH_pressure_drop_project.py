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
# import numpy as np
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
        if Re_1 > 1e6:
            print('Reynolds number out of range for McAdams correlation: Re_C = ', Re_1)
    dp_friction = f_1 * (L_heated / D_pipe_1) * (rho_subcooled * u_1**2 / 2) / 1e5 # [bar]
    dp_friction_list_1.append(dp_friction)
#===Plot the friction pressure drop===
plt.plot(mass_flow_rate_list, dp_friction_list_1, label='D_pipe_1 = 2e-2 m, McAdams correlation')
plt.xlabel('Mass Flow Rate [kg/s]')
plt.ylabel('Friction Pressure Drop [bar]')
plt.title('Friction Pressure Drop vs Mass Flow Rate')
plt.legend()
plt.grid()
plt.show()

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
        if Re_2 > 1e6:
            print('Reynolds number out of range for McAdams correlation: Re_C = ', Re_2)
    dp_friction = f_2 * (L_heated / Hydraulic_diameter) * (rho_subcooled * u_2**2 / 2) / 1e5 # [bar]
    dp_friction_list_2.append(dp_friction)
#===Plot the friction pressure drop===
plt.plot(mass_flow_rate_list, dp_friction_list_2, label='BWR subchannel, McAdams correlation')
plt.xlabel('Mass Flow Rate [kg/s]')
plt.ylabel('Friction Pressure Drop [bar]')
plt.title('Friction Pressure Drop vs Mass Flow Rate')
plt.legend()
plt.grid()
plt.show()

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
plt.plot(mass_flow_rate_list, dp_friction_list_3, label='BWR subchannel, Colebrook-White equation')
plt.xlabel('Mass Flow Rate [kg/s]')
plt.ylabel('Friction Pressure Drop [bar]')
plt.title('Friction Pressure Drop vs Mass Flow Rate')
plt.legend()
plt.grid()
plt.show()


#===Second exercise: Calculate two phase pressure drop through the channel ===
#=== D : First assume the inlet equilibrium quality is 0.15 and the channel is adiabatic (no heating) ===
x_in_D = 0.15
rho_L_D = steamTable.rhoL_p(p_in) # [kg/m³]
rho_G_D = steamTable.rhoV_p(p_in) # [kg/m³]
void_fraction_D = 1/(1 + ((1 - x_in_D)/x_in_D) * (rho_G_D/rho_L_D)) # [-]
rho_m_D = void_fraction_D * rho_G_D + (1 - void_fraction_D) * rho_L_D # [kg/m³]
print(f'Void fraction at inlet quality of 0.15: {void_fraction_D}')
print(f'Mixture density at inlet quality of 0.15: {rho_m_D} kg/m³')
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
        if Re_L0_D > 1e6:
            print('Reynolds number out of range for McAdams correlation: Re_C = ', Re_L0_D)
    friction_factor_L0_list_D.append(friction_factor_L0_D)
    friction_factor_TP_D = friction_factor_L0_D * (mu_ratio_D)**0.2
    dp_fric = (friction_factor_TP_D * (L_heated / Hydraulic_diameter) * (G_m_D**2 / (2 * rho_m_D))) / 1e5 # [bar]
    dp_fric_list_D.append(dp_fric)
    dp_tot = dp_fric + dp_grav_D + dp_acc__D
    dp_tot_list_D.append(dp_tot)
#===Plot friction factor for liquid only at inlet conditions===
plt.plot(mass_flow_rate_list, friction_factor_L0_list_D, label='Friction factor for liquid only at inlet conditions')
plt.xlabel('Mass Flow Rate [kg/s]')
plt.ylabel('Friction Factor [-]')
plt.ylim(0, 0.03)
plt.title('Friction Factor vs Mass Flow Rate')
plt.legend()
plt.grid()
plt.show()
#===Plot the two phase pressure drop===
plt.plot(mass_flow_rate_list, dp_tot_list_D, label='Inlet quality = 0.15, Adiabatic channel')
plt.xlabel('Mass Flow Rate [kg/s]')
plt.ylabel('Two Phase Pressure Drop [bar]')
plt.title('Two Phase Pressure Drop vs Mass Flow Rate')
plt.legend()
plt.grid()
plt.show()


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
        if Re_L0_E > 1e6:
            print('Reynolds number out of range for McAdams correlation: Re_C = ', Re_L0_E)
    #===Integration over the heated length===
    dp_grav_E_dz = lambda z: 9.81 * (void_fraction_E(z) * rho_G_E + (1 - void_fraction_E(z)) *rho_L_E) # [Pa/m]
    dp_grav_E, _  = quad(dp_grav_E_dz, 0, L_heated)
    dp_grav_list_E.append(dp_grav_E / 1e5) # [bar]

    phi_2_L0_E = lambda z: mu_ratio_E(z)**0.2 * (rho_L_E / rho_m_E(z))
    dp_fric_E_dz = lambda z: (friction_factor_L0_E * (G_m_E**2) / (Hydraulic_diameter * 2 * rho_L_E)) * phi_2_L0_E(z) # [Pa/m]
    dp_fric_E, _ = quad(dp_fric_E_dz, 0, L_heated)
    dp_fric_list_E.append(dp_fric_E / 1e5) # [bar]

    dp_acc_E_dz = G_m_E**2 * (1/rho_L_E + (1/rho_G_E - 1/rho_L_E) * (x_out_E - x_in_E)/L_heated) # [Pa/m]
    dp_acc_E = dp_acc_E_dz * L_heated
    dp_acc_list_E.append(dp_acc_E / 1e5) # [bar] # too high value compared to others, need to check
    
    print(f'Mass Flow Rate: {mass_flow_rate} kg/s, dp_grav: {dp_grav_E/1e5} bar, dp_acc: {dp_acc_E/1e5} bar, dp_fric: {dp_fric_E/1e5} bar')
    
    dp_tot_E = (dp_grav_E + dp_acc_E + dp_fric_E) / 1e5 # [bar]
    dp_tot_list_E.append(dp_tot_E)
    
    
#===Plot the two phase pressure drop===
plt.plot(mass_flow_rate_list, dp_tot_list_E, label='Inlet quality = 0.0, Outlet quality = 0.15')
plt.xlabel('Mass Flow Rate [kg/s]')
plt.ylabel('Two Phase Pressure Drop [bar]')
plt.title('Two Phase Pressure Drop vs Mass Flow Rate')
plt.legend()
plt.grid()
plt.show()


#===Third exercise: Axial heat flux distribution===
#===F: Find axial position where the equilibrium quality is zero for a mass flow rate of 0,1kg/s===
#===This bulk boiling point (zB). Assume cosine axial heat distribution neglecting pressure drop for material properties===

#===Cosinusoïdal axial heat flux distribution===
#=== z = 0 at the middle of the heated length ===
#=== Integral of cos(z) from -pi/2 to pi/2 = 2 ===
# Therefore, maximum heat flux is the average heat flux  divided by 2/pi
mass_flow_rate_F = 0.1 # [kg/s]
heat_flux_pin_max = heat_flux_pin_avg / (2 / math.pi) # [W/m²]
heat_flux_pin = lambda z: heat_flux_pin_max * (math.cos(math.pi * z / L_heated)) # [W/m²], z = 0 at the middle of the heated length

#=== Graphical representation of the axial heat flux distribution ===
n_axial = 100
z_values = [i * 1/n_axial for i in range(-int((L_heated/2)*n_axial), int((L_heated/2)*n_axial))] # [m]
heat_flux_values = [heat_flux_pin(z) for z in z_values]
plt.plot(heat_flux_values, z_values)
plt.ylabel('Axial Position [m]')
plt.xlabel('Heat Flux [W/m²]')
plt.title('Axial Heat Flux Distribution')
plt.grid()
plt.show()

# #===Verification of the average heat flux calculation===
# integral_heat_flux = sum(heat_flux_values) * 1/n_axial # [W/m²]
# average_heat_flux = integral_heat_flux / L_heated # [W/m²]
# print(f'Calculated Average Heat Flux: {average_heat_flux} W/m²')

#===Find axial position of bulk boiling point (zB)===
def bulk_boiling_point_eq(z):
    #===Calculate inlet properties===
    h_in = steamTable.h_pt(p_in, T_in) # [kJ/kg]
    rho_in = steamTable.rho_ph(p_in, h_in) # [kg/m³]
    
    #===Calculate enthalpy at position z===
    q_dz = lambda z: heat_flux_pin(z) * (math.pi * D_rod ) # [W]
    h_z = h_in + (quad(q_dz, -L_heated/2, z)[0]/ mass_flow_rate_F)  / 1000 # [kJ/kg]

    T_z = steamTable.t_ph(p_in, h_z) # [°C]
    #===Check if the temperature at position z is equal to the saturation temperature===
    T_sat = steamTable.tsat_p(p_in) # [°C]
    return T_z - T_sat

zB_solution = fsolve(bulk_boiling_point_eq, 0.01) # Initial guess at the middle of the heated length
zB = zB_solution[0]
print(f'Axial position of bulk boiling point (zB): {zB} m')



   


#===B Find axial position where the flow (dynamic quality) is one▪This is werewe have saturated vapour (zV)===


# #===Plot the total pressure drop and its individual components ===
# #=== in function of the mass flow rate through the subchannel. ===
# mass_flow_rate_list = []
# dp_total_list = []
# dp_friction_list = []
# dp_elevation_list = []
# dp_acceleration_list = []

# for mdot in range(0, 2.5, 10): # [kg/s]
#     mass_flow_rate_list.append(mdot)
    
#     #===Calculate inlet properties===
#     h_in = steamTable.h_pt(p_in, T_in) # [kJ/kg]
#     rho_in = steamTable.rho_ph(p_in, h_in) # [kg/m³]
#     v1 = 1 / rho_in # [m³/kg]
    
#     #===Calculate outlet properties===
#     G = mdot / A_flow # [kg/m²/s]
#     q_total = heat_flux_pin_avg * (3.1416 * D_rod * L_heated) # [W]
#     h2 = h_in + q_total / mdot / 1000 # [kJ/kg]
#     p2 = p_in # [bar], assuming negligible pressure drop for outlet enthalpy calculation
#     T2 = steamTable.t_ph(p2, h2) # [°C]
#     rho2 = steamTable.rho_ph(p2, h2) # [kg/m³]
#     v2 = 1 / rho2 # [m³/kg]
    
#     #===Calculate pressure drop components===
#     # Frictional pressure drop
#     Re_C = G * D_rod / steamTable.my_ph(p_in, h_in) # Reynolds number
#     f = steamTable.friction_factor(Re_C, relative_wall_roughness / D_rod) # Darcy friction factor
#     dp_friction = f * (L_heated / D_rod) * (G**2 / (2 * rho_in)) / 100000 # [bar]
    
#     # Elevation pressure drop
#     g = 9.81 # [m/s²]
#     dp_elevation = rho_in * g * L_heated / 100000 # [bar]
    
#     # Acceleration pressure drop
#     dp_acceleration = G**2 * (v2 - v1) / 100000 # [bar]
    
#     # Total pressure drop
#     dp_total = dp_friction + dp_elevation + dp_acceleration
    
#     dp_total_list.append(dp_total)
#     dp_friction_list.append(dp_friction)
#     dp_elevation_list.append(dp_elevation)
#     dp_acceleration_list.append(dp_acceleration)