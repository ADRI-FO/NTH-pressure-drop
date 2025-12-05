#===IMPORT PACKAGES===
from pyXSteam.XSteam import XSteam
steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS) # m/kg/sec/°C/bar/W
import matplotlib.pyplot as plt
import math as math
import tkinter as tk
from tkinter import messagebox
from iapws import IAPWS97
from CoolProp.CoolProp import PropsSI
from scipy.optimize import root
from scipy.optimize import fsolve
from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy.integrate import cumulative_trapezoid
import numpy as np
# import pandas as pd
# import six
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

# #===First exercise: Calculate friction monophase pressure drop through channel without heating===
# #=== For a mass flow rate ranging from 0 to 2.5 kg/s ===
# rho_subcooled = steamTable.rho_pt(p_in, T_in) # [kg/m³]
# n_mass_flow_rate = 50
# mass_flow_rate_list = [i * 2.5 / n_mass_flow_rate for i in range(1, n_mass_flow_rate+1)] # [kg/s]
# dp_friction_list_1 = []
# #===A : First for a pipe diameter equal to 2e-2 [m] and with McAdams correlation
# #=== f = 0,184 * Re_C^(-0,2) if 3e4 < Re_C < 1e6 (turbulent) ===
# D_pipe_1 = 2e-2 # [m]
# A_pipe = math.pi * (D_pipe_1/2)**2 # [m²]

# for mass_flow_rate in mass_flow_rate_list:
#     u_1 = mass_flow_rate/(rho_subcooled * A_pipe)
#     Re_1 = rho_subcooled * u_1 * D_pipe_1 / steamTable.my_pt(p_in, T_in)
#     if Re_1 < 30000:
#         f_1 = 64 / Re_1  # Laminar flow
#     if Re_1 > 30000:
#         f_1 = 0.184 * Re_1**(-0.2)
#         # if Re_1 > 1e6:
#         #     print('Reynolds number out of range for McAdams correlation: Re_C = ', Re_1)
#     dp_friction = f_1 * (L_heated / D_pipe_1) * (rho_subcooled * u_1**2 / 2) / 1e5 # [bar]
#     dp_friction_list_1.append(dp_friction)
# # #===Plot the friction pressure drop===
# # plt.plot(mass_flow_rate_list, dp_friction_list_1, label='D_pipe_1 = 2e-2 m, McAdams correlation')
# # plt.xlabel('Mass Flow Rate [kg/s]')
# # plt.ylabel('Friction Pressure Drop [bar]')
# # plt.title('Friction Pressure Drop vs Mass Flow Rate')
# # plt.legend()
# # plt.grid()
# # plt.show()

#===B : Then for the actual BWR subchannel still with McAdams correlation===
A_flow = (Pitch**2 - math.pi * (D_rod/2)**2) # [m²]
Perimeter_heated = math.pi * D_rod # [m]
Hydraulic_diameter = 4 * A_flow / Perimeter_heated # [m] 
# dp_friction_list_2 = []
# for mass_flow_rate in mass_flow_rate_list:
#     u_2 = mass_flow_rate/(rho_subcooled * A_flow)
#     Re_2 = rho_subcooled * u_2 * Hydraulic_diameter / steamTable.my_pt(p_in, T_in)
#     if Re_2 < 30000:
#         f_2 = 64 / Re_2  # Laminar flow
#     if Re_2 > 30000:
#         f_2 = 0.184 * Re_2**(-0.2)
#         # if Re_2 > 1e6:
#         #     print('Reynolds number out of range for McAdams correlation: Re_C = ', Re_2)
#     dp_friction = f_2 * (L_heated / Hydraulic_diameter) * (rho_subcooled * u_2**2 / 2) / 1e5 # [bar]
#     dp_friction_list_2.append(dp_friction)
# # #===Plot the friction pressure drop===
# # plt.plot(mass_flow_rate_list, dp_friction_list_2, label='BWR subchannel, McAdams correlation')
# # plt.xlabel('Mass Flow Rate [kg/s]')
# # plt.ylabel('Friction Pressure Drop [bar]')
# # plt.title('Friction Pressure Drop vs Mass Flow Rate')
# # plt.legend()
# # plt.grid()
# # plt.show()

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

# dp_friction_list_3 = []
# for mass_flow_rate in mass_flow_rate_list:
#     u_C = mass_flow_rate / (rho_subcooled * A_flow)
#     Re_C = rho_subcooled * u_C * Hydraulic_diameter / steamTable.my_pt(p_in, T_in)
#     f = 0.184 * Re_C**(-0.2)   # Initial guess
#     rel_rough = relative_wall_roughness   # already ε/D
#     for _ in range(50):
#         f_old = f
#         f = 1.0 / ( -2.0 * math.log10(rel_rough/3.7 + 2.51/(Re_C*math.sqrt(f_old))) )**2
#         if abs(f - f_old) < 1e-8:
#             break
#     friction_factor_C = f
#     # Alternatively, use the root-finding function:
#     #friction_factor_C = colebrook_root(Re_C, relative_wall_roughness)
#     dp_C = friction_factor_C * (L_heated / Hydraulic_diameter) * (rho_subcooled * u_C**2 / 2) / 1e5
#     dp_friction_list_3.append(dp_C)
# #===Plot the friction pressure drop===
# # plt.plot(mass_flow_rate_list, dp_friction_list_3, label='BWR subchannel, Colebrook-White equation')
# # plt.xlabel('Mass Flow Rate [kg/s]')
# # plt.ylabel('Friction Pressure Drop [bar]')
# # plt.title('Friction Pressure Drop vs Mass Flow Rate')
# # plt.legend()
# # plt.grid()
# # plt.show()


# #===Second exercise: Calculate two phase pressure drop through the channel ===
# #=== D : First assume the inlet equilibrium quality is 0.15 and the channel is adiabatic (no heating) ===
# x_in_D = 0.15
# rho_L_D = steamTable.rhoL_p(p_in) # [kg/m³]
# rho_G_D = steamTable.rhoV_p(p_in) # [kg/m³]
# void_fraction_D = 1/(1 + ((1 - x_in_D)/x_in_D) * (rho_G_D/rho_L_D)) # [-]
# rho_m_D = void_fraction_D * rho_G_D + (1 - void_fraction_D) * rho_L_D # [kg/m³]
# # print(f'Void fraction at inlet quality of 0.15: {void_fraction_D}')
# # print(f'Mixture density at inlet quality of 0.15: {rho_m_D} kg/m³')
# dp_acc__D = 0 # No acceleration in adiabatic case without change in quality
# friction_factor_L0_list_D = []
# dp_grav_D = rho_m_D * 9.81 * L_heated / 1e5 # [bar]
# dp_tot_list_D = []
# dp_fric_list_D = []
# for mass_flow_rate in mass_flow_rate_list:
#     G_m_D = mass_flow_rate / A_flow # [kg/m²/s]
#     mu_L_D = PropsSI('VISCOSITY', 'P', p_in * 1e5, 'Q', 0, 'Water') # [Pa.s]
#     mu_G_D = PropsSI('VISCOSITY', 'P', p_in * 1e5, 'Q', 1, 'Water') # [Pa.s]
#     mu_ratio_D = 1/(1 + x_in_D * (mu_L_D/mu_G_D - 1)) # [-]
#     u_L0_D = mass_flow_rate/(rho_L_D * A_flow) # [m/s]   # Liquid only velocity at inlet conditions (all the mass flow rate is considered to be liquid)
#     Re_L0_D = rho_L_D * u_L0_D * Hydraulic_diameter / mu_L_D
#     if Re_L0_D < 30000:
#         friction_factor_L0_D = 64 / Re_L0_D  # Laminar flow
#     if Re_L0_D > 30000:
#         friction_factor_L0_D = 0.184 * Re_L0_D**(-0.2)
#         # if Re_L0_D > 1e6:
#         #     print('Reynolds number out of range for McAdams correlation: Re_C = ', Re_L0_D)
#     friction_factor_L0_list_D.append(friction_factor_L0_D)
#     friction_factor_TP_D = friction_factor_L0_D * (mu_ratio_D)**0.2
#     dp_fric = (friction_factor_TP_D * (L_heated / Hydraulic_diameter) * (G_m_D**2 / (2 * rho_m_D))) / 1e5 # [bar]
#     dp_fric_list_D.append(dp_fric)
#     dp_tot = dp_fric + dp_grav_D + dp_acc__D
#     dp_tot_list_D.append(dp_tot)
# # #===Plot friction factor for liquid only at inlet conditions===
# # plt.plot(mass_flow_rate_list, friction_factor_L0_list_D, label='Friction factor for liquid only at inlet conditions')
# # plt.xlabel('Mass Flow Rate [kg/s]')
# # plt.ylabel('Friction Factor [-]')
# # plt.ylim(0, 0.03)
# # plt.title('Friction Factor vs Mass Flow Rate')
# # plt.legend()
# # plt.grid()
# # plt.show()
# # #===Plot the two phase pressure drop===
# # plt.plot(mass_flow_rate_list, dp_tot_list_D, label='Inlet quality = 0.15, Adiabatic channel')
# # plt.xlabel('Mass Flow Rate [kg/s]')
# # plt.ylabel('Two Phase Pressure Drop [bar]')
# # plt.title('Two Phase Pressure Drop vs Mass Flow Rate')
# # plt.legend()
# # plt.grid()
# # plt.show()

def rho_m(x, p):
    rho_L = steamTable.rhoL_p(p)
    rho_G = steamTable.rhoV_p(p)
    return 1 / (x / rho_G + (1 - x) / rho_L)


# #===E : Then assume the inlet and outlet equilibrium quality are 0 and 0.15 respectively.===
# #=== Thus instead of using a heat flux as input it is assumed that the quality evolves linearly from inlet to outlet===
# x_in_E = 0.0
# x_out_E = 0.15
# x_E = lambda z: x_in_E + (x_out_E - x_in_E) * (z / L_heated) # [-], linear evolution of quality along the heated length
# rho_L_E = steamTable.rhoL_p(p_in) # [kg/m³]
# rho_G_E = steamTable.rhoV_p(p_in) # [kg/m³]

# dp_tot_list_E = []
# dp_acc_list_E = []
# dp_fric_list_E = []
# dp_grav_list_E = []
# # Not a good idea to descretize the lenght. Better to integrate directly over the lenght.
# for mass_flow_rate in mass_flow_rate_list:
#     void_fraction_E = lambda z: 1/(1 + ((1 - x_E(z))/x_E(z)) * (rho_G_E/rho_L_E)) if x_E(z) != 0 else 0.0  # [-]
#     rho_m_E = lambda z: void_fraction_E(z) * rho_G_E + (1 - void_fraction_E(z)) * rho_L_E  # [kg/m³]
#     G_m_E = mass_flow_rate / A_flow # [kg/m²/s]
#     mu_L_E = PropsSI('VISCOSITY', 'P', p_in * 1e5, 'Q', 0, 'Water') # [Pa.s]
#     mu_G_E = PropsSI('VISCOSITY', 'P', p_in * 1e5, 'Q', 1, 'Water') # [Pa.s]
#     mu_ratio_E = lambda z: 1/(1 + x_E(z) * (mu_L_E/mu_G_E - 1)) # [-]
#     u_L0_E = mass_flow_rate /(rho_L_E * A_flow ) # [m/s]
#     Re_L0_E = rho_L_E * u_L0_E * Hydraulic_diameter / mu_L_E
#     if Re_L0_E < 30000:
#         friction_factor_L0_E = 64 / Re_L0_E  # Laminar flow
#     if Re_L0_E > 30000:
#         friction_factor_L0_E = 0.184 * Re_L0_E**(-0.2)
#         # if Re_L0_E > 1e6:
#         #     print('Reynolds number out of range for McAdams correlation: Re_C = ', Re_L0_E)
#     #===Integration over the heated length===
#     dp_grav_E_dz = lambda z: 9.81 * (void_fraction_E(z) * rho_G_E + (1 - void_fraction_E(z)) *rho_L_E) # [Pa/m]
#     dp_grav_E, _  = quad(dp_grav_E_dz, 0, L_heated)
#     dp_grav_list_E.append(dp_grav_E / 1e5) # [bar]

#     phi_2_L0_E = lambda z: mu_ratio_E(z)**0.2 * (rho_L_E / rho_m_E(z))
#     dp_fric_E_dz = lambda z: (friction_factor_L0_E * (G_m_E**2) / (Hydraulic_diameter * 2 * rho_L_E)) * phi_2_L0_E(z) # [Pa/m]
#     dp_fric_E, _ = quad(dp_fric_E_dz, 0, L_heated)
#     dp_fric_list_E.append(dp_fric_E / 1e5) # [bar] 

#     dp_acc_E = G_m_E**2 * (1/rho_m(x_out_E, p_in) - 1/rho_m(x_in_E, p_in)) # [Pa]
#     dp_acc_list_E.append(dp_acc_E / 1e5) # [bar]

#     # print(f'Mass Flow Rate: {mass_flow_rate} kg/s, dp_grav: {dp_grav_E/1e5} bar, dp_acc: {dp_acc_E/1e5} bar, dp_fric: {dp_fric_E/1e5} bar}')

#     dp_tot_E = (dp_grav_E + dp_acc_E + dp_fric_E) / 1e5 # [bar]
#     dp_tot_list_E.append(dp_tot_E)
        
# # #===Plot the two phase pressure drop===
# # plt.plot(mass_flow_rate_list, dp_tot_list_E, label='Inlet quality = 0.0, Outlet quality = 0.15')
# # plt.xlabel('Mass Flow Rate [kg/s]')
# # plt.ylabel('Two Phase Pressure Drop [bar]')
# # plt.title('Two Phase Pressure Drop vs Mass Flow Rate')
# # plt.legend()
# # plt.grid()
# # plt.show()



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
    z = np.asarray(z)  # accepte scalaires et tableaux
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


# ---- Build integration grid ONCE ----
N = 800  # resolution
z_grid = np.linspace(-L_heated/2, L_heated/2, N)

# Heat flux per unit length (W/m) on the grid
q_grid = heat_flux_pin(z_grid) * (math.pi * D_rod)

# Precompute integral Q(z)
Q_grid = cumulative_trapezoid(q_grid, z_grid, initial=0)

# Build interpolation function
Q_interp = interp1d(z_grid, Q_grid, kind='cubic')

def h_m_z(z, mass_flow_rate):
    h_in = steamTable.h_pt(p_in, T_in)  # kJ/kg
    Q = float(Q_interp(z))              # W (integrated)
    return h_in + Q / (mass_flow_rate * 1000)  # kJ/kg


def x_e_z(z, mass_flow_rate, p):
    #===Calculate dynamic quality at position z===
    h_l_sat = steamTable.hL_p(p) # [kJ/kg]
    h_g_sat = steamTable.hV_p(p) # [kJ/kg]
    x_e = (h_m_z(z, mass_flow_rate) - h_l_sat) / (h_g_sat - h_l_sat) # [-]
    return x_e

#===Find axial position of bulk boiling point (zB)===
def bulk_boiling_point_eq(z, mass_flow_rate):
    #===Check if the temperature at position z is equal to the saturation temperature===
    h_l_sat = steamTable.hL_p(p_in) # [kJ/kg] #This is not perfect, normally should be at the pressure at z but need an iterative process
    return h_m_z(z, mass_flow_rate) - h_l_sat

def zB(mass_flow_rate):
    h_l_sat_top = steamTable.hL_p(p_in) #This is not perfect, normally should be at the pressure at L_heated/2 but need an iterative process

    # enthalpy at outlet
    h_at_outlet = h_m_z(L_heated/2, mass_flow_rate)

    # Check if boiling is reached
    if h_at_outlet < h_l_sat_top:
        # No boiling in whole channel
        return L_heated/2
    
    # Otherwise solve for zB
    sol = fsolve(bulk_boiling_point_eq, 0.0, args=(mass_flow_rate,))
    z_val = float(sol[0])

    # Safety clamp
    return min(max(z_val, -L_heated/2), L_heated/2)


zB_F = zB(mass_flow_rate_F)
#print(f'Axial position of bulk boiling point from the channel inlet (zB): {zB_F + (L_heated/2)} m')


#===Find axial position where the flow (dynamic quality) is one. This is were we have saturated vapour (zV)===
def saturated_vapour_point_eq(z, mass_flow_rate, p):
    #===Calculate dynamic quality at position z===
    return x_e_z(z, mass_flow_rate, p) - 1.0

def zV(mass_flow_rate, p):

    # Quality at channel outlet
    x_out = x_e_z(L_heated/2, mass_flow_rate, p)

    # If we never reach x=1, return outlet
    if x_out < 1.0:
        return L_heated/2

    # Otherwise solve for saturated vapor point
    sol = fsolve(saturated_vapour_point_eq, 0.0, args=(mass_flow_rate, p))
    z_val = float(sol[0])

    # Safety clamp
    return min(max(z_val, -L_heated/2), L_heated/2)


zV_F = zV(mass_flow_rate_F, p_in)
#print(f'Axial position of saturated vapour point from the channel inlet (zV): {zV_F + (L_heated/2)} m')


#===Find axial position where the flow (dynamic quality) is zero: bubble detachment (zD) (start of subcool boiling)===
#Using Saha and Zuber correlation or Stanton for bubble detachment (low flow rates), and also do you have to conscider the local or the average heat flux? (check in the book)
def bubble_detachment_point_eq(z, mass_flow_rate, p):
    # Clamp z pour rester dans le canal
    z = max(min(z, L_heated/2), -L_heated/2)

    # Paramètres de flux
    G_m = mass_flow_rate / A_flow # [kg/m²/s]
    cp_L = steamTable.CpL_p(p) * 1e3
    k_L  = steamTable.tcL_p(p)
    Dh = Hydraulic_diameter
    q = heat_flux_pin(z)
    Tsat = steamTable.tsat_p(p)

    # Peclet
    Pe = G_m * Dh * cp_L / k_L

    # T_bulk
    if Pe < 7e4:
        T_bulk = Tsat - q*Dh/(455*k_L)
    else:
        T_bulk = Tsat - q/(6.5e-3*G_m*cp_L)

    # Clamp T_bulk pour SteamTable
    T_bulk = max(Tsat-50, min(T_bulk, Tsat))

    # Si encore hors domaine → retourner grand résidu pour fsolve
    if T_bulk <= 0 or T_bulk > 2*Tsat:
        return 1e6

    h_bulk = steamTable.h_pt(p, T_bulk)
    return h_m_z(z, mass_flow_rate) - h_bulk


def zD(mass_flow_rate, p):
    """Compute bubble detachment point robustly."""
    # --- Enthalpy at inlet/outlet ---
    h_in  = h_m_z(-L_heated/2, mass_flow_rate)
    h_out = h_m_z( L_heated/2, mass_flow_rate)

    # --- Reference h_bulk at center ---
    Tsat = steamTable.tsat_p(p)
    h_bulk_center = steamTable.h_pt(p, Tsat - heat_flux_pin(0) * Hydraulic_diameter / (455 * steamTable.tcL_p(p)))

    # --- Check if detachment exists inside the channel ---
    if h_out < h_bulk_center:
        return L_heated/2     # no detachment
    if h_in > h_bulk_center:
        return -L_heated/2    # detachment occurs immediately

    # --- Solve for root with safe initial guess -1 ---
    sol = fsolve(bubble_detachment_point_eq, x0=0.0, args=(mass_flow_rate, p))
    z_val = float(sol[0])

    # --- Clamp final solution to domain ---
    return max(min(z_val, L_heated/2), -L_heated/2)


# zD_sol = fsolve(bubble_detachment_point_eq, -1, args=(mass_flow_rate_F, p_in))[0]
zD_new = zD(mass_flow_rate_F, p_in)
# print(f'value of zD (from middle of heated length): {zD_new} m')
# print(f'Axial position of bubble detachment point from the channel inlet (zD_new): {zD_new + (L_heated/2)} m')

def x_flow_z(z, mass_flow_rate, p):
    zD_x = zD(mass_flow_rate, p)
    xe_zD = x_e_z(zD_x, mass_flow_rate, p)
    epsilon = xe_zD * math.exp((x_e_z(z, mass_flow_rate, p)/xe_zD) - 1)
    if abs(epsilon) < 1e-4:
        return x_e_z(z, mass_flow_rate, p)
    else:
        return x_e_z(z, mass_flow_rate, p) - epsilon



def void_fraction_HEM(x, p):
    rho_L = steamTable.rhoL_p(p)
    rho_G = steamTable.rhoV_p(p)
    return x/(x +  (rho_G/rho_L) * (1 - x)) # [-]


def rho_m_plus(x,p): #velocotiy of each phase assumed to be radially uniform
    alpha = void_fraction_HEM(x, p)
    rho_L = steamTable.rhoL_p(p)
    rho_G = steamTable.rhoV_p(p)
    return 1 / (x**2/ (alpha * rho_G) + (1 - x)**2 / ((1 - alpha) * (rho_L)))


def friction_factor_1phase(Re, rel_rough):
    if Re < 2100:
        return 64 / Re  # Laminar flow
    elif 2100 <= Re < 3000:
        # Transitional flow, use linear interpolation between laminar and turbulent
        f_laminar = 64 / Re
        f_turbulent = colebrook_root(Re, rel_rough)
        return f_laminar + (f_turbulent - f_laminar) * ((Re - 2100) / (3000 - 2100))
    else:
        return colebrook_root(Re, rel_rough)
    
def McAdams_factor_1phase(Re):
    if Re < 2100:
        return 64 / Re
    elif 2100 <= Re < 3000:
        f_laminar = 64/Re
        f_turbulent_blach = 0.316*Re**(-0.25)
        return f_laminar + (f_turbulent_blach-f_laminar)*((Re - 2100) / (3000 - 2100))
    elif 3000<=Re<30000:
        f_turbulent_blach = 0.316*Re**(-0.25)
        return f_turbulent_blach
    else:
        f_turbulent_McAdams = 0.184*Re**(-0.2)
        return f_turbulent_McAdams
    
def phi_2_L0_McAdams(Re, x, p):
    mu_L = PropsSI('VISCOSITY', 'P', p * 1e5, 'Q', 0, 'Water')  # [Pa.s]
    mu_G = PropsSI('VISCOSITY', 'P', p * 1e5, 'Q', 1, 'Water')  # [Pa.s]
    mu_ratio = 1 / (1 + x * (mu_L / mu_G - 1))  # [-]
    rho_L = steamTable.rhoL_p(p)  # [kg/m³]
    if Re<=30000:
        return mu_ratio ** 0.25 * (rho_L / rho_m(x, p))  # [-]
    else:
        return mu_ratio ** 0.2 * (rho_L / rho_m(x, p))  # [-]

def omega_Jones(G_m, p):
    Gm_omega = G_m * 737.5621 #[lb/hr/ft²]
    p_omega = p * 14.5038 #[psi]
    if Gm_omega/1e6 <= 0.7:
        return 1.36 + 5e-4 * p_omega + 0.1 * Gm_omega/1e6 - 7.14e-4 * p_omega * (Gm_omega/1e6)
    else:
        return 1.26 + 4e-4 * p_omega + 0.119 * 1e6/Gm_omega + 2.8e-4 * p_omega * (1e6/Gm_omega)
    
def phi_2_L0_Jones(G_m, x, p):
    rho_L = steamTable.rhoL_p(p)  # [kg/m³]
    rho_V = steamTable.rhoV_p(p)  # [kg/m³]
    omega = omega_Jones(G_m, p)
    return omega * (1.2 * ((rho_L / rho_V) - 1) * x**0.824) + 1 # [-]

def Re(mass_flow_rate, mu, A_flow, Hydraulic_diameter):
    return mass_flow_rate * Hydraulic_diameter / (A_flow * mu)

# --- Functions called Tkinter ---
def choose_zB():
    global separation_point
    separation_point = "zB"
    messagebox.showinfo("Selection", "You have chosen zB as the separation point.")
    root_tk.destroy()

def choose_zD():
    global separation_point
    separation_point = "zD"
    messagebox.showinfo("Selection", "You have chosen zD as the separation point.")
    root_tk.destroy()

def ask_constant_properties():
    def set_choice(value):
        nonlocal choice
        choice = value
        window.destroy()

    choice = None
    window = tk.Tk()
    window.title("Property Selection")
    window.geometry("340x150")

    label = tk.Label(window, text="Constant properties along each region ?", font=("Arial", 12))
    label.pack(pady=15)

    # Buttons
    btn_yes = tk.Button(window, text="YES", width=12, command=lambda: set_choice("YES"))
    btn_yes.pack(side="left", padx=25, pady=10)

    btn_no = tk.Button(window, text="NO", width=12, command=lambda: set_choice("NO"))
    btn_no.pack(side="right", padx=25, pady=10)

    window.mainloop()
    return choice

#=== FINAL EXAM CHALLENGE : Plot the total pressure drop through the subchanneland its individual components in function of the mass flow rate ===
#=== for mass flow rates ranging from 0 to 2.5 kg/s.===
#=== Use McAdams correlation for the two phase pressure drop and Coolebrook for the monophase pressure drop===
#=== Assume constant material properties and use the inlet pressure to calculate those for the subcooled region===
#=== E.g. Four graphs (liquid, two phase, vapour, total) with each three components (acceleration, gravity and friction) ===

#=== ask the user to choose between zB and zD ===

# --- Tkinter windows---
root_tk = tk.Tk()
root_tk.title("Choose Separation Point")

label = tk.Label(root_tk, text="Choose the point to separate the liquid and two-phase region:")
label.pack(pady=10)

button_zB = tk.Button(root_tk, text="zB (Bulk Boiling Point)", command=choose_zB, width=30)
button_zB.pack(pady=5)

button_zD = tk.Button(root_tk, text="zD (Bubble Detachment)", command=choose_zD, width=30)
button_zD.pack(pady=5)

root_tk.mainloop()

# After closure of the window
if separation_point is None:
    print("No choice made. Defaulting to zD.")
    separation_point = "zD"

print("Separation point chosen:", separation_point)

#=== ask the user to choose if constant properties are used along each region ===

constant_properties_choice = ask_constant_properties()
print("Constant properties selected:", constant_properties_choice)


def choose_mc_adams():
    global heat_transfer_correlation_choice
    heat_transfer_correlation_choice = "McAdams"
    root_correlation.destroy()

def choose_jones():
    global heat_transfer_correlation_choice
    heat_transfer_correlation_choice = "Jones"
    root_correlation.destroy()

#--- Tkinter window for heat transfer correlation selection ---
# Create main window
root_correlation = tk.Tk()
root_correlation.title("Select Heat Transfer Correlation")

# Instructions
label = tk.Label(root_correlation, text="Choose the heat transfer correlation:")
label.pack(pady=10)

# Buttons
btn_mc_adams = tk.Button(root_correlation, text="Mc Adams correlation", width=25, command=choose_mc_adams)
btn_mc_adams.pack(pady=5)

btn_jones = tk.Button(root_correlation, text="Jones' correlation", width=25, command=choose_jones)
btn_jones.pack(pady=5)

# Start GUI
root_correlation.mainloop()

# Print the result (optional)
print("Correlation chosen:", heat_transfer_correlation_choice)

n_mass_flow_rate = 90 # max 90
mass_flow_rate_list_final = [i * 2.5 / n_mass_flow_rate for i in range(1, n_mass_flow_rate+1)] # [kg/s]

#===First graph: Liquid region (from 0 to zD (Zsc))===

dp_acc_list_liquid = []
dp_fric_list_liquid = []
dp_grav_list_liquid = []
dp_tot_list_liquid = []

for mass_flow_rate in mass_flow_rate_list_final:

    # Separation point
    if separation_point == "zB":
        z_sep = zB(mass_flow_rate)
    else:
        z_sep = zD(mass_flow_rate, p_in)

    L_liquid = z_sep + (L_heated/2)

    # -----------------------
    # INITIAL VALUES (inlet)
    # -----------------------
    rho_liquid_in = steamTable.rho_pt(p_in, T_in)      # [kg/m³]
    mu_liquid_in  = steamTable.my_pt(p_in, T_in)       # [Pa.s]

    # initial guess for pressure at outlet of liquid region
    p_guess = p_in

    if constant_properties_choice == "NO":
        # Non-constant properties: we will update them in the loop
        
        # -----------------------
        # ITERATION LOOP
        # -----------------------
        for k in range(50):

            # enthalpy at separation point
            h_liquid_out = h_m_z(z_sep, mass_flow_rate)

            # OUTLET PROPERTIES
            T_liquid_out  = steamTable.t_ph(p_guess, h_liquid_out)
            rho_liquid_out = steamTable.rho_ph(p_guess, h_liquid_out)
            mu_liquid_out  = PropsSI('VISCOSITY', 'P', p_guess*1e5, 'H', h_liquid_out*1e3, 'Water')

            # MEAN PROPERTIES
            rho_liquid_mean = 0.5 * (rho_liquid_in + rho_liquid_out)
            mu_liquid_mean  = 0.5 * (mu_liquid_in  + mu_liquid_out)

            # VELOCITY
            u_liquid = mass_flow_rate / (rho_liquid_mean * A_flow)

            # REYNOLDS WITH MEAN PROPERTIES
            Re_liquid = rho_liquid_mean * u_liquid * Hydraulic_diameter / mu_liquid_mean

            # FRICTION FACTOR (Colebrook)
            f_liquid = friction_factor_1phase(Re_liquid, relative_wall_roughness)

            # PRESSURE DROPS
            dp_acc_liquid = (mass_flow_rate/A_flow)**2 * (1/rho_liquid_out - 1/rho_liquid_in) / 1e5
            dp_friction_liquid = f_liquid * (L_liquid / Hydraulic_diameter) * (rho_liquid_mean * u_liquid**2 / 2) / 1e5
            dp_gravity_liquid = rho_liquid_mean * g * L_liquid / 1e5

            dp_total_iter = dp_acc_liquid + dp_friction_liquid + dp_gravity_liquid

            # NEW PRESSURE AT OUTLET
            p_new = p_in - dp_total_iter

            # convergence
            if abs(p_new - p_guess) < 1e-5:
                break
            
            p_guess = p_new #0.5*p_guess + 0.5*p_new # relaxation (stable)

    else:
        #Assuming constant material properties and use the inlet pressure to calculate those for the subcooled region:
        u_liquid_in = mass_flow_rate / (rho_liquid_in * A_flow) # [m/s]
        Re_liquid_in = Re(mass_flow_rate, mu_liquid_in, A_flow, Hydraulic_diameter)
        # print(f'Mass Flow Rate: {mass_flow_rate} kg/s, Liquid Reynolds number: {Re_liquid}')
        
        #===Pressure drops===
        dp_acc_liquid = 0
        f_liquid = friction_factor_1phase(Re_liquid_in, relative_wall_roughness)
        dp_friction_liquid = f_liquid * L_liquid / Hydraulic_diameter * (rho_liquid_in * u_liquid_in**2 / 2) / 1e5
        dp_gravity_liquid = rho_liquid_in * g * L_liquid/ 1e5
        dp_total_iter = dp_acc_liquid + dp_friction_liquid + dp_gravity_liquid

    # -----------------------
    # STORE FINAL VALUES
    # -----------------------
    dp_acc_list_liquid.append(dp_acc_liquid)
    dp_fric_list_liquid.append(dp_friction_liquid)
    dp_grav_list_liquid.append(dp_gravity_liquid)
    dp_tot_list_liquid.append(dp_total_iter)

    # Optionally print final state (debug)
    # print(f"mass={mass_flow_rate}: rho_mean={rho_liquid_mean}, mu_mean={mu_liquid_mean}, p_end={p_new}")


    # #===Total pressure drop===
    # dp_total_liquid = dp_acc_liquid + dp_friction_liquid + dp_gravity_liquid
    # dp_tot_list_liquid.append(dp_total_liquid)

#comparer la température en zb avec celle de l'enthalpie calculée et celle de la formule slide 16 cours 14

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

for idx, mass_flow_rate in enumerate(mass_flow_rate_list_final):
    # inlet pressure to two-phase region (bar)
    p_in_TP = p_in - dp_tot_list_liquid[idx]
    # print(f'For mass flow rate {mass_flow_rate} kg/s, pressure at two phase region inlet: {p_in_TP} bar')
    G_m_TP = mass_flow_rate / A_flow # [kg/m²/s]  

    if separation_point == "zB":
        zB_sol = zB(mass_flow_rate)
        if zB_sol == L_heated/2:
            # No two phase region
            dp_acc_list_two_phase.append(0.0)
            dp_fric_list_two_phase.append(0.0)
            dp_grav_list_two_phase.append(0.0)
            dp_tot_list_two_phase.append(0.0)
            continue
    else:
        zD_sol = zD(mass_flow_rate, p_in_TP)
        if zD_sol == L_heated/2:
            # No two phase region
            dp_acc_list_two_phase.append(0.0)
            dp_fric_list_two_phase.append(0.0)
            dp_grav_list_two_phase.append(0.0)
            dp_tot_list_two_phase.append(0.0)
            continue
    if constant_properties_choice == "YES":    
        p_out_TP = p_in_TP # Initial guess, can be improved with iteration
        zV_sol = zV(mass_flow_rate, p_in_TP)
        # print(f'Mass Flow Rate: {mass_flow_rate} kg/s, zB: {zB_sol + (L_heated/2)} m, zV: {zV_sol + (L_heated/2)} m}')
        if separation_point == "zB":
            x_in_TP = x_e_z(zB_sol, mass_flow_rate, p_in_TP) # [-]
            x_out_TP = x_e_z(zV_sol, mass_flow_rate, p_out_TP) # [-]
        else:
            #===Levy correlation to find x===
            x_in_TP = x_flow_z(zD_sol, mass_flow_rate, p_in_TP) # [-]
            x_out_TP = x_flow_z(zV_sol, mass_flow_rate, p_out_TP) # [-]

        # #===Acceleration pressure drop===
        dp_acc_TP = G_m_TP**2 * (1/rho_m(x_out_TP, p_out_TP) - 1/rho_m(x_in_TP, p_in_TP)) / 1e5  # [bar]

        #===Gravity pressure drop===
        if separation_point == "zB":
            dp_gravity_TP = quad(lambda z: rho_m(x_e_z(z, mass_flow_rate, p_in_TP), p_in_TP) * g, zB_sol, zV_sol)[0] / 1e5 # [bar]
        else:
            dp_gravity_TP = quad(lambda z: rho_m(x_flow_z(z, mass_flow_rate, p_in_TP), p_in_TP) * g, zD_sol, zV_sol)[0] / 1e5 # [bar]
        
        # #===Friction pressure drop (using McAdams correlation)===
        mu_LO = PropsSI('VISCOSITY', 'P', p_in_TP * 1e5, 'Q', 0, 'Water') # [Pa.s]
        Re_LO = Re(mass_flow_rate, mu_LO, A_flow, Hydraulic_diameter)
        friction_factor_L0 = McAdams_factor_1phase(Re_LO)  # Using the liquid only Reynolds number at the inlet of the two phase region
        #Have to change to correct Mcadams one phase and not colebrook 
        if separation_point == "zB":
            if heat_transfer_correlation_choice == "McAdams":
                dp_friction_TP = quad(lambda z: (friction_factor_L0 * (G_m_TP**2) / (Hydraulic_diameter * 2 * steamTable.rhoL_p(p_in_TP))) * phi_2_L0_McAdams(Re_LO, x_e_z(z, mass_flow_rate, p_in_TP), p_in_TP), zB_sol, zV_sol)[0] / 1e5 # [bar]
            else:
                dp_friction_TP = quad(lambda z: (friction_factor_L0 * (G_m_TP**2) / (Hydraulic_diameter * 2 * steamTable.rhoL_p(p_in_TP))) * phi_2_L0_Jones(G_m_TP, x_e_z(z, mass_flow_rate, p_in_TP), p_in_TP), zB_sol, zV_sol)[0] / 1e5 # [bar] with Jones
        else:
            if heat_transfer_correlation_choice == "McAdams":
                dp_friction_TP = quad(lambda z: (friction_factor_L0 * (G_m_TP**2) / (Hydraulic_diameter * 2 * steamTable.rhoL_p(p_in_TP))) * phi_2_L0_McAdams(Re_LO, x_flow_z(z, mass_flow_rate, p_in_TP), p_in_TP), zD_sol, zV_sol)[0] / 1e5 # [bar]
            else:
                dp_friction_TP = quad(lambda z: (friction_factor_L0 * (G_m_TP**2) / (Hydraulic_diameter * 2 * steamTable.rhoL_p(p_in_TP))) * phi_2_L0_Jones(G_m_TP, x_flow_z(z, mass_flow_rate, p_in_TP), p_in_TP), zD_sol, zV_sol)[0] / 1e5 # [bar] with Jones

    else:
        #===Iterative process to find outlet pressure with non-constant properties===
        p_guess = p_in_TP - 0.0  # initial guess for outlet pressure (bar)

        for it in range(50):
            zV_sol = zV(mass_flow_rate, p_guess)
            if separation_point == "zB":
                x_in_TP = x_e_z(zB_sol, mass_flow_rate, p_in_TP) # [-]
                x_out_TP = x_e_z(zV_sol, mass_flow_rate, p_guess) # [-]
            else:
                x_in_TP = x_flow_z(zD_sol, mass_flow_rate, p_in_TP) # [-]
                x_out_TP = x_flow_z(zV_sol, mass_flow_rate, p_guess) # [-]

            #===Acceleration pressure drop===
            dp_acc_TP = G_m_TP**2 * (1/rho_m(x_out_TP, p_guess) - 1/rho_m(x_in_TP, p_in_TP)) / 1e5  # [bar]

            p_mid_TP = (p_in_TP + p_guess) / 2  # Midpoint pressure for gravity and friction calculations
            #===Gravity pressure drop===
            if separation_point == "zB":
                dp_gravity_TP = quad(lambda z: rho_m(x_e_z(z, mass_flow_rate, p_mid_TP), p_mid_TP) * g, zB_sol, zV_sol)[0] / 1e5 # [bar]
            else:
                dp_gravity_TP = quad(lambda z: rho_m(x_flow_z(z, mass_flow_rate, p_mid_TP), p_mid_TP) * g, zD_sol, zV_sol)[0] / 1e5 # [bar]

            #===Friction pressure drop (using McAdams correlation)===
            mu_LO = PropsSI('VISCOSITY', 'P', p_mid_TP * 1e5, 'Q', 0, 'Water') # [Pa.s]
            Re_LO = Re(mass_flow_rate, mu_LO, A_flow, Hydraulic_diameter)
            friction_factor_L0 = McAdams_factor_1phase(Re_LO)  # Using the liquid only Reynolds number at the inlet of the two phase region

            if separation_point == "zB":
                if heat_transfer_correlation_choice == "McAdams":
                    dp_friction_TP = quad(lambda z: (friction_factor_L0 * (G_m_TP**2) / (Hydraulic_diameter * 2 * steamTable.rhoL_p(p_mid_TP))) * phi_2_L0_McAdams(Re_LO, x_e_z(z, mass_flow_rate, p_mid_TP), p_mid_TP), zB_sol, zV_sol)[0] / 1e5 # [bar]
                else:
                    dp_friction_TP = quad(lambda z: (friction_factor_L0 * (G_m_TP**2) / (Hydraulic_diameter * 2 * steamTable.rhoL_p(p_mid_TP))) * phi_2_L0_Jones(G_m_TP, x_e_z(z, mass_flow_rate, p_mid_TP), p_mid_TP), zB_sol, zV_sol)[0] / 1e5 # [bar] with Jones
            else:   
                if heat_transfer_correlation_choice == "McAdams":
                    dp_friction_TP = quad(lambda z: (friction_factor_L0 * (G_m_TP**2) / (Hydraulic_diameter * 2 * steamTable.rhoL_p(p_mid_TP))) * phi_2_L0_McAdams(Re_LO, x_flow_z(z, mass_flow_rate, p_mid_TP), p_mid_TP), zD_sol, zV_sol)[0] / 1e5 # [bar]
                else:
                    dp_friction_TP = quad(lambda z: (friction_factor_L0 * (G_m_TP**2) / (Hydraulic_diameter * 2 * steamTable.rhoL_p(p_mid_TP))) * phi_2_L0_Jones(G_m_TP, x_flow_z(z, mass_flow_rate, p_mid_TP), p_mid_TP), zD_sol, zV_sol)[0] / 1e5 # [bar] with Jones

            #===Total pressure drop===
            dp_total_TP = dp_acc_TP + dp_gravity_TP + dp_friction_TP   
            #===New outlet pressure===
            p_new = p_in_TP - dp_total_TP
            #===Convergence check===
            if abs(p_new - p_guess) < 1e-4:
                break
            p_guess = p_new #0.5*p_guess + 0.5*p_new # relaxation (stable)

    #===Store final values===

    dp_acc_list_two_phase.append(dp_acc_TP)
    dp_grav_list_two_phase.append(dp_gravity_TP)
    dp_fric_list_two_phase.append(dp_friction_TP)
    
    #===Total pressure drop===
    dp_total_TP = dp_acc_TP + dp_gravity_TP + dp_friction_TP
    dp_tot_list_two_phase.append(dp_total_TP)

#===Plot the two phase region pressure drop components===   
plt.plot(mass_flow_rate_list_final, dp_acc_list_two_phase, label='Acceleration Pressure Drop')
plt.plot(mass_flow_rate_list_final, dp_fric_list_two_phase, label='Friction Pressure Drop')
plt.plot(mass_flow_rate_list_final, dp_grav_list_two_phase, label='Gravity Pressure Drop')
plt.plot(mass_flow_rate_list_final, dp_tot_list_two_phase, label='Total Pressure Drop')
plt.xlabel('Mass Flow Rate [kg/s]')
plt.ylabel('Pressure Drop [bar]')
plt.title('Pressure Drop Components in Two Phase Region vs Mass Flow Rate')
plt.legend()
plt.grid()
plt.show()

#===Third graph: Vapour region (from zV to L_heated/2)===
#===Still have to change to change with constant properties or not===
dp_acc_list_vapour = []
dp_fric_list_vapour = []
dp_grav_list_vapour = []
dp_tot_list_vapour = []

for mass_flow_rate in mass_flow_rate_list_final:
    zV_sol = zV(mass_flow_rate, p_in)
    p_in_vapour = p_in - dp_tot_list_liquid[mass_flow_rate_list_final.index(mass_flow_rate)] - dp_tot_list_two_phase[mass_flow_rate_list_final.index(mass_flow_rate)]  # [bar]
    # print(f'For mass flow rate {mass_flow_rate} kg/s, pressure at vapour region inlet: {p_in_vapour} bar')
    if zV_sol == L_heated/2:
        # No vapour region
        dp_acc_list_vapour.append(0.0)
        dp_fric_list_vapour.append(0.0)
        dp_grav_list_vapour.append(0.0)
        dp_tot_list_vapour.append(0.0)
        continue    
    else:
        h_m_vapour_in = h_m_z(zV_sol , mass_flow_rate) # [kJ/kg]
        rho_vapour_in = steamTable.rho_ph(p_in_vapour, h_m_vapour_in) # [kg/m³]
        #mu_vapour_in = steamTable.my_ph(p_in_vapour, h_m_vapour_in) # [Pa.s]
        mu_vapour_in = PropsSI('VISCOSITY', 'P', p_in_vapour * 1e5, 'H', h_m_vapour_in * 1e3, 'Water') # [Pa.s]
        u_vapour_in = mass_flow_rate / (rho_vapour_in * A_flow)
        Re_vapour_in = rho_vapour_in * u_vapour_in * Hydraulic_diameter / mu_vapour_in
        L_vap = L_heated/2 - zV_sol
        

        if constant_properties_choice == "NO":
            # initial guess for outlet pressure (bar)
            p_guess = p_in_vapour - 0.0

            # iteration to get rho_out -> rho_mean -> dp -> p_out
            for it in range(50):
                # enthalpy at region outlet (fixed by heating)
                h_vap_out = h_m_z(L_heated/2, mass_flow_rate)   # [kJ/kg]
                # compute outlet properties at current p_guess
                # protect calls with try/except if steamTable can fail (optional)
                try:
                    rho_vapour_out = steamTable.rho_ph(p_guess, h_vap_out)
                    mu_vapour_out  = PropsSI('VISCOSITY', 'P', p_guess * 1e5, 'H', h_vap_out * 1e3, 'Water')
                except Exception:
                    # fallback: clamp to inlet props to avoid crash
                    rho_vapour_out = rho_vapour_in
                    mu_vapour_out  = mu_vapour_in

                # mean properties (user requested mean usage)
                rho_vapour_mean = (rho_vapour_in + rho_vapour_out)/2
                mu_vapour_mean  = (mu_vapour_in  + mu_vapour_out)/2

                # velocity and Reynolds based on mean properties
                u_vapour = mass_flow_rate / (rho_vapour_mean * A_flow)
                Re_vapour = rho_vapour_mean * u_vapour * Hydraulic_diameter / mu_vapour_mean

                # friction factor (uses your function)
                f_vapour = friction_factor_1phase(Re_vapour, relative_wall_roughness)

                # local pressure drop contributions using mean properties
                dp_acc_vapour = (mass_flow_rate / A_flow)**2 * (1.0 / rho_vapour_out - 1.0 / rho_vapour_in) / 1e5
                dp_friction_vapour = f_vapour * (L_vap / Hydraulic_diameter) * (rho_vapour_mean * u_vapour**2 / 2.0) / 1e5
                dp_gravity_vapour = rho_vapour_mean * g * L_vap / 1e5

                dp_total_iter = dp_acc_vapour + dp_friction_vapour + dp_gravity_vapour

                # updated outlet pressure (bar)
                p_new = p_in_vapour - dp_total_iter

                # convergence test on outlet pressure
                if abs(p_new - p_guess) < 1e-5:
                    p_guess = p_new
                    break
                
                p_guess = p_new #0.5 * p_guess + 0.5 * p_new # relaxation to stabilize

            dp_acc_vapour = dp_acc_vapour
            dp_friction_vapour = dp_friction_vapour
            dp_gravity_vapour = dp_gravity_vapour
            dp_total_vapour = dp_total_iter

        else:
            #===pressure drops===
            dp_acc_vapour = 0
            f_vapour = friction_factor_1phase(Re_vapour_in, relative_wall_roughness)
            dp_friction_vapour = f_vapour * (L_vap / Hydraulic_diameter) * (rho_vapour_in * u_vapour_in**2 / 2) / 1e5
            dp_gravity_vapour = rho_vapour_in * g * L_vap / 1e5
            dp_total_vapour = dp_acc_vapour + dp_friction_vapour + dp_gravity_vapour

        dp_acc_list_vapour.append(dp_acc_vapour)
        dp_fric_list_vapour.append(dp_friction_vapour)
        dp_grav_list_vapour.append(dp_gravity_vapour)
        dp_tot_list_vapour.append(dp_total_vapour)
#===Plot the vapour region pressure drop components===
plt.plot(mass_flow_rate_list_final, dp_acc_list_vapour, label='Acceleration Pressure Drop')
plt.plot(mass_flow_rate_list_final, dp_fric_list_vapour, label='Friction Pressure Drop')
plt.plot(mass_flow_rate_list_final, dp_grav_list_vapour, label='Gravity Pressure Drop')
plt.plot(mass_flow_rate_list_final, dp_tot_list_vapour, label='Total Pressure Drop')
plt.xlabel('Mass Flow Rate [kg/s]')
plt.ylabel('Pressure Drop [bar]')
plt.title('Pressure Drop Components in Vapour Region vs Mass Flow Rate')
plt.legend()
plt.grid()
plt.show()

#===Total pressure drop graph===
#===Sum over the three regions for each component===
dp_acc_list_total = [dp_acc_list_liquid[i] + dp_acc_list_two_phase[i] + dp_acc_list_vapour[i] for i in range(len(mass_flow_rate_list_final))]
dp_fric_list_total = [dp_fric_list_liquid[i] + dp_fric_list_two_phase[i] + dp_fric_list_vapour[i] for i in range(len(mass_flow_rate_list_final))]
dp_grav_list_total = [dp_grav_list_liquid[i] + dp_grav_list_two_phase[i] + dp_grav_list_vapour[i] for i in range(len(mass_flow_rate_list_final))]
dp_tot_list_total = [dp_tot_list_liquid[i] + dp_tot_list_two_phase[i] + dp_tot_list_vapour[i] for i in range(len(mass_flow_rate_list_final))]   
#===Plot the total pressure drop components===
plt.plot(mass_flow_rate_list_final, dp_acc_list_total, label='Acceleration Pressure Drop')
plt.plot(mass_flow_rate_list_final, dp_fric_list_total, label='Friction Pressure Drop')
plt.plot(mass_flow_rate_list_final, dp_grav_list_total, label='Gravity Pressure Drop')
plt.plot(mass_flow_rate_list_final, dp_tot_list_total, label='Total Pressure Drop')
plt.xlabel('Mass Flow Rate [kg/s]')
plt.ylabel('Pressure Drop [bar]')
plt.title('Total Pressure Drop Components vs Mass Flow Rate')
plt.legend()
plt.grid()
plt.show()  



#=== Plot zB and zV in terms of mass flow rate ===
zB_list = []
zV_list = []
for mass_flow_rate in mass_flow_rate_list_final:
    zB_sol = zB(mass_flow_rate)
    zV_sol = zV(mass_flow_rate, p_in) #TO IMPROVE: pass p after pressure drop calculation
    zB_list.append(zB_sol)
    zV_list.append(zV_sol)
plt.plot(mass_flow_rate_list_final, zB_list, label='Bulk Boiling Point zB')
plt.plot(mass_flow_rate_list_final, zV_list, label='Saturated Vapour Point zV')
# Add a horizontal dotted line at z = L_heated/2 and z = -L_heated/2
plt.axhline(y=L_heated/2, color='r', linestyle='--', label='Channel Outlet')
plt.axhline(y=-L_heated/2, color='g', linestyle='--', label='Channel Inlet')
plt.xlabel('Mass Flow Rate [kg/s]')
plt.ylabel('Axial Position [m]')
plt.title('Axial Positions zB and zV vs Mass Flow Rate')
plt.legend()
plt.grid()
plt.show()

#===plot zD in terms of mass flow rate===
zD_list = []
for mass_flow_rate in mass_flow_rate_list_final:
    zD_sol = zD(mass_flow_rate, p_in) #TO IMPROVE: pass p after pressure drop calculation
    zD_list.append(zD_sol)
plt.plot(mass_flow_rate_list_final, zD_list, label='Bubble Detachment Point zD')
# Add a horizontal dotted line at z = L_heated/2 and z = -L_heated/2
plt.axhline(y=L_heated/2, color='r', linestyle='--', label='Channel Outlet')
plt.axhline(y=-L_heated/2, color='g', linestyle='--', label='Channel Inlet')
plt.xlabel('Mass Flow Rate [kg/s]')
plt.ylabel('Axial Position [m]')
plt.title('Axial Position zD vs Mass Flow Rate')
plt.legend()
plt.grid()
plt.show()

# Plot equilibrium and flow quality in terms of height for a mass flow of 0,1 kg/s
#===Graphs of h and xe  along z for mass flow rate = 0.1 kg/s===

mass_flow_rate_graphs = 0.1 # [kg/s]
n_axial_G = 20
z_values_G = [i * L_heated / n_axial_G - (L_heated/2) for i in range(n_axial_G +1)] # [m]
h_values_G = [h_m_z(z, mass_flow_rate_graphs) for z in z_values_G]
plt.plot(h_values_G, z_values_G, label='Enthalpy along z')
plt.xlabel('Enthalpy [kJ/kg]')
plt.ylabel('Axial Position z [m]')
plt.title(f'Enthalpy  vs Axial Position z (Mass Flow Rate = {mass_flow_rate_graphs} kg/s)')
plt.legend()
plt.grid()
plt.show()

#=== In the same graph, plot also the flow quality along z between ZD and zV===
z_values_G_flow = [z for z in z_values_G if z >= zD(mass_flow_rate_graphs, p_in) and z <= zV(mass_flow_rate_graphs, p_in)]
x_e_values_G = [x_e_z(z, mass_flow_rate_graphs, p_in) for z in z_values_G]  
x_values_G = [x_flow_z(z, mass_flow_rate_graphs, p_in) for z in z_values_G_flow]
plt.plot(x_e_values_G, z_values_G, label='Equilibrium Quality')
plt.plot(x_values_G, z_values_G_flow, label='Flow Quality')
plt.xlabel(' Equilibrium Quality and Flow Quality [-]')
plt.ylabel('Axial Position z [m]')
plt.title(f'Equilibrium Quality and Flow Quality along Axial Position z (Mass Flow Rate = {mass_flow_rate_graphs} kg/s)')
plt.legend()
plt.grid()
plt.show()

#===Have to plot the graph of the flow quality along z too after using Levy correlation===