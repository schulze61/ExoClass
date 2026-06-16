# -*- coding: utf-8 -*-
"""
Created on Sun Jan  7 12:03:02 2024

@author: jgsch
"""

import numpy as np
import scipy.optimize as so
from spec_eos import *
import pandas as pd

#maybe add some sort of function for making initial guesses. Might speed it up.
def rhofind_Fe(P, EOS = Smith_sFe):
    rho = np.array([float(so.fsolve(EOS, 9000.0 , args = (P))) for P in P])
    return rho

def rhofind_mant(P, EOS = EP_mantle):
    if EOS == EP_mantle:
        rho = EP_mantle(P)
    else:
        rho = np.array([float(so.fsolve(EOS, 5000.0 , args = (P))) for P in P])
    return rho

def rhofind_h2o(P, EOS):
    rho = np.array([so.fsolve(EOS, 2000.0, args = (P))[0] for P in P])
    return rho

elike = pd.read_csv('./data_files/elike.csv')
elike = si.interp1d(elike['mass'], elike['density'], fill_value = 'extrapolate')

G = 6.6743*(10**(-11))
me = 5.972*(10**24)
re = 6371000.0


def guesses(mass):
    rho_guess = 1000*elike(mass)
    rad_guess = pow(mass*me/((4.0/3.0)*np.pi*rho_guess), 1.0/3.0)/re
    central_pressure_guess = ((rho_guess/(5510.0))**2)*((rad_guess**2))*364.0
    
    return central_pressure_guess, rad_guess



def initialize_planet(mass, cmf, core_eos = Anderson_lFe, mant_eos = EP_mantle, num_core_layers = 300, num_mant_layers = 300, get_shell_masses = False):
    if cmf == 0:
        num_core_layers = 0
    if cmf == 1:
        num_mant_layers = 0
        
    mass_array = np.zeros(num_core_layers + num_mant_layers)
    pressure_array = np.zeros(num_core_layers + num_mant_layers)
    radius_array = np.zeros(num_core_layers + num_mant_layers)
    density_array = np.zeros(num_core_layers + num_mant_layers)
    
    
    if cmf == 0:
        mass_array[:] = mass*(1.0-cmf)/num_mant_layers
    elif cmf == 1:
        mass_array[:] = mass*cmf/num_core_layers
    else:
        mass_array[:num_core_layers] = mass*cmf/num_core_layers
        mass_array[num_core_layers:] = mass*(1.0-cmf)/num_mant_layers
    
    
    cum_mass = np.array([sum(mass_array[:i+1]) for i in range(0, len(mass_array))])
    central_pressure_guess, Radius_planet_guess = guesses(mass)
    CMB_P_guess = (262.*(Radius_planet_guess)-550.*pow(Radius_planet_guess,2.) + 432.*pow(Radius_planet_guess,3.))

    
    if cmf == 0:
        dP_dr = CMB_P_guess/num_mant_layers
        pressure_array = np.linspace(num_mant_layers, 0, num_mant_layers)*dP_dr
        density_array = rhofind_mant(pressure_array, mant_eos)
    elif cmf == 1:
        dP_dr = (central_pressure_guess - CMB_P_guess)/(num_core_layers)
        pressure_array[:num_core_layers] = np.linspace(num_core_layers, 0, num_core_layers)*dP_dr + CMB_P_guess
        density_array = rhofind_Fe(pressure_array, core_eos)
    else:   
        dP_dr = (central_pressure_guess - CMB_P_guess)/(num_core_layers)
        pressure_array[:num_core_layers] = np.linspace(num_core_layers, 0, num_core_layers)*dP_dr + CMB_P_guess
    
        dP_dr = CMB_P_guess/num_mant_layers
        pressure_array[num_core_layers:] = np.linspace(num_mant_layers, 0, num_mant_layers)*dP_dr
    
        density_array[:num_core_layers] = rhofind_Fe(pressure_array[:num_core_layers], core_eos)
        density_array[num_core_layers:] = rhofind_mant(pressure_array[num_core_layers:], mant_eos)
        
    
    for i in range(1, len(radius_array)):
        radius_array[i] = pow((3.0*(cum_mass[i] - cum_mass[i-1])*me/(4.0*np.pi*(density_array[i] + density_array[i-1])/2.0)) + (radius_array[i-1])**3, 1.0/3.0)
    
    if get_shell_masses:
      return pressure_array, density_array, radius_array, cum_mass, mass_array
    else:  
        return pressure_array, density_array, radius_array, cum_mass


def build_rocky_planet(mass, cmf, core_eos = Anderson_lFe, mant_eos = EP_mantle, num_core_layers = 600, num_mant_layers = 500, get_shell_masses = False):
    if get_shell_masses:
        pressure_array, density_array, radius_array, mass_array, ms = initialize_planet(mass, cmf, core_eos, mant_eos, num_core_layers, num_mant_layers, get_shell_masses = True)
    else:
        pressure_array, density_array, radius_array, mass_array = initialize_planet(mass, cmf, core_eos, mant_eos, num_core_layers, num_mant_layers)
    
    pressure_array[-1] = 101325.0
    dold = 100000*density_array
    
    check = False
    icount = 0
    
    #while check == False and icount < 100:
    while check == False:

        for i in range(1, len(pressure_array))[::-1]:
            pressure_array[i-1] = pressure_array[i] + (G*mass_array[i]*me*density_array[i]/(radius_array[i]**2))*np.diff(radius_array)[i-1]
        density_array = np.zeros(len(dold))
        if cmf == 0:
            density_array = rhofind_mant(pressure_array/1e9, mant_eos)
        elif cmf == 1:
            density_array = rhofind_Fe(pressure_array/1e9, core_eos)
        else:
            density_array[:num_core_layers] = rhofind_Fe(pressure_array[:num_core_layers]/1e9, core_eos)
            density_array[num_core_layers:] = rhofind_mant(pressure_array[num_core_layers:]/1e9, mant_eos)
        
        for i in range(1, len(radius_array)):
            radius_array[i] = pow((3.0*(mass_array[i] - mass_array[i-1])*me/(4.0*np.pi*(density_array[i] + density_array[i-1])/2.0)) + (radius_array[i-1])**3, 1.0/3.0)
        
        check = check_converge(dold, density_array, icount)
        dold = density_array
        
        icount += 1
        
        #if icount == 30:
        #    print('skipped')
        #    if get_shell_masses:
        #        return np.zeros(1), np.zeros(1), np.zeros(1), np.zeros(1), np.zeros(1), np.zeros(1), np.zeros(1)
        #    else:
        #        return np.zeros(1), np.zeros(1), np.zeros(1), np.zeros(1)
    
    #print(icount)
    
    if get_shell_masses:
        return pressure_array, density_array, radius_array, mass_array, ms, num_core_layers, num_mant_layers
    else:
        return pressure_array, density_array, radius_array, mass_array
    



def check_converge(dold, dnew, icount = 0):
    t1 = abs(1.0 - dnew/dold)
    ci = np.where(t1 > 10**(-3))[0]
    check = len(np.where(t1 > 10**(-3))[0])
    if icount <= 15:
        if check > 0:
            return False
        else:
            return True
    else:
        #print(t1[ci], ci)
        if check > 2:
            return False
        else:
            return True

    
    
    
    
def initialize_planet_wat(mass, cmf, core_eos = Smith_sFe, mant_eos = EP_mantle, wt_frac_water = 0.1, h2o_eos = EP_h2o, num_core_layers = 600, num_mant_layers = 500, num_water_layers = 700):
        
    mass_array = np.zeros(num_core_layers + num_mant_layers + num_water_layers)
    pressure_array = np.zeros(num_core_layers + num_mant_layers + num_water_layers)
    radius_array = np.zeros(num_core_layers + num_mant_layers + num_water_layers)
    density_array = np.zeros(num_core_layers + num_mant_layers + num_water_layers)
    

    mass_array[:num_core_layers] = mass*cmf*(1-wt_frac_water)/num_core_layers
    mass_array[num_core_layers:num_mant_layers+num_core_layers] = mass*(1.0-cmf)*(1-wt_frac_water)/num_mant_layers
    mass_array[num_core_layers+num_mant_layers:] = mass*(wt_frac_water)/num_water_layers
    
    
    cum_mass = np.array([sum(mass_array[:i]) for i in range(0, len(mass_array))])
    central_pressure_guess, Radius_planet_guess = guesses(mass)
    CMB_P_guess = (262.*(Radius_planet_guess)-550.*pow(Radius_planet_guess,2.) + 432.*pow(Radius_planet_guess,3.))

  
    dP_dr = (central_pressure_guess - CMB_P_guess)/(num_core_layers)
    pressure_array[:num_core_layers] = np.linspace(num_core_layers, 0, num_core_layers)*dP_dr + CMB_P_guess
    
    #return pressure_array
    
    WMBP = 1.0
    dP_dr = (CMB_P_guess-WMBP)/num_mant_layers
    pressure_array[num_core_layers:num_mant_layers+num_core_layers] = pressure_array[num_core_layers-1] - np.arange(num_mant_layers)*dP_dr 
    
    
    dP_dr = WMBP/num_water_layers
    pressure_array[num_core_layers+num_mant_layers:] = pressure_array[num_core_layers + num_mant_layers - 1] - np.arange(num_water_layers)*dP_dr

    
    density_array[:num_core_layers] = rhofind_Fe(pressure_array[:num_core_layers], core_eos)
    density_array[num_core_layers:num_mant_layers+num_core_layers] = rhofind_mant(pressure_array[num_core_layers:num_mant_layers+num_core_layers], mant_eos)
    density_array[num_mant_layers+num_core_layers:] = rhofind_h2o(pressure_array[num_mant_layers+num_core_layers:], h2o_eos)
        
    
    for i in range(1, len(radius_array)):
        radius_array[i] = pow((3.0*(cum_mass[i] - cum_mass[i-1])*me/(4.0*np.pi*(density_array[i] + density_array[i-1])/2.0)) + (radius_array[i-1])**3, 1.0/3.0)
    
    
    
    return pressure_array, density_array, radius_array, cum_mass


def build_wet_planet(mass, cmf, wt_frac_water, core_eos, mant_eos, h2o_eos, num_core_layers = 600, num_mant_layers = 500, num_water_layers = 700):
    pressure_array, density_array, radius_array, mass_array = initialize_planet_wat(mass, cmf, core_eos, mant_eos, wt_frac_water, h2o_eos, num_core_layers, num_mant_layers, num_water_layers)
    
    pressure_array[-1] = 101325.0
    dold = 100000*density_array
    
    check = False
    
    while check == False:
        for i in range(1, len(pressure_array))[::-1]:
            pressure_array[i-1] = pressure_array[i] + (G*mass_array[i]*me*density_array[i]/(radius_array[i]**2))*np.diff(radius_array)[i-1]

        density_array[:num_core_layers] = rhofind_Fe(pressure_array[:num_core_layers]/1e9, core_eos)
        density_array[num_core_layers:num_mant_layers+num_core_layers] = rhofind_mant(pressure_array[num_core_layers:num_mant_layers+num_core_layers]/1e9, mant_eos)
        density_array[num_core_layers + num_mant_layers:] = rhofind_h2o(pressure_array[num_core_layers+num_mant_layers:]/1e9, h2o_eos)
        
        for i in range(1, len(radius_array)):
            radius_array[i] = pow((3.0*(mass_array[i] - mass_array[i-1])*me/(4.0*np.pi*(density_array[i] + density_array[i-1])/2.0)) + (radius_array[i-1])**3, 1.0/3.0)
    
        
        check = check_converge(dold, density_array)
        
        dold = density_array
        
    return pressure_array, density_array, radius_array, mass_array


    




    
