# -*- coding: utf-8 -*-
"""
Created on Sun Jan 21 14:23:36 2024

@author: jgsch
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings("ignore")
 
h = 6.62607015e-34 #plancks constant SI
c = 2.998e8 #speed of light SI
kb =  1.380649e-23 # Boltzman constant SI
rsun = 6.957e8 # Solar radius in SI
au2m = 1.496e11 # au to meters
re = 6371000.0 # radius of Earth SI
nm2m = 1e-9 #nanometers to meters
sol_2_earth_mass = 333060.402
lyr_2_au = 63241.1
pc_2_lyr = 3.26156
arcsec_2_rad = (10**(-6))*0.000004848137


def thermal_intensity(params, star = True):
    if star:
        R = params['star radius']*rsun
        T = params['star temperature']
    else:
        R = params['planet radius']*re
        T = equilibrium_temperature(params)
    lam = params['wavelength']*nm2m
    numerator = (2.0*h*(c**2))
    denom = (lam**5)*(np.exp((h*c)/(lam*kb*T)) - 1.0)
    
    return 4.0*np.pi*(R**2)*numerator/denom

def equilibrium_temperature(params):
    Tstar = params['star temperature']
    Rstar = params['star radius']*rsun
    a = params['planet semimajor axis']*au2m
    A = params['planet albedo']
    return Tstar*((1-A)**(1.0/4.0))*np.sqrt(Rstar/(2.0*a))


def thermal_contrast(params):
    return thermal_intensity(params, star = False)/thermal_intensity(params)

def reflected_contrast(params):
    Rplanet = params['planet radius']*re
    a = params['planet semimajor axis']*au2m
    return params['planet albedo']*(Rplanet**2)/(4.0*(a**2))

def min_planet_mass_astrometry(Mstar, a, d, theta_min = 400):
    return theta_min*arcsec_2_rad*Mstar*sol_2_earth_mass*((d*lyr_2_au)/a)


def direct_imaging_plot(params):
    pldata = pd.read_csv('./data/PS_2026.07.24_07.42.56.csv')
    pldata['sy_dist_lyr'] = pldata['sy_dist']*3.2615637769

    #double check this
    params['planet wavelength peak'] = 0.002897/equilibrium_temperature(params)
    params['theta min'] = 1.22*params['planet wavelength peak']/params['telescope diameter']

    from matplotlib.ticker import MultipleLocator
    fig, ax = plt.subplots(2,2, figsize = (10,10))

    #-----------------top left panel----------------------#
    dist_array = np.linspace(0.0, 100, 100)
    ax[0,0].fill_between(dist_array, dist_array*lyr_2_au*params['theta min'], alpha = 0.3, color = 'gray')
    ax[0,0].set_xlabel('Distance to Star (lyr)', fontsize = 24)
    ax[0,0].set_ylabel('Semi-major axis (au)', fontsize = 24)
    ax[0,0].plot(params['distance to star'], params['planet semimajor axis'], 'ro', markersize = 10)
    ax[0,0].grid( linestyle = '--', linewidth = 0.5, zorder = 100)

    ax[0,0].set_xlim(0, max(dist_array))
    ax[0,0].set_ylim(0, max(dist_array*lyr_2_au*params['theta min']))

    ax[0,0].xaxis.set_major_locator(MultipleLocator(10))
    ax[0,0].xaxis.set_minor_locator(MultipleLocator(2))
    ax[0,0].yaxis.set_major_locator(MultipleLocator(2.))
    ax[0,0].yaxis.set_minor_locator(MultipleLocator(0.5))

    ax[0,0].tick_params(which = 'major', direction = 'in', top = True, right = True, length = 10, width = 2, 
                      labelsize = 16)
    ax[0,0].tick_params(which = 'minor', direction = 'in', top = True, right = True, length = 7.5, width = 1)

    ax[0,0].annotate('Not resolvable' + '\n' + 'from host star', 
                     (max(dist_array)/2, max(dist_array*lyr_2_au*params['theta min'])/4), 
                     fontsize = 12)

    ax[0,0].patch.set_edgecolor('black')  
    ax[0,0].patch.set_linewidth(2) 
    #----------------------------------------------------#

    #-----------------top right panel--------------------#
    ax[0,1].hist(pldata['sy_dist_lyr'], bins = 30, range = (0, 100), cumulative = True, alpha = 0.3, color = 'dodgerblue')
    ax[0,1].axvline(params['distance to star'], color = 'r', lw = 2)
    ax[0,1].set_xlabel(r'Distance to Star (lyr)', fontsize = 24)
    ax[0,1].set_ylabel(r'N systems w/in d', fontsize = 24)

    ax[0,1].xaxis.set_major_locator(MultipleLocator(10))
    ax[0,1].xaxis.set_minor_locator(MultipleLocator(1))
    ax[0,1].yaxis.set_major_locator(MultipleLocator(50))
    ax[0,1].yaxis.set_minor_locator(MultipleLocator(10))

    ax[0,1].tick_params(which = 'major', direction = 'in', top = True, right = True, length = 10, width = 2, 
                      labelsize = 16)
    ax[0,1].tick_params(which = 'minor', direction = 'in', top = True, right = True, length = 7.5, width = 1)

    ax[0,1].patch.set_edgecolor('black')  
    ax[0,1].patch.set_linewidth(2) 
    ax[0,1].grid( linestyle = '--', linewidth = 0.5)

    ax[0,1].set_xlim(0, 100)
    ax[0,1].set_ylim(0, 450)
    #----------------------------------------------------#

    #-----------------bottom left panel------------------#
    params['wavelength'] = (10**np.linspace(-1, 2, 1000))*1000

    thermal_star = thermal_intensity(params)/max(thermal_intensity(params))
    thermal_planet = thermal_intensity(params, star = False)/max(thermal_intensity(params))
    reflected = reflected_contrast(params)*thermal_intensity(params)/max(thermal_intensity(params))

    ax[1,0].plot(params['wavelength']/1000.0, np.log10(thermal_star), 
                 '-', color = 'darkorange', lw = 3, label = 'Star')

    ax[1,0].plot(params['wavelength']/1000.0, np.log10(reflected+thermal_planet), 
                 ls = (0, (1, 1)), color = 'k', lw = 3, label = 'Planet', zorder = 10)

    ax[1,0].plot(params['wavelength']/1000.0, np.log10(thermal_planet), 
                 '-', color = 'olivedrab', lw = 2, label = 'Emitted')
    ax[1,0].plot(params['wavelength']/1000.0, np.log10(reflected), 
                 '-', color = 'mediumpurple', lw = 2, label = 'Reflected')

    ax[1,0].set_ylim(-14.9, 2)
    ax[1,0].set_xlim(10**(-1), 100)

    ax[1,0].set_xscale('log')
    ax[1,0].legend(frameon = False)

    ax[1,0].tick_params(which = 'major', direction = 'in', top = True, right = True, length = 10, width = 2, 
                      labelsize = 16)
    ax[1,0].tick_params(which = 'minor', direction = 'in', top = True, right = True, length = 7.5, width = 1)

    ax[1,0].patch.set_edgecolor('black')  
    ax[1,0].patch.set_linewidth(2) 

    ax[1,0].yaxis.set_major_locator(MultipleLocator(5))
    ax[1,0].yaxis.set_minor_locator(MultipleLocator(1))

    ax[1,0].set_xlabel(r'Wavelength ($\mu$m)', fontsize = 24)
    ax[1,0].set_ylabel(r'Normalized Intensity', fontsize = 24)
    #----------------------------------------------------#

    #-----------------bottom right panel-----------------#
    ax[1,1].plot(params['wavelength']/1000.0, thermal_contrast(params),
                '-', color = 'olivedrab', lw = 2, label = 'Emitted')
    ax[1,1].plot(params['wavelength']/1000.0, np.ones(len(params['wavelength']/1000.0))*reflected_contrast(params), 
                 '-', color = 'mediumpurple', lw = 2, label = 'Reflected')
    ax[1,1].plot(params['wavelength']/1000.0, np.ones(len(params['wavelength']/1000.0))*reflected_contrast(params) 
               + thermal_contrast(params), 
                ls = (0, (1, 1)), color = 'k', lw = 3, label = 'Both', zorder = 10)

    ax[1,1].legend(frameon = False, fontsize = 12)
    ax[1,1].fill_between([min(params['wavelength']/1000.0), max(params['wavelength']/1000.0)], 
                         [10**(-7), 10**(-7)], 
                         color = 'gray', 
                         alpha = 0.3)
    ax[1,1].set_xlim([min(params['wavelength']/1000.0), max(params['wavelength']/1000.0)])

    ax[1,1].set_ylim(10**(-8), 1.1*max([max(thermal_contrast(params)), reflected_contrast(params)]))
    ax[1,1].set_xlabel(r'$\lambda$ (nm)', fontsize = 24)
    ax[1,1].set_ylabel(r'Star-planet contrast', fontsize = 24)
    ax[1,1].set_xscale('log')
    ax[1,1].set_yscale('log')

    ax[1,1].set_xscale('log')
    ax[1,1].legend(frameon = False)

    ax[1,1].tick_params(which = 'major', direction = 'in', top = True, right = True, length = 10, width = 2, 
                      labelsize = 16)
    ax[1,1].tick_params(which = 'minor', direction = 'in', top = True, right = True, length = 7.5, width = 1)

    ax[1,1].patch.set_edgecolor('black')  
    ax[1,1].patch.set_linewidth(2) 
    ax[1,1].annotate('Not detectable',
                     (0.2, 10**(-7.3)),
                     fontsize = 12)

    ax[1,1].set_ylabel(r'Star-planet contrast', fontsize = 24)
    ax[1,1].set_xlabel(r'Wavelength ($\mu$m)', fontsize = 24)


    plt.tight_layout()
    
    
    
    
    
