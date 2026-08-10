# -*- coding: utf-8 -*-
"""
Created on Sun Jan 21 14:23:36 2024

@author: jgsch
"""

import numpy as np
import matplotlib.pyplot as plt




h = 6.62607015e-27
c = 2.998e10
kb =  1.380649e-16
rsun_cm = 6.96e10
au2cm = 1.496e13
re = 637100000.0
nm = 1e-7


def chen_kipping(mass_array):
    rad = np.zeros(len(mass_array))
    for i in range(0, len(mass_array)):
        M = np.log10(mass_array[i])
        C1 = 0
        C2 = C1 + ((0.2790-0.589)*np.log10(2.04))
        C3 = C2 + ((0.589+0.044)*np.log10(318.0*0.414))
        if mass_array[i] <= 2.04:
            rad[i] = C1 + M*0.2790
        elif mass_array[i] > 2.04 and mass_array[i] <= 318.0*0.414:
            rad[i] = C2 + M*0.589
        else:
            rad[i] = C3 + M*(-0.044)
            
    return rad

def animate_through_time(pldata):
    pldata.insert(len(pldata.columns), 'rho', 5.512*pldata['pl_bmasse']/(pldata['pl_rade']**3))
    colors = {'Transit': 'mediumorchid',
              'Radial Velocity': 'mediumseagreen',
              'Microlensing': 'dodgerblue',
              'Imaging': 'goldenrod',
              'Transit Timing Variations': 'gray',
              'Eclipse Timing Variations': 'gray',
              'Orbital Brightness Modulation': 'gray',
              'Pulsar Timing': 'gray',
              'Astrometry': 'violet',
              'Pulsation Timing Variations':'gray',
              'Disk Kinematics': 'gray'}

    markers = {'Transit': 'o',
              'Radial Velocity': 's',
              'Microlensing': 'D',
              'Imaging': 'P',
              'Transit Timing Variations': '>',
              'Eclipse Timing Variations': '>',
              'Orbital Brightness Modulation': '>',
              'Pulsar Timing': '>',
              'Astrometry': '*',
              'Pulsation Timing Variations':'>',
              'Disk Kinematics': '>'}

    labels = {'Transit': 'Transit',
              'Radial Velocity': 'Radial Velocity',
              'Microlensing': 'Microlensing',
              'Imaging': 'Imaging',
              'Transit Timing Variations': 'other',
              'Eclipse Timing Variations': None,
              'Orbital Brightness Modulation': None,
              'Pulsar Timing': None,
              'Astrometry': 'Astrometry',
              'Pulsation Timing Variations':None,
              'Disk Kinematics': None}

    years = np.arange(1999, 2027)
    fig, ax = plt.subplots(2,2, figsize = (10,10))
    def animate(i):
        ax[0,0].cla()
        ax[0,1].cla()
        ax[1,0].cla()
        ax[1,1].cla()
        for disc in np.unique(pldata['discoverymethod']):
            dhold = pldata.loc[(pldata['discoverymethod'] == disc) & (pldata['disc_year'] <= years[i])]
            if len(dhold) > 0:
                if len(np.where(np.isnan(dhold['pl_rade']) == False)[0]) > 0:
                    ax[0,0].plot(dhold['pl_orbper'], dhold['pl_rade'], 
                               lw = 0, 
                               fillstyle = 'none', 
                               marker = markers[disc], 
                               color = colors[disc], 
                               label = labels[disc])
                    ax[0,0].legend(loc = 'lower right')
                if len(np.where(np.isnan(dhold['pl_rade']) == False)[0]) > 0: 
                    ax[0,1].plot(dhold['pl_orbper'], dhold['pl_bmasse'], 
                               lw = 0, 
                               fillstyle = 'none', 
                               marker = markers[disc], 
                               color = colors[disc])
                if len(np.where((np.isnan(dhold['pl_rade']) == False) & (np.isnan(dhold['pl_bmasse']) == False))[0]) > 0: 
                    ax[1,0].plot(dhold['pl_bmasse'], dhold['pl_rade'], 
                               lw = 0, 
                               fillstyle = 'none', 
                               marker = markers[disc], 
                               color = colors[disc])
                if len(np.where((np.isnan(dhold['rho']) == False) & (np.isnan(dhold['pl_bmasse']) == False))[0]) > 0: 
                    ax[1,1].plot(dhold['pl_bmasse'], dhold['rho'], 
                               lw = 0, 
                               fillstyle = 'none', 
                               marker = markers[disc], 
                               color = colors[disc])
                
        
        ax[0,0].set_yscale('log')
        ax[0,0].set_xscale('log')
        ax[0,0].set_ylim(0.2, 70)
        ax[0,0].set_xlim(0.01, 10**8)
        ax[0,0].set_ylabel(r'Radius [R$_\oplus$]', fontsize = 18)
        ax[0,0].set_xlabel(r'Orbital Period [days]', fontsize = 18)
        
        ax[0,1].set_yscale('log')
        ax[0,1].set_xscale('log')
        ax[0,1].set_ylim(0.1, 10000)
        ax[0,1].set_xlim(0.01, 10**8)
        ax[0,1].set_ylabel(r'Mass [M$_\oplus$]', fontsize = 18)
        ax[0,1].set_xlabel(r'Orbital Period [days]', fontsize = 18)
        
        ax[1,0].set_yscale('log')
        ax[1,0].set_xscale('log')
        ax[1,0].set_ylim(0.2, 70)
        ax[1,0].set_xlim(0.1, 10000)
        ax[1,0].set_xlabel(r'Mass [M$_\oplus$]', fontsize = 18)
        ax[1,0].set_ylabel(r'Radius [R$_\oplus$]', fontsize = 18)
        
        ax[1,1].set_yscale('log')
        ax[1,1].set_xscale('log')
        ax[1,1].set_ylim(0.2, 70)
        ax[1,1].set_xlim(0.1, 10000)
        ax[1,1].set_xlabel(r'Mass [M$_\oplus$]', fontsize = 18)
        ax[1,1].set_ylabel(r'Density [g/cc]', fontsize = 18)


        ax[0,0].annotate('Year = ' + str(years[i]), (10**4, 2), fontsize = 12)

                
        ax[0,0].tick_params(which = 'major', direction = 'in', top = True, right = True, length = 10, width = 2, 
                          labelsize = 16)
        ax[0,0].tick_params(which = 'minor', direction = 'in', top = True, right = True, length = 7.5, width = 1)

        ax[0,1].tick_params(which = 'major', direction = 'in', top = True, right = True, length = 10, width = 2, 
                          labelsize = 16)
        ax[0,1].tick_params(which = 'minor', direction = 'in', top = True, right = True, length = 7.5, width = 1)
        
        ax[1,0].tick_params(which = 'major', direction = 'in', top = True, right = True, length = 10, width = 2, 
                          labelsize = 16)
        ax[1,0].tick_params(which = 'minor', direction = 'in', top = True, right = True, length = 7.5, width = 1)
        
        ax[1,1].tick_params(which = 'major', direction = 'in', top = True, right = True, length = 10, width = 2, 
                          labelsize = 16)
        ax[1,1].tick_params(which = 'minor', direction = 'in', top = True, right = True, length = 7.5, width = 1)
        
        
        ax[0,0].patch.set_edgecolor('black')  
        ax[0,0].patch.set_linewidth(2) 
    
        ax[0,1].patch.set_edgecolor('black')  
        ax[0,1].patch.set_linewidth(2) 
        
        ax[1,0].patch.set_edgecolor('black')  
        ax[1,0].patch.set_linewidth(2) 
        
        ax[1,1].patch.set_edgecolor('black')  
        ax[1,1].patch.set_linewidth(2)
        
        
        plt.tight_layout()

    import matplotlib.animation

    ani = matplotlib.animation.FuncAnimation(fig, animate, frames=len(years))


    writer = matplotlib.animation.PillowWriter(fps=1,
                                     metadata=dict(artist='Me'))

    directory = './figures/'
    gifname = 'through_time'
    ani.save(directory + gifname + '.gif', writer=writer)
    
        
        
        
def plot_year(pldata, year):
    pldata.insert(len(pldata.columns), 'rho', 5.512*pldata['pl_bmasse']/(pldata['pl_rade']**3))
    colors = {'Transit': 'mediumorchid',
              'Radial Velocity': 'mediumseagreen',
              'Microlensing': 'dodgerblue',
              'Imaging': 'goldenrod',
              'Transit Timing Variations': 'gray',
              'Eclipse Timing Variations': 'gray',
              'Orbital Brightness Modulation': 'gray',
              'Pulsar Timing': 'gray',
              'Astrometry': 'violet',
              'Pulsation Timing Variations':'gray',
              'Disk Kinematics': 'gray'}

    markers = {'Transit': 'o',
              'Radial Velocity': 's',
              'Microlensing': 'D',
              'Imaging': 'P',
              'Transit Timing Variations': '>',
              'Eclipse Timing Variations': '>',
              'Orbital Brightness Modulation': '>',
              'Pulsar Timing': '>',
              'Astrometry': '*',
              'Pulsation Timing Variations':'>',
              'Disk Kinematics': '>'}

    labels = {'Transit': 'Transit',
              'Radial Velocity': 'Radial Velocity',
              'Microlensing': 'Microlensing',
              'Imaging': 'Imaging',
              'Transit Timing Variations': 'other',
              'Eclipse Timing Variations': None,
              'Orbital Brightness Modulation': None,
              'Pulsar Timing': None,
              'Astrometry': 'Astrometry',
              'Pulsation Timing Variations':None,
              'Disk Kinematics': None}


    fig, ax = plt.subplots(2,2, figsize = (10,10))

    for disc in np.unique(pldata['discoverymethod']):
        dhold = pldata.loc[(pldata['discoverymethod'] == disc) & (pldata['disc_year'] <= int(year))]
        
        if len(dhold) > 0:
            if len(np.where(np.isnan(dhold['pl_rade']) == False)[0]) > 0:
                ax[0,0].plot(dhold['pl_orbper'], dhold['pl_rade'], 
                           lw = 0, 
                           fillstyle = 'none', 
                           marker = markers[disc], 
                           color = colors[disc], 
                           label = labels[disc])
                ax[0,0].legend(loc = 'lower right')
            if len(np.where(np.isnan(dhold['pl_rade']) == False)[0]) > 0: 
                ax[0,1].plot(dhold['pl_orbper'], dhold['pl_bmasse'], 
                           lw = 0, 
                           fillstyle = 'none', 
                           marker = markers[disc], 
                           color = colors[disc])
            if len(np.where((np.isnan(dhold['pl_rade']) == False) & (np.isnan(dhold['pl_bmasse']) == False))[0]) > 0: 
                ax[1,0].plot(dhold['pl_bmasse'], dhold['pl_rade'], 
                           lw = 0, 
                           fillstyle = 'none', 
                           marker = markers[disc], 
                           color = colors[disc])
            if len(np.where((np.isnan(dhold['rho']) == False) & (np.isnan(dhold['pl_bmasse']) == False))[0]) > 0: 
                ax[1,1].plot(dhold['pl_bmasse'], dhold['rho'], 
                           lw = 0, 
                           fillstyle = 'none', 
                           marker = markers[disc], 
                           color = colors[disc])
        
            
        
    ax[0,0].set_yscale('log')
    ax[0,0].set_xscale('log')
    ax[0,0].set_ylim(0.2, 70)
    ax[0,0].set_xlim(0.01, 10**8)
    ax[0,0].set_ylabel(r'Radius [R$_\oplus$]', fontsize = 18)
    ax[0,0].set_xlabel(r'Orbital Period [days]', fontsize = 18)
    
    ax[0,1].set_yscale('log')
    ax[0,1].set_xscale('log')
    ax[0,1].set_ylim(0.1, 10000)
    ax[0,1].set_xlim(0.01, 10**8)
    ax[0,1].set_ylabel(r'Mass [M$_\oplus$]', fontsize = 18)
    ax[0,1].set_xlabel(r'Orbital Period [days]', fontsize = 18)
    
    ax[1,0].set_yscale('log')
    ax[1,0].set_xscale('log')
    ax[1,0].set_ylim(0.2, 70)
    ax[1,0].set_xlim(0.1, 10000)
    ax[1,0].set_xlabel(r'Mass [M$_\oplus$]', fontsize = 18)
    ax[1,0].set_ylabel(r'Radius [R$_\oplus$]', fontsize = 18)
    
    ax[1,1].set_yscale('log')
    ax[1,1].set_xscale('log')
    ax[1,1].set_ylim(0.2, 70)
    ax[1,1].set_xlim(0.1, 10000)
    ax[1,1].set_xlabel(r'Mass [M$_\oplus$]', fontsize = 18)
    ax[1,1].set_ylabel(r'Density [g/cc]', fontsize = 18)


    ax[0,0].annotate('Year = ' + str(year), (10**4, 2), fontsize = 12)

            
    ax[0,0].tick_params(which = 'major', direction = 'in', top = True, right = True, length = 10, width = 2, 
                      labelsize = 16)
    ax[0,0].tick_params(which = 'minor', direction = 'in', top = True, right = True, length = 7.5, width = 1)

    ax[0,1].tick_params(which = 'major', direction = 'in', top = True, right = True, length = 10, width = 2, 
                      labelsize = 16)
    ax[0,1].tick_params(which = 'minor', direction = 'in', top = True, right = True, length = 7.5, width = 1)
    
    ax[1,0].tick_params(which = 'major', direction = 'in', top = True, right = True, length = 10, width = 2, 
                      labelsize = 16)
    ax[1,0].tick_params(which = 'minor', direction = 'in', top = True, right = True, length = 7.5, width = 1)
    
    ax[1,1].tick_params(which = 'major', direction = 'in', top = True, right = True, length = 10, width = 2, 
                      labelsize = 16)
    ax[1,1].tick_params(which = 'minor', direction = 'in', top = True, right = True, length = 7.5, width = 1)
    
    
    ax[0,0].patch.set_edgecolor('black')  
    ax[0,0].patch.set_linewidth(2) 

    ax[0,1].patch.set_edgecolor('black')  
    ax[0,1].patch.set_linewidth(2) 
    
    ax[1,0].patch.set_edgecolor('black')  
    ax[1,0].patch.set_linewidth(2) 
    
    ax[1,1].patch.set_edgecolor('black')  
    ax[1,1].patch.set_linewidth(2)
    
    
    plt.tight_layout()
    plt.show()
        

        
    
    
    


        
    
    
    
