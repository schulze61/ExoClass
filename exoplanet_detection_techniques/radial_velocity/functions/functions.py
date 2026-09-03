# -*- coding: utf-8 -*-
"""
Created on Thu Jan  4 10:49:44 2024

@author: jgsch
"""
import numpy as np
from scipy.optimize import newton
import matplotlib.pyplot as plt
import scipy.stats as sp
import pandas as pd


#some constants
G = 6.67430*(10**(-11)) #SI
Msun = 2*(10**30) #SI
au2m = 1.496*(10**11)
Me = 5.972*(10**24)
c = 2.998e+8 #speed of light in m/s
halpha_rest = 6562.8 #angstrom


def transcendental_equ(E, M, e):
    return M - (E - e*np.sin(E))

def build_rv_curve(Mp = 1.0, a = 1.0, i = np.pi/2.0, e = 0.0, V0 = 0.0, w = 0.0,  T = 0.5, numpts = 200, n_orbits = 1):
    a = a*au2m
    P = 2.0*np.pi*np.sqrt((a**3)/(G*Msun))
    t = P*np.linspace(0, n_orbits, numpts)
    T = T*P
    Mp = Mp*Me
    K = ((2.0*np.pi*G/P)**(1.0/3.0))*(Mp*np.sin(i)/((Msun+Mp)**(2.0/3.0)))*(1.0/np.sqrt(1.0-(e**2)))
    #K = ((2.0*np.pi*G/P)**(1.0/3.0))*(Mp*np.sin(i)/((Msun)**(2.0/3.0)))*(1.0/np.sqrt(1.0-(e**2)))

    
    M = 2.0*np.pi*(t-T)/P
    E = np.array([newton(transcendental_equ, m-e, args = (m, e)) for m in M])
    if e!=0:
        beta = (1.0 - np.sqrt(1.0 - (e**2)))/e
        v = E + 2.0*np.arctan((beta*np.sin(E))/(1.0 - beta*np.cos(E)))
    else:
        v = np.arccos((np.cos(E) - e) / (1 - e*np.cos(E)))
        
    V = V0 + K*(np.cos(v+w) + e*np.cos(w))
    
    dlambda = np.sqrt((1.0 + V/c) / (1.0 - V/c)) - 1.0
    halpha_obs = dlambda*halpha_rest + halpha_rest
    
    return t/P, V, E, halpha_obs, K


def animate_rv_curve(Mp = 333030/2.0, a = 0.05, e = 0.8, inc = np.pi/2.0, w = 0, numpts = 100,
                     arrow = False, T=0.5,
                     directory = './animations/',
                     static = False,
                     orb_phase = 0.5):
    phase, vel, E, ha_obs, K = build_rv_curve(Mp = Mp, a = a, i = inc, e = e, w = w, numpts = numpts, T=T)
    vel = vel/1000
    c = 299792
    spect = pd.read_csv('./data/mock_spectra.csv')

    b = a*np.sqrt(1.0 - e**2) #impact parameter
    focci = e*a #ellipse focus for plotting purposes
    x = (a*np.cos(E) - focci) 
    y = b*np.sin(E)

    astar = a*(Mp*Me)/Msun
    bstar = astar*np.sqrt(1.0 - e**2)
    xstar = (astar*np.cos(E + np.pi) + e*astar)
    ystar =  bstar*np.sin(E + np.pi)


    xmax = max([max(abs(xstar)), max(abs(x))])
    xmax = xmax + 0.05*xmax

    ymax = max([max(abs(ystar)), max(abs(y))])
    ymax = ymax + 0.05*ymax

    pause = min(abs(vel))/abs(vel)

    fig, hold = plt.subplots(1,1,figsize = (10,10))
    hold.set_axis_off()

    ax = plt.subplot2grid((4,2), loc = (0,0), colspan = 1, rowspan = 1)
    ax4 = plt.subplot2grid((4,2), loc = (1,0), colspan = 2, rowspan = 1)
    ax2 = plt.subplot2grid((4,2), loc = (0,1), colspan = 1, rowspan = 1)
    ax3 = plt.subplot2grid((4,2), loc = (2,0), colspan = 2, rowspan = 2)

    #ax.set_axis_off()

    ax.plot(x*np.cos(inc - np.pi/2.0),y, 'r-')
    ax.plot(xstar*np.cos(inc - np.pi/2.0), ystar, 'c-')

    ax2.plot(xstar*np.cos(inc - np.pi/2.0), -xstar*np.sin(inc - np.pi/2.0), 'c-')
    ax2.plot(x*np.cos(inc - np.pi/2.0), -x*np.sin(inc - np.pi/2.0), 'r-', alpha = 1.0)

    ax.grid(color='lightgray',linestyle='--');
    ax2.grid(color='lightgray',linestyle='--');
    ax3.grid(color='lightgray',linestyle='--');
    ax4.grid(color='lightgray',linestyle='--', which = 'both', axis = 'x');

    from matplotlib.ticker import MultipleLocator

    if xmax > ymax:
        ax.set_xlim(-xmax, xmax)
        ax.set_ylim(-xmax, xmax)
        if arrow:
            dx = (xmax/6)*np.sin(w+np.pi)
            dy = (xmax/6)*np.cos(w+np.pi)
            #ax.arrow(0, 0, dx = dx, dy = dy, 
            #         width = 0.001*xmax, head_width = 0.01*xmax, zorder = 0, facecolor = 'k')
            ax.annotate("", xytext=(0, 0), xy=(dx, dy), size=25, color = 'k',
                        arrowprops=dict(arrowstyle="simple"))
            ax.plot([-1000, 1000], [-1000*dy/dx, 1000*dy/dx], 'k--')
        ax2.set_xlim(-xmax, xmax)
        ax2.set_ylim(-xmax, xmax)
        xmaj = 2.0*xmax/3
        
    else:
        ax.set_xlim(-ymax, ymax)
        ax.set_ylim(-ymax, ymax)
        if arrow:
            dx = (ymax/6)*np.sin(w+np.pi)
            dy = (ymax/6)*np.cos(w+np.pi)
            ax.arrow(0, 0, dx = dx, dy = dy, 
                     width = 0.001, head_width = 0.01, zorder=0, facecolor = 'k')
            ax.plot([-1000, 1000], [-1000*dy/dx, 1000*dy/dx], 'k--')
        ax2.set_xlim(-ymax, ymax)
        ax2.set_ylim(-ymax, ymax)
        xmaj = 2.0*ymax/3

    xmaj = round(xmaj, -int(np.floor(np.log10(abs(xmaj)))))
    ax.xaxis.set_major_locator(MultipleLocator(xmaj))
    ax.xaxis.set_minor_locator(MultipleLocator(xmaj/5))
    ax.yaxis.set_major_locator(MultipleLocator(xmaj))
    ax.yaxis.set_minor_locator(MultipleLocator(xmaj/5))
    ax2.xaxis.set_major_locator(MultipleLocator(xmaj))
    ax2.xaxis.set_minor_locator(MultipleLocator(xmaj/5))
    ax2.yaxis.set_major_locator(MultipleLocator(xmaj))
    ax2.yaxis.set_minor_locator(MultipleLocator(xmaj/5))
    ax3.xaxis.set_major_locator(MultipleLocator(0.2))
    ax3.xaxis.set_minor_locator(MultipleLocator(0.05))
    xmaj = round(((max(vel)-min(vel))/5), 
                    -int(np.floor(np.log10(abs(((max(vel)-min(vel))/5))))))
    ax3.yaxis.set_major_locator(MultipleLocator(xmaj))
    ax3.yaxis.set_minor_locator(MultipleLocator(xmaj/5))
    ax3.set_xlim(0,1)


    ax.set_aspect('equal', 'box')
    ax2.set_aspect('equal', 'box') 

    ax3.plot(phase, vel, 'c-', alpha = 1.0);
    ax3.set_ylim(-1.1*max(abs(vel)), 1.1*max(abs(vel)))

    ax3.annotate(r'Max $v_r$ = ' + str(round(max(vel), 2)) + ' km/s', (0.15, 0.45), xycoords = 'figure fraction', fontsize = 14)
    ax3.annotate(r'Min $v_r$ = ' + str(round(min(vel), 2)) + ' km/s', (0.15, 0.425), xycoords = 'figure fraction', fontsize = 14)
    ax3.annotate(r'$K$ = ' + str(round(K/1000, 2)) + ' km/s', (0.15, 0.4), xycoords = 'figure fraction', fontsize = 14)

    ax4.plot(spect['lambda_rest'], spect['flux'], 'k-', alpha = 0.5)

    ax.set_xlabel('x [au]', fontsize = 24)
    ax.set_ylabel('y [au]', fontsize = 24)

    ax2.set_xlabel('x [au]', fontsize = 24)
    ax2.set_ylabel('z [au]', fontsize = 24)


    ax3.set_ylabel('Radial Velocity [km/s]', fontsize = 24)
    ax3.set_xlabel('Phase', fontsize = 24)


    ax.tick_params(which = 'major', direction = 'in', top = True, right = True, length = 10, width = 2, labelsize = 12)
    ax.tick_params(which = 'minor', direction = 'in', top = True, right = True, length = 7.5, width = 1)
    ax2.tick_params(which = 'major', direction = 'in', top = True, right = True, length = 10, width = 2, labelsize = 12)
    ax2.tick_params(which = 'minor', direction = 'in', top = True, right = True, length = 7.5, width = 1)
    ax3.tick_params(which = 'major', direction = 'in', top = True, right = True, length = 10, width = 2, labelsize = 16)
    ax3.tick_params(which = 'minor', direction = 'in', top = True, right = True, length = 7.5, width = 1)


    from matplotlib.animation import FuncAnimation
    import matplotlib.animation

    l, = ax.plot([],[], 'ko', markersize = 10);
    lstar, = ax.plot([],[], 'k*', markersize = 20);

    l2, = ax3.plot([],[], 'k*', markersize = 20);
    l3, = ax2.plot([],[], 'k*', markersize = 20);
    l4, = ax2.plot([],[], 'ko', markersize = 10);

    l5, = ax4.plot([],[], 'k-', markersize = 10);

    ax.patch.set_edgecolor('black')  
    ax.patch.set_linewidth(2) 
    ax2.patch.set_edgecolor('black')  
    ax2.patch.set_linewidth(2)
    ax3.patch.set_edgecolor('black')  
    ax3.patch.set_linewidth(2)

    ax4.set_xlim(4850, 4880)
    ax4.set_ylim(1.25, 2.0)


    plt.tight_layout()

    def animate(i):
        l.set_data([x[i]*np.cos(inc - np.pi/2.0)], [y[i]])
        lstar.set_data([xstar[i]*np.cos(inc - np.pi/2.0)], [ystar[i]])
        l2.set_data([phase[i]], [vel[i]])
        l3.set_data([xstar[i]*np.cos(inc - np.pi/2.0)], [-xstar[i]*np.sin(inc - np.pi/2.0)])
        l4.set_data([x[i]*np.cos(inc - np.pi/2.0)], [-x[i]*np.sin(inc - np.pi/2.0)])
        ax4.clear()
        ax4.plot(spect['lambda_rest'], spect['flux'], 'k-', alpha = 0.5)
        ax4.plot(spect['lambda_rest']*np.sqrt((c+vel[i])/(c-vel[i])), spect['flux'], 'c-')
        ax4.grid(color='lightgray',linestyle='--', which = 'both', axis = 'x');
        ax4.set_xlim(4850, 4880)
        ax4.set_ylim(1.25, 2.0)

        ax4.set_xlabel(r'Observed Wavelength [$\AA$]', fontsize = 24)
        ax4.set_ylabel(r'Flux', fontsize = 24)
        ax4.xaxis.set_major_locator(MultipleLocator(5))
        ax4.xaxis.set_minor_locator(MultipleLocator(1))
        ax4.tick_params(which = 'major', direction = 'in', top = True, right = True, length = 10, width = 2, labelsize = 16)
        ax4.tick_params(which = 'minor', direction = 'in', top = True, right = True, length = 7.5, width = 1)
        ax4.patch.set_edgecolor('black')  
        ax4.patch.set_linewidth(2)
        
        ax4.set_yticks([])
        ax4.get_xaxis().get_major_formatter().set_useOffset(False)
        plt.pause(pause[i])
        

    if static:
        i = np.argmin(abs(orb_phase-phase))
        l.set_data([x[i]*np.cos(inc - np.pi/2.0)], [y[i]])
        lstar.set_data([xstar[i]*np.cos(inc - np.pi/2.0)], [ystar[i]])
        l2.set_data([phase[i]], [vel[i]])
        l3.set_data([xstar[i]*np.cos(inc - np.pi/2.0)], [-xstar[i]*np.sin(inc - np.pi/2.0)])
        l4.set_data([x[i]*np.cos(inc - np.pi/2.0)], [-x[i]*np.sin(inc - np.pi/2.0)])
        ax4.clear()
        ax4.plot(spect['lambda_rest'], spect['flux'], 'k-', alpha = 0.5)
        ax4.plot(spect['lambda_rest']*np.sqrt((c+vel[i])/(c-vel[i])), spect['flux'], 'c-')
        ax4.grid(color='lightgray',linestyle='--', which = 'both', axis = 'x');
        ax4.set_xlim(4850, 4880)
        ax4.set_ylim(1.25, 2.0)

        ax4.set_xlabel(r'Observed Wavelength [$\AA$]', fontsize = 24)
        ax4.set_ylabel(r'Flux', fontsize = 24)
        ax4.xaxis.set_major_locator(MultipleLocator(5))
        ax4.xaxis.set_minor_locator(MultipleLocator(1))
        ax4.tick_params(which = 'major', direction = 'in', top = True, right = True, length = 10, width = 2, labelsize = 16)
        ax4.tick_params(which = 'minor', direction = 'in', top = True, right = True, length = 7.5, width = 1)
        ax4.patch.set_edgecolor('black')  
        ax4.patch.set_linewidth(2)
        
        ax4.set_yticks([])
        ax4.get_xaxis().get_major_formatter().set_useOffset(False)
        gifname = 'Mp=' + str(round(Mp, 2)) + 'Me' + '_a=' + str(round(a, 2)) + 'AU_e=' + str(round(e, 2)) + '_i=' + str(round(inc, 2)) + 'rad' + '_w=' + str(round(w,2))+'rad'

        plt.savefig(directory + gifname + '.png', dpi = 300)
        
    else:
        ani = matplotlib.animation.FuncAnimation(fig, animate, frames=len(x))


        from IPython.display import HTML
        HTML(ani.to_jshtml())


        writer = matplotlib.animation.PillowWriter(fps=15,
                                         metadata=dict(artist='Me'))
        gifname = 'Mp=' + str(round(Mp, 2)) + 'Me' + '_a=' + str(round(a, 2)) + 'AU_e=' + str(round(e, 2)) + '_i=' + str(round(inc, 2)) + 'rad' + '_w=' + str(round(w,2))+'rad'
        ani.save(directory + gifname + '.gif', writer=writer)
