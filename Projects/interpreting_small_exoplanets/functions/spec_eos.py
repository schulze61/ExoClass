# -*- coding: utf-8 -*-
"""
Created on Sun Jan  7 11:57:27 2024

@author: jgsch
"""

import pandas as pd
import numpy as np
import scipy.interpolate as si
import scipy.optimize as so

from general_eos import *



mfe = 55.845

#------------------------------iron---------------------------------------#
hfe = pd.read_csv('./data_files/hakim_data.csv')
def hakim_dft(P):
    l = 7
    ind = np.argmin(abs(P-hfe['P']))
    parr = np.array(hfe.loc[ind-l:ind+l, 'P'])
    varr = np.array(hfe.loc[ind-l:ind+l, 'Volume'])
    rho = 1000*mfe/(varr)
    p = so.curve_fit(BM3, rho, parr, p0 = [9000.0, 170.0, 5.0])
    return p[0]

def Anderson_lFe(rho, P):
    rho0 = 7019.0
    k = 109.7
    kp = 4.66
    kpp = -0.043
        
    return P - BM4(rho, rho0, k, kp, kpp)


def Anderson_sFe(rho, P):
    K0 = 156.2
    Kp0 = 6.08
    rho0 = 8.30*1000.0
    return P - vinet(rho, rho0, K0, Kp0)


def Z16_eFe(rho, P):
    rho0 = 7050.0
    K0 = 201
    return P - BM2(rho, rho0, K0)


def Sotin_core(rho, P):
    #rho0 = (0.87 + 0.13)/(0.87/8340 + 0.13/4900)
    rho0 = 8340.0
    K0 = 135.0
    Kp0 = 6.0
    return P - BM3(rho, rho0, K0, Kp0)


def Belonoshko_Fe(rho, P):
    rho0 = 1000*55.845/6.695
    K0 = 173.98
    Kp0 = 5.297
    return P - BM3(rho, rho0, K0, Kp0)


def Smith_sFe(rho, P):
    rho0 = 1000.0*8.43
    K0 = 177.7
    Kp0 = 5.64
    return P - vinet(rho, rho0, K0, Kp0)

def hakim_sFe(rho, P):
    if P < 234.4:
        rho0, K0, Kp0 = hakim_dft(P)
        return P - BM3(rho, rho0, K0, Kp0)
    else:
        rho0 = 1000*mfe/4.28575
        P0 = 234.4
        K0 = 1145.7
        c0 = 3.19
        c2 = -2.40
        return P - holzapfel_hakim(rho, rho0, K0, c0, c2, P0)
    
    
def bouchet_hcpFe(rho, P):
    rho0 = 1000.0*55.845/6.290
    K0 = 253.844
    Kp0 = 4.719
    P0 = 1003.6*((26.0/6.290)**(5.0/3.0))
    c0 = -np.log((3.0*K0)/P0)
    c2 = ((3.0/2.0)*(Kp0-3.0)) - c0
    return P - holzapfel(rho, rho0, K0, c0, c2)

def valencia_Fe(rho, P):
    rho0 = 8300.0
    K0 = 160.2
    Kp0 = 5.82
    return P - vinet(rho, rho0, K0, Kp0)
    
    
#----------------------------MgSiO3------------------------------------#

def Karki_pv(rho, P):
    K0 = 247.0
    Kp0 = 3.97
    Kpp0 = -0.016
    rho0 = 4.10*1000.0
    return P - BM4(rho, rho0, K0, Kp0, Kpp0)


dz16 = pd.read_csv('./data_files/SD_data.csv')
dz16 = dz16.sort_values(by = 'rho')

def Z16_mant(rho, P):
    if P > 23.83:
        rho0 = 3980.0
        K0 = 206.0
        return P - BM2(rho, rho0, K0)
    else:
        t = si.interp1d(dz16['rho'], dz16['P'], fill_value = 'extrapolate' )
        #print(t(rho), float(P), float(P)-t(rho))
        #print(float(P))
        return round(P,8) - round(float(t(rho)),8)

def Sotin_mant(rho, P):
    if P < 23.0:
        rho0 = 3215.0
        K0 = 111.0
        Kp0 = 7.0
    else:
        rho0 = 4108.0
        K0 = 263.0
        Kp0 = 3.9
    #print(BM3(rho, rho0, K0, Kp0))
    return P - BM3(rho, rho0, K0, Kp0)

def ZS13_mant(rho, P):
    if P > 122.0:
        V0 = 40.80 #A3/ 3 O atoms
        rho0 = (100.39/1000.0/(6.0221408e+23))/(V0*1e-30) #CHECK
        #print(rho0)
        K0 = 203.0
        Kp0 = 4.19
        
    else:
        V0 = 40.78 #A3/ 3 O atoms
        rho0 = (100.39/1000.0/(6.0221408e+23))/(V0*1e-30) #CHECK
        #print(rho0)
        K0 = 232.0
        Kp0 = 3.86        
    return P - BM3(rho, rho0, K0, Kp0)

def magrathea_mant(rho, P):
    if P < 115.4079944467753:
        #brg Oganov & Ono 04 V0 = 25.206
        rho0 = 1000*100.3875/25.206
        K0 = 230.05
        Kp0 = 4.142
        return P - vinet(rho, rho0, K0, Kp0)
    else:
        #Sakai 2016 -- Keane v0 = 24.73
        rho0 = 1000*100.39/24.73
        K0 = 203.0
        Kp0 = 5.35
        ginf = 0.93
        #print(keane(rho, rho0, K0, Kp0, ginf))
        return P  - (keane(rho, rho0, K0, Kp0, ginf))

LM = pd.read_csv('./data_files/0.07CaMg_0.00FeMg_0.09AlMg_0.9SiMg_0.0NaMg_0.00Fe_LM_results.txt')
UM = pd.read_csv('./data_files/0.07CaMg_0.00FeMg_0.09AlMg_0.9SiMg_0.0NaMg_0.00Fe_UM_results.txt')


LM.columns = ['P[bar]', ' T[K]', ' rho[g/cm3]', ' log10(alpha)[1/K]',
       ' cp[J/(kg*K)]', ' c2/c', ' fc2/c', ' per', ' wus', ' aperov', ' perov',
       ' fperov', ' ab', ' an', ' sp', ' herc', ' fo', ' fa', ' wad', ' fwad',
       ' ring', ' fring', ' odi', ' en', ' fs', ' ts', ' jd', ' di', ' hed',
       ' cen', ' cts', ' cor', ' aki', ' faki', ' gr', ' alm', ' maj', ' py',
       ' jmaj', ' cmaj', ' fmaj', ' appv', ' ppv', ' fppv', ' mfer', ' ffer',
       ' nfer', ' ca-pv', ' cfs', ' ky', ' neph', ' coe', ' seif', ' q', ' s',
       ' bad']

UM.columns = ['P[bar]', ' T[K]', ' rho[g/cm3]', ' log10(alpha)[1/K]',
       ' cp[J/(kg*K)]', ' c2/c', ' fc2/c', ' per', ' wus', ' aperov', ' perov',
       ' fperov', ' ab', ' an', ' sp', ' herc', ' fo', ' fa', ' wad', ' fwad',
       ' ring', ' fring', ' odi', ' en', ' fs', ' ts', ' jd', ' di', ' hed',
       ' cen', ' cts', ' cor', ' aki', ' faki', ' gr', ' alm', ' maj', ' py',
       ' jmaj', ' cmaj', ' fmaj', ' appv', ' ppv', ' fppv', ' mfer', ' ffer',
       ' nfer', ' ca-pv', ' cfs', ' ky', ' neph', ' coe', ' seif', ' q', ' s',
       ' bad']




UM['P_GPa'] = 0.0001*UM['P[bar]']
LM['P_GPa'] = 0.0001*LM['P[bar]']

UM['rho_kgm3'] = 1000*UM[' rho[g/cm3]']
LM['rho_kgm3'] = 1000*LM[' rho[g/cm3]'] 

UM = UM.iloc[np.where(UM[' T[K]'] == min(UM[' T[K]']))]
LM = LM.iloc[np.where(LM[' T[K]'] == min(LM[' T[K]']))]

UM = UM.sort_values(by = 'P_GPa')
LM = LM.sort_values(by = 'P_GPa')

mantle_grid = pd.concat((UM, LM))

EP_mantle_EOS = si.interp1d(mantle_grid['P_GPa'], mantle_grid['rho_kgm3'], fill_value = 'extrapolate')

def EP_mantle(P):
    rho = np.array([EP_mantle_EOS(P) for P in P])
    return rho


def valencia_mant(rho, P):
    if P <= 11.4675:
        rho0 = 3347.0
        K0 = 126.8
        Kp0 = 4.274
    elif (P > 11.4675) and (P <= 26.546):
        rho0 = 3644.0
        K0 = 174.5
        Kp0 = 4.247
    elif (P > 26.546) and (P <= 106.578947368):
        rho0 = 4152.0
        K0 = 223.6
        Kp0 = 4.274
    elif P > 106.578947368:
        rho0 = 4270.0
        K0 = 233.6
        Kp0 = 4.524
        
    return P - vinet(rho, rho0, K0, Kp0)


#---------------------------- wat/ice ----------------------------------------#
seag_dft = np.loadtxt('./data_files/seag_DFT.txt')
interp_seag_dft = si.interp1d(seag_dft[:, 1]*1000.0, seag_dft[:, 2], fill_value = 'extrapolate')

sf_wat = pd.read_csv('./data_files/sf_wat.csv')

def seag_h2o(rho, P):
    if P <=44.3:
        return P - BM3(rho,  1460.0, 23.7, 4.15)
    else:
        return P - interp_seag_dft(rho)
    
    
def magrathea_h2o(rho, P):
    mh2o = 18.01528
    if P <=  0.993:
        rho0 = (mh2o/18.047)*1000.0
        K0 = 2.18
        Kp0 = 4.0
    if (P > 0.993) and (P <= 3.0):
        rho0 = (mh2o/14.17)*1000.0
        K0 = 14.05
        Kp0 = 4.0
    if P > 3.0:
        rho0 = (mh2o/12.49)*1000.0
        K0 = 20.15
        Kp0 = 4.0
    
    return P - BM3(rho, rho0, K0, Kp0) 

def EP_h2o(rho, P):
    mh2o = 18.01528

    
    if P <= 0.993:
        interp_sf = si.interp1d(sf_wat['rho_kgm3'], sf_wat['p_GPa'], fill_value = 'extrapolate')
        return P - float(interp_sf(rho))
    if (P > 0.993) and (P <= 3.0):
        #out = sf.seafreeze(PT, 'VI')
        #return rho-out.rho[0]
        rho0 = (mh2o/13.62)*1000.0
        K0 = 15.2
        Kp0 = 6.5
        return P - BM3(rho, rho0, K0, Kp0)   
    if P > 3.0:
        rho0 = (mh2o/12.49)*1000.0
        K0 = 20.15
        Kp0 = 4.0
        return P - BM3(rho, rho0, K0, Kp0)  
    
    
def sotin_h2o(rho, P):   
    if P <= 0.993:
        rho0 = 1000.0
        K0 = 2.2
        Kp0 = 4.0
    else:
        rho0 = 1460.0
        K0 = 23.9
        Kp0 = 4.2
    return P - BM3(rho, rho0, K0, Kp0)  


def valencia_h2o(rho, P):
    mh2o = 18.01528
    if P <=  0.993:
        rho0 = 998.23
        K0 = 2.18
        Kp0 = 4.0
        return P - BM3(rho, rho0, K0, Kp0)
    else:
        rho0 = 1463.0
        K0 = 2.308*10
        Kp0 = 4.532
        return P - vinet(rho, rho0, K0, Kp0)
        
    
        





