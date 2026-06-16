# -*- coding: utf-8 -*-
"""
Created on Sun Jan  7 11:53:41 2024

@author: jgsch
"""

import numpy as np

def BM4(rho, rho0, k, kp, kpp):
    x = (rho/rho0)**(1.0/3.0)
    ep1 = 3.0*(4.0 - kp)/4.0
    ep2 = (3.0/8.0)*k*kpp + (3.0/8.0)*kp*(kp-7.0) + 143.0/24.0
    t1 = (3.0/2.0)*k*(x**7 - x**5)
    t2 = (1.0 + ep1 - ep1*(x**2) + ep2*((x**2 - 1)**2))
    return t1*t2

def BM3(rho, rho0, K0, Kp0):
    x = (rho/rho0)**(1.0/3.0)
    t1 = (3.0/2.0)*K0*((x**7) - (x**5))
    t2 = 1.0 - (3.0/4.0)*(4.0 - Kp0)*((x**2) - 1)
    return t1*t2

def vinet(rho, rho0, K0, Kp0):
    x = rho/rho0
    return 3.0*K0*(x**(2.0/3.0))*(1.0 - x**(-1.0/3.0))*np.exp(((3.0/2.0)*(Kp0 - 1.0))*(1.0 - (x**(-1.0/3.0))))
      
def BM2(rho, rho0, K0):
    x = rho/rho0
    return (3.0/2.0)*K0*((x**(7.0/3.0)) - (x**(5.0/3.0)))

def holzapfel(rho, rho0, K0, c0, c2):
    x = rho0/rho
    t1 = 3.0*K0*(x**(-5.0/3.0))
    t2 = (1.0 - (x**(1.0/3.0)))
    return t1*t2*(1.0 + c2*(x**(1.0/3.0))*t2)*np.exp(c0*t2)


def holzapfel_hakim(rho, rho0, K0, c0, c2, P0):
    x = rho0/rho
    t1 = 3.0*K0*(x**(-5.0/3.0))
    t2 = (1.0 - (x**(1.0/3.0)))
    return P0 + t1*t2*(1.0 + c2*(x**(1.0/3.0))*t2)*np.exp(c0*t2)


def keane(rho, rho0, K0, Kp0, ginf):
    y = rho/rho0
    kpinf = 2.0*(ginf + (1.0/6.0))
    t1 = ((Kp0*K0)/(kpinf**2))*((y**(kpinf)) - 1.0)
    t2 = (Kp0 - kpinf)*(K0/kpinf)*np.log(y)
    return (t1 - t2)