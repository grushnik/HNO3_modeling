# -*- coding: utf-8 -*-
"""
Created on Fri May 26 12:58:15 2023

@author: QiRao
"""

import numpy as np

def column(vars):
    '''
    Packing column parameters
    
    Parameters:
        T - absolute temperature [K]
        pT - total pressure [kPa]
    Returns:
        area - effective interfacial area [m^2/m^3]
        ap - specific dry surface area [m^2/m^3]
        Re_L - Reynolds number
        Pm - power consumption per unit mass of the gas [W/kg]
        kL - liquid-side mass transfer coefficient [m/s]
        kG_star - gas-side mass transfer coefficient without multiplying with yT [(kmol/s)/(m^2*(kN/m^2))]
    '''
    T, pT, vG, vL, rho_G, mu_G, D_G, rho_L, mu_L, D_L, mu_w, para_packing, dp, epsilon = vars
    D_NO, D_NO2, D_N2O3, D_N2O4, D_HNO2, D_HNO3, D_H2O = D_G
    R = 8.314                      
    f, m, n, l, alpha, beta, Z, s, eta = para_packing
    gamma = n
    
    area = alpha*vL**beta*dp**gamma/epsilon**3
    ap = m*dp**n
    Re_L = rho_L*vL*dp/mu_L/epsilon
    epsilon_L = (1.53e-4 + 2.9e-5*Re_L**0.66*(mu_L/mu_w)**0.75)*dp**(-1.2)
    Pm = f*vG**3*ap/6/(epsilon - epsilon_L)**4
    kL = D_L/l*eta*Re_L*np.sqrt(mu_L/rho_L/D_L)
    
    kG_NO_star = D_NO/(R*T*l)*0.553*((Pm*l)**(1/3)*l*rho_G/mu_G)**0.62*(mu_G/rho_G/D_NO)**(1/3)
    kG_NO2_star = D_NO2/(R*T*l)*0.553*((Pm*l)**(1/3)*l*rho_G/mu_G)**0.62*(mu_G/rho_G/D_NO2)**(1/3)
    kG_N2O3_star = D_N2O3/(R*T*l)*0.553*((Pm*l)**(1/3)*l*rho_G/mu_G)**0.62*(mu_G/rho_G/D_N2O3)**(1/3)
    kG_N2O4_star = D_N2O4/(R*T*l)*0.553*((Pm*l)**(1/3)*l*rho_G/mu_G)**0.62*(mu_G/rho_G/D_N2O4)**(1/3)
    kG_HNO2_star = D_HNO2/(R*T*l)*0.553*((Pm*l)**(1/3)*l*rho_G/mu_G)**0.62*(mu_G/rho_G/D_HNO2)**(1/3)
    kG_HNO3_star = D_HNO3/(R*T*l)*0.553*((Pm*l)**(1/3)*l*rho_G/mu_G)**0.62*(mu_G/rho_G/D_HNO3)**(1/3)
    kG_H2O_star = D_H2O/(R*T*l)*0.553*((Pm*l)**(1/3)*l*rho_G/mu_G)**0.62*(mu_G/rho_G/D_H2O)**(1/3)
    kG_star = [kG_NO_star, kG_NO2_star, kG_N2O3_star, kG_N2O4_star, kG_HNO2_star, kG_HNO3_star, kG_H2O_star]
    
    return area, ap, Re_L, epsilon_L, Pm, kL, kG_star