# -*- coding: utf-8 -*-
"""
Created on Tue May 30 11:50:16 2023

@author: QiRao
"""

import numpy as np
import matplotlib.pyplot as plt

def psatPure(T):
    '''
    Equilibrium vapor pressures for pure nitric acid and water
    
    Parameters:
        T - absolute temperature [K]
    Returns:
        p0_H2O, p0_HNO3 - vapor pressures [Pa]
    '''
    Torr = 101325/760
    mbar = 1e-3*100000
    p0_HNO3 = 10**(7.61628 - 1486.238/(T-43)) * Torr
    p0_H2O = 10**(8.42926609 - 1827.17843/T - 71208.271/T**2) * mbar
    return p0_H2O, p0_HNO3  

def psatMix(x_H2O,x_HNO3,T):
    '''
    Calculates partial vapor pressures over 
    a solution of HNO3 in water
    
    Parameters:
        x_H2O, x_HNO3 - molar fractions in the solution
        T - absolute temperature [K]
    Returns:
        p_H2O, p_HNO3 - partial pressures [Pa]
    '''
    p0_H2O, p0_HNO3 = psatPure(T)
    
    A1 = -391.43 - 7.44e4/T
    A2 = -627.739- 1.406e5/T
    B1 = 0.5695
    B2 = 1/B1
    
    x1 = x_H2O
    x2 = x_HNO3
    
    c1 = A1*x2**2/(x2 + B1*x1)**2
    c2 = A2*x1**2/(x1 + B2*x2)**2
    
    gamma1 = 10**(c1/T)
    gamma2 = 10**(c2/T)
    
    p_H2O = gamma1 * x1 * p0_H2O
    p_HNO3 = gamma2 * x2 * p0_HNO3
    
    return p_H2O, p_HNO3