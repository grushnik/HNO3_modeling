# -*- coding: utf-8 -*-
"""
Created on Fri May 26 13:35:25 2023

@author: QiRao
"""
import numpy as np

def HkD(vars):
    T, pT, W, rho = vars
    W = W/100
    Ma = 63
    
    factor = 10**(1.872 - 1.424*(W/(1-W)**1.83))*np.exp(-1285.3/T)
    HkD_NO2 = 2.2e-7*factor
    HkD_N2O4 = 8.18e-6*factor
    HkD_N2O3 = 1.57e-5*factor
    
    Ca = rho*W/Ma    # [kmol/m^3]  HNO3 concentration
    Hw_HNO2 = 0.484  # [kmol/m3/kN/m2]
    ks = -0.111 + 0.323 - 0.000517
    H_HNO2 = Hw_HNO2*10**(-1*ks*Ca)
    
    return HkD_NO2, HkD_N2O4, HkD_N2O3, H_HNO2


    
    
    