# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 10:08:40 2023

@author: QiRao
"""
import numpy as np

def K6(T,w,K2):
    if w < 5:
        KH = np.exp(31.96 - 0.0693*T + (3.27e-4*T - 0.4193)*w)
    else:
        KH = np.exp(30.086 - 0.0693*T - (0.197 - 3.27e-4*T)*w + 1.227*np.log10((100 - w)/w**2))
    K6 = KH/K2**1.5/101.33**2        # [kPa^-0.5]
    return K6

def K6_Joshi(T,w,K2):
    K6 = np.exp(7.412 - 20.28921*w/100 + 32.47322*(w/100)**2 - 30.87*(w/100)**3)/101.33**0.5  # [kPa^-0.5]
    return K6

def K6_Suchak(T,w,K2):
    lnKH_Suchak = 2.188e7/T**2.58 -  4.571e4/T**1.424 * w/100
    KH_Suchak = np.exp(lnKH_Suchak)  # atm^-2 
    K6_Suchak = KH_Suchak/K2**1.5/101.33**2    # kPa^(-1/2)
    return K6_Suchak