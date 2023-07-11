# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 09:45:54 2023

@author: QiRao
"""

def DGas(vars):
    '''
    Diffusivities for gas species
    
    Parameters:
        T - absolute temperature [K]
        pT - total pressure [kPa]
    Returns:
        D_NO - diffusivity of NO [m^2/s]
        D_NO2 - diffusivity of NO2 [m^2/s]
        D_N2O3 - diffusivity of N2O3 [m^2/s]
        D_N2O4 - diffusivity of N2O4 [m^2/s]
        D_HNO2 - diffusivity of HNO2 [m^2/s]
        D_HNO3 - diffusivity of HNO3 [m^2/s]
        D_H2O - diffusivity of H2O [m^2/s]
    '''
    
    T, pT = vars
    pT = pT
    
    D0_NO = 2.3e-5    #[m2/s]  at 296 K
    D0_NO2 = 1.4e-5   #[m2/s]  at 296 K
    D0_N2O3 = 0.81e-5  #[m2/s]  at 273 K
    D0_N2O4 = 0.84e-5  #[m2/s]  at 273 K
    D0_HNO2 = 1.3e-5  #[m2/s]  at 296 K
    D0_HNO3 = 1.1e-5  #[m2/s]  at 296 K
    D0_H2O = 2.49e-5   #[m2/s]  at 293 K
    
    D_NO = D0_NO*101.33/pT*(T/296)**1.75
    D_NO2 = D0_NO2*101.33/pT*(T/296)**1.75
    D_N2O3 = D0_N2O3*101.33/pT*(T/273)**1.75
    D_N2O4 = D0_N2O4*101.33/pT*(T/273)**1.75
    D_HNO2 = D0_HNO2*101.33/pT*(T/296)**1.75
    D_HNO3 = D0_HNO3*101.33/pT*(T/296)**1.75
    D_H2O = D0_H2O*101.33/pT*(T/293)**1.75
    
    return D_NO, D_NO2, D_N2O3, D_N2O4, D_HNO2, D_HNO3, D_H2O