# -*- coding: utf-8 -*-
"""
Created on Tue May 30 10:51:32 2023

@author: QiRao
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

def mu_HNO3(wPer,Tc):
    '''
    Viscosity of nitric-acid solution
    
    Parameters
        wPer - concentration [% w/w]
        Tc - temperature [C]
    Return:
        mu - dynamic viscosity [mPa * s]
        
    Note: 0 < T < 100C, 0 < wPer < 60
    '''
    
    # Pure water
    
    T = Tc + 273.16
    mu0 = np.exp(-3.7188 + 578.919/(T-137.546))
    
    # 60% HNO3
    
    Tp = np.array([0,20,50,100])
    mup = np.array([3.5,2.5,1.4,0.7])
    mu_interp = interp1d(Tp,mup,'cubic')
    mu60 = mu_interp(Tc)
    
    # wPer HNO3
    
    mu = mu0 + wPer/60 * (mu60 - mu0)
    return mu

if __name__ == "__main__":
    Tc = np.linspace(0,100,101)

    wPer = 20   # percentage concentration HNO3 [% w/w]
    
    plt.figure(figsize=(12,8))
    plt.plot(Tc,mu_HNO3(0,Tc),'b',label = 'water')
    plt.plot(Tc,mu_HNO3(60,Tc),'r', label='60% HNO3')
    plt.plot(Tc,mu_HNO3(wPer,Tc),'--k', label='%d%% HNO3' %wPer)
    
    plt.xlabel('Temperature - T[C]')
    plt.ylabel(r'Dynamic viscosity - $\mu$ [mPa*s]')
    plt.legend()
    plt.axis([0,100,0,4])
    plt.grid()
    plt.draw()