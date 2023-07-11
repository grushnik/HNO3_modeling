# -*- coding: utf-8 -*-
"""
Created on Mon May 22 14:51:04 2023

@author: QiRao
"""

import numpy as np
from scipy.optimize import fsolve

# Equations (1.7)-(1.11)

def eqConst(T):
    """
    Calculates NOx equilibrium constants in gas phase
    
    Parameters:
    T - absolute temperature [K]
    
    Returns:
    k1 - rate constant for 2NO + O2 <-> NO2 reaction [1/(kPa^2*s)]
    K2 - equilibrium constant for 2NO2 <-> N2O4 reaction [1/kPa]
    K3 - equilibrium constant for NO + NO2 <-> N2O3 reaction [1/kPa]
    K4 - equilibrium constant for NO + NO2 + H2O <-> 2HNO2 reaction [1/kPa]
    K5 - equilibrium constant for 3NO2 + H2O <-> 2HNO3 + NO reaction [1/kPa]
    """
    logk1 = 652.1/T - 4.747
    logK2 = 2993/T - 11.232
    logK3 = 2072/T - 9.2397
    logK4 = 2051.17/T - 8.7385
    logK5 = 2003.8/T - 10.763
    
    k1 = 10**logk1
    K2 = 10**logK2
    K3 = 10**logK3
    K4 = 10**logK4
    K5 = 10**logK5

    return k1,K2,K3,K4,K5
    
def gas(vars):
    T, pT, yO2, yN_star, yNO_star, yH2O_star, K = vars
    
    # Dependent variables from equations (3.6)-(3.9)

    def depVars(vars,*K):
        yT, yNO, yNO2, yH2O= vars
        K2, K3, K4, K5, pT = K
        yN2O4 = K2*pT*yNO2**2/yT
        yN2O3 = K3*pT*yNO*yNO2/yT
        yHNO2 = np.sqrt(K4*pT*yNO*yNO2*yH2O/yT)
        yHNO3 = np.sqrt(K5*pT*yNO2**3*yH2O/(yNO*yT))
#        if K4*pT*yNO*yNO2*yH2O/yT < 0:
#            print(yNO,yNO2,yH2O)
        return yN2O4, yN2O3, yHNO2, yHNO3
    
    # System of nonlinear equations (2.2), (2.4), (2.5), (2.6)
    
    def eqns(vars,*K):
        yT, yNO, yNO2, yH2O = vars
        K2, K3, K4, K5, pT = K
        yN2O4, yN2O3, yHNO2, yHNO3 = depVars(vars,*K)
        eq1 = -yT + 1 + yO2 + yH2O + yNO + yNO2 + yN2O3 + yN2O4 + yHNO2 + yHNO3
        eq2 = -yN_star + yNO + yNO2 + 2*yN2O3 + 2*yN2O4 + yHNO2 + yHNO3
        eq3 = -yNO_star + yNO + yN2O3 + 0.5*yHNO2 - 0.5*yHNO3 
        eq4 = -yH2O_star + yH2O + 0.5*yHNO2 + 0.5*yHNO3
        return eq1,eq2,eq3,eq4
    
    # Intial guess

    yT_0 = 1+yO2+yH2O_star+yN_star
    yNO_0 = yNO_star
    yNO2_0 = yN_star-yNO_star
    yH2O_0 = yH2O_star
    
    vars_0 = np.array([yT_0, yNO_0, yNO2_0, yH2O_0])
    
    # Solve

    sol = fsolve(eqns,vars_0,args=K)
    
    #    print('Max residual: %e' %(np.max(np.abs(eqns(sol)))))

    # Extract independent and dependent variables

    yT, yNO, yNO2, yH2O = sol     # independent
    yN2O4, yN2O3, yHNO2, yHNO3 = depVars(sol,*K)  # dependent
    results = [yT, yNO, yNO2, yH2O, yN2O4, yN2O3, yHNO2, yHNO3]
    return results


if __name__ == "__main__":
    # Input data

    Tc = 31 # [C]
    pT = 101.33  # [kPa]
    yO2 = 0.267
    yN_star = 0.541
    yNO_star = 0.72*yN_star
    yH2O_star = 0.350
    
    T = 273.16 + Tc    # [K]

    k1, K2,K3,K4,K5 = eqConst(T)
    K = (K2, K3, K4, K5, pT)
    print('K2 = %f atm^-1\nK3 = %f atm^-1\nK4 = %f atm^-1\nK5 = %f atm^-1' %(K2,K3,K4,K5))
    vars = [T, pT, yO2, yN_star, yNO_star, yH2O_star, K]
    yT, yNO, yNO2, yH2O, yN2O4, yN2O3, yHNO2, yHNO3 = gas(vars)

    # Report

    print('yN2   = %f     xN2   = %5.2f %%\nyO2   = %f     xO2   = %5.2f %%\n'\
          'yH2O  = %f     xH2O  = %5.2f %%\nyNO   = %f     xNO   = %5.2f %%\n'\
          'yNO2  = %f     xNO2  = %5.2f %%\n'\
          'yN2O3 = %f     xN2O3 = %5.2f %%\nyN2O4 = %f     xN2O4 = %5.2f %%\n'\
          'yHNO2 = %f     xHNO2 = %5.2f %%\nyHNO3 = %f     xHNO3 = %5.2f %%' \
           %(1, 1/yT*100, yO2, yO2/yT*100, yH2O, yH2O/yT*100, yNO, yNO/yT*100, \
             yNO2, yNO2/yT*100, yN2O3, yN2O3/yT*100, yN2O4, yN2O4/yT*100, \
             yHNO2, yHNO2/yT*100, yHNO3, yHNO3/yT*100))
    print('----------------     ----------------')
    print('yT    = %f     xT    = %.2f %%' %(yT, 100))