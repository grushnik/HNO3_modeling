# -*- coding: utf-8 -*-
"""
Created on Tue May 23 15:07:01 2023

@author: QiRao
"""

import numpy as np
from scipy.optimize import fsolve

def interface(vars):
    T, pT, yO2, yN_star, yNO_star, yH2O_star, K = vars
    K2, K3, K4, K6, a, HkD_NO2, HkD_N2O4, HkD_N2O3, H_HNO2, pH2O_i, pNO_o, pNO2_o, pN2O4_o, pN2O3_o, pHNO2_o, kGa_NO2, kGa_N2O4, kGa_N2O3, kGa_HNO2, kGa_NO, kL = K
    def eqns(vars,*K):
        pNO_i, pNO2_i = vars
        K2, K3, K4, K6, a, HkD_NO2, HkD_N2O4, HkD_N2O3, H_HNO2, pH2O_i, pNO_o, pNO2_o, pN2O4_o, pN2O3_o, pHNO2_o, kGa_NO2, kGa_N2O4, kGa_N2O3, kGa_HNO2, kGa_NO, kL = K
        
        A1 = a*HkD_NO2/3
        A2 = 1/(K2**0.5*K6**(1/3))
        A3 = 2/3*K2*a*HkD_N2O4
        A4 = 1/(K2*K6**(2/3))
        A5 = 1/3*K3*a*HkD_N2O3
        A7 = 1/6*K4**0.5*pH2O_i**0.5*kL*a*H_HNO2
        A8 = 1/(K2**(1/4)*K6**(1/6))
        A9 = -1*K3*kGa_N2O3
        A10 = -0.5*K4**0.5*pH2O_i**0.5*kGa_HNO2
        A11 = -1*kGa_NO
        A12 = kGa_NO*pNO_o + kGa_N2O3*pN2O3_o + 0.5*kGa_HNO2*pHNO2_o
        
        B1 = 3*A1
        B2 = A2
        B3 = kGa_NO2
        B4 = 3*A3
        B5 = A4
        B6 = 2*K2*kGa_N2O4
        B7 = 3*A5
        B8 = -1*A9
        B9 = 3*A7
        B10 = -1*A10
        B11 = -1*kGa_NO2*pNO2_o - kGa_N2O3*pN2O3_o - 2*kGa_N2O4*pN2O4_o - 0.5*kGa_HNO2*pHNO2_o
        B12 = A8
        
        if pNO2_i - A2*pNO_i**(1/3) >= 0:
        
            eq1 = A1*(pNO2_i - A2*pNO_i**(1/3))**1.5 + A3*(pNO2_i**2 - A4*pNO_i**(2/3)) \
                  + A5*(pNO_i*pNO2_i - A2*pNO_i**(4/3)) + A7*(pNO_i**0.5*pNO2_i**0.5 - A8*pNO_i**(2/3)) \
                  + A9*pNO_i*pNO2_i + A10*pNO_i**0.5*pNO2_i**0.5 + A11*pNO_i + A12
            eq2 = B1*(pNO2_i - B2*pNO_i**(1/3))**1.5 + B3*pNO2_i + B4*(pNO2_i**2- B5*pNO_i**(2/3)) \
                  + B6*pNO2_i**2 + B7*(pNO_i*pNO2_i - B2*pNO_i**(4/3)) + B9*(pNO_i**0.5*pNO2_i**0.5 - B12*pNO_i**(2/3)) \
                  + B8*pNO_i*pNO2_i + B10*pNO_i**0.5*pNO2_i**0.5 + B11
        else:
            eq1 = A9*pNO_i*pNO2_i + A10*pNO_i**0.5*pNO2_i**0.5 + A11*pNO_i + A12
            eq2 = B3*pNO2_i + B6*pNO2_i**2 + B8*pNO_i*pNO2_i + B10*pNO_i**0.5*pNO2_i**0.5 + B11

        return eq1,eq2
    
    # Intial guess
    
    A2 = 1/(K2**0.5*K6**(1/3))
    A4 = 1/(K2*K6**(2/3))
    A8 = 1/(K2**(1/4)*K6**(1/6))
    
    g1 = (pNO2_o/A2)**3
    g2 = (pNO2_o**2/A4)**1.5
    g3 = (pNO2_o**0.5/A8)**6
    g4 = (A2*pNO_o*1.08)**(1/3)
#    pNO_i0 = min([pNO_o, g1]) #pNO_o
#    pNO2_i0 = pNO2_o*1.005*2
    
    pNO_i0 = pNO_o*1.001
    pNO2_i0 = max([pNO2_o*0.99, g4])

    vars_0 = np.array([pNO_i0, pNO2_i0])
    # Solve

    sol = fsolve(eqns,vars_0,args=K)
    
#     N = 50
#     best_obj = np.abs(eqns(sol,*K)[0]) +  np.abs(eqns(sol,*K)[1])
#     best_x = sol
# #    print(np.abs(eqns(sol,*K)[0]) +  np.abs(eqns(sol,*K)[1]))
#     x = np.linspace(pNO_o*0, pNO_o*2,N)
#     y = np.linspace(0, pNO2_o*2,N)
#     for x0 in x:
#         for y0 in y:
#             if x0<=(y0/A2)**3:
#                 sol_new = fsolve(eqns, np.array([x0, y0]),args=K)
#                 fun = np.abs(eqns(sol_new,*K)[0]) +  np.abs(eqns(sol_new,*K)[1])
#                 if fun < best_obj:
#                     best_obj = fun
#                     best_x = sol_new
    # Extract independent and dependent variables
    pNO_i, pNO2_i = sol #best_x #sol
    
    pN2O4_i = K2*pNO2_i**2
    pN2O3_i = K3*pNO_i*pNO2_i
    pHNO2_i = K4**0.5*pH2O_i**0.5*pNO_i**0.5*pNO2_i**0.5
    
    pN2O4_b = pNO_i**(2/3)/(K6**(2/3))
    pNO2_b = pNO_i**(1/3)/(K2**0.5*K6**(1/3))
    pN2O3_b = K3*pNO_i**(4/3)/(K2**0.5*K6**(1/3))
    pHNO2_b = K4**0.5*pH2O_i**0.5*pNO_i**(2/3)/(K2**0.25*K6**(1/6))
    
    results = [pNO_i, pNO2_i, pN2O4_i, pN2O3_i, pHNO2_i, pN2O4_b, pNO2_b, pN2O3_b, pHNO2_b]
    return results

if __name__ == "__main__":
    
    Tc = 25 # [C]
    pT = 1  # [atm]
    yO2 = 0.267
    yN_star = 0.541
    yNO_star = 0.72*yN_star
    yH2O_star = 0.350
    
    T = 273.16 + Tc    # [K]
    
    K = (6.494880861116165,\
 0.5194240482149102,\
 1.4023547514847117,\
 129.15948645746076,\
 124.9262565069449,\
 2.0733546121789872e-07,\
 7.709109421647324e-06,\
 1.4796212459640954e-05,\
 0.21415985298976148,\
 2.493300492301609,\
 0.5374590322676245,\
 5.895653406720646,\
 2.2279068818464474,\
 0.0162428158331849,\
 0.19737915176305432,\
 0.00021291067715863195,\
 0.0001664493803892933,\
 0.00016246232770682095,\
 0.00020264737647945273,\
 0.00029643579319598313,\
 4.41596633964505e-05)
    
    vars = [T, pT, yO2, yN_star, yNO_star, yH2O_star, K]
    pNO_i, pNO2_i, pN2O4_i, pN2O3_i, pHNO2_i, pN2O4_b, pNO2_b, pN2O3_b, pHNO2_b = interface(vars)
        
    