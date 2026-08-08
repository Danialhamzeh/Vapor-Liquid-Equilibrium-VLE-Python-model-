# components are put in the order of CO2, N2, C2H5OH, H2O
import math
from Rachford_Rice import Rachford
from Vapor import Vap_func
from Liquid import Liq_func

T = 50+273.15 ; P = 50e5 ; R = 8.314         # T [K] , P [Pa] , and R [J/(mol.K)]
zi = [0.25,0.05,0.4,0.3]

K_hysys = [36.3268191442449,2481.83778564272,1.79072824666026e-002,7.28087519013496e-003]

W = [0.238940000534058,3.99998016655445e-002,0.644370019435883,0.344000011682510]
Pc = [7370.00000000000e3,3394.37011718750e3,6147.00000000000e3,22120.0000000000e3]                                       # Pc [Pa]
Tc = [30.9500061035156+273.15,-146.955999755859+273.15,240.750024414063+273.15,374.149011230469+273.15]                  # Tc [K] 


# K_val = [20,2300,0.1,5e-3]                                           # initial guess for K_values
K_val = [2,2,0.01,0.01]
while True:
    suum = 0 
    Beta, xi, yi = Rachford (zi, K_val)
    fi_v, phi_sat, Hr_mix = Vap_func (R,Tc,Pc,W,T,P,yi)
    fi_L = Liq_func (R,Tc,Pc,W,T,P,xi,phi_sat)
    for i in range(0,4):
        suum += (1-fi_L[i]/fi_v[i])**2 
    if suum < 1e-12:
        break
    else:
        for i in range(0,4):
            K_val[i] = fi_L[i]/fi_v[i] * K_val[i]
    
# print(Beta,xi,yi)    
# print(fi_v,phi_sat)
# print(fi_L)
# print(K_val)
# print(K_hysys)
# print(Hr_mix)

print(f"Vapor Fraction (Beta) : {Beta}")
print(f"Liquid Mole Fractions (x) : CO2 = {xi[0]}, N2 = {xi[1]},  EtOH = {xi[2]},   H2O = {xi[3]}")
print(f"Vapor Mole Fractions (y) : CO2 = {yi[0]}, N2 = {yi[1]},  EtOH = {yi[2]},   H2O = {yi[3]}")
print(f"K-values (Python) : CO2 = {K_val[0]},  N2 = {K_val[1]},  EtOH = {K_val[2]},  H2O = {K_val[3]}")
print(f"K-values (HYSYS) : CO2 = {K_hysys[0]},  N2 = {K_hysys[1]},  EtOH = {K_hysys[2]},  H2O = {K_hysys[3]}")
print(f"Liquid Fugacity (fi_L) : CO2 = {fi_L[0]},  N2 = {fi_L[1]},  EtOH = {fi_L[2]},  H2O = {fi_L[3]}")
print(f"Vapor Fugacity (fi_L) : CO2 = {fi_v[0]},  N2 = {fi_v[1]},  EtOH = {fi_v[2]},  H2O = {fi_v[3]}")






