    
def Liq_func(R,Tc,Pc,W,T,P,xi,phi_sat):
    
    import math 
    import numpy as np

    x_EtOH = xi[2]/(xi[2]+xi[3]) 
    x_H2O = xi[3]/(xi[2]+xi[3])

    c1 = [1.336e2 , 3.541e1 , 8.649e1 , 6.593e1] ; c2 = [-4.735e3 , -9.662e2 , -7.931e3 , -7.228e3]
    c3 = [0,0,0,0]  ; c4 = [-2.127e1 , -4.318 , -1.025e1 , -7.177]
    c5 = [4.091e-2 , 7.932e-5 , 6.389e-6 , 4.031e-6] ; c6 = [1 , 2 , 2 , 2]
    Pi_sat = []
    for i in range(0,4):
        ln_Pi_sat = c1[i] + c2[i]/(T+c3[i]) + c4[i]*math.log(T) + c5[i]*(T)**c6[i]
        Pi_sat.append(math.exp(ln_Pi_sat))
        Pi_sat[i] = Pi_sat[i]*1e3
    # print(Pi_sat)

    # find Henry parameters for CO2 and N2
    B_CO2 = [-3232.22827148438,676.278015136719]
    A_CO2 = [69.6759109497070,-183.691406250000]
    D_CO2 = [1.06596446130425e-003,-9.77071821689606e-002]
    C_CO2 = [-8.40765953063965,39.0684318542480]

    B_N2 =[1280.52441406250,7.05463600158691]
    A_N2 = [20.0323772430420,30.1560802459717]
    D_N2 = [-3.66481882520020e-003,-3.69720801245421e-004]
    C_N2 = [-1.15200400352478,-2.33697700500488] 

    H_CO2_list = [] 
    H_N2_list = []
    for i in range(0,2):
        H_CO2_list.append(math.exp(A_CO2[i]+ B_CO2[i]/T + C_CO2[i]*math.log(T) + D_CO2[i]*T))             # H [kPa]
        H_N2_list.append(math.exp(A_N2[i]+ B_N2[i]/T + C_N2[i]*math.log(T) + D_N2[i]*T))
    # print(H_CO2_list,H_N2_list)

    # calc H_CO2_overall and H_N2_overall
    H_CO2 = H_CO2_list[0]**x_EtOH * H_CO2_list[1]**x_H2O
    H_N2 = H_N2_list[0]**x_EtOH * H_N2_list[1]**x_H2O
    # print(H_CO2,H_N2)
    H = [H_CO2,H_N2]
    # print(H)


    ## find partial molar volume for all components

    # find a and b for each component as a function of their critical properties based on SRK EOS
    nu =[]
    m = [] ; alpha = [] ; ac = [] ; a = [] ; b = [] ; Tr = [] ; A = [] ; B = [] ; a_ij = [] ; a_mix = 0 ; b_mix = 0
    for i in range(0,4):
        m.append(0.48 + 1.574*W[i] - 0.176*(W[i])**2)
        Tr.append(T/Tc[i])
        alpha.append((1+m[i]*(1-math.sqrt(Tr[i])))**2)
        ac.append(0.42748*(R*Tc[i])**2/Pc[i])
        a.append(ac[i]*alpha[i])
        b.append(0.08664*R*Tc[i]/Pc[i])
    K_ij = [[0,0,0,0],[0,0,0,-0],[0,0,0,0],[0,0,0,0]] # binary interaction coefficients 

    for i in range (0,4):
        b_mix += xi[i]*b[i]
        for j in range(0,4):
            a_ij.append((math.sqrt(a[i]*a[j]))*(1-K_ij[i][j]))
            a_mix += xi[i]*xi[j]*a_ij[i+j+3*i]
            #print(a_mix)
    # print(a_ij)
    # print(b_mix,a_mix)
    

    A = a_mix*P/(R*T)**2
    B = b_mix*P/(R*T)
    # print(a,b,A,B)

    # find Real roots of an cubic polynomial
    z1 = 0.4                                         # initial guess for Z-factor
    while True:
        fz = z1**3-z1**2+(A-B-B**2)*z1-A*B
        fzD = 3*z1**2-2*z1+(A-B-B**2)
        z1 = z1-fz/fzD
        if abs(fz/fzD) < 1e-6:
            break
    # print(f"the root is {z1}")
    Sum = 1-z1
    Pro = A*B/z1
    if Sum**2-4*Pro < 0:
        Z = z1
        # print(Z)
    else:
        k1 = math.sqrt(Sum**2-4*Pro)
        k2 = Sum
        z2 = (k1 + k2)/2
        z3 = Sum - z2
        Z = min(z1,z2,z3)
        # print(z1,z2,z3,Z)
    nu = Z*R*T/P
    # print(nu)
    #print(Z)                                                       
    # print(Z,R,T,P,nu)

    dP_dn1 = [0,0,0,0] ; dP_dV = [0,0,0,0] ; dV_dn1 = [0,0,0,0] ; Summ = 0
    for j in range(0,4):
        for i in range(0,4):
            Summ += 2*xi[i]*a_ij[i+j+3*j]
            if i == 3:
                dP_dn1[j] = R*T/(nu-b_mix) + b[j]*R*T/(nu-b_mix)**2 - Summ/(nu*(nu+b_mix)) + a_mix*b[j]/(nu*(nu+b_mix)**2)
                dP_dV[j] = -R*T/(nu-b_mix)**2 + a_mix*(2*nu+b_mix)/(nu*(nu+b_mix))**2
                dV_dn1[j] = -dP_dn1[j]/dP_dV[j]
                Summ = 0
    # print(dV_dn1)    
    nu_partial_SRK = dV_dn1                                                


    nu_hysys = [] 
    rho_molar = [24.8241657423117,24.8241849078693,16.7593951554272,54.9467049787379]     # rho [kmol/m^3]
    for i in range(0,4):
        nu_hysys.append(1/rho_molar[i]*1e-3)
    # print(nu_hysys)

    # find pressure correction factor (hysys)
    PF1 = []
    Pi_sat_solv = x_EtOH*Pi_sat[2] + x_H2O*Pi_sat[3]
    # print(Pi_sat_solv)
    for i in range(0,2):                                             # CO2 and N2
        PF1.append(math.exp(nu_hysys[i]*(P-Pi_sat_solv)/(R*T)))
    for i in range(2,4):                                             # EtOH and H2O
        PF1.append(math.exp(nu_hysys[i]*(P-Pi_sat[i])/(R*T)))
    # print(PF1)                                                                                   # PF1 : using nu_hysys

    # find pressure correction factor (SRK)                                                          # PF2 : using nu_partial_SRK
    PF2 = []
    Pi_sat_solv = x_EtOH*Pi_sat[2] + x_H2O*Pi_sat[3]
    # print(Pi_sat_solv)
    for i in range(0,2):                                              # CO2 and N2
        PF2.append(math.exp(nu_partial_SRK[i]*(P-Pi_sat_solv)/(R*T)))
    for i in range(2,4):                                              # EtOH and H2O
        PF2.append(math.exp(nu_hysys[i]*(P-Pi_sat[i])/(R*T)))
    # print(PF2)


    # calc activity coefficients for EtOH and H2O

    A = [1332.31201171875,-109.633903503418]              # 1 : EtOH-H2O 2: H2O-EtOH
    B = [1.00000000317108e-029,1.00000000317108e-029]
    C = [1.00000000317108e-029,1.00000000317108e-029] 
    alpha_ij1 = 0.303099006414413
    alpha_ij2 = 1.00000000317108e-029
    tu_ij = [0,0] ; t = 50  ; G_ij =[0,0] ; R_new = 1.98721  # t [°C] ; R_new [cal/(gmol K)]           
    for i in range(0,2):
        tu_ij[i] = (A[i] + B[i]*t + C[i]/T)/(R_new*T)
    a_ij = alpha_ij1 + alpha_ij2
    for i in range(0,2):
        G_ij[i] = math.exp(-tu_ij[i]*a_ij)
    Gamma = [0,0]       # 1: EtOH    2 : H2O
    x_prime = [x_EtOH,x_H2O]
    Gamma_EtOH = math.exp((x_prime[1])**2*(tu_ij[1]*(G_ij[1]/(x_prime[0]+x_prime[1]*G_ij[1]))**2 + tu_ij[0]*G_ij[0]/(x_prime[1]+x_prime[0]*G_ij[0])**2))
    Gamma_H2O = math.exp((x_prime[0])**2*(tu_ij[0]*(G_ij[0]/(x_prime[1]+x_prime[0]*G_ij[0]))**2 + tu_ij[1]*G_ij[1]/(x_prime[0]+x_prime[1]*G_ij[1])**2))
    # print(Gamma_EtOH,Gamma_H2O)
    Gamma = [Gamma_EtOH,Gamma_H2O]


    fi_L = [0,0,0,0]
    for i in range(0,4):
        if i == 0 or i == 1:
            fi_L[i] = xi[i]*H[i]*1000*PF2[i]
        else:
            fi_L[i] = xi[i]*Gamma[i-2]*phi_sat[i-2]*Pi_sat[i]*PF2[i]
            
    
    #print(fi_L)
    #print(Gamma)
    #print(phi_sat)
    #print(PF2)
    #print(nu_partial_SRK)
    #print(nu_hysys) 
    #print(a_mix,b_mix)
    # print(H)
    return fi_L


        
