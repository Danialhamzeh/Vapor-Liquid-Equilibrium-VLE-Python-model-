# components are put in the order of CO2, N2, C2H5OH, H2O

def Vap_func(R,Tc,Pc,W,T,P,yi):
    
    import math


## part one: estimate the fugacity coefficient of each component in mixture in vapour phase (Φ_i_v)


    m = [] ; alpha = [] ; ac = [] ; a = [] ; b = [] ; Tr = [] ; A = [] ; B = [] ; a_ij = [] ; a_mix = 0 ; b_mix = 0

# find a and b for each component as a function of their critical properties based on SRK EOS

    for i in range(0,4):
        m.append(0.48 + 1.574*W[i] - 0.176*(W[i])**2)
        Tr.append(T/Tc[i])
        alpha.append((1+m[i]*(1-math.sqrt(Tr[i])))**2)
        ac.append(0.42748*(R*Tc[i])**2/Pc[i])
        a.append(ac[i]*alpha[i])
        b.append(0.08664*R*Tc[i]/Pc[i])    
    

# find a_mix and b_mix based on mole fractions(y_i), a, b, and binary interaction coefficients (K_ij)

    K_ij = [[0,0,0,0],[0,0,0,-0],[0,0,0,0],[0,0,0,0]]                                              # binary interaction coefficients 
    for i in range (0,4):
        b_mix += yi[i]*b[i]
        for j in range(0,4):
            a_ij.append((math.sqrt(a[i]*a[j]))*(1-K_ij[i][j]))
            a_mix += yi[i]*yi[j]*a_ij[i+j+3*i]
            #print(a_mix)
    # print(b_mix,a_mix)
    
    

# find A and B which are function of a_mix and b_mix and are used in Cubic polynomial in terms of Z

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
        Z = max(z1,z2,z3)
        # print(z1,z2,z3,Z)
    V = Z*R*T/P
    #print(Z)                                                        ## first verifying parameter 
    # print(Z,R,T,P,V)


# last step: find Φ_i_v based on its formula for SRK EOS


    phi_v = [0,0,0,0]
    parameter_phi= [0,0,0,0]
    for i in range(0,4):
        for j in range(0,4):
            parameter_phi[i] += yi[j]*math.sqrt(a[j]) * (1-K_ij[i][j])
        #print(parameter_phi) 
        phi_v[i] = math.exp(-math.log(Z-P*b_mix/(R*T)) + (Z-1)*b[i]/b_mix - a_mix/(b_mix*R*T)*((1/a_mix)*(2*math.sqrt(a[i])*parameter_phi[i])-b[i]/b_mix) * math.log(1+b_mix/V))
        #print(phi_v)
    #print(phi_v)
    
    fi_v = []
    for i in range(0,4):
        fi_v.append(yi[i]*P*phi_v[i])
    






    D_alphaT = [0,0,0,0]; D_aT = [0,0,0,0];  XX = [] ; D_a_mix = 0
    
    # find Hr_mix 
    for i in range(0,4):
        D_alphaT[i] = -m[i]/math.sqrt(T*Tc[i])*(1+m[i]*(1-math.sqrt(Tr[i])))
        D_aT[i] = ac[i] * D_alphaT[i] 
    for i in range(0,4):
        for j in range(0,4):
            XX.append( 0.5 / (math.sqrt(a[i]*a[j])) + ( a[j]*D_aT[i] + a[j]*D_aT[i] ) )
            D_a_mix +=  yi[i]*yi[j]*XX[i+j+3*i]
    Hr_mix = R*T* ( Z-1 - 1/(b_mix*R*T) * (a_mix - T*D_a_mix) * math.log( 1 + b_mix/V) ) 












## part two: estimate the pure saturation fugacity coefficient of each component (Φ_i_sat_pure)

# phi_sat_pure is only needed for C2H5OH and H2O in liquid phase. 






# find Pi_sat from Antoine equation 

# saturation vapor pressure is only meaningful for C2H5OH and H2O ....>
# which have critical temperatures more than operating temperature(T = 50°C)

    c1 = [1.336e2 , 3.541e1 , 8.649e1 , 6.593e1] ; c2 = [-4.735e3 , -9.662e2 , -7.931e3 , -7.228e3]
    c3 = [0,0,0,0]  ; c4 = [-2.127e1 , -4.318 , -1.025e1 , -7.177]
    c5 = [4.091e-2 , 7.932e-5 , 6.389e-6 , 4.031e-6] ; c6 = [1 , 2 , 2 , 2]
    Pi_sat = []
    for i in range(0,4):
        ln_Pi_sat = c1[i] + c2[i]/(T+c3[i]) + c4[i]*math.log(T) + c5[i]*(T)**c6[i]
        Pi_sat.append(math.exp(ln_Pi_sat))
        Pi_sat[i] = Pi_sat[i]*1e3
    # print(Pi_sat)

# first step: calc Z for each of pure components
    A_pure = [0,0]
    B_pure = [0,0]
    Z_pure = [0,0]
    V_pure = [0,0]

    for i in range(0,2):
        A_pure[i] = a[i+2]*Pi_sat[i+2]/(R*T)**2             
        B_pure[i] = b[i+2]*Pi_sat[i+2]/(R*T)
        #print(A_pure,B_pure)
        z1 = 0.4
        while True:
            fz = z1**3-z1**2+(A_pure[i]-B_pure[i]-(B_pure[i])**2)*z1-A_pure[i]*B_pure[i]
            fzD = 3*z1**2-2*z1+(A_pure[i]-B_pure[i]-(B_pure[i])**2)
            z1 = z1-fz/fzD
            if abs(fz/fzD) < 1e-6:
                break
        #print(f"the root is {z1}")
        Sum = 1-z1
        Pro = A_pure[i]*B_pure[i]/z1
        if Sum**2-4*Pro < 0:
            Z_pure[i] = z1
        else:
            k1 = math.sqrt(Sum**2-4*Pro)
            k2 = Sum
            z2 = (k1 + k2)/2
            z3 = Sum - z2
            Z_pure[i] = max(z1,z2,z3)
        # print(Z_pure[i])                                 ## second and third verifying parameters 
    #print(Z_pure)
    #print(Pi_sat)

    V_pure = [0,0]
    for j in range(0,2):
        # print(Pi_sat[j+2])
        V_pure[j] = Z_pure[j]*R*T/Pi_sat[j+2]
    #print(V_pure)



# second step: calc phi_sat_pure for EtOH and Water

    parameter_phi_sat = [0,0,0,0]
    phi_sat = [0,0]
    for i in range(0,2):
        for j in range(0,4):
            parameter_phi_sat[i] += yi[j]*math.sqrt(a[j]) * (1-K_ij[i][j])
        #print(parameter_phi) 
        phi_sat[i] = math.exp(-math.log(Z_pure[i]-Pi_sat[i+2]*b[i+2]/(R*T)) + (Z_pure[i]-1)*1 - a[i+2]/(b[i+2]*R*T)*((1/a[i+2])*(2*math.sqrt(a[i+2])*parameter_phi_sat[i])-1) * math.log(1+b[i+2]/V_pure[i]))
        #print(phi_sat)
    #print(phi_sat)
    # print(phi_v,phi_sat)
    # print(phi_v)
    # print(Pi_sat)
    # print(Z,Z_pure)
    return fi_v, phi_sat, Hr_mix