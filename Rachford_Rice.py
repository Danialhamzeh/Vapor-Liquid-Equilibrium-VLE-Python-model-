def Rachford(zi, K_val): 


    # Find Vapor Volume Fraction (Beta)
    
    Beta = 0.1                                     # initial guess for vapor volume fraction
    while True:
        fRR = 0
        fRRD = 0
        for i in range(0,4):
            fRR = fRR + (zi[i]*(K_val[i]-1))/(1+Beta*(K_val[i]-1))
            fRRD = fRRD - (zi[i]*(K_val[i]-1)**2)/(1+Beta*(K_val[i]-1))**2
        Beta = Beta - fRR/fRRD
            # print(- fRR/fRRD)
        if abs(- fRR/fRRD) < 1e-10:
            break
            # print(Beta)


        
        

    # find mole fractions in liquid and vapor phases
    # xi : mole fractions in liquid phase      yi = mole fractions in vapor phase        
    xi = [0,0,0,0]
    yi = [0,0,0,0]
    for i in range (0,4):
        xi[i] = zi[i]/(1+Beta*(K_val[i]-1))
        yi[i] = xi[i]*K_val[i]
        # print(xi,yi,Beta)
    return Beta, xi, yi 





#zi = [0.25,0.05,0.4,0.3]
#K_val = [36.33,2482,1.791e-3,7.281e-003]
#[Beta,xi,yi] = Rachford(zi, K_val)
#print(yi)

