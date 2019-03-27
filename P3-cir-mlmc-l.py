def P3_CIR_anti_mlmc(x0, s0, mu1, mu2, rho, lambd, gamma, N0, method, mlmc_l, M, eps, L):
#
# Initialisation
#
    theta = 0.75
    Nl = np.zeros(L+1)
    suml = np.zeros((2,L+1))
    dNl = np.ones(L+1)*N0
    
    while sum(dNl) > 0.02*sum(Nl):
#
# Update sample sums
#
        for l in range(L+1):
            if dNl[l] > 0:
                sums = P3_CIR_anti_mlmc_l(x0,s0,mu1,mu2,rho,lambd,l,gamma,method,M,dNl[l])
                Nl[l] += dNl[l]
                suml[0,l] += sums[0]
                suml[1,l] += sums[1]
# Compute absolute average and variance     
#
        ml = abs(suml[0,:]/Nl)
        Vl = np.array([max(0,val) for val in (suml[1,:]/Nl - ml**2)])
#
# Set optimal number of additional samples
#
        Ns = np.ceil( np.sqrt(Vl/M**l)*sum(np.sqrt(Vl*M**l)) / ((1-theta)*eps**2) ) 
        dNl = Ns-Nl
        dNl = np.array([max(0,val) for val in dNl])
#
# Finally, evaluate the multilevel estimator
#
    C = sum(Nl*M**np.arange(len(Nl)))
    P = sum(suml[0,:]/Nl) # Final expected value
    V = sum(Vl/Nl)        # Variance of estimator
    
    return P, V