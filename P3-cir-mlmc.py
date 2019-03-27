def P3_CIR_anti_mlmc_l(x0,s0,mu1,mu2,rho,lambd,l,gamma,method,M,N0):
    """ x0: initial price
        s0: initial volatility level
        N: number of MC runs/simulations """
    N = int(N0/2)
    
    # Lowest level:
    if l == 0:
        dt = 5/(5*(M**l))
        L = 5*(M**l) + 1
        x = np.zeros((N,L))
        x_anti = np.zeros((N,L))
        s = np.zeros((N,L))
        s_anti = np.zeros((N,L))
        x[:,0] = x_anti[:,0] = x0
        s[:,0] = s_anti[:,0] = s0
        mil_cor = 0 if method == 'euler-maruyama' else 1 # Milstein method correction term
        dW1 = np.random.normal(0,np.sqrt(dt),(N,L-1))
        dW2 = np.random.normal(0,np.sqrt(dt),(N,L-1))

        for k in range (1,L):
            s[:,k] = s[:,k-1] + mu2*(lambd-s[:,k-1])*dt + gamma*np.sqrt(abs(s[:,k-1]))*dW2[:,k-1] + (0.5**2)*mil_cor*(gamma**2)*(dW2[:,k-1]**2 - dt)
            x[:,k] = x[:,k-1] + mu1*x[:,k-1]*dt + s[:,k-1]*x[:,k-1]*(rho*dW1[:,k-1] + (1-rho)*dW2[:,k-1]) + 0.5*mil_cor*(s[:,k-1]**2)*x[:,k-1]*(rho**2)*(dW1[:,k-1]**2 - dt) + 0.5*mil_cor*(s[:,k-1]**2)*x[:,k-1]*((1-rho)**2)*(dW2[:,k-1]**2 - dt)
            s_anti[:,k] = s_anti[:,k-1] + mu2*(lambd-s_anti[:,k-1])*dt - gamma*np.sqrt(abs(s_anti[:,k-1]))*dW2[:,k-1] + (0.5**2)*mil_cor*(gamma**2)*(dW2[:,k-1]**2 - dt)
            x_anti[:,k] = x_anti[:,k-1] + mu1*x_anti[:,k-1]*dt + s_anti[:,k-1]*x_anti[:,k-1]*(-rho*dW1[:,k-1] - (1-rho)*dW2[:,k-1]) + 0.5*mil_cor*(s_anti[:,k-1]**2)*x_anti[:,k-1]*(rho**2)*(dW1[:,k-1]**2 - dt) + 0.5*mil_cor*(s_anti[:,k-1]**2)*x_anti[:,k-1]*((1-rho)**2)*(dW2[:,k-1]**2 - dt)

        y = np.min(np.take(x,[3*(M**l),4*(M**l),-1],1),1)
        y_anti = np.min(np.take(x_anti,[3*(M**l),4*(M**l),-1],1),1)
    
        y[y<x0] = 0
        y[y>=x0] = 1
        y_anti[y_anti<x0] = 0
        y_anti[y_anti>=x0] = 1
        
        xs = 0.35*y
        xs_anti = 0.35*y_anti
            
    else:
        dtf = 5/(5*(M**l))
        dtc = 5/(5*(M**(l-1)))
        Lf = 5*(M**l) + 1
        Lc = 5*(M**(l-1)) + 1
        xf = np.zeros((N,Lf))
        xc = np.zeros((N,Lc))
        sf = np.zeros((N,Lf))
        sc = np.zeros((N,Lc))
        xf_anti = np.zeros((N,Lf))
        xc_anti = np.zeros((N,Lc))
        sf_anti = np.zeros((N,Lf))
        sc_anti = np.zeros((N,Lc))
        xf[:,0] = xc[:,0] = xf_anti[:,0] = xc_anti[:,0] = x0
        sf[:,0] = sc[:,0] = xf_anti[:,0] = xc_anti[:,0] = s0
        mil_cor = 0 if method == 'euler-maruyama' else 1 # Milstein method correction term

        dW1f = np.random.normal(0,np.sqrt(dtf),(N,Lf-1))
        dW1c = np.zeros((N,int(dW1f.shape[1]/M)))
        dW2f = np.random.normal(0,np.sqrt(dtf),(N,Lf-1))
        dW2c = np.zeros((N,int(dW2f.shape[1]/M)))

        for j in range(N):
            dW1c[j,:] = [sum(dW1f[j,i-M+1:i+1]) for i in range(M-1,dW1f.shape[1],M)]
            dW2c[j,:] = [sum(dW2f[j,i-M+1:i+1]) for i in range(M-1,dW2f.shape[1],M)]

        for k in range(1,Lf):
            sf[:,k] = sf[:,k-1] + mu2*(lambd-sf[:,k-1])*dtf + gamma*np.sqrt(abs(sf[:,k-1]))*dW2f[:,k-1] + (0.5**2)*mil_cor*(gamma**2)*(dW2f[:,k-1]**2 - dtf)
            xf[:,k] = xf[:,k-1] + mu1*xf[:,k-1]*dtf + sf[:,k-1]*xf[:,k-1]*(rho*dW1f[:,k-1] + (1-rho)*dW2f[:,k-1]) + 0.5*mil_cor*(sf[:,k-1]**2)*xf[:,k-1]*(rho**2)*(dW1f[:,k-1]**2 - dtf) + 0.5*mil_cor*(sf[:,k-1]**2)*xf[:,k-1]*((1-rho)**2)*(dW2f[:,k-1]**2 - dtf)
            sf_anti[:,k] = sf_anti[:,k-1] + mu2*(lambd-sf_anti[:,k-1])*dtf - gamma*np.sqrt(abs(sf_anti[:,k-1]))*dW2f[:,k-1] + (0.5**2)*mil_cor*(gamma**2)*(dW2f[:,k-1]**2 - dtf)
            xf_anti[:,k] = xf_anti[:,k-1] + mu1*xf_anti[:,k-1]*dtf + sf_anti[:,k-1]*xf_anti[:,k-1]*(-rho*dW1f[:,k-1] - (1-rho)*dW2f[:,k-1]) + 0.5*mil_cor*(sf_anti[:,k-1]**2)*xf_anti[:,k-1]*(rho**2)*(dW1f[:,k-1]**2 - dtf) + 0.5*mil_cor*(sf_anti[:,k-1]**2)*xf_anti[:,k-1]*((1-rho)**2)*(dW2f[:,k-1]**2 - dtf)
            

        for m in range(1,Lc):
            sc[:,m] = sc[:,m-1] + mu2*(lambd-sc[:,m-1])*dtc + gamma*np.sqrt(abs(sc[:,m-1]))*dW2c[:,m-1] + (0.5**2)*mil_cor*(gamma**2)*(dW2c[:,m-1]**2 - dtc)
            xc[:,m] = xc[:,m-1] + mu1*xc[:,m-1]*dtc + sc[:,m-1]*xc[:,m-1]*(rho*dW1c[:,m-1] + (1-rho)*dW2c[:,m-1]) + 0.5*mil_cor*(sc[:,m-1]**2)*xc[:,m-1]*(rho**2)*(dW1c[:,m-1]**2 - dtc) + 0.5*mil_cor*(sc[:,m-1]**2)*xc[:,m-1]*((1-rho)**2)*(dW2c[:,m-1]**2 - dtc)
            sc_anti[:,m] = sc_anti[:,m-1] + mu2*(lambd-sc_anti[:,m-1])*dtc - gamma*np.sqrt(abs(sc_anti[:,m-1]))*dW2c[:,m-1] + (0.5**2)*mil_cor*(gamma**2)*(dW2c[:,m-1]**2 - dtc)
            xc_anti[:,m] = xc_anti[:,m-1] + mu1*xc_anti[:,m-1]*dtc + sc_anti[:,m-1]*xc_anti[:,m-1]*(-rho*dW1c[:,m-1] - (1-rho)*dW2c[:,m-1]) + 0.5*mil_cor*(sc_anti[:,m-1]**2)*xc_anti[:,m-1]*(rho**2)*(dW1c[:,m-1]**2 - dtc) + 0.5*mil_cor*(sc_anti[:,m-1]**2)*xc_anti[:,m-1]*((1-rho)**2)*(dW2c[:,m-1]**2 - dtc)

        yf = np.min(np.take(xf,[3*(M**l),4*(M**l),-1],1),1)
        yc = np.min(np.take(xc,[3*(M**(l-1)),4*(M**(l-1)),5*(M**(l-1))],1),1)
        yf_anti = np.min(np.take(xf_anti,[3*(M**l),4*(M**l),-1],1),1)
        yc_anti = np.min(np.take(xc_anti,[3*(M**(l-1)),4*(M**(l-1)),5*(M**(l-1))],1),1)
        
        yf[yf<x0] = 0
        yf[yf>=x0] = 1
        yc[yc<x0] = 0
        yc[yc>=x0] = 1
        
        yf_anti[yf_anti<x0] = 0
        yf_anti[yf_anti>=x0] = 1
        yc_anti[yc_anti<x0] = 0
        yc_anti[yc_anti>=x0] = 1
        
        Pf = 0.35*yf
        Pc = 0.35*yc
        Pf_anti = 0.35*yf_anti
        Pc_anti = 0.35*yc_anti
        
        xs = Pf-Pc
        xs_anti = Pf_anti-Pc_anti
        
    return np.array([sum(xs+xs_anti), sum(xs**2 + xs_anti**2)])