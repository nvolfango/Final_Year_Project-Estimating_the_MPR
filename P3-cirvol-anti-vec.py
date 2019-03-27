def P3_CIRvol_anti_vec(x0,s0,mu1,mu2,rho,lambd,gamma,M,N):
    """ x0: initial price
        s0: initial volatility level
        M: number of steps per day
        N: number of MC runs/simulations """
    Dt = 1/(253*M);
    L = 5*253*M+1; # Total no. of steps in the period
    x = np.zeros((N,L));  x_anti = np.zeros((N,L));
    y = np.zeros(N);      y_anti = np.zeros(N);
    s = np.zeros((N,L));  s_anti = np.zeros((N,L));
    x[:,0] = x_anti[:,0] = x0;
    s[:,0] = s_anti[:,0] = s0;
    
    for k in range(1,L):
        dz1 = np.random.normal(size=N);
        dz2 = np.random.normal(size=N);
        s[:,k] = s[:,k-1] + mu2*(lambd-s[:,k-1])*Dt + gamma*np.sqrt(abs(s[:,k-1]))*np.sqrt(Dt)*dz2;
        x[:,k] = x[:,k-1] + mu1*x[:,k-1]*Dt + s[:,k-1]*x[:,k-1]*(rho*dz1 + (1-rho)*dz2)*np.sqrt(Dt);
        s_anti[:,k] = s_anti[:,k-1] + mu2*(lambd-s_anti[:,k-1])*Dt - gamma*np.sqrt(abs(s_anti[:,k-1]))*np.sqrt(Dt)*dz2;
        x_anti[:,k] = x_anti[:,k-1] + mu1*x_anti[:,k-1]*Dt - s_anti[:,k-1]*x_anti[:,k-1]*(rho*dz1 + (1-rho)*dz2)*np.sqrt(Dt);
    
    y = np.min(np.take(x,[3*M*253-1,4*M*253-1,-1],1),1);
    y_anti = np.min(np.take(x_anti,[3*M*253-1,4*M*253-1,-1],1),1);

    y[y<x0] = 0;  y_anti[y_anti<x0] = 0;
    y[y>=x0] = 1; y_anti[y_anti>=x0] = 1;
    return(0.35*np.mean(0.5*(y+y_anti)))