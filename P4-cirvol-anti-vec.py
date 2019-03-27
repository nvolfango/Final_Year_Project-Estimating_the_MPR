def P4_CIRvol_anti_vec(s0,mu1,mu2,rho,lambd,gamma,r,M,N):
    """ x0: initial price
        s0: initial volatility level
        M: number of steps per day
        N: number of MC runs/simulations """
    Dt = 1/(253*M);
    L = 6*253*M+1; # Total no. of steps in the period
    x = np.zeros((N,L)); x_anti = np.zeros((N,L));
    s = np.zeros((N,L)); s_anti = np.zeros((N,L));
    x[:,0] = x_anti[:,0] = 1;
    s[:,0] = s_anti[:,0] = s0;
    
    pay = np.zeros(N); pay_anti = np.zeros(N);

    for k in range(1,L):
        dz1 = np.random.normal(size=N);
        dz2 = np.random.normal(size=N);
        s[:,k] = s[:,k-1] + mu2*(lambd-s[:,k-1])*Dt + gamma*np.sqrt(abs(s[:,k-1]))*np.sqrt(Dt)*dz2;
        x[:,k] = x[:,k-1] + mu1*x[:,k-1]*Dt + s[:,k-1]*x[:,k-1]*(rho*dz1 + (1-rho)*dz2)*np.sqrt(Dt);
        s_anti[:,k] = s_anti[:,k-1] + mu2*(lambd-s_anti[:,k-1])*Dt - gamma*np.sqrt(abs(s_anti[:,k-1]))*np.sqrt(Dt)*dz2;
        x_anti[:,k] = x_anti[:,k-1] + mu1*x_anti[:,k-1]*Dt - s_anti[:,k-1]*x_anti[:,k-1]*(rho*dz1 + (1-rho)*dz2)*np.sqrt(Dt);

    x = np.cumsum(x,1); x_anti = np.cumsum(x_anti,1);
    PIL = (np.take(x,[1*253*M,2*253*M,3*253*M,4*253*M,5*253*M,6*253*M],1) - np.take(x,[1*253*M-5,2*253*M-5,3*253*M-5,4*253*M-5,5*253*M-5,6*253*M-5],1)) / 5;
    PIL_anti = (np.take(x_anti,[1*253*M,2*253*M,3*253*M,4*253*M,5*253*M,6*253*M],1) - np.take(x_anti,[1*253*M-5,2*253*M-5,3*253*M-5,4*253*M-5,5*253*M-5,6*253*M-5],1)) / 5;
    
    cyr = np.zeros(N); cyr_anti = np.zeros(N);
    for i in range(6):
        pay = np.array([pay[count] + 0.05*(cyr[count]+1)*(1+r)**(5-i) if (pil>0.9) else pay[count] for count,pil in enumerate(PIL[:,i])]);
        cyr = [0 if (pil>0.9) else cyr[count]+1 for count,pil in enumerate(PIL[:,i])];
        pay_anti = np.array([pay_anti[count] + 0.05*(cyr_anti[count]+1)*(1+r)**(5-i) if (pil>0.9) else pay_anti[count] for count,pil in enumerate(PIL_anti[:,i])]);
        cyr_anti = [0 if (pil>0.9) else cyr_anti[count]+1 for count,pil in enumerate(PIL_anti[:,i])];

    return(np.mean(0.5*(pay+pay_anti)));