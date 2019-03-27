def P5_CIRvol_vec(x0,s0,mu1,mu2,rho,lambd,gamma,r,rp,M,N):
    """ x0: initial price
        s0: initial volatility level
        M: number of steps per week
        N: number of MC runs/simulations """
    wkspryr = 304;
    h = 6/wkspryr;
    rw = (1+r)**h-1;
    rqp = rp/4;
    qw1 = [12,13,13,13];
    qw2 = [12,13,12,13];
    qw = np.concatenate((qw1,qw1,qw2,qw1,qw1,qw2));
    q = np.cumsum(qw);
    intwk = wkspryr-q;
    
    income = np.zeros(N);
    loss = np.zeros(N);
    
    Dt = 6/(wkspryr*M);
    L = wkspryr*M+1;
    x = np.zeros((N,L));
    s = np.zeros((N,L));
    xs = np.zeros((N,wkspryr));
    x[:,0] = x0;
    s[:,0] = s0;
    
    for k in range(1,L):
        dz1 = np.random.normal(size=N);
        dz2 = np.random.normal(size=N);
        s[:,k] = s[:,k-1] + mu2*(lambd-s[:,k-1])*Dt + gamma*np.sqrt(abs(s[:,k-1]))*np.sqrt(Dt)*dz2;
        x[:,k] = x[:,k-1] + mu1*x[:,k-1]*Dt + s[:,k-1]*x[:,k-1]*(rho*dz1 + (1-rho)*dz2)*np.sqrt(Dt);
        
    xs = np.take(x,np.arange(M-1,L-1,M),1);
    PI = np.array([[0 if val<4500 or val>9000 else 1 for val in xs[row,:]] for row in range(N)]);
    PI = np.cumsum(PI,1);
    p = (np.take(PI,q-1,1) - np.append(np.zeros((N,1)),np.take(PI,q[0:-1]-1,1),1)) / qw;
    
    income = np.sum(rqp*p*(1+rw)**intwk,1);
    loss = [(x0-last)/x0 if last<4250 else 0 for last in xs[:,-1]];
        
    meanincome = np.mean(income);
    meanloss = np.mean(loss);
    P5_mpr = meanincome-meanloss;
    
    return(P5_mpr);