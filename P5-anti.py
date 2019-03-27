def P5_anti(x0,m,s,r,rp,N):
    h = 6/304;
    rw = (1+r)**h-1;
    rqp = rp/4;
    b = m-0.5*s**2;
    qw1 = [12,13,13,13];
    qw2 = [12,13,12,13];
    qw = np.concatenate((qw1,qw1,qw2,qw1,qw1,qw2));
    q = np.cumsum(qw);
    intwk = 304-q;
    income = np.zeros(N);
    loss = np.zeros(N);
    income_anti = np.zeros(N);
    loss_anti = np.zeros(N);
   
    for i in range(N):
        Z = np.random.normal(size=304);
        
        x = np.exp(b*h + s*np.sqrt(h)*Z);
        x = x0*np.cumprod(x);
        x_anti = np.exp(b*h - s*np.sqrt(h)*Z);
        x_anti = x0*np.cumprod(x_anti);
        
        PI = np.array(list(map(lambda val: 0 if val<4500 or val>9000 else 1,x)));
        PI = np.cumsum(PI);
        PI_anti = np.array(list(map(lambda val: 0 if val<4500 or val>9000 else 1,x_anti)));
        PI_anti = np.cumsum(PI_anti);
        
        p = (np.take(PI,q-1) - np.append([0],np.take(PI,q[0:-1]-1))) / qw;
        p_anti = (np.take(PI_anti,q-1) - np.append([0],np.take(PI_anti,q[0:-1]-1))) / qw;
        
        income[i] = np.sum(rqp*p*(1+rw)**intwk);
        income_anti[i] = np.sum(rqp*p_anti*(1+rw)**intwk);
        
        I = 1 if (x[-1]<4250) else 0;
        loss[i] = I*(x0-x[-1])/x0;
        
        I_anti = 1 if (x_anti[-1]<4250) else 0;
        loss_anti[i] = I*(x0-x_anti[-1])/x0;
    
    
    P5_pr = income-loss;
    P5_pr_anti = income_anti-loss_anti;
    
    pr = 0.5*(P5_pr + P5_pr_anti)
    prs = mpr**2
    
    return(0.5*(P5_pr+P5_pr_anti));