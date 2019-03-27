def P5_vec(x0,m,s,r,rp,N):
    wkspryr = 304;
    h = 6/wkspryr;
    rw = (1+r)**h-1;
    rqp = rp/4;
    b = m-0.5*s**2;
    qw1 = [12,13,13,13];
    qw2 = [12,13,12,13];
    qw = np.concatenate((qw1,qw1,qw2,qw1,qw1,qw2));
    q = np.cumsum(qw);
    intwk = wkspryr-q;
    income = np.zeros(N);
    loss = np.zeros(N);
    
    
    Z = np.random.normal(size=(N,wkspryr));
    x = np.exp(b*h + s*np.sqrt(h)*Z);
    x = x0*np.cumprod(x,1);

    PI = np.array([[0 if (val<4500 or val>9000) else 1 for val in x[row,:]] for row in range(N)]);
    PI = np.cumsum(PI,1);
    p = (np.take(PI,q-1,1) - np.append(np.zeros((N,1)),np.take(PI,q[0:-1]-1,1),1)) / qw;
    
    income = np.sum(rqp*p*(1+rw)**intwk,1);
    loss = [(x0-last)/x0 if last<4250 else 0 for last in x[:,-1]];

    pay = income-loss
    mpr = np.mean(pay)
    vpr = np.var(pay)
    
    return(mpr, vpr);