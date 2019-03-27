def P5(x0,m,s,r,rp,N):
    h = 6/304
    rw = (1+r)**h-1
    rqp = rp/4
    b = m-0.5*s**2
    qw1 = [12,13,13,13]
    qw2 = [12,13,12,13]
    qw = np.concatenate((qw1,qw1,qw2,qw1,qw1,qw2))
    q = np.cumsum(qw)
    intwk = 304-q
    income = np.zeros(N)
    loss = np.zeros(N)
    
    for i in range(N):
        Z = np.random.normal(size=304);
        x = np.exp(b*h + s*np.sqrt(h)*Z)
        x = x0*np.cumprod(x)
        PI = np.array(list(map(lambda val: 0 if val<4500 or val>9000 else 1,x)))
        PI = np.cumsum(PI)
        p = (np.take(PI,q-1) - np.append([0],np.take(PI,q[0:-1]-1))) / qw
        income[i] = np.sum(rqp*p*(1+rw)**intwk)
        
        I = 1 if (x[-1]<4250) else 0
        loss[i] = I*(x0-x[-1])/x0
    
    meanincome = np.mean(income)
    meanloss = np.mean(loss)
    P5_mpr = meanincome-meanloss
    
    return(P5_mpr)