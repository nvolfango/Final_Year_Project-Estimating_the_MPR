def P4_vec(m,s,r,N):
    h = 1/253;
    b = m-0.5*s**2;
    pay = np.zeros(N);
    
    Z = np.random.normal(size=(N,30));
    t = np.tile([1-4*h,h,h,h,h],(N,6));
    x = np.exp(b*t + s*np.sqrt(t)*Z);
    x = np.cumprod(x, axis=1);
    x = np.cumsum(x, axis=1);
    PIL = (np.take(x,[4,9,14,19,24,29],1) - np.append(np.zeros((N,1)),np.take(x,[4,9,14,19,24],1),1)) / 5;

    cyr = np.zeros(N);
    for i in range(6):
        pay = [pay[count] + 0.05*(cyr[count]+1)*(1+r)**(5-i) if (pil>0.9) else pay[count] for count,pil in enumerate(PIL[:,i])];
        cyr = [0 if (pil>0.9) else cyr[count]+1 for count,pil in enumerate(PIL[:,i])];
    
    mpr = np.mean(pay)
    vpr = np.var(pay)
    return(mpr, vpr);