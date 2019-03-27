def P3_anti_vec(x0,mu,s,N):
    y = np.zeros(N); y_anti = np.zeros(N);
    x = x0*np.ones((N,3)); x_anti = x0*np.ones((N,3))
    z = np.random.normal(size=(N,3));
    
    x[:,0] = x0*np.exp((mu-0.5*s**2)*3 + s*np.sqrt(3)*z[:,0]);
    x[:,1] = x[:,0]*np.exp((mu-0.5*s**2) + s*z[:,1]);
    x[:,2] = x[:,1]*np.exp((mu-0.5*s**2) + s*z[:,2]);

    x_anti[:,0] = x0*np.exp((mu-0.5*s**2)*3 - s*np.sqrt(3)*z[:,0]);
    x_anti[:,1] = x_anti[:,0]*np.exp((mu-0.5*s**2) - s*z[:,1]);
    x_anti[:,2] = x_anti[:,1]*np.exp((mu-0.5*s**2) - s*z[:,2]);

    y = np.min(x, axis=1); y_anti = np.min(x_anti, axis=1);
    y[y<x0] = 0;
    y[y>=x0] = 1;
    y_anti[y_anti<x0] = 0;
    y_anti[y_anti>=x0] = 1;
    
    pr = 0.35*0.5*(y+y_anti);
    mpr = np.mean(pr)
    vpr = np.var(pr)
    
    return(mpr, vpr);