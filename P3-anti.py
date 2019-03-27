def P3_anti(x0,mu,s,N):
    y = np.zeros(N)
    y_anti = np.zeros(N)
    x = np.arange(0,3)
    x_anti = np.arange(0,3)

    for j in range(N):
        samples = np.random.normal(0,1,3)
        
        x[0] = x0*np.exp((mu-0.5*s**2)*3 + s*np.sqrt(3)*samples[0])
        x[1] = x[0]*np.exp((mu-0.5*s**2) + s*samples[1])
        x[2] = x[1]*np.exp((mu-0.5*s**2) + s*samples[2])
        
        x_anti[0] = x0*np.exp((mu-0.5*s**2)*3 - s*np.sqrt(3)*samples[0])
        x_anti[1] = x_anti[0]*np.exp((mu-0.5*s**2) - s*samples[1])
        x_anti[2] = x_anti[1]*np.exp((mu-0.5*s**2) - s*samples[2])
    
        y[j] = np.min(x)
        y_anti[j] = np.min(x_anti)

    y[y<x0] = 0
    y[y>=x0] = 1
    
    y_anti[y_anti<x0] = 0
    y_anti[y_anti>=x0] = 1
    
    pr = 0.35*0.5*(y+y_anti)
    prs = pr**2
    
    return(sum(pr), sum(prs))