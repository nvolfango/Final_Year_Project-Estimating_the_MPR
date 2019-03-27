def P3(x0,mu,s,N):
    y = np.zeros(N)
    x = np.arange(0,3)

    for j in range(N):
        x[0] = x0*np.exp((mu-0.5*s**2)*3 + s*np.sqrt(3)*np.random.normal())
        x[1] = x[0]*np.exp((mu-0.5*s**2) + s*np.random.normal())
        x[2] = x[1]*np.exp((mu-0.5*s**2) + s*np.random.normal())
        
        y[j] = np.min(x)

    y[y<x0] = 0
    y[y>=x0] = 1
    return(0.35*np.mean(y))