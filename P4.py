def P4(m,s,r,N):
    h = 1/253
    b = m-0.5*s**2
    pay = np.zeros(N)
    
    for i in range(N):
        Z = np.random.normal(size=30)    # 6 years, 5 business days
        t = np.tile([1-4*h,h,h,h,h],6)
        x = np.exp(b*t + s*np.sqrt(t)*Z)
        x = np.cumprod(x)               # Index path levels
        x = np.cumsum(x)
        PIL = (np.take(x,[4,9,14,19,24,29]) - np.append([0],np.take(x,[4,9,14,19,24]))) / 5
        
        cyr = 0
        for j in range(6):
            if PIL[j] > 0.9:
                pay[i] = pay[i] + 0.05*(cyr+1)*(1+r)**(5-j)
                cyr = 0
            else:
                cyr = cyr+1
    
    return(np.mean(pay))