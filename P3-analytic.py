def analyticP3(steps):
    # Integral limits
    x3a = -0.8227
    x3b = 4.
    ua = -3.525
    ub = 0.
    
    xdt = (x3b-x3a)/steps
    udt = (ub-ua)/steps

    x3Int = np.zeros(steps)
    x3Seq = np.linspace(x3a,x3b,steps)

    for m in range(steps):
        x3 = x3Seq[m]

        uSeq = np.linspace(ua,ub,steps)
        uInt = (1/np.sqrt(2*np.pi))*np.exp(-0.5*(-np.sqrt(3)*x3-1.9-uSeq)**2)*N(uSeq-0.475)

        x3Int[m] = (1/np.sqrt(2*np.pi))*np.exp(-0.5*x3**2)*(N(-np.sqrt(3)*x3-1.9) + sum(uInt)*udt)

    MPR = 0.2781312 - 0.35*sum(x3Int)*xdt
    return(MPR)