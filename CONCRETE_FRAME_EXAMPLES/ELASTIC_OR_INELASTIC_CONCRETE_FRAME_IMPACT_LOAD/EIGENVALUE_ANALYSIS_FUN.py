def EIGENVALUE_ANALYSIS(N_DOF, PLOT=True):
    import openseespy.opensees as ops
    import numpy as np
    lambda_vals = ops.eigen('-fullGenLapack', N_DOF)
    #Lambda_vals = ops.eigen('-genBandArpack', N_DOF)
    
    omega1 = np.sqrt(lambda_vals[0])
    PERIOD = []
    if PLOT == True:
        print('\n')
        print('+--------------------------------------------+')
        print("|  lambda   |  omega   |  period | frequency |")
        print('+--------------------------------------------+')
        for lam in lambda_vals:
            omega = np.sqrt(lam)
            period = (2.0 * np.pi) / omega
            PERIOD.append(period)
            frequ = 1.0 / period
            print("| %5.3e | %8.4f | %7.4f | %9.4f |" % (lam, omega, period, frequ))
        print('+--------------------------------------------+')
        print('\n')
    #return lam, omega, period, frequ, np.max(PERIOD)
    return np.min(PERIOD), np.max(PERIOD)