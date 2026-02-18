def EIGENVALUE_AND_EIGENVECTOR_ANALYSIS_FUN(N_DOF, node, mode, dof, PLOT=True):
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/eigen.html
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/nodeEigenvector.html
    import openseespy.opensees as ops
    import numpy as np
    lambda_vals = ops.eigen('-fullGenLapack', N_DOF)
    #Lambda_vals = ops.eigen('-genBandArpack', N_DOF)
    print('Eigenvector:')
    print(ops.nodeEigenvector(node, mode, dof))     
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
        
    if PLOT == False:  
        for lam in lambda_vals:
            omega = np.sqrt(lam)
            period = (2.0 * np.pi) / omega
            PERIOD.append(period)
            frequ = 1.0 / period
    #return lam, omega, period, frequ, np.max(PERIOD)
    return np.min(PERIOD), np.max(PERIOD)