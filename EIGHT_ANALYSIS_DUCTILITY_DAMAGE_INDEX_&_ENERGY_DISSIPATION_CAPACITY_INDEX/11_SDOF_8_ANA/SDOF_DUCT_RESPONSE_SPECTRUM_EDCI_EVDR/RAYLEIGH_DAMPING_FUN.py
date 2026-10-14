def RAYLEIGH_DAMPING(Nmodes, zeta_i, zeta_j, target_mode_i, target_mode_j): 
    import openseespy.opensees as ops
    import numpy as np
    
    # Define the number of modes to consider in the eigenvalue analysis
    #Nmodes = 10  # This can be adjusted based on your structural model
    
    # Perform eigenvalue analysis to get the squared frequencies (w2) of the first Nmodes
    Lambda = ops.eigen('-fullGenLapack', Nmodes)
    #Lambda = ops.eigen('-genBandArpack',Nmodes)
    
    # Select two target modes and their desired damping ratios
    # Mode 1: 3% damping
    #target_mode_i = 0       # First mode (0-based index)
    Omega_I = np.sqrt(Lambda[target_mode_i])  # Natural frequency of mode 1 (rad/sec)
    #zeta_i = 0.03           # 3% damping ratio for mode 1
    
    # Mode 3: 2% damping (you can choose any other mode)
    #target_mode_j = 2       # Third mode (0-based index)
    Omega_J = np.sqrt(Lambda[target_mode_j])  # Natural frequency of mode 3 (rad/sec)
    #zeta_j = 0.02           # 2% damping ratio for mode 3
    
    # Construct the coefficient matrix for the Rayleigh damping equation
    # The system relates alpha and beta coefficients to damping ratios at two frequencies
    A = np.array([
        [1/Omega_I, Omega_I],  # First row for mode i
        [1/Omega_J, Omega_J]   # Second row for mode j
    ])
    
    # Right-hand side vector (2*zeta for each mode)
    b = np.array([2*zeta_i, 2*zeta_j])
    
    # Solve the linear system to find alpha (x[0]) and beta (x[1]) coefficients
    x = np.linalg.solve(A, b)
    
    # Apply Rayleigh damping to the model using the computed coefficients
    # Arguments are: alpha_M (mass proportional), beta_K (stiffness proportional)
    # beta_Kinit and beta_Kcomm are set to 0.0 here (not used in basic Rayleigh damping)
    ops.rayleigh(x[0], 0.0, 0.0, x[1])
    
    PERIOD_01 = (np.pi * 2) / Omega_I # Structure First Period
    PERIOD_02 = (np.pi * 2) / Omega_J # Structure Second Period
    print('Structure First Period:  ', PERIOD_01)
    print('Structure Second Period: ', PERIOD_02) 
    return PERIOD_01, PERIOD_02