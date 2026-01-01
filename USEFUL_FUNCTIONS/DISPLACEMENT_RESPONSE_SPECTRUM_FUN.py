#%% Displacement Response Spectrum
def DISPLACEMENT_RESPONSE_SPECTRUM(SPRING_KIND, acc, dt, periods, ductility_ratio, damping_ratio=0.05):
    import openseespy.opensees as ops
    import numpy as np
    import ANALYSIS_FUNCTION as S01
    import OPENSEEES_HYSTERETICSM_FORCE_DISP_FUN as S02
    # DEFINE ANALYSIS PROPERTIES
    MAX_ITERATIONS = 20000    # Convergence iteration for test
    MAX_TOLERANCE = 1.0e-10   # Convergence tolerance for test
    
    spectral_disps = []
    spectral_accels = []
    
    acc_list = list(acc)
    n_steps = len(acc_list)
    
    for T in periods:
        ops.wipe()
        
        # Create model
        ops.model('basic', '-ndm', 1, '-ndf', 1)
        
        # Nodes
        ops.node(1, 0.0)
        ops.node(2, 0.0)
        
        # Fix base
        ops.fix(1, 1)
        
        # Calculate parameters
        omega = 2 * np.pi / T
        mass = 1.0
        stiffness = mass * omega**2
        print('stiffness', stiffness)
        print('Period', T)
        
        # Define  Structural Properties
        FY = 8.5                                         # [N] Yield Force of Structure
        FU = 1.5 * FY                                    # [N] Ultimate Force of Structure
        #Ke = 4500000.0                                   
        Ke = stiffness                                   # [N/m] Spring Elastic Stiffness
        DY = FY / Ke                                     # [m] Yield Displacement
        print('DY', DY)
                                
        DSU = DY * ductility_ratio                       # [m] Ultimate Displacement
        Ksh = (FU - FY) / (DSU - DY)                     # [N/m] Displacement Hardening Modulus
        Kp = FU / DSU                                    # [N/m] Spring Plastic Stiffness
        b = Ksh / Ke                                     # Displacement Hardening Ratio
        # Define material properties
        MatTag_S = 1 # SPRING TAG
        # FORCE-DISPLACEMENT RELATIONSHIP OF LATERAL SPRING AND PLOT 
        DP = [0, 0, 0, 0]
        FP = [0, 0, 0, 0]
        DN = [0, 0, 0, 0]
        FN = [0, 0, 0, 0]
        #print(DSU,"------------" ,FU)
        DP[0], FP[0] = DY, FY
        DP[1], FP[1] = DSU, FU 
        DP[2], FP[2] = 1.1*DSU, 0.20*FU
        DP[3], FP[3] = 1.25*DSU, 0.10*FU
        DN[0], FN[0] = -DY, -FY 
        DN[1], FN[1] = -DSU, -FU
        DN[2], FN[2] = -1.1*DSU, -0.20*FU   
        DN[3], FN[3] = -1.25*DSU, -0.10*FU
        
        c = 2 * damping_ratio * omega * mass  # Damping coefficient
        
        # Material and element
        if SPRING_KIND == 'ELASTIC':
            ops.uniaxialMaterial('Elastic', MatTag_S, stiffness) # INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Elastic_Uniaxial_Material
        if SPRING_KIND == 'INELASTIC':
            S02.OPENSEEES_HYSTERETICSM_FORCE_DISP_FUN(MatTag_S, DP, FP, DN, FN, PLOT = False, X_LABEL='Displacement (mm)', Y_LABEL='Force [N]', TITLE='FORCE-DISPLACEMENT CURVE')
        ops.uniaxialMaterial('Viscous', 2, c, 1)      # INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Viscous_Material
        ops.element('zeroLength', 1, 1, 2, '-mat', 1, 2, '-dir', 1, 1)
        
        # Mass
        ops.mass(2, mass)
        
        # Analysis setup
        ops.constraints('Plain')
        ops.numberer('Plain')
        ops.system('BandGeneral')
        ops.test('NormDispIncr', 1.0e-8, 10)
        ops.algorithm('Newton')
        
        # Create recorder to track response
        ops.recorder('Node', '-file', f'nodeDisp_T_{T:.3f}.txt', '-node', 2, '-dof', 1, 'disp')
        ops.recorder('Node', '-file', f'nodeAccel_T_{T:.3f}.txt', '-node', 2, '-dof', 1, 'accel')
        
        # Time series and pattern
        ops.timeSeries('Path', 1, '-dt', dt, '-values', *acc_list)
        ops.pattern('UniformExcitation', 1, 1, '-accel', 1)
        
        # Damping
        ops.rayleigh(0.0, 0.0, 0.0, 2*damping_ratio/omega)
        
        # Set up analysis
        ops.integrator('Newmark', 0.5, 0.25)
        ops.analysis('Transient')
        
        # Initialize
        u_max = 0.0
        a_max = 0.0
        
        # Perform analysis
        current_time = 0.0
        for i in range(n_steps):
            OK = ops.analyze(1, dt)
            S01.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS
            current_time += dt
            
            # Get responses
            u = ops.nodeDisp(2, 1)
            
            # For acceleration in support excitation:
            # We can approximate total acceleration from equilibrium
            ops.reactions()
            R = ops.nodeReaction(2, 1)
            a_relative = -R / mass
            a_total = a_relative - acc_list[i]
            
            # Track maxima
            u_max = max(u_max, abs(u))
            a_max = max(a_max, abs(a_total))
        
        spectral_disps.append(u_max)
        spectral_accels.append(a_max)
        """
        # Remove temporary files
        import os
        if os.path.exists(f'nodeDisp_T_{T:.3f}.txt'):
            os.remove(f'nodeDisp_T_{T:.3f}.txt')
        if os.path.exists(f'nodeAccel_T_{T:.3f}.txt'):
            os.remove(f'nodeAccel_T_{T:.3f}.txt')
        """
    return np.array(spectral_disps), np.array(spectral_accels)