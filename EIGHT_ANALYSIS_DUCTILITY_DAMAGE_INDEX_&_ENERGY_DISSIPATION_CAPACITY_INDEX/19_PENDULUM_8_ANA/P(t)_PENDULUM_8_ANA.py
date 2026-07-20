###########################################################################################################
#                   >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                      #
#         COMPREHENSIVE NONLINEAR SEISMIC ASSESSMENT OF A PENDULUM STRUCTURE : AN OPENSEES FRAMEWORK      #
#    FOR STATIC PUSHOVER, CYCLIC DEGRADATION, STATIC TIME-HISTORY AND DYNAMIC TIME-HISTORY ANALYSIS       #
#---------------------------------------------------------------------------------------------------------#
#                                                  P(t) = P0 sin(wt)                                      #
#                                           P(t) = P0 exp(-0.05wt) sin(wt)                                #
#---------------------------------------------------------------------------------------------------------#
#         ASSESSMENT OF DUCTILITY DAMAGE INDICES FOR STRUCTURAL ELEMENTS AND SYSTEMS AND EVALUATION       #
#                                        OF ENERGY DISSIPATION CAPACITY INDEX                             #
#---------------------------------------------------------------------------------------------------------#
#                    THIS PYTHON SCRIPT WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                     #
#                                   EMAIL: salar.d.ghashghaei@gmail.com                                   #
###########################################################################################################
"""
Nonlinear Seismic Performance Assessment of a PENDULUM:
An OpenSeesPy Framework for Material and Geometric Nonlinearity Under Static, Cyclic, and Earthquake Loading

This OpenSeesPy script performs rigorous nonlinear static and dynamic analysis of a PENDULUM
 for performance-based earthquake engineering. The 2D model incorporates both material nonlinearity
 (elastic-perfectly plastic or hysteretic steel with strain hardening)
 and geometric nonlinearity (corotational truss formulation) to capture P-delta effects
 and large displacements—essential for collapse assessment.
 
Eight analysis protocols are implemented:    
(1) [PERIOD] : Structural Period
(2) [STATIC] : Gravity load analysis establishing dead load state
(3) [PUSHOVER] : Displacement-controlled pushover generating full capacity curves
 and plastic mechanism identification
(4) [CYCLIC_DISPLACEMENT] : Symmetric cyclic displacement protocol capturing hysteresis,
 pinching behavior, and energy dissipation degradation
(5) [STATIC_EXTERNAL_TIME-DEPENDENT_LOADING] : Static Analysis of External time-dependent loading P(t) = P0 sin(wt) or P(t) = P0 exp(-0.05wt) sin(wt)  
(6) [DYNAMIC_EXTERNAL_TIME-DEPENDENT_LOADING] : Dynamic Analysis of External time-dependent loading P(t) = P0 sin(wt) or P(t) = P0 exp(-0.05wt) sin(wt)  
(7) [FREE-VIBRATION] : Free-vibration with initial conditions extracting damping ratios
 via logarithmic decrement
(8) [SEISMIC] : Multi-directional seismic excitation with Rayleigh damping (3% ratio)
 and uniform acceleration patterns.

The code continuously records displacements, enabling member-level demand-to-capacity
 ratio tracking. Period monitoring during inelastic response reveals softening
 and potential period shift into resonant ground motion frequency ranges.
 This framework enables vulnerability curve development, seismic fragility assessment,
 and retrofit prioritization—directly supporting next-generation bridge seismic design codes
 and risk-informed asset management decisions.
 
------------------------------------------------ 
Energy Dissipation Capacity Index (EDCI):    
The Energy Dissipation Capacity Index is a quantitative
 measure used in structural engineering to evaluate how
 effectively a structural element (e.g., a beam, column, shear wall, or connection)
 can absorb and dissipate energy during seismic loading compared to its
 performance under controlled cyclic displacement loading.
It compares the actual energy absorbed during an earthquake
 with the maximum energy dissipation capacity that the component
 demonstrates in a laboratory‑style cyclic test.  
 
Why This Index Is Important:
During an earthquake, structures undergo repeated cycles of deformation.
 A system with high energy dissipation capacity can withstand more damage
 without collapsing because it can convert seismic input energy into hysteretic energy, not elastic rebound.

The EDCI helps engineers understand:
[1] Ductility performance
[2] Hysteretic behavior
[3] Damage tolerance
[4] Collapse prevention capability
It is especially used in performance‑based seismic evaluation and retrofit design. 
"""
#%%----------------------------------------------------
import numpy as np
import openseespy.opensees as ops
import matplotlib.pyplot as plt 
import PLOT_2D as S01
import ANALYSIS_FUNCTION as S02
import PERIOD_FUN as S03
import DAMPING_RATIO_FUN as S04
import EIGENVALUE_ANALYSIS_FUN as S05
import RAYLEIGH_DAMPING_FUN as S06
import BILINEAR_CURVE as S07
import OPENSEEES_HYSTERETICSM_FORCE_DISP_FUN as S08
#%%----------------------------------------------------
def PENDULUM(SPRING_TYPE, MAT_TYPE, TOTAL_MASS):
    # Initialize model
    ops.wipe()
    # Create a 2D model with 3 DOF per node
    ops.model('Basic', '-ndm', 2, '-ndf', 3)
    
    MAX_ITERATIONS = 5000   # Maximum number of iterations
    MAX_TOLERANCE = 1.0e-6  # Specified tolerance for convergence
    L = 1.0                 # [m] Element Length
    GMfact = 9.81           # [m/s^2] standard acceleration of gravity or standard acceleration 
    DR = 0.05               # Damping ratio
      
    # Moment-Rotation Relation for Rotational Spring 
    MY = 12000.0                      # [N.m] Yield Moment of Structure
    MU = 1.05 * MY                    # [N.m] Ultimate Moment of Structure
    KeI = 820000.0                    # [N.m/rad] Spring Elastic Stiffness
    TY = MY / KeI                     # [rad] Yield Rotation
    TSU = 0.36                        # [rad] Ultimate Rotation
    KshI = (MU - MY) / (TSU - TY)     # [N.m/rad] Rotation Hardening Modulus
    KpI = MU / TSU                    # [N.m/rad] Spring Plastic Stiffness
    b = KshI / KeI                    # Rotation Hardening Ratio
    # Positive branch points
    pos_disp = [0, TY, TSU, 1.1*TSU, 1.25*TSU]
    pos_force = [0, MY, MU, 0.2*MU, 0.1*MU]
    KPI = np.array([MY, TY, MU, TSU, 0.2*MU, 1.1*TSU, 0.1*MU, 1.25*TSU])
    
    # Negative branch points
    neg_disp = [0, -TY, -TSU, -1.1*TSU, -1.25*TSU]
    neg_force = [0, -MY, -MU, -0.2*MU, -0.1*MU]
    KNI = np.array([-MY, -TY, -MU, -TSU, -0.2*MU, -1.1*TSU, -0.1*MU, -1.25*TSU])
    
    # Plot
    plt.plot(pos_disp, pos_force, marker='o', color='red')
    plt.plot(neg_disp, neg_force, marker='o', color='black')
    
    plt.xlabel("Rotation [rad]")
    plt.ylabel("Moment [N.m]")
    plt.title("Moment–Rotation Curve for Rotational Spring")
    plt.grid(True)
    plt.axhline(0, linewidth=0.5)
    plt.axvline(0, linewidth=0.5)
    plt.show()

    # Add nodes
    ops.node(1, 0.0, 0.0)
    ops.node(2, 0.0, -L)
    ops.node(3, L, -L)#-1.0e-8
    
    ops.node(11, 0.0, 0.0)#-1.0e-8

    
    # Define mass to node 2
    ops.mass(2, TOTAL_MASS, TOTAL_MASS, TOTAL_MASS)
    
    # Fix node 1
    ops.fix(1, 1, 1, 0)
    ops.fix(3, 1, 1, 1) 
    
    # Fix spring node
    ops.fix(11, 1, 1, 1)
        
    # Elastic Element
    B = 0.05                # [m] Element Section Width
    H = 0.05                # [m] Element Section Height
    AREA = B*H              # [m^2] Element Section Area
    Iz = B*(H**3) / 12
    Es = 200000.0/(1000)**-2 # [N/m^2] Element Young's Modulus
    DENSITY_STEEL = 7850.0   # [kg/m^3] Steel Material Density
    
    # Geometric transformation
    transfTag = 1
    #ops.geomTransf('Linear', transfTag)
    #ops.geomTransf('PDelta', transfTag)
    ops.geomTransf('Corotational', transfTag)
    ops.element('elasticBeamColumn', 100, 1, 2, AREA, Es, Iz, transfTag, '-mass', AREA*DENSITY_STEEL)
    #INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/elasticBeamColumn.html
    
    # Define material and element for Rotational Spring Kθ(M-θ)
    omega = np.sqrt(KeI/TOTAL_MASS)   # [N/m] KeI: Rotational Elastic Stiffness
    C = 2 * DR * omega * TOTAL_MASS   # Damping Coefficient
    if SPRING_TYPE == 'ELASTIC': 
        ops.uniaxialMaterial('Elastic', 1, KeI, C)
        ops.element('zeroLength', 10, 11, 1,'-mat', 1,'-dir', 3) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/ZeroLength.html
    if SPRING_TYPE == 'INELASTIC':
        ops.uniaxialMaterial('HystereticSM', 1, '-posEnv', *KPI.flatten(), '-negEnv', *KNI.flatten(), '-pinch', 1, 1)
        #INFO LINK: https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/HystereticSM.html
        ops.uniaxialMaterial('Viscous', 2, C, 1.0)  # Material for C (alpha=1.0 for linear)
        ops.element('zeroLength', 10, 11, 1,'-mat', 1, 2,'-dir', 3, 3)
        
    # Define material and element for Lateral Spring K∆(P-∆)
    Ele_Tag_03 = 300
    MAT_TAG = 5000
    if MAT_TYPE == 'ELASTIC':
        E = 200000.0/10**6                                       # [N/m²] Modulus of steel
        #ops.uniaxialMaterial('Elastic', MAT_TAG, E)             # TESNSION AND COMPRESSION IS SAME VALUES
        ops.uniaxialMaterial('Elastic', MAT_TAG, E ,0.0, 0.5*E)  # TESNSION AND COMPRESSION IS NOT SAME VALUES
        # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/ElasticUni.html
    # Material (inelastic steel)    
    if MAT_TYPE == 'INELASTIC':
        Fy = 0.240			    # [N/m²] Steel yield stress
        E = 200000.0/10**6      # [N/m²] Modulus of steel
        ey = Fy/Es			    # [m/m] Steel yield strain
        Fu = 1.1818*Fy          # [N/m²] Steel Ultimate Strength
        esu = 0.12              # [m/m] Steel Ultimate Strain
        Esh = (Fu - Fy)/(esu - ey)
        Bs = Esh / Es           # strain-hardening ratio 
        pinchX = 0.8            # Pinching factor in X direction
        pinchY = 0.5            # Pinching factor in Y direction
        damage1 = 0.0           # Damage due to ductility
        damage2 = 0.0           # Damage due to energy
        beta = 0.1              # Stiffness degradation parameter
        ops.uniaxialMaterial('Hysteretic', MAT_TAG, Fy, ey, Fu, esu, 0.2*Fu, 1.1*esu, -Fy, -ey, -Fu, -esu, -0.2*Fu, -1.1*esu, pinchX, pinchY, damage1, damage2, beta)
        # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/Hysteretic.html
    
    # INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Corotational_Truss_Element
    # element corotTruss $eleTag $iNode $jNode $A $matTag <-rho $rho> <-cMass $cFlag> <-doRayleigh $rFlag>
    ops.element('corotTruss', Ele_Tag_03, 2, 3, AREA, MAT_TAG, '-rho',AREA*DENSITY_STEEL , '-doRayleigh', 1)
    
    center_node = 2 # Target Node
    
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(2, 1.0, 0.0, 0.0)
    
    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.system('BandGeneral')
    
    if MAT_TYPE == 'ELASTIC':
        ops.algorithm('Linear')
    if MAT_TYPE == 'INELASTIC':
        ops.algorithm('Newton')
        

    time = []
    dispX, dispY, dispZ = [], [], []
    velX, velY, velZ = [], [], []
    accelX, accelY, accelZ = [], [], []
    reactionX, reactionY, reactionZ = [], [], []
    stiffnessX, stiffnessY, stiffnessZ = [], [], []
    OMEGA, PERIOD = [], []
    PERIOD_MIN, PERIOD_MAX = [], []
    
    if ANAL_TYPE == 'PERIOD': 
      #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(1, 0.5*DR, DR, 0, 1)
      PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(1, PLOT=True) 
      # Compute modal properties
      ops.modalProperties("-print", "-file", f"SALAR_ModalReport_{ANAL_TYPE}.txt", "-unorm")
      return PERIODmin, PERIODmax
  
    if ANAL_TYPE == 'STATIC': 
        ops.integrator('LoadControl', 1.0)
        ops.analysis('Static')
        OK = ops.analyze(1)
        S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS)
        ops.reactions()
        reactionX.append(ops.nodeReaction(1, 1))           # BASE REACTION IN X DIR.
        reactionY.append(ops.nodeReaction(1, 2))           # BASE REACTION IN Y DIR.
        reactionZ.append(ops.nodeReaction(1, 3))           # BASE REACTION IN Z DIR.
        dispX.append(ops.nodeDisp(center_node, 1))         # DISPLACEMENT IN X DIR.  
        dispY.append(ops.nodeDisp(center_node, 2))         # DISPLACEMENT IN Y DIR.
        dispZ.append(ops.nodeDisp(center_node, 3))         # ROTATION IN Z DIR.
        stiffnessX.append(np.abs(reactionX[-1]) / np.abs(dispX[-1]))
        stiffnessY.append(np.abs(reactionY[-1]) / np.abs(dispY[-1]))
        stiffnessZ.append(np.abs(reactionZ[-1]) / np.abs(dispZ[-1]))
        print('\n\nSTATIC ANALYSIS DONE.\n\n')  
            
        DATA = (reactionX, dispX, velX, accelX,
                reactionY, dispY, velY, accelY,
                reactionZ, dispZ, velZ, accelZ,
                stiffnessX, stiffnessY, stiffnessZ)
        
        return  DATA
    
    if ANAL_TYPE == 'PUSHOVER': # STATIC TIME-HISTORY ANALYSIS
        IDctrlDOF = 1   # 1: Horizental Dispalcement - 2: Vertical Dispalcement
        DINCR = -0.001  # [m] Incremental Vertical Displacement
        DMAX = -2.35    # [m] Max. Displacement
        ops.integrator('DisplacementControl', center_node, IDctrlDOF, DINCR)
        ops.analysis('Static')
        Nsteps =  int(np.abs(DMAX/ DINCR)) 
        STEP = 0.0
        for step in range(Nsteps):
            OK = ops.analyze(1)
            S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS)
            ops.reactions()
            reactionX.append(ops.nodeReaction(1, 1))           # BASE REACTION IN X DIR.
            reactionY.append(ops.nodeReaction(1, 2))           # BASE REACTION IN Y DIR.
            reactionZ.append(ops.nodeReaction(1, 3))           # BASE REACTION IN Z DIR.
            dispX.append(ops.nodeDisp(center_node, 1))         # DISPLACEMENT IN X DIR.  
            dispY.append(ops.nodeDisp(center_node, 2))         # DISPLACEMENT IN Y DIR.
            dispZ.append(ops.nodeDisp(center_node, 3))         # ROTATION IN Z DIR.
            velX.append(ops.nodeVel(center_node, 1))           # VELOCITY IN X DIR.
            velY.append(ops.nodeVel(center_node, 2))           # VELOCITY IN Y DIR.
            velZ.append(ops.nodeVel(center_node, 3))           # VELOCITY IN Z DIR.
            accelX.append(ops.nodeAccel(center_node, 1))       # ACCELERATION IN X DIR.
            accelY.append(ops.nodeAccel(center_node, 2))       # ACCELERATION IN Y DIR.
            accelZ.append(ops.nodeAccel(center_node, 3))       # ACCELERATION IN Z DIR.
            stiffnessX.append(np.abs(reactionX[-1]) / np.abs(dispX[-1]))
            stiffnessY.append(np.abs(reactionY[-1]) / np.abs(dispY[-1]))
            stiffnessZ.append(np.abs(reactionZ[-1]) / np.abs(dispZ[-1]))      
            # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
            #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(1, 0.5*DR, DR, 0, 1)
            PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(1, PLOT=True)
            PERIOD_MIN.append(PERIODmin)
            PERIOD_MAX.append(PERIODmax) 
            STEP += 1
            #print(STEP, disp[-1], reaction[-1])
            print(f"Step: {STEP}, Displacement: {dispX[-1]:.4f} mm, Reaction: {reactionX[-1]:.2f} N")     
        else:
            print('\n\nPUSHOVER ANALYSIS DONE.\n\n')
            
        DATA = (reactionX, dispX, velX, accelX,
                reactionY, dispY, velY, accelY,
                reactionZ, dispZ, velZ, accelZ,
                stiffnessX, stiffnessY, stiffnessZ,
                np.array(PERIOD_MIN), np.array(PERIOD_MAX))
    
        return  DATA
    
    if ANAL_TYPE == 'CYCLIC_DISPLACEMENT': # STATIC TIME-HISTORY ANALYSIS
        IDctrlDOF = 1   # 1: Horizental Dispalcement - 2: Vertical Dispalcement
        DMAX = -2.35    # [m] Max. Displacement
        n_points = 10000
        # 4. CYCLIC DISPALCEMENT PROTOCOL
        # Key strain points (same logic as your protocol)
        key_disp = np.array([
             0.0,
             0.1*DMAX,   -0.1*DMAX,
             0.5*DMAX,   -0.5*DMAX,
             0.8*DMAX,   -0.8*DMAX,
             DMAX,       -DMAX,
             0.1*DMAX,   -0.1*DMAX,
             0.5*DMAX,   -0.5*DMAX,
             0.8*DMAX,   -0.8*DMAX,
             DMAX,       -DMAX,
             0.0
        ])
        
        # Generate 1000-point displacement protocol
        disp_protocol = np.interp(
            np.linspace(0, len(key_disp) - 1, n_points),
            np.arange(len(key_disp)),
            key_disp
        )
        
        ops.analysis('Static')
        STEP = 0.0
        for target_disp in disp_protocol:
            current_disp = ops.nodeDisp(center_node, 1) # DISPALCEMENT APPLIED IN MIDDLE NOD IN X DIR.
            dU = target_disp - current_disp
            ops.integrator('DisplacementControl', center_node, IDctrlDOF, dU)
            OK = ops.analyze(1)
            S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS)
            ops.reactions()
            reactionX.append(ops.nodeReaction(1, 1))           # BASE REACTION IN X DIR.
            reactionY.append(ops.nodeReaction(1, 2))           # BASE REACTION IN Y DIR.
            reactionZ.append(ops.nodeReaction(1, 3))           # BASE REACTION IN Z DIR.
            dispX.append(ops.nodeDisp(center_node, 1))         # DISPLACEMENT IN X DIR.  
            dispY.append(ops.nodeDisp(center_node, 2))         # DISPLACEMENT IN Y DIR.
            dispZ.append(ops.nodeDisp(center_node, 3))         # ROTATION IN Z DIR.
            stiffnessX.append(np.abs(reactionX[-1]) / np.abs(dispX[-1]))
            stiffnessY.append(np.abs(reactionY[-1]) / np.abs(dispY[-1]))
            stiffnessZ.append(np.abs(reactionZ[-1]) / np.abs(dispZ[-1]))      
            # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
            #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(1, 0.5*DR, DR, 0, 1)
            PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(1, PLOT=True)
            PERIOD_MIN.append(PERIODmin)
            PERIOD_MAX.append(PERIODmax)
            STEP += 1
            #print(STEP, disp[-1], reaction[-1])
            print(f"Step: {STEP}, Displacement: {dispX[-1]:.4f} mm, Reaction: {reactionX[-1]:.2f} N")     
        else:
            print('\n\nCYCLIC DISPLAEMENT ANALYSIS DONE.\n\n')
            
        DATA = (reactionX, dispX, velX, accelX,
                reactionY, dispY, velY, accelY,
                reactionZ, dispZ, velZ, accelZ,
                stiffnessX, stiffnessY, stiffnessZ,
                np.array(PERIOD_MIN), np.array(PERIOD_MAX))
    
        return  DATA 

    if ANAL_TYPE == 'STATIC_EXTERNAL_TIME-DEPENDENT_LOADING': # STATIC TIME-HISTORY ANALYSIS
        #%% DEFINE EXTERNAL TIME-DEPENDENT LOADING PROPERTIES
        # IN HERE ANALYSIS TIME AND DURATION ARE LOAD STEPS
        duration = 20.0             # [s] Analysis duration
        dt = 0.01                   # [s] Time step
        DT = dt                     # [s] Time step
        DT_time = 5.0               # [s] Total external Load Analysis Durations [*******]
        force_amplitude = 960.0     # [N] Amplitude Force
        omega_DT = 5.0715           # [rad/s] Natural angular frequency
        time_steps = int(duration/dt)
        
        # Check function
        def CHECK_FUN(DT_time ,duration):
            if DT_time > duration:
                print('\n\nAnalysis Duration Must be greater than External Load Duration\n\n')        
                exit()
            return -1

        CHECK_FUN(DT_time ,duration)
        def EXTERNAL_TIME_DEPENDENT(force_amplitude, omega_DT, DT, DT_time): # P(t) = P0 sin(wt)
            import numpy as np
            import matplotlib.pyplot as plt
            # External Load Durations
            num_steps = int(DT_time / DT)
            load_time = np.linspace(0, DT_time, num_steps) 
            target_frequency = 1.0 * omega_DT  # Target excitation frequency
            DT_load = force_amplitude * np.sin(target_frequency * load_time)
            # Plot External Time-dependent Loading
            plt.figure(figsize=(10, 6))
            plt.plot(load_time, DT_load, label=f'External Loading - Max: {np.max(DT_load):.3f}', linewidth=5)
            plt.title('External Time-dependent Loading Over Time')
            plt.xlabel('Time (s)')
            plt.ylabel('Force (N)')
            plt.grid(True)
            plt.legend()
            plt.show()
            return DT_load

        #DT_load = EXTERNAL_TIME_DEPENDENT(force_amplitude, omega_DT, DT, DT_time)

        def EXTERNAL_TIME_DEPENDENT_02(force_amplitude, omega_DT, DT, DT_time): #  P(t) = P0 exp(-0.05wt) sin(wt) 
            import numpy as np
            import matplotlib.pyplot as plt
            # External Load Durations
            num_steps = int(DT_time / DT)
            load_time = np.linspace(0, DT_time, num_steps) 
            target_frequency = 1.0 * omega_DT  # Target excitation frequency
            DT_load = force_amplitude * np.exp(-0.05*target_frequency * load_time) * np.sin(target_frequency * load_time)
            # Plot External Time-dependent Loading
            plt.figure(figsize=(10, 6))
            plt.plot(load_time, DT_load, label=f'External Loading - Max: {np.max(DT_load):.3f}', linewidth=5)
            plt.title('External Time-dependent Loading Over Time')
            plt.xlabel('Time (s)')
            plt.ylabel('Force (N)')
            plt.grid(True)
            plt.legend()
            plt.show()
            #print(load_time, DT_load)
            return DT_load

        DT_load02 = EXTERNAL_TIME_DEPENDENT_02(force_amplitude, omega_DT, DT, DT_time)
        # Static Time-depenent External loading analysis
        TS_TAG = 3
        PATT_TAG = 3
        # Apply time-dependent explosion loading
        ops.timeSeries('Path', TS_TAG, '-dt', dt, '-values', *DT_load02)
        ops.pattern('Plain', TS_TAG, PATT_TAG)
        ops.load(center_node, 1.0, 0.0, 0.0)
        
        ops.system('BandGeneral')
        ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/normDispIncr.html
        ops.algorithm('Newton')  # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/algorithm.html
        ops.integrator('LoadControl', dt)
        #ops.integrator('DisplacementControl', center_node, 0.001)
        ops.analysis('Static') # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/analysis.html
        
        STEP = 0.0
        stable = 0
        
        for JJ in range(time_steps):
            stable = ops.analyze(1)
            S02.ANALYSIS(stable, 1, MAX_TOLERANCE, MAX_ITERATIONS)
            ops.reactions()
            reactionX.append(ops.nodeReaction(1, 1))           # BASE REACTION IN X DIR.
            reactionY.append(ops.nodeReaction(1, 2))           # BASE REACTION IN Y DIR.
            reactionZ.append(ops.nodeReaction(1, 3))           # BASE REACTION IN Z DIR.
            dispX.append(ops.nodeDisp(center_node, 1))         # DISPLACEMENT IN X DIR.  
            dispY.append(ops.nodeDisp(center_node, 2))         # DISPLACEMENT IN Y DIR.
            dispZ.append(ops.nodeDisp(center_node, 3))         # ROTATION IN Z DIR.
            stiffnessX.append(np.abs(reactionX[-1]) / np.abs(dispX[-1]))
            stiffnessY.append(np.abs(reactionY[-1]) / np.abs(dispY[-1]))
            stiffnessZ.append(np.abs(reactionZ[-1]) / np.abs(dispZ[-1]))     
            # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
            #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(1, 0.5*DR, DR, 0, 1)
            PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(1, PLOT=True)
            PERIOD_MIN.append(PERIODmin)
            PERIOD_MAX.append(PERIODmax)
            STEP += 1
            #print(STEP, disp[-1], reaction[-1])
            print(f"Step: {STEP}, Displacement: {dispX[-1]:.4f} mm, Reaction: {reactionX[-1]:.2f} N")     
        else:
            print('\n\nSTATIC EXTERNAL TIME-DEPENDENT LOADING ANALYSIS DONE.\n\n')
            
        DATA = (reactionX, dispX, velX, accelX,
                reactionY, dispY, velY, accelY,
                reactionZ, dispZ, velZ, accelZ,
                stiffnessX, stiffnessY, stiffnessZ,
                np.array(PERIOD_MIN), np.array(PERIOD_MAX))
    
        return  DATA
    
    if ANAL_TYPE == 'DYNAMIC_EXTERNAL_TIME-DEPENDENT_LOADING': # DYNAMIC TIME-HISTORY ANALYSIS
        #%% DEFINE EXTERNAL TIME-DEPENDENT LOADING PROPERTIES
        duration = 20.0             # [s] Analysis duration
        dt = 0.01                   # [s] Time step
        DT = dt                     # [s] Time step
        DT_time = 5.0               # [s] Total external Load Analysis Durations [*******]
        force_amplitude = 960.0     # [N] Amplitude Force
        omega_DT = 5.0715           # [rad/s] Natural angular frequency
        DR = 0.03                   # Damping Ratio

        # Check function
        def CHECK_FUN(DT_time ,duration):
            if DT_time > duration:
                print('\n\nAnalysis Duration Must be greater than External Load Duration\n\n')        
                exit()
            return -1

        CHECK_FUN(DT_time ,duration)
        def EXTERNAL_TIME_DEPENDENT(force_amplitude, omega_DT, DT, DT_time): # P(t) = P0 sin(wt)
            import numpy as np
            import matplotlib.pyplot as plt
            # External Load Durations
            num_steps = int(DT_time / DT)
            load_time = np.linspace(0, DT_time, num_steps) 
            target_frequency = 1.0 * omega_DT  # Target excitation frequency
            DT_load = force_amplitude * np.sin(target_frequency * load_time)
            # Plot External Time-dependent Loading
            plt.figure(figsize=(10, 6))
            plt.plot(load_time, DT_load, label=f'External Loading - Max: {np.max(DT_load):.3f}', linewidth=5)
            plt.title('External Time-dependent Loading Over Time')
            plt.xlabel('Time (s)')
            plt.ylabel('Force (N)')
            plt.grid(True)
            plt.legend()
            plt.show()
            return DT_load

        #DT_load = EXTERNAL_TIME_DEPENDENT(force_amplitude, omega_DT, DT, DT_time)

        def EXTERNAL_TIME_DEPENDENT_02(force_amplitude, omega_DT, DT, DT_time): #  P(t) = P0 exp(-0.05wt) sin(wt) 
            import numpy as np
            import matplotlib.pyplot as plt
            # External Load Durations
            num_steps = int(DT_time / DT)
            load_time = np.linspace(0, DT_time, num_steps) 
            target_frequency = 1.0 * omega_DT  # Target excitation frequency
            DT_load = force_amplitude * np.exp(-0.05*target_frequency * load_time) * np.sin(target_frequency * load_time)
            # Plot External Time-dependent Loading
            plt.figure(figsize=(10, 6))
            plt.plot(load_time, DT_load, label=f'External Loading - Max: {np.max(DT_load):.3f}', linewidth=5)
            plt.title('External Time-dependent Loading Over Time')
            plt.xlabel('Time (s)')
            plt.ylabel('Force (N)')
            plt.grid(True)
            plt.legend()
            plt.show()
            #print(load_time, DT_load)
            return DT_load

        DT_load02 = EXTERNAL_TIME_DEPENDENT_02(force_amplitude, omega_DT, DT, DT_time)
        # Static Time-depenent External loading analysis
        TS_TAG = 3
        PATT_TAG = 3
        # Apply time-dependent explosion loading
        ops.timeSeries('Path', TS_TAG, '-dt', dt, '-values', *DT_load02)
        ops.pattern('Plain', TS_TAG, PATT_TAG)
        ops.load(center_node, 1.0, 0.0, 0.0)
        
        ops.constraints('Plain')
        ops.numberer('Plain')
        ops.system('BandGeneral')
        ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/normDispIncr.html
        #ops.integrator('CentralDifference')  # JUST FOR LINEAR ANALYSIS - INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/centralDifference.html
        alpha=0.5; beta=0.25;
        ops.integrator('Newmark', alpha, beta) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/newmark.html
        #alpha=2/3;gamma=1.5-alpha; gamma=1.5-alpha;beta=(2-alpha)**2/4;
        #ops.integrator('HHT', alpha, gamma, beta) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/hht.html
        ops.algorithm('Newton')  # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/algorithm.html
        ops.analysis('Transient') # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/analysis.html
        
        stable = 0
        current_time = 0.0
        
        while stable == 0 and current_time < duration:
            ops.analyze(1, dt)
            S02.ANALYSIS(stable, 1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS
            current_time = ops.getTime()
            time.append(current_time)
            ops.reactions()
            reactionX.append(ops.nodeReaction(1, 1))           # BASE REACTION IN X DIR.
            reactionY.append(ops.nodeReaction(1, 2))           # BASE REACTION IN Y DIR.
            reactionZ.append(ops.nodeReaction(1, 3))           # BASE REACTION IN Z DIR.
            dispX.append(ops.nodeDisp(center_node, 1))         # DISPLACEMENT IN X DIR.  
            dispY.append(ops.nodeDisp(center_node, 2))         # DISPLACEMENT IN Y DIR.
            dispZ.append(ops.nodeDisp(center_node, 3))         # ROTATION IN Z DIR.
            velX.append(ops.nodeVel(center_node, 1))           # VELOCITY IN X DIR.
            velY.append(ops.nodeVel(center_node, 2))           # VELOCITY IN Y DIR.
            velZ.append(ops.nodeVel(center_node, 3))           # VELOCITY IN Z DIR.
            accelX.append(ops.nodeAccel(center_node, 1))       # ACCELERATION IN X DIR.
            accelY.append(ops.nodeAccel(center_node, 2))       # ACCELERATION IN Y DIR.
            accelZ.append(ops.nodeAccel(center_node, 3))       # ACCELERATION IN Z DIR.
            stiffnessX.append(np.abs(reactionX[-1]) / np.abs(dispX[-1]))
            stiffnessY.append(np.abs(reactionY[-1]) / np.abs(dispY[-1]))
            stiffnessZ.append(np.abs(reactionZ[-1]) / np.abs(dispZ[-1]))
            OMEGA.append(np.sqrt(stiffnessZ[-1]/TOTAL_MASS))
            PERIOD.append((np.pi * 2) / OMEGA[-1])   
            # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
            #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(1, 0.5*DR, DR, 0, 1)
            PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(1, PLOT=True)
            PERIOD_MIN.append(PERIODmin)
            PERIOD_MAX.append(PERIODmax)
            #print(time[-1], disp[-1], veloY[-1])
            print(f"Time: {time[-1]:.4f}, Displacement: {dispX[-1]:.4f} mm, Reaction: {reactionX[-1]:.2f} N")      
        else:
            print('\n\nDYNAMIC EXTERNAL TIME-DEPENDENT LOADING ANALYSIS DONE.\n\n')  
        # Calculating Damping Ratio and Period Using Logarithmic Decrement Analysis 
        damping_ratio = S04.DAMPING_RATIO(dispX)  
        
        # Compute modal properties
        ops.modalProperties("-print", "-file", f"SALAR_ModalReport_{ANAL_TYPE}.txt", "-unorm") 
        
        DATA = (time, reactionX, dispX, velX, accelX,
                reactionY, dispY, velY, accelY,
                reactionZ, dispZ, velZ, accelZ,
                stiffnessX, stiffnessY, stiffnessZ,
                PERIOD, damping_ratio,
                np.array(PERIOD_MIN), np.array(PERIOD_MAX))
        
        return  DATA
        
    if ANAL_TYPE == 'FREE-VIBRATION': # DYNAMIC TIME-HISTORY ANALYSIS
        #%% DEFINE PARAMETERS FOR FREE-VIBRATION ANALYSIS
        u0 = -0.10                         # [m] Initial displacement
        v0 = 0.15                          # [m/s] Initial velocity
        a0 = 0.065                         # [m/s^2] Initial acceleration
        IU = True                          # Free Vibration with Initial Displacement
        IV = True                          # Free Vibration with Initial Velocity
        IA = True                          # Free Vibration with Initial Acceleration
        duration = 20.0                    # [s] Analysis duration
        dt = 0.001                         # [s] Time step
        DR = 0.03                          # Damping Ratio
        
        if IU == True:
            # Define initial displacment
            ops.setNodeDisp(center_node, 1, u0, '-commit') # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/setNodeDisp.html
        if IV == True:
            # Define initial velocity
            ops.setNodeVel(center_node, 1, v0, '-commit')  # INFO LINK: https://openseespydoc.readthedocs.io/en/stable/src/setNodeVel.html
        if IA == True:
            # Define initial  acceleration
            ops.setNodeAccel(center_node, 1, a0, '-commit') # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/setNodeAccel.html
            
        ops.constraints('Plain')
        ops.numberer('Plain')
        ops.system('BandGeneral')
        ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/normDispIncr.html
        #ops.integrator('CentralDifference')  # JUST FOR LINEAR ANALYSIS - INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/centralDifference.html
        alpha=0.5; beta=0.25;
        ops.integrator('Newmark', alpha, beta) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/newmark.html
        #alpha=2/3;gamma=1.5-alpha; gamma=1.5-alpha;beta=(2-alpha)**2/4;
        #ops.integrator('HHT', alpha, gamma, beta) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/hht.html
        ops.algorithm('Newton')  # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/algorithm.html
        ops.analysis('Transient') # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/analysis.html
        
        stable = 0
        current_time = 0.0
        while stable == 0 and current_time < duration:
            ops.analyze(1, dt)
            S02.ANALYSIS(stable, 1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS
            current_time = ops.getTime()
            time.append(current_time)
            ops.reactions()
            reactionX.append(ops.nodeReaction(1, 1))           # BASE REACTION IN X DIR.
            reactionY.append(ops.nodeReaction(1, 2))           # BASE REACTION IN Y DIR.
            reactionZ.append(ops.nodeReaction(1, 3))           # BASE REACTION IN Z DIR.
            dispX.append(ops.nodeDisp(center_node, 1))         # DISPLACEMENT IN X DIR.  
            dispY.append(ops.nodeDisp(center_node, 2))         # DISPLACEMENT IN Y DIR.
            dispZ.append(ops.nodeDisp(center_node, 3))         # ROTATION IN Z DIR.
            velX.append(ops.nodeVel(center_node, 1))           # VELOCITY IN X DIR.
            velY.append(ops.nodeVel(center_node, 2))           # VELOCITY IN Y DIR.
            velZ.append(ops.nodeVel(center_node, 3))           # VELOCITY IN Z DIR.
            accelX.append(ops.nodeAccel(center_node, 1))       # ACCELERATION IN X DIR.
            accelY.append(ops.nodeAccel(center_node, 2))       # ACCELERATION IN Y DIR.
            accelZ.append(ops.nodeAccel(center_node, 3))       # ACCELERATION IN Z DIR.
            stiffnessX.append(np.abs(reactionX[-1]) / np.abs(dispX[-1]))
            stiffnessY.append(np.abs(reactionY[-1]) / np.abs(dispY[-1]))
            stiffnessZ.append(np.abs(reactionZ[-1]) / np.abs(dispZ[-1]))
            OMEGA.append(np.sqrt(stiffnessZ[-1]/TOTAL_MASS))
            PERIOD.append((np.pi * 2) / OMEGA[-1]) 
            # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
            #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(1, 0.5*DR, DR, 0, 1)
            PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(1, PLOT=True)
            PERIOD_MIN.append(PERIODmin)
            PERIOD_MAX.append(PERIODmax)
            #print(time[-1], disp[-1], veloY[-1])
            print(f"Time: {time[-1]:.4f}, Displacement: {dispX[-1]:.4f} mm, Reaction: {reactionX[-1]:.2f} N")      
        else:
            print('\n\nFREE-VIBRATION ANALYSIS DONE.\n\n')    
        # Calculating Damping Ratio and Period Using Logarithmic Decrement Analysis 
        damping_ratio = S04.DAMPING_RATIO(dispX)  
        
        # Compute modal properties
        ops.modalProperties("-print", "-file", f"SALAR_ModalReport_{ANAL_TYPE}.txt", "-unorm") 
        
        DATA = (time, reactionX, dispX, velX, accelX,
                reactionY, dispY, velY, accelY,
                reactionZ, dispZ, velZ, accelZ,
                stiffnessX, stiffnessY, stiffnessZ,
                PERIOD, damping_ratio,
                np.array(PERIOD_MIN), np.array(PERIOD_MAX))
        
        return DATA
    if ANAL_TYPE == 'SEISMIC': # DYNAMIC TIME-HISTORY ANALYSIS
        #%% DEFINE PARAMETERS FOR SEISMIC ANALYSIS
        duration = 20.0                    # [s] Analysis duration
        dt = 0.01                          # [s] Time step
        GMfact = 9810                      # [mm/s²]standard acceleration of gravity or standard acceleration
        SSF_X = 0.10                       # Seismic Acceleration Scale Factor in X Direction
        SSF_Y = 0.10                       # Seismic Acceleration Scale Factor in Y Direction
        iv0_X = 0.0005                     # [mm/s] Initial velocity applied to the node  in X Direction
        iv0_Y = 0.0005                     # [mm/s] Initial velocity applied to the node  in Y Direction
        st_iv0 = 0.0                       # [s] Initial velocity applied starting time
        SEI = 'X'                          # Seismic Direction
        DR = 0.03                          # Damping ratio
        
        # Define time series for input motion (Acceleration time history)
        if SEI == 'X':
            SEISMIC_TAG_01 = 100
            ops.timeSeries('Path', SEISMIC_TAG_01, '-dt', dt, '-filePath', 'Ground_Acceleration_X.txt', '-factor', GMfact, '-startTime', st_iv0) # SEISMIC-X
            # Define load patterns
            # pattern UniformExcitation $patternTag $dof -accel $tsTag <-vel0 $vel0> <-fact $cFact>
            ops.pattern('UniformExcitation', SEISMIC_TAG_01, 1, '-accel', SEISMIC_TAG_01, '-vel0', iv0_X, '-fact', SSF_X) # SEISMIC-X
        if SEI == 'Y':
            SEISMIC_TAG_02 = 200
            ops.timeSeries('Path', SEISMIC_TAG_02, '-dt', dt, '-filePath', 'Ground_Acceleration_Y.txt', '-factor', GMfact) # SEISMIC-Z
            ops.pattern('UniformExcitation', SEISMIC_TAG_02, 2, '-accel', SEISMIC_TAG_02, '-vel0', iv0_Y, '-fact', SSF_Y) 
        if SEI == 'XY':
            SEISMIC_TAG_01 = 100
            ops.timeSeries('Path', SEISMIC_TAG_01, '-dt', dt, '-filePath', 'Ground_Acceleration_X.txt', '-factor', GMfact, '-startTime', st_iv0) # SEISMIC-X
            # Define load patterns
            # pattern UniformExcitation $patternTag $dof -accel $tsTag <-vel0 $vel0> <-fact $cFact>
            ops.pattern('UniformExcitation', SEISMIC_TAG_01, 1, '-accel', SEISMIC_TAG_01, '-vel0', iv0_X, '-fact', SSF_X) # SEISMIC-X 
            SEISMIC_TAG_02 = 200
            ops.timeSeries('Path', SEISMIC_TAG_02, '-dt', dt, '-filePath', 'Ground_Acceleration_Y.txt', '-factor', GMfact) # SEISMIC-Z
            ops.pattern('UniformExcitation', SEISMIC_TAG_02, 2, '-accel', SEISMIC_TAG_02, '-vel0', iv0_Y, '-fact', SSF_Y)  # SEISMIC-Z
        print('Seismic Defined Done.')
        
        ops.constraints('Plain')
        ops.numberer('Plain')
        ops.system('BandGeneral')
        ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/normDispIncr.html
        #ops.integrator('CentralDifference')  # JUST FOR LINEAR ANALYSIS - INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/centralDifference.html
        alpha=0.5; beta=0.25;
        ops.integrator('Newmark', alpha, beta) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/newmark.html
        #alpha=2/3;gamma=1.5-alpha; gamma=1.5-alpha;beta=(2-alpha)**2/4;
        #ops.integrator('HHT', alpha, gamma, beta) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/hht.html
        ops.algorithm('Newton')  # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/algorithm.html
        ops.analysis('Transient') # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/analysis.html
        
        stable = 0
        current_time = 0.0
        while stable == 0 and current_time < duration:
            ops.analyze(1, dt)
            S02.ANALYSIS(stable, 1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS
            current_time = ops.getTime()
            time.append(current_time)
            ops.reactions()
            reactionX.append(ops.nodeReaction(1, 1))           # BASE REACTION IN X DIR.
            reactionY.append(ops.nodeReaction(1, 2))           # BASE REACTION IN Y DIR.
            reactionZ.append(ops.nodeReaction(1, 3))           # BASE REACTION IN Z DIR.
            dispX.append(ops.nodeDisp(center_node, 1))         # DISPLACEMENT IN X DIR.  
            dispY.append(ops.nodeDisp(center_node, 2))         # DISPLACEMENT IN Y DIR.
            dispZ.append(ops.nodeDisp(center_node, 3))         # ROTATION IN Z DIR.
            velX.append(ops.nodeVel(center_node, 1))           # VELOCITY IN X DIR.
            velY.append(ops.nodeVel(center_node, 2))           # VELOCITY IN Y DIR.
            velZ.append(ops.nodeVel(center_node, 3))           # VELOCITY IN Z DIR.
            accelX.append(ops.nodeAccel(center_node, 1))       # ACCELERATION IN X DIR.
            accelY.append(ops.nodeAccel(center_node, 2))       # ACCELERATION IN Y DIR.
            accelZ.append(ops.nodeAccel(center_node, 3))       # ACCELERATION IN Z DIR.
            stiffnessX.append(np.abs(reactionX[-1]) / np.abs(dispX[-1]))
            stiffnessY.append(np.abs(reactionY[-1]) / np.abs(dispY[-1]))
            stiffnessZ.append(np.abs(reactionZ[-1]) / np.abs(dispZ[-1]))
            OMEGA.append(np.sqrt(stiffnessZ[-1]/TOTAL_MASS))
            PERIOD.append((np.pi * 2) / OMEGA[-1])
            # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
            #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(1, 0.5*DR, DR, 0, 1)
            PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(1, PLOT=True)
            PERIOD_MIN.append(PERIODmin)
            PERIOD_MAX.append(PERIODmax)
            #print(time[-1], disp[-1], veloY[-1])
            print(f"Time: {time[-1]:.4f}, Displacement: {dispX[-1]:.4f} mm, Reaction: {reactionX[-1]:.2f} N")

        else:
            print('\n\nSEISMIC ANALYSIS DONE.\n\n')     
        # Calculating Damping Ratio and Period Using Logarithmic Decrement Analysis 
        damping_ratio = S04.DAMPING_RATIO(dispX)  
        
        # Compute modal properties
        ops.modalProperties("-print", "-file", f"SALAR_ModalReport_{ANAL_TYPE}.txt", "-unorm") 
        
        DATA = (time, reactionX, dispX, velX, accelX,
                reactionY, dispY, velY, accelY,
                reactionZ, dispZ, velZ, accelZ,
                stiffnessX, stiffnessY, stiffnessZ,
                PERIOD, damping_ratio,
                np.array(PERIOD_MIN), np.array(PERIOD_MAX))
        
        return DATA    

#%%-------------------------------------------------------
def PLOT_TIME_HISTORY(time, reaction, disp, velo, acc):

    import matplotlib.pyplot as plt

    data_dict = {
        "Reaction [N]": reaction,
        "Displacement [m]": disp,
        "Velocity [m/s]": velo,
        "Acceleration [m/s²]": acc,
    }

    fig, axes = plt.subplots(4, 2, figsize=(15, 12))
    axes = axes.flatten()

    for i, (label, data) in enumerate(data_dict.items()):
        axes[i].plot(time, data, color='black', linewidth=2)
        axes[i].set_title(label, fontsize=9)
        axes[i].set_xlabel("Time (s)")
        axes[i].grid(True)

  
    for j in range(len(data_dict), len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout()
    plt.show()


def PLOT(XDATA, YDATA, TITLE, XLABEL, YLABEL, COLOR, SEMILOGY):
    import matplotlib.pyplot as plt
    fig = plt.figure(-1, figsize=(12, 8))
    plt.plot(XDATA, YDATA, color=COLOR, linewidth=2)
    plt.title(TITLE)
    plt.ylabel(YLABEL)
    plt.xlabel(XLABEL)
    if SEMILOGY == True:
        plt.semilogy()
    plt.grid()
    
# Plotting Nodal Displacements
def PLOT_DISPLAEMENTS(time_steps, displacements_dict, TITLE):
    plt.figure(figsize=(10, 6))
    
    for node_id, disp_values in displacements_dict.items():
        plt.plot(time_steps, disp_values, label=f'Node {node_id} - MAX. ABS. : {np.max(np.abs(disp_values)): 0.4e}', linewidth=2)
    
    plt.xlabel('Time [s]')
    plt.ylabel('Displacement [m]')
    plt.title(TITLE)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()  
    
# Plotting Element Forces
def PLOT_FORCES(time_steps, forces_dict, YLABEL, TITLE):
    plt.figure(figsize=(10, 6))
    
    for node_id, force_values in forces_dict.items():
        plt.plot(time_steps, force_values, label=f'Ele. {node_id} - MAX. ABS. : {np.max(np.abs(force_values)): 0.4e}', linewidth=2)
    
    plt.xlabel('Time [s]')
    plt.ylabel(YLABEL)
    plt.title(TITLE)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()  
    
def PLOT_2D(X, Y, Xfit, Yfit, X2, Y2, XLABEL, YLABEL, TITLE, LEGEND01, LEGEND02, LEGEND03, COLOR, Z):
    import matplotlib.pyplot as plt
    plt.figure(figsize=(12, 8))
    if Z == 1:
        # Plot 1 line
        plt.plot(X, Y,color=COLOR)
        plt.xlabel(XLABEL)
        plt.ylabel(YLABEL)
        plt.title(TITLE)
        plt.grid(True)
        plt.show()
    if Z == 2:
        # Plot 2 lines
        plt.plot(X, Y, Xfit, Yfit, 'r--', linewidth=3)
        plt.title(TITLE)
        plt.xlabel(XLABEL)
        plt.ylabel(YLABEL)
        plt.legend([LEGEND01, LEGEND02], loc='lower right')
        plt.grid(True)
        plt.show()
    if Z == 3:
        # Plot 3 lines
        plt.plot(X, Y, Xfit, Yfit, 'r--', X2, Y2, 'g-*', linewidth=3)
        plt.title(TITLE)
        plt.xlabel(XLABEL)
        plt.ylabel(YLABEL)
        plt.legend([LEGEND01, LEGEND02, LEGEND03], loc='lower right')
        plt.grid(True)
        plt.show() 
#%%----------------------------------------------------
TOTAL_MASS = 5000.0   # [kg] Total Mass of Structure
#%%----------------------------------------------------
# PERIOD ANALYSIS
SPRING_TYPE = 'ELASTIC' # Kθ(M-θ) - 'ELASTIC' OR 'INELASTIC' 
MAT_TYPE = 'ELASTIC'    # K∆(P-∆) 'ELASTIC' OR 'INELASTIC'
ANAL_TYPE = 'PERIOD'

DATA = PENDULUM(SPRING_TYPE, MAT_TYPE, TOTAL_MASS)
(PERIOD_MIN_X, PERIOD_MAX_X) = DATA
print('Structure First Period:  ', PERIOD_MIN_X)
print('Structure Second Period: ', PERIOD_MAX_X) 


#%%----------------------------------------------------
# STATIC ANALYSIS
SPRING_TYPE = 'ELASTIC' # Kθ(M-θ) - 'ELASTIC' OR 'INELASTIC' 
MAT_TYPE = 'ELASTIC'    # K∆(P-∆) 'ELASTIC' OR 'INELASTIC'
ANAL_TYPE = 'STATIC'

DATA = PENDULUM(SPRING_TYPE, MAT_TYPE, TOTAL_MASS)
(reactionX, dispX, velX, accelX,
        reactionY, dispY, velY, accelY,
        reactionZ, dispZ, velZ, accelZ,
        stiffnessX, stiffnessY, stiffnessZ) = DATA

S01.PLOT_2D_FRAME(deformed_scale=1)  # Adjust scale factor as needed
#%%----------------------------------------------------
# PUSHOVER ANALYSIS (STATIC TIME-HISTORY ANALYSIS)
SPRING_TYPE = 'ELASTIC' # Kθ(M-θ) - 'ELASTIC' OR 'INELASTIC' 
MAT_TYPE = 'INELASTIC'    # K∆(P-∆) 'ELASTIC' OR 'INELASTIC'
ANAL_TYPE = 'PUSHOVER'

DATA = PENDULUM(SPRING_TYPE, MAT_TYPE, TOTAL_MASS)
(reactionX_PUSH, dispX_PUSH, velX_PUSH, accelX_PUSH,
        reactionY_PUSH, dispY_PUSH, velY_PUSH, accelY_PUSH,
        reactionZ_PUSH, dispZ_PUSH, velZ_PUSH, accelZ_PUSH,
        stiffnessX_PUSH, stiffnessY_PUSH, stiffnessZ_PUSH,
 PERIOD_MIN_PUSH, PERIOD_MAX_PUSH) = DATA


XDATA = dispX_PUSH
YDATA = reactionX_PUSH
XLABEL = 'Displacement [m]'
YLABEL = 'Base Reaction [N]'
TITLE = 'Base Reaction and Displacement of Structure During Pushover Analysis'
COLOR = 'black'
SEMILOGY = False
PLOT(XDATA, YDATA, TITLE, XLABEL, YLABEL, COLOR, SEMILOGY)

DATA = S07.BILNEAR_CURVE(np.abs(dispX_PUSH), np.abs(reactionX_PUSH), SLOPE_NODE=10)
(X_PUSH, Y_PUSH, Elastic_ST, Plastic_ST, Tangent_ST, Ductility_Rito, Over_Strength_Factor) = DATA

# PLOT STRUCTURAL PERIOD DURING THE ANALYSIS
plt.figure(0, figsize=(12, 8))
plt.plot(dispX_PUSH, PERIOD_MIN_PUSH, linewidth=3)
plt.plot(dispX_PUSH, PERIOD_MAX_PUSH, linewidth=3)
plt.title('Period of Structure During Pushover Analysis')
plt.ylabel('Structural Period [s]')
plt.xlabel('Displacement [m]')
#plt.semilogy()
plt.grid()
plt.legend([f'PERIOD - MIN VALUES: Min: {np.min(PERIOD_MIN_PUSH):.3f} (s) - Mean: {np.mean(PERIOD_MIN_PUSH):.3f} (s) - Max: {np.max(PERIOD_MIN_PUSH):.3f} (s)', 
            f'PERIOD - MAX VALUES:  Min: {np.min(PERIOD_MAX_PUSH):.3f} (s) - Mean: {np.mean(PERIOD_MAX_PUSH):.3f} (s) - Max: {np.max(PERIOD_MAX_PUSH):.3f} (s)',
            ])
plt.show()

S01.PLOT_2D_FRAME(deformed_scale=1)  # Adjust scale factor as needed
#%%----------------------------------------------------
# CYCLIC DISPLACEMENT ANALYSIS (STATIC TIME-HISTORY ANALYSIS)
MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'
ANAL_TYPE = 'CYCLIC_DISPLACEMENT'

DATA = PENDULUM(SPRING_TYPE, MAT_TYPE, TOTAL_MASS)
(reactionX_CP, dispX_CP, velX_CP_CP, accelX_CP,
        reactionY_CP, dispY_CP, velY_CP, accelY_CP,
        reactionZ_CP, dispZ_CP, velZ_CP, accelZ_CP,
        stiffnessX_CP, stiffnessY_CP, stiffnessZ_CP,
 PERIOD_MIN_CP, PERIOD_MAX_CP) = DATA


XDATA = dispX_CP
YDATA = reactionX_CP
XLABEL = 'Displacement [m]'
YLABEL = 'Base Reaction [N]'
TITLE = 'Base Reaction and Displacement of Structure During Cyclic-Displacement Analysis'
COLOR = 'black'
SEMILOGY = False
PLOT(XDATA, YDATA, TITLE, XLABEL, YLABEL, COLOR, SEMILOGY)

# PLOT STRUCTURAL PERIOD DURING THE ANALYSIS
plt.figure(0, figsize=(12, 8))
plt.plot(dispX_CP, PERIOD_MIN_CP, linewidth=3)
plt.plot(dispX_CP, PERIOD_MAX_CP, linewidth=3)
plt.title('Period of Structure During Cyclic Displacement Analysis')
plt.ylabel('Structural Period [s]')
plt.xlabel('Displacement [m]')
#plt.semilogy()
plt.grid()
plt.legend([f'PERIOD - MIN VALUES: Min: {np.min(PERIOD_MIN_CP):.3f} (s) - Mean: {np.mean(PERIOD_MIN_CP):.3f} (s) - Max: {np.max(PERIOD_MIN_CP):.3f} (s)', 
            f'PERIOD - MAX VALUES:  Min: {np.min(PERIOD_MAX_CP):.3f} (s) - Mean: {np.mean(PERIOD_MAX_CP):.3f} (s) - Max: {np.max(PERIOD_MAX_CP):.3f} (s)',
            ])
plt.show()

S01.PLOT_2D_FRAME(deformed_scale=1)  # Adjust scale factor as needed
#%%----------------------------------------------------
# EXTERNAL TIME-DEPENDENT LOADING ANALYSIS (STATIC TIME-HISTORY ANALYSIS)
MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'
ANAL_TYPE = 'STATIC_EXTERNAL_TIME-DEPENDENT_LOADING'

DATA = PENDULUM(SPRING_TYPE, MAT_TYPE, TOTAL_MASS)

(reactionX_ETDLS, dispX_ETDLS, velX_ETDLS, accelX_ETDLS,
        reactionY_ETDLS, dispY_ETDLS, velY_ETDLS, accelY_ETDLS,
        reactionZ_ETDLS, dispZ_ETDLS, velZ_ETDLS, accelZ_ETDLS,
        stiffnessX_ETDLS, stiffnessY_ETDLS, stiffnessZ_ETDLS,
 PERIOD_MIN_ETDLS, PERIOD_MAX_ETDLS) = DATA


XDATA = dispX_ETDLS
YDATA = reactionX_ETDLS
XLABEL = 'Displacement [m]'
YLABEL = 'Base Reaction [N]'
TITLE = 'Base Reaction and Displacement of Structure During Static External Time-dependent Loading Analysis'
COLOR = 'black'
SEMILOGY = False
PLOT(XDATA, YDATA, TITLE, XLABEL, YLABEL, COLOR, SEMILOGY)

# PLOT STRUCTURAL PERIOD DURING THE ANALYSIS
plt.figure(0, figsize=(12, 8))
plt.plot(dispX_ETDLS, PERIOD_MIN_ETDLS, linewidth=3)
plt.plot(dispX_ETDLS, PERIOD_MAX_ETDLS, linewidth=3)
plt.title('Period of Structure During Static External Time-dependent Loading Analysis')
plt.ylabel('Structural Period [s]')
plt.xlabel('Displacement [m]')
#plt.semilogy()
plt.grid()
plt.legend([f'PERIOD - MIN VALUES: Min: {np.min(PERIOD_MIN_ETDLS):.3f} (s) - Mean: {np.mean(PERIOD_MIN_ETDLS):.3f} (s) - Max: {np.max(PERIOD_MIN_ETDLS):.3f} (s)', 
            f'PERIOD - MAX VALUES:  Min: {np.min(PERIOD_MAX_ETDLS):.3f} (s) - Mean: {np.mean(PERIOD_MAX_ETDLS):.3f} (s) - Max: {np.max(PERIOD_MAX_ETDLS):.3f} (s)',
            ])
plt.show()

S01.PLOT_2D_FRAME(deformed_scale=1)  # Adjust scale factor as needed
#%%----------------------------------------------------
# EXTERNAL TIME-DEPENDENT LOADING ANALYSIS (DYNAMIC TIME-HISTORY ANALYSIS)
SPRING_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'
MAT_TYPE = 'INELASTIC'      # 'ELASTIC' OR 'INELASTIC'
ANAL_TYPE = 'DYNAMIC_EXTERNAL_TIME-DEPENDENT_LOADING'

DATA = PENDULUM(SPRING_TYPE, MAT_TYPE, TOTAL_MASS)

(time_ETDLD, reactionX_ETDLD, dispX_ETDLD, velX_ETDLD, accelX_ETDLD,
        reactionY_ETDLD, dispY_ETDLD, velY_ETDLD, accelY_ETDLD,
        reactionZ_ETDLD, dispZ_ETDLD, velZ_ETDLD, accelZ_ETDLD,
        stiffnessX_ETDLD, stiffnessY_ETDLD, stiffnessZ_ETDLD,
PERIOD_ETDLD, damping_ratio_ETDLD,
PERIOD_MIN_ETDLD, PERIOD_MAX_ETDLD) = DATA



XDATA = dispX_ETDLD
YDATA = reactionX_ETDLD
XLABEL = 'Displacement [m]'
YLABEL = 'Base Reaction [N]'
TITLE = 'Base Reaction and Displacement of Structure During Dynamic External Time-dependent Loading Analysis'
COLOR = 'black'
SEMILOGY = False
PLOT(XDATA, YDATA, TITLE, XLABEL, YLABEL, COLOR, SEMILOGY)

PLOT_TIME_HISTORY(time_ETDLD, reactionX_ETDLD, dispX_ETDLD, velX_ETDLD, accelX_ETDLD)

# PLOT STRUCTURAL PERIOD DURING THE ANALYSIS
plt.figure(0, figsize=(12, 8))
plt.plot(dispX_ETDLD, PERIOD_MIN_ETDLD, linewidth=3)
plt.plot(dispX_ETDLD, PERIOD_MAX_ETDLD, linewidth=3)
plt.title('Period of Structure During Dynamic External Time-dependent Loading Analysis')
plt.ylabel('Structural Period [s]')
plt.xlabel('Displacement [m]')
#plt.semilogy()
plt.grid()
plt.legend([f'PERIOD - MIN VALUES: Min: {np.min(PERIOD_MIN_ETDLD):.3f} (s) - Mean: {np.mean(PERIOD_MIN_ETDLD):.3f} (s) - Max: {np.max(PERIOD_MIN_ETDLD):.3f} (s)', 
            f'PERIOD - MAX VALUES:  Min: {np.min(PERIOD_MAX_ETDLD):.3f} (s) - Mean: {np.mean(PERIOD_MAX_ETDLD):.3f} (s) - Max: {np.max(PERIOD_MAX_ETDLD):.3f} (s)',
            ])
plt.show()

S01.PLOT_2D_FRAME(deformed_scale=1)  # Adjust scale factor as needed
#%%----------------------------------------------------
# FREE-VIBRATION ANALYSIS (DYNAMIC TIME-HISTORY ANALYSIS)
MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'
ANAL_TYPE = 'FREE-VIBRATION'
DATA = PENDULUM(SPRING_TYPE, MAT_TYPE, TOTAL_MASS)

(time_FV, reactionX_FV, dispX_FV, velX_FV, accelX_FV,
        reactionY_FV, dispY_FV, velY_FV, accelY_FV,
        reactionZ_FV, dispZ_FV, velZ_FV, accelZ_FV,
        stiffnessX_FV, stiffnessY_FV, stiffnessZ_FV,
PERIOD_FV, damping_ratio_FV,
PERIOD_MIN_FV, PERIOD_MAX_FV) = DATA

XDATA = dispX_FV
YDATA = reactionX_FV
XLABEL = 'Displacement [m]'
YLABEL = 'Base Reaction [N]'
TITLE = 'Base Reaction and Dispalcement of Structure During Free-vibration Analysis'
COLOR = 'black'
SEMILOGY = False
PLOT(XDATA, YDATA, TITLE, XLABEL, YLABEL, COLOR, SEMILOGY)

PLOT_TIME_HISTORY(time_FV, reactionX_FV, dispX_FV, velX_FV, accelX_FV)

# PLOT STRUCTURAL PERIOD DURING THE ANALYSIS
plt.figure(0, figsize=(12, 8))
plt.plot(time_FV, PERIOD_MIN_FV, linewidth=3)
plt.plot(time_FV, PERIOD_MAX_FV, linewidth=3)
plt.title('Period of Structure During Free-vibration Analysis')
plt.ylabel('Structural Period [s]')
plt.xlabel('Time [s]')
#plt.semilogy()
plt.grid()
plt.legend([f'PERIOD - MIN VALUES: Min: {np.min(PERIOD_MIN_FV):.3f} (s) - Mean: {np.mean(PERIOD_MIN_FV):.3f} (s) - Max: {np.max(PERIOD_MIN_FV):.3f} (s)', 
            f'PERIOD - MAX VALUES:  Min: {np.min(PERIOD_MAX_FV):.3f} (s) - Mean: {np.mean(PERIOD_MAX_FV):.3f} (s) - Max: {np.max(PERIOD_MAX_FV):.3f} (s)',
            ])
plt.show()

S01.PLOT_2D_FRAME(deformed_scale=1)  # Adjust scale factor as needed

#%% EQULIVALENT VISCOUS DAMPING RATIO
def EQULIVALENT_VISCOUS_DAMPING_RATIO_FUN(displacement, base_shear, title="Hysteresis Curve"):
    """
    Compute the equivalent viscous damping ratio from a hysteresis loop.

    The dissipated energy per cycle (E_d) is the area enclosed by the loop,
    computed using the Shoelace formula. The maximum strain energy (E_s) is
    estimated as 0.5 * F_max * u_max, where u_max is the maximum absolute
    displacement and F_max is the absolute base shear at that same displacement.

    Parameters
    ----------
    displacement : array-like
        Displacement values (assumed to form a closed hysteresis loop).
    base_shear : array-like
        Corresponding base shear values.
    title : str, optional
        Title for the plot.

    Returns
    -------
    zeta_eq : float
        Equivalent viscous damping ratio.
    fig : matplotlib.figure.Figure
        Figure object with the plotted hysteresis loop and energy representation.
    """
    import numpy as np
    import matplotlib.pyplot as plt

    # Data preparation
    disp = np.asarray(displacement, dtype=float)
    shear = np.asarray(base_shear, dtype=float)

    if disp.size != shear.size:
        raise ValueError("Displacement and base shear arrays must have the same length.")
    if disp.size < 3:
        raise ValueError("At least 3 points are required to form a closed loop.")

    # Close the loop if not already closed (important for shoelace)
    if not (disp[0] == disp[-1] and shear[0] == shear[-1]):
        disp = np.append(disp, disp[0])
        shear = np.append(shear, shear[0])

    # Dissipated energy (E_d) via Shoelace formula ---
    x = disp
    y = shear
    E_d = 0.5 * np.abs(np.dot(x[:-1], y[1:]) - np.dot(y[:-1], x[1:]))

    # Maximum strain energy (E_s)
    # Find the point of maximum absolute displacement
    idx_max = np.argmax(np.abs(disp))
    u_max = np.abs(disp[idx_max])
    F_at_max_u = np.abs(shear[idx_max])

    E_s = 0.5 * F_at_max_u * u_max

    # Equivalent viscous damping ratio
    zeta_eq = 100 * E_d / (4.0 * np.pi * E_s) if E_s != 0 else 0.0

    # Plotting
    fig, ax = plt.subplots(figsize=(7, 6))

    # Hysteresis curve
    ax.plot(disp, shear, 'k-', linewidth=1.2, label="Hysteresis Loop")
    ax.scatter(disp[idx_max], shear[idx_max], color='blue', s=80,
               zorder=5, label=r"$(u_{\rm max}, F_{\rm max})$")

    # Fill the enclosed area
    ax.fill(disp, shear, color='red', alpha=0.25, label=f"E$_d$ = {E_d:.3f} N·m")

    # Mark origin
    ax.axhline(0, color='gray', linestyle='--', linewidth=0.8)
    ax.axvline(0, color='gray', linestyle='--', linewidth=0.8)

    # Annotations and text box
    textstr = (f"$\\xi_{{\\rm eq}} = {zeta_eq:.4f}$ [%]\n"
               f"$E_d = {E_d:.3f}$ N·m\n"
               f"$E_s = {E_s:.3f}$ N·m")
    props = dict(boxstyle='round', facecolor='white', alpha=0.8)
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes,
            fontsize=11, verticalalignment='top', bbox=props)

    ax.set_title(title)
    ax.set_xlabel("Displacement [m]")
    ax.set_ylabel("Base Shear [N]")
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.legend(loc='lower right')
    fig.tight_layout()

    return zeta_eq, fig

u, F = dispX_FV, reactionX_FV
zeta, fig = EQULIVALENT_VISCOUS_DAMPING_RATIO_FUN(u, F, title="Free-Vibration Hysteresis")
fig.show()
print(f"Equivalent viscous damping ratio = {zeta:.4f}")
#%%----------------------------------------------------
# SEISMIC ANALYSIS (DYNAMIC TIME-HISTORY ANALYSIS)
MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'
ANAL_TYPE = 'SEISMIC'

DATA = PENDULUM(SPRING_TYPE, MAT_TYPE, TOTAL_MASS)

(time_SEI, reactionX_SEI, dispX_SEI, velX_SEI, accelX_SEI,
        reactionY_SEI, dispY_SEI, velY_SEI, accelY_SEI,
        reactionZ_SEI, dispZ_SEI, velZ_SEI, accelZ_SEI,
        stiffnessX_SEI, stiffnessY_SEI, stiffnessZ_SEI,
 PERIOD_SEI, damping_ratio_SEI,
 PERIOD_MIN_SEI, PERIOD_MAX_SEI) = DATA


XDATA = dispX_SEI
YDATA = reactionX_SEI
XLABEL = 'Displacement [m]'
YLABEL = 'Base Reaction [N]'
TITLE = 'Base Reaction and Dispalcement of Structure During Seismic Analysis'
COLOR = 'black'
SEMILOGY = False
PLOT(XDATA, YDATA, TITLE, XLABEL, YLABEL, COLOR, SEMILOGY)

PLOT_TIME_HISTORY(time_SEI, reactionX_SEI, dispX_SEI, velX_SEI, accelX_SEI)

# PLOT STRUCTURAL PERIOD DURING THE ANALYSIS
plt.figure(0, figsize=(12, 8))
plt.plot(time_SEI, PERIOD_MIN_SEI, linewidth=3)
plt.plot(time_SEI, PERIOD_MAX_SEI, linewidth=3)
plt.title('Period of Structure During Seismic Analysis')
plt.ylabel('Structural Period [s]')
plt.xlabel('Time [s]')
#plt.semilogy()
plt.grid()
plt.legend([f'PERIOD - MIN VALUES: Min: {np.min(PERIOD_MIN_SEI):.3f} (s) - Mean: {np.mean(PERIOD_MIN_SEI):.3f} (s) - Max: {np.max(PERIOD_MIN_SEI):.3f} (s)', 
            f'PERIOD - MAX VALUES:  Min: {np.min(PERIOD_MAX_SEI):.3f} (s) - Mean: {np.mean(PERIOD_MAX_SEI):.3f} (s) - Max: {np.max(PERIOD_MAX_SEI):.3f} (s)',
            ])
plt.show()

S01.PLOT_2D_FRAME(deformed_scale=1)  # Adjust scale factor as needed
#%%----------------------------------------------------
# --------------------------------------
#  Plot BaseAxial-Displacement Analysis 
# --------------------------------------
XX = np.abs(dispX_PUSH); YY = np.abs(reactionX_PUSH); # ABSOLUTE VALUE
SLOPE_NODE = 10

DATA = S07.BILNEAR_CURVE(XX, YY, SLOPE_NODE)
X, Y, Elastic_ST, Plastic_ST, Tangent_ST, Ductility_Rito, Over_Strength_Factor = DATA

XLABEL = 'Displacement [m]'
YLABEL = 'Base Reaction [N]'
LEGEND01 = 'Curve'
LEGEND02 = 'Bilinear Fitted'
LEGEND03 = 'Undefined'
TITLE = f'Last Data of BaseShear-Displacement Analysis - Ductility Ratio: {X[2]/X[1]:.4f} - Over Strength Factor: {Y[2]/Y[1]:.4f}'
COLOR = 'black'
PLOT_2D(np.abs(dispX_PUSH), np.abs(reactionX_PUSH), X, Y, X, Y, XLABEL, YLABEL, TITLE, LEGEND01, LEGEND02, LEGEND03, COLOR='black', Z=2) 
#print(f'\t\t Ductility Ratio: {Y[2]/Y[1]:.4f}')

# Calculate Over Strength Coefficient (Ω0)
Omega_0 = Y[2] / Y[1]
# Calculate Displacement Ductility Ratio (μ)
mu = X[2] / X[1]
# Calculate Ductility Coefficient (Rμ)
#R_mu = 1
#R_mu = (2 * mu - 1) ** 0.5
R_mu = mu
# Calculate Structural Behavior Coefficient (R)
R = Omega_0 * R_mu
print(f'Over Strength Coefficient (Ω0):      {Omega_0:.4f}')
print(f'Displacement Ductility Ratio (μ):    {mu:.4f}')
print(f'Ductility Coefficient (Rμ):          {R_mu:.4f}')
print(f'Structural Behavior Coefficient (R): {R:.4f}')
Dd = np.max(np.abs(dispX_SEI))
DIy = (Dd - X[1]) /(X[2] - X[1])
print(f'Structural Ductility Damage Index in X Direction: {100*DIy:.4f} (%)')
#%%----------------------------------------------------
# EVALUATION OF DISSIPATED ENERGY CAPACITY INDEX
def DISSIPATED_ENERGY_FUN_WITH_PLOT(displacement, base_shear, method, title="Hysteresis Curve"):
    if method == 1:
        """
        Compute dissipated energy using convex hull and plot the hysteresis curve
        with the outer hull area shaded.
    
        Parameters
        ----------
        displacement : array-like
        base_shear  : array-like
        title       : str
    
        Returns
        -------
        float
            Area of convex hull (dissipated energy)
        """
        import numpy as np
        from scipy.spatial import ConvexHull
        import matplotlib.pyplot as plt
        displacement = np.asarray(displacement)
        base_shear  = np.asarray(base_shear)
    
        if displacement.size != base_shear.size:
            raise ValueError("Displacement and base shear arrays must have equal lengths.")
    
        points = np.column_stack((displacement, base_shear))
        hull = ConvexHull(points)
        area = hull.volume   # 2D hull → area
    
        fig, ax = plt.subplots(figsize=(7, 6))
    
        # Plot full hysteresis
        ax.plot(displacement, base_shear, 'k-', linewidth=1, label="Hysteresis Curve")
    
        # Plot convex hull edges
        hull_pts = points[hull.vertices]
        ax.plot(hull_pts[:, 0], hull_pts[:, 1], 'r--', lw=2, label="Convex Hull")
    
        # Shade hull area
        ax.fill(hull_pts[:, 0], hull_pts[:, 1], color='red', alpha=0.25, label="Hull Area (Energy)")
    
        # Labels and style
        ax.set_title(f"{title} - (Convex Hull)")
        ax.set_xlabel("Displacement (m)")
        ax.set_ylabel("Base Shear (N)")
        ax.grid(True, linestyle='--', alpha=0.5)
        ax.legend()
        
    if method == 2:  
        import numpy as np
        import matplotlib.pyplot as plt
        # Data preparation
        disp = np.asarray(displacement, dtype=float)
        shear = np.asarray(base_shear, dtype=float)

        if disp.size != shear.size:
            raise ValueError("Displacement and base shear arrays must have the same length.")
        if disp.size < 3:
            raise ValueError("At least 3 points are required to form a closed loop.")

        # Close the loop if not already closed (important for shoelace)
        if not (disp[0] == disp[-1] and shear[0] == shear[-1]):
            disp = np.append(disp, disp[0])
            shear = np.append(shear, shear[0])

        # Dissipated energy (E_d) via Shoelace formula
        x = disp
        y = shear
        area = 0.5 * np.abs(np.dot(x[:-1], y[1:]) - np.dot(y[:-1], x[1:]))
        
        # Plotting
        fig, ax = plt.subplots(figsize=(7, 6))
    
        # Hysteresis curve
        idx_max = np.argmax(np.abs(disp))
        ax.plot(disp, shear, 'k-', linewidth=1.2, label="Hysteresis Loop")
        ax.scatter(disp[idx_max], shear[idx_max], color='blue', s=80,
                   zorder=5, label=r"$(u_{\rm max}, F_{\rm max})$")
    
        # Fill the enclosed area
        ax.fill(disp, shear, color='red', alpha=0.25, label=f"E$_d$ = {area:.3f} N·m")
    
    
        ax.set_title(title)
        ax.set_xlabel("Displacement (m)")
        ax.set_ylabel("Base Shear (N)")
        ax.grid(True, linestyle='--', alpha=0.5)
        ax.legend(loc='lower right')
        fig.tight_layout()
    
        return area, fig

Ed_SEI, fig_SEI = DISSIPATED_ENERGY_FUN_WITH_PLOT(
    dispX_SEI, reactionX_SEI, method = 2, 
    title="Earthquake Response – Dissipated Energy"
)
fig_SEI.show()

print(f"Dissipated Energy from Earthquake= {Ed_SEI:.2f} N·m")

Ed_CP, fig_CP = DISSIPATED_ENERGY_FUN_WITH_PLOT(
    dispX_CP, reactionX_CP, method = 2,
    title="Cyclic Loading – Dissipated Energy"
)
fig_CP.show()

print(f"Dissipated Energy from Cyclic Displacement= {Ed_CP:.2f} N·m")


DECI = 100 * Ed_SEI / Ed_CP

if DECI <= 100:
    print(f'\n\tDISSIPATED ENERGY CAPACITY INDEX: {DECI:.3f} [%]\n')
else:
    print('\n\tFOR EVALUATION OF DISSIPATED ENERGY CAPACITY INDEX:')
    print('\n\tCHECK THE CYCLIC DISPLACEMENT ANALYSIS AND IF IT IS POSSIBLE')
    print('\t\t\tINCREASE THE DISPLACEMENT.\n')
    
if DECI <= 0:
    print("\n\tZONE 0: NO DAMAGE\n")
elif DECI > 0 and DECI <= 10:
    print("\n\tZONE 1: VERY MINOR DAMAGE\n")
elif DECI > 10 and DECI <= 20:
    print("\n\tZONE 2: MINOR DAMAGE\n")
elif DECI > 20 and DECI <= 30:
    print("\n\tZONE 3: MODERATE–LOW DAMAGE\n")
elif DECI > 30 and DECI <= 40:
    print("\n\tZONE 4: MODERATE DAMAGE\n")
elif DECI > 40 and DECI <= 50:
    print("\n\tZONE 5: MODERATE–HIGH DAMAGE\n")
elif DECI > 50 and DECI <= 60:
    print("\n\tZONE 6: SEVERE–LOW DAMAGE\n")
elif DECI > 60 and DECI <= 70:
    print("\n\tZONE 7: SEVERE–MEDIUM DAMAGE\n")    
elif DECI > 70 and DECI <= 80:
    print("\n\tZONE 8: SEVERE–HIGH DAMAGE\n")
elif DECI > 80 and DECI <= 90:
    print("\n\tZONE 9: VERY SEVERE DAMAGE\n")
elif DECI > 90 and DECI <= 100:
    print("\n\tZONE 10: FAILURE DAMAGE\n")      
#%%----------------------------------------------------
# %% FRAGILITY ANALYSIS
#exec(open("FRAGILITY_FUN.py").read())
import FRAGILITY_FUN as FF

damage_states = {
    'Minor Damage Level': (20, 5),# Median DI=20%, β=5%
    'Moderate Damage Level': (40, 5),
    'Severe Damage Level': (60, 5),
    'Failure Level': (100, 5)
}

# Intensity Measure (IM) values from 0.0 to 100.0
Ddi = np.abs(dispX_SEI)
DIyi = 100*(Ddi - X[1]) /(X[2] - X[1])
for KK in range(len(DIyi)):
    if DIyi[KK] <= 0.0: # DAMAGE INDEX -> 0.0 <= DI <=100 
        DIyi[KK] = 0.0
    if DIyi[KK] >= 100: 
        DIyi[KK] = 100.0    
        
im_values = DIyi # Structural Ductility Damage Index
TITLE = 'Structural Ductility Damage Index (%)  [IM]'
FF.FRAGILITY_ANALYSIS(damage_states, im_values, TITLE, SCATTER='True', SEMI_LOG='False')

# Add interpretation
results = FF.INTERPRET_FRAGILITY(damage_states, im_values)
print(f"\nProbability of Failure at max displacement: {results['probabilities']['Failure Level']*100:.1f}%")
#%%----------------------------------------------------
# EXCEL OUTPUT
import pandas as pd

# Create DataFrame function
def create_df(reaction, disp, PERIOD_MIN, PERIOD_MAX):
    df = pd.DataFrame({
        "reaction": reaction,
        "disp": disp,
        "PERIOD_MIN": PERIOD_MIN,
        "PERIOD_MAX": PERIOD_MAX,        
    })
    return df


# Save to Excel
with pd.ExcelWriter("PENDULUM_OUTPUT.xlsx", engine='openpyxl') as writer:
    
    # PUSHOVER
    df1 = create_df(reactionX_PUSH, dispX_PUSH, PERIOD_MIN_PUSH, PERIOD_MAX_PUSH)
    df1.to_excel(writer, sheet_name="PUSHOVER", index=False)
                 
    # CYCLIC DISPLACEMENT
    df1 = create_df(reactionX_CP, dispX_CP, PERIOD_MIN_CP, PERIOD_MAX_CP)
    df1.to_excel(writer, sheet_name="CYCLIC_DISPLACEMENT", index=False)
    
    # STATIC EXTERNAL TIME-DEPENDENT LOADING
    df2 = create_df(reactionX_ETDLS, dispX_ETDLS, PERIOD_MIN_ETDLS, PERIOD_MAX_ETDLS)
    df2.to_excel(writer, sheet_name="STATIC_EXTERNAL_TIME-DEPENDENT_LOADING", index=False)

    # DYNAMIC EXTERNAL TIME-DEPENDENT LOADING
    df3 = create_df(reactionX_ETDLD, dispX_ETDLD, PERIOD_MIN_ETDLD, PERIOD_MAX_ETDLD)
    df3.to_excel(writer, sheet_name="DYNAMIC_EXTERNAL_TIME-DEPENDENT_LOADING", index=False)
    
    # FREE-VIBRATION
    df3 = create_df(reactionX_FV, dispX_FV, PERIOD_MIN_FV, PERIOD_MAX_FV)
    df3.to_excel(writer, sheet_name="FREE-VIBRATION", index=False)

    # SEISMIC
    df4 = create_df(reactionX_SEI, dispX_SEI, PERIOD_MIN_SEI, PERIOD_MAX_SEI)
    df4.to_excel(writer, sheet_name="SEISMIC", index=False)
#%%----------------------------------------------------
