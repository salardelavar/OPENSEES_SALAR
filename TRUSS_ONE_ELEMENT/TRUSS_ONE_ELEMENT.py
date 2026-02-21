###########################################################################################################
#                   >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                      #
#      COMPREHENSIVE NONLINEAR SEISMIC ASSESSMENT OF STEEL A TRUSS ELEMENT : AN OPENSEES FRAMEWORK FOR    #
#       STATIC PUSHOVER, CYCLIC DEGRADATION, STATIC TIME-HISTORY AND DYNAMIC TIME-HISTORY ANALYSIS        #
#---------------------------------------------------------------------------------------------------------#
#         ASSESSMENT OF DUCTILITY DAMAGE INDICES FOR STRUCTURAL ELEMENTS AND SYSTEMS AND EVALUATION       #
#                                        OF ENERGY DISSIPATION CAPACITY INDEX                             #
#---------------------------------------------------------------------------------------------------------#
#                          THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                     #
#                                   EMAIL: salar.d.ghashghaei@gmail.com                                   #
###########################################################################################################
"""
Nonlinear Seismic Performance Assessment of a truss element:
An OpenSeesPy Framework for Material and Geometric Nonlinearity Under Static, Cyclic, and Earthquake Loading

This OpenSeesPy script performs rigorous nonlinear static and dynamic analysis of a truss element
 for performance-based earthquake engineering. The 2D model incorporates both material nonlinearity
 (elastic-perfectly plastic or hysteretic steel with strain hardening)
 and geometric nonlinearity (corotational truss formulation) to capture P-delta effects
 and large displacements—essential for collapse assessment.

The bridge spans 43m with 10-panel Warren configuration, assigning distinct cross-sectional
 areas for chords and diagonals. Eigenvalue analysis tracks period elongation throughout
 loading, directly quantifying stiffness degradation and damage progression.
 
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

The code continuously records element axial forces, stresses, strains,
 and full-field nodal displacements, enabling member-level demand-to-capacity
 ratio tracking. Period monitoring during inelastic response reveals softening
 and potential period shift into resonant ground motion frequency ranges.
 This framework enables vulnerability curve development, seismic fragility assessment,
 and retrofit prioritization—directly supporting next-generation bridge seismic design codes
 and risk-informed asset management decisions.
"""
#%%----------------------------------------------------
# WIKIPEDIA: Truss bridge
'https://en.wikipedia.org/w/index.php?title=Truss_bridge'
#%%----------------------------------------------------
import numpy as np
import openseespy.opensees as ops
import matplotlib.pyplot as plt 
import PLOT_2D_TRUSS as S01
import ANALYSIS_FUNCTION as S02
import PERIOD_FUN as S03
import DAMPING_RATIO_FUN as S04
import EIGENVALUE_ANALYSIS_FUN as S05
import RAYLEIGH_DAMPING_FUN as S06
import BILINEAR_CURVE as S07
#%%----------------------------------------------------
def TRUSS_ONE_ELEMENT(LENGTH, AREA, MAT_TYPE,
                        ELE_TYPE, ANAL_TYPE, TOTAL_MASS):
    # Initialize model
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 2)
    
    MAX_ITERATIONS = 5000   # Maximum number of iterations
    MAX_TOLERANCE = 1.0e-6  # Specified tolerance for convergence
    
    # Material (elastic steel)
    MAT_TAG = 1
    if MAT_TYPE == 'ELASTIC':
        E = 200000.0  # [N/mm²] Modulus of steel
        #ops.uniaxialMaterial('Elastic', MAT_TAG, E)             # TESNSION AND COMPRESSION IS SAME VALUES
        ops.uniaxialMaterial('Elastic', MAT_TAG, E ,0.0, 0.5*E) # TESNSION AND COMPRESSION IS NOT SAME VALUES
        # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/ElasticUni.html
    # Material (inelastic steel)    
    if MAT_TYPE == 'INELASTIC':
        Fy = 240.0			    # [N/mm²] Steel yield stress
        Es = 200000.0		    # [N/mm²] Modulus of steel
        ey = Fy/Es			    # [mm/mm] Steel yield strain
        Fu = 1.1818*Fy          # [N/mm²] Steel Ultimate Strength
        esu = 0.12              # [mm/mm] Steel Ultimate Strain
        Esh = (Fu - Fy)/(esu - ey)
        Bs = Esh / Es           # strain-hardening ratio 
        pinchX = 0.8            # Pinching factor in X direction
        pinchY = 0.5            # Pinching factor in Y direction
        damage1 = 0.0           # Damage due to ductility
        damage2 = 0.0           # Damage due to energy
        beta = 0.1              # Stiffness degradation parameter
        ops.uniaxialMaterial('Hysteretic', MAT_TAG, Fy, ey, Fu, esu, 0.2*Fu, 1.1*esu, -Fy, -ey, -Fu, -esu, -0.2*Fu, -1.1*esu, pinchX, pinchY, damage1, damage2, beta)
        # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/Hysteretic.html
    #%% GEOMETRY DEFINITION
    ops.node(1, 0.0, 0.0)
    ops.node(2, LENGTH, 0.0)
    #%% SUPPORTS
    ops.fix(1, 1, 1)      # left  - pinned
    ops.fix(2, 0, 1)      # right - rolled
    
    #%% CREATE ELEMENTS
    ele_tag = 1

    if ELE_TYPE == 'corotTruss':
        # INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Corotational_Truss_Element
        # element corotTruss $eleTag $iNode $jNode $A $matTag <-rho $rho> <-cMass $cFlag> <-doRayleigh $rFlag>
        ops.element('corotTruss', ele_tag, 1, 2, AREA, MAT_TAG, '-rho',AREA*0.00785 , '-doRayleigh', 1)
    if ELE_TYPE == 'Truss': 
        # INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Truss_Element
        # element truss $eleTag $iNode $jNode $A $matTag <-rho $rho> <-cMass $cFlag> <-doRayleigh $rFlag>
        ops.element('Truss', ele_tag, 1, 2, AREA, MAT_TAG, '-rho', AREA*0.00785 , '-doRayleigh', 1)
        ele_tag += 1
    
    #%%--------------------------------------
    # Recorder for element axial forces
    ops.recorder(
        'Element',
        '-file', f'ElementForces_{ANAL_TYPE}.txt',
        '-time',
        '-ele', *range(1, ele_tag),
        'axialForce'
    )
    # Recorder for element axial strains
    ops.recorder(
        'Element',
        '-file', f'ElementStrain_{ANAL_TYPE}.txt',
        '-time',
        '-ele', *range(1, ele_tag),
        'strain'
    )
    # Recorder for element axial stress
    ops.recorder(
        'Element',
        '-file', f'ElementStress_{ANAL_TYPE}.txt',
        '-time',
        '-ele', *range(1, ele_tag),
        'stress'
    )
    #%%--------------------------------------
    tag = 2
    GMfact = 9810.0           # [mm/s²] standard acceleration of gravity or standard acceleration
    TOTAL_WEIGHT = TOTAL_MASS * GMfact
    ENW = TOTAL_WEIGHT / tag  # EACH NODE WEIGHT
    ENM = TOTAL_MASS / tag    # EACH NODE MASS
    #%% LOAD & MASS & ANALYSIS
    center_node = 2
    TS_TAG = 1
    PATT_TAG = 1
    ops.timeSeries('Linear', TS_TAG)
    ops.pattern('Plain', PATT_TAG, TS_TAG)
    for ii in range(1, tag):
        ops.load(ii, 0.0, -ENW)   # [N] downward
        ops.mass(ii, ENW, ENW)
        
    TS_TAG = 2
    PATT_TAG = 2
    ops.timeSeries('Linear', TS_TAG)
    ops.pattern('Plain', PATT_TAG, TS_TAG)    
    #ops.timeSeries('Linear', 2)
    #ops.pattern('Plain', 2, 1)    
    ops.load(center_node, 1.0, 0.0)   # [N] downward
    
    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.system('BandGeneral')
    
    if ELE_TYPE == 'ELASTIC':
        ops.algorithm('Linear')
    if ELE_TYPE == 'INELASTIC':
        ops.algorithm('Newton')
        
    time = []
    force = []  
    stress, strain = [], []
    dispX, dispY = [], []
    veloX, veloY = [], []
    accX, accY = [], []
    reaction, disp_mid = [], []
    stiffness = []
    OMEGA, PERIOD = [], []
    PERIOD_MIN, PERIOD_MAX = [], []
    # Initialize lists for each element's axiaal force, stress and strain
    ele_axialforce = {
        1: [],  # AXIALFORCE-01
        }  
    ele_strain = {
        1: [],  # STRAIN-01
        }
    ele_di = {
        1: [],  # DAMAGE_INDEX-01
        }
    ele_stress = {
        1: [],  # STRESS-01
        }    
    # Initialize lists for each node's displacement
    node_displacementsX = {
        1: [],  # DISP01
        2: [],  # DISP02
        }
    node_displacementsY = {
        1: [],  # DISP01
        2: [],  # DISP02
        }
    if ANAL_TYPE == 'PERIOD': 
      #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(2, 0.5*DR, DR, 0, 1)
      PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(1, PLOT=True) 
      return PERIODmin, PERIODmax
  
    if ANAL_TYPE == 'STATIC': 
        ops.integrator('LoadControl', 1.0)
        ops.analysis('Static')
        OK = ops.analyze(1)
        S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS)
        ops.reactions()
        reaction.append(ops.nodeReaction(1, 1))                         # AXIAL BASE REACTION
        disp_mid.append(ops.nodeDisp(center_node, 1))                   # DISPLACEMENT NODE 02 
        dispX.append(ops.nodeDisp(center_node, 1))                      # DISPLACEMENT NODE 02 IN X DIR 
        dispY.append(ops.nodeDisp(center_node, 2))                      # DISPLACEMENT NODE 02 IN Y DIR 
        # Store axial forces, strain and stress
        for ele_id in ele_axialforce.keys(): 
            ele_axialforce[ele_id].append(ops.eleResponse(ele_id, 'axialForce'))     # [N] ELEMENT FORCE
            ele_stress[ele_id].append(ops.eleResponse(ele_id, 'material', 'stress')) # [N/mm^2] ELEMENT STRESS 
            ele_strain[ele_id].append(ops.eleResponse(ele_id, 'material', 'strain')) # [mm/mm] ELEMENT STRAIN
            if ELE_TYPE == 'INELASTIC':
                ele_di[ele_id].append(100*(abs(ele_strain[ele_id]) - ey) - (esu - ey))
            if ELE_TYPE == 'ELASTIC':
                ele_di[ele_id].append(0.0)                
        # Store displacements
        for node_id in node_displacementsX.keys():    
            node_displacementsX[node_id].append(ops.nodeDisp(node_id, 1))
            node_displacementsY[node_id].append(ops.nodeDisp(node_id, 2))
        else:
            print('\n\nSTATIC ANALYSIS DONE.\n\n')  
            
        DATA = (reaction, disp_mid,
                ele_axialforce, ele_stress, ele_strain, ele_di,
                node_displacementsX, node_displacementsY,
                dispX, dispY)
        
        return  DATA
    
    if ANAL_TYPE == 'PUSHOVER': # STATIC TIME-HISTORY ANALYSIS
        IDctrlDOF = 1   # 1: Horizental Dispalcement - 2: Vertical Dispalcement
        DINCR = -0.01   # [mm] Incremental Vertical Displacement
        DMAX = -50.0    # [mm] Max. Vertical Displacement
        ops.integrator('DisplacementControl', center_node, IDctrlDOF, DINCR)
        ops.analysis('Static')
        Nsteps =  int(np.abs(DMAX/ DINCR)) 
        STEP = 0.0
        for step in range(Nsteps):
            OK = ops.analyze(1)
            S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS)
            ops.reactions()
            reaction.append(ops.nodeReaction(1, 1))                         # AXIAL BASE REACTION
            disp_mid.append(ops.nodeDisp(center_node, 1))                   # DISPLACEMENT NODE 02 
            dispX.append(ops.nodeDisp(center_node, 1))                      # DISPLACEMENT NODE 02 IN X DIR 
            dispY.append(ops.nodeDisp(center_node, 2))                      # DISPLACEMENT NODE 02 IN Y DIR 
            # Store elements force, strain and stress
            for ele_id in ele_axialforce.keys(): 
                ele_axialforce[ele_id].append(ops.eleResponse(ele_id, 'axialForce'))     # [N] ELEMENT FORCE
                ele_stress[ele_id].append(ops.eleResponse(ele_id, 'material', 'stress')) # [N/mm^2] ELEMENT STRESS 
                ele_strain[ele_id].append(ops.eleResponse(ele_id, 'material', 'strain')) # [mm/mm] ELEMENT STRAIN      
                if ELE_TYPE == 'INELASTIC':
                    ele_di[ele_id].append(100*(abs(ele_strain[ele_id]) - ey) - (esu - ey))
                if ELE_TYPE == 'ELASTIC':
                    ele_di[ele_id].append(0.0) 
            # Store displacements
            for node_id in node_displacementsX.keys():    
                node_displacementsX[node_id].append(ops.nodeDisp(node_id, 1))
                node_displacementsY[node_id].append(ops.nodeDisp(node_id, 2))      
            # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
            #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(1, 0.5*DR, DR, 0, 1)
            PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(1, PLOT=True)
            PERIOD_MIN.append(PERIODmin)
            PERIOD_MAX.append(PERIODmax) 
            STEP += 1
            #print(STEP, dispY[-1], reaction[-1])
            print(f"Step: {STEP}, Displacement: {dispX[-1]:.4f} mm, Reaction: {reaction[-1]:.2f} N")     
        else:
            print('\n\nPUSHOVER ANALYSIS DONE.\n\n')
            
        DATA = (reaction, disp_mid,
                ele_axialforce, ele_stress, ele_strain, ele_di,
                node_displacementsX, node_displacementsY,
                dispX, dispY,
                np.array(PERIOD_MIN), np.array(PERIOD_MAX))
    
        return  DATA
    
    if ANAL_TYPE == 'CYCLIC_DISPLACEMENT': # STATIC TIME-HISTORY ANALYSIS
        IDctrlDOF = 1   # 1: Horizental Dispalcement - 2: Vertical Dispalcement
        DMAX = -50.0    # [mm] Max. Vertical Displacement
        n_points = 10000
        # 4. CYCLIC DISPALCEMENT PROTOCOL
        # Key strain points (same logic as your protocol)
        key_disp = np.array([
             0.0,
             0.1*DMAX,   -0.1*DMAX,
             0.5*DMAX,   -0.5*DMAX,
             0.8*DMAX,   -0.8*DMAX,
             DMAX,       -DMAX,
             0.1*DMAX,     -0.1*DMAX,
             0.5*DMAX,     -0.5*DMAX,
             0.8*DMAX,     -0.8*DMAX,
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
            reaction.append(ops.nodeReaction(1, 1))                         # AXIAL BASE REACTION
            disp_mid.append(ops.nodeDisp(center_node, 1))                   # DISPLACEMENT NODE 02 
            dispX.append(ops.nodeDisp(center_node, 1))                      # DISPLACEMENT NODE 02 IN X DIR 
            dispY.append(ops.nodeDisp(center_node, 2))                      # DISPLACEMENT NODE 02 IN Y DIR 
            # Store elements force, strain and stress
            for ele_id in ele_axialforce.keys(): 
                ele_axialforce[ele_id].append(ops.eleResponse(ele_id, 'axialForce'))     # [N] ELEMENT FORCE
                ele_stress[ele_id].append(ops.eleResponse(ele_id, 'material', 'stress')) # [N/mm^2] ELEMENT STRESS 
                ele_strain[ele_id].append(ops.eleResponse(ele_id, 'material', 'strain')) # [mm/mm] ELEMENT STRAIN      
                if ELE_TYPE == 'INELASTIC':
                    ele_di[ele_id].append(100*(abs(ele_strain[ele_id]) - ey) - (esu - ey))
                if ELE_TYPE == 'ELASTIC':
                    ele_di[ele_id].append(0.0) 
            # Store displacements
            for node_id in node_displacementsX.keys():    
                node_displacementsX[node_id].append(ops.nodeDisp(node_id, 1))
                node_displacementsY[node_id].append(ops.nodeDisp(node_id, 2))      
            # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
            #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(1, 0.5*DR, DR, 0, 1)
            PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(1, PLOT=True)
            PERIOD_MIN.append(PERIODmin)
            PERIOD_MAX.append(PERIODmax) 
            STEP += 1
            #print(STEP, dispY[-1], reaction[-1])
            print(f"Step: {STEP}, Displacement: {dispX[-1]:.4f} mm, Reaction: {reaction[-1]:.2f} N")     
        else:
            print('\n\nCYCLIC DISPLAEMENT ANALYSIS DONE.\n\n')
            
        DATA = (reaction, disp_mid,
                ele_axialforce, ele_stress, ele_strain, ele_di,
                node_displacementsX, node_displacementsY,
                dispX, dispY,
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
        ops.load(center_node, 1.0, 0.0)
        
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
            reaction.append(ops.nodeReaction(1, 1))                         # AXIAL BASE REACTION
            disp_mid.append(ops.nodeDisp(center_node, 1))                   # DISPLACEMENT NODE 02 
            dispX.append(ops.nodeDisp(center_node, 1))                      # DISPLACEMENT NODE 02 IN X DIR 
            dispY.append(ops.nodeDisp(center_node, 2))                      # DISPLACEMENT NODE 02 IN Y DIR 
            # Store elements force, strain and stress
            for ele_id in ele_axialforce.keys(): 
                ele_axialforce[ele_id].append(ops.eleResponse(ele_id, 'axialForce'))     # [N] ELEMENT FORCE
                ele_stress[ele_id].append(ops.eleResponse(ele_id, 'material', 'stress')) # [N/mm^2] ELEMENT STRESS 
                ele_strain[ele_id].append(ops.eleResponse(ele_id, 'material', 'strain')) # [mm/mm] ELEMENT STRAIN      
                if ELE_TYPE == 'INELASTIC':
                    ele_di[ele_id].append(100*(abs(ele_strain[ele_id]) - ey) - (esu - ey))
                if ELE_TYPE == 'ELASTIC':
                    ele_di[ele_id].append(0.0) 
            # Store displacements
            for node_id in node_displacementsX.keys():    
                node_displacementsX[node_id].append(ops.nodeDisp(node_id, 1))
                node_displacementsY[node_id].append(ops.nodeDisp(node_id, 2))      
            # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
            #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(1, 0.5*DR, DR, 0, 1)
            PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(1, PLOT=True)
            PERIOD_MIN.append(PERIODmin)
            PERIOD_MAX.append(PERIODmax) 
            STEP += 1
            #print(STEP, dispY[-1], reaction[-1])
            print(f"Step: {STEP}, Displacement: {dispX[-1]:.4f} mm, Reaction: {reaction[-1]:.2f} N")     
        else:
            print('\n\nSTATIC EXTERNAL TIME-DEPENDENT LOADING ANALYSIS DONE.\n\n')
            
        DATA = (reaction, disp_mid,
                ele_axialforce, ele_stress, ele_strain, ele_strain,
                node_displacementsX, node_displacementsY,
                dispX, dispY,
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
        ops.load(center_node, 1.0, 0.0)
        
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
            reaction.append(ops.nodeReaction(1, 1))               # AXIAL BASE REACTION
            disp_mid.append(ops.nodeDisp(center_node, 1))         # DISPLACEMENT NODE 02 
            dispX.append(ops.nodeDisp(center_node, 1))            # DISPLACEMENT NODE 02 IN X DIR 
            dispY.append(ops.nodeDisp(center_node, 2))            # DISPLACEMENT NODE 02 IN Y DIR 
            veloX.append(ops.nodeVel(center_node, 1))             # VELOCITY NODE 02
            veloY.append(ops.nodeVel(center_node, 2))             # VELOCITY NODE 02
            accX.append(ops.nodeAccel(center_node, 1))            # ACCELERATION NODE 02
            accY.append(ops.nodeAccel(center_node, 2))            # ACCELERATION NODE 02
            stiffness.append(np.abs(reaction[-1] / disp_mid[-1]))
            OMEGA.append(np.sqrt(stiffness[-1]/TOTAL_MASS))
            PERIOD.append((np.pi * 2) / OMEGA[-1])
            # Store elements force, strain and stress
            for ele_id in ele_axialforce.keys(): 
                ele_axialforce[ele_id].append(ops.eleResponse(ele_id, 'axialForce'))     # [N] ELEMENT FORCE
                ele_stress[ele_id].append(ops.eleResponse(ele_id, 'material', 'stress')) # [N/mm^2] ELEMENT STRESS 
                ele_strain[ele_id].append(ops.eleResponse(ele_id, 'material', 'strain')) # [mm/mm] ELEMENT STRAIN 
                if ELE_TYPE == 'INELASTIC':
                    ele_di[ele_id].append(100*(abs(ele_strain[ele_id]) - ey) - (esu - ey))
                if ELE_TYPE == 'ELASTIC':
                    ele_di[ele_id].append(0.0) 
            # Store displacements
            for node_id in node_displacementsX.keys():    
                node_displacementsX[node_id].append(ops.nodeDisp(node_id, 1))
                node_displacementsY[node_id].append(ops.nodeDisp(node_id, 2))   
            # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
            #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(1, 0.5*DR, DR, 0, 1)
            PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(1, PLOT=True)
            PERIOD_MIN.append(PERIODmin)
            PERIOD_MAX.append(PERIODmax)
            #print(time[-1], dispY[-1], veloY[-1])
            print(f"Time: {time[-1]:.4f}, Displacement: {dispX[-1]:.4f} mm, Reaction: {reaction[-1]:.2f} N")      
        else:
            print('\n\nDYNAMIC EXTERNAL TIME-DEPENDENT LOADING ANALYSIS DONE.\n\n')  
        # Calculating Damping Ratio and Period Using Logarithmic Decrement Analysis 
        damping_ratio = S04.DAMPING_RATIO(dispY)  
        
        # Compute modal properties
        ops.modalProperties("-print", "-file", f"SALAR_ModalReport_{ANAL_TYPE}.txt", "-unorm") 
        
        DATA = (time, reaction, disp_mid,
                ele_axialforce, ele_stress, ele_strain, ele_di,
                node_displacementsX, node_displacementsY,
                dispX, dispY,
                veloX, veloY,
                accX, accY,
                stiffness, PERIOD, damping_ratio,
                np.array(PERIOD_MIN), np.array(PERIOD_MAX))
        
        return  DATA
        
    if ANAL_TYPE == 'FREE-VIBRATION': # DYNAMIC TIME-HISTORY ANALYSIS
        #%% DEFINE PARAMETERS FOR FREE-VIBRATION ANALYSIS
        u0 = -5.35                         # [mm] Initial displacement
        v0 = 0.15                          # [mm/s] Initial velocity
        a0 = 0.065                         # [mm/s^2] Initial acceleration
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
            reaction.append(ops.nodeReaction(1, 1))               # AXIAL BASE REACTION
            disp_mid.append(ops.nodeDisp(center_node, 1))         # DISPLACEMENT NODE 02 
            dispX.append(ops.nodeDisp(center_node, 1))            # DISPLACEMENT NODE 02 IN X DIR 
            dispY.append(ops.nodeDisp(center_node, 2))            # DISPLACEMENT NODE 02 IN Y DIR 
            veloX.append(ops.nodeVel(center_node, 1))             # VELOCITY NODE 02
            veloY.append(ops.nodeVel(center_node, 2))             # VELOCITY NODE 02
            accX.append(ops.nodeAccel(center_node, 1))            # ACCELERATION NODE 02
            accY.append(ops.nodeAccel(center_node, 2))            # ACCELERATION NODE 02
            stiffness.append(np.abs(reaction[-1] / disp_mid[-1]))
            OMEGA.append(np.sqrt(stiffness[-1]/TOTAL_MASS))
            PERIOD.append((np.pi * 2) / OMEGA[-1])
            # Store elements force, strain and stress
            for ele_id in ele_axialforce.keys(): 
                ele_axialforce[ele_id].append(ops.eleResponse(ele_id, 'axialForce'))     # [N] ELEMENT FORCE
                ele_stress[ele_id].append(ops.eleResponse(ele_id, 'material', 'stress')) # [N/mm^2] ELEMENT STRESS 
                ele_strain[ele_id].append(ops.eleResponse(ele_id, 'material', 'strain')) # [mm/mm] ELEMENT STRAIN 
                if ELE_TYPE == 'INELASTIC':
                    ele_di[ele_id].append(100*(abs(ele_strain[ele_id]) - ey) - (esu - ey))
                if ELE_TYPE == 'ELASTIC':
                    ele_di[ele_id].append(0.0) 
            # Store displacements
            for node_id in node_displacementsX.keys():    
                node_displacementsX[node_id].append(ops.nodeDisp(node_id, 1))
                node_displacementsY[node_id].append(ops.nodeDisp(node_id, 2))   
            # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
            #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(1, 0.5*DR, DR, 0, 1)
            PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(1, PLOT=True)
            PERIOD_MIN.append(PERIODmin)
            PERIOD_MAX.append(PERIODmax)
            #print(time[-1], dispY[-1], veloY[-1])
            print(f"Time: {time[-1]:.4f}, Displacement: {dispX[-1]:.4f} mm, Reaction: {reaction[-1]:.2f} N")      
        else:
            print('\n\nFREE-VIBRATION ANALYSIS DONE.\n\n')    
        # Calculating Damping Ratio and Period Using Logarithmic Decrement Analysis 
        damping_ratio = S04.DAMPING_RATIO(dispY)  
        
        # Compute modal properties
        ops.modalProperties("-print", "-file", f"SALAR_ModalReport_{ANAL_TYPE}.txt", "-unorm") 
        
        DATA = (time, reaction, disp_mid,
                ele_axialforce, ele_stress, ele_strain, ele_di,
                node_displacementsX, node_displacementsY,
                dispX, dispY,
                veloX, veloY,
                accX, accY,
                stiffness, PERIOD, damping_ratio,
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
            reaction.append(ops.nodeReaction(1, 1))               # AXIAL BASE REACTION
            disp_mid.append(ops.nodeDisp(center_node, 1))         # DISPLACEMENT NODE 02 
            dispX.append(ops.nodeDisp(center_node, 1))            # DISPLACEMENT NODE 02 IN X DIR 
            dispY.append(ops.nodeDisp(center_node, 2))            # DISPLACEMENT NODE 02 IN Y DIR 
            veloX.append(ops.nodeVel(center_node, 1))             # VELOCITY NODE 02
            veloY.append(ops.nodeVel(center_node, 2))             # VELOCITY NODE 02
            accX.append(ops.nodeAccel(center_node, 1))            # ACCELERATION NODE 02
            accY.append(ops.nodeAccel(center_node, 2))            # ACCELERATION NODE 02
            stiffness.append(np.abs(reaction[-1] / disp_mid[-1]))
            OMEGA.append(np.sqrt(stiffness[-1]/TOTAL_MASS))
            PERIOD.append((np.pi * 2) / OMEGA[-1])
            # Store elements force, strain and stress
            for ele_id in ele_axialforce.keys(): 
                ele_axialforce[ele_id].append(ops.eleResponse(ele_id, 'axialForce'))     # [N] ELEMENT FORCE
                ele_stress[ele_id].append(ops.eleResponse(ele_id, 'material', 'stress')) # [N/mm^2] ELEMENT STRESS 
                ele_strain[ele_id].append(ops.eleResponse(ele_id, 'material', 'strain')) # [mm/mm] ELEMENT STRAIN 
                if ELE_TYPE == 'INELASTIC':
                    ele_di[ele_id].append(100*(abs(ele_strain[ele_id]) - ey) - (esu - ey))
                if ELE_TYPE == 'ELASTIC':
                    ele_di[ele_id].append(0.0) 
            # Store displacements
            for node_id in node_displacementsX.keys():    
                node_displacementsX[node_id].append(ops.nodeDisp(node_id, 1))
                node_displacementsY[node_id].append(ops.nodeDisp(node_id, 2))
            # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
            #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(1, 0.5*DR, DR, 0, 1)
            PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(1, PLOT=True)
            PERIOD_MIN.append(PERIODmin)
            PERIOD_MAX.append(PERIODmax)
            #print(time[-1], dispY[-1], veloY[-1])
            print(f"Time: {time[-1]:.4f}, Displacement: {dispX[-1]:.4f} mm, Reaction: {reaction[-1]:.2f} N")

        else:
            print('\n\nSEISMIC ANALYSIS DONE.\n\n')     
        # Calculating Damping Ratio and Period Using Logarithmic Decrement Analysis 
        damping_ratio = S04.DAMPING_RATIO(dispY)  
        
        # Compute modal properties
        ops.modalProperties("-print", "-file", f"SALAR_ModalReport_{ANAL_TYPE}.txt", "-unorm") 
        
        DATA = (time, reaction, disp_mid,
                ele_axialforce, ele_stress, ele_strain, ele_di,
                node_displacementsX, node_displacementsY,
                dispX, dispY,
                veloX, veloY,
                accX, accY,
                stiffness, PERIOD, damping_ratio,
                np.array(PERIOD_MIN), np.array(PERIOD_MAX))
        
        return DATA    

#%%-------------------------------------------------------
def PLOT_TIME_HISTORY(time,
                      reaction, disp_mid,
                      dispX, dispY,
                      veloX, veloY,
                      accX, accY):

    import matplotlib.pyplot as plt

    data_dict = {
        "Reaction [N]": reaction,
        "Mid Displacement [mm]": disp_mid,
        "Disp X [mm]": dispX,
        "Disp Y [mm]": dispY,
        "Velocity X [mm/s]": veloX,
        "Velocity Y [mm/s]": veloY,
        "Acceleration X [mm/s²]": accX,
        "Acceleration Y [mm/s²]": accY,
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
    plt.ylabel('Displacement [mm]')
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
LENGTH = 3000.0       # [mm] Bridge Length
AREA  = 2*2040.0      # [mm^2] Section Area for diagonal elements - 2 x UNP140
TOTAL_MASS = 0.0      # [kg] Total Mass of Structure
#%%----------------------------------------------------
# PERIOD ANALYSIS
ELE_TYPE = 'Truss'      # MATERIAL NONLINEARITY
#ELE_TYPE = 'corotTruss'  # MATERIAL AND GEOMETRIC NONLINEARITIES
MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'
ANAL_TYPE = 'PERIOD'

DATA = TRUSS_ONE_ELEMENT(LENGTH, AREA, MAT_TYPE, ELE_TYPE, ANAL_TYPE, TOTAL_MASS)
(PERIOD_MIN_X, PERIOD_MAX_X) = DATA
print('Structure First Period:  ', PERIOD_MIN_X)
print('Structure Second Period: ', PERIOD_MAX_X) 

S01.PLOT_2D_FRAME_TRUSS(deformed_scale=1.0)
#%%----------------------------------------------------
# STATIC ANALYSIS
#ELE_TYPE = 'Truss'      # MATERIAL NONLINEARITY
ELE_TYPE = 'corotTruss'  # MATERIAL AND GEOMETRIC NONLINEARITIES
MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'
ANAL_TYPE = 'STATIC'

DATA = TRUSS_ONE_ELEMENT(LENGTH, AREA, MAT_TYPE, ELE_TYPE, ANAL_TYPE, TOTAL_MASS)
(reaction, disp_mid, 
 ele_force, ele_stress, ele_strain, ele_di,
 node_displacementsX, node_displacementsY,
 dispX, dispY) = DATA

S01.PLOT_2D_FRAME_TRUSS(deformed_scale=1.0)
#%%----------------------------------------------------
# PUSHOVER ANALYSIS (STATIC TIME-HISTORY ANALYSIS)
#ELE_TYPE = 'Truss'      # MATERIAL NONLINEARITY
ELE_TYPE = 'corotTruss'  # MATERIAL AND GEOMETRIC NONLINEARITIES
MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'
ANAL_TYPE = 'PUSHOVER'

DATA = TRUSS_ONE_ELEMENT(LENGTH, AREA, MAT_TYPE, ELE_TYPE, ANAL_TYPE, TOTAL_MASS)
(reaction_PUSH, disp_mid_PUSH, 
 ele_force_PUSH, ele_stress_PUSH, ele_strain_PUSH, ele_di_PUSH,
 node_displacementsX_PUSH, node_displacementsY_PUSH,
  dispX_PUSH, dispY_PUSH,
 PERIOD_MIN_PUSH, PERIOD_MAX_PUSH) = DATA


XDATA = disp_mid_PUSH
YDATA = reaction_PUSH
XLABEL = 'Displacement in Middle Span [mm]'
YLABEL = 'Base Reaction [N]'
TITLE = 'Base Reaction and Displacement of Structure During Pushover Analysis'
COLOR = 'black'
SEMILOGY = False
PLOT(XDATA, YDATA, TITLE, XLABEL, YLABEL, COLOR, SEMILOGY)

DATA = S07.BILNEAR_CURVE(np.abs(disp_mid_PUSH), np.abs(reaction_PUSH), SLOPE_NODE=10)
(X_PUSH, Y_PUSH, Elastic_ST, Plastic_ST, Tangent_ST, Ductility_Rito, Over_Strength_Factor) = DATA

# PLOT STRUCTURAL PERIOD DURING THE ANALYSIS
plt.figure(0, figsize=(12, 8))
plt.plot(disp_mid_PUSH, PERIOD_MIN_PUSH, linewidth=3)
plt.plot(disp_mid_PUSH, PERIOD_MAX_PUSH, linewidth=3)
plt.title('Period of Structure During Pushover Analysis')
plt.ylabel('Structural Period [s]')
plt.xlabel('Displacement [mm]')
#plt.semilogy()
plt.grid()
plt.legend([f'PERIOD - MIN VALUES: Min: {np.min(PERIOD_MIN_PUSH):.3f} (s) - Mean: {np.mean(PERIOD_MIN_PUSH):.3f} (s) - Max: {np.max(PERIOD_MIN_PUSH):.3f} (s)', 
            f'PERIOD - MAX VALUES:  Min: {np.min(PERIOD_MAX_PUSH):.3f} (s) - Mean: {np.mean(PERIOD_MAX_PUSH):.3f} (s) - Max: {np.max(PERIOD_MAX_PUSH):.3f} (s)',
            ])
plt.show()

S01.PLOT_2D_FRAME_TRUSS(deformed_scale=1.0)
#%%----------------------------------------------------
# CYCLIC DISPLACEMENT ANALYSIS (STATIC TIME-HISTORY ANALYSIS)
#ELE_TYPE = 'Truss'      # MATERIAL NONLINEARITY
ELE_TYPE = 'corotTruss'  # MATERIAL AND GEOMETRIC NONLINEARITIES
MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'
ANAL_TYPE = 'CYCLIC_DISPLACEMENT'

DATA = TRUSS_ONE_ELEMENT(LENGTH, AREA, MAT_TYPE, ELE_TYPE, ANAL_TYPE, TOTAL_MASS)
(reaction_CP, disp_mid_CP, 
 ele_force_CP, ele_stress_CP, ele_strain_CP, ele_di_CP,
 node_displacementsX_CP, node_displacementsY_CP,
  dispX_CP, dispY_CP,
 PERIOD_MIN_CP, PERIOD_MAX_CP) = DATA


XDATA = disp_mid_CP
YDATA = reaction_CP
XLABEL = 'Displacement in Middle Span [mm]'
YLABEL = 'Base Reaction [N]'
TITLE = 'Base Reaction and Displacement of Structure During Cyclic-Displacement Analysis'
COLOR = 'black'
SEMILOGY = False
PLOT(XDATA, YDATA, TITLE, XLABEL, YLABEL, COLOR, SEMILOGY)

# PLOT STRUCTURAL PERIOD DURING THE ANALYSIS
plt.figure(0, figsize=(12, 8))
plt.plot(disp_mid_CP, PERIOD_MIN_CP, linewidth=3)
plt.plot(disp_mid_CP, PERIOD_MAX_CP, linewidth=3)
plt.title('Period of Structure During Cyclic Displacement Analysis')
plt.ylabel('Structural Period [s]')
plt.xlabel('Displacement [mm]')
#plt.semilogy()
plt.grid()
plt.legend([f'PERIOD - MIN VALUES: Min: {np.min(PERIOD_MIN_CP):.3f} (s) - Mean: {np.mean(PERIOD_MIN_CP):.3f} (s) - Max: {np.max(PERIOD_MIN_CP):.3f} (s)', 
            f'PERIOD - MAX VALUES:  Min: {np.min(PERIOD_MAX_CP):.3f} (s) - Mean: {np.mean(PERIOD_MAX_CP):.3f} (s) - Max: {np.max(PERIOD_MAX_CP):.3f} (s)',
            ])
plt.show()

S01.PLOT_2D_FRAME_TRUSS(deformed_scale=1.0)
#%%----------------------------------------------------
# EXTERNAL TIME-DEPENDENT LOADING ANALYSIS (STATIC TIME-HISTORY ANALYSIS)
#ELE_TYPE = 'Truss'      # MATERIAL NONLINEARITY
ELE_TYPE = 'corotTruss'  # MATERIAL AND GEOMETRIC NONLINEARITIES
MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'
ANAL_TYPE = 'STATIC_EXTERNAL_TIME-DEPENDENT_LOADING'

DATA = TRUSS_ONE_ELEMENT(LENGTH, AREA, MAT_TYPE, ELE_TYPE, ANAL_TYPE, TOTAL_MASS)

(reaction_ETDLS, disp_mid_ETDLS, 
 ele_force_ETDLS, ele_stress_ETDLS, ele_strain_ETDLS, ele_di_ETDLS,
 node_displacementsX_ETDLS, node_displacementsY_ETDLS,
  dispX_ETDLS, dispY_ETDLS,
 PERIOD_MIN_ETDLS, PERIOD_MAX_ETDLS) = DATA


XDATA = disp_mid_ETDLS
YDATA = reaction_ETDLS
XLABEL = 'Displacement in Middle Span [mm]'
YLABEL = 'Base Reaction [N]'
TITLE = 'Base Reaction and Displacement of Structure During Static External Time-dependent Loading Analysis'
COLOR = 'black'
SEMILOGY = False
PLOT(XDATA, YDATA, TITLE, XLABEL, YLABEL, COLOR, SEMILOGY)

# PLOT STRUCTURAL PERIOD DURING THE ANALYSIS
plt.figure(0, figsize=(12, 8))
plt.plot(disp_mid_ETDLS, PERIOD_MIN_ETDLS, linewidth=3)
plt.plot(disp_mid_ETDLS, PERIOD_MAX_ETDLS, linewidth=3)
plt.title('Period of Structure During Static External Time-dependent Loading Analysis')
plt.ylabel('Structural Period [s]')
plt.xlabel('Displacement [mm]')
#plt.semilogy()
plt.grid()
plt.legend([f'PERIOD - MIN VALUES: Min: {np.min(PERIOD_MIN_ETDLS):.3f} (s) - Mean: {np.mean(PERIOD_MIN_ETDLS):.3f} (s) - Max: {np.max(PERIOD_MIN_ETDLS):.3f} (s)', 
            f'PERIOD - MAX VALUES:  Min: {np.min(PERIOD_MAX_ETDLS):.3f} (s) - Mean: {np.mean(PERIOD_MAX_ETDLS):.3f} (s) - Max: {np.max(PERIOD_MAX_ETDLS):.3f} (s)',
            ])
plt.show()

S01.PLOT_2D_FRAME_TRUSS(deformed_scale=1.0)
#%%----------------------------------------------------
# EXTERNAL TIME-DEPENDENT LOADING ANALYSIS (DYNAMIC TIME-HISTORY ANALYSIS)
#ELE_TYPE = 'Truss'      # MATERIAL NONLINEARITY
ELE_TYPE = 'corotTruss'  # MATERIAL AND GEOMETRIC NONLINEARITIES
MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'
ANAL_TYPE = 'DYNAMIC_EXTERNAL_TIME-DEPENDENT_LOADING'

DATA = TRUSS_ONE_ELEMENT(LENGTH, AREA, MAT_TYPE, ELE_TYPE, ANAL_TYPE, TOTAL_MASS)

(time_ETDLD, reaction_ETDLD, disp_mid_ETDLD,
ele_force_ETDLD, ele_stress_ETDLD, ele_strain_ETDLD, ele_di_ETDLD,
node_displacementsX_ETDLD, node_displacementsY_ETDLD,
dispX_ETDLD, dispY_ETDLD,
veloX_ETDLD, veloY_ETDLD,
accX_ETDLD, accY_ETDLD,
stiffness_ETDLD, PERIOD_ETDLD, damping_ratio_ETDLD,
PERIOD_MIN_ETDLD, PERIOD_MAX_ETDLD) = DATA


XDATA = disp_mid_ETDLD
YDATA = reaction_ETDLD
XLABEL = 'Displacement in Middle Span [mm]'
YLABEL = 'Base Reaction [N]'
TITLE = 'Base Reaction and Displacement of Structure During Dynamic External Time-dependent Loading Analysis'
COLOR = 'black'
SEMILOGY = False
PLOT(XDATA, YDATA, TITLE, XLABEL, YLABEL, COLOR, SEMILOGY)

PLOT_TIME_HISTORY(time_ETDLD, reaction_ETDLD, disp_mid_ETDLD,
                          dispX_ETDLD, dispY_ETDLD,
                          veloX_ETDLD, veloY_ETDLD,
                          accX_ETDLD, accY_ETDLD)

# PLOT ELEMENTS AXIAL FORCE
YLABEL = 'Element Axial Force [N]'  
TITLE = "elements  force in X Dir. vs Time for Truss Element During Dynamic External Time-dependent Loading Analysis"
PLOT_FORCES(time_ETDLD, ele_force_ETDLD, YLABEL, TITLE) # TRUSS ELEMENT - ELEMENTS AXIAL FORCE
# PLOT ELEMENTS STRESS
YLABEL = 'Element Stress [N/mm^2]'  
TITLE = "elements  Stress in X Dir. vs Time for Truss Element During Dynamic External Time-dependent Loading Analysis"
PLOT_FORCES(time_ETDLD, ele_stress_ETDLD, YLABEL, TITLE) # TRUSS ELEMENT - ELEMENTS STRESS
# PLOT ELEMENTS STRAIN
YLABEL = 'Element Strain [mm]'  
TITLE = "elements  Strain in X Dir. vs Time for Truss Element During Dynamic External Time-dependent Loading Analysis"
PLOT_FORCES(time_ETDLD, ele_strain_ETDLD, YLABEL, TITLE) # TRUSS ELEMENT - ELEMENTS AXIAL FORCE
# PLOT ELEMENTS DAMAGE-INDEX
YLABEL = 'Element Damage-index [%]'  
TITLE = "elements Damage-index [%] in X Dir. vs Time for Truss Element During Dynamic External Time-dependent Loading Analysis"
PLOT_FORCES(time_ETDLD, ele_strain_ETDLD, YLABEL, TITLE) # TRUSS ELEMENT - ELEMENTS AXIAL FORCE

# PLOT NODAL DISPALEMENTS
PLOT_DISPLAEMENTS(time_ETDLD, node_displacementsX_ETDLD, TITLE = "Node Displacements in X Dir. vs Time for Truss Element During Dynamic External Time-dependent Loading Analysis") # TRUSS ELEMENT IN X DIR.
PLOT_DISPLAEMENTS(time_ETDLD, node_displacementsY_ETDLD, TITLE = "Node Displacements in Y Dir. vs Time for Truss Element During Dynamic External Time-dependent Loading Analysis") # TRUSS ELEMENT IN Y DIR.


# PLOT STRUCTURAL PERIOD DURING THE ANALYSIS
plt.figure(0, figsize=(12, 8))
plt.plot(disp_mid_ETDLD, PERIOD_MIN_ETDLD, linewidth=3)
plt.plot(disp_mid_ETDLD, PERIOD_MAX_ETDLD, linewidth=3)
plt.title('Period of Structure During Dynamic External Time-dependent Loading Analysis')
plt.ylabel('Structural Period [s]')
plt.xlabel('Displacement [mm]')
#plt.semilogy()
plt.grid()
plt.legend([f'PERIOD - MIN VALUES: Min: {np.min(PERIOD_MIN_ETDLD):.3f} (s) - Mean: {np.mean(PERIOD_MIN_ETDLD):.3f} (s) - Max: {np.max(PERIOD_MIN_ETDLD):.3f} (s)', 
            f'PERIOD - MAX VALUES:  Min: {np.min(PERIOD_MAX_ETDLD):.3f} (s) - Mean: {np.mean(PERIOD_MAX_ETDLD):.3f} (s) - Max: {np.max(PERIOD_MAX_ETDLD):.3f} (s)',
            ])
plt.show()

S01.PLOT_2D_FRAME_TRUSS(deformed_scale=1.0)
#%%----------------------------------------------------
# FREE-VIBRATION ANALYSIS (DYNAMIC TIME-HISTORY ANALYSIS)
#ELE_TYPE = 'Truss'      # MATERIAL NONLINEARITY
ELE_TYPE = 'corotTruss'  # MATERIAL AND GEOMETRIC NONLINEARITY
MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'
ANAL_TYPE = 'FREE-VIBRATION'
DATA = TRUSS_ONE_ELEMENT(LENGTH, AREA, MAT_TYPE, ELE_TYPE, ANAL_TYPE, TOTAL_MASS)

(time_FV, reaction_FV, disp_mid_FV,
ele_force_FV, ele_stress_FV, ele_strain_FV, ele_di_FV,
node_displacementsX_FV, node_displacementsY_FV,
dispX_FV, dispY_FV,
veloX_FV, veloY_FV,
accX_FV, accY_FV,
stiffness_FV, PERIOD_FV, damping_ratio_FV,
PERIOD_MIN_FV, PERIOD_MAX_FV) = DATA

XDATA = disp_mid_FV
YDATA = reaction_FV
XLABEL = 'Displacement in Middle Span [mm]'
YLABEL = 'Base Reaction [N]'
TITLE = 'Base Reaction and Dispalcement of Structure During Free-vibration Analysis'
COLOR = 'black'
SEMILOGY = False
PLOT(XDATA, YDATA, TITLE, XLABEL, YLABEL, COLOR, SEMILOGY)

PLOT_TIME_HISTORY(time_FV, reaction_FV, disp_mid_FV,
                          dispX_FV, dispY_FV,
                          veloX_FV, veloY_FV,
                          accX_FV, accY_FV)

# PLOT ELEMENTS AXIAL FORCE
YLABEL = 'Element Axial Force [N]'
TITLE = "elements  force in X Dir. vs Time for Truss Element During Free-vibration Analysis"
PLOT_FORCES(time_FV, ele_force_FV, YLABEL, TITLE) # TRUSS ELEMENT - ELEMENTS AXIAL FORCE
# PLOT ELEMENTS STRESS
YLABEL = 'Element Stress [N/mm^2]'  
TITLE = "elements  Stress in X Dir. vs Time for Truss Element During Free-vibration Analysis"
PLOT_FORCES(time_FV, ele_stress_FV, YLABEL, TITLE) # TRUSS ELEMENT - ELEMENTS STRESS
# PLOT ELEMENTS STRAIN
YLABEL = 'Element Strain [mm]'  
TITLE = "elements  Strain in X Dir. vs Time for Truss Element During Free-vibration Analysis"
PLOT_FORCES(time_FV, ele_strain_FV, YLABEL, TITLE) # TRUSS ELEMENT - ELEMENTS AXIAL FORCE
# PLOT ELEMENTS DAMAGE-INDEX
YLABEL = 'Element Damage-index [%]'  
TITLE = "elements Damage-index [%] in X Dir. vs Time for Truss Element During Free-vibration Analysis"
PLOT_FORCES(time_FV, ele_strain_FV, YLABEL, TITLE) # TRUSS ELEMENT - ELEMENTS AXIAL FORCE

# PLOT NODAL DISPALEMENTS
PLOT_DISPLAEMENTS(time_FV, node_displacementsX_FV, TITLE = "Node Displacements in X Dir. vs Time for TRUSS ELEMENT During Free-vibration Analysis") # TRUSS ELEMENT IN X DIR.
PLOT_DISPLAEMENTS(time_FV, node_displacementsY_FV, TITLE = "Node Displacements in Y Dir. vs Time for TRUSS ELEMENT During Free-vibration Analysis") # TRUSS ELEMENT IN Y DIR.

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

S01.PLOT_2D_FRAME_TRUSS(deformed_scale=1.0)
#%%----------------------------------------------------
# SEISMIC ANALYSIS (DYNAMIC TIME-HISTORY ANALYSIS)
#ELE_TYPE = 'Truss'      # MATERIAL NONLINEARITY
ELE_TYPE = 'corotTruss'  # MATERIAL AND GEOMETRIC NONLINEARITY
MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'
ANAL_TYPE = 'SEISMIC'

DATA = TRUSS_ONE_ELEMENT(LENGTH, AREA, MAT_TYPE, ELE_TYPE, ANAL_TYPE, TOTAL_MASS)

(time_SEI, reaction_SEI, disp_mid_SEI,
 ele_axialforce_SEI, ele_stress_SEI, ele_strain_SEI, ele_di_SEI,
 node_displacementsX_SEI, node_displacementsY_SEI,
 dispX_SEI, dispY_SEI,
 veloX_SEI, veloY_SEI,
 accX_SEI, accY_SEI,
 stiffness_SEI, PERIOD_SEI, damping_ratio_SEI,
 PERIOD_MIN_SEI, PERIOD_MAX_SEI) = DATA


XDATA = disp_mid_SEI
YDATA = reaction_SEI
XLABEL = 'Displacement in Middle Span [mm]'
YLABEL = 'Base Reaction [N]'
TITLE = 'Base Reaction and Dispalcement of Structure During Seismic Analysis'
COLOR = 'black'
SEMILOGY = False
PLOT(XDATA, YDATA, TITLE, XLABEL, YLABEL, COLOR, SEMILOGY)

PLOT_TIME_HISTORY(time_SEI, reaction_SEI, disp_mid_SEI,
                          dispX_SEI, dispY_SEI,
                          veloX_SEI, veloY_SEI,
                          accX_SEI, accY_SEI)

# PLOT ELEMENTS AXIAL FORCE
YLABEL = 'Element Axial Force [N]'
TITLE = "elements  force in X Dir. vs Time for Truss Element During Seismic Analysis"
PLOT_FORCES(time_SEI, ele_axialforce_SEI, YLABEL, TITLE) # TRUSS ELEMENT - ELEMENTS AXIAL FORCE
# PLOT ELEMENTS STRESS
YLABEL = 'Element Stress [N/mm^2]'  
TITLE = "elements  Stress in X Dir. vs Time for Truss Element During Seismic Analysis"
PLOT_FORCES(time_SEI, ele_stress_SEI, YLABEL, TITLE) # TRUSS ELEMENT - ELEMENTS STRESS
# PLOT ELEMENTS STRAIN
YLABEL = 'Element Strain [mm]'  
TITLE = "elements  Strain in X Dir. vs Time for Truss Element During Seismic Analysis"
PLOT_FORCES(time_SEI, ele_strain_SEI, YLABEL, TITLE) # TRUSS ELEMENT - ELEMENTS AXIAL FORCE
# PLOT ELEMENTS DAMAGE-INDEX
YLABEL = 'Element Damage-index [%]'  
TITLE = "elements Damage-index [%] in X Dir. vs Time for Truss Element During Seismic Analysis"
PLOT_FORCES(time_SEI, ele_strain_SEI, YLABEL, TITLE) # TRUSS ELEMENT - ELEMENTS AXIAL FORCE

# PLOT NODAL DISPALEMENTS
PLOT_DISPLAEMENTS(time_SEI, node_displacementsX_SEI, TITLE = "Node Displacements in X Dir. vs Time for TRUSS ELEMENT During Seismic Analysis") # TRUSS ELEMENT IN X DIR.
PLOT_DISPLAEMENTS(time_SEI, node_displacementsY_SEI, TITLE = "Node Displacements in Y Dir. vs Time for TRUSS ELEMENT During Seismic Analysis") # TRUSS ELEMENT IN Y DIR.

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

S01.PLOT_2D_FRAME_TRUSS(deformed_scale=1.0)
#%%----------------------------------------------------
# --------------------------------------
#  Plot BaseAxial-Displacement Analysis 
# --------------------------------------
XX = np.abs(disp_mid_PUSH); YY = np.abs(reaction_PUSH); # ABSOLUTE VALUE
SLOPE_NODE = 10

DATA = S07.BILNEAR_CURVE(XX, YY, SLOPE_NODE)
X, Y, Elastic_ST, Plastic_ST, Tangent_ST, Ductility_Rito, Over_Strength_Factor = DATA

XLABEL = 'Displacement in X [mm]'
YLABEL = 'Base-Axial Reaction [N]'
LEGEND01 = 'Curve'
LEGEND02 = 'Bilinear Fitted'
LEGEND03 = 'Undefined'
TITLE = f'Last Data of BaseShear-Displacement Analysis - Ductility Ratio: {X[2]/X[1]:.4f} - Over Strength Factor: {Y[2]/Y[1]:.4f}'
COLOR = 'black'
PLOT_2D(np.abs(disp_mid_PUSH), np.abs(reaction_PUSH), X, Y, X, Y, XLABEL, YLABEL, TITLE, LEGEND01, LEGEND02, LEGEND03, COLOR='black', Z=2) 
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
Dd = np.max(np.abs(disp_mid_SEI))
DIy = (Dd - X[1]) /(X[2] - X[1])
print(f'Structural Ductility Damage Index in X Direction: {100*DIy:.4f} (%)')
#%%----------------------------------------------------
# EXCEL OUTPUT
import pandas as pd

# Create DataFrame function
def create_df(reaction, disp_mid, dispX, dispY, PERIOD_MIN, PERIOD_MAX):
    df = pd.DataFrame({
        "reaction": reaction,
        "disp_mid": disp_mid,
        "dispX": dispX,
        "dispY": dispY,
        "PERIOD_MIN": PERIOD_MIN,
        "PERIOD_MAX": PERIOD_MAX,        
    })
    return df


# Save to Excel
with pd.ExcelWriter("TRUSS_ONE_ELEMENT_OUTPUT.xlsx", engine='openpyxl') as writer:
    
    # PUSHOVER
    df1 = create_df(reaction_PUSH, disp_mid_PUSH, dispX_PUSH, dispY_PUSH, PERIOD_MIN_PUSH, PERIOD_MAX_PUSH)
    df1.to_excel(writer, sheet_name="PUSHOVER", index=False)
                 
    # CYCLIC DISPLACEMENT
    df1 = create_df(reaction_CP, disp_mid_CP, dispX_CP, dispY_CP, PERIOD_MIN_CP, PERIOD_MAX_CP)
    df1.to_excel(writer, sheet_name="CYCLIC_DISPLACEMENT", index=False)
    
    # STATIC EXTERNAL TIME-DEPENDENT LOADING
    df2 = create_df(reaction_ETDLS, disp_mid_ETDLS, dispX_ETDLS, dispY_ETDLS, PERIOD_MIN_ETDLS, PERIOD_MAX_ETDLS)
    df2.to_excel(writer, sheet_name="STATIC_EXTERNAL_TIME-DEPENDENT_LOADING", index=False)

    # DYNAMIC EXTERNAL TIME-DEPENDENT LOADING
    df3 = create_df(reaction_ETDLD, disp_mid_ETDLD, dispX_ETDLD, dispY_ETDLD, PERIOD_MIN_ETDLD, PERIOD_MAX_ETDLD)
    df3.to_excel(writer, sheet_name="DYNAMIC_EXTERNAL_TIME-DEPENDENT_LOADING", index=False)
    
    # FREE-VIBRATION
    df3 = create_df(reaction_FV, disp_mid_FV, dispX_FV, dispY_FV, PERIOD_MIN_FV, PERIOD_MAX_FV)
    df3.to_excel(writer, sheet_name="FREE-VIBRATION", index=False)

    # SEISMIC
    df4 = create_df(reaction_SEI, disp_mid_SEI, dispX_SEI, dispY_SEI, PERIOD_MIN_SEI, PERIOD_MAX_SEI)
    df4.to_excel(writer, sheet_name="SEISMIC", index=False)
#%%----------------------------------------------------