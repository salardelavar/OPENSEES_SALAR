###########################################################################################################
#                   >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                      #
# PUSHOVER ANALYSIS OF A MULTI-DEGREE-OF-FREEDOM STRUCTURE VIA EVALUATION OF A MULTILINEAR FITTING CURVE  #
#---------------------------------------------------------------------------------------------------------#
# EQUIVALENT SDOF SYSTEM DERIVATION VIA DISPLACEMENT-BASED SEISMIC DESIGN PROCEDURE WITH PUSHOVER ANALYSIS#
#---------------------------------------------------------------------------------------------------------#
#                  THIS PYTHON SCRIPT IS WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                    #
#                                   EMAIL: salar.d.ghashghaei@gmail.com                                   #
###########################################################################################################
"""
This Python script performs nonlinear pushover analysis of a 10-story MDOF
 building frame using OpenSeesPy to derive an equivalent SDOF system for
 displacement-based seismic design.
1. It models a 10-DOF lumped-mass shear building with 4 parallel zeroLength 
springs per floor, each combining a hysteretic uniaxial material (elastic or inelastic)
 and a viscous damper.
 
2. Spring properties (yield force, ultimate force, elastic stiffness, ultimate displacement) 
are defined per column, from which yield displacement, hardening ratio, and damping coefficient are computed.

3. Multiple analysis types are supported: STATIC, PUSHOVER, CYCLIC_DISPLACEMENT, FREE-VIBRATION,
 SEISMIC, and external time-dependent loading (static or dynamic).

4. In the PUSHOVER case, displacement control is applied at the top node with small increments
 up to twice the ultimate displacement, recording base reaction, nodal displacement, element
 forces, and damage index at each step.

5. At every step, eigenvalue analysis (via `EIGENVALUE_ANALYSIS_FUN`) recomputes the structure's
 period (min/max) to track stiffness degradation.

6. The damage index is computed using a ductility-based function (`DAMAGE_INDEX_FUN`)
 comparing current displacement to yield and ultimate displacements.

7. Results are visualized: base shear vs. displacement, period evolution, element forces,
 nodal displacements, and damage index.

8. A `PLOT_1D_SPRING` routine optionally renders a virtual 1D spring representation of the
 deformed structure.

9. The pushover curve is then fitted with multilinear idealizations — bilinear, trilinear,
 quadrilinear, and pentalinear — using an area-preserving piecewise-linear fitting algorithm
 (`MULTILINEAR_CURVE`).

10. Fitting outputs include elastic stiffness, plastic stiffness, tangent stiffness, ductility ratio,
 and over-strength factor, plotted alongside the original pushover curve for comparison.
 
------------------------------------------------ 
EQUIVALENT SDOF SYSTEM DERIVATION VIA DISPLACEMENT-BASED PUSHOVER ANALYSIS:
Change MDOF to SDOF System with Displacement Based Design Concept

A displacement-based pushover transformation,
 converting a multi-degree-of-freedom (MDOF) system into an equivalent
 single-degree-of-freedom (SDOF) system for seismic assessment.
 It calculates effective modal properties—displacement, mass, and
 stiffness—by weighting element forces and nodal displacements according
 to a presumed deformed shape.
 The derived effective period provides a simplified dynamic characteristic
 for performance-based engineering. The visualizations effectively track the
 evolution of these equivalent parameters throughout the nonlinear analysis steps. 
"""
# BOOK: Dynamics of Structures in SI Units -  Anil Kumar Chopra
'https://share.google/TDN5O4eWmw5pH8zUH'
# BOOK: Differential Equations for Engineers-Wei-Chau Xie-CAMBRIDGE-2010
# BOOK: Structural dynamics THEORY AND COMPUTATION MARIO BAZ-5 EDITION
# PAPER: Displacement-based seismic design of buildings theory - M.S. Medhekar, D.J.L. Kennedy - 1999 Elsevier
# YOUTUBE:
'https://www.youtube.com/watch?v=MZUhSHmIUdI'
#%%----------------------------------------------------
import numpy as np
import openseespy.opensees as ops
import matplotlib.pyplot as plt 
import PLOT_1D_SPRING as S01
import ANALYSIS_FUNCTION as S02
import PERIOD_FUN as S03
import DAMPING_RATIO_FUN as S04
import EIGENVALUE_ANALYSIS_FUN as S05
import RAYLEIGH_DAMPING_FUN as S06
import BILINEAR_CURVE as S07
import OPENSEEES_HYSTERETICSM_FORCE_DISP_FUN as S08
import SALAR_MATH as S09
import EQULIVALENT_VISCOUS_DAMPING_RATIO_FUN as S055
import DAMAGE_INDEX_FUN as S066
#%%----------------------------------------------------
def MDOF(MAT_TYPE, TOTAL_MASS, ANAL_TYPE):
    # Initialize model
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)
    
    MAX_ITERATIONS = 5000   # Maximum number of iterations
    MAX_TOLERANCE = 1.0e-6  # Specified tolerance for convergence
    GMfact = 9.81           # [m/s^2] standard acceleration of gravity or standard acceleration 
    NUM_FLOOR = 10          # Number of Floors (if this variable is changed, you must rewrite other parameters)
    NUM_COL = 4             # Number of Columns (if this variable is changed, you must rewrite other parameters)
    
    # EACH COLUMN IN EACH FLOOR PROPERTIES
    # YILED STRENGTH [N]
    FYi = [85000.0,      # COLUMN 01
           95000.0,      # COLUMN 02
           75000.0,      # COLUMN 03
           65000.0]      # COLUMN 04
    
    # ULTIMATE STRENGTH [N]    
    FUi = [1.18 * FYi[0],   # COLUMN 01
           1.12 * FYi[1],   # COLUMN 02
           1.20 * FYi[2],   # COLUMN 03
           1.10 * FYi[3]]   # COLUMN 04
    
    # ELASTIC STIFFNESS [N/m]
    Kei = [4500000.0,       # COLUMN 01
           4100000.0,       # COLUMN 02
           4600000.0,       # COLUMN 03
           4300000.0]       # COLUMN 04
    
    # ULTIMATE DISPLACEMENT [m]
    DSUi = [0.32,           # COLUMN 01
            0.34,           # COLUMN 02
            0.30,           # COLUMN 03
            0.36]           # COLUMN 04
    
    # EACH FLOOR MASS [kg]
    MASS = [0.10 * TOTAL_MASS, # DOF 01
            0.10 * TOTAL_MASS, # DOF 02
            0.10 * TOTAL_MASS, # DOF 03
            0.10 * TOTAL_MASS, # DOF 04
            0.10 * TOTAL_MASS, # DOF 05
            0.10 * TOTAL_MASS, # DOF 06
            0.10 * TOTAL_MASS, # DOF 07
            0.10 * TOTAL_MASS, # DOF 08
            0.10 * TOTAL_MASS, # DOF 09
            0.10 * TOTAL_MASS] # DOF 10
    
    
    # EACH FLOOR DAMPING RATIO
    DRi = [0.05, # DOF 01 
           0.01, # DOF 02
           0.02, # DOF 03
           0.03, # DOF 04
           0.01, # DOF 05
           0.02, # DOF 06
           0.03, # DOF 07 
           0.01, # DOF 08
           0.02, # DOF 09
           0.03] # DOF 10   
    
    # Define nodes
    ops.node(1, 0.0)  # Fixed base
        
    # Define boundary conditions
    ops.fix(1, 1)

    # Define masses
    for JJ in range(1, NUM_FLOOR+1):
        ops.node(JJ+1, 0.0)
        ops.mass(JJ, MASS[JJ-1])
    
    #%% FOUR COLUMNS AND FOUR DEGREES OF FREEDOM STRUCTURE, CALCULATE LATERAL STIFFNESS AND DAMPING
    ele_tags = []
    for II in range(0, NUM_FLOOR):     # FOR EACH FLOOR 
        for JJ in range(0, NUM_COL):   # FOR EACH COLUMN
            FY = FYi[JJ]                                     # [N] Yield Force of Structure
            FU = FUi[JJ]                                     # [N] Ultimate Force of Structure
            Ke = Kei[JJ]                                     # [N/m] Spring Elastic Stiffness
            DY = FY / Ke                                     # [m] Yield Displacement
            DSU = DSUi[JJ]                                   # [m] Ultimate Displacement
            Ksh = (FU - FY) / (DSU - DY)                     # [N/m] Displacement Hardening Modulus
            Kp = FU / DSU                                    # [N/m] Spring Plastic Stiffness
            b = Ksh / Ke                                     # Displacement Hardening Ratio
            """
            # Positive branch points
            pos_disp = [0, DY, DSU, 1.1*DSU, 1.25*DSU]
            pos_force = [0, FY, FU, 0.2*FU, 0.1*FU]
            KP = np.array([FY, DY, FU, DSU, 0.2*FU, 1.1*DSU, 0.1*FU, 1.25*DSU])
            
            # Negative branch points
            neg_disp = [0, -DY, -DSU, -1.1*DSU, -1.25*DSU]
            neg_force = [0, -FY, -FU, -0.2*FU, -0.1*FU]
            KN = np.array([-FY, -DY, -FU, -DSU, -0.2*FU, -1.1*DSU, -0.1*FU, -1.25*DSU])
    
            # Plot
            plt.figure(0, figsize=(12, 10))
            plt.plot(pos_disp, pos_force, marker='o', color='red')
            plt.plot(neg_disp, neg_force, marker='o', color='black')
            
            plt.xlabel("Displacement [m]")
            plt.ylabel("Force [N]")
            plt.title(f"Force–Displacement Curve for Element {II+1}")
            plt.grid(True)
            plt.axhline(0, linewidth=0.5)
            plt.axvline(0, linewidth=0.5)
            plt.show()
            """
            # Define material properties
            MAT_TAG = 1000 + (II * 4 + JJ) # SPRING TAG
            # FORCE-DISPLACEMENT RELATIONSHIP OF LATERAL SPRING AND PLOT 
            DP = [0, 0, 0, 0]
            FP = [0, 0, 0, 0]
            DN = [0, 0, 0, 0]
            FN = [0, 0, 0, 0]
            #DSU = DY * duct # IN EACH STEP IT WILL BE CHNAGED
            #FU = FY * osf   # IN EACH STEP IT WILL BE CHNAGED
            #print(DSU,"------------" ,FU)
            DP[0], FP[0] = DY, FY
            DP[1], FP[1] = DSU, FU 
            DP[2], FP[2] = 1.1*DSU, 0.20*FU
            DP[3], FP[3] = 1.25*DSU, 0.10*FU
            DN[0], FN[0] = -DY, -FY 
            DN[1], FN[1] = -DSU, -FU
            DN[2], FN[2] = -1.1*DSU, -0.20*FU   
            DN[3], FN[3] = -1.25*DSU, -0.10*FU
            #print(DP, FP)
            #print(DN, FN)
            
            if MAT_TYPE == 'INELASTIC':
                S08.OPENSEEES_HYSTERETICSM_FORCE_DISP_FUN(MAT_TAG, DP, FP, DN, FN, PLOT = False, X_LABEL='Displacement (m)', Y_LABEL='Force [N]', TITLE='FORCE-DISPLACEMENT CURVE')
            if MAT_TYPE == 'ELASTIC':
                #ops.uniaxialMaterial('Elastic', MAT_TAG, Ke)             # TESNSION AND COMPRESSION IS SAME VALUES
                ops.uniaxialMaterial('Elastic', MAT_TAG, Ke ,0.0, 0.5*Ke) # TESNSION AND COMPRESSION IS NOT SAME VALUES
                # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/ElasticUni.html
            
            MAT_TAG_C = 2000 + (II * 4 + JJ) # SPRING DAMPER
            print(MAT_TAG, MAT_TAG_C)
            alpha = 1.0    # velocity exponent (usually 0.3–1.0)
            omega = np.sqrt(Ke / MASS[II])
            Cd = 2 * DRi[II] * omega * MASS[II]  # [N·s/m] Damping coefficient 
            ops.uniaxialMaterial('Viscous', MAT_TAG_C, Cd, alpha)  # Material for C (alpha=1.0 for linear)
            
            # Define element
            eleTag = (II * NUM_FLOOR) + (JJ + 1); ele_tags.append(eleTag);
            ops.element('zeroLength', eleTag, II+1, II+2, '-mat',  MAT_TAG, MAT_TAG_C, '-dir', 1, 1)  # DOF LATERAL SPRING

    center_node = NUM_FLOOR + 1
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    for JJ in range(1, NUM_FLOOR+1):
        ops.load(JJ+1, 1.0)

    
    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.system('BandGeneral')
    
    if MAT_TYPE == 'ELASTIC':
        ops.algorithm('Linear')
    if MAT_TYPE == 'INELASTIC':
        ops.algorithm('Newton')
        
    time = []
    disp = []
    velo = []
    acc = []
    reaction = []
    stiffness = []
    OMEGA, PERIOD = [], []
    PERIOD_MIN, PERIOD_MAX = [], []
    DI = []
    
    # Initialize lists for each element's force
    ele_force = {tag: [] for tag in ele_tags}
    
    # Initialize lists for each node's displacement
    node_displacements = {ZZ: [] for ZZ in range(1, NUM_FLOOR+2)}
    
    if ANAL_TYPE == 'PERIOD': 
      #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(5, 0.5*DR, DR, 0, 1)
      PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(5, PLOT=True) 
      # Compute modal properties
      ops.modalProperties("-print", "-file", f"SALAR_ModalReport_{ANAL_TYPE}.txt", "-unorm")
      return PERIODmin, PERIODmax
  
    if ANAL_TYPE == 'STATIC': 
        ops.integrator('LoadControl', 1.0)
        ops.analysis('Static')
        OK = ops.analyze(1)
        S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS)
        ops.reactions()
        reaction.append(ops.nodeReaction(1, 1))                         # BASE REACTION
        disp.append(ops.nodeDisp(center_node, 1))                       # DISPLACEMENT NODE 11 IN X DIR                
        # EVALUATION OF DUCTILITY DAMAGE INDEX
        if MAT_TYPE == 'INELASTIC':
            di = S066.DAMAGE_INDEX_FUN(disp[-1], DY, DSU)
            DI.append(di)                   # DAMAGE INDEX
        if MAT_TYPE == 'ELASTIC':
            DI.append(0.0)
        print('\n\nSTATIC ANALYSIS DONE.\n\n')  
            
        DATA = (reaction, disp, DI)
        
        return  DATA
    
    if ANAL_TYPE == 'PUSHOVER': # STATIC TIME-HISTORY ANALYSIS
        IDctrlDOF = 1   # 1: Horizental Dispalcement - 2: Vertical Dispalcement
        DINCR = -0.001  # [m] Incremental Vertical Displacement
        DMAX = -2.0*DSU     # [m] Max. Displacement
        ops.integrator('DisplacementControl', center_node, IDctrlDOF, DINCR)
        ops.analysis('Static')
        Nsteps =  int(np.abs(DMAX/ DINCR)) 
        STEP = 0.0
        for step in range(Nsteps):
            OK = ops.analyze(1)
            S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS)
            ops.reactions()
            reaction.append(ops.nodeReaction(1, 1))                         # BASE REACTION
            disp.append(ops.nodeDisp(center_node, 1))                      # DISPLACEMENT NODE 11 IN X DIR       
            # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
            #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(5, 0.5*DR, DR, 0, 1)
            PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(5, PLOT=True)
            PERIOD_MIN.append(PERIODmin)
            PERIOD_MAX.append(PERIODmax) 
            # EVALUATION OF DUCTILITY DAMAGE INDEX
            if MAT_TYPE == 'INELASTIC':
                di = S066.DAMAGE_INDEX_FUN(disp[-1], DY, DSU)
                DI.append(di)                   # DAMAGE INDEX
            if MAT_TYPE == 'ELASTIC':
                DI.append(0.0)
            # Store forces and displacements
            for ele_id in ele_force.keys(): 
                ele_force[ele_id].append(ops.eleResponse(ele_id, 'force')[0])        # [N] ELEMENT AXIAL FORCE       # [N] ELEMENT MOMENT FORCE                
            # Store displacements
            for node_id in node_displacements.keys():    
                node_displacements[node_id].append(ops.nodeDisp(node_id, 1))
            STEP += 1
            #print(STEP, disp[-1], reaction[-1])
            print(f"Step: {STEP}, Displacement: {disp[-1]:.4f} m, Reaction: {reaction[-1]:.2f} N")     
        else:
            print('\n\nPUSHOVER ANALYSIS DONE.\n\n')
            
        DATA = (reaction, disp, DI,
                ele_force, node_displacements,
                np.array(PERIOD_MIN), np.array(PERIOD_MAX))
    
        # Run the file loading effective properties
        exec(open("COMPUTE_EFFECTIVE_PROPERTIES_FUN_PUSHOVER.py").read())
        
        return  DATA
    
    if ANAL_TYPE == 'CYCLIC_DISPLACEMENT': # STATIC TIME-HISTORY ANALYSIS
        IDctrlDOF = 1   # 1: Horizental Dispalcement - 2: Vertical Dispalcement
        DMAX = -1.0*DSU     # [m] Max. Displacement
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
            reaction.append(ops.nodeReaction(1, 1))                         # BASE REACTION
            disp.append(ops.nodeDisp(center_node, 1))                       # DISPLACEMENT NODE 11 IN X DIR       
            # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
            #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(5, 0.5*DR, DR, 0, 1)
            PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(5, PLOT=True)
            PERIOD_MIN.append(PERIODmin)
            PERIOD_MAX.append(PERIODmax)
            # EVALUATION OF DUCTILITY DAMAGE INDEX
            if MAT_TYPE == 'INELASTIC':
                di = S066.DAMAGE_INDEX_FUN(disp[-1], DY, DSU)
                DI.append(di)                   # DAMAGE INDEX
            if MAT_TYPE == 'ELASTIC':
                DI.append(0.0)
            # Store forces and displacements
            for ele_id in ele_force.keys(): 
                ele_force[ele_id].append(ops.eleResponse(ele_id, 'force')[0])        # [N] ELEMENT AXIAL FORCE       # [N] ELEMENT MOMENT FORCE                
            # Store displacements
            for node_id in node_displacements.keys():    
                node_displacements[node_id].append(ops.nodeDisp(node_id, 1))
            STEP += 1
            #print(STEP, disp[-1], reaction[-1])
            print(f"Step: {STEP}, Displacement: {disp[-1]:.4f} m, Reaction: {reaction[-1]:.2f} N")     
        else:
            print('\n\nCYCLIC DISPLAEMENT ANALYSIS DONE.\n\n')
            
        DATA = (reaction, disp, DI,
                ele_force, node_displacements,
                np.array(PERIOD_MIN), np.array(PERIOD_MAX))
    
        return  DATA 

    if ANAL_TYPE == 'STATIC_EXTERNAL_TIME-DEPENDENT_LOADING': # STATIC TIME-HISTORY ANALYSIS
        #%% DEFINE EXTERNAL TIME-DEPENDENT LOADING PROPERTIES
        # IN HERE ANALYSIS TIME AND DURATION ARE LOAD STEPS
        duration = 20.0             # [s] Analysis duration
        dt = 0.01                   # [s] Time step
        DT = dt                     # [s] Time step
        DT_time = 5.0               # [s] Total external Load Analysis Durations [*******]
        force_amplitude = 5960.0    # [N] Amplitude Force
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
        ops.load(center_node, 1.0)
        
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
            reaction.append(ops.nodeReaction(1, 1))                         # BASE REACTION
            disp.append(ops.nodeDisp(center_node, 1))                       # DISPLACEMENT NODE 11 IN X DIR      
            # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
            #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(5, 0.5*DR, DR, 0, 1)
            PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(5, PLOT=True)
            PERIOD_MIN.append(PERIODmin)
            PERIOD_MAX.append(PERIODmax)
            # EVALUATION OF DUCTILITY DAMAGE INDEX
            if MAT_TYPE == 'INELASTIC':
                di = S066.DAMAGE_INDEX_FUN(disp[-1], DY, DSU)
                DI.append(di)                   # DAMAGE INDEX
            if MAT_TYPE == 'ELASTIC':
                DI.append(0.0)
            # Store forces and displacements
            for ele_id in ele_force.keys(): 
                ele_force[ele_id].append(ops.eleResponse(ele_id, 'force')[0])        # [N] ELEMENT AXIAL FORCE       # [N] ELEMENT MOMENT FORCE                
            # Store displacements
            for node_id in node_displacements.keys():    
                node_displacements[node_id].append(ops.nodeDisp(node_id, 1))
            STEP += 1
            #print(STEP, disp[-1], reaction[-1])
            print(f"Step: {STEP}, Displacement: {disp[-1]:.4f} m, Reaction: {reaction[-1]:.2f} N")     
        else:
            print('\n\nSTATIC EXTERNAL TIME-DEPENDENT LOADING ANALYSIS DONE.\n\n')
            
        DATA = (reaction, disp, DI,
                ele_force, node_displacements,
                np.array(PERIOD_MIN), np.array(PERIOD_MAX))
    
        return  DATA
    
    if ANAL_TYPE == 'DYNAMIC_EXTERNAL_TIME-DEPENDENT_LOADING': # DYNAMIC TIME-HISTORY ANALYSIS
        #%% DEFINE EXTERNAL TIME-DEPENDENT LOADING PROPERTIES
        duration = 20.0             # [s] Analysis duration
        dt = 0.01                   # [s] Time step
        DT = dt                     # [s] Time step
        DT_time = 5.0               # [s] Total external Load Analysis Durations [*******]
        force_amplitude = 5960.0    # [N] Amplitude Force
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
        ops.load(center_node, 1.0)
        
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
            reaction.append(ops.nodeReaction(1, 1))               # BASE REACTION
            disp.append(ops.nodeDisp(center_node, 1))             # DISPLACEMENT NODE 11 IN X DIR 
            velo.append(ops.nodeVel(center_node, 1))              # VELOCITY NODE 11
            acc.append(ops.nodeAccel(center_node, 1))             # ACCELERATION NODE 11
            stiffness.append(np.abs(reaction[-1] / disp[-1]))
            OMEGA.append(np.sqrt(stiffness[-1]/TOTAL_MASS))
            PERIOD.append((np.pi * 2) / OMEGA[-1])   
            # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
            #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(5, 0.5*DR, DR, 0, 1)
            PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(5, PLOT=True)
            PERIOD_MIN.append(PERIODmin)
            PERIOD_MAX.append(PERIODmax)
            # EVALUATION OF DUCTILITY DAMAGE INDEX
            if MAT_TYPE == 'INELASTIC':
                di = S066.DAMAGE_INDEX_FUN(disp[-1], DY, DSU)
                DI.append(di)                   # DAMAGE INDEX
            if MAT_TYPE == 'ELASTIC':
                DI.append(0.0)
            # Store forces and displacements
            for ele_id in ele_force.keys(): 
                ele_force[ele_id].append(ops.eleResponse(ele_id, 'force')[0])        # [N] ELEMENT AXIAL FORCE            
            # Store displacements
            for node_id in node_displacements.keys():    
                node_displacements[node_id].append(ops.nodeDisp(node_id, 1))
            #print(time[-1], disp[-1], velo[-1])
            print(f"Time: {time[-1]:.4f}, Displacement: {disp[-1]:.4f} m, Reaction: {reaction[-1]:.2f} N")      
        else:
            print('\n\nDYNAMIC EXTERNAL TIME-DEPENDENT LOADING ANALYSIS DONE.\n\n')  
        # Calculating Damping Ratio and Period Using Logarithmic Decrement Analysis 
        damping_ratio = S04.DAMPING_RATIO(disp)  
        
        # Compute modal properties
        ops.modalProperties("-print", "-file", f"SALAR_ModalReport_{ANAL_TYPE}.txt", "-unorm") 
        
        DATA = (time, reaction, disp, velo, acc, DI,
                ele_force, node_displacements,
                stiffness, PERIOD, damping_ratio,
                np.array(PERIOD_MIN), np.array(PERIOD_MAX))
        
        return  DATA
        
    if ANAL_TYPE == 'FREE-VIBRATION': # DYNAMIC TIME-HISTORY ANALYSIS
        #%% DEFINE PARAMETERS FOR FREE-VIBRATION ANALYSIS
        u0 = -0.010                       # [m] Initial displacement
        v0 = 0.0015                       # [m/s] Initial velocity
        a0 = 0.0065                       # [m/s^2] Initial acceleration
        IU = True                          # Free Vibration with Initial Displacement
        IV = True                          # Free Vibration with Initial Velocity
        IA = True                          # Free Vibration with Initial Acceleration
        duration = 20.0                    # [s] Analysis duration
        dt = 0.001                         # [s] Time step
        DR = 0.03                          # Damping Ratio
        
        disp_02, disp_03, disp_04, disp_05 = [], [], [], []
        disp_06, disp_07, disp_08, disp_09, disp_10, disp_11 = [], [], [], [], [], []
        velo_02, velo_03, velo_04, velo_05 = [], [], [], []
        velo_06, velo_07, velo_08, velo_09, velo_10, velo_11 = [], [], [], [], [], []
        accel_02, accel_03, accel_04, accel_05 = [], [], [], []
        accel_06, accel_07, accel_08, accel_09, accel_10, accel_11 = [], [], [], [], [], []
        
        if IU == True:
            # Define initial displacment
            for JJ in range(1, NUM_FLOOR+1):
                ops.setNodeDisp(JJ+1, 1, u0, '-commit')

            # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/setNodeDisp.html
        if IV == True:
            # Define initial velocity
            for JJ in range(1, NUM_FLOOR+1):
                ops.setNodeVel(JJ+1, 1, v0, '-commit')

            # INFO LINK: https://openseespydoc.readthedocs.io/en/stable/src/setNodeVel.html
        if IA == True:
            # Define initial  acceleration
            for JJ in range(1, NUM_FLOOR+1):
                ops.setNodeAccel(JJ+1, 1, a0, '-commit')

            # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/setNodeAccel.html
            
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
            reaction.append(ops.nodeReaction(1, 1))               # BASE REACTION
            disp.append(ops.nodeDisp(center_node, 1))             # DISPLACEMENT NODE 11 IN X DIR 
            velo.append(ops.nodeVel(center_node, 1))              # VELOCITY NODE 11
            acc.append(ops.nodeAccel(center_node, 1))             # ACCELERATION NODE 11
            disp_02.append(ops.nodeDisp(2, 1)); velo_02.append(ops.nodeVel(2, 1)); accel_02.append(ops.nodeAccel(2, 1));
            disp_03.append(ops.nodeDisp(3, 1)); velo_03.append(ops.nodeVel(3, 1)); accel_03.append(ops.nodeAccel(3, 1));
            disp_04.append(ops.nodeDisp(4, 1)); velo_04.append(ops.nodeVel(4, 1)); accel_04.append(ops.nodeAccel(4, 1));
            disp_05.append(ops.nodeDisp(5, 1)); velo_05.append(ops.nodeVel(5, 1)); accel_05.append(ops.nodeAccel(5, 1));
            disp_06.append(ops.nodeDisp(6, 1)); velo_06.append(ops.nodeVel(6, 1)); accel_06.append(ops.nodeAccel(6, 1));
            disp_07.append(ops.nodeDisp(7, 1)); velo_07.append(ops.nodeVel(7, 1)); accel_07.append(ops.nodeAccel(7, 1));
            disp_08.append(ops.nodeDisp(8, 1)); velo_08.append(ops.nodeVel(8, 1)); accel_08.append(ops.nodeAccel(8, 1));
            disp_09.append(ops.nodeDisp(9, 1)); velo_09.append(ops.nodeVel(9, 1)); accel_09.append(ops.nodeAccel(9, 1));
            disp_10.append(ops.nodeDisp(10, 1)); velo_10.append(ops.nodeVel(10, 1)); accel_10.append(ops.nodeAccel(10, 1));
            disp_11.append(ops.nodeDisp(11, 1)); velo_11.append(ops.nodeVel(11, 1)); accel_11.append(ops.nodeAccel(11, 1));
            stiffness.append(np.abs(reaction[-1] / disp[-1]))
            OMEGA.append(np.sqrt(stiffness[-1]/TOTAL_MASS))
            PERIOD.append((np.pi * 2) / OMEGA[-1]) 
            # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
            #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(5, 0.5*DR, DR, 0, 1)
            PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(5, PLOT=True)
            PERIOD_MIN.append(PERIODmin)
            PERIOD_MAX.append(PERIODmax)
            # EVALUATION OF DUCTILITY DAMAGE INDEX
            if MAT_TYPE == 'INELASTIC':
                di = S066.DAMAGE_INDEX_FUN(disp[-1], DY, DSU)
                DI.append(di)                   # DAMAGE INDEX
            if MAT_TYPE == 'ELASTIC':
                DI.append(0.0)
            # Store forces and displacements
            for ele_id in ele_force.keys(): 
                ele_force[ele_id].append(ops.eleResponse(ele_id, 'force')[0])        # [N] ELEMENT AXIAL FORCE             
            # Store displacements
            for node_id in node_displacements.keys():    
                node_displacements[node_id].append(ops.nodeDisp(node_id, 1))
            #print(time[-1], disp[-1], velo[-1])
            print(f"Time: {time[-1]:.4f}, Displacement: {disp[-1]:.4f} m, Reaction: {reaction[-1]:.2f} N")      
        else:
            print('\n\nFREE-VIBRATION ANALYSIS DONE.\n\n')    
        # Calculating Damping Ratio and Period Using Logarithmic Decrement Analysis  
        damping_ratio = S04.DAMPING_RATIO(disp)       # DAMAPING RATIO FROM DOF 05
        
        damping_ratio_02 = S04.DAMPING_RATIO(disp_02) # DAMAPING RATIO FROM DOF 02
        damping_ratio_03 = S04.DAMPING_RATIO(disp_03) # DAMAPING RATIO FROM DOF 03
        damping_ratio_04 = S04.DAMPING_RATIO(disp_04) # DAMAPING RATIO FROM DOF 04
        damping_ratio_05 = S04.DAMPING_RATIO(disp_05) # DAMAPING RATIO FROM DOF 05
        damping_ratio_06 = S04.DAMPING_RATIO(disp_06) # DAMAPING RATIO FROM DOF 06
        damping_ratio_07 = S04.DAMPING_RATIO(disp_07) # DAMAPING RATIO FROM DOF 07
        damping_ratio_08 = S04.DAMPING_RATIO(disp_08) # DAMAPING RATIO FROM DOF 08
        damping_ratio_09 = S04.DAMPING_RATIO(disp_09) # DAMAPING RATIO FROM DOF 09
        damping_ratio_10 = S04.DAMPING_RATIO(disp_10) # DAMAPING RATIO FROM DOF 10
        damping_ratio_11 = S04.DAMPING_RATIO(disp_11) # DAMAPING RATIO FROM DOF 11
        #print(damping_ratio_02, damping_ratio_03, damping_ratio_04, damping_ratio_05)
        
        # Compute modal properties
        ops.modalProperties("-print", "-file", f"SALAR_ModalReport_{ANAL_TYPE}.txt", "-unorm") 
        
        # Run the file loading effective properties
        exec(open("COMPUTE_EFFECTIVE_PROPERTIES_FUN_FREE_VIBRATION.py").read())
        
        DATA = (time, reaction, disp, velo, acc, DI,
                ele_force, node_displacements,
                stiffness, PERIOD, damping_ratio,
                np.array(PERIOD_MIN), np.array(PERIOD_MAX))
        
        return DATA
    if ANAL_TYPE == 'SEISMIC': # DYNAMIC TIME-HISTORY ANALYSIS
        #%% DEFINE PARAMETERS FOR SEISMIC ANALYSIS
        duration = 20.0                    # [s] Analysis duration
        dt = 0.01                          # [s] Time step
        GMfact = 9.810                     # [m/s²]standard acceleration of gravity or standard acceleration
        SSF_X = 5.5                        # Seismic Acceleration Scale Factor in X Direction
        SSF_Y = 5.5                        # Seismic Acceleration Scale Factor in Y Direction
        iv0_X = 0.0005                     # [m/s] Initial velocity applied to the node  in X Direction
        iv0_Y = 0.0005                     # [m/s] Initial velocity applied to the node  in Y Direction
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
            reaction.append(ops.nodeReaction(1, 1))               # BASE REACTION
            disp.append(ops.nodeDisp(center_node, 1))             # DISPLACEMENT NODE 11 IN X DIR 
            velo.append(ops.nodeVel(center_node, 1))              # VELOCITY NODE 11
            acc.append(ops.nodeAccel(center_node, 1))             # ACCELERATION NODE 11
            stiffness.append(np.abs(reaction[-1] / disp[-1]))
            OMEGA.append(np.sqrt(stiffness[-1]/TOTAL_MASS))
            PERIOD.append((np.pi * 2) / OMEGA[-1])
            # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
            #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(5, 0.5*DR, DR, 0, 1)
            PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(5, PLOT=True)
            PERIOD_MIN.append(PERIODmin)
            PERIOD_MAX.append(PERIODmax)
            # EVALUATION OF DUCTILITY DAMAGE INDEX
            if MAT_TYPE == 'INELASTIC':
                di = S066.DAMAGE_INDEX_FUN(disp[-1], DY, DSU)
                DI.append(di)                   # DAMAGE INDEX
            if MAT_TYPE == 'ELASTIC':
                DI.append(0.0)
            # Store forces and displacements
            for ele_id in ele_force.keys(): 
                ele_force[ele_id].append(ops.eleResponse(ele_id, 'force')[0])        # [N] ELEMENT AXIAL FORCE       # [N] ELEMENT MOMENT FORCE                
            # Store displacements
            for node_id in node_displacements.keys():    
                node_displacements[node_id].append(ops.nodeDisp(node_id, 1))
            #print(time[-1], disp[-1], velo[-1])
            print(f"Time: {time[-1]:.4f}, Displacement: {disp[-1]:.4f} m, Reaction: {reaction[-1]:.2f} N")

        else:
            print('\n\nSEISMIC ANALYSIS DONE.\n\n')     
        # Calculating Damping Ratio and Period Using Logarithmic Decrement Analysis 
        damping_ratio = S04.DAMPING_RATIO(disp)  
        
        # Compute modal properties
        ops.modalProperties("-print", "-file", f"SALAR_ModalReport_{ANAL_TYPE}.txt", "-unorm") 
        
        DATA = (time, reaction, disp, velo, acc, DI,
                ele_force, node_displacements,
                stiffness, PERIOD, damping_ratio,
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
TOTAL_MASS = 500000.0   # [kg] Total Mass of Structure
#%%----------------------------------------------------
# PUSHOVER ANALYSIS (STATIC TIME-HISTORY ANALYSIS)
MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'
ANAL_TYPE = 'PUSHOVER'

DATA = MDOF(MAT_TYPE, TOTAL_MASS, ANAL_TYPE)
(reaction_PUSH, disp_PUSH, DI_PUSH,
 ele_force_PUSH, node_displacements_PUSH,
 PERIOD_MIN_PUSH, PERIOD_MAX_PUSH) = DATA


XDATA = disp_PUSH
YDATA = reaction_PUSH
XLABEL = 'Displacement [m]'
YLABEL = 'Base Reaction [N]'
TITLE = 'Base Reaction and Displacement of Structure During Pushover Analysis'
COLOR = 'black'
SEMILOGY = False
PLOT(XDATA, YDATA, TITLE, XLABEL, YLABEL, COLOR, SEMILOGY)

DATA = S07.BILNEAR_CURVE(np.abs(disp_PUSH), np.abs(reaction_PUSH), SLOPE_NODE=10)
(X_PUSH, Y_PUSH, Elastic_ST, Plastic_ST, Tangent_ST, Ductility_Rito, Over_Strength_Factor) = DATA

# PLOT STRUCTURAL PERIOD DURING THE ANALYSIS
plt.figure(0, figsize=(12, 8))
plt.plot(disp_PUSH, PERIOD_MIN_PUSH, linewidth=3)
plt.plot(disp_PUSH, PERIOD_MAX_PUSH, linewidth=3)
plt.title('Period of Structure During Pushover Analysis')
plt.ylabel('Structural Period [s]')
plt.xlabel('Displacement [m]')
#plt.semilogy()
plt.grid()
plt.legend([f'PERIOD - MIN VALUES: Min: {np.min(PERIOD_MIN_PUSH):.3f} (s) - Mean: {np.mean(PERIOD_MIN_PUSH):.3f} (s) - Max: {np.max(PERIOD_MIN_PUSH):.3f} (s)', 
            f'PERIOD - MAX VALUES:  Min: {np.min(PERIOD_MAX_PUSH):.3f} (s) - Mean: {np.mean(PERIOD_MAX_PUSH):.3f} (s) - Max: {np.max(PERIOD_MAX_PUSH):.3f} (s)',
            ])
plt.show()

# PLOT FORCES
YLABEL = 'Force [N]'
TITLE = "elements  force vs Time During Pushover Analysis"
STOP = len(disp_PUSH)
NUM = STOP
step_PUSH = np.linspace(1, STOP, NUM)
PLOT_FORCES(step_PUSH, ele_force_PUSH, YLABEL, TITLE) # ELEMENTS FORCE

# PLOT NODAL DISPALEMENTS
STOP = len(disp_PUSH)
NUM = STOP
step_PUSH = np.linspace(1, STOP, NUM)
PLOT_DISPLAEMENTS(step_PUSH, node_displacements_PUSH, TITLE = "Node Displacements vs Time for During Pushover Analysis") 

plt.figure(-1, figsize=(12, 8))
plt.plot(disp_PUSH, DI_PUSH, color='black', linewidth=2)
plt.xlabel('Displacement [m]')
plt.ylabel('Structural Damage Index [%]')
plt.title(f'Displacement vs Structural Damage Index')
plt.grid()
plt.show()

S01.PLOT_1D_SPRING(deformed_scale=1.0, virtual_spring_length=1.0, virtual_spring_angle=0.0)
#%%----------------------------------------------------
# --------------------------------------
#  Plot BaseAxial-Displacement Analysis 
# --------------------------------------

#%% Core fitting functions
def _area_under_curve(Cur, Mom):
    Cur = np.asarray(Cur, dtype=float)
    Mom = np.asarray(Mom, dtype=float)
    return np.sum((Mom[:-1] + Mom[1:]) * 0.5 * np.diff(Cur))


def _piecewise_linear_fit(Cur, Mom, fixed_nodes, SLOPE_NODE):
    """
    Fit a piecewise linear curve while preserving the area under the original curve.

    fixed_nodes : indices of intermediate breakpoints after yield.
                  Example: [] for bilinear, [8] for trilinear,
                           [7, 11] for quadrilinear, etc.
    SLOPE_NODE  : index used to compute the initial elastic slope k0.
    """
    Cur = np.asarray(Cur, dtype=float)
    Mom = np.asarray(Mom, dtype=float)

    if len(Cur) != len(Mom):
        raise ValueError("Cur and Mom must have the same length")
    if not (0 <= SLOPE_NODE < len(Cur)):
        raise ValueError("SLOPE_NODE out of range")
    if Cur[SLOPE_NODE] == 0:
        raise ValueError("Cur[SLOPE_NODE] must be non-zero")

    k0 = Mom[SLOPE_NODE] / Cur[SLOPE_NODE]
    Area = _area_under_curve(Cur, Mom)

    Cu = Cur[-1]
    Mu = Mom[-1]

    fixed_nodes = list(fixed_nodes)
    for n in fixed_nodes:
        if not (0 <= n < len(Cur) - 1):
            raise ValueError("fixed node must be a valid index before the last point")
    if any(fixed_nodes[i] >= fixed_nodes[i + 1] for i in range(len(fixed_nodes) - 1)):
        raise ValueError("fixed_nodes must be strictly increasing")

    X_fixed = [Cur[i] for i in fixed_nodes] + [Cu]
    Y_fixed = [Mom[i] for i in fixed_nodes] + [Mu]

    x2, y2 = X_fixed[0], Y_fixed[0]

    # Area from the first fixed breakpoint to the end
    A_rest = 0.0
    for j in range(len(X_fixed) - 1):
        A_rest += 0.5 * (Y_fixed[j] + Y_fixed[j + 1]) * (X_fixed[j + 1] - X_fixed[j])

    denom = k0 * x2 - y2
    if abs(denom) < 1e-12:
        raise ValueError("Cannot solve for first yield point: denominator near zero")

    x1 = (2.0 * Area - 2.0 * A_rest - y2 * x2) / denom
    y1 = k0 * x1

    if not (0.0 < x1 < x2):
        print(f"Warning: computed x1 = {x1:.4f} is not between 0 and x2 = {x2:.4f}")

    X = np.concatenate(([0.0, x1], X_fixed))
    Y = np.concatenate(([0.0, y1], Y_fixed))
    return X, Y, k0


def _postprocess(X, Y):
    Elastic_ST = Y[1] / X[1]
    Plastic_ST = Y[-1] / X[-1]
    Tangent_ST = np.diff(Y) / np.diff(X)
    Ductility_Rito = X[-1] / X[1]
    Over_Strength_Factor = Y[-1] / Y[1]

    print('+==========================+')
    print('=   Analysis curve fitted =')
    print('     Disp        Base Shear')
    print('----------------------------')
    print(np.column_stack((X, Y)))
    print('+==========================+')
    print('+----------------------------------------------------+')
    print(f' Structure Elastic Stiffness :     {Elastic_ST:.2f}')
    print(f' Structure Plastic Stiffness :     {Plastic_ST:.2f}')
    print(f' Structure Tangent Stiffness :     {np.array2string(Tangent_ST, precision=2)}')
    print(f' Structure Ductility Ratio :       {Ductility_Rito:.2f}')
    print(f' Structure Over Strength Factor:   {Over_Strength_Factor:.2f}')
    print('+----------------------------------------------------+')

    return X, Y, Elastic_ST, Plastic_ST, Tangent_ST, Ductility_Rito, Over_Strength_Factor


def MULTILINEAR_CURVE(Cur, Mom, SLOPE_NODE, FIXED_NODES):
    X, Y, k0 = _piecewise_linear_fit(Cur, Mom, FIXED_NODES, SLOPE_NODE)
    return _postprocess(X, Y)


def BILINEAR_CURVE(Cur, Mom, SLOPE_NODE):
    return MULTILINEAR_CURVE(Cur, Mom, SLOPE_NODE, [])


def TRILINEAR_CURVE(Cur, Mom, SLOPE_NODE, NODE2):
    return MULTILINEAR_CURVE(Cur, Mom, SLOPE_NODE, [NODE2])


def QUADRILINEAR_CURVE(Cur, Mom, SLOPE_NODE, NODE2, NODE3):
    return MULTILINEAR_CURVE(Cur, Mom, SLOPE_NODE, [NODE2, NODE3])


def PENTALINEAR_CURVE(Cur, Mom, SLOPE_NODE, NODE2, NODE3, NODE4):
    return MULTILINEAR_CURVE(Cur, Mom, SLOPE_NODE, [NODE2, NODE3, NODE4])


SLOPE_NODE = 1

# Choose breakpoint indices from the original curve
# %%
TRI_NODE = int(0.6 * len(disp_PUSH))

QUAD_NODES = (int(0.6 * (len(disp_PUSH))),
              int(0.65 * (len(disp_PUSH))))

PENTA_NODES = (int(0.6 * (len(disp_PUSH))),
               int(0.65 * (len(disp_PUSH))),
               int(0.9 * (len(disp_PUSH))))

#%% Fit curves
disp_PUSH, reaction_PUSH = np.abs(disp_PUSH), np.abs(reaction_PUSH)
X_bi, Y_bi, *_ = BILINEAR_CURVE(disp_PUSH, reaction_PUSH, SLOPE_NODE)
X_tri, Y_tri, *_ = TRILINEAR_CURVE(disp_PUSH, reaction_PUSH, SLOPE_NODE, TRI_NODE)
X_quad, Y_quad, *_ = QUADRILINEAR_CURVE(disp_PUSH, reaction_PUSH, SLOPE_NODE, *QUAD_NODES)
X_penta, Y_penta, *_ = PENTALINEAR_CURVE(disp_PUSH, reaction_PUSH, SLOPE_NODE, *PENTA_NODES)

#%% Plot
plt.figure(figsize=(10, 6))

plt.plot(disp_PUSH, reaction_PUSH, 'ko-', linewidth=2, markersize=5, label='Original curve')

plt.plot(X_bi, Y_bi, 'b--', linewidth=2, marker='s', markersize=5,
         label='Bilinear')

plt.plot(X_tri, Y_tri, 'g-.', linewidth=2, marker='^', markersize=5,
         label='Trilinear')

plt.plot(X_quad, Y_quad, 'm:', linewidth=2, marker='D', markersize=5,
         label='Quadrilinear')

plt.plot(X_penta, Y_penta, 'r-', linewidth=2, marker='o', markersize=5,
         label='Pentalinear')

plt.xlabel('Displacement [m]')
plt.ylabel('Base Shear [N]')
plt.title('Multilinear Fitting of Pushover Curve')
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()

#%%----------------------------------------------------
