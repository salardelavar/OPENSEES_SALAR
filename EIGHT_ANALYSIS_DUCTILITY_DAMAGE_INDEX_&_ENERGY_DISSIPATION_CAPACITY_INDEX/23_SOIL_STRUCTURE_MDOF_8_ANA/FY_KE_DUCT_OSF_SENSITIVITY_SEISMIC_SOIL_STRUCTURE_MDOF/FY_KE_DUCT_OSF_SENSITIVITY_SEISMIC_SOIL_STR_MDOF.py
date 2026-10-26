###########################################################################################################
#                   >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                      #
#            SENSITIVITY ANALYSIS OF SOIL'S DAMPING RATIO, YIELD STRENGTH, ELASTIC STIFFNESS,             #
#                 DUCTILITY RATIO AND OVER-STRENGTH FACTOR ON THE SEISMIC PERFORMANCE OF                  #
#                   SOIL–STRUCTURE MDOF SYSTEMS USING OPENSEES: CYCLIC DISPLACEMENT AND                   #
#                    SEISMIC TIME-HISTORY ANALYSES WITH EVALUATION OF PERIOD ELONGATION,                  #
#                    ENERGY DISSIPATION CAPACITY INDEX, EQUIVALENT VISCOUS DAMPING RATIO,                 #
#                                 DUCTILITY DAMAGE INDEX AND RECONSTRUCTION COST                          #
#---------------------------------------------------------------------------------------------------------#
#                  THIS PYTHON SCRIPT IS WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                    #
#                                   EMAIL: salar.d.ghashghaei@gmail.com                                   #
###########################################################################################################
"""
This OpenSeesPy-based Python framework performs a comprehensive multi-parameter
 sensitivity analysis of a two-degree-of-freedom nonlinear soil–structure system subjected
 to both cyclic displacement protocols and selected ground-motion time histories.
 By systematically varying soil yield strength, elastic stiffness, ductility ratio,
 over-strength factor, and viscous damping ratios of both soil and structure across a
 discrete design space, the code assembles a high-fidelity hysteretic model using multi-linear
 OpenSees hysteretic materials and parallel viscous dampers, extracts instantaneous eigenvalue
 periods at every analysis step to quantify period elongation, and computes a suite of performance
 metrics—including peak displacement, velocity and acceleration demands, energy-dissipation capacity 
 index, equivalent viscous damping ratio from hysteresis loops, ductility-based damage index, and 
 reconstruction cost. Extensive post-processing then generates correlation scatter plots,
 three-dimensional response surfaces, Pearson sensitivity heatmaps, Random-Forest safety
 classification, and both main-effect and interaction ANOVA tables, thereby ranking the relative
 influence of each soil and structural parameter on seismic performance and providing a rigorous
 quantitative basis for performance-based design and risk assessment of soil–structure systems.
"""
# BOOK: 
'https://www.fema.gov/sites/default/files/documents/fema-p-2091-soil-structure-interaction.pdf'    
# BOOK: Soil-Foundation-Structure Interaction
'https://www.routledge.com/Soil-Foundation-Structure-Interaction/Orense-Chouw-Pender/p/book/9781138114128'
# REPORT: Soil-Structure Interaction for Building Structures (NEHRP Consultants)
'https://www.nehrp.gov/pdf/nistgcr12-917-21.pdf'    
# BOOK: Dynamics of Structures in SI Units -  Anil Kumar Chopra
'https://share.google/TDN5O4eWmw5pH8zUH'
# BOOK: Differential Equations for Engineers-Wei-Chau Xie-CAMBRIDGE-2010
# BOOK: Structural dynamics THEORY AND COMPUTATION MARIO BAZ-5 EDITION
# BOOK: Nonlinear analysis in soil mechanics - WF Chen, WF Chen, E Mizuno - Elsevier Science Pub Co - 1990
# BOOK: Limit analysis in soil mechanics - WF Chen, XL Liu - Elsevier - 2012
# BOOK: Limit analysis and soil plasticity - WF Chen - Elsevier - 2013
# BOOK: Dynamics of Structure and Foundation; A Unified Approach-1. Fundamentals-Indrajit Chowdhury-CRC Press-2009
# BOOK: Dynamics of structure and foundation; A Uniﬁed Approach-2. Applications-Indrajit Chowdhury-CRC Press-2009
# BOOK: Foundation Analysis and Design-Joseph-E.-Bowles-1996
# REPORT: A Practical Guide o Soil-Structure Interaction_FEMA P-2091 
# PAPER: Displacement-based seismic design of buildings theory - M.S. Medhekar, D.J.L. Kennedy - 1999 Elsevier
# PAPER: Modelling of soil-structure interaction in OpenSees: A practical approach for performance-based seismic design
'https://www.sciencedirect.com/science/article/abs/pii/S2352012421000060?utm_source=chatgpt.com'
# YOUTUBE:
'https://www.youtube.com/watch?v=MZUhSHmIUdI'
#%%----------------------------------------------------
import numpy as np
import pandas as pd
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
import DISSIPATED_ENERGY_WITH_PLOT_FUN as S10
import EQULIVALENT_VISCOUS_DAMPING_RATIO_FUN as S055
import DAMAGE_INDEX_FUN as S066
import PLOT_CONTOUR_3D_2D_FUN as S11
import RECONSTRUCTION_COST_DI_FUN as S088
#%%----------------------------------------------------

def SOIL_STRUCTURE_MDOF(XIz_soil, XIz_str, FYz, KEz, DUCTz, OSFz, MAT_TYPE, TOTAL_MASS, ANAL_TYPEm, i):
    # Initialize model
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)
    
    MAX_ITERATIONS = 5000   # Maximum number of iterations
    MAX_TOLERANCE = 1.0e-6  # Specified tolerance for convergence
    GMfact = 9.81 # [m/s^2] standard acceleration of gravity or standard acceleration 
    
    # [kg] Mass
    MASS = [0.90 * TOTAL_MASS, # SOIL
            0.10 * TOTAL_MASS] # STRUCTURE  
    # Damping ratio
    DRi = [XIz_soil,  # SOIL
           XIz_str]   # STRUCTURE     
    
    # Define nodes
    ops.node(1, 0.0)  # Fixed base
    ops.node(2, 0.0)  # Mass node
    ops.node(3, 0.0)  # Mass node
    
        
    # Define boundary conditions
    ops.fix(1, 1)

    # Define masses
    for JJ in range(2, 4):
        ops.mass(JJ, MASS[JJ-2])
    
    #%% TWO DEGREES OF FREEDOM STRUCTURE, CALCULATE LATERAL STIFFNESS AND DAMPING 
    for II in range(0, 2):
        if II == 0: # SOIL
            FY = FYz                                         # [N] Yield Force of Soil 
            FU = OSFz * FY                                   # [N] Ultimate Force of Soil 
            Ke = KEz                                         # [N/m] Spring Elastic Stiffness
            DY = FY / Ke                                     # [m] Yield Displacement
            DSU = DUCTz * DY                                 # [m] Ultimate Displacement
            Ksh = (FU - FY) / (DSU - DY)                     # [N/m] Displacement Hardening Modulus
            Kp = FU / DSU                                    # [N/m] Spring Plastic Stiffness
            b = Ksh / Ke                                     # Displacement Hardening Ratio
            
        if II == 1: # STRUCTURE
            FY = 850_000.0                                   # [N] Yield Force of Structure
            FU = 1.18 * FY                                   # [N] Ultimate Force of Structure
            Ke = 4_500_000.0                                 # [N/m] Spring Elastic Stiffness
            DY = FY / Ke                                     # [m] Yield Displacement
            DSU = 0.36                                       # [m] Ultimate Displacement
            Ksh = (FU - FY) / (DSU - DY)                     # [N/m] Displacement Hardening Modulus
            Kp = FU / DSU                                    # [N/m] Spring Plastic Stiffness
            b = Ksh / Ke                                     # Displacement Hardening Ratio
            
        # Define material properties
        MAT_TAG = 1000 + II # SPRING TAG
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
            S08.OPENSEEES_HYSTERETICSM_FORCE_DISP_FUN(MAT_TAG, DP, FP, DN, FN, PLOT = False, X_LABEL='Displacement (mm)', Y_LABEL='Force [N]', TITLE='FORCE-DISPLACEMENT CURVE')
        if MAT_TYPE == 'ELASTIC':
            #ops.uniaxialMaterial('Elastic', MAT_TAG, Ke)             # TESNSION AND COMPRESSION IS SAME VALUES
            ops.uniaxialMaterial('Elastic', MAT_TAG, Ke ,0.0, 0.5*Ke) # TESNSION AND COMPRESSION IS NOT SAME VALUES
            # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/ElasticUni.html
        
        MatTag_C = 2000 + II # SPRING DAMPER
        alpha = 1.0    # velocity exponent (usually 0.3–1.0)
        omega = np.sqrt(Ke / MASS[II])
        Cd = 2 * DRi[II] * omega * MASS[II]  # [N·s/m] Damping coefficient 
        ops.uniaxialMaterial('Viscous', MatTag_C, Cd, alpha)  # Material for C (alpha=1.0 for linear)
        
        # Define element
        ops.element('zeroLength', II+1, II+1, II+2, '-mat',  MAT_TAG, MatTag_C, '-dir', 1, 1)  # DOF LATERAL SPRING

    center_node = 3
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(2, 1.0)
    ops.load(center_node, 1.0)
    
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
    DI, REC = [], []
    ele_force = {
        1: [],  # AXIALFORCE-01
        2: [],  # AXIALFORCE-02
        }
    # Initialize lists for each node's displacement
    node_displacements = {
        1: [],  # DISP01
        2: [],  # DISP02
        3: [],  # DISP03
        }
    if ANAL_TYPE == 'PERIOD': 
      #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(2, 0.5*DRi[0], DRi[0], 0, 1)
      PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(2, PLOT=True) 
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
        disp.append(ops.nodeDisp(center_node, 1))                       # DISPLACEMENT NODE 05 IN X DIR                
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
        DMAX = -DSU     # [m] Max. Displacement
        ops.integrator('DisplacementControl', center_node, IDctrlDOF, DINCR)
        ops.analysis('Static')
        Nsteps =  int(np.abs(DMAX/ DINCR)) 
        STEP = 0.0
        for step in range(Nsteps):
            OK = ops.analyze(1)
            S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS)
            ops.reactions()
            reaction.append(ops.nodeReaction(1, 1))                         # BASE REACTION
            disp.append(ops.nodeDisp(center_node, 1))                      # DISPLACEMENT NODE 05 IN X DIR       
            # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
            #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(2, 0.5*DR, DR, 0, 1)
            PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(2, PLOT=True)
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
        exec(open("COMPUTE_EFFECTIVE_PROPERTIES_FUN.py").read())
        
        return  DATA
    
    if ANAL_TYPE == 'CYCLIC_DISPLACEMENT': # STATIC TIME-HISTORY ANALYSIS
        IDctrlDOF = 1   # 1: Horizental Dispalcement - 2: Vertical Dispalcement
        DMAX = -DSU     # [m] Max. Displacement
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
            reaction.append(ops.nodeReaction(1, 1))                         # BASE REACTION
            disp.append(ops.nodeDisp(center_node, 1))                       # DISPLACEMENT NODE 05 IN X DIR       
            # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
            #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(2, 0.5*DR, DR, 0, 1)
            PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(2, PLOT=True)
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
            disp.append(ops.nodeDisp(center_node, 1))                       # DISPLACEMENT NODE 05 IN X DIR      
            # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
            #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(2, 0.5*DR, DR, 0, 1)
            PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(2, PLOT=True)
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
            disp.append(ops.nodeDisp(center_node, 1))             # DISPLACEMENT NODE 05 IN X DIR 
            velo.append(ops.nodeVel(center_node, 1))              # VELOCITY NODE 05
            acc.append(ops.nodeAccel(center_node, 1))             # ACCELERATION NODE 05
            stiffness.append(np.abs(reaction[-1] / disp[-1]))
            OMEGA.append(np.sqrt(stiffness[-1]/TOTAL_MASS))
            PERIOD.append((np.pi * 2) / OMEGA[-1])   
            # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
            #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(2, 0.5*DR, DR, 0, 1)
            PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(2, PLOT=True)
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
        u0 = -0.010                        # [m] Initial displacement
        v0 = 0.000015                      # [m/s] Initial velocity
        a0 = 0.000065                      # [m/s^2] Initial acceleration
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
            reaction.append(ops.nodeReaction(1, 1))               # BASE REACTION
            disp.append(ops.nodeDisp(center_node, 1))             # DISPLACEMENT NODE 05 IN X DIR 
            velo.append(ops.nodeVel(center_node, 1))              # VELOCITY NODE 05
            acc.append(ops.nodeAccel(center_node, 1))             # ACCELERATION NODE 05
            stiffness.append(np.abs(reaction[-1] / disp[-1]))
            OMEGA.append(np.sqrt(stiffness[-1]/TOTAL_MASS))
            PERIOD.append((np.pi * 2) / OMEGA[-1]) 
            # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
            #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(2, 0.5*DR, DR, 0, 1)
            PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(2, PLOT=True)
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
            print('\n\nFREE-VIBRATION ANALYSIS DONE.\n\n')    
        # Calculating Damping Ratio and Period Using Logarithmic Decrement Analysis 
        damping_ratio = S04.DAMPING_RATIO(disp)  
        
        # Compute modal properties
        ops.modalProperties("-print", "-file", f"SALAR_ModalReport_{ANAL_TYPE}.txt", "-unorm") 
        
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
        SSF_X = 15.0                       # Seismic Acceleration Scale Factor in X Direction
        SSF_Y = 15.0                       # Seismic Acceleration Scale Factor in Y Direction
        iv0_X = 0.0005                     # [m/s] Initial velocity applied to the node  in X Direction
        iv0_Y = 0.0005                     # [m/s] Initial velocity applied to the node  in Y Direction
        st_iv0 = 0.0                       # [s] Initial velocity applied starting time
        SEI = 'X'                          # Seismic Direction
        DR = 0.03                          # Damping ratio
        
        # Define time series for input motion (Acceleration time history)
        if SEI == 'X':
            SEISMIC_TAG_01 = 100
            gm_accels = np.loadtxt(f'Ground_Acceleration_{i+1}.txt')  # Assumes acceleration in m/s²
            ops.timeSeries('Path', SEISMIC_TAG_01, '-dt', dt, '-values', *gm_accels.tolist(), '-factor', GMfact) # SEISMIC-X
            # Define load patterns
            # pattern UniformExcitation $patternTag $dof -accel $tsTag <-vel0 $vel0> <-fact $cFact>
            ops.pattern('UniformExcitation', SEISMIC_TAG_01, 1, '-accel', SEISMIC_TAG_01, '-vel0', iv0_X, '-fact', SSF_X) # SEISMIC-X
        if SEI == 'Y':
            SEISMIC_TAG_02 = 200
            gm_accels = np.loadtxt(f'Ground_Acceleration_{i+1}.txt')  # Assumes acceleration in m/s²
            ops.timeSeries('Path', SEISMIC_TAG_02, '-dt', dt, '-values', *gm_accels.tolist(), '-factor', GMfact, '-vel0', iv0_Y) # SEISMIC-Y
            # Define load patterns
            # pattern UniformExcitation $patternTag $dof -accel $tsTag <-vel0 $vel0> <-fact $cFact>
            ops.pattern('UniformExcitation', SEISMIC_TAG_02, 2, '-accel', SEISMIC_TAG_02, '-vel0', iv0_Y, '-fact', SSF_Y) # SEISMIC-Y
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
            disp.append(ops.nodeDisp(center_node, 1))             # DISPLACEMENT NODE 05 IN X DIR 
            velo.append(ops.nodeVel(center_node, 1))              # VELOCITY NODE 05
            acc.append(ops.nodeAccel(center_node, 1))             # ACCELERATION NODE 05
            stiffness.append(np.abs(reaction[-1] / disp[-1]))
            OMEGA.append(np.sqrt(stiffness[-1]/TOTAL_MASS))
            PERIOD.append((np.pi * 2) / OMEGA[-1])
            # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
            #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(2, 0.5*DR, DR, 0, 1)
            PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(2, PLOT=True)
            PERIOD_MIN.append(PERIODmin)
            PERIOD_MAX.append(PERIODmax)
            # EVALUATION OF DUCTILITY DAMAGE INDEX
            if MAT_TYPE == 'INELASTIC':
                di = S066.DAMAGE_INDEX_FUN(disp[-1], DY, DSU)
                DI.append(di)                   # DAMAGE INDEX
                reci = S088.RECONSTRUCTION_COST_DI_FUN(disp[-1], DY, DSU, 1_000_000) # EVALUATION OF RECONSTRUCTION COST: 1_000_1000
                REC.append(reci)
            if MAT_TYPE == 'ELASTIC':
                DI.append(0.0)
                REC.append(0.0)
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
        
        DATA = (time, reaction, disp, velo, acc, DI, REC,
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
TOTAL_MASS = 1_000_000.0   # [kg] Total Mass of Structure
#%%----------------------------------------------------
MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'

# --------------------------------------------------------------
# SENSITIVITY ANALYSIS BY CHANGING SOIL PARAMETERS
# --------------------------------------------------------------

import time as TI

FY_MIN = 0.80 * 13_320_000       # [N] MIN. SOIL'S YIELD STRENGTH 
FY_MAX = 1.20 * 13_320_000       # [N] MAX. SOIL'S YIELD STRENGTH 

KE_MIN = 0.80 * 74_000_000.0     # [N/m] MIN. SOIL'S ELASTIC STIFFNESS  
KE_MAX = 1.20 * 74_000_000.0     # [N/m] MAX. SOIL'S ELASTIC STIFFNESS  

DUCT_MIN = 20.0                  # [m/m] MIN. SOIL'S DUCTILITY RATIO
DUCT_MAX = 50.0                  # [m/m] MAX. SOIL'S DUCTILITY RATIO

OSF_MIN = 1.05                   # [N/N] MIN. SOIL'S OVER-STRENGTH FACTOR
OSF_MAX = 1.20                   # [N/N] MAX. SOIL'S OVER-STRENGTH FACTOR
 
XI_soil_MIN =  0.001             # MIN. SOIL'S DAMPING RATIO
XI_soil_MAX =  0.10              # MAX. SOIL'S DAMPING RATIO

XI_str_MIN = 0.01                # MIN. STRUCTURE'S DAMPING RATIO
XI_str_MAX = 0.05                # MAX. STRUCTURE'S DAMPING RATIO

NUM_SIM = 3

III = int(25)        # Select Seismic  Ground Motion Number 26 (II+1)

SOIL_DUCT = []   
SOIL_FY = []
SOIL_KE = []
SOIL_OSF = []
SOIL_XI = []
STR_XI = []
min_T = []   # max period from the seismic analysis
max_T = []   # max period from the seismic analysis
max_edci = []
max_zeta = []
max_di, max_rec = [], []
max_disp, max_velo, max_acc =  [], [], []

STEP = 0
# Analysis Durations:
starttime = TI.process_time() # start the CPU clock



for AA in range(0, NUM_SIM):                                   # SOIL DAMPING RATIO
    # linear interpolation of the soil damping ratio between XI_soil_MIN and XI_soil_MAX
    XIz_soil = XI_soil_MIN + (XI_soil_MAX - XI_soil_MIN) * (AA / (NUM_SIM-1))
    
    for BB in range(0, NUM_SIM):                               # STRUCTURE DAMPING RATIO
        # linear interpolation of the structure damping ratio between XI_str_MIN and XI_str_MAX
        XIz_str = XI_str_MIN + (XI_str_MAX - XI_str_MIN) * (BB / (NUM_SIM-1))    
    
        for XX in range(0, NUM_SIM):                               # SOIL YIELD STRENGTH [N]
            # linear interpolation of the soil yield strength between FY_MIN and FY_MAX
            FYz = FY_MIN + (FY_MAX - FY_MIN) * (XX / 5)
            
            for YY in range(0, NUM_SIM):                           # SOIL ELASTIC STIFFNESS [N/m]
                # linear interpolation of the soil elastic stiffness between KE_MIN and KE_MAX
                KEz = KE_MIN + (KE_MAX - KE_MIN) * (YY / (NUM_SIM-1))
                
                for ZZ in range(0, NUM_SIM):                       # SOIL OVER-STRENGTH FACTOR [N/N]
                    # linear interpolation of the soil over-strength factor between OSF_MIN and OSF_MAX
                    OSFz = OSF_MIN + (OSF_MAX - OSF_MIN) * (ZZ / (NUM_SIM-1))   
                    
                    for JJ in range(0, NUM_SIM):                   # SOIL DUCTILITY RATIO [m/m]
                        # linear interpolation of the ductility ratio between DUCT_MIN and DUCT_MAX
                        DUCTz = DUCT_MIN + (DUCT_MAX - DUCT_MIN) * (JJ / (NUM_SIM-1))
                        
                        SOIL_XI.append(XIz_soil)             # store the current variable value
                        STR_XI.append(XIz_str)               # store the current variable value
                        SOIL_DUCT.append(DUCTz)              # store the current variable value
                        SOIL_FY.append(FYz)                  # store the current variable value
                        SOIL_KE.append(KEz)                  # store the current variable value
                        SOIL_OSF.append(OSFz)                # store the current variable value
                    
                        
                        # CYCLIC DISPLACEMENT ANALYSIS (STATIC TIME-HISTORY ANALYSIS)
                        MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'
                        ANAL_TYPE = 'CYCLIC_DISPLACEMENT'

                        DATA = SOIL_STRUCTURE_MDOF(XIz_soil, XIz_str, FYz, KEz, DUCTz, OSFz, MAT_TYPE, TOTAL_MASS, ANAL_TYPE, 0)         # call nonlinear SOIL_STRUCTURE_MDOF seissmic solver
                        (reaction_CP, disp_CP, DI_CP,
                         ele_force_CP, node_displacements_CP,
                         PERIOD_MIN_CP, PERIOD_MAX_CP) = DATA

                        # SEISMIC ANALYSIS
                        ANAL_TYPE = 'SEISMIC'
                        DATA = SOIL_STRUCTURE_MDOF(XIz_soil, XIz_str, FYz, KEz, DUCTz, OSFz, MAT_TYPE, TOTAL_MASS, ANAL_TYPE, III)   # call nonlinear SOIL_STRUCTURE_MDOF seissmic solver
                        
                        (time_SEI, reaction_SEI, disp_SEI, velo_SEI, acc_SEI, DI_SEI, REC_SEI,
                                ele_force_SEI, node_displacements_SEI,
                                stiffness_SEI, PERIOD_SEI, damping_ratio_SEI,
                                PERIOD_MIN_SEI, PERIOD_MAX_SEI) = DATA    # unpack solver output
                                             
                        
                        # ENERGY DISSIPATION CAPACITY INDEX
                        edci = S10.ENERGY_DISSIPATION_CAPACITY_INDEX(
                                   disp_SEI, reaction_SEI, disp_CP, reaction_CP)
                        
                            
                        #%% START - Equivalent viscous damping ratio
                        XLABEL = "Displacement [m]"
                        YLABEL = "Base Reaction [N]"
                        TITLE = "Equivalent viscous damping ratio - Dynamic Analysis Hysteresis"
                        method = 2 # 
                        zeta = S055.EQULIVALENT_VISCOUS_DAMPING_RATIO_FUN(disp_SEI, reaction_SEI, method, XLABEL, YLABEL, TITLE)
                        print(f"Equivalent viscous damping ratio = {zeta:.4f}")
                        #%% END - Equivalent viscous damping ratio
                        
                        # PLOT THE SPRING
                        #S01.PLOT_1D_SPRING(deformed_scale=1.0, virtual_spring_length=1.0, virtual_spring_angle=0.0)
                        
                    
                        # store the maximum (final) period of the seismic curve
                        min_T.append(np.max(PERIOD_MIN_SEI)) 
                        max_T.append(np.max(PERIOD_MAX_SEI)) 
                        
                        max_disp.append(np.max(np.abs(disp_SEI)))
                        max_velo.append(np.max(np.abs(velo_SEI)))
                        max_acc.append(np.max(np.abs(acc_SEI)))
                        max_di.append(np.max(np.abs(DI_SEI)))
                        max_rec.append(np.max(np.abs(REC_SEI)))
                        max_edci.append(edci)
                        max_zeta.append(zeta)
                        

                        STEP = STEP + 1
                        #if STEP == 100:
                        #    break
                        print(f'\t\t\t STEP {STEP} DONE. \n\n')


#exit()
# TIMER REPORT
totaltime = TI.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n')
#%%----------------------------------------------------
SOIL_XI = np.array(SOIL_XI) * 100
STR_XI = np.array(STR_XI) * 100

XDATA = max_disp
YDATA = min_T
XLABEL = 'Max. Structure Displacent [m]'
YLABEL = 'MIN. STRUCTURAL PERIOD [s]' 
TITLE = f'{YLABEL} and {XLABEL} DURING SEISMIC ANALYSIS'
COLOR = 'brown'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = max_disp
YDATA = max_T
XLABEL = 'Max. Structure Displacent [m]'
YLABEL = 'MAX. STRUCTURAL PERIOD [s]' 
TITLE = f'{YLABEL} and {XLABEL} DURING SEISMIC ANALYSIS'
COLOR = 'brown'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = max_velo
YDATA = max_T
XLABEL = 'Max. Structure Velocity [m/s]'
YLABEL = 'MAX. STRUCTURAL PERIOD [s]' 
TITLE = f'{YLABEL} and {XLABEL} DURING SEISMIC ANALYSIS'
COLOR = 'brown'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = max_acc
YDATA = max_T
XLABEL = 'Max. Structure Acceleration [m/s^2]'
YLABEL = 'MAX. STRUCTURAL PERIOD [s]' 
TITLE = f'{YLABEL} and {XLABEL} DURING SEISMIC ANALYSIS'
COLOR = 'brown'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = SOIL_XI
YDATA = max_T
XLABEL = 'Soil Damping Ratio [N]'
YLABEL = 'MAX. STRUCTURAL PERIOD [s]' 
TITLE = f'{YLABEL} and {XLABEL} DURING SEISMIC ANALYSIS'
COLOR = 'brown'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)
                
XDATA = SOIL_FY
YDATA = max_T
XLABEL = 'Soil Yield Strength [N]'
YLABEL = 'MAX. STRUCTURAL PERIOD [s]' 
TITLE = f'{YLABEL} and {XLABEL} DURING SEISMIC ANALYSIS'
COLOR = 'brown'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = SOIL_KE
YDATA = max_T
XLABEL = 'Soil Elastic Stiffness [N/m]'
YLABEL = 'MAX. STRUCTURAL PERIOD [s]' 
TITLE = f'{YLABEL} and {XLABEL} DURING SEISMIC ANALYSIS'
COLOR = 'brown'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = SOIL_DUCT
YDATA = max_T
XLABEL = 'Soil ductility Ratio [m/m]'
YLABEL = 'MAX. STRUCTURAL PERIOD [s]' 
TITLE = f'{YLABEL} and {XLABEL} DURING SEISMIC ANALYSIS'
COLOR = 'brown'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = SOIL_OSF
YDATA = max_T
XLABEL = 'Soil Over-strength Factor [N/N]'
YLABEL = 'MAX. STRUCTURAL PERIOD [s]' 
TITLE = f'{YLABEL} and {XLABEL} DURING SEISMIC ANALYSIS'
COLOR = 'brown'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = STR_XI
YDATA = max_T
XLABEL = 'Structure Damping Ratio [N]'
YLABEL = 'MAX. STRUCTURAL PERIOD [s]' 
TITLE = f'{YLABEL} and {XLABEL} DURING SEISMIC ANALYSIS'
COLOR = 'brown'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = max_disp
YDATA = max_edci
XLABEL = 'Max. Structure Displacent [m]'
YLABEL = 'Soil-Structure Energy Dissipation Capacity Index [%]' 
TITLE = f'{YLABEL} and {XLABEL} DURING SEISMIC ANALYSIS'
COLOR = 'blue'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = max_velo
YDATA = max_edci
XLABEL = 'Max. Structure Velocity [m/s]'
YLABEL = 'Soil-Structure Energy Dissipation Capacity Index [%]' 
TITLE = f'{YLABEL} and {XLABEL} DURING SEISMIC ANALYSIS'
COLOR = 'blue'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = max_acc
YDATA = max_edci
XLABEL = 'Max. Structure Acceleration [m/s^2]'
YLABEL = 'Soil-Structure Energy Dissipation Capacity Index [%]' 
TITLE = f'{YLABEL} and {XLABEL} DURING SEISMIC ANALYSIS'
COLOR = 'blue'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = SOIL_XI
YDATA = max_edci
XLABEL = 'Soil Damping Ratio [N]'
YLABEL = 'Soil-Structure Energy Dissipation Capacity Index [%]' 
TITLE = f'{YLABEL} and {XLABEL} DURING SEISMIC ANALYSIS'
COLOR = 'blue'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = SOIL_FY
YDATA = max_edci
XLABEL = 'Soil Yield Strength [N]'
YLABEL = 'Soil-Structure Energy Dissipation Capacity Index [%]' 
TITLE = f'{YLABEL} and {XLABEL} DURING SEISMIC ANALYSIS'
COLOR = 'blue'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = SOIL_KE
YDATA = max_edci
XLABEL = 'Soil Elastic Stiffness [N/m]'
YLABEL = 'Soil-Structure Energy Dissipation Capacity Index [%]' 
TITLE = f'{YLABEL} and {XLABEL} DURING SEISMIC ANALYSIS'
COLOR = 'blue'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = SOIL_OSF
YDATA = max_edci
XLABEL = 'Soil Over-strength Factor [N/N]'
YLABEL = 'Soil-Structure Energy Dissipation Capacity Index [%]' 
TITLE = f'{YLABEL} and {XLABEL} DURING SEISMIC ANALYSIS'
COLOR = 'blue'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = SOIL_DUCT
YDATA = max_edci
XLABEL = 'Soil ductility Ratio [m/m]'
YLABEL = 'Soil-Structure Energy Dissipation Capacity Index [%]' 
TITLE = f'{YLABEL} and {XLABEL} DURING SEISMIC ANALYSIS'
COLOR = 'blue'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = STR_XI
YDATA = max_edci
XLABEL = 'Structure Damping Ratio [%]'
YLABEL = 'Soil-Structure Energy Dissipation Capacity Index [%]'
TITLE = f'{YLABEL} and {XLABEL} DURING SEISMIC ANALYSIS'
COLOR = 'blue'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = max_disp
YDATA = max_zeta
XLABEL = 'Max. Structure Displacent [m]'
YLABEL = 'Soil-Structure Equivalent viscous damping ratio [%]'  
TITLE = f'{YLABEL} and {XLABEL} DURING SEISMIC ANALYSIS'
COLOR = 'blue'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = max_velo
YDATA = max_zeta
XLABEL = 'Max. Structure Velocity [m/s]'
YLABEL = 'Soil-Structure Equivalent viscous damping ratio [%]'  
TITLE = f'{YLABEL} and {XLABEL} DURING SEISMIC ANALYSIS'
COLOR = 'blue'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = max_acc
YDATA = max_zeta
XLABEL = 'Max. Structure Acceleration [m/s^2]'
YLABEL = 'Soil-Structure Equivalent viscous damping ratio [%]'  
TITLE = f'{YLABEL} and {XLABEL} DURING SEISMIC ANALYSIS'
COLOR = 'blue'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = SOIL_XI
YDATA = max_zeta
XLABEL = 'Soil Damping Ratio [N]'
YLABEL = 'Soil-Structure Equivalent viscous damping ratio [%]'  
TITLE = f'{YLABEL} and {XLABEL} DURING SEISMIC ANALYSIS'
COLOR = 'red'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = SOIL_FY
YDATA = max_zeta
XLABEL = 'Soil Yield Strength [N]'
YLABEL = 'Soil-Structure Equivalent viscous damping ratio [%]'  
TITLE = f'{YLABEL} and {XLABEL} DURING SEISMIC ANALYSIS'
COLOR = 'red'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = SOIL_KE
YDATA = max_zeta
XLABEL = 'Soil Elastic Stiffness [N/m]'
YLABEL = 'Soil-Structure Equivalent viscous damping ratio [%]'  
TITLE = f'{YLABEL} and {XLABEL} DURING SEISMIC ANALYSIS'
COLOR = 'red'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = SOIL_OSF
YDATA = max_zeta
XLABEL = 'Soil Over-strength Factor [N/N]'
YLABEL = 'Soil-Structure Equivalent viscous damping ratio [%]' 
TITLE = f'{YLABEL} and {XLABEL} DURING SEISMIC ANALYSIS'
COLOR = 'red'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = SOIL_DUCT
YDATA = max_zeta
XLABEL = 'Soil ductility Ratio [m/m]'
YLABEL = 'Soil-Structure Equivalent viscous damping ratio [%]' 
TITLE = f'{YLABEL} and {XLABEL} DURING SEISMIC ANALYSIS'
COLOR = 'red'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = STR_XI
YDATA = max_zeta
XLABEL = 'Soil ductility Ratio [m/m]'
YLABEL = 'Soil-Structure Equivalent viscous damping ratio [%]' 
TITLE = f'{YLABEL} and {XLABEL} DURING SEISMIC ANALYSIS'
COLOR = 'red'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)


XDATA = max_disp
YDATA = max_di
XLABEL = 'Max. Structure Displacent [m]'
YLABEL = 'Structure Ductility Damage Index [%]'   
TITLE = f'{YLABEL} and {XLABEL} DURING SEISMIC ANALYSIS'
COLOR = 'green'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = max_velo
YDATA = max_di
XLABEL = 'Max. Structure Velocity [m/s]'
YLABEL = 'Structure Ductility Damage Index [%]'   
TITLE = f'{YLABEL} and {XLABEL} DURING SEISMIC ANALYSIS'
COLOR = 'green'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = max_acc
YDATA = max_di
XLABEL = 'Max. Structure Acceleration [m/s^2]'
YLABEL = 'Structure Ductility Damage Index [%]'   
TITLE = f'{YLABEL} and {XLABEL} DURING SEISMIC ANALYSIS'
COLOR = 'green'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = SOIL_XI
YDATA = max_di
XLABEL = 'Soil Damping Ratio [N]'
YLABEL = 'Structure Ductility Damage Index [%]'    
TITLE = f'{YLABEL} and {XLABEL} DURING SEISMIC ANALYSIS'
COLOR = 'green'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = SOIL_FY
YDATA = max_di
XLABEL = 'Soil Yield Strength [N]'
YLABEL = 'Structure Ductility Damage Index [%]'  
TITLE = f'{YLABEL} and {XLABEL} DURING SEISMIC ANALYSIS'
COLOR = 'green'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = SOIL_KE
YDATA = max_di
XLABEL = 'Soil Elastic Stiffness [N/m]'
YLABEL = 'Structure Ductility Damage Index [%]' 
TITLE = f'{YLABEL} and {XLABEL} DURING SEISMIC ANALYSIS'
COLOR = 'green'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = SOIL_OSF
YDATA = max_di
XLABEL = 'Soil Over-strength Factor [N/N]'
YLABEL = 'Structure Ductility Damage Index [%]' 
TITLE = f'{YLABEL} and {XLABEL} DURING SEISMIC ANALYSIS'
COLOR = 'green'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = SOIL_DUCT
YDATA = max_di
XLABEL = 'Soil ductility Ratio [m/m]'
YLABEL = 'Structure Ductility Damage Index [%]' 
TITLE = f'{YLABEL} and {XLABEL} DURING SEISMIC ANALYSIS'
COLOR = 'green'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = STR_XI
YDATA = max_di
XLABEL = 'Soil ductility Ratio [m/m]'
YLABEL = 'Structure Ductility Damage Index [%]' 
TITLE = f'{YLABEL} and {XLABEL} DURING SEISMIC ANALYSIS'
COLOR = 'red'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)



XDATA = max_disp
YDATA = max_rec
XLABEL = 'Max. Structure Displacent [m]'
YLABEL = 'Structure Reconstruction Cost [$]'   
TITLE = f'{YLABEL} and {XLABEL} DURING SEISMIC ANALYSIS'
COLOR = 'green'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = max_velo
YDATA = max_rec
XLABEL = 'Max. Structure Velocity [m/s]'
YLABEL = 'Structure Reconstruction Cost [$]'   
TITLE = f'{YLABEL} and {XLABEL} DURING SEISMIC ANALYSIS'
COLOR = 'green'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = max_acc
YDATA = max_rec
XLABEL = 'Max. Structure Acceleration [m/s^2]'
YLABEL = 'Structure Reconstruction Cost [$]'   
TITLE = f'{YLABEL} and {XLABEL} DURING SEISMIC ANALYSIS'
COLOR = 'green'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = SOIL_XI
YDATA = max_rec
XLABEL = 'Soil Damping Ratio [N]'
YLABEL = 'Structure Reconstruction Cost [$]'    
TITLE = f'{YLABEL} and {XLABEL} DURING SEISMIC ANALYSIS'
COLOR = 'green'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = SOIL_FY
YDATA = max_rec
YLABEL = 'Structure Reconstruction Cost [$]'    
TITLE = f'{YLABEL} and {XLABEL} DURING SEISMIC ANALYSIS'
COLOR = 'green'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = SOIL_KE
YDATA = max_rec
XLABEL = 'Soil Elastic Stiffness [N/m]'
YLABEL = 'Structure Reconstruction Cost [$]'   
TITLE = f'{YLABEL} and {XLABEL} DURING SEISMIC ANALYSIS'
COLOR = 'green'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = SOIL_OSF
YDATA = max_rec
XLABEL = 'Soil Over-strength Factor [N/N]'
YLABEL = 'Structure Reconstruction Cost [$]'    
TITLE = f'{YLABEL} and {XLABEL} DURING SEISMIC ANALYSIS'
COLOR = 'green'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = SOIL_DUCT
YDATA = max_rec
XLABEL = 'Soil ductility Ratio [m/m]'
YLABEL = 'Structure Reconstruction Cost [$]'   
TITLE = f'{YLABEL} and {XLABEL} DURING SEISMIC ANALYSIS'
COLOR = 'green'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = STR_XI
YDATA = max_rec
XLABEL = 'Soil ductility Ratio [m/m]'
YLABEL = 'Structure Reconstruction Cost [$]'    
TITLE = f'{YLABEL} and {XLABEL} DURING SEISMIC ANALYSIS'
COLOR = 'red'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)



# 3D PLOT
X, Y, Z = max_disp, SOIL_KE, min_T
XLABEL, YLABEL, ZLABEL = 'MAx. Structure Displacement [m]', 'Soil Stiffness [N/m]', 'Structure Min. Period [s]',               
S11.PLOT_CONTOUR_3D_2D_FUN(120, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = max_disp, SOIL_KE, max_T
XLABEL, YLABEL, ZLABEL = 'MAx. Structure Displacement [m]', 'Soil Stiffness [N/m]', 'Structure Max. Period [s]',               
S11.PLOT_CONTOUR_3D_2D_FUN(120, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = max_velo, SOIL_KE, max_T
XLABEL, YLABEL, ZLABEL = 'MAx. Structure velocity [m/s]', 'Soil Stiffness [N/m]', 'Structure Max. Period [s]',               
S11.PLOT_CONTOUR_3D_2D_FUN(120, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = max_acc, SOIL_KE, max_T
XLABEL, YLABEL, ZLABEL = 'MAx. Structure Acceleration [m/s^2]', 'Soil Stiffness [N/m]', 'Structure Max. Period [s]',               
S11.PLOT_CONTOUR_3D_2D_FUN(120, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = SOIL_FY, SOIL_KE, max_T
XLABEL, YLABEL, ZLABEL = 'Soil Yield Strength [N]', 'Soil Stiffness [N/m]', 'Structure Max. Period [s]',               
S11.PLOT_CONTOUR_3D_2D_FUN(120, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = SOIL_XI, STR_XI, max_T
XLABEL, YLABEL, ZLABEL = 'Soil Damping Ratio [%]', 'Structure Damping Ratio [%]', 'Structure Max. Period [s]',               
S11.PLOT_CONTOUR_3D_2D_FUN(120, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = SOIL_FY, SOIL_KE, max_zeta
XLABEL, YLABEL, ZLABEL = 'Soil Elastic Stiffness [N/m]', 'Soil Stiffness [N/m]', 'Soil-Structure Equivalent viscous damping ratio [%]',               
S11.PLOT_CONTOUR_3D_2D_FUN(120, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = SOIL_XI, STR_XI, max_zeta
XLABEL, YLABEL, ZLABEL = 'Soil Damping Ratio [%]', 'Structure Damping Ratio [%]', 'Soil-Structure Equivalent viscous damping ratio [%]',               
S11.PLOT_CONTOUR_3D_2D_FUN(120, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = SOIL_FY, SOIL_KE, max_edci
XLABEL, YLABEL, ZLABEL = 'Soil Elastic Stiffness [N/m]', 'Soil Stiffness [N/m]', 'Soil-Structure Energy Dissipation Capacity Index [%]', 
S11.PLOT_CONTOUR_3D_2D_FUN(120, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = SOIL_XI, STR_XI, max_edci
XLABEL, YLABEL, ZLABEL = 'Soil Damping Ratio [%]', 'Structure Damping Ratio [%]', 'Soil-Structure Energy Dissipation Capacity Index [%]',               
S11.PLOT_CONTOUR_3D_2D_FUN(120, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = SOIL_FY, SOIL_KE, max_di
XLABEL, YLABEL, ZLABEL = 'Soil Elastic Stiffness [N/m]', 'Soil Stiffness [N/m]', 'Structure Damage Index [%]', 
S11.PLOT_CONTOUR_3D_2D_FUN(120, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = SOIL_XI, STR_XI, max_di
XLABEL, YLABEL, ZLABEL = 'Soil Damping Ratio [%]', 'Structure Damping Ratio [%]', 'Structure Damage Index [%]',               
S11.PLOT_CONTOUR_3D_2D_FUN(120, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = SOIL_XI, STR_XI, max_disp
XLABEL, YLABEL, ZLABEL = 'Soil Damping Ratio [%]', 'Structure Damping Ratio [%]', 'Max. Structure Displacement [m]',               
S11.PLOT_CONTOUR_3D_2D_FUN(120, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = SOIL_XI, STR_XI, max_velo
XLABEL, YLABEL, ZLABEL = 'Soil Damping Ratio [%]', 'Structure Damping Ratio [%]', 'Max. Structure Velocity [m/s]',               
S11.PLOT_CONTOUR_3D_2D_FUN(120, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = SOIL_XI, STR_XI, max_acc
XLABEL, YLABEL, ZLABEL = 'Soil Damping Ratio [%]', 'Structure Damping Ratio [%]', 'Max. Structure Acceleration [m]',               
S11.PLOT_CONTOUR_3D_2D_FUN(120, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = SOIL_XI, STR_XI, min_T
XLABEL, YLABEL, ZLABEL = 'Soil Damping Ratio [%]', 'Structure Damping Ratio [%]', 'Min. Structure Period [s]',               
S11.PLOT_CONTOUR_3D_2D_FUN(120, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = SOIL_XI, STR_XI, max_rec
XLABEL, YLABEL, ZLABEL = 'Soil Damping Ratio [%]', 'Structure Damping Ratio [%]', 'Reconsruction Cost [$]',               
S11.PLOT_CONTOUR_3D_2D_FUN(120, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = SOIL_FY, SOIL_KE, max_rec
XLABEL, YLABEL, ZLABEL = 'Soil Elastic Stiffness [N/m]', 'Soil Stiffness [N/m]', 'Reconsruction Cost [$]',               
S11.PLOT_CONTOUR_3D_2D_FUN(120, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = max_zeta, max_edci, max_rec
XLABEL, YLABEL, ZLABEL = 'Soil-Structure Equivalent viscous damping ratio [%]', 'Soil-Structure Energy Dissipation Capacity Index [%]', 'Reconsruction Cost [$]',               
S11.PLOT_CONTOUR_3D_2D_FUN(120, X, Y, Z, XLABEL, YLABEL, ZLABEL)

#%%------------------------------------------------------
# RANDOM FOREST ANALYSIS
"""
This code predicts the seismic safety of a structure using simulation data by training a Random Forest Classifier to
 classify whether the system is "safe" or "unsafe" based on features like maximum displacement, velocity, acceleration,
 and base reaction. A regression model is also trained to estimate safety likelihood. It evaluates model performance using
 metrics like classification accuracy, mean squared error, and R² score. Additionally, it identifies key features influencing
 safety through feature importance analysis. The tool aids in seismic risk assessment, structural optimization, and understanding
 critical safety parameters.
"""

data = {
    "SOIL_DR":                            SOIL_XI,
    "STR_DR":                             STR_XI,
    "SOIL_FY":                            SOIL_FY,
    "SOIL_KE":                            SOIL_KE,        
    "SOIL_DUCT":                          SOIL_DUCT,
    "SOIL_OSF":                           SOIL_OSF,
    "DISPLACEMENT":                       max_disp,
    "VELOCITY":                           max_velo,
    "ACCELERION":                         max_acc,
    "MIN-PEIOD":                          min_T,
    "MAX-PEIOD":                          max_T,
    "EDCI":                               max_edci,
    "EDVR":                               max_zeta,
}


# Convert to DataFrame
df = pd.DataFrame(data)
#%%------------------------------------------------------
# PLOT HEATMAP FOR CORRELATION 
S09.PLOT_HEATMAP(df)
#print(df)


threshold_ductility_ratio = 35.0
S09.RANDOM_FOREST(df, threshold_ductility_ratio)
#%%------------------------------------------------------
# MULTIPLE REGRESSION MODEL
#S01.MULTIPLE_REGRESSION(df) 
#%%-------------------------------------------------------------------
# Plots a heatmap of sensitivity coefficients (correlation or SRC) between inputs X and outputs Y.
import SENSITIVITY_HEATMAP_FUN as S099
X = np.column_stack([SOIL_FY, SOIL_KE, SOIL_DUCT, SOIL_OSF, SOIL_XI])          
Y = np.column_stack([max_disp, max_acc, max_edci, max_zeta, max_T])  

X_LABELS = ['Soil Yield Strength [N]',
            'Soil Elastic Stiffness [N/m]',
            'Soil Ductility Ratio [m/m]',
            'Soil Over-strength Factor [N/N]',
            'Soil Damping Ratio [%]']

Y_LABELS = ['Max. Disp. [m]',
            'Max. Acc. [m]',
            'Soil-Structure Energy Dissipation Capacity Index [%]',
            'Soil-Structure Equivalent viscous damping ratio [%]',
            'Soil-Structure Max. Period [s]']

coeffs = S099.SENSITIVITY_HEATMAP_FUN(X, Y, X_LABELS, Y_LABELS, method='pearson')
#%%-------------------------------------------------------------------
# ANOVA with automatic binning of continuous predictors.
import ANOVA_SENSITIVITY_FUN as S100

df_sens = pd.DataFrame({
    "SOIL_DR":                            SOIL_XI,
    "STR_DR":                             STR_XI,
    "SOIL_FY":                            SOIL_FY,
    "SOIL_KE":                            SOIL_KE,        
    "SOIL_DUCT":                          SOIL_DUCT,
    "SOIL_OSF":                           SOIL_OSF,
    "DISPLACEMENT":                       max_disp,
    "VELOCITY":                           max_velo,
    "ACCELERION":                         max_acc,
    "EDVR":                               max_zeta,
})

param_list = ["SOIL_DR",
              "STR_DR", 
              "SOIL_FY",
              "SOIL_KE",
              "SOIL_DUCT",
              "SOIL_OSF",
              "DISPLACEMENT",
              "VELOCITY",
              "ACCELERION"]

# ANOVA – main effects only
anova_table, _ = S100.ANOVA_SENSITIVITY_FUN(
    df_sens, output_col="EDVR",
    param_cols=param_list,
    n_bins=4,
    include_interactions=False,
    plot=True,
)
plt.show()
print("\nANOVA table (main effects):")
print(anova_table.round(4))

# ANOVA – with two-way interactions
anova_table_inter, _ = S100.ANOVA_SENSITIVITY_FUN(
    df_sens, output_col="EDVR",
    param_cols=param_list,
    n_bins=4,
    include_interactions=True,
    plot=True,
)
plt.show()
print("\nANOVA table (with interactions):")
print(anova_table_inter.round(4))

#%%----------------------------------------------------
