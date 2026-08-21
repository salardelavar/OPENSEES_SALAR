#############################################################################################################
#                          >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                 #
#SENSITIVITY ANALYSIS OF POST-FIRE SINGLE-DEGREE-FREEOM (SDOF) STRUCTURES USING INCREMENTAL DYNAMIC ANALYSIS#
#         EFFECTS OF THERMAL, STRUCTURAL MASS, ELASTIC STIFFNESS, STRUCTURAL DUCTILITY RATIO                #
#       AND OVER-STRENGTH FACTOR ON OUTPUT KEY PARAMETERS SUCH AS EQUIVALENT VISCOUS DAMPING RATIO          #
#                        FROM NONLINEAR DYNAMIC ANALYSES USING PYTHON AND OPENSEES                          #
#-----------------------------------------------------------------------------------------------------------#
#                   EQUIVALENT VISCOUS DAMPING RATIO: xi_eq = 100 * E_d / (4 * pi * E_s)                    #
#-----------------------------------------------------------------------------------------------------------#
#                    THIS PYTHON SCRIPT IS WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                    #
#                                   EMAIL: salar.d.ghashghaei@gmail.com                                     #
#############################################################################################################
"""
1. This script performs a "sensitivity analysis of a post-fire single-degree-of-freedom (SDOF) structure" using "Incremental Dynamic" in OpenSeesPy.
2. The objective is to study the effects of "initial displacement, sructural elastic stiffness, ductility ratio, damping ratio, and over-strength factor" on nonlinear dynamic response.
3. Structural properties such as "yield force, ultimate force, elastic stiffness, hardening ratio, and damping ratio" are defined first.
4. Elastic and plastic periods of the system are computed and reported.
5. Key seismic performance parameters are calculated, including "over-strength (Ω₀)", "ductility (μ)", and "behavior factor (R)".
6. A nonlinear "hysteretic spring model" is used to represent structural behavior with strength degradation.
7. A "viscous damper" is added to simulate energy dissipation.
8. The analysis considers "Incremental Dynamic due to initial displacement" (no external ground motion).
9. The function `ANALYSIS_SDOF` builds the OpenSees model, applies initial conditions, and runs a transient analysis.
10. "Newmark time integration" and "Newton–Raphson iteration" are used for nonlinear solution.
11. "Rayleigh damping coefficients" are computed from the first eigenvalue.
12. Time histories of "displacement, velocity, acceleration, base reaction, stiffness, and damage index" are recorded.
13. A "ductility-based damage index" is calculated at each time step.
14. The effective stiffness degradation is monitored during the response.
15. Multiple simulations are executed by looping over ranges of "damping ratio, initial displacement, ductility, and over-strength".
16. For each simulation, "maximum response values" are extracted.
17. Damage index values are limited between "0% and 100%".
18. All results are stored for post-processing and statistical analysis.
19. The script reports completion of each simulation and total runtime.
20. Overall, the code provides a "parametric nonlinear dynamic assessment" of SDOF systems under Incremental Dynamic conditions.

This Python Script encodes the temperature-dependent degradation of structural steel per the experimental data.
It takes the steel temperature T (°C) and returns the reduction factors for elastic modulus (kE), yield strength (kFy), ultimate strength (kFu), and the associated ultimate strain (esu).
The stepwise logic reflects the distinct stages of material deterioration: a gradual decline up to 400 °C, a precipitous drop in stiffness and strength in the 400–500 °C range, and severe residual capacity beyond 700 °C.
These factors are critical for nonlinear thermo-mechanical analysis, enabling accurate prediction of load-bearing capacity, deflections, and failure modes during fire exposure.
Implementation follows the tabulated values in the referenced paper, providing a simple yet experimentally grounded fire-response material model for structural steel. 
"""
# BOOK: Dynamics of Structure and Foundation; A Unified Approach-1. Fundamentals-Indrajit Chowdhury-CRC Press-2009
# BOOK: Dynamics of structure and foundation; A Uniﬁed Approach-2. Applications-Indrajit Chowdhury-CRC Press-2009
# BOOK: Plastic Design and Second-Order Analysis of Steel Frames-W. F. Chen, I. Sohal-Springer New York-1995
# BOOK: Structural Stability-W. F. Chen'
# BOOK: Differential Equations for Engineers-Wei-Chau Xie-CAMBRIDGE-2010
# BOOK: Structural dynamics THEORY AND COMPUTATION MARIO BAZ-5 EDITION   
# BOOK: Dynamics of Structures in SI Units -  Anil Kumar Chopra
'https://share.google/TDN5O4eWmw5pH8zUH' 
# PAPER: Mechanical properties of structural steel at elevated temperatures and after cooling down
'https://www.researchgate.net/publication/227780727_Mechanical_properties_of_structural_steel_at_elevated_temperatures_and_after_cooling_down'
# PAPER: Post-fire structural material behaviour of high- and ultra-high-strength steels – Experimental study and aids for assessment for further use and potential reuse after fires
'https://www.sciencedirect.com/science/article/pii/S0950061824028848'    
# PAPER: POST-FIRE MECHANICAL PROPERTIES OF STRUCTURAL STEEL
'https://www.researchgate.net/publication/264545278_POST-FIRE_MECHANICAL_PROPERTIES_OF_STRUCTURAL_STEEL'
# WEBSITE:
'https://www.constructionspecifier.com/post-fire-structural-steel-evaluation/'
# WEBSITE:
'https://bews.es/services/structural-engineering-consultants/post-fire-analysis/'
#%%----------------------------------------------------
import time as TI
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
import OPENSEEES_HYSTERETICSM_FORCE_DISP_FUN as S08
import SALAR_MATH as S09
import PLOT_CONTOUR_3D_2D_FUN as S10
import EQUIVALENT_VISCOUS_DAMPING_RATIO_FUN as S055
import DAMAGE_INDEX_FUN as S066
#%%----------------------------------------------------
def POST_FIRE_SDOF_FOUR_ELEMENTS(THI, MI, KI, OSF, DPR, JJ, MAT_TYPE, TOTAL_MASS, ANAL_TYPE):
    # Initialize model
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)
    
    MAX_ITERATIONS = 5000   # Maximum number of iterations
    MAX_TOLERANCE = 1.0e-6  # Specified tolerance for convergence
    GMfact = 9.81 # [m/s^2] standard acceleration of gravity or standard acceleration 
    
    # [kg] Mass
    MASS = [0.40 * MI, # ELEMENT 01
            0.20 * MI, # ELEMENT 02
            0.30 * MI, # ELEMENT 03
            0.10 * MI] # ELEMENT 04  
    
    # Damping ratio
    DRi = [DPR, # ELEMENT 01
           DPR, # ELEMENT 02
           DPR, # ELEMENT 03
           DPR] # ELEMENT 04      
    
    # Define nodes
    ops.node(1, 0.0)  # Fixed base
    ops.node(2, 0.0)  # Mass node
        
    # Define boundary conditions
    ops.fix(1, 1)

    # Define masses
    ops.mass(2, np.sum(MASS))
    
    def STEEL_FIRE_PROPERTIES(T):
        # THERMAL EFFECTS OF STEEL MATERIAL STRESS-STRAIN
        # WRITTEN BY SALAR DElAVAR GHASHGHAEI (QASHQAI)
        # PAPER: Mechanical properties of structural steel at elevated temperatures and after cooling down
        'https://www.researchgate.net/publication/227780727_Mechanical_properties_of_structural_steel_at_elevated_temperatures_and_after_cooling_down'    
        if T <= 20:
            kE  = 1.00
            kFy = 1.00
            kFu = 1.00
            esu = 0.25
            
        elif 20 < T <= 100:
            kE  = 0.90
            kFy = 0.90
            kFu = 0.90
            esu = 0.277
            
        elif 100 < T <= 200:
            kE  = 0.80
            kFy = 0.80
            kFu = 0.80
            esu = 0.3125
            
        elif 200 < T <= 300:
            kE  = 0.70
            kFy = 0.70
            kFu = 0.70
            esu = 0.3571
            
        elif 300 < T <= 400:
            kE  = 0.60
            kFy = 0.60
            kFu = 0.60
            esu = 0.3814 
            
        elif 400 < T <= 500:
            kE  = 0.30
            kFy = 0.30
            kFu = 0.30
            esu = 0.4012              
        
        elif 500 < T <= 700:
            kE  = 0.13
            kFy = 0.23
            kFu = 0.35
            esu = 0.4124
    
        elif T > 700:
            kE  = 0.09
            kFy = 0.11
            kFu = 0.20
            esu = 0.45
    
        return kE, kFy, kFu, esu
    #%% FOUR DEGREES OF FREEDOM STRUCTURE, CALCULATE LATERAL STIFFNESS AND DAMPING 
    for II in range(0, 4):
    # THERMAL EFFECTS OF STEEL MATERIAL FORCE-DISPLACEMENT   
        # FIRE EFFECTS ON ALL COLUMNS
        kE, kFy, kFu, esu = STEEL_FIRE_PROPERTIES(THI)   # THERMAL EFFECTS [°C] 
        EsI = kE*KI                                      # [N/m] Spring Elastic Stiffness
        fyI = kFy*85000.0                                # [N] Yield Force of Structure
        fuI = kFu*(OSF*85000.0)                          # [N] Ultimate Force of Structure
        
        FY = fyI                                         # [N] Yield Force of Structure
        FU = fuI                                         # [N] Ultimate Force of Structure
        Ke = EsI                                         # [N/m] Spring Elastic Stiffness
        DY = fyI / Ke                                    # [m] Yield Displacement
        DSU = esu                                        # [m] Ultimate Displacement
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
            S08.OPENSEEES_HYSTERETICSM_FORCE_DISP_FUN(MAT_TAG, DP, FP, DN, FN, PLOT = False, X_LABEL='Displacement (m)', Y_LABEL='Force [N]', TITLE='FORCE-DISPLACEMENT CURVE')
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
        ops.element('zeroLength', II+1, 1, 2, '-mat',  MAT_TAG, MatTag_C, '-dir', 1, 1)  # DOF LATERAL SPRING

    center_node = 2
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
    DI = []
    ele_force = {
        1: [],  # AXIALFORCE-01
        }
    # Initialize lists for each node's displacement
    node_displacements = {
        1: [],  # DISP01
        2: [],  # DISP02
        }
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
            #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(1, 0.5*DR, DR, 0, 1)
            PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(1, PLOT=True)
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
            disp.append(ops.nodeDisp(center_node, 1))                       # DISPLACEMENT NODE 05 IN X DIR       
            # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
            #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(1, 0.5*DR, DR, 0, 1)
            PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(1, PLOT=True)
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
            #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(1, 0.5*DR, DR, 0, 1)
            PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(1, PLOT=True)
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
            #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(1, 0.5*DR, DR, 0, 1)
            PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(1, PLOT=True)
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
        u0 = -0.00010                      # [m] Initial displacement
        v0 = 0.000015                      # [m/s] Initial velocity
        a0 = 0.000065                      # [m/s^2] Initial acceleration
        IU = True                          # Free Vibration with Initial Displacement
        IV = True                          # Free Vibration with Initial Velocity
        IA = True                          # Free Vibration with Initial Acceleration
        duration = 20.0                    # [s] Analysis duration
        dt = 0.001                         # [s] Time step
        DR = 0.03                          # Damping Ratio
        
        disp_02, velo_02, accel_02 = [], [], []
        
        if IU == True:
            # Define initial displacment
            ops.setNodeDisp(center_node, 1, u0, '-commit') 
            # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/setNodeDisp.html
        if IV == True:
            # Define initial velocity
            ops.setNodeVel(center_node, 1, v0, '-commit')
            # INFO LINK: https://openseespydoc.readthedocs.io/en/stable/src/setNodeVel.html
        if IA == True:
            # Define initial  acceleration
            ops.setNodeAccel(center_node, 1, a0, '-commit')
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
            disp.append(ops.nodeDisp(center_node, 1))             # DISPLACEMENT NODE 05 IN X DIR 
            velo.append(ops.nodeVel(center_node, 1))              # VELOCITY NODE 05
            acc.append(ops.nodeAccel(center_node, 1))             # ACCELERATION NODE 05
            disp_02.append(ops.nodeDisp(2, 1)); velo_02.append(ops.nodeVel(2, 1)); accel_02.append(ops.nodeAccel(2, 1));
            stiffness.append(np.abs(reaction[-1] / disp[-1]))
            OMEGA.append(np.sqrt(stiffness[-1]/TOTAL_MASS))
            PERIOD.append((np.pi * 2) / OMEGA[-1]) 
            # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
            #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(1, 0.5*DR, DR, 0, 1)
            PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(1, PLOT=True)
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
        damping_ratio = S04.DAMPING_RATIO(disp)       # DAMAPING RATIO FROM DOF 05
        
        damping_ratio_02 = S04.DAMPING_RATIO(disp_02) # DAMAPING RATIO FROM DOF 02
        
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
        GMfact = JJ * 2 * 9.810            # [m/s²] standard acceleration of gravity or standard acceleration
        SSF_X = 200.0                      # Seismic Acceleration Scale Factor in X Direction
        SSF_Y = 200.0                      # Seismic Acceleration Scale Factor in Y Direction
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
            disp.append(ops.nodeDisp(center_node, 1))             # DISPLACEMENT NODE 05 IN X DIR 
            velo.append(ops.nodeVel(center_node, 1))              # VELOCITY NODE 05
            acc.append(ops.nodeAccel(center_node, 1))             # ACCELERATION NODE 05
            stiffness.append(np.abs(reaction[-1] / disp[-1]))
            OMEGA.append(np.sqrt(stiffness[-1]/TOTAL_MASS))
            PERIOD.append((np.pi * 2) / OMEGA[-1])
            # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
            #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(1, 0.5*DR, DR, 0, 1)
            PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(1, PLOT=True)
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
# STRUCTURAL SENSITIVITY ANALYSIS OF STRUCTURE

# Analysis Durations:
starttime = TI.process_time()

THERMAL_MAX, MASS_MAX, COE_MAX, KI_MAX, OSF_MAX, DPR_MAX = [], [], [], [], [], [] 
reaction_SEI_MAX, disp_SEI_MAX, velo_SEI_MAX, acc_SEI_MAX = [], [], [], []
DI_SEI_MAX, damping_ratio_SEI_MAX, PERIOD_MIN_SEI_MAX, PERIOD_MAX_SEI_MAX  = [], [], [], []
zeta_MAX = []

# STRUCTURAL THERMAL
THI = [50.0, 100.0, 200.0, 300.0, 400.0]
# STRUCTURAL MASS
TOTAL_MASS = 5000.0   # [kg] Total Mass of Structure
MI = [0.8*TOTAL_MASS, 0.9*TOTAL_MASS, 1.0*TOTAL_MASS, 1.1*TOTAL_MASS, 1.2*TOTAL_MASS]
# STRUCTURAL ELASTIC STIFFNEES
KI = [0.8*4500000.0, 0.9*4500000.0, 1.0*4500000.0, 1.1*4500000.0, 1.2*4500000.0]
#KI = [0.9*4500000.0, 1.0*4500000.0, 1.1*4500000.0]
# STRUCTURAL DUCTILITY RATIO
#DUCT = [2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0]
DUCT = [6.0, 6.5, 7.0, 7.5, 8.0]
# STRUCTURAL OVER-STRENGTH FACTOR
OSF = [1.05, 1.10, 1.15, 1.20, 1.25]
# STRUCTURAL DAMPING RATIO
#DPR = [0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10]
DPR = [0.01, 0.02, 0.03, 0.04, 0.05]
# STANDARD ACCELERATION OF GRAVITY OR STANDARD ACCELERATION SCALE FACTOR
#COE = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00]
COE = [0.10, 0.25, 0.50, 0.75, 0.90, 1.00]

MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'
ANAL_TYPE = 'SEISMIC'

II = 0
for ti in THI:
    for mi in MI:
        for ki in KI:
            for osf in OSF:
                for dpr in DPR:
                    for coe in COE:
                        DATA = POST_FIRE_SDOF_FOUR_ELEMENTS(ti, mi, ki, osf, dpr, coe, MAT_TYPE, TOTAL_MASS, ANAL_TYPE)
                                
                        (time_SEI, reaction_SEI, disp_SEI, velo_SEI, acc_SEI, DI_SEI,
                         ele_force_SEI, node_displacements_SEI,
                         stiffness_SEI, PERIOD_SEI, damping_ratio_SEI,
                         PERIOD_MIN_SEI, PERIOD_MAX_SEI) = DATA
                                
                        XLABEL = "Displacement [m]"
                        YLABEL = "Base Reaction [N]"
                        TITLE = "Equivalent viscous damping ratio - Seismic Hysteresis"
                        method = 1 # 
                        PLOT = False
                        zeta = S055.EQUIVALENT_VISCOUS_DAMPING_RATIO_FUN(disp_SEI, reaction_SEI, method, XLABEL, YLABEL, TITLE, PLOT)
                        print(f"Equivalent viscous damping ratio = {zeta:.4f}")
                            
                        THERMAL_MAX.append(ti)
                        MASS_MAX.append(mi)
                        COE_MAX.append(coe)
                        KI_MAX.append(ki)
                        OSF_MAX.append(osf)
                                
                        reaction_SEI_MAX.append(np.max(np.abs(reaction_SEI))) 
                        disp_SEI_MAX.append(np.max(np.abs(disp_SEI))) 
                        velo_SEI_MAX.append(np.max(np.abs(acc_SEI)))  
                        acc_SEI_MAX.append(np.max(np.abs(reaction_SEI))) 
                        DI_SEI_MAX.append(np.max(np.abs(DI_SEI))) 
                        PERIOD_MAX_SEI_MAX.append(np.max(np.abs(PERIOD_MAX_SEI)))
                        zeta_MAX.append(zeta)
                        II = II + 1
                        print(f'\n STEP: {II} \n COE: {coe} - KI: {ki} - DPR: {dpr} - OSF: {osf} \n')
else:
    print('\n\n Analysis completed successfully')    

totaltime = TI.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n') 
#%%-------------------------------------------------------------------
# PLOT 2D AND 3D CONTOUR 
"""        
X, Y, Z = COE_MAX, DUCT_MAX, OSF_MAX 
XLABEL, YLABEL, ZLABEL = 'GRAVITY COEFFICIENT', 'DUCTILITY RATIO [m/m]', 'OVER-STRENGTH FACTOR [N/N]'               
S10.PLOT_CONTOUR_3D_2D_FUN(111, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = COE_MAX, DUCT_MAX, KI_MAX 
XLABEL, YLABEL, ZLABEL = 'GRAVITY COEFFICIENT', 'DUCTILITY RATIO [m/m]', 'ELASTIC STIFFNESS [N/m]'               
S10.PLOT_CONTOUR_3D_2D_FUN(112, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = COE_MAX, DUCT_MAX, DPR_MAX 
XLABEL, YLABEL, ZLABEL = 'GRAVITY COEFFICIENT', 'DUCTILITY RATIO [m/m]', 'DAMPING RATIO [%]'               
S10.PLOT_CONTOUR_3D_2D_FUN(113, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = DUCT_MAX, OSF_MAX, DPR_MAX 
XLABEL, YLABEL, ZLABEL = 'DUCTILITY RATIO [m/m]', 'OVER-STRENGTH FACTOR [N/N]', 'DAMPING RATIO [%]'               
S10.PLOT_CONTOUR_3D_2D_FUN(114, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = DUCT_MAX, OSF_MAX, DI_SEI_MAX 
XLABEL, YLABEL, ZLABEL = 'DUCTILITY RATIO [m/m]', 'OVER-STRENGTH FACTOR [N/N]', 'STRUCTURAL DAMAGE INDEX [%]'               
S10.PLOT_CONTOUR_3D_2D_FUN(118, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = KI_MAX, DPR_MAX, disp_SEI_MAX 
XLABEL, YLABEL, ZLABEL = 'ELASTIC STIFFNESS [N/m]', 'DAMPING RATIO [%]', 'MAX. DISP. [m]'               
S10.PLOT_CONTOUR_3D_2D_FUN(119, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = KI_MAX, DPR_MAX, velo_SEI_MAX 
XLABEL, YLABEL, ZLABEL = 'ELASTIC STIFFNESS [N/m]', 'DAMPING RATIO [%]', 'MAX. VELOCITY [m/s]'               
S10.PLOT_CONTOUR_3D_2D_FUN(120, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = KI_MAX, DPR_MAX, acc_SEI_MAX 
XLABEL, YLABEL, ZLABEL = 'ELASTIC STIFFNESS [N/m]', 'DAMPING RATIO [%]', 'MAX. ACCELERATION [m/s^2]'               
S10.PLOT_CONTOUR_3D_2D_FUN(121, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = KI_MAX, DPR_MAX, DI_SEI_MAX 
XLABEL, YLABEL, ZLABEL = 'ELASTIC STIFFNESS [N/m]', 'DAMPING RATIO [%]', 'STRUCTURAL DAMAGE INDEX [%]'               
S10.PLOT_CONTOUR_3D_2D_FUN(122, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = KI_MAX, DPR_MAX, PERIOD_MAX_SEI_MAX
XLABEL, YLABEL, ZLABEL = 'ELASTIC STIFFNESS [N/m]', 'DAMPING RATIO [%]', 'STRUCTURAL PERIOD [s]'               
S10.PLOT_CONTOUR_3D_2D_FUN(123, X, Y, Z, XLABEL, YLABEL, ZLABEL)
"""
X, Y, Z = KI_MAX, acc_SEI_MAX, zeta_MAX
XLABEL, YLABEL, ZLABEL = 'ELASTIC STIFFNESS [N/m]', 'MAX. ACCELERATION [m/s^2]', 'Equivalent viscous damping ratio [%]'               
S10.PLOT_CONTOUR_3D_2D_FUN(124, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = KI_MAX, disp_SEI_MAX, zeta_MAX
XLABEL, YLABEL, ZLABEL = 'ELASTIC STIFFNESS [N/m]', 'MAX. VELOCITY [m/s]', 'Equivalent viscous damping ratio [%]'               
S10.PLOT_CONTOUR_3D_2D_FUN(125, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = KI_MAX, velo_SEI_MAX , zeta_MAX
XLABEL, YLABEL, ZLABEL = 'ELASTIC STIFFNESS [N/m]', 'MAX. DISP. [m]', 'Equivalent viscous damping ratio [%]'               
S10.PLOT_CONTOUR_3D_2D_FUN(126, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = KI_MAX, acc_SEI_MAX, PERIOD_MAX_SEI_MAX
XLABEL, YLABEL, ZLABEL = 'ELASTIC STIFFNESS [N/m]', 'MAX. ACCELERATION [m/s^2]', 'STRUCTURAL PERIOD [s]'               
S10.PLOT_CONTOUR_3D_2D_FUN(124, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = KI_MAX, disp_SEI_MAX, PERIOD_MAX_SEI_MAX
XLABEL, YLABEL, ZLABEL = 'ELASTIC STIFFNESS [N/m]', 'MAX. VELOCITY [m/s]', 'STRUCTURAL PERIOD [s]'               
S10.PLOT_CONTOUR_3D_2D_FUN(125, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = KI_MAX, velo_SEI_MAX , PERIOD_MAX_SEI_MAX
XLABEL, YLABEL, ZLABEL = 'ELASTIC STIFFNESS [N/m]', 'MAX. DISP. [m]', 'STRUCTURAL PERIOD [s]'               
S10.PLOT_CONTOUR_3D_2D_FUN(126, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = KI_MAX, MASS_MAX , PERIOD_MAX_SEI_MAX
XLABEL, YLABEL, ZLABEL = 'ELASTIC STIFFNESS [N/m]', 'MASS [kg]', 'STRUCTURAL PERIOD [s]'               
S10.PLOT_CONTOUR_3D_2D_FUN(128, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = KI_MAX, MASS_MAX , zeta_MAX 
XLABEL, YLABEL, ZLABEL = 'ELASTIC STIFFNESS [N/m]', 'MASS [kg]', 'Equivalent viscous damping ratio [%]'               
S10.PLOT_CONTOUR_3D_2D_FUN(129, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = KI_MAX, MASS_MAX , disp_SEI_MAX 
XLABEL, YLABEL, ZLABEL = 'ELASTIC STIFFNESS [N/m]', 'MASS [kg]', 'MAX. DISP. [m]'               
S10.PLOT_CONTOUR_3D_2D_FUN(130, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = KI_MAX, MASS_MAX , velo_SEI_MAX 
XLABEL, YLABEL, ZLABEL = 'ELASTIC STIFFNESS [N/m]', 'MASS [kg]', 'MAX. VELOCITY [m/s]'               
S10.PLOT_CONTOUR_3D_2D_FUN(131, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = KI_MAX, MASS_MAX , acc_SEI_MAX 
XLABEL, YLABEL, ZLABEL = 'ELASTIC STIFFNESS [N/m]', 'MASS [kg]', 'MAX. ACCELERATION [m/s^2]'               
S10.PLOT_CONTOUR_3D_2D_FUN(132, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = COE_MAX, MASS_MAX , PERIOD_MAX_SEI_MAX
XLABEL, YLABEL, ZLABEL = 'GRAVITY COEFFICIENT', 'MASS [kg]', 'STRUCTURAL PERIOD [s]'               
S10.PLOT_CONTOUR_3D_2D_FUN(133, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = COE_MAX, MASS_MAX , zeta_MAX 
XLABEL, YLABEL, ZLABEL = 'GRAVITY COEFFICIENT', 'MASS [kg]', 'Equivalent viscous damping ratio [%]'               
S10.PLOT_CONTOUR_3D_2D_FUN(134, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = COE_MAX, MASS_MAX , disp_SEI_MAX 
XLABEL, YLABEL, ZLABEL = 'GRAVITY COEFFICIENT', 'MASS [kg]', 'MAX. DISP. [m]'               
S10.PLOT_CONTOUR_3D_2D_FUN(135, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = COE_MAX, MASS_MAX , velo_SEI_MAX 
XLABEL, YLABEL, ZLABEL = 'GRAVITY COEFFICIENT', 'MASS [kg]', 'MAX. VELOCITY [m/s]'               
S10.PLOT_CONTOUR_3D_2D_FUN(136, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = COE_MAX, MASS_MAX , acc_SEI_MAX 
XLABEL, YLABEL, ZLABEL = 'GRAVITY COEFFICIENT', 'MASS [kg]', 'MAX. ACCELERATION [m/s^2]'               
S10.PLOT_CONTOUR_3D_2D_FUN(137, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = THERMAL_MAX, zeta_MAX, disp_SEI_MAX 
XLABEL, YLABEL, ZLABEL = 'THERMAL EFFECTS [°C]', 'Equivalent viscous damping ratio [%]', 'MAX. DISP. [m]'               
S10.PLOT_CONTOUR_3D_2D_FUN(138, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = THERMAL_MAX, zeta_MAX, velo_SEI_MAX 
XLABEL, YLABEL, ZLABEL = 'THERMAL EFFECTS [°C]', 'Equivalent viscous damping ratio [%]', 'MAX. VELOCITY [m/s]'               
S10.PLOT_CONTOUR_3D_2D_FUN(139, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = THERMAL_MAX, zeta_MAX, acc_SEI_MAX 
XLABEL, YLABEL, ZLABEL = 'TTHERMAL EFFECTS [°C]', 'Equivalent viscous damping ratio [%]', 'MAX. ACCELERATION [m/s^2]'               
S10.PLOT_CONTOUR_3D_2D_FUN(140, X, Y, Z, XLABEL, YLABEL, ZLABEL)
#%%-------------------------------------------------------------------
X, HISTO_COLOR, LABEL = disp_SEI_MAX, 'lime', 'DISPLACEMENT FROM IDA ANALYSIS'
S09.HISROGRAM_BOXPLOT(X, HISTO_COLOR, LABEL)
X, HISTO_COLOR, LABEL = velo_SEI_MAX, 'pink', 'VELOCITY FROM IDA ANALYSIS'
S09.HISROGRAM_BOXPLOT(X, HISTO_COLOR, LABEL)
X, HISTO_COLOR, LABEL = acc_SEI_MAX, 'cyan', 'ACCELERATION FROM IDA ANALYSIS'
S09.HISROGRAM_BOXPLOT(X, HISTO_COLOR, LABEL)
X, HISTO_COLOR, LABEL = reaction_SEI_MAX, 'purple', 'REACTION FROM IDA ANALYSIS'
S09.HISROGRAM_BOXPLOT(X, HISTO_COLOR, LABEL)
X, HISTO_COLOR, LABEL = PERIOD_MAX_SEI_MAX , 'green', 'STRUCTURAL PERIOD [s] FROM IDA ANALYSIS'
S09.HISROGRAM_BOXPLOT(X, HISTO_COLOR, LABEL)
X, HISTO_COLOR, LABEL = DI_SEI_MAX, 'orange', 'DUCTILITY DAMAGE INDEX FROM IDA ANALYSIS'
S09.HISROGRAM_BOXPLOT(X, HISTO_COLOR, LABEL)
X, HISTO_COLOR, LABEL = zeta_MAX , 'green', 'EQUIVALENT VISCOUS DAMPING RATIO [%] FROM IDA ANALYSIS'
S09.HISROGRAM_BOXPLOT(X, HISTO_COLOR, LABEL)
X, HISTO_COLOR, LABEL = THERMAL_MAX , 'brown', 'THERMAL EFFFECTS [°C] FROM IDA ANALYSIS'
S09.HISROGRAM_BOXPLOT(X, HISTO_COLOR, LABEL)
#%%-------------------------------------------------------------------
import pandas as pd
DATA_TOTAL = { 
    'DISPLACEMENT': disp_SEI_MAX,
    'VELOCITY': velo_SEI_MAX,
    'ACCELERATION ': acc_SEI_MAX,
    'REACTION': reaction_SEI_MAX, 
    #'DUCTILITY DAMAGE INDEX': DI_SEI_MAX,
    'EQULIVALENT VISCOUS DAMPING RATIO': zeta_MAX,
}
# Convert to DataFrame
results_df = pd.DataFrame(DATA_TOTAL)
# CORRELATION MATRIX
S09.PLOT_HEATMAP(results_df)
#%%-------------------------------------------------------------------
# Plots a heatmap of sensitivity coefficients (correlation or SRC) between inputs X and outputs Y.
import SENSITIVITY_HEATMAP_FUN as S099
X = np.column_stack([THERMAL_MAX, KI_MAX, MASS_MAX, OSF_MAX])          
Y = np.column_stack([velo_SEI_MAX, acc_SEI_MAX, PERIOD_MAX_SEI_MAX, zeta_MAX])  
X_LABELS = ['THERMAL EFFECTS [°C]', 'ELASTIC STIFFNESS [N/m]', 'STRUCTURAL MASS [kg]', 'OVER-STRENGTH FACTOR [N/N]']
Y_LABELS = ['VELOCITY [m/s]', 'ACCELERATION [m/s^2]', 'STRUCTURAL PERIOD [s]', 'EQUIVALENT VISCOUS DAMPING RATIO [%]']
coeffs = S099.SENSITIVITY_HEATMAP_FUN(X, Y, X_LABELS, Y_LABELS, method='pearson')
#%%-------------------------------------------------------------------
# ANOVA with automatic binning of continuous predictors.
import ANOVA_SENSITIVITY_FUN as S100

df_sens = pd.DataFrame({
    "THERMAL": THERMAL_MAX,
    "MASS":    MASS_MAX,
    "COE":     COE_MAX,
    "KI":      KI_MAX,
    "OSF":     OSF_MAX,
    "ZETA":    zeta_MAX,
    "DISP":    disp_SEI_MAX,
    "VEL":     velo_SEI_MAX,
    "ACC":     acc_SEI_MAX,
    "PERIOD":  PERIOD_MAX_SEI_MAX,
})

param_list = ["THERMAL", "KI", "OSF", "COE", "DISP", "VEL", "ACC", "PERIOD"]

# ANOVA – main effects only
anova_table, _ = S100.ANOVA_SENSITIVITY_FUN(
    df_sens, output_col="ZETA",
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
    df_sens, output_col="ZETA",
    param_cols=param_list,
    n_bins=4,
    include_interactions=True,
    plot=True,
)
plt.show()
print("\nANOVA table (with interactions):")
print(anova_table_inter.round(4))
#%%-------------------------------------------------------------------
