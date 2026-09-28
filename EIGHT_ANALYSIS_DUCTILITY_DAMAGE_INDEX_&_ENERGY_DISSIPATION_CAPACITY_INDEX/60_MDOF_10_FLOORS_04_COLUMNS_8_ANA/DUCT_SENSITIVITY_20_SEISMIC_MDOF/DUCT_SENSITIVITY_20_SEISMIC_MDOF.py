###########################################################################################################
#                   >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                      #
#        SENSITIVITY ANALYSIS OF COLUMN DUCTILITY RATIO  AND OPENSEES VIA 20 SEISMIC GROUND MOTIONS       #
#---------------------------------------------------------------------------------------------------------#
#                EQUIVALENT VISCOUS DAMPING RATIO: xi_eq = 100 * E_d / (4 * pi * E_s)                     #
#---------------------------------------------------------------------------------------------------------#
#          ENERGY DISSIPATION CAPACITY INDEX = 100 * E_d(earthquake) / E_d(cyclic displacement)           #
#---------------------------------------------------------------------------------------------------------#
#                                                  P(t) = P0 sin(wt)                                      #
#                                           P(t) = P0 exp(-0.05wt) sin(wt)                                #
#---------------------------------------------------------------------------------------------------------#
#         ASSESSMENT OF DUCTILITY DAMAGE INDICES FOR STRUCTURAL ELEMENTS AND SYSTEMS AND EVALUATION       #
#                                        OF ENERGY DISSIPATION CAPACITY INDEX                             #
#---------------------------------------------------------------------------------------------------------#
#                  THIS PYTHON SCRIPT IS WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                    #
#                                   EMAIL: salar.d.ghashghaei@gmail.com                                   #
###########################################################################################################
"""
Assume that a multi-degree-of-freedom structure is subjected to nonlinear dynamic analysis under
 twenty different ground motion records. To investigate the influence of column ductility on
 seismic performance, a sensitivity analysis is performed by systematically varying each column’s
 ductility ratio from 20 to 50 in 21 steps. For each ductility level, a pushover analysis provides
 the capacity curve and equivalent SDOF parameters, a cyclic displacement analysis establishes a
 reference hysteretic energy, and twenty nonlinear time-history analyses are conducted using
 different ground motions. From these analyses, the Energy Dissipation Capacity Index (EDCI)
 and other response metrics such as maximum displacement, velocity, acceleration, damage index,
 over-strength factor, ductility ratio, and equivalent viscous damping ratio are computed.
 The median response over the twenty ground motions is then evaluated for each ductility level.
 Finally, the results are assessed through trend plots, 3D contour surfaces, correlation heatmaps,
 Random Forest regression, and ANOVA sensitivity analysis to identify the optimal column ductility
 ratio that maximizes energy dissipation capacity while maintaining structural safety.

 
Nonlinear Seismic Performance Assessment of a MDOF:
An OpenSeesPy Framework for Material and Geometric Nonlinearity Under Static, Cyclic, and Earthquake Loading

This OpenSeesPy script performs rigorous nonlinear static and dynamic analysis of a MDOF
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
"""

==========================================================
SENSITIVITY ANALYSIS BY CHANGING COLUMN DUCTILITY RATIO
==========================================================

1. Sets MAT_TYPE = 'INELASTIC' and sweeps the COLUMN DUCTILITY RATIO from 20.0 to 50.0
   in 21 linear steps (DUCT_MIN -> DUCT_MAX).  COL_DUCT stores each swept value.
2. Declares ~20 accumulator lists: raw per-record responses, plus *_MED lists that
   will hold one MEDIAN value per ductility level.
3. Starts a CPU timer.  Outer loop `for JJ in range(0, 21)` computes DUCT and resets
   the per-ductility raw accumulators (DISP, VELO, ACC, EDCI, DII, OMEGA, MU, RR, EDVR).
4. STEP 1 -- PUSHOVER: MDOF(DUCT, MAT_TYPE, TOTAL_MASS, 'PUSHOVER', II) returns the
   capacity curve (disp_PUSH, reaction_PUSH) and the equivalent-SDOF system
   (mass, stiffness, period) via displacement-based pushover reduction.
5. Stores SDOF_ef_* arrays and PERIOD_PUSH = max(PERIOD_MAX_PUSH)  (secant period).
6. STEP 2 -- CYCLIC_DISPLACEMENT: benchmark constant-amplitude hysteresis loop.
   S055.EQULIVALENT_VISCOUS_DAMPING_RATIO_FUN returns zeta_CP (Jacobsen area method).
7. STEP 3 -- SEISMIC LOOP `for II in range(40, 60)`: 20 nonlinear time-history runs.
8. Each run unpacks time/reaction/disp/velo/acc/DI/stiffness/period arrays from MDOF().
9. S12.DUCTILITY_DAMAGE_INDEX_FUN fits the pushover curve and returns Park-Ang style
   DIx, over-strength Omega_0, ductility mu, and behaviour coefficient R.
10. S10.ENERGY_DISSIPATION_CAPACITY_INDEX compares seismic vs cyclic dissipated energy
   -> EDCI (%).  Peak DISP/VELO/ACC also appended.
11. Ratio metrics appended: PERIOD_RATIO, DI_RATIO, REACTION_RATIO, DISP_RATIO
   (each = SEISMIC / PUSHOVER).
12. S055 again computes zeta_SEI for the seismic loop -> EDVR and EDVR_RATIO = zeta_SEI/zeta_CP.
13. STEP 4 -- after 20 records, np.median() reduces each raw list to ONE robust value
   (DISP_MED, DI_RATIO_MED, EDCI_MED, EDVR_MED, ...), the design-point response.
14. Outer loop ends after 21 design points.  Total CPU time is printed.
15. PLOTTING: ~17 calls to S01.PLOT_SCATTER(COL_DUCT, Y, ...) draw scatter + polynomial
   fit of order 1/3/7 for every response vs. the ductility ratio.
16. 5 calls to S11.PLOT_CONTOUR_3D_2D_FUN build 2D contours + 3D surfaces pairing two
   responses (e.g. DISP_RATIO vs DI_RATIO) against the ductility ratio axis.
17. MACHINE LEARNING: builds a pandas DataFrame with features DISP, VEL, ACC,
   REACTION_RATIO, EDVR_RATIO, DI_RATIO, then calls S01.RANDOM_FOREST(df) to train a
   RandomForestClassifier (safe/unsafe) + Regressor (safety likelihood) and print
   accuracy, MSE, R2, and feature importances.
18. S01.PLOT_HEATMAP(df) draws the Pearson correlation matrix of the six features.
19. S099.SENSITIVITY_HEATMAP_FUN(X=[DISP,VEL,ACC], Y=[EDVR_RATIO,DUCT,DI_RATIO],
   method='pearson') returns a 3x3 sensitivity-coefficient heatmap.
20. S100.ANOVA_SENSITIVITY_FUN(df_sens, output_col="DISP", param_cols=[VEL,ACC,
   EDVR_RATIO,DI_RATIO], n_bins=4) bins continuous predictors and runs classical ANOVA.
21. Two ANOVA runs: (a) MAIN EFFECTS only, (b) MAIN + TWO-WAY INTERACTIONS; each prints
   a SumSq / df / F / p-value table and shows the diagnostic plot.

22. Engineering reading: plots reveal whether higher column ductility reduces damage
   (DII), damping (EDVR), energy index (EDCI) and peak demand -- the core question of
   performance-based seismic design.

23. ML + correlation + ANOVA together rank which EDP (disp, vel or acc) drives damage.
24. The whole file is the FORWARD problem of a PBEE workflow; the companion Newton-
   Raphson script solves the INVERSE problem (find DUCT for target EDCI).
=====================================================================================

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
import pandas as pd
import openseespy.opensees as ops
import matplotlib.pyplot as plt 
import SALAR_MATH as S01
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
import DISSIPATED_ENERGY_WITH_PLOT_FUN as S10
import PLOT_CONTOUR_3D_2D_FUN as S11
import DUCTILITY_DAMAGE_INDEX_FUN as S12
#%%----------------------------------------------------
def MDOF(DUCT, MAT_TYPE, TOTAL_MASS, ANAL_TYPE, i):
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
    FUi = [1.25 * FYi[0],   # COLUMN 01
           1.12 * FYi[1],   # COLUMN 02
           1.20 * FYi[2],   # COLUMN 03
           1.10 * FYi[3]]   # COLUMN 04
    
    # ELASTIC STIFFNESS [N/m]
    Ke = 4500000.0          # [N/m] Elastic Stiffness
    Kei = [Ke * 1.20,       # COLUMN 01
           Ke * 0.90,       # COLUMN 02
           Ke * 1.00,       # COLUMN 03
           Ke * 0.70]       # COLUMN 04
    """ 
    # ULTIMATE DISPLACEMENT [m]
    DSUi = [0.32,           # COLUMN 01
            0.34,           # COLUMN 02
            0.30,           # COLUMN 03
            0.36]           # COLUMN 04
    """
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
    DRi = [0.03, # DOF 01 
           0.03, # DOF 02
           0.03, # DOF 03
           0.03, # DOF 04
           0.03, # DOF 05
           0.03, # DOF 06
           0.03, # DOF 07 
           0.03, # DOF 08
           0.03, # DOF 09
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
            DSU = DY * DUCT                                  # [m] Ultimate Displacement
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
                ops.uniaxialMaterial('Elastic', MAT_TAG, Ke)             # TESNSION AND COMPRESSION IS SAME VALUES
                #ops.uniaxialMaterial('Elastic', MAT_TAG, Ke ,0.0, 0.5*Ke) # TESNSION AND COMPRESSION IS NOT SAME VALUES
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
        DMAX = -1.0 * DSU    # [m] Max. Displacement
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
            
        #%%------------------------------------------------------- 
        # Run the file loading effective properties
        #exec(open("COMPUTE_EFFECTIVE_PROPERTIES_FUN_PUSHOVER.py").read())
        #% EQUIVALENT SDOF SYSTEM DERIVATION VIA DISPLACEMENT-BASED SEISMIC DESIGN PROCEDURE WITH PUSHOVER ANALYSIS
        # Change MDOF to SDOF System with Displacement Based Design Concept
        # THIS PYTHON SCRIPT IS WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
        """
        This script implements a displacement-based pushover transformation,
         converting a multi-degree-of-freedom (MDOF) system into an equivalent
         single-degree-of-freedom (SDOF) system for seismic assessment.
         It calculates effective modal properties—displacement, mass, and
         stiffness—by weighting element forces and nodal displacements according
         to a presumed deformed shape.
         The derived effective period provides a simplified dynamic characteristic
         for performance-based engineering. The visualizations effectively track the
         evolution of these equivalent parameters throughout the nonlinear analysis steps.
        """
        #--------------------------------------------------------------------------- 
        displacement_X_1 = np.array(list(node_displacements.values())[1])   # DOF 02    
        displacement_X_2 = np.array(list(node_displacements.values())[2])   # DOF 03 
        displacement_X_3 = np.array(list(node_displacements.values())[3])   # DOF 04 
        displacement_X_4 = np.array(list(node_displacements.values())[4])   # DOF 05 
        displacement_X_5 = np.array(list(node_displacements.values())[5])   # DOF 06    
        displacement_X_6 = np.array(list(node_displacements.values())[6])   # DOF 07 
        displacement_X_7 = np.array(list(node_displacements.values())[7])   # DOF 08 
        displacement_X_8 = np.array(list(node_displacements.values())[8])   # DOF 09 
        displacement_X_9 = np.array(list(node_displacements.values())[9])   # DOF 10 
        displacement_X_10 = np.array(list(node_displacements.values())[10])   # DOF 11 

        ele_force_01 = np.array(list(ele_force.values())[0])   # ELEMENT 01 
        ele_force_02 = np.array(list(ele_force.values())[1])   # ELEMENT 02 
        ele_force_03 = np.array(list(ele_force.values())[2])   # ELEMENT 03 
        ele_force_04 = np.array(list(ele_force.values())[3])   # ELEMENT 04 
        ele_force_05 = np.array(list(ele_force.values())[4])   # ELEMENT 05 
        ele_force_06 = np.array(list(ele_force.values())[5])   # ELEMENT 06 
        ele_force_07 = np.array(list(ele_force.values())[6])   # ELEMENT 07 
        ele_force_08 = np.array(list(ele_force.values())[7])   # ELEMENT 08 
        ele_force_09 = np.array(list(ele_force.values())[8])   # ELEMENT 09 
        ele_force_10 = np.array(list(ele_force.values())[9])   # ELEMENT 10 

        STIFF_X_1 = np.abs(ele_force_01 / displacement_X_1)
        STIFF_X_2 = np.abs(ele_force_02 / displacement_X_2)
        STIFF_X_3 = np.abs(ele_force_03 / displacement_X_3)
        STIFF_X_4 = np.abs(ele_force_04 / displacement_X_4)
        STIFF_X_5 = np.abs(ele_force_05 / displacement_X_5)
        STIFF_X_6 = np.abs(ele_force_06 / displacement_X_6)
        STIFF_X_7 = np.abs(ele_force_07 / displacement_X_7)
        STIFF_X_8 = np.abs(ele_force_08 / displacement_X_8)
        STIFF_X_9 = np.abs(ele_force_09 / displacement_X_9)
        STIFF_X_10 = np.abs(ele_force_10 / displacement_X_10)

        MX2 = (MASS[0] * np.square(displacement_X_1) + 
               MASS[1] * np.square(displacement_X_2) + 
               MASS[2] * np.square(displacement_X_3) +
               MASS[3] * np.square(displacement_X_4) +
               MASS[4] * np.square(displacement_X_5) + 
               MASS[5] * np.square(displacement_X_6) +
               MASS[6] * np.square(displacement_X_7) + 
               MASS[7] * np.square(displacement_X_8) + 
               MASS[8] * np.square(displacement_X_9) +
               MASS[9] * np.square(displacement_X_10))

        MX = (MASS[0] * np.array(displacement_X_1) + 
              MASS[1] * np.array(displacement_X_2) + 
              MASS[2] * np.array(displacement_X_3) + 
              MASS[3] * np.array(displacement_X_4) +
              MASS[4] * np.array(displacement_X_5) + 
              MASS[5] * np.array(displacement_X_6) + 
              MASS[6] * np.array(displacement_X_7) +
              MASS[7] * np.array(displacement_X_8) + 
              MASS[8] * np.array(displacement_X_9) + 
              MASS[9] * np.array(displacement_X_10))
        EFFECTIVE_DISP_X = MX2 / np.abs(MX)

        EFFECTIVE_MASS_X = np.abs(MX) / EFFECTIVE_DISP_X


        # Effective Stiffness
        KX = (np.array(STIFF_X_1) * np.array(displacement_X_1) + 
              np.array(STIFF_X_2) * np.array(displacement_X_2) + 
              np.array(STIFF_X_3) * np.array(displacement_X_3) + 
              np.array(STIFF_X_4) * np.array(displacement_X_4) + 
              np.array(STIFF_X_5) * np.array(displacement_X_5) + 
              np.array(STIFF_X_6) * np.array(displacement_X_6) + 
              np.array(STIFF_X_7) * np.array(displacement_X_7) + 
              np.array(STIFF_X_8) * np.array(displacement_X_8) + 
              np.array(STIFF_X_9) * np.array(displacement_X_9) + 
              np.array(STIFF_X_10) * np.array(displacement_X_10))

        EFFECTIVE_STIFF_X = np.abs(KX) / EFFECTIVE_DISP_X


        # Effective Period
        EFFECTIVE_PERIOD_X = 2 * np.pi / np.sqrt(EFFECTIVE_STIFF_X/EFFECTIVE_MASS_X)


        MED_DISP = np.median(EFFECTIVE_DISP_X)
        MED_MASS = np.median(EFFECTIVE_MASS_X)
        MED_STIFF = np.median(EFFECTIVE_STIFF_X)
        MED_PERIOD = np.median(EFFECTIVE_PERIOD_X)

        print('Median Effective Displacement:           ', MED_DISP)
        print('Median Effective Mass:                   ', MED_MASS)
        print('Median Effective Stiffness:              ', MED_STIFF)
        print('Median Effective Period:                 ', MED_PERIOD)

        # Create a figure with two subplots (Effective Displacement and Effective Mass)
        fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(12, 10))

        ax1.plot(EFFECTIVE_DISP_X, color='green', linewidth=3)
        ax1.set_title(f'Effective Displacement - Median: {np.median(EFFECTIVE_DISP_X): .5f}')
        ax1.set_xlabel('Step')
        ax1.set_ylabel('Effective Displacement')
        #ax1.legend(loc='upper right')
        ax1.grid(True)

        ax2.plot(EFFECTIVE_MASS_X, color='magenta', linewidth=3)
        ax2.set_title(f'Effective Mass - Median: {np.median(EFFECTIVE_MASS_X): .5f}')
        ax2.set_xlabel('Step')
        ax2.set_ylabel('Effective Mass')
        #ax2.legend(loc='upper right')
        ax2.grid(True)

        ax3.plot(EFFECTIVE_STIFF_X, color='cyan', linewidth=3)
        ax3.set_title(f'Effective Stiffness - Median: {np.median(EFFECTIVE_STIFF_X): .5f}')
        ax3.set_xlabel('Step')
        ax3.set_ylabel('Effective Stiffness')
        #ax3.legend(loc='upper right')
        ax3.grid(True)

        ax4.plot(EFFECTIVE_PERIOD_X, color='black', linewidth=3)
        ax4.set_title(f'Effective Period - Median: {np.median(EFFECTIVE_PERIOD_X): .5f}')
        ax4.set_xlabel('Step')
        ax4.set_ylabel('Effective Period')
        #ax4.legend(loc='upper right')
        ax4.semilogy()
        ax4.grid(True)

        plt.tight_layout()
        plt.show()
        
        SDOF_EFFE_DISP = MED_DISP 
        SDOF_EFFE_MASS = MED_MASS
        SDOF_EFFE_STIFF = MED_STIFF 
        SDOF_EFFE_PERIOD = MED_PERIOD
        #%%-------------------------------------------------------
        
        DATA = (reaction, disp, DI,
                ele_force, node_displacements,
                np.array(PERIOD_MIN), np.array(PERIOD_MAX),
                SDOF_EFFE_DISP, SDOF_EFFE_MASS, SDOF_EFFE_STIFF, SDOF_EFFE_PERIOD)
        
        return  DATA
    
    if ANAL_TYPE == 'CYCLIC_DISPLACEMENT': # STATIC TIME-HISTORY ANALYSIS
        IDctrlDOF = 1   # 1: Horizental Dispalcement - 2: Vertical Dispalcement
        DMAX = -1.0 * DSU     # [m] Max. Displacement
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
        duration = 10.0                    # [s] Analysis duration
        dt = 0.01                          # [s] Time step
        GMfact = 9.810                     # [m/s²]standard acceleration of gravity or standard acceleration
        SSF_X = 50.10                      # Seismic Acceleration Scale Factor in X Direction
        SSF_Y = 50.10                      # Seismic Acceleration Scale Factor in Y Direction
        iv0_X = 0.0005                     # [m/s] Initial velocity applied to the node  in X Direction
        iv0_Y = 0.0005                     # [m/s] Initial velocity applied to the node  in Y Direction
        st_iv0 = 0.0                       # [s] Initial velocity applied starting time
        SEI = 'X'                          # Seismic Direction
        DR = np.mean(DRi)                  # Damping ratio
        
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
TOTAL_MASS = 5_000_000.0   # [kg] Total Mass of Structure
#%%----------------------------------------------------
MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'

# --------------------------------------------------------------------------------------
# SENSITIVITY ANALYSIS BY CHANGING EACH COLUMN DUCTILITY RATIO
# --------------------------------------------------------------------------------------


import time as TI
import numpy as np

DUCT_MIN = 20.0          # MIN. COLUMN'S DUCTILITY RATIO
DUCT_MAX = 50.0          # MAX. COLUMN'S DUCTILITY RATIO

COL_DUCT       = []   # the swept variable itself  (21 values)
PERIOD_RATIO   = []   # T_max(seismic) / T_max(pushover)
DI_RATIO       = []   # damage index seismic / damage index pushover
REACTION_RATIO = []   # base reaction ratio
DISP_RATIO     = []   # displacement ratio
EDVR_RATIO     = []   # equivalent viscous damping ratio  (seismic / cyclic)
DISP, VELO, ACC = [], [], []   # max displacement / velocity / acceleration
EDCI           = []   # Energy Dissipation Capacity Index
DII, OMEGA, MU, RR = [], [], [], []   # damage index, over-strength, ductility, behaviour coeff.
EDVR           = []   # equivalent viscous damping ratio (seismic)
PERIOD_PUSH    = []   # max period from the pushover analysis
PERIOD_DYN     = []   # max period from the dynamic (seismic) analysis

# SDOF-equivalent quantities derived from the pushover curve (displacement-based pushover)
SDOF_ef_DISP_PUSH, SDOF_ef_MASS_PUSH, SDOF_ef_STIFF_PUSH, SDOF_ef_PERIOD_PUSH = [], [], [], []

# MEDIAN containers: one value per DUCT value 
PERIOD_RATIO_MED   = []
DI_RATIO_MED       = []
REACTION_RATIO_MED = []
DISP_RATIO_MED     = []
EDVR_RATIO_MED     = []
DISP_MED, VELO_MED, ACC_MED = [], [], []
EDCI_MED  = []
DII_MED, OMEGA_MED, MU_MED, RR_MED = [], [], [], []
EDVR_MED  = []
PERIOD_DYN_MED = []

# Analysis Durations:
starttime = TI.process_time() # start the CPU clock

for JJ in range(0, 21):                # 21 points -> 20 intervals of the sweep
    # linear interpolation of the ductility ratio between DUCT_MIN and DUCT_MAX
    DUCT = DUCT_MIN + (DUCT_MAX - DUCT_MIN) * (JJ / 20)
    COL_DUCT.append(DUCT)              # store the current design variable value

    II = 0                             # ground-motion index, will be overwritten later

    # per-ductility reset of the raw (non-median) accumulators
    DISP, VELO, ACC = [], [], []
    EDCI = []
    DII, OMEGA, MU, RR = [], [], [], []
    EDVR = []

    # STEP 1 : PUSHOVER ANALYSIS
    ANAL_TYPE = 'PUSHOVER'
    DATA = MDOF(DUCT, MAT_TYPE, TOTAL_MASS, ANAL_TYPE, II)   # call nonlinear MDOF solver
    (reaction_PUSH, disp_PUSH, DI_PUSH,
     ele_force_PUSH, node_displacements_PUSH,
     PERIOD_MIN_PUSH, PERIOD_MAX_PUSH,
     SDOF_EFFE_DISP_PUSH, SDOF_EFFE_MASS_PUSH, SDOF_EFFE_STIFF_PUSH,
     SDOF_EFFE_PERIOD_PUSH) = DATA                           # unpack solver output

    # store the MDOF -> SDOF equivalent system parameters
    SDOF_ef_DISP_PUSH.append(SDOF_EFFE_DISP_PUSH)
    SDOF_ef_MASS_PUSH.append(SDOF_EFFE_MASS_PUSH)
    SDOF_ef_STIFF_PUSH.append(SDOF_EFFE_STIFF_PUSH)
    SDOF_ef_PERIOD_PUSH.append(SDOF_EFFE_PERIOD_PUSH)

    # store the maximum (final) period of the pushover curve
    PERIOD_PUSH.append(max(PERIOD_MAX_PUSH))

    # STEP 2 : CYCLIC_DISPLACEMENT ANALYSIS  (reference / benchmark hysteresis)
    ANAL_TYPE = 'CYCLIC_DISPLACEMENT'
    DATA = MDOF(DUCT, MAT_TYPE, TOTAL_MASS, ANAL_TYPE, II)
    (reaction_CP, disp_CP, DI_CP,
     ele_force_CP, node_displacements_CP,
     PERIOD_MIN_CP, PERIOD_MAX_CP) = DATA

    # compute the equivalent viscous damping ratio of the cyclic hysteresis loop
    XLABEL = "Displacement [m]"
    YLABEL = "Base Reaction [N]"
    TITLE  = "Equivalent viscous damping ratio - Cyclic Displaecment Hysteresis"
    method = 1
    zeta_CP = S055.EQULIVALENT_VISCOUS_DAMPING_RATIO_FUN(
                  disp_CP, reaction_CP, method, XLABEL, YLABEL, TITLE)
    print(f"Equivalent viscous damping ratio = {zeta_CP:.4f}")

    # STEP 3 : SEISMIC LOOP  -> 20 ground motions
    for II in range(40, 60):           # 20 ground-motion records

        ANAL_TYPE = 'SEISMIC'
        DATA = MDOF(DUCT, MAT_TYPE, TOTAL_MASS, ANAL_TYPE, II)
        (time_SEI, reaction_SEI, disp_SEI, velo_SEI, acc_SEI, DI_SEI,
         ele_force_SEI, node_displacements_SEI,
         stiffness_SEI, PERIOD_SEI, damping_ratio_SEI,
         PERIOD_MIN_SEI, PERIOD_MAX_SEI) = DATA

        PERIOD_DYN.append(max(PERIOD_MAX_SEI))

        # DUCTILITY DAMAGE INDEX
        SLOPE_NODE = 10
        PLOT = True
        DIx, Omega_0, mu, R_mu, R = S12.DUCTILITY_DAMAGE_INDEX_FUN(
                                        disp_PUSH, reaction_PUSH,
                                        SLOPE_NODE, disp_SEI, PLOT)
        DII.append(DIx); OMEGA.append(Omega_0); MU.append(mu); RR.append(R)

        # ENERGY DISSIPATION CAPACITY INDEX
        edci = S10.ENERGY_DISSIPATION_CAPACITY_INDEX(
                   disp_SEI, reaction_SEI, disp_CP, reaction_CP)
        EDCI.append(edci)

        # peak response quantities
        DISP.append(max(disp_SEI))
        VELO.append(max(velo_SEI))
        ACC.append(max(acc_SEI))

        # ratio quantities (seismic vs. pushover)
        PERIOD_RATIO.append(  max(PERIOD_MAX_SEI) / max(PERIOD_MAX_PUSH))
        DI_RATIO.append(      max(DI_SEI)         / max(DI_PUSH))
        REACTION_RATIO.append(max(reaction_SEI)   / max(reaction_PUSH))
        DISP_RATIO.append(    max(disp_SEI)       / max(disp_PUSH))

        # equivalent viscous damping of the seismic hysteresis
        TITLE = "Equivalent viscous damping ratio - Seismic Hysteresis"
        zeta_SEI = S055.EQULIVALENT_VISCOUS_DAMPING_RATIO_FUN(
                       disp_SEI, reaction_SEI, method, XLABEL, YLABEL, TITLE)
        print(f"Equivalent viscous damping ratio = {zeta_SEI:.4f}")

        EDVR.append(zeta_SEI)
        EDVR_RATIO.append(zeta_SEI / zeta_CP)   # seismic damping / cyclic damping

    # STEP 4 : REDUCE 20 GROUND MOTIONS -> 1 MEDIAN VALUE PER DUCTILITY
    DISP_MED.append(np.median(DISP)); VELO_MED.append(np.median(VELO)); ACC_MED.append(np.median(ACC))
    PERIOD_RATIO_MED.append(  np.median(PERIOD_RATIO))
    DI_RATIO_MED.append(      np.median(DI_RATIO))
    REACTION_RATIO_MED.append(np.median(REACTION_RATIO))
    DISP_RATIO_MED.append(    np.median(DISP_RATIO))
    EDCI_MED.append(          np.median(EDCI))
    EDVR_RATIO_MED.append(    np.median(EDVR_RATIO))
    EDVR_MED.append(          np.median(EDVR))
    DII_MED.append(np.median(DII)); OMEGA_MED.append(np.median(OMEGA))
    MU_MED.append(np.median(MU));  RR_MED.append(np.median(RR))
    PERIOD_DYN_MED.append(np.median(PERIOD_DYN))


# TIMER REPORT
totaltime = TI.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n')

#%% ---------------------------------

XDATA = COL_DUCT
YDATA = SDOF_ef_DISP_PUSH
XLABEL = 'EACH COL. DUCTILITY RATIO'
YLABEL = 'EFFECTIVE DISPLACEMENT [m]' 
TITLE = f'{YLABEL} and {XLABEL} DURING PERIOD ANALYSIS \n EQUIVALENT SDOF SYSTEM DERIVATION VIA DISPLACEMENT-BASED PUSHOVER ANALYSIS'
COLOR = 'purple'
S01.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = COL_DUCT
YDATA = SDOF_ef_MASS_PUSH
XLABEL = 'EACH COL. DUCTILITY RATIO'
YLABEL = 'EFFECTIVE MASS [kg]' 
TITLE = f'{YLABEL} and {XLABEL} DURING PERIOD ANALYSIS \n EQUIVALENT SDOF SYSTEM DERIVATION VIA DISPLACEMENT-BASED PUSHOVER ANALYSIS'
COLOR = 'purple'
S01.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 3)

XDATA = COL_DUCT
YDATA = SDOF_ef_STIFF_PUSH
XLABEL = 'EACH COL. DUCTILITY RATIO'
YLABEL = 'EFFECTIVE STIFFNESS [N/m]' 
TITLE = f'{YLABEL} and {XLABEL} DURING PERIOD ANALYSIS \n EQUIVALENT SDOF SYSTEM DERIVATION VIA DISPLACEMENT-BASED PUSHOVER ANALYSIS'
COLOR = 'purple'
S01.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 3)

XDATA = COL_DUCT
YDATA = SDOF_ef_PERIOD_PUSH
XLABEL = 'EACH COL. DUCTILITY RATIO'
YLABEL = 'EFFECTIVE PERIOD [s]' 
TITLE = f'{YLABEL} and {XLABEL} DURING PERIOD ANALYSIS \n EQUIVALENT SDOF SYSTEM DERIVATION VIA DISPLACEMENT-BASED PUSHOVER ANALYSIS'
COLOR = 'purple'
S01.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = COL_DUCT
YDATA = DII_MED
XLABEL = 'EACH COL. DUCTILITY RATIO'
YLABEL = 'STRUCTURAL DAMAGE INDEX [%]' 
TITLE = f'{YLABEL} and {XLABEL} DURING PERIOD ANALYSIS WITH 20 GROUND MOTIONS'
COLOR = 'purple'
S01.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 3)

XDATA = COL_DUCT
YDATA = OMEGA_MED
XLABEL = 'EACH COL. DUCTILITY RATIO'
YLABEL = 'OVER-STRENGTH FACTOR' 
TITLE = f'{YLABEL} and {XLABEL} DURING PERIOD ANALYSIS WITH 20 GROUND MOTIONS'
COLOR = 'blue'
S01.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 3)

XDATA = COL_DUCT
YDATA = MU_MED
XLABEL = 'EACH COL. DUCTILITY RATIO'
YLABEL = 'STRUCTURAL DUCTILITY RATIO' 
TITLE = f'{YLABEL} and {XLABEL} DURING PERIOD ANALYSIS WITH 20 GROUND MOTIONS'
COLOR = 'blue'
S01.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = COL_DUCT
YDATA = RR_MED
XLABEL = 'EACH COL. DUCTILITY RATIO'
YLABEL = 'STRUCTURAL BEHAVIOR COEFFICIENT' 
TITLE = f'{YLABEL} and {XLABEL} DURING PERIOD ANALYSIS WITH 20 GROUND MOTIONS'
COLOR = 'blue'
S01.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = COL_DUCT
YDATA = EDCI_MED
XLABEL = 'EACH COL. DUCTILITY RATIO'
YLABEL = 'ENERGY DISSIPATION CAPACITY INDEX [%]' 
TITLE = f'{YLABEL} and {XLABEL} DURING PERIOD ANALYSIS WITH 20 GROUND MOTIONS'
COLOR = 'red'
S01.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 3)

XDATA = COL_DUCT
YDATA = EDVR_MED
XLABEL = 'EACH COL. DUCTILITY RATIO'
YLABEL = 'EQUIVALENT VISCOUS DAMPING RATIO [%]' 
TITLE = f'{YLABEL} and {XLABEL} DURING PERIOD ANALYSIS WITH 20 GROUND MOTIONS'
COLOR = 'brown'
S01.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 7)


XDATA = COL_DUCT
YDATA = DISP_MED
XLABEL = 'EACH COL. DUCTILITY RATIO'
YLABEL = 'MAX. DISPLACEMENT [m]' 
TITLE = f'{YLABEL} and {XLABEL} DURING PERIOD ANALYSIS WITH 20 GROUND MOTIONS'
COLOR = 'blue'
S01.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 3)

XDATA = COL_DUCT
YDATA = VELO_MED
XLABEL = 'EACH COL. DUCTILITY RATIO'
YLABEL = 'MAX. VELOCITY [m/s]' 
TITLE = f'{YLABEL} and {XLABEL} DURING PERIOD ANALYSIS WITH 20 GROUND MOTIONS'
COLOR = 'blue'
S01.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 3)

XDATA = COL_DUCT
YDATA = ACC_MED
XLABEL = 'EACH COL. DUCTILITY RATIO'
YLABEL = 'MAX. ACCELERATION [m/s^2]' 
TITLE = f'{YLABEL} and {XLABEL} DURING PERIOD ANALYSIS WITH 20 GROUND MOTIONS'
COLOR = 'blue'
S01.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 3)

XDATA = COL_DUCT
YDATA = PERIOD_PUSH
XLABEL = 'EACH COL. DUCTILITY RATIO'
YLABEL = 'MAX. PERIOD DURING PUSHOVER ANAL. [s]' 
TITLE = f'{YLABEL} and {XLABEL} DURING PERIOD ANALYSIS WITH 20 GROUND MOTIONS'
COLOR = 'green'
S01.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = COL_DUCT
YDATA = PERIOD_DYN_MED
XLABEL = 'EACH COL. DUCTILITY RATIO'
YLABEL = 'MAX. PERIOD DURING DYN. ANAL. [s]' 
TITLE = f'{YLABEL} and {XLABEL} DURING PERIOD ANALYSIS WITH 20 GROUND MOTIONS'
COLOR = 'pink'
S01.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = COL_DUCT
YDATA = PERIOD_RATIO_MED
XLABEL = 'EACH COL. DUCTILITY RATIO'
YLABEL = 'Tmax(SEI) / Tmax(PUSH)' # 'PERIOD RATIO (SEISMIC DIVDED BY CYCLIC DISPLACEMENT)'
TITLE = f'{YLABEL} and {XLABEL} DURING PERIOD ANALYSIS WITH 20 GROUND MOTIONS'
COLOR = 'purple'
S01.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 7)

XDATA = COL_DUCT
YDATA = REACTION_RATIO_MED
XLABEL = 'EACH COL. DUCTILITY RATIO'
YLABEL = 'REACTIONmax(SEI) / REACTIONmax(PUSH)' # 'REACTION RATIO (SEISMIC DIVDED BY CYCLIC DISPLACEMENT)'
TITLE = f'{YLABEL} and {XLABEL} DURING PERIOD ANALYSIS WITH 20 GROUND MOTIONS'
COLOR = 'red'
S01.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 3)

XDATA = COL_DUCT
YDATA = DISP_RATIO_MED
XLABEL = 'EACH COL. DUCTILITY RATIO'
YLABEL = 'DISPmax(SEI) / DISPmax(PUSH)' # 'DISPLACEMENT RATIO (SEISMIC DIVDED BY CYCLIC DISPLACEMENT)'
TITLE = f'{YLABEL} and {XLABEL} DURING PERIOD ANALYSIS WITH 20 GROUND MOTIONS'
COLOR = 'blue'
S01.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 3)

XDATA = COL_DUCT
YDATA = DI_RATIO_MED
XLABEL = 'EACH COL. DUCTILITY RATIO'
YLABEL = 'DImax(SEI) / DImax(PUSH)' # 'DUCTILITY DAMAGE INDEX (SEISMIC DIVDED BY CYCLIC DISPLACEMENT)'
TITLE = f'{YLABEL} and {XLABEL} DURING PERIOD ANALYSIS WITH 20 GROUND MOTIONS'
COLOR = 'green'
S01.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 3)

XDATA = COL_DUCT
YDATA = EDVR_RATIO_MED
XLABEL = 'EACH COL. DUCTILITY RATIO'
YLABEL = 'EDVR(SEI) / EDVR(PUSH)' # 'EDVR (SEISMIC DIVDED BY CYCLIC DISPLACEMENT)'
TITLE = f'{YLABEL} and {XLABEL} DURING PERIOD ANALYSIS WITH 20 GROUND MOTIONS'
COLOR = 'brown'
S01.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 3)


# 3D PLOT
X, Y, Z = DISP_RATIO_MED, DI_RATIO_MED, COL_DUCT
XLABEL, YLABEL, ZLABEL = 'DISPmax(SEI) / DISPmax(PUSH)', 'DImax(SEI) / DImax(PUSH)', 'EACH COL. DUCTILITY RATIO'               
S11.PLOT_CONTOUR_3D_2D_FUN(120, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = EDVR_RATIO_MED, DI_RATIO_MED, COL_DUCT
XLABEL, YLABEL, ZLABEL = 'EDVR(SEI) / EDVR(CP)', 'DImax(SEI) / DImax(CP)', 'EACH COL. DUCTILITY RATIO'               
S11.PLOT_CONTOUR_3D_2D_FUN(121, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = DISP_MED, DI_RATIO_MED, COL_DUCT
XLABEL, YLABEL, ZLABEL = 'MAX. DISPLACEMENT [m]', 'DImax(SEI) / DImax(PUSH)', 'EACH COL. DUCTILITY RATIO'               
S11.PLOT_CONTOUR_3D_2D_FUN(120, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = VELO_MED, DI_RATIO_MED, COL_DUCT
XLABEL, YLABEL, ZLABEL = 'MAX. VELOCITY  [m/s]', 'DImax(SEI) / DImax(PUSH)', 'EACH COL. DUCTILITY RATIO'               
S11.PLOT_CONTOUR_3D_2D_FUN(121, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = ACC_MED, DI_RATIO_MED, COL_DUCT
XLABEL, YLABEL, ZLABEL = 'MAX. ACCELERATION [m/s^2]', 'DImax(SEI) / DImax(PUSH)', 'EACH COL. DUCTILITY RATIO'               
S11.PLOT_CONTOUR_3D_2D_FUN(121, X, Y, Z, XLABEL, YLABEL, ZLABEL)

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
    "DISP":              DISP_MED,
    "VEL":               VELO_MED,
    "ACC":               ACC_MED,
    "REACTION_RATIO":    REACTION_RATIO_MED,
    "EDVR_RATIO":        EDVR_RATIO_MED,
    "DI_RATIO":          DI_RATIO_MED,
}


# Convert to DataFrame
df = pd.DataFrame(data)
#print(df)
S01.RANDOM_FOREST(df)
#%%------------------------------------------------------
# PLOT HEATMAP FOR CORRELATION 
S01.PLOT_HEATMAP(df)
#%%------------------------------------------------------
# MULTIPLE REGRESSION MODEL
#S01.MULTIPLE_REGRESSION(df) 
#%%-------------------------------------------------------------------
# Plots a heatmap of sensitivity coefficients (correlation or SRC) between inputs X and outputs Y.
import SENSITIVITY_HEATMAP_FUN as S099
X = np.column_stack([DISP_MED, VELO_MED, ACC_MED])          
Y = np.column_stack([EDVR_RATIO_MED, COL_DUCT, DI_RATIO_MED])  
X_LABELS = ['Max. Displacement [m]', 'Max. Velocity [m/s]', 'Max. Acceleration [m/s^2]']
Y_LABELS = ['EDVR RATIO', 'EACH COL. DUCTILITY RATIO', 'DI RATIO']
coeffs = S099.SENSITIVITY_HEATMAP_FUN(X, Y, X_LABELS, Y_LABELS, method='pearson')
#%%-------------------------------------------------------------------
# ANOVA with automatic binning of continuous predictors.
import ANOVA_SENSITIVITY_FUN as S100

df_sens = pd.DataFrame({
    "DISP":              DISP_MED,
    "VEL":               VELO_MED,
    "ACC":               ACC_MED,
    "REACTION_RATIO":    REACTION_RATIO_MED,
    "EDVR_RATIO":        EDVR_RATIO_MED,
    "DI_RATIO":          DI_RATIO_MED,
})

param_list = ["VEL", "ACC", "EDVR_RATIO", "DI_RATIO"]

# ANOVA – main effects only
anova_table, _ = S100.ANOVA_SENSITIVITY_FUN(
    df_sens, output_col="DISP",
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
    df_sens, output_col="DISP",
    param_cols=param_list,
    n_bins=4,
    include_interactions=True,
    plot=True,
)
plt.show()
print("\nANOVA table (with interactions):")
print(anova_table_inter.round(4))
#%%-------------------------------------------------------------------
exit()
#%%----------------------------------------------------
# PERIOD ANALYSIS
MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'
ANAL_TYPE = 'PERIOD'

DATA = MDOF(DUCT, MAT_TYPE, TOTAL_MASS, ANAL_TYPE, II)
(PERIOD_MIN_X, PERIOD_MAX_X) = DATA
print('Structure First Period:  ', PERIOD_MIN_X)
print('Structure Second Period: ', PERIOD_MAX_X) 


#%%----------------------------------------------------
# STATIC ANALYSIS
MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'
ANAL_TYPE = 'STATIC'

DATA = MDOF(DUCT, MAT_TYPE, TOTAL_MASS, ANAL_TYPE, II)
(reaction, disp, DI) = DATA

#%%----------------------------------------------------
# PUSHOVER ANALYSIS (STATIC TIME-HISTORY ANALYSIS)
MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'
ANAL_TYPE = 'PUSHOVER'

DATA = MDOF(DUCT, MAT_TYPE, TOTAL_MASS, ANAL_TYPE, II)
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
#%%----------------------------------------------------
# CYCLIC DISPLACEMENT ANALYSIS (STATIC TIME-HISTORY ANALYSIS)
MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'
ANAL_TYPE = 'CYCLIC_DISPLACEMENT'

DATA = MDOF(DUCT, MAT_TYPE, TOTAL_MASS, ANAL_TYPE, II)
(reaction_CP, disp_CP, DI_CP,
 ele_force_CP, node_displacements_CP,
 PERIOD_MIN_CP, PERIOD_MAX_CP) = DATA


XDATA = disp_CP
YDATA = reaction_CP
XLABEL = 'Displacement [m]'
YLABEL = 'Base Reaction [N]'
TITLE = 'Base Reaction and Displacement of Structure During Cyclic-Displacement Analysis'
COLOR = 'black'
SEMILOGY = False
PLOT(XDATA, YDATA, TITLE, XLABEL, YLABEL, COLOR, SEMILOGY)

# PLOT STRUCTURAL PERIOD DURING THE ANALYSIS
plt.figure(0, figsize=(12, 8))
plt.plot(disp_CP, PERIOD_MIN_CP, linewidth=3)
plt.plot(disp_CP, PERIOD_MAX_CP, linewidth=3)
plt.title('Period of Structure During Cyclic Displacement Analysis')
plt.ylabel('Structural Period [s]')
plt.xlabel('Displacement [m]')
#plt.semilogy()
plt.grid()
plt.legend([f'PERIOD - MIN VALUES: Min: {np.min(PERIOD_MIN_CP):.3f} (s) - Mean: {np.mean(PERIOD_MIN_CP):.3f} (s) - Max: {np.max(PERIOD_MIN_CP):.3f} (s)', 
            f'PERIOD - MAX VALUES:  Min: {np.min(PERIOD_MAX_CP):.3f} (s) - Mean: {np.mean(PERIOD_MAX_CP):.3f} (s) - Max: {np.max(PERIOD_MAX_CP):.3f} (s)',
            ])
plt.show()


# PLOT FORCES
YLABEL = 'Force [N]'
TITLE = "elements  force vs Time During Cyclic Displacement Analysis"
STOP = len(disp_CP)
NUM = STOP
step_CP = np.linspace(1, STOP, NUM)
PLOT_FORCES(step_CP, ele_force_CP, YLABEL, TITLE) # ELEMENTS FORCE

# PLOT NODAL DISPALEMENTS
STOP = len(disp_CP)
NUM = STOP
step_CP = np.linspace(1, STOP, NUM)
PLOT_DISPLAEMENTS(step_CP, node_displacements_CP, TITLE = "Node Displacements vs Time for During Cyclic Displacement Analysis") 


plt.figure(-1, figsize=(12, 8))
plt.plot(disp_CP, DI_CP, color='black', linewidth=2)
plt.xlabel('Displacement [m]')
plt.ylabel('Structural Damage Index [%]')
plt.title(f'Displacement vs Structural Damage Index')
plt.grid()
plt.show()

XLABEL = "Displacement [m]"
YLABEL = "Base Reaction [N]"
TITLE = "Equivalent viscous damping ratio - Cyclic Displacement Hysteresis"
method = 1 # 
zeta = S055.EQULIVALENT_VISCOUS_DAMPING_RATIO_FUN(disp_CP, reaction_CP, method, XLABEL, YLABEL, TITLE)
print(f"Equivalent viscous damping ratio = {zeta:.4f}")
#%%----------------------------------------------------
# EXTERNAL TIME-DEPENDENT LOADING ANALYSIS (STATIC TIME-HISTORY ANALYSIS)
MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'
ANAL_TYPE = 'STATIC_EXTERNAL_TIME-DEPENDENT_LOADING'

DATA = MDOF(DUCT, MAT_TYPE, TOTAL_MASS, ANAL_TYPE, II)

(reaction_ETDLS, disp_ETDLS, DI_ETDLS,
 ele_force_ETDLS, node_displacements_ETDLS,
 PERIOD_MIN_ETDLS, PERIOD_MAX_ETDLS) = DATA


XDATA = disp_ETDLS
YDATA = reaction_ETDLS
XLABEL = 'Displacement [m]'
YLABEL = 'Base Reaction [N]'
TITLE = 'Base Reaction and Displacement of Structure During Static External Time-dependent Loading Analysis'
COLOR = 'black'
SEMILOGY = False
PLOT(XDATA, YDATA, TITLE, XLABEL, YLABEL, COLOR, SEMILOGY)

# PLOT STRUCTURAL PERIOD DURING THE ANALYSIS
plt.figure(0, figsize=(12, 8))
plt.plot(disp_ETDLS, PERIOD_MIN_ETDLS, linewidth=3)
plt.plot(disp_ETDLS, PERIOD_MAX_ETDLS, linewidth=3)
plt.title('Period of Structure During Static External Time-dependent Loading Analysis')
plt.ylabel('Structural Period [s]')
plt.xlabel('Displacement [m]')
#plt.semilogy()
plt.grid()
plt.legend([f'PERIOD - MIN VALUES: Min: {np.min(PERIOD_MIN_ETDLS):.3f} (s) - Mean: {np.mean(PERIOD_MIN_ETDLS):.3f} (s) - Max: {np.max(PERIOD_MIN_ETDLS):.3f} (s)', 
            f'PERIOD - MAX VALUES:  Min: {np.min(PERIOD_MAX_ETDLS):.3f} (s) - Mean: {np.mean(PERIOD_MAX_ETDLS):.3f} (s) - Max: {np.max(PERIOD_MAX_ETDLS):.3f} (s)',
            ])
plt.show()


# PLOT FORCES
YLABEL = 'Force [N]'
TITLE = "elements  force vs Time During Static External Time-dependent Loading Analysis"
STOP = len(disp_ETDLS)
NUM = STOP
step_ETDLS = np.linspace(1, STOP, NUM)
PLOT_FORCES(step_ETDLS, ele_force_ETDLS, YLABEL, TITLE) # ELEMENTS FORCE

# PLOT NODAL DISPALEMENTS
STOP = len(disp_ETDLS)
NUM = STOP
step_ETDLS = np.linspace(1, STOP, NUM)
PLOT_DISPLAEMENTS(step_ETDLS, node_displacements_ETDLS, TITLE = "Node Displacements vs Time for During Static External Time-dependent Loading Analysis") 


plt.figure(-1, figsize=(12, 8))
plt.plot(disp_ETDLS, DI_ETDLS, color='black', linewidth=2)
plt.xlabel('Displacement [m]')
plt.ylabel('Structural Damage Index [%]')
plt.title(f'Displacement vs Structural Damage Index')
plt.grid()
plt.show()

XLABEL = "Displacement [m]"
YLABEL = "Base Reaction [N]"
TITLE = "Equivalent viscous damping ratio - Static External Time-dependent Loading Hysteresis"
method = 2 # 
zeta = S055.EQULIVALENT_VISCOUS_DAMPING_RATIO_FUN(disp_ETDLS, reaction_ETDLS, method, XLABEL, YLABEL, TITLE)
print(f"Equivalent viscous damping ratio = {zeta:.4f}")
#%%----------------------------------------------------
# EXTERNAL TIME-DEPENDENT LOADING ANALYSIS (DYNAMIC TIME-HISTORY ANALYSIS)
MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'
ANAL_TYPE = 'DYNAMIC_EXTERNAL_TIME-DEPENDENT_LOADING'

DATA = MDOF(DUCT, MAT_TYPE, TOTAL_MASS, ANAL_TYPE, II)

(time_ETDLD, reaction_ETDLD, disp_ETDLD, velo_ETDLD, acc_ETDLD,  DI_ETDLD,
 ele_force_ETDLD, node_displacements_ETDLD,
stiffness, PERIOD, damping_ratio,
PERIOD_MIN_ETDLD, PERIOD_MAX_ETDLD) = DATA



XDATA = disp_ETDLD
YDATA = reaction_ETDLD
XLABEL = 'Displacement [m]'
YLABEL = 'Base Reaction [N]'
TITLE = 'Base Reaction and Displacement of Structure During Dynamic External Time-dependent Loading Analysis'
COLOR = 'black'
SEMILOGY = False
PLOT(XDATA, YDATA, TITLE, XLABEL, YLABEL, COLOR, SEMILOGY)

PLOT_TIME_HISTORY(time_ETDLD, reaction_ETDLD, disp_ETDLD, velo_ETDLD, acc_ETDLD)

# PLOT STRUCTURAL PERIOD DURING THE ANALYSIS
plt.figure(0, figsize=(12, 8))
plt.plot(disp_ETDLD, PERIOD_MIN_ETDLD, linewidth=3)
plt.plot(disp_ETDLD, PERIOD_MAX_ETDLD, linewidth=3)
plt.title('Period of Structure During Dynamic External Time-dependent Loading Analysis')
plt.ylabel('Structural Period [s]')
plt.xlabel('Displacement [m]')
#plt.semilogy()
plt.grid()
plt.legend([f'PERIOD - MIN VALUES: Min: {np.min(PERIOD_MIN_ETDLD):.3f} (s) - Mean: {np.mean(PERIOD_MIN_ETDLD):.3f} (s) - Max: {np.max(PERIOD_MIN_ETDLD):.3f} (s)', 
            f'PERIOD - MAX VALUES:  Min: {np.min(PERIOD_MAX_ETDLD):.3f} (s) - Mean: {np.mean(PERIOD_MAX_ETDLD):.3f} (s) - Max: {np.max(PERIOD_MAX_ETDLD):.3f} (s)',
            ])
plt.show()

# PLOT FORCES
YLABEL = 'Force [N]'
TITLE = "elements  force vs Time During Dynamic External Time-dependent Loading Analysis"
PLOT_FORCES(time_ETDLD, ele_force_ETDLD, YLABEL, TITLE) # ELEMENTS FORCE

# PLOT NODAL DISPALEMENTS
PLOT_DISPLAEMENTS(time_ETDLD, node_displacements_ETDLD, TITLE = "Node Displacements vs Time for During Dynamic External Time-dependent Loading Analysis") 

plt.figure(-1, figsize=(12, 8))
plt.plot(disp_ETDLD, DI_ETDLS, color='black', linewidth=2)
plt.xlabel('Displacement [m]')
plt.ylabel('Structural Damage Index [%]')
plt.title(f'Displacement vs Structural Damage Index')
plt.grid()
plt.show()

XLABEL = "Displacement [m]"
YLABEL = "Base Reaction [N]"
TITLE = "Equivalent viscous damping ratio - Dynamic External Time-dependent Loading Hysteresis"
method = 2 # 
zeta = S055.EQULIVALENT_VISCOUS_DAMPING_RATIO_FUN(disp_ETDLD, reaction_ETDLD, method, XLABEL, YLABEL, TITLE)
print(f"Equivalent viscous damping ratio = {zeta:.4f}")
#%%----------------------------------------------------
# FREE-VIBRATION ANALYSIS (DYNAMIC TIME-HISTORY ANALYSIS)
MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'
ANAL_TYPE = 'FREE-VIBRATION'
DATA = MDOF(DUCT, MAT_TYPE, TOTAL_MASS, ANAL_TYPE, II)

(time_FV, reaction_FV, disp_FV, velo_FV, acc_FV, DI_FV,
 ele_force_FV, node_displacements_FV,
stiffness_FV, PERIOD_FV, damping_ratio_FV,
PERIOD_MIN_FV, PERIOD_MAX_FV) = DATA

XDATA = disp_FV
YDATA = reaction_FV
XLABEL = 'Displacement [m]'
YLABEL = 'Base Reaction [N]'
TITLE = 'Base Reaction and Dispalcement of Structure During Free-vibration Analysis'
COLOR = 'black'
SEMILOGY = False
PLOT(XDATA, YDATA, TITLE, XLABEL, YLABEL, COLOR, SEMILOGY)

PLOT_TIME_HISTORY(time_FV, reaction_FV, disp_FV, velo_FV, acc_FV)

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

# PLOT FORCES
YLABEL = 'Force [N]'
TITLE = "elements  force vs Time During Free-vibration Analysis"
PLOT_FORCES(time_FV, ele_force_FV, YLABEL, TITLE) # ELEMENTS FORCE

# PLOT NODAL DISPALEMENTS
PLOT_DISPLAEMENTS(time_FV, node_displacements_FV, TITLE = "Node Displacements vs Time for During Free-vibration Analysis")

plt.figure(-1, figsize=(12, 8))
plt.plot(disp_FV, DI_FV, color='black', linewidth=2)
plt.xlabel('Displacement [m]')
plt.ylabel('Structural Damage Index [%]')
plt.title(f'Displacement vs Structural Damage Index')
plt.grid()
plt.show()

X, HISTO_COLOR, LABEL = disp_FV, 'cyan', 'DISPLACEMENT FROM FREE-VIBRATION ANALYSIS'
S09.HISROGRAM_BOXPLOT(X, HISTO_COLOR, LABEL)
X, HISTO_COLOR, LABEL = reaction_FV, 'purple', 'REACTION FROM FREE-VIBRATION ANALYSIS'
S09.HISROGRAM_BOXPLOT(X, HISTO_COLOR, LABEL)
X, HISTO_COLOR, LABEL = DI_FV, 'orange', 'DUCTILITY DAMAGE INDEX FROM FREE-VIBRATION ANALYSIS'
S09.HISROGRAM_BOXPLOT(X, HISTO_COLOR, LABEL)

XLABEL = "Displacement [m]"
YLABEL = "Base Reaction [N]"
TITLE = "Equivalent viscous damping ratio - Free-vibration Hysteresis"
method = 2 # 
zeta = S055.EQULIVALENT_VISCOUS_DAMPING_RATIO_FUN(disp_FV, reaction_FV, method, XLABEL, YLABEL, TITLE)
print(f"Equivalent viscous damping ratio = {zeta:.4f}")
#%%----------------------------------------------------
# SEISMIC ANALYSIS (DYNAMIC TIME-HISTORY ANALYSIS)
MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'
ANAL_TYPE = 'SEISMIC'

DATA = MDOF(DUCT, MAT_TYPE, TOTAL_MASS, ANAL_TYPE, II)

(time_SEI, reaction_SEI, disp_SEI, velo_SEI, acc_SEI, DI_SEI,
 ele_force_SEI, node_displacements_SEI,
 stiffness_SEI, PERIOD_SEI, damping_ratio_SEI,
 PERIOD_MIN_SEI, PERIOD_MAX_SEI) = DATA


XDATA = disp_SEI
YDATA = reaction_SEI
XLABEL = 'Displacement [m]'
YLABEL = 'Base Reaction [N]'
TITLE = 'Base Reaction and Dispalcement of Structure During Seismic Analysis'
COLOR = 'black'
SEMILOGY = False
PLOT(XDATA, YDATA, TITLE, XLABEL, YLABEL, COLOR, SEMILOGY)

PLOT_TIME_HISTORY(time_SEI, reaction_SEI, disp_SEI, velo_SEI, acc_SEI)

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

# PLOT FORCES
YLABEL = 'Force [N]'
TITLE = "elements  force vs Time During Seismic Analysis"
PLOT_FORCES(time_SEI, ele_force_SEI, YLABEL, TITLE) # ELEMENTS FORCE

# PLOT NODAL DISPALEMENTS
PLOT_DISPLAEMENTS(time_SEI, node_displacements_SEI, TITLE = "Node Displacements vs Time for During Seismic Analysis")

plt.figure(-1, figsize=(12, 8))
plt.plot(disp_SEI, DI_SEI, color='black', linewidth=2)
plt.xlabel('Displacement [m]')
plt.ylabel('Structural Damage Index [%]')
plt.title(f'Displacement vs Structural Damage Index')
plt.grid()
plt.show()

X, HISTO_COLOR, LABEL = disp_SEI, 'cyan', 'DISPLACEMENT FROM SEISMIC ANALYSIS'
S09.HISROGRAM_BOXPLOT(X, HISTO_COLOR, LABEL)
X, HISTO_COLOR, LABEL = reaction_SEI, 'purple', 'REACTION FROM SEISMIC ANALYSIS'
S09.HISROGRAM_BOXPLOT(X, HISTO_COLOR, LABEL)
X, HISTO_COLOR, LABEL = DI_SEI, 'orange', 'DUCTILITY DAMAGE INDEX FROM SEISMIC ANALYSIS'
S09.HISROGRAM_BOXPLOT(X, HISTO_COLOR, LABEL)

XLABEL = "Displacement [m]"
YLABEL = "Base Reaction [N]"
TITLE = "Equivalent viscous damping ratio - Seismic Hysteresis"
method = 1 # 
zeta = S055.EQULIVALENT_VISCOUS_DAMPING_RATIO_FUN(disp_SEI, reaction_SEI, method, XLABEL, YLABEL, TITLE)
print(f"Equivalent viscous damping ratio = {zeta:.4f}")
#%%----------------------------------------------------
# --------------------------------------
#  Plot BaseAxial-Displacement Analysis 
# --------------------------------------
XX = np.abs(disp_PUSH); YY = np.abs(reaction_PUSH); # ABSOLUTE VALUE
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
PLOT_2D(np.abs(disp_PUSH), np.abs(reaction_PUSH), X, Y, X, Y, XLABEL, YLABEL, TITLE, LEGEND01, LEGEND02, LEGEND03, COLOR='black', Z=2) 
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
Dd = np.max(np.abs(disp_SEI))
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
    disp_SEI, reaction_SEI, method = 2, 
    title="Earthquake Response – Dissipated Energy"
)
fig_SEI.show()

print(f"Dissipated Energy from Earthquake= {Ed_SEI:.2f} N·m")

Ed_CP, fig_CP = DISSIPATED_ENERGY_FUN_WITH_PLOT(
    disp_CP, reaction_CP, method = 2,
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
import FRAGILITY_CURVE_FUN as FF

DAMAGE_STATES = {
    'Minor Damage Level': 20.0, # Median Value = 20%
    'Moderate Damage Level': 40.0,
    'Severe Damage Level': 70.0,
    'Failure Level': 100.0
}


A = 50.0      # Scale factor (median EDP at IM = 1)
B = 1.2       # Power‑law exponent (elastic or inelastic response trend)
BETA = 0.4    # Dispersion (logarithmic standard deviation of EDP)

IMs, XLABEL, SCATTER, SEMILOGY = DI_ETDLD, 'DYN. TIME-HISTORY ANA. DUCTILITY DAMAGE INDEX [%]', False, False
FF.FRAGILITY_CURVE_FUN_02(IMs, DAMAGE_STATES, XLABEL, SCATTER, SEMILOGY, A, B, BETA)

IMs, XLABEL, SCATTER, SEMILOGY = DI_FV, 'FREE-VIBRTION ANA. DUCTILITY DAMAGE INDEX [%]', False, False
FF.FRAGILITY_CURVE_FUN_02(IMs, DAMAGE_STATES, XLABEL, SCATTER, SEMILOGY, A, B, BETA)
    
IMs, XLABEL, SCATTER, SEMILOGY = DI_SEI, 'SEISMIC ANA. DUCTILITY DAMAGE INDEX [%]', False, False
FF.FRAGILITY_CURVE_FUN_02(IMs, DAMAGE_STATES, XLABEL, SCATTER, SEMILOGY, A, B, BETA)
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
with pd.ExcelWriter("MDOF_10_FLOORS_04_COLUMNS_8_ANA_OUTPUT.xlsx", engine='openpyxl') as writer:
    
    # PUSHOVER
    df1 = create_df(reaction_PUSH, disp_PUSH, PERIOD_MIN_PUSH, PERIOD_MAX_PUSH)
    df1.to_excel(writer, sheet_name="PUSHOVER", index=False)
                 
    # CYCLIC DISPLACEMENT
    df1 = create_df(reaction_CP, disp_CP, PERIOD_MIN_CP, PERIOD_MAX_CP)
    df1.to_excel(writer, sheet_name="CYCLIC_DISPLACEMENT", index=False)
    
    # STATIC EXTERNAL TIME-DEPENDENT LOADING
    df2 = create_df(reaction_ETDLS, disp_ETDLS, PERIOD_MIN_ETDLS, PERIOD_MAX_ETDLS)
    df2.to_excel(writer, sheet_name="STATIC_EXTERNAL_TIME-DEPENDENT_LOADING", index=False)

    # DYNAMIC EXTERNAL TIME-DEPENDENT LOADING
    df3 = create_df(reaction_ETDLD, disp_ETDLD, PERIOD_MIN_ETDLD, PERIOD_MAX_ETDLD)
    df3.to_excel(writer, sheet_name="DYNAMIC_EXTERNAL_TIME-DEPENDENT_LOADING", index=False)
    
    # FREE-VIBRATION
    df3 = create_df(reaction_FV, disp_FV, PERIOD_MIN_FV, PERIOD_MAX_FV)
    df3.to_excel(writer, sheet_name="FREE-VIBRATION", index=False)

    # SEISMIC
    df4 = create_df(reaction_SEI, disp_SEI, PERIOD_MIN_SEI, PERIOD_MAX_SEI)
    df4.to_excel(writer, sheet_name="SEISMIC", index=False)
#%%----------------------------------------------------