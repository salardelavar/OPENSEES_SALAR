###########################################################################################################
#                   >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                      #
#          SENSITIVITY ANALYSIS OF COLUMN YIELD STRENGTH, ELASTIC STIFFNESS, DUCTILITY RATIO,             #
#    OVER-STRENGTH FACTOR AND OPENSEES VIA PUSHOVER ANALYSIS OF A MULTI-DEGREE-OF-FREEDOM STRUCTURE       #
#           AND EVALUATION OF A MULTILINEAR FITTING CURVE FOR STRUCTURAL ELASTIC STIFFNESS,               #
#                    PLASTIC STIFFNESS, DUCTILITY RATIO, AND OVER-STRENGTH FACTOR                         #
#---------------------------------------------------------------------------------------------------------#
# EQUIVALENT SDOF SYSTEM DERIVATION VIA DISPLACEMENT-BASED SEISMIC DESIGN PROCEDURE WITH PUSHOVER ANALYSIS#
#---------------------------------------------------------------------------------------------------------#
#                  THIS PYTHON SCRIPT IS WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                    #
#                                   EMAIL: salar.d.ghashghaei@gmail.com                                   #
###########################################################################################################
"""
# ==========================================================================================
#   SENSITIVITY ANALYSIS OF A MULTI-DEGREE-OF-FREEDOM (MDOF) BUILDING FRAME
#   ------------------------------------------------------------------------------
#   Effect of Column Yield Strength, Elastic Stiffness, Ductility Ratio, Over-stength Fctor, on the Equivalent SDOF System and on the
#   Multilinear Idealization of the Pushover Curve
#
#   METHODOLOGY
#   -----------
#   1.  A 10-story MDOF shear-building frame is modelled in OpenSeesPy
#       (10 lumped masses, 4 parallel zeroLength springs per floor).
#
#   2.  The COLUMN DUCTILITY RATIO is chosen as the single design parameter and is
#       swept over 21 equally spaced values between DUCT_MIN and DUCT_MAX.
#
#   3.  For every value of the swept parameter, a nonlinear PUSHOVER ANALYSIS is
#       carried out under displacement control applied at the top node.
#
#   4.  The resulting base-shear / displacement curve is post-processed to derive the
#       EQUIVALENT SDOF SYSTEM (effective displacement, mass, stiffness, and period)
#       for displacement-based seismic design.
#
#   5.  The same pushover curve is then idealised with four multilinear fits —
#       BILINEAR, TRILINEAR, QUADRILINEAR, and PENTALINEAR — using an area-preserving
#       piecewise-linear algorithm.
#
#   6.  The following structural metrics are extracted from each fit:
#           • Elastic Stiffness
#           • Plastic Stiffness
#           • Ductility Ratio
#           • Over-Strength Factor
#
#   7.  The whole set of metrics — SDOF-equivalent and multilinear — is then plotted
#       against the swept column ductility ratio to expose trends, sensitivities,
#       and the influence of the idealization order (bi → penta).
#
#   PURPOSE
#   -------
#   Quantify how the column ductility ratio drives (a) the equivalent SDOF
#   parameters used in displacement-based design and (b) the parameters obtained
#   from multilinear fitting of the pushover curve — and to what extent the choice
#   of idealization order (bi/tri/quad/penta-linear) changes those conclusions.
# ==========================================================================================

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
import EQULIVALENT_VISCOUS_DAMPING_RATIO_FUN as S055
import DAMAGE_INDEX_FUN as S066
import PLOT_CONTOUR_3D_2D_FUN as S11
#%%----------------------------------------------------
def MDOF(FYz, KEz, DUCTz, OSFz, MAT_TYPE, TOTAL_MASS, ANAL_TYPE):
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
    FYi = [FYz,      # COLUMN 01
           FYz,      # COLUMN 02
           FYz,      # COLUMN 03
           FYz]      # COLUMN 04
    
    # ULTIMATE STRENGTH [N]    
    FUi = [OSFz * FYi[0],   # COLUMN 01
           OSFz * FYi[1],   # COLUMN 02
           OSFz * FYi[2],   # COLUMN 03
           OSFz * FYi[3]]   # COLUMN 04
    
    # ELASTIC STIFFNESS [N/m]
    Kei = [KEz,       # COLUMN 01
           KEz,       # COLUMN 02
           KEz,       # COLUMN 03
           KEz]       # COLUMN 04
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
            DSU = DUCTz * DY                                 # [m] Ultimate Displacement
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
        DMAX = -1.0*DSU     # [m] Max. Displacement
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

        ax1.plot(abs(displacement_X_10), EFFECTIVE_DISP_X, color='green', linewidth=3)
        ax1.set_title(f'Effective Displacement - Median: {np.median(EFFECTIVE_DISP_X): .5f}')
        ax1.set_xlabel('Abs. Top Displacement [m]')
        ax1.set_ylabel('Effective Displacement [m]')
        #ax1.legend(loc='upper right')
        ax1.grid(True)

        ax2.plot(abs(displacement_X_10), EFFECTIVE_MASS_X, color='magenta', linewidth=3)
        ax2.set_title(f'Effective Mass - Median: {np.median(EFFECTIVE_MASS_X): .5f}')
        ax2.set_xlabel('Abs. Top Displacement [m]')
        ax2.set_ylabel('Effective Mass [kg]')
        #ax2.legend(loc='upper right')
        ax2.grid(True)

        ax3.plot(abs(displacement_X_10), EFFECTIVE_STIFF_X, color='cyan', linewidth=3)
        ax3.set_title(f'Effective Stiffness - Median: {np.median(EFFECTIVE_STIFF_X): .5f}')
        ax3.set_xlabel('Step')
        ax3.set_ylabel('Effective Stiffness')
        #ax3.legend(loc='upper right')
        ax3.grid(True)

        ax4.plot(abs(displacement_X_10), EFFECTIVE_PERIOD_X, color='black', linewidth=3)
        ax4.set_title(f'Effective Period - Median: {np.median(EFFECTIVE_PERIOD_X): .5f}')
        ax4.set_xlabel('Abs. Top Displacement [m]')
        ax4.set_ylabel('Effective Period [s]')
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
#%%----------------------------------------------------
TOTAL_MASS = 5_000_000.0   # [kg] Total Mass of Structure
#%%----------------------------------------------------
MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'

# --------------------------------------------------------------
# SENSITIVITY ANALYSIS BY CHANGING EACH COLUMN DUCTILITY RATIO
# --------------------------------------------------------------

import time as TI

FY_MIN = 0.80 * 85000.0     # MIN. COLUMN'S YIELD STRENGTH 
FY_MAX = 1.20 * 85000.0     # MAX. COLUMN'S YIELD STRENGTH 

KE_MIN = 0.80 * 4500000.0   # MIN. COLUMN'S ELASTIC STIFFNESS  
KE_MAX = 1.20 * 4500000.0   # MAX. COLUMN'S ELASTIC STIFFNESS  

DUCT_MIN = 20.0             # MIN. COLUMN'S DUCTILITY RATIO
DUCT_MAX = 50.0             # MAX. COLUMN'S DUCTILITY RATIO

OSF_MIN = 1.05              # MIN. COLUMN'S OVER-STRENGTH FACTOR
OSF_MAX = 1.20              # MAX. COLUMN'S OVER-STRENGTH FACTOR

COL_DUCT = []   
COL_FY = []
COL_KE = []
COL_OSF = []
PERIOD_PUSH = []   # max period from the pushover analysis

# SDOF-equivalent quantities derived from the pushover curve (displacement-based pushover)
SDOF_ef_DISP_PUSH, SDOF_ef_MASS_PUSH, SDOF_ef_STIFF_PUSH, SDOF_ef_PERIOD_PUSH = [], [], [], []

C_Elastic_ST_bi, C_Plastic_ST_bi, C_Ductility_Rito_bi, C_Over_Strength_Factor_bi = [], [], [], []
C_Elastic_ST_tri, C_Plastic_ST_tri, C_Ductility_Rito_tri, C_Over_Strength_Factor_tri = [], [], [], []
C_Elastic_ST_quad, C_Plastic_ST_quad, C_Ductility_Rito_quad, C_Over_Strength_Factor_quad = [], [], [], []
C_Elastic_ST_penta, C_Plastic_ST_penta, C_Ductility_Rito_penta, C_Over_Strength_Factor_penta = [], [], [], []

STEP = 0
# Analysis Durations:
starttime = TI.process_time() # start the CPU clock

for XX in range(0, 6):                               # COLUMN YIELD STRENGTH [N]
    # linear interpolation of the ductility ratio between FY_MIN and FY_MAX
    FYz = FY_MIN + (FY_MAX - FY_MIN) * (XX / 5)
    
    for YY in range(0, 6):                           # COLUMN ELASTIC STIFFNESS [N/m]
        # linear interpolation of the ductility ratio between KE_MIN and KE_MAX
        KEz = KE_MIN + (KE_MAX - KE_MIN) * (YY / 5)
        
        for ZZ in range(0, 6):                       # COLUMN OVER-STRENGTH FACTOR [N/N]
            # linear interpolation of the ductility ratio between OSF_MIN and OSF_MAX
            OSFz = OSF_MIN + (OSF_MAX - OSF_MIN) * (ZZ / 5)   
            
            for JJ in range(0, 6):                   # COLUMN DUCTILITY RATIO [m/m]
                # linear interpolation of the ductility ratio between DUCT_MIN and DUCT_MAX
                DUCTz = DUCT_MIN + (DUCT_MAX - DUCT_MIN) * (JJ / 5)
                
                COL_DUCT.append(DUCTz)              # store the current variable value
                COL_FY.append(FYz)                  # store the current variable value
                COL_KE.append(KEz)                  # store the current variable value
                COL_OSF.append(OSFz)                # store the current variable value
            
                # PUSHOVER ANALYSIS
                ANAL_TYPE = 'PUSHOVER'
                DATA = MDOF(FYz, KEz, DUCTz, OSFz, MAT_TYPE, TOTAL_MASS, ANAL_TYPE)   # call nonlinear MDOF solver
                
                (reaction_PUSH, disp_PUSH, DI_PUSH,
                 ele_force_PUSH, node_displacements_PUSH,
                 PERIOD_MIN_PUSH, PERIOD_MAX_PUSH,
                 SDOF_EFFE_DISP_PUSH, SDOF_EFFE_MASS_PUSH, SDOF_EFFE_STIFF_PUSH,
                 SDOF_EFFE_PERIOD_PUSH) = DATA                           # unpack solver output
            
                # PLOT THE SPRING
                #S01.PLOT_1D_SPRING(deformed_scale=1.0, virtual_spring_length=1.0, virtual_spring_angle=0.0)
                
                # store the MDOF -> SDOF equivalent system parameters
                SDOF_ef_DISP_PUSH.append(SDOF_EFFE_DISP_PUSH)
                SDOF_ef_MASS_PUSH.append(SDOF_EFFE_MASS_PUSH)
                SDOF_ef_STIFF_PUSH.append(SDOF_EFFE_STIFF_PUSH)
                SDOF_ef_PERIOD_PUSH.append(SDOF_EFFE_PERIOD_PUSH)
            
                # store the maximum (final) period of the pushover curve
                PERIOD_PUSH.append(max(PERIOD_MAX_PUSH))
                
                SLOPE_NODE = 1
            
                # Choose breakpoint indices from the original curve
                # %%
                TRI_NODE = int(0.6 * len(disp_PUSH))
            
                QUAD_NODES = (int(0.6 * (len(disp_PUSH))),
                              int(0.65 * (len(disp_PUSH))))
            
                PENTA_NODES = (int(0.6 * (len(disp_PUSH))),
                               int(0.65 * (len(disp_PUSH))),
                               int(0.9 * (len(disp_PUSH))))
            
                #%% Bilinear Curve Fitted
                disp_PUSH, reaction_PUSH = np.abs(disp_PUSH), np.abs(reaction_PUSH)
                zdata = BILINEAR_CURVE(disp_PUSH, reaction_PUSH, SLOPE_NODE)
                (X_bi, Y_bi,
                 Elastic_ST_bi, Plastic_ST_bi,
                 Tangent_ST_bi, Ductility_Rito_bi,
                 Over_Strength_Factor_bi) = zdata
                
                C_Elastic_ST_bi.append(Elastic_ST_bi)
                C_Plastic_ST_bi.append(Plastic_ST_bi)
                C_Ductility_Rito_bi.append(Ductility_Rito_bi)
                C_Over_Strength_Factor_bi.append(Over_Strength_Factor_bi)
                #%% Trilinear Curve Fitted
                zdata = TRILINEAR_CURVE(disp_PUSH, reaction_PUSH, SLOPE_NODE, TRI_NODE)
                (X_tri, Y_tri,
                 Elastic_ST_tri, Plastic_ST_tri,
                 Tangent_ST_tri, Ductility_Rito_tri,
                 Over_Strength_Factor_tri) = zdata
                
                C_Elastic_ST_tri.append(Elastic_ST_tri)
                C_Plastic_ST_tri.append(Plastic_ST_tri)
                C_Ductility_Rito_tri.append(Ductility_Rito_tri)
                C_Over_Strength_Factor_tri.append(Over_Strength_Factor_tri)
                #%% Quadrilinear Curve Fitted
                zdata = QUADRILINEAR_CURVE(disp_PUSH, reaction_PUSH, SLOPE_NODE, *QUAD_NODES)
                (X_quad, Y_quad,
                 Elastic_ST_quad, Plastic_ST_quad,
                 Tangent_ST_quad, Ductility_Rito_quad,
                 Over_Strength_Factor_quad) = zdata
                
                C_Elastic_ST_quad.append(Elastic_ST_quad)
                C_Plastic_ST_quad.append(Plastic_ST_quad)
                C_Ductility_Rito_quad.append(Ductility_Rito_quad)
                C_Over_Strength_Factor_quad.append(Over_Strength_Factor_quad)
                #%% Pentalinear Curve Fitted
                zdata = PENTALINEAR_CURVE(disp_PUSH, reaction_PUSH, SLOPE_NODE, *PENTA_NODES)
                (X_penta, Y_penta,
                 Elastic_ST_penta, Plastic_ST_penta,
                 Tangent_ST_penta, Ductility_Rito_penta,
                 Over_Strength_Factor_penta) = zdata
            
                C_Elastic_ST_penta.append(Elastic_ST_penta)
                C_Plastic_ST_penta.append(Plastic_ST_penta)
                C_Ductility_Rito_penta.append(Ductility_Rito_penta)
                C_Over_Strength_Factor_penta.append(Over_Strength_Factor_penta)
                """
                #%% Plot
                # --------------------------------------
                #  Plot BaseAxial-Displacement Analysis 
                # --------------------------------------
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
                """
                STEP = STEP + 1
                #if STEP == 100:
                #    break
                print(f'\t\t\t STEP {STEP} DONE. \n\n')


#exit()
# TIMER REPORT
totaltime = TI.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n')
#%%----------------------------------------------------
# --------------------------------------
#  Plot BaseAxial-Displacement Analysis 
# --------------------------------------
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
                
XDATA = COL_DUCT
YDATA = SDOF_ef_DISP_PUSH
XLABEL = 'EACH COL. DUCTILITY RATIO'
YLABEL = 'EFFECTIVE DISPLACEMENT [m]' 
TITLE = f'{YLABEL} and {XLABEL} DURING PERIOD ANALYSIS \n EQUIVALENT SDOF SYSTEM DERIVATION VIA DISPLACEMENT-BASED PUSHOVER ANALYSIS'
COLOR = 'purple'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = COL_DUCT
YDATA = SDOF_ef_MASS_PUSH
XLABEL = 'EACH COL. DUCTILITY RATIO'
YLABEL = 'EFFECTIVE MASS [kg]' 
TITLE = f'{YLABEL} and {XLABEL} DURING PERIOD ANALYSIS \n EQUIVALENT SDOF SYSTEM DERIVATION VIA DISPLACEMENT-BASED PUSHOVER ANALYSIS'
COLOR = 'purple'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = COL_DUCT
YDATA = SDOF_ef_STIFF_PUSH
XLABEL = 'EACH COL. DUCTILITY RATIO'
YLABEL = 'EFFECTIVE STIFFNESS [N/m]' 
TITLE = f'{YLABEL} and {XLABEL} DURING PERIOD ANALYSIS \n EQUIVALENT SDOF SYSTEM DERIVATION VIA DISPLACEMENT-BASED PUSHOVER ANALYSIS'
COLOR = 'purple'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = COL_DUCT
YDATA = SDOF_ef_PERIOD_PUSH
XLABEL = 'EACH COL. DUCTILITY RATIO'
YLABEL = 'EFFECTIVE PERIOD [s]' 
TITLE = f'{YLABEL} and {XLABEL} DURING PERIOD ANALYSIS \n EQUIVALENT SDOF SYSTEM DERIVATION VIA DISPLACEMENT-BASED PUSHOVER ANALYSIS'
COLOR = 'purple'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)


XDATA = COL_FY
YDATA = PERIOD_PUSH
XLABEL = 'Column Yield Strength [N]'
YLABEL = 'MAX. STRUCTURAL PERIOD [s]' 
TITLE = f'{YLABEL} and {XLABEL} DURING PERIOD ANALYSIS'
COLOR = 'brown'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = COL_KE
YDATA = PERIOD_PUSH
XLABEL = 'Column Elastic Stiffness [N/m]'
YLABEL = 'MAX. STRUCTURAL PERIOD [s]' 
TITLE = f'{YLABEL} and {XLABEL} DURING PERIOD ANALYSIS'
COLOR = 'brown'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = COL_DUCT
YDATA = PERIOD_PUSH
XLABEL = 'EACH COL. DUCTILITY RATIO [m/m]'
YLABEL = 'MAX. STRUCTURAL PERIOD [s]' 
TITLE = f'{YLABEL} and {XLABEL} DURING PERIOD ANALYSIS'
COLOR = 'brown'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = COL_OSF
YDATA = PERIOD_PUSH
XLABEL = 'Column Over-strength Factor [N/N]'
YLABEL = 'MAX. STRUCTURAL PERIOD [s]' 
TITLE = f'{YLABEL} and {XLABEL} DURING PERIOD ANALYSIS'
COLOR = 'brown'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)


XDATA = COL_FY
YDATA = C_Ductility_Rito_penta
XLABEL = 'Column Yield Strength [N]'
YLABEL = 'STRUCTURAL DUCTILITY RATIO [m/m]' 
TITLE = f'{YLABEL} and {XLABEL} DURING PERIOD ANALYSIS'
COLOR = 'blue'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = COL_KE
YDATA = C_Ductility_Rito_penta
XLABEL = 'Column Elastic Stiffness [N/m]'
YLABEL = 'STRUCTURAL DUCTILITY RATIO [m/m]' 
TITLE = f'{YLABEL} and {XLABEL} DURING PERIOD ANALYSIS'
COLOR = 'blue'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = COL_OSF
YDATA = C_Ductility_Rito_penta
XLABEL = 'Column Over-strength Factor [N/N]'
YLABEL = 'STRUCTURAL DUCTILITY RATIO [m/m]' 
TITLE = f'{YLABEL} and {XLABEL} DURING PERIOD ANALYSIS'
COLOR = 'blue'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = COL_DUCT
YDATA = C_Ductility_Rito_penta
XLABEL = 'Column uctility Ratio [m/m]'
YLABEL = 'STRUCTURAL DUCTILITY RATIO [m/m]' 
TITLE = f'{YLABEL} and {XLABEL} DURING PERIOD ANALYSIS'
COLOR = 'blue'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)


XDATA = COL_FY
YDATA = C_Over_Strength_Factor_penta
XLABEL = 'Column Yield Strength [N]'
YLABEL = 'OVER-STRENGTH FACTOR [N/N]'  
TITLE = f'{YLABEL} and {XLABEL} DURING PERIOD ANALYSIS'
COLOR = 'blue'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = COL_KE
YDATA = C_Over_Strength_Factor_penta
XLABEL = 'Column Elastic Stiffness [N/m]'
YLABEL = 'OVER-STRENGTH FACTOR [N/N]'  
TITLE = f'{YLABEL} and {XLABEL} DURING PERIOD ANALYSIS'
COLOR = 'blue'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = COL_OSF
YDATA = C_Over_Strength_Factor_penta
XLABEL = 'Column Over-strength Factor [N/N]'
YLABEL = 'OVER-STRENGTH FACTOR [N/N]' 
TITLE = f'{YLABEL} and {XLABEL} DURING PERIOD ANALYSIS'
COLOR = 'blue'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

XDATA = COL_DUCT
YDATA = C_Over_Strength_Factor_penta
XLABEL = 'Column uctility Ratio [m/m]'
YLABEL = 'OVER-STRENGTH FACTOR [N/N]' 
TITLE = f'{YLABEL} and {XLABEL} DURING PERIOD ANALYSIS'
COLOR = 'blue'
S09.PLOT_SCATTER(XDATA, YDATA , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)



# 3D PLOT
X, Y, Z = COL_FY, SDOF_ef_DISP_PUSH, SDOF_ef_MASS_PUSH
XLABEL, YLABEL, ZLABEL = 'Column Yield Strength [N]', 'Effective Displacement [m]', 'Effective Mass [kg]',               
S11.PLOT_CONTOUR_3D_2D_FUN(120, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = COL_KE, SDOF_ef_DISP_PUSH, SDOF_ef_MASS_PUSH
XLABEL, YLABEL, ZLABEL = 'Column Elastic Stiffness [N/m]', 'Effective Displacement [m]', 'Effective Mass [kg]',               
S11.PLOT_CONTOUR_3D_2D_FUN(120, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = COL_DUCT, SDOF_ef_DISP_PUSH, SDOF_ef_MASS_PUSH
XLABEL, YLABEL, ZLABEL = 'Column Ductility Ratio [m/m]', 'Effective Displacement [m]', 'Effective Mass [kg]',               
S11.PLOT_CONTOUR_3D_2D_FUN(120, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = COL_OSF, SDOF_ef_DISP_PUSH, SDOF_ef_MASS_PUSH
XLABEL, YLABEL, ZLABEL = 'Column Over-strength Factor [N/N]', 'Effective Displacement [m]', 'Effective Mass [kg]',               
S11.PLOT_CONTOUR_3D_2D_FUN(120, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = COL_DUCT, SDOF_ef_DISP_PUSH, SDOF_ef_STIFF_PUSH
XLABEL, YLABEL, ZLABEL = 'Column Ductility Ratio [m/m]', 'Effective Displacement [m]', 'Effective Stiffness [N/m]',               
S11.PLOT_CONTOUR_3D_2D_FUN(121, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = COL_DUCT, SDOF_ef_DISP_PUSH, SDOF_ef_PERIOD_PUSH
XLABEL, YLABEL, ZLABEL = 'Column Ductility Ratio [m/m]', 'Effective Displacement [m]', 'Effective Period [s]',               
S11.PLOT_CONTOUR_3D_2D_FUN(122, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = COL_DUCT, SDOF_ef_DISP_PUSH, C_Ductility_Rito_penta
XLABEL, YLABEL, ZLABEL = 'Column Ductility Ratio [m/m]', 'Effective Displacement [m]', 'structural Ductility Ratio [m/m]',               
S11.PLOT_CONTOUR_3D_2D_FUN(123, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = COL_DUCT, SDOF_ef_DISP_PUSH, C_Over_Strength_Factor_penta
XLABEL, YLABEL, ZLABEL = 'Column Ductility Ratio [m/m]', 'Effective Displacement [m]', 'structural Over-strength Factor [N/N]',               
S11.PLOT_CONTOUR_3D_2D_FUN(124, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X, Y, Z = COL_DUCT, SDOF_ef_DISP_PUSH, PERIOD_PUSH
XLABEL, YLABEL, ZLABEL = 'Column Ductility Ratio [m/m]', 'Effective Displacement [m]', 'structural Max. Period [s]',               
S11.PLOT_CONTOUR_3D_2D_FUN(124, X, Y, Z, XLABEL, YLABEL, ZLABEL)

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
    "COL_FY":                            COL_FY,
    "COL_KE":                            COL_KE,        
    "COL_DUCT":                          COL_DUCT,
    "COL_OSF":                           COL_OSF,
    "SDOF_ef_DISP_PUSH":                 SDOF_ef_DISP_PUSH,
    "SDOF_ef_MASS_PUSH":                 SDOF_ef_MASS_PUSH,
    "SDOF_ef_STIFF_PUSH":                SDOF_ef_STIFF_PUSH,
    "SDOF_ef_PERIOD_PUSH":               SDOF_ef_PERIOD_PUSH,
    "C_Ductility_Rito_penta":            C_Ductility_Rito_penta,
    "C_Over_Strength_Factor_penta":      C_Over_Strength_Factor_penta,
}


# Convert to DataFrame
df = pd.DataFrame(data)
#print(df)
threshold_ductility_ratio = 35.0
S09.RANDOM_FOREST(df, threshold_ductility_ratio)
#%%------------------------------------------------------
# PLOT HEATMAP FOR CORRELATION 
S09.PLOT_HEATMAP(df)
#%%------------------------------------------------------
# MULTIPLE REGRESSION MODEL
#S01.MULTIPLE_REGRESSION(df) 
#%%-------------------------------------------------------------------
# Plots a heatmap of sensitivity coefficients (correlation or SRC) between inputs X and outputs Y.
import SENSITIVITY_HEATMAP_FUN as S099
X = np.column_stack([COL_FY, COL_KE, COL_DUCT, COL_OSF])          
Y = np.column_stack([SDOF_ef_MASS_PUSH, SDOF_ef_PERIOD_PUSH, C_Ductility_Rito_penta, C_Over_Strength_Factor_penta])  
X_LABELS = ['Column Yield Strength [N]','Column Elastic Stiffness [N/m]','Column Ductility Ratio [m/m]', 'Column Over-srength Factor [N/N]']
Y_LABELS = ['Effective Mass [kg]','Effective Period [s]', 'Structural Ductility Ratio [m/m]', 'Structural Over-srength Factor [N/N]']
coeffs = S099.SENSITIVITY_HEATMAP_FUN(X, Y, X_LABELS, Y_LABELS, method='pearson')
#%%-------------------------------------------------------------------
# ANOVA with automatic binning of continuous predictors.
import ANOVA_SENSITIVITY_FUN as S100

df_sens = pd.DataFrame({
    "COL_FY":                            COL_FY,
    "COL_KE":                            COL_KE,        
    "COL_DUCT":                          COL_DUCT,
    "COL_OSF":                           COL_OSF,
    "SDOF_ef_DISP_PUSH":                 SDOF_ef_DISP_PUSH,
    "SDOF_ef_MASS_PUSH":                 SDOF_ef_MASS_PUSH,
    "SDOF_ef_STIFF_PUSH":                SDOF_ef_STIFF_PUSH,
    "SDOF_ef_PERIOD_PUSH":               SDOF_ef_PERIOD_PUSH,
    "C_Ductility_Rito_penta":            C_Ductility_Rito_penta,
    "C_Over_Strength_Factor_penta":      C_Over_Strength_Factor_penta,
})

param_list = ["COL_FY", "COL_KE", "COL_OSF", "SDOF_ef_MASS_PUSH", "SDOF_ef_PERIOD_PUSH", "C_Ductility_Rito_penta", "C_Over_Strength_Factor_penta"]

# ANOVA – main effects only
anova_table, _ = S100.ANOVA_SENSITIVITY_FUN(
    df_sens, output_col="COL_DUCT",
    param_cols=param_list,
    n_bins=3,
    include_interactions=False,
    plot=True,
)
plt.show()
print("\nANOVA table (main effects):")
print(anova_table.round(4))

# ANOVA – with two-way interactions
anova_table_inter, _ = S100.ANOVA_SENSITIVITY_FUN(
    df_sens, output_col="COL_DUCT",
    param_cols=param_list,
    n_bins=4,
    include_interactions=True,
    plot=True,
)
plt.show()
print("\nANOVA table (with interactions):")
print(anova_table_inter.round(4))

#%%----------------------------------------------------
# PLOT THE EQUIVALENT SDOF SYSTEM DERIVATION VIA DISPLACEMENT-BASED SEISMIC DESIGN PROCEDURE WITH PUSHOVER ANALYSIS
# Prepare data
x = np.asarray(COL_DUCT).ravel()

def aligned(vals, x_ref):
    """
    Truncate `vals` to the common length with `x_ref`, then reorder it
    with the same sort key used for `x_ref`.
    """
    v = np.asarray(vals).ravel()
    m = min(len(v), len(x_ref))          # safe common length
    x_cut = x_ref[:m]
    v_cut = v[:m]
    order = np.argsort(x_cut)            # reorder indices stay within bounds
    return x_cut[order], v_cut[order]

# Rebuild x_s together with each series, using the same truncation rule
x_s, y_disp = aligned(SDOF_ef_DISP_PUSH,  x)
_,   y_mass = aligned(SDOF_ef_MASS_PUSH,  x)
_,   y_stiff= aligned(SDOF_ef_STIFF_PUSH, x)
_,   y_per  = aligned(SDOF_ef_PERIOD_PUSH,x)

series = {
    'Effective Displacement [m]' : y_disp,
    'Effective Mass [kg]'        : y_mass,
    'Effective Stiffness [N/m]'  : y_stiff,
    'Effective Period [s]'       : y_per,
}

# ----------------------------------------------------------------------
# Plot: one subplot per quantity
# ----------------------------------------------------------------------
n = len(series)
fig, axes = plt.subplots(1, n, figsize=(5 * n, 4.5), sharex=True)
if n == 1:
    axes = [axes]

colors = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3']

for ax, (name, y), c in zip(axes, series.items(), colors):
    ax.plot(x_s, y, '-o', color=c, lw=1.8, ms=6, mfc='white', mew=1.5)

    # annotate each point with its value
    for xi, yi in zip(x_s, y):
        ax.annotate(f'{yi:.3g}',
                    xy=(xi, yi), xytext=(0, 6),
                    textcoords='offset points',
                    ha='center', va='bottom', fontsize=8, color='black')

    ax.set_title(name, fontsize=11)
    ax.set_xlabel('Column Ductility Ratio [m/m]')
    ax.grid(True, which='both', ls=':', alpha=0.4)

    # headroom so labels aren't clipped
    ymin, ymax = np.min(y), np.max(y)
    span = ymax - ymin if ymax > ymin else abs(ymax) or 1
    ax.set_ylim(ymin - 0.08 * span, ymax + 0.15 * span)

axes[0].set_ylabel('Value')

#fig.suptitle('EQUIVALENT SDOF SYSTEM DERIVATION VIA DISPLACEMENT-BASED SEISMIC DESIGN PROCEDURE WITH PUSHOVER ANALYSIS \n SDOF Effective Parameters vs Column Ductility Ratio',fontsize=13, y=1.03)
plt.tight_layout()
plt.show()
#%%----------------------------------------------------
# PLOT THE SENSITIVITY ANALYSIS BY CHANGING EACH COLUMN DUCTILITY RATIO
x = np.asarray(COL_DUCT)

# --- Style dictionary so each idealization is consistently colored ---
style = {
    'Bilinear':     dict(color='#e41a1c', marker='o', ls='-',  lw=1.8, ms=6),
    'Trilinear':    dict(color='#377eb8', marker='s', ls='--', lw=1.8, ms=6),
    'Quadrilinear': dict(color='#4daf4a', marker='^', ls='-.', lw=1.8, ms=6),
    'Pentalinear':  dict(color='#984ea3', marker='D', ls=':',  lw=1.8, ms=6),
}

# Each metric family -> {label: list-of-values}
metric_families = {
    'Structure Elastic Stiffness': [
        ('Bilinear',     C_Elastic_ST_bi),
        ('Trilinear',    C_Elastic_ST_tri),
        ('Quadrilinear', C_Elastic_ST_quad),
        ('Pentalinear',  C_Elastic_ST_penta),
    ],
    'Structure Plastic Stiffness': [
        ('Bilinear',     C_Plastic_ST_bi),
        ('Trilinear',    C_Plastic_ST_tri),
        ('Quadrilinear', C_Plastic_ST_quad),
        ('Pentalinear',  C_Plastic_ST_penta),
    ],
    'Structure Ductility Ratio': [
        ('Bilinear',     C_Ductility_Rito_bi),
        ('Trilinear',    C_Ductility_Rito_tri),
        ('Quadrilinear', C_Ductility_Rito_quad),
        ('Pentalinear',  C_Ductility_Rito_penta),
    ],
    'Structure Over-Strength Factor': [
        ('Bilinear',     C_Over_Strength_Factor_bi),
        ('Trilinear',    C_Over_Strength_Factor_tri),
        ('Quadrilinear', C_Over_Strength_Factor_quad),
        ('Pentalinear',  C_Over_Strength_Factor_penta),
    ],
}
"""
# ---------------------------------------
# Layout: one subplot per metrix family
# ---------------------------------------
n = len(metric_families)
fig, axes = plt.subplots(1, n, figsize=(4.2 * n, 4.6), sharex=True)
if n == 1:
    axes = [axes]

for ax, (title, series_list) in zip(axes, metric_families.items()):
    for name, vals in series_list:
        vals = np.asarray(vals)
        # guard: values list may be shorter/longer than COL_DUCT
        m = min(len(x), len(vals))
        ax.plot(x[:m], vals[:m], label=name, **style[name])

    ax.set_title(title, fontsize=11)
    ax.set_xlabel('Column Ductility Ratio [m/m]')
    ax.grid(True, which='both', ls=':', alpha=0.4)
    #ax.semilogy()

axes[0].set_ylabel('Value')
axes[-1].legend(loc='best', fontsize=9, framealpha=0.9)

fig.suptitle('Statistical Metrics vs Column Ductility Ratio', fontsize=13, y=1.02)
plt.tight_layout()
plt.show()
"""
"""
# ---------------------------------------
# Layout: one subplot per metric family
# ---------------------------------------
n = len(metric_families)
fig, axes = plt.subplots(1, n, figsize=(4.2 * n, 4.6), sharex=True)
if n == 1:
    axes = [axes]

for ax, (title, series_list) in zip(axes, metric_families.items()):
    for name, vals in series_list:
        vals = np.asarray(vals)
        # Guard: values list may be shorter/longer than x
        m = min(len(x), len(vals))
        
        # Plot only points using scatter
        ax.scatter(x[:m], vals[:m],
                   label=name,
                   s=25,                # marker size
                   **{k: v for k, v in style[name].items()
                      if k not in ['linestyle', 'marker', 'markersize', 'linewidth']})
    
    ax.set_title(title, fontsize=11)
    ax.set_xlabel('Column Ductility Ratio [m/m]')
    ax.grid(True, which='both', ls=':', alpha=0.4)
    # ax.semilogy()

axes[0].set_ylabel('Value')
axes[-1].legend(loc='best', fontsize=9, framealpha=0.9)
fig.suptitle('Statistical Metrics vs Column Ductility Ratio', fontsize=13, y=1.02)
plt.tight_layout()
plt.show()
"""
#%%----------------------------------------------------


