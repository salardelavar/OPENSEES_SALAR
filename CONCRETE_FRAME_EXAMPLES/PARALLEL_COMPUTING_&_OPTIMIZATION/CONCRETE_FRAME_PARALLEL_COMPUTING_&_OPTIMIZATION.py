######################################################################################################################
#                               IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL                           #
#                     CONCRETE COLUMN SECTION REBAR OPTIMIZATION BASED ON DEMAND BASE-SHEAR REACTION.                #
# UTILIZING PARALLEL PROCESSING PROCEDURES FOR THE SIMULTANEOUS EXECUTION OF NONLINEAR STATIC AND DYNAMIC CONCRETE   #
#                                             STRUCTURAL ANALYSIS, USING OPENSEES.                                   #
#--------------------------------------------------------------------------------------------------------------------#
#                           OPTIMIZATION ALOGORITHM: NEWTON-RAPHSON METHOD PARALLEL COMPUTING                        #
#--------------------------------------------------------------------------------------------------------------------#
# PARALLEL COMPUTING IS A METHOD OF PERFORMING MULTIPLE CALCULATIONS OR PROCESSES SIMULTANEOUSLY BY DIVIDING A TASK  #
# INTO SMALLER SUB-TASKS. THESE SUB-TASKS RUN CONCURRENTLY ON MULTIPLE PROCESSORS OR CORES TO SPEED UP COMPUTATION.  #
# IT'S COMMONLY USED IN HIGH-PERFORMANCE TASKS LIKE SIMULATIONS, DATA ANALYSIS, AND MACHINE LEARNING.                #
#--------------------------------------------------------------------------------------------------------------------#
#                              THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                            #
#                                       EMAIL: salar.d.ghashghaei@gmail.com                                          #
######################################################################################################################
"""
[1] Nonlinear Frame Modeling: 2D RC frame with distributed plasticity (fiber sections) using `nonlinearBeamColumn` elements.  
[2] Material Laws:  
   - *Concrete*: `Concrete01` with confined (core) and unconfined (cover) properties.  
   - *Steel*: `Hysteretic` model with pinching, hardening, and cyclic degradation.  
[3] Seismic Loads:  
   - Pushover: Displacement-controlled lateral loading to failure.  
   - Dynamic: Uniform excitation with user-defined ground motions (X/Y components).  
[4] Damping: Rayleigh damping (a0, a1) calibrated via eigenvalue analysis (modes 1–2).  
[5] Performance Metrics:  
   - Ductility Ratio (μ): Derived from bilinearized pushover curves.  
   - Overstrength (Ω₀): Yield vs. ultimate capacity.  
[6] Advanced Solver: HHT-α integrator (unconditionally stable) with Newton-Raphson iterations.  
[7] Outputs:  
   - Hysteretic responses (P-M, V-Δ, M-θ).  
   - Time-history plots (displacement, base shear).  
   - Stiffness degradation tracking.  
[8] Validation: Logarithmic decrement method for damping ratio verification.  

[9] The Python multiprocessing library enables you to execute your programs in parallel
     and utilize your system's full processor capacity. Essentially, this library allows
     you to run multiple Python processes simultaneously, as if executing separate Python
     programs independently.
     
Key Innovations:  
- Combines static (pushover) and dynamic analyses for comprehensive fragility assessment.  
- Automated bilinear curve fitting for rapid performance evaluation.  
- Exportable results for FEMA P-58 compliance checks. 
"""
#%%-----------------------------------------------------
# YOUTUBE: Stanford CS149 I Parallel Computing I 2023 I Lecture 1 - Why Parallelism? Why Efficiency?
'https://www.youtube.com/watch?v=V1tINV2-9p4'    
# YOUTUBE: Parallel Computing Explained In 3 Minutes
'https://www.youtube.com/watch?v=q7sgzDH1cR8'  
# YOUTUBE: Parallel Programming with Python
'https://www.youtube.com/watch?v=AG1soUh4-nU'
#%%-----------------------------------------------------  
import openseespy.opensees as ops
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import ANALYSIS_FUNCTION as S02
import CONCRETE_SECTION_FUN as S03
import PLOT_2D as S04
import SALAR_PLOT_FUN as S05
import EXCEL_EXPORT_FUN as S06

# Define materials for nonlinear columns and beam
# Define parameters (units: mm, N)
# ------------------------------------------
# CONCRETE                  tag   f'c        ec0   f'cu        ecu
# Core concrete (confined)
fcC = -27.6         # [N/mm²] Concrete Compressive Strength
ec0C = -0.0045      # [mm/mm] Concrete Compressive Strain
fcUC = -21          # [N/mm²] Concrete Compressive Ultimate Strength
ecuC = -0.015       # [mm/mm] Concrete Compressive Ultimate Strain

# Cover concrete (unconfined)
fcU = -18           # [N/mm²] Concrete Compressive Strength
ec0U = -0.0025      # [mm/mm] Concrete Compressive Strain
fcUU = -2           # [N/mm²] Concrete Compressive Ultimate Strength
ecuU = -0.008       # [mm/mm] Concrete Compressive Ultimate Strain
 
# STEEL
# Reinforcing steel
fy = 400          # [N/mm²] Steel Rebar Yield Strength   
Es = 2e5          # [N/mm²] Modulus of Elasticity
ey = fy/Es        # [mm/mm] Steel Rebar Yield Strain
fu = 1.1818*fy    # [N/mm²] Steel Rebar Ultimate Strength
esu = ey*75.2     # [mm/mm] Steel Rebar Ultimate Strain
Esh = (fu - fy)/(esu - ey)
Bs = Esh / Es

# Column Section
Bc = 500                 # [mm] Depth of the Section 
Hc = 500                 # [mm] Height of the Section  
coverC = 50              # [mm] Concrete Section Cover
DIAc = 25               # [mm] Column Rebar Size Diameter
AsC = np.pi*(DIAc**2)/4  # [mm²] Column Area of Rebar


# Beam Section
Bb = 500                 # [mm] Depth of the Section 
Hb = 300                 # [mm] Height of the Section  
coverB = 50              # [mm] Concrete Section Cover
DIAb = DIAc - 4          # [mm] Beam Rebar Size Diameter
AsB = np.pi*(DIAb**2)/4  # [mm²] Beam Area of Rebar

#%% DEFINE THE ELEMENTS LENGTH
LENGTH_COL = 3000        # [mm] Column Length 
LENGTH_BM = 7000         # [mm] Beam Length 

#%% DEFINE PARAMEETRS FOR NONLINEAR DYNAMIC ANALYSIS
GMfact = 9810    # [mm/s²]standard acceleration of gravity or standard acceleration
SSF_X = 0.0001   # Seismic Acceleration Scale Factor in X Direction
SSF_Y = 0.0001   # Seismic Acceleration Scale Factor in Y Direction
iv0_X = 0.00005  # [mm/s] Initial velocity applied to the node  in X Direction
iv0_Y = 0.00005  # [mm/s] Initial velocity applied to the node  in Y Direction
st_iv0 = 0.0     # [s] Initial velocity applied starting time
SEI = 'X'        # Seismic Direction
DR = 0.05        # Intial Guess for Damping ratio
duration = 15.0  # [s] Total simulation duration
dt = 0.01        # [s] Time step
MASS = 12000     # [kg] Mass on the each column

#%% DEFINE PARAMEETRS FOR NONLINEAR STATIC ANALYSIS 
DMAX = 175       # [mm] Maximum Displacement
DINCR = 0.05     # [mm] Incremental Displacement

#%% DEFINE ANALYSIS PROPERTIES
MAX_ITERATIONS = 20000    # Convergence iteration for test
MAX_TOLERANCE = 1.0e-10   # Convergence tolerance for test
#STEEL_KIND: 1 -> WITHOUT HARDENING AND ULTIMATE STRAIN
#STEEL_KIND: 2 -> WITH HARDENING AND ULTIMATE STRAIN

#%%------------------------------------------------------------------------------
def PD_ANALYSIS(DIAc, STEEL_KIND, ANA_KIND):
    
    #DIAc = 25               # [mm] Column Rebar Size Diameter
    AsC = np.pi*(DIAc**2)/4  # [mm²] Column Area of Rebar
    
    DIAb = DIAc - 4          # [mm] Beam Rebar Size Diameter
    AsB = np.pi*(DIAb**2)/4  # [mm²] Beam Area of Rebar
    
    # Reset model
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)

    # CORNER NODES
    ops.node(1, 0.0, 0.0)
    ops.node(2, LENGTH_BM, 0.0)
    ops.node(3, 0.0, LENGTH_COL)
    ops.node(4, LENGTH_BM, LENGTH_COL)

    # BEAM
    ops.node(11, 0.25*LENGTH_BM, LENGTH_COL)
    ops.node(12, 0.50*LENGTH_BM, LENGTH_COL)
    ops.node(13, 0.75*LENGTH_BM, LENGTH_COL)
    
    
    
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 1, 1, 1)

    secTagC = 10
    secTagB = 20
    coreTag = 1
    coverTag = 2
    steelTag = 3
    
    if STEEL_KIND == 1:# WITHOUT HARDENING AND ULTIMATE STRAIN
        ops.uniaxialMaterial('Steel01', steelTag, fy, Es, 0.0) 
    if STEEL_KIND == 2:# WITH HARDENING AND ULTIMATE STRAIN    
        pinchX = 0.8   # Pinching factor in X direction
        pinchY = 0.5   # Pinching factor in Y direction
        damage1 = 0.0  # Damage due to ductility
        damage2 = 0.0  # Damage due to energy
        beta = 0.1 # Stiffness degradation parameter
        ops.uniaxialMaterial('Hysteretic', steelTag, fy, ey, fu, esu, 0.2*fu, 1.1*esu, -fy, -ey, -fu, -esu, -0.2*fu, -1.1*esu, pinchX, pinchY, damage1, damage2, beta)
        # INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Hysteretic_Material


    ops.uniaxialMaterial('Concrete01', coreTag, fcC, ec0C, fcUC, ecuC)  # Core concrete (confined)
    ops.uniaxialMaterial('Concrete01', coverTag, fcU, ec0U, fcUU, ecuU) # Cover concrete (unconfined)
    
    # COLUMN SECTION
    S03.CONFINED_CONCRETE_SECTION(secTagC, Hc, Bc, coverC, AsC, coreTag, coverTag, steelTag, COL=True)
    # BEAM SECTION
    S03.CONFINED_CONCRETE_SECTION(secTagB, Hb, Bb, coverB, AsB, coreTag, coverTag, steelTag, COL=False)
    
    # Define geometric transformation (corotational for large displacements)
    transfTag = 1
    #ops.geomTransf('Linear', transfTag)
    #ops.geomTransf('PDelta', transfTag)
    ops.geomTransf('Corotational', transfTag)
    numIntgrPts = 5
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/nonlinearBeamColumn.html
    # LEFTH COLUMN
    ops.element('nonlinearBeamColumn', 1, 1, 3, numIntgrPts, secTagC, transfTag, '-mass', Bc*Hc*0.000025) # COLUMN 01
    # RIGHT COLUMN
    ops.element('nonlinearBeamColumn', 2, 2, 4, numIntgrPts, secTagC, transfTag, '-mass', Bc*Hc*0.000025) # COLUMN 02
    # BEAM
    ops.element('nonlinearBeamColumn', 9, 3, 11, numIntgrPts, secTagB, transfTag, '-mass', Bb*Hb*0.000025)   # BEAM 01
    ops.element('nonlinearBeamColumn', 10, 11, 12, numIntgrPts, secTagB, transfTag, '-mass', Bb*Hb*0.000025) # BEAM 01
    ops.element('nonlinearBeamColumn', 11, 12, 13, numIntgrPts, secTagB, transfTag, '-mass', Bb*Hb*0.000025) # BEAM 01
    ops.element('nonlinearBeamColumn', 12, 13, 4, numIntgrPts, secTagB, transfTag, '-mass', Bb*Hb*0.000025)  # BEAM 01

    
    if ANA_KIND == 'PUSHOVER':
        WEIGHT = MASS * GMfact
        # Data storage
        FORCE_S, FORCE_A, MOMENT = [], [], []
        DISP_X, DISP_Y, ROT, DRIFT = [], [], [], []
        KA, KS, KI, STEP = [], [], [], []
    
        # Define time series and load pattern
        ops.timeSeries('Linear', 1)
        ops.pattern('Plain', 1, 1)
        #ops.load(3, 1.0, -WEIGHT, 0.0)
        #ops.load(4, 1.0, -WEIGHT, 0.0)
        ops.load(3, 1.0, -1, 0.0)
        ops.load(4, 1.0, -1, 0.0)
        
        # Total steps
        steps = int(np.abs(DMAX)/np.abs(DINCR))
    
        ops.constraints('Plain')
        ops.numberer('Plain')
        ops.system('BandGeneral')
        ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/normDispIncr.html
        ops.algorithm('Newton')
        ops.analysis('Static')
        
        for step in range(steps):
            
            ops.integrator('DisplacementControl', 3, 1, DINCR) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/displacementControl.html
            ops.integrator('DisplacementControl', 4, 1, DINCR) 
            OK = ops.analyze(1)
            S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS
            # Record results
            ops.reactions()
            S = ops.nodeReaction(1, 1) + ops.nodeReaction(2, 1) # SHEAR BASE REACTION
            A = ops.nodeReaction(1, 2) + ops.nodeReaction(2, 2) # AXIAL BASE REACTION
            M = ops.nodeReaction(1, 3) + ops.nodeReaction(2, 3) # MOMENT BASE REACTION
            #print(rot, M)
            disp_X = ops.nodeDisp(3, 1) # LATERAL DISPLACEMENT IN X FOR NODE 3
            disp_Y = ops.nodeDisp(3, 2) # LATERAL DISPLACEMENT IN Y FOR NODE 3
            rot = ops.nodeDisp(3, 3)    # ROTATION IN Z FOR NODE 3
            FORCE_S.append(S)
            FORCE_A.append(A)
            MOMENT.append(M)
            DISP_X.append(disp_X)
            DISP_Y.append(disp_Y)
            DRIFT.append(disp_X / LENGTH_COL) # STRUCTURAL LATERAL RIFT
            ROT.append(rot)
            KS.append(np.abs(S/disp_X)) # LATERAL STIFFNESS IN X
            KA.append(np.abs(A/disp_Y)) # LATERAL STIFFNESS IN Y
            KI.append(np.abs(M/rot))    # ROTATIONAL STIFFNESS IN Z
            STEP.append(step+1)
            #print(step+1, disp_X, S)
            
        #ops.wipe()    
        print('\n Nonlinear Static Analysis Done.\n')  
        #%% PLOT THE RESULTS FOR NONLINEAR STATIC ANALYSIS    
        #S05.PLOT_PUSHOVER(FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y, ROT, KA, KS, KI, STEP)
        #%% PLOT 2D FRAME SHAPES FOR NONLINEAR STATIC ANALYSIS
        #S04.PLOT_2D_FRAME(deformed_scale=10)  # Adjust scale factor as needed
        #%% EXPORT DATA TO EXCEL FOR DYNAMIC ANALYSIS
        #S06.EXPORT_DATA_STATIC(FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y, ROT, KA, KS, KI, STEP)
        #%% OUTPUT MODEL TO JSON FILE
        # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/printModel.html
        #ops.printModel("-JSON", "-file", "CONCRETE_FRAME_PARALLEL_COMPUTING_PUSHOVER.json")
        #ops.printModel() # Print the Analysis Data
        return FORCE_S # RETURN THE SUPPLY VALUE FOR YOUR OPTIMIZATION PROBLEM
#%%----------------------------------------------------------------------------------------------
               
    if ANA_KIND == 'DYNAMIC':
        # Define mass
        ops.mass(3, MASS, MASS, 0.0)
        ops.mass(4, MASS, MASS, 0.0)
        
        # Static analysis
        ops.timeSeries('Linear', 1)
        ops.pattern('Plain', 1, 1)
        if SEI == 'X':
            ops.load(3, 1.0, 0.0, 0.0)
            ops.load(4, 1.0, 0.0, 0.0)
        if SEI == 'Y': 
            ops.load(3, 0.0, 1.0, 0.0)
            ops.load(4, 0.0, 1.0, 0.0)
        if SEI == 'XY':
            ops.load(3, 1.0, 1.0, 0.0)
            ops.load(4, 1.0, 1.0, 0.0)
        
        # Dynamic analysis
        ops.constraints('Plain')
        ops.numberer('Plain')
        ops.system('BandGeneral')
        ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/normDispIncr.html
        #ops.integrator('Newmark', 0.5, 0.25) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/newmark.html
        alpha = 1;gamma=1.5-alpha; beta=((2-alpha)**2)/4;
        ops.integrator('HHT', alpha, gamma, beta) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/hht.html
        ops.algorithm('Newton') # INFO LINK: https://openseespydoc.readthedocs.io/en/stable/src/algorithm.html
        ops.analysis('Transient')
        
        # Calculate Rayleigh damping factors
        Lambda01 = ops.eigen('-fullGenLapack', 2)  # eigenvalue mode 2
        #Lambda01 = ops.eigen('-genBandArpack', 2) # eigenvalue mode 2
        Omega01 = np.power(max(Lambda01), 0.5)
        Omega02 = np.power(min(Lambda01), 0.5)
        a0 = (2 * Omega01 * Omega02 * DR) / (Omega01 + Omega02) # c = a0 * m : Mass-proportional damping
        a1 = (DR * 2) / (Omega01 + Omega02)   # c = a1 * k : Stiffness-proportional damping
        # Apply Rayleigh damping
        ops.rayleigh(a0, a1, 0, 0)   # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/reyleigh.html
        #ops.rayleigh(0, 0, 2 * DR * Omega01, 0) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/reyleigh.html
        PERIOD_01 = (np.pi * 2) / Omega01 # Structure First Period
        PERIOD_02 = (np.pi * 2) / Omega02 # Structure Second Period
        #print('Structure First Period:  ', PERIOD_01)
        #print('Structure Second Period: ', PERIOD_02) 
        
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
        #print('Seismic Defined Done.')
            
        # Data storage
        FORCE_S, FORCE_A, MOMENT = [], [], []
        DISP_X, DISP_Y, ROT, DRIFT = [], [], [], []
        KA, KS, KI, STEP = [], [], [], []
        time = []
        displacement = []
        velocity_X, velocity_Y = [], []
        acceleration_X, acceleration_Y = [], []
            
        stable = 0
        current_time = 0.0
        
        while stable == 0 and current_time < duration:
            stable = ops.analyze(1, dt)
            S02.ANALYSIS(stable, 1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS
            current_time = ops.getTime()
            time.append(current_time)
            # Record results
            ops.reactions()
            S = ops.nodeReaction(1, 1) + ops.nodeReaction(2, 1) # SHEAR BASE REACTION
            A = ops.nodeReaction(1, 2) + ops.nodeReaction(2, 2) # AXIAL BASE REACTION
            M = ops.nodeReaction(1, 3) + ops.nodeReaction(2, 3) # MOMENT BASE REACTION
            #print(rot, M)
            disp_X = ops.nodeDisp(3, 1) # LATERAL DISPLACEMENT IN X FOR NODE 3
            disp_Y = ops.nodeDisp(3, 2) # LATERAL DISPLACEMENT IN Y FOR NODE 3
            rot = ops.nodeDisp(3, 3)    # ROTATION IN Z FOR NODE 3
            DRIFT.append(disp_X / LENGTH_COL)          # STRUCTURAL LATERAL RIFT
            velocity_X.append(ops.nodeVel(3, 1))       # LATERAL VELOCITY IN X FOR NODE 3
            acceleration_X.append(ops.nodeAccel(3, 1)) # LATERAL ACCELERATION IN X FOR NODE 3
            velocity_Y.append(ops.nodeVel(3, 2))       # LATERAL VELOCITY IN Y FOR NODE 3
            acceleration_Y.append(ops.nodeAccel(3, 2)) # LATERAL ACCELERATION IN Y FOR NODE 3
            FORCE_S.append(S)
            FORCE_A.append(A)
            MOMENT.append(M)
            DISP_X.append(disp_X)
            DISP_Y.append(disp_Y)
            ROT.append(rot)
            KS.append(np.abs(S/disp_X)) # LATERAL STIFFNESS IN X
            KA.append(np.abs(A/disp_Y)) # LATERAL STIFFNESS IN Y
            KI.append(np.abs(M/rot))    # ROTATIONAL STIFFNESS IN Z
            #print(current_time, disp_X, S)
        
        #ops.wipe() 
        print('\n Nonlinear Dynamic Analysis Done.\n')
        # Calculating Damping Ratio Using Logarithmic Decrement Analysis 
        displacement = np.array(DISP_X)
        peaks = np.array([displacement[i] for i in range(1, len(displacement)-1) if displacement[i] > displacement[i-1] and displacement[i] > displacement[i+1]])
        # Natural logarithm
        delta = np.log(peaks[:-1] / peaks[1:])    
        
        #%% PLOT THE RESULTS FOR NONLINEAR DYNAMIC ANALYSIS
        #S05.PLOT_DYNAMIC(time, DISP_X, DISP_Y, velocity_X, velocity_Y, acceleration_X, acceleration_Y, FORCE_S, FORCE_A, MOMENT, ROT, delta) 
        #%% PLOT 2D FRAME SHAPES FOR NONLINEAR DYNAMIC ANALYSIS
        #S04.PLOT_2D_FRAME(deformed_scale=100000)  # Adjust scale factor as needed
        #%% EXPORT DATA TO EXCEL FOR DYNAMIC ANALYSIS
        #S06.EXPORT_DATA_DYNAMIC(time, DISP_X, DISP_Y, velocity_X, velocity_Y, acceleration_X, acceleration_Y, FORCE_S, FORCE_A, MOMENT, ROT, delta, KA, KS, KI, PERIOD_01, PERIOD_02)
        #%% OUTPUT MODEL TO JSON FILE
        # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/printModel.html
        #ops.printModel("-JSON", "-file", "CONCRETE_FRAME_PARALLEL_COMPUTING_DYNAMIC.json")
        #ops.printModel() # Print the Analysis Data
        return FORCE_S # RETURN THE SUPPLY VALUE FOR YOUR OPTIMIZATION PROBLEM
#%%------------------------------------------------------------------------------
   
# ----------------------------------------------------------------------------------------------------
# FIND THE OPTIMUM VALUE (NEWTON-RAPHSON SOLVER FOR OPTIMAL REBAR DIAMETER) WITH PARALLEL COMPUTING
# ----------------------------------------------------------------------------------------------------

import concurrent.futures
from multiprocessing import freeze_support
import time as TI  # Import time module with alias

ESP = 1e-2        # Finite difference derivative Convergence Tolerance
TOLERANCE = 1e-6  # Convergence Tolerance
RESIDUAL = 100    # Convergence Residual 
IT = 0            # Intial Iteration
ITMAX = 100000    # Max. Iteration
DEMAND = 215000.0 # [N] Demand Base Shear Reaction
#DEMAND_NODE = 4  # Demand Displacement node

# Helper function must be defined outside main guard
def Run_Analysis(X):
    """Helper function for parallel execution"""
    FORCE_S = PD_ANALYSIS(X, STEEL_KIND=2, ANA_KIND='PUSHOVER') # -> IF YOU WANT TO RUN NONLINEAR STATIC ANALYSIS FOR YOUR OPTIMIZATION PROBLEM
    #FORCE_S = PD_ANALYSIS(X, STEEL_KIND=2, ANA_KIND='DYNAMIC') # -> IF YOU WANT TO RUN NONLINEAR DYNAMIC ANALYSIS FOR YOUR OPTIMIZATION PROBLEM
    SUPPLY = np.max(np.abs(FORCE_S))
    print(f'SUPPLY: {SUPPLY:.5f}')
    return SUPPLY  # Return SUPPLY value
    
def Optimize_Rebar_Diameter():
    X = DIAc          # [mm] Intial Guess for Column Rebar Diameter
    ESP = 1e-2        # Finite difference derivative Convergence Tolerance
    TOLERANCE = 1e-6  # Convergence Tolerance
    RESIDUAL = 100    # Convergence Residual 
    IT = 0            # Intial Iteration
    ITMAX = 100000    # Max. Iteration
    #DEMAND = 215000.0   # [N] Demand Base Shear Reaction
    

    with concurrent.futures.ProcessPoolExecutor() as executor:
        while RESIDUAL > TOLERANCE and IT < ITMAX:
            params = [X, X - ESP, X + ESP]
            
            futures = [executor.submit(Run_Analysis, param) for param in params]
            SUPPLY, SUPPLYmin, SUPPLYmax = [f.result() for f in futures]

            F = SUPPLY - DEMAND
            Fmin = SUPPLYmin - DEMAND
            Fmax = SUPPLYmax - DEMAND

            DF = (Fmax - Fmin) / (2 * ESP)
            print(DF)
            DX = F / DF if abs(DF) > 1e-12 else 0
            RESIDUAL = abs(DX)
            X -= DX
            IT += 1

            print(f'IT: {IT} - RESIDUAL: {RESIDUAL:.6e} - COLUMN SECTION REBAR DIAMETER {X:.6e}')
            
            if RESIDUAL < TOLERANCE:
                print(f'\t\t Optimum Section Rebar Diameter:   {X:.4f}')
                print(f'\t\t Iteration Counts:                 {IT}')
                print(f'\t\t Convergence Residual:             {RESIDUAL:.10e}')
                break

if __name__ == '__main__':
    current_time = TI.strftime("%H:%M:%S", TI.localtime())
    print(f"Start time (HH:MM:SS): {current_time}\n\n")
    

    freeze_support()  # Required for Windows support
    Optimize_Rebar_Diameter()

    
    current_time = TI.strftime("%H:%M:%S", TI.localtime())
    print(f"Finish time (HH:MM:SS): {current_time}\n\n")    

#%%--------------------------------------------------------------------



    
