################################################################################################################
#                          >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                    #
#  SENSITIVITY ANALYSIS OF CONCRETE STRUCTURES USING PUSHOVER & FREE VIBRATION: EFFECTS OF COLUMN HEIGHT,      #
# REBAR DIAMETER, AND CONCRETE CONFINEMENT IN RC FRAMES ON OUTPUT KEY PARAMETERS FROM NONLINEAR STATIC AND     #
#                                  DYNAMIC ANALYSES USING PYTHON AND OPENSEES                                  #
#--------------------------------------------------------------------------------------------------------------#
#                  FREE VIBRATION ANALYSIS USING INITIAL DISPLACEMENT, VELOCITY, AND ACCELERATION              #
#--------------------------------------------------------------------------------------------------------------#
#                                             PARALLEL COMPUTING VERSION                                       #
#--------------------------------------------------------------------------------------------------------------#
# PARALLEL PROCESSING MEANS RUNNING SEVERAL TASKS AT THE SAME TIME INSTEAD OF ONE AFTER ANOTHER.               #
# IN THE CODE, EACH STEP ANALYSIS WAS CALCULATED IN SEQUENCE,                                                  #
# SO THE CPU WORKED ON ONLY ONE MODE AT ANY MOMENT. IN THE REWRITTEN VERSION, THE JOBLIB LIBRARY ALLOWS        #
# ALL FOUR MODES TO RUN SIMULTANEOUSLY ON DIFFERENT CPU CORES. EACH CORE PROCESSES ONE MODE INDEPENDENTLY,     #
# SO THE TOTAL COMPUTATION TIME BECOMES MUCH SHORTER.                                                          #
#                                                                                                              #
# MODERN COMPUTERS USUALLY HAVE MULTIPLE CORES, FOR EXAMPLE 4, 8, OR EVEN MORE. WHEN WE USE PARALLEL           #
# PROCESSING, WE DIVIDE THE WORKLOAD ACROSS THESE CORES. BECAUSE EACH MODE IS A SEPARATE AND INDEPENDENT       #
# ANALYSIS, THEY ARE PERFECT FOR PARALLEL EXECUTION. INSTEAD OF WAITING FOR MODE 1 TO FINISH BEFORE            #
# STARTING MODE 2, ALL MODES START TOGETHER AND FINISH ALMOST TOGETHER.                                        #
#                                                                                                              #
# IN PRACTICE, THE SPEED IMPROVEMENT DEPENDS ON HOW MANY CORES YOUR CPU HAS. IF YOUR COMPUTER HAS 4 CORES,     #
# THE RUNTIME CAN BE UP TO FOUR TIMES FASTER. IN MANY CASES THE SPEEDUP IS AROUND 3–4 TIMES,                   #
# BECAUSE THERE IS A SMALL OVERHEAD WHEN STARTING PARALLEL TASKS. THE REWRITTEN CODE USES PARALLEL             #
# AND DELAYED TO AUTOMATICALLY SEND EACH MODE TO A DIFFERENT CORE AND THEN COLLECT ALL RESULTS                 #
# IN THE CORRECT ORDER. THIS MAKES THE ANALYSIS MORE EFFICIENT WITHOUT CHANGING THE ENGINEERING RESULTS.       #
#                                                                                                              #
# PARALLEL PROCESSING IS ESPECIALLY HELPFUL IN STRUCTURAL ENGINEERING SIMULATIONS WHERE EACH ANALYSIS          #
# REQUIRES HEAVY NUMERICAL CALCULATION. BY USING ALL AVAILABLE CPU POWER,                                      #
# YOU FINISH THE WORK FASTER AND CAN TEST MORE CASES OR MORE MODELS IN THE SAME AMOUNT OF TIME.                #
#--------------------------------------------------------------------------------------------------------------#
#                           PROGRAM DEVELOPED BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                            #
#                                   CONTACT: salar.d.ghashghaei@gmail.com                                      #
################################################################################################################

"""
1. Objective:
   This code performs a sensitivity analysis on a 2D reinforced concrete frame under free-vibration
   conditions. It evaluates the effects of varying column rebar diameters, section depths,
   and confinement enhancement ratios on structural behavior coefficients using OpenSees.

2. Model Setup:
   A nonlinear 2D frame is modeled with columns and beams using distributed plasticity elements.
   Corotational transformation is applied to account for geometric nonlinearity due to large displacements.

3. Material Modeling:
   Confined and unconfined concrete are modeled using the Modified Kent-Scott-Park formulation.
   Steel reinforcement behavior follows either a bilinear or strain hardening model, depending
   on input parameters.

4. Analysis Types:
   Both static (pushover) and dynamic (free-vibration) analyses are supported. Rayleigh damping
   is calibrated based on modal frequencies for dynamic analysis.

5. Key Outputs:
   The model outputs base shear, top displacement, initial and tangent stiffness, ductility ratios
   , overstrength factors, and seismic response modification coefficients (R-factors).

6. Sensitivity Parameters:
   Column rebar diameters (20–32 mm), confinement enhancement ratios (1.15–1.35), and section
   depths are varied systematically to assess their influence on seismic performance.

7. Bilinear Fitting:
   Pushover curves are processed using a bilinear approximation algorithm to extract elastic/plastic
   stiffness, yield points, ductility capacity, and R-factors.

8. Visualization:
   2D graphs and 3D contour plots visualize relationships between rebar size, confinement level,
   and seismic performance indicators such as stiffness and R-factors.

9. Validation:
   Eigenvalue analysis validates dynamic properties (e.g., periods, mode shapes). Convergence
   checks are embedded to ensure numerical stability and robustness.

10. Applications:
   This tool supports performance-based design, seismic retrofitting, reinforcement optimization,
   and code compliance studies (e.g., FEMA 356, Eurocode 8, Standard No. 2800).

This code combines nonlinear modeling, parametric analysis, and post-processing to quantify
 how key structural and material variables influence the seismic response of RC frames.
"""
"""
Advantages and Disadvantages of Using Free Vibration Analysis Instead of Earthquake Ground Motion
 Analysis – From a Structural and Earthquake Engineering Perspective

Advantages of Free Vibration Analysis:

1. Simpler Modeling and Analysis:
   No need to select, scale, or apply earthquake ground motion records. The analysis setup is faster and easier.

2. Insight into Natural Dynamic Behavior:
   Helps to understand the inherent modal and dynamic characteristics of the structure without the complexity of ground motion.

3. Better for Parametric Studies:
   When analyzing the impact of design variables (e.g., confinement ratio, rebar diameter, or column size), free vibration isolates the structure’s response from external excitations.

4. Avoids Uncertainties Related to Ground Motions:
   Eliminates issues such as selecting appropriate records, scaling methods, and site-specific ground motion properties.

5. Useful for Education and Research:
   Ideal for conceptual studies, teaching structural dynamics, and investigating fundamental system behavior.

Disadvantages Compared to Earthquake Ground Motion Analysis:

1. Unrealistic Excitation:
   Free vibration uses initial displacement or velocity, which does not represent actual earthquake-induced ground accelerations.

2. Cannot Evaluate Real Performance:
   Real design must consider how the structure behaves under real seismic loading, which free vibration cannot simulate accurately.

3. No Accumulated or Multi-Pulse Effects:
   Earthquakes may include multiple strong pulses or long durations, which can cause cumulative damage – not captured in free vibration.

4. Ignores Nonlinear and Hysteretic Behavior:
   Free vibration is typically linear or mildly damped; it cannot capture yielding, stiffness degradation, or energy dissipation from plasticity.

5. Not Recognized for Final Design:
   Building codes (e.g., ASCE 7, Eurocode 8, Iranian 2800) require seismic ground motion analysis for structural design and performance evaluation.

6. Highly Dependent on Initial Conditions:
   Results vary significantly based on chosen initial displacement or velocity, which can be arbitrary or unrealistic.

Summary Table:

| Item                         | Free Vibration Analysis       | Earthquake Ground Motion Analysis  |
|------------------------------|-------------------------------|------------------------------------|
| Purpose                      | Academic, conceptual studies  | Real seismic performance evaluation|
| Suitable for Final Design?   |    No                         |   Yes                              |
| Realistic Seismic Input?     |    No                         |   Yes                              |
| Captures Nonlinear Effects?  |    Limited                    |   Yes                              |
| Parametric Sensitivity       |    High                       |   Moderate                         |

Free vibration analysis is valuable for understanding structural dynamics and performing sensitivity studies.
However, for any realistic seismic assessment, performance-based evaluation, or code-compliant design,
 earthquake ground motion analysis is essential and irreplaceable.

"""
import openseespy.opensees as ops
import matplotlib.pyplot as plt
import numpy as np
import time as TI
import ANALYSIS_FUNCTION as S02
import CONCRETE_SECTION_FUN_TWO as S03
import PLOT_2D as S04
import DAMPING_RATIO_FUN as S05
import RAYLEIGH_DAMPING_FUN as S06
import BILINEAR_CURVE as BC
from joblib import Parallel, delayed
import numpy as np

# Define materials for nonlinear columns and beam
# Define parameters (units: mm, N)
# ------------------------------------------
fc = 25                  # [N/mm²] Unconfined concrete compressive strength

# Column Section
Bc = 500                 # [mm] Depth of the Section 
Hc = 500                 # [mm] Height of the Section  
coverC = 50              # [mm] Concrete Section Cover
DIAc = 25                # [mm] Rebar Size Diameter
AsC = np.pi*(DIAc**2)/4  # [mm²] Area of Rebar

# Beam Section
Bb = 500                 # [mm] Depth of the Section 
Hb = 300                 # [mm] Height of the Section  
coverB = 50              # [mm] Concrete Section Cover
DIAb = 18                # [mm] Rebar Size Diameter
AsB = np.pi*(DIAb**2)/4  # [mm²] Area of Rebar

#%% DEFINE THE ELEMENTS LENGTH
LENGTH_COL = 3000        # [mm] Column Length 
LENGTH_BM = 7000         # [mm] Beam Length 

#%% DEFINE PARAMETERS FOR NONLINEAR DYNAMIC ANALYSIS
GMfact = 9810            # [mm/s²]standard acceleration of gravity or standard acceleration
FV = 'X'                 # Free-Vibration Direction
u0 = 0.031               # [mm] Initial displacement applied to the node
v0 = 0.00015             # [mm/s] Initial velocity
a0 = 0.000065            # [mm/s^2] Initial acceleration
DR = 0.05                # Damping ratio
duration = 20.0          # [s] Total simulation duration
dt = 0.001               # [s] Time step
MASS = 12000             # [kg] Mass on the each column

#%% DEFINE PARAMETERS FOR NONLINEAR STATIC ANALYSIS 
DMAX = 675               # [mm] Maximum Displacement
DINCR = 0.05             # [mm] Incremental Displacement

#%% DEFINE ANALYSIS PROPERTIES
MAX_ITERATIONS = 5000    # Convergence iteration for test
MAX_TOLERANCE = 1.0e-6   # Convergence tolerance for test
#STEEL_KIND: 1 -> WITHOUT HARDENING AND ULTIMATE STRAIN
#STEEL_KIND: 2 -> WITH HARDENING AND ULTIMATE STRAIN
#%%------------------------------------------------------------------------------
def PD_ANALYSIS(Hc, DIAc, Kfc, u0, v0, a0, IA, IU, IV, STEEL_KIND, ANA_KIND):
    
    #DIAc = 25                    # [mm] # Rebar Size Diameter
    AsC = np.pi*(DIAc**2)/4       # [mm²] Area of Rebar
    RO_COL = 100*(8*AsC)/(Hc*Bc)  # Column Rebar Ratio
    
    DIAb = DIAc-4                 # [mm] # Rebar Size Diameter
    AsB = np.pi*(DIAb**2)/4       # [mm²] Area of Rebar
    RO_BE = 100*(6*AsB)/(Hb*Bb)   # Beam Rebar Ratio
    
    # Reset model
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)

    # CORNER NODES
    ops.node(1, 0.0, 0.0)
    ops.node(2, LENGTH_BM, 0.0)
    ops.node(3, 0.0, LENGTH_COL)
    ops.node(4, LENGTH_BM, LENGTH_COL)
    
    # LEFT COLUMN MIDDLE NODES
    ops.node(5, 0.0, (1/3)*LENGTH_COL)
    ops.node(6, 0.0, (2/3)*LENGTH_COL)
    
    # RIGHT COLUMN MIDDLE NODES
    ops.node(7, LENGTH_BM, (1/3)*LENGTH_COL)
    ops.node(8, LENGTH_BM, (2/3)*LENGTH_COL)
    
    # BEAM MIDDLE NODES
    ops.node(9, (1/7)*LENGTH_BM, LENGTH_COL)
    ops.node(10, (6/7)*LENGTH_BM, LENGTH_COL)
    
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 1, 1, 1)

    secTagC01, secTagC02 = 10, 20
    secTagB01, secTagB02 = 30, 40

    # COLUMN SECTION FOR ELEMENT 01, 03, 04, 06
    #Kc = 1.25        # Confinement Enhancement Ratio
    S03.CONFINED_CONCRETE_SECTION(secTagC01, Hc, Bc, coverC, AsC, STEEL_KIND, fc, Kfc, COL=True)
    # COLUMN SECTION FOR ELEMENT 02, 05
    S03.CONFINED_CONCRETE_SECTION(secTagC02, Hc, Bc, coverC, AsC, STEEL_KIND, fc, Kfc-0.05, COL=True)
    # BEAM SECTION FOR ELEMENT 07, 09
    S03.CONFINED_CONCRETE_SECTION(secTagB01, Hb, Bb, coverB, AsB, STEEL_KIND, fc, Kfc-0.10, COL=False)
    # BEAM SECTION FOR ELEMENT 08
    S03.CONFINED_CONCRETE_SECTION(secTagB02, Hb, Bb, coverB, AsB, STEEL_KIND, fc, Kfc-0.15, COL=False)
    
    # Define geometric transformation (corotational for large displacements)
    transfTag = 1
    #ops.geomTransf('Linear', transfTag)
    #ops.geomTransf('PDelta', transfTag)
    ops.geomTransf('Corotational', transfTag)
    numIntgrPts = 10
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/nonlinearBeamColumn.html
    
    # LEFT COLUMN
    ops.element('nonlinearBeamColumn', 1, 1, 5, numIntgrPts, secTagC01, transfTag, '-mass', Bc*Hc*0.000025) # COLUMN - ELEMENT 01
    ops.element('nonlinearBeamColumn', 2, 5, 6, numIntgrPts, secTagC02, transfTag, '-mass', Bc*Hc*0.000025) # COLUMN - ELEMENT 02
    ops.element('nonlinearBeamColumn', 3, 6, 3, numIntgrPts, secTagC01, transfTag, '-mass', Bc*Hc*0.000025) # COLUMN - ELEMENT 03
    # RIGHT COLUMN
    ops.element('nonlinearBeamColumn', 4, 2, 7, numIntgrPts, secTagC01, transfTag, '-mass', Bc*Hc*0.000025) # COLUMN - ELEMENT 04
    ops.element('nonlinearBeamColumn', 5, 7, 8, numIntgrPts, secTagC02, transfTag, '-mass', Bc*Hc*0.000025) # COLUMN - ELEMENT 05
    ops.element('nonlinearBeamColumn', 6, 8, 4, numIntgrPts, secTagC01, transfTag, '-mass', Bc*Hc*0.000025) # COLUMN - ELEMENT 06
    # BEAM
    ops.element('nonlinearBeamColumn', 7, 3, 9, numIntgrPts, secTagB01, transfTag, '-mass', Bb*Hb*0.000025)  # BEAM - ELEMENT 07
    ops.element('nonlinearBeamColumn', 8, 9, 10, numIntgrPts, secTagB02, transfTag, '-mass', Bb*Hb*0.000025) # BEAM - ELEMENT 08
    ops.element('nonlinearBeamColumn', 9, 10, 4, numIntgrPts, secTagB01, transfTag, '-mass', Bb*Hb*0.000025) # BEAM - ELEMENT 09

    if ANA_KIND == 'PUSHOVER':
        WEIGHT = MASS * GMfact
        # Data storage
        FORCE_S, FORCE_A, MOMENT = [], [], []
        DISP_X, DISP_Y, ROT = [], [], []
        KA, KS, KI, STEP = [], [], [], []
        DRIFT = [] 
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
        ops.algorithm('Newton')
        ops.analysis('Static')
        
        for step in range(steps):
            
            ops.integrator('DisplacementControl', 3, 1, DINCR) 
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
            ROT.append(rot)
            KS.append(np.abs(S)/np.abs(disp_X)) # LATERAL STIFFNESS IN X
            KA.append(np.abs(A)/np.abs(disp_Y)) # LATERAL STIFFNESS IN Y
            KI.append(np.abs(M)/np.abs(rot))    # ROTATIONAL STIFFNESS IN Z
            DRIFT.append(100*disp_X/LENGTH_COL)
            STEP.append(step)
            #print(step+1, disp_X, S)
        
        # ---------------------------------------
        #  Plot BaseShear-Displacement Analysis
        # ---------------------------------------
        XX = np.abs(DISP_X); YY = np.abs(FORCE_S); # ABSOLUTE VALUE
        SLOPE_NODE = 10

        DATA = BC.BILNEAR_CURVE(XX, YY, SLOPE_NODE)
        X, Y, Elastic_ST, Plastic_ST, Tangent_ST, Ductility_Rito, Over_Strength_Factor = DATA
        """
        XLABEL = 'Displacement in X [mm]'
        YLABEL = 'Base-Shear [N]'
        LEGEND01 = 'Curve'
        LEGEND02 = 'Bilinear Fitted'
        LEGEND03 = 'Undefined'
        TITLE = f'Last Data of BaseShear-Displacement Analysis - Ductility Ratio: {X[2]/X[1]:.4f} - Over Strength Factor: {Y[2]/Y[1]:.4f}'
        COLOR = 'black'
        BC.PLOT_2D(np.abs(DISP_X), np.abs(FORCE_S), X, Y, X, Y, XLABEL, YLABEL, TITLE, LEGEND01, LEGEND02, LEGEND03, COLOR='black', Z=2) 
        print(f'\t\t Ductility Ratio: {Y[2]/Y[1]:.4f}')
        """

        # Calculate Over Strength Coefficient (Ω0)
        Omega_0 = Y[2] / Y[1]
        # Calculate Displacement Ductility Ratio (μ)
        mu = X[2] / X[1]
        # Calculate Ductility Coefficient (Rμ)
        #R_mu = (2 * mu - 1) ** 0.5
        #R_mu = 1
        R_mu = mu
        # Calculate Structural Behavior Coefficient (R)
        R = Omega_0 * R_mu
        """
        print(f'Over Strength Coefficient (Ω0):      {Omega_0:.4f}')
        print(f'Displacement Ductility Ratio (μ):    {mu:.4f}')
        print(f'Ductility Coefficient (Rμ):          {R_mu:.4f}')
        print(f'Structural Behavior Coefficient (R): {R:.4f}') 
        """
        
        return RO_COL, RO_BE, FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y, ROT, KA, KS, KI, Elastic_ST, Plastic_ST, Tangent_ST, Ductility_Rito, Over_Strength_Factor, R, DRIFT, STEP        
    
    if ANA_KIND == 'DYNAMIC':
        # Define mass
        ops.mass(3, MASS, MASS, 0.0)
        ops.mass(4, MASS, MASS, 0.0)
        
        # Dynamic analysis to apply initial displacement
        ops.timeSeries('Linear', 1)
        ops.pattern('Plain', 1, 1)
        if FV == 'X':
            ops.load(3, 1.0, 0.0, 0.0)
            ops.load(4, 1.0, 0.0, 0.0)
        if FV == 'Y': 
            ops.load(3, 0.0, 1.0, 0.0)
            ops.load(4, 0.0, 1.0, 0.0)
        
        if IU == True and FV == 'X':
            # Define initial displacment
            ops.setNodeDisp(3, 1, u0, '-commit') # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/setNodeDisp.html
        if IV == True and FV == 'X':
            # Define initial velocity
            ops.setNodeVel(3, 1, v0, '-commit')  # INFO LINK: https://openseespydoc.readthedocs.io/en/stable/src/setNodeVel.html
        if IA == True and FV == 'X':
            # Define initial  acceleration
            ops.setNodeAccel(3, 1, a0, '-commit') # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/setNodeAccel.html
            
        if IU == True and FV == 'Y':
            # Define initial displacment
            ops.setNodeDisp(3, 2, u0, '-commit') # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/setNodeDisp.html
        if IV == True and FV == 'Y':
            # Define initial velocity
            ops.setNodeVel(3, 2, v0, '-commit')  # INFO LINK: https://openseespydoc.readthedocs.io/en/stable/src/setNodeVel.html
        if IA == True and FV == 'Y':
            # Define initial  acceleration
            ops.setNodeAccel(3, 2, a0, '-commit') # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/setNodeAccel.html    
    
        
        # Dynamic analysis
        ops.constraints('Plain')
        ops.numberer('Plain')
        ops.system('BandGeneral')
        ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS)
        ops.integrator('Newmark', 0.5, 0.25) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/newmark.html
        #alpha = 1;gamma=1.5-alpha; beta=((2-alpha)**2)/4;
        #ops.integrator('HHT', alpha, gamma, beta) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/hht.html
        ops.algorithm('Newton') # INFO LINK: https://openseespydoc.readthedocs.io/en/stable/src/algorithm.html
        ops.analysis('Transient')
        
        # Calculate Rayleigh damping factors
        PERIOD_01, PERIOD_02 = S06.RAYLEIGH_DAMPING(2, DR, 0.6*DR, 0, 1) 
                    
        # Data storage
        FORCE_S, FORCE_A, MOMENT = [], [], []
        DISP_X, DISP_Y, ROT = [], [], []
        KA, KS, KI, STEP = [], [], [], []
        time = []
        velocity_X, velocity_Y = [], []
        acceleration_X, acceleration_Y = [], []
        DRIFT = [] 
        
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
            disp_X = ops.nodeDisp(3, 1)                # LATERAL DISPLACEMENT IN X FOR NODE 3
            disp_Y = ops.nodeDisp(3, 2)                # LATERAL DISPLACEMENT IN Y FOR NODE 3
            rot = ops.nodeDisp(3, 3)                   # ROTATION IN Z FOR NODE 3
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
            DRIFT.append(100*disp_X/LENGTH_COL)
            #print(current_time, disp_X, S)
           
        damping_ratio = S05.DAMPING_RATIO(DISP_X)
        #ops.wipe()  
        return RO_COL, RO_BE, FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y, ROT, KA, KS, KI, time, velocity_X, velocity_Y, acceleration_X, acceleration_Y, PERIOD_01, PERIOD_02, DRIFT, damping_ratio

#%%--------------------------------------------------------------------------------------------------
# SENSTIVITY ANALYSIS BY CHANGING REBAR DIAMETER AND CONFINEMENT ENHANCEMENT RATIO AND COLUMN DEPTH 
IU = True        # Free Vibration with Initial Displacement
IV = True        # Free Vibration with Initial Velocity
IA = True        # Free Vibration with Initial Acceleration

# STORAGE ARRAYS
RO_COL_MAX, RO_BE_MAX = [], []
FORCE_S_MAX, FORCE_A_MAX, MOMENT_MAX = [], [], []
DISP_X_MAX, DISP_Y_MAX, ROT_MAX = [], [], []
KA_MAX, KS_MAX, KI_MAX, R_MAX = [], [], [], []
Elastic_ST_MAX, Plastic_ST_MAX, Tangent_ST_MAX = [], [], []
Ductility_Rito_MAX, Over_Strength_Factor_MAX = [], []
HCC_MAX, DIAc_MAX, Kfc_MAX = [], [], []
DRIFT_MAX = []

FORCE_S_MAX_d, FORCE_A_MAX_d, MOMENT_MAX_d = [], [], []
DISP_X_MAX_d, DISP_Y_MAX_d, ROT_MAX_d = [], [], []
velocity_X_MAX_d, velocity_Y_MAX_d = [], []
acceleration_X_MAX_d, acceleration_Y_MAX_d = [], []
KA_MAX_d, KS_MAX_d, KI_MAX_d = [], [], []
PERIOD_01_MAX, PERIOD_02_MAX, DAMPING_RATIO_MAX = [], [], []
DRIFT_MAX_d = []

def run_analysis(hc, dia, kf, u0, v0, a0, IA, IU, IV):
    
    # ---------- PUSHOVER ----------
    DATA = PD_ANALYSIS(
        hc, dia, kf, u0, v0, a0,
        IA, IU, IV,
        STEEL_KIND=2,
        ANA_KIND='PUSHOVER'
    )

    (RO_COL, RO_BE, FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y,
     ROT, KA, KS, KI, Elastic_ST, Plastic_ST,
     Tangent_ST, Ductility_Rito, Over_Strength_Factor, R,
     DRIFT, STEP) = DATA

    pushover_results = {
        "RO_COL": RO_COL,
        "RO_BE": RO_BE,
        "FORCE_S": np.max(np.abs(FORCE_S)),
        "FORCE_A": np.max(np.abs(FORCE_A)),
        "MOMENT": np.max(np.abs(MOMENT)),
        "DISP_X": np.max(np.abs(DISP_X)),
        "DISP_Y": np.max(np.abs(DISP_Y)),
        "ROT": np.max(np.abs(ROT)),
        "KA": np.max(np.abs(KA)),
        "KS": np.max(np.abs(KS)),
        "KI": np.max(np.abs(KI)),
        "DRIFT": np.max(np.abs(DRIFT)),
        "Elastic_ST": Elastic_ST,
        "Plastic_ST": Plastic_ST,
        "Tangent_ST": Tangent_ST,
        "Ductility_Rito": Ductility_Rito,
        "Over_Strength_Factor": Over_Strength_Factor,
        "R": R
    }

    # ---------- DYNAMIC ----------
    DATA = PD_ANALYSIS(
        hc, dia, kf, u0, v0, a0,
        IA, IU, IV,
        STEEL_KIND=2,
        ANA_KIND='DYNAMIC'
    )

    (_, _, FORCE_Sd, FORCE_Ad, MOMENTd, DISP_Xd, DISP_Yd, ROTd,
     KAd, KSd, KId, timed,
     velocity_Xd, velocity_Yd,
     acceleration_Xd, acceleration_Yd,
     PERIOD_01, PERIOD_02, DRIFT_d, damping_ratio) = DATA

    dynamic_results = {
        "FORCE_S_d": np.max(np.abs(FORCE_Sd)),
        "FORCE_A_d": np.max(np.abs(FORCE_Ad)),
        "MOMENT_d": np.max(np.abs(MOMENTd)),
        "DISP_X_d": np.max(np.abs(DISP_Xd)),
        "DISP_Y_d": np.max(np.abs(DISP_Yd)),
        "ROT_d": np.max(np.abs(ROTd)),
        "VEL_X_d": np.max(np.abs(velocity_Xd)),
        "VEL_Y_d": np.max(np.abs(velocity_Yd)),
        "ACC_X_d": np.max(np.abs(acceleration_Xd)),
        "ACC_Y_d": np.max(np.abs(acceleration_Yd)),
        "KA_d": np.max(np.abs(KAd)),
        "KS_d": np.max(np.abs(KSd)),
        "KI_d": np.max(np.abs(KId)),
        "DRIFT_d": np.max(np.abs(DRIFT_d)),
        "PERIOD_01": PERIOD_01,
        "PERIOD_02": PERIOD_02,
        "DAMPING_RATIO": damping_ratio
    }

    return hc, dia, kf, pushover_results, dynamic_results


# COLUMN SECTION DEPTH
HCC = [350.0, 400.0, 450.0, 500.0, 550.0, 600.0]
#HCC = [400.0, 500.0, 500.0]
# REBAR DIAMETER FOR COLUMN SECTION
#DIAc = [8, 10, 12, 14, 16, 18, 20, 22, 25, 28, 30, 32]
DIAc = [28, 30, 32]
# CONFINEMENT ENHANCEMENT FACTOR (K)
#Kfc = [1.15, 1.20, 1.25, 1.30, 1.35]
Kfc = [1.25]

params = [(hc, dia, kf) for hc in HCC for dia in DIAc for kf in Kfc]

# ----------------------  PARALLEL PROCESSING  ----------------------
# Analysis Durations:
current_time = TI.strftime("%H:%M:%S", TI.localtime())
print("Start Time:", current_time)

results = Parallel(n_jobs=-1, verbose=10)(
delayed(run_analysis)(hc, dia, kf, u0, v0, a0, IA, IU, IV)
    for hc, dia, kf in params
)


for hc, dia, kf, pushover, dynamic in results:
    
    HCC_MAX.append(hc)
    DIAc_MAX.append(dia)
    Kfc_MAX.append(kf)

    RO_COL_MAX.append(pushover["RO_COL"])
    RO_BE_MAX.append(pushover["RO_BE"])

    FORCE_S_MAX.append(pushover["FORCE_S"])
    FORCE_A_MAX.append(pushover["FORCE_A"])
    MOMENT_MAX.append(pushover["MOMENT"])
    DISP_X_MAX.append(pushover["DISP_X"])
    DISP_Y_MAX.append(pushover["DISP_Y"])
    ROT_MAX.append(pushover["ROT"])

    KA_MAX.append(pushover["KA"])
    KS_MAX.append(pushover["KS"])
    KI_MAX.append(pushover["KI"])
    DRIFT_MAX.append(pushover["DRIFT"])

    Elastic_ST_MAX.append(pushover["Elastic_ST"])
    Plastic_ST_MAX.append(pushover["Plastic_ST"])
    Tangent_ST_MAX.append(pushover["Tangent_ST"])
    Ductility_Rito_MAX.append(pushover["Ductility_Rito"])
    Over_Strength_Factor_MAX.append(pushover["Over_Strength_Factor"])
    R_MAX.append(pushover["R"])

    FORCE_S_MAX_d.append(dynamic["FORCE_S_d"])
    FORCE_A_MAX_d.append(dynamic["FORCE_A_d"])
    MOMENT_MAX_d.append(dynamic["MOMENT_d"])
    DISP_X_MAX_d.append(dynamic["DISP_X_d"])
    DISP_Y_MAX_d.append(dynamic["DISP_Y_d"])
    ROT_MAX_d.append(dynamic["ROT_d"])

    velocity_X_MAX_d.append(dynamic["VEL_X_d"])
    velocity_Y_MAX_d.append(dynamic["VEL_Y_d"])
    acceleration_X_MAX_d.append(dynamic["ACC_X_d"])
    acceleration_Y_MAX_d.append(dynamic["ACC_Y_d"])

    KA_MAX_d.append(dynamic["KA_d"])
    KS_MAX_d.append(dynamic["KS_d"])
    KI_MAX_d.append(dynamic["KI_d"])
    DRIFT_MAX_d.append(dynamic["DRIFT_d"])

    PERIOD_01_MAX.append(dynamic["PERIOD_01"])
    PERIOD_02_MAX.append(dynamic["PERIOD_02"])
    DAMPING_RATIO_MAX.append(dynamic["DAMPING_RATIO"])

current_time = TI.strftime("%H:%M:%S", TI.localtime())
print("Finish Time:", current_time)    

#%%------------------------------------------------------------------------------
def PLOT_3D(TAG, X, Y, Z, XLABEL, YLABEL, ZLABEL):
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.interpolate import griddata
    # Convert to NumPy arrays
    X = np.array(X)
    Y = np.array(Y)
    Z = np.array(Z)
    
    # Create grid for contour plot
    xi = np.linspace(min(X), max(X), 100)
    yi = np.linspace(min(Y), max(Y), 100)
    xi, yi = np.meshgrid(xi, yi)
    
    # Interpolate Z values on the grid
    zi = griddata((X, Y), Z, (xi, yi), method='cubic')
    
    # Plot 3D contour
    fig = plt.figure(TAG, figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')
    contour = ax.plot_surface(xi, yi, zi, cmap='viridis', edgecolor='none')
    
    ax.set_xlabel(XLABEL)
    ax.set_ylabel(YLABEL)
    ax.set_zlabel(ZLABEL,)
    
    fig.colorbar(contour, ax=ax, shrink=0.5, aspect=5)
    plt.title(f'3D Contour Plot of {ZLABEL}')
    plt.show()

#%%------------------------------------------------------------------------------
X = RO_COL_MAX
Y = RO_BE_MAX
Z = Kfc_MAX
XLABEL = 'Column Steel Reinforcement Ratio [%]'   
YLABEL = 'Beam Steel Reinforcement Ratio [%]' 
ZLABEL = 'Confinement Enhancement Ratio' 
PLOT_3D(0, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = RO_COL_MAX
Y = RO_BE_MAX
Z = FORCE_S_MAX
XLABEL = 'Column Steel Reinforcement Ratio [%]'   
YLABEL = 'Beam Steel Reinforcement Ratio [%]' 
ZLABEL = 'Structure Shear Reaction in X Dir. [N] from Pushover Analysis'
PLOT_3D(1, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = RO_COL_MAX
Y = RO_BE_MAX
Z = FORCE_S_MAX_d
XLABEL = 'Column Steel Reinforcement Ratio [%]'   
YLABEL = 'Beam Steel Reinforcement Ratio [%]' 
ZLABEL = 'Structure Shear Reaction in X Dir. [N] from Dynamic Analysis'
PLOT_3D(2, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = RO_COL_MAX
Y = RO_BE_MAX
Z = np.divide(FORCE_S_MAX_d, FORCE_S_MAX)
XLABEL = 'Column Steel Reinforcement Ratio [%]'   
YLABEL = 'Beam Steel Reinforcement Ratio [%]' 
ZLABEL = 'Structure Shear Reaction Ratio in X Dir. '
PLOT_3D(3, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = RO_COL_MAX
Y = RO_BE_MAX
Z =  FORCE_A_MAX
XLABEL = 'Column Steel Reinforcement Ratio [%]'   
YLABEL = 'Beam Steel Reinforcement Ratio [%]' 
ZLABEL = 'Structure Axial Reaction in Y Dir. [N] from Pushover Analysis'
PLOT_3D(4, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = RO_COL_MAX
Y = RO_BE_MAX
Z =  FORCE_A_MAX_d
XLABEL = 'Column Steel Reinforcement Ratio [%]'   
YLABEL = 'Beam Steel Reinforcement Ratio [%]' 
ZLABEL = 'Structure Axial Reaction in Y Dir. [N] from Dynamic Analysis'
PLOT_3D(5, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = RO_COL_MAX
Y = RO_BE_MAX
Z =  np.divide(FORCE_A_MAX_d, FORCE_A_MAX)
XLABEL = 'Column Steel Reinforcement Ratio [%]'   
YLABEL = 'Beam Steel Reinforcement Ratio [%]' 
ZLABEL = 'Structure Axial Reaction Ratio in Y Dir.'
PLOT_3D(6, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = RO_COL_MAX
Y = RO_BE_MAX
Z =  MOMENT_MAX
XLABEL = 'Column Steel Reinforcement Ratio [%]'   
YLABEL = 'Beam Steel Reinforcement Ratio [%]' 
ZLABEL = 'Structure Moment Reaction in Y Dir. [N.mm] from Pushover Analysis'
PLOT_3D(7, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = RO_COL_MAX
Y = RO_BE_MAX
Z =  MOMENT_MAX_d
XLABEL = 'Column Steel Reinforcement Ratio [%]'   
YLABEL = 'Beam Steel Reinforcement Ratio [%]' 
ZLABEL = 'Structure Moment Reaction in Y Dir. [N.mm] from Dynamic Analysis'
PLOT_3D(8, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = RO_COL_MAX
Y = RO_BE_MAX
Z =  np.divide(MOMENT_MAX_d, MOMENT_MAX)
XLABEL = 'Column Steel Reinforcement Ratio [%]'   
YLABEL = 'Beam Steel Reinforcement Ratio [%]' 
ZLABEL = 'Structure Moment Reaction Ratio in Y Dir.'
PLOT_3D(9, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = RO_COL_MAX
Y = RO_BE_MAX
Z =  DISP_X_MAX
XLABEL = 'Column Steel Reinforcement Ratio [%]'   
YLABEL = 'Beam Steel Reinforcement Ratio [%]'  
ZLABEL = 'Structure Displacement in X Dir. [mm] from Pushover Analysis'
PLOT_3D(10, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = RO_COL_MAX
Y = RO_BE_MAX
Z =  DISP_X_MAX_d
XLABEL = 'Column Steel Reinforcement Ratio [%]'   
YLABEL = 'Beam Steel Reinforcement Ratio [%]' 
ZLABEL = 'Structure Displacement in X Dir. [mm] from Dynamic Analysis'
PLOT_3D(11, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = RO_COL_MAX
Y = RO_BE_MAX
Z =  np.divide(DISP_X_MAX_d, DISP_X_MAX)
XLABEL = 'Column Steel Reinforcement Ratio [%]'   
YLABEL = 'Beam Steel Reinforcement Ratio [%]' 
ZLABEL = 'Structure Displacement Ratio in X Dir.'
PLOT_3D(12, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = RO_COL_MAX
Y = RO_BE_MAX
Z =  DISP_Y_MAX
XLABEL = 'Column Steel Reinforcement Ratio [%]'   
YLABEL = 'Beam Steel Reinforcement Ratio [%]' 
ZLABEL = 'Structure Displacement in Y Dir. [mm] from Pushover Analysis'
PLOT_3D(13, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = RO_COL_MAX
Y = RO_BE_MAX
Z =  DISP_Y_MAX_d
XLABEL = 'Column Steel Reinforcement Ratio [%]'   
YLABEL = 'Beam Steel Reinforcement Ratio [%]' 
ZLABEL = 'Structure Displacement in Y Dir. [mm] from Dynamic Analysis'
PLOT_3D(14, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = RO_COL_MAX
Y = RO_BE_MAX
Z =  np.divide(DISP_Y_MAX_d, DISP_Y_MAX)
XLABEL = 'Column Steel Reinforcement Ratio [%]'   
YLABEL = 'Beam Steel Reinforcement Ratio [%]' 
ZLABEL = 'Structure Displacement Ratio in Y Dir.'
PLOT_3D(15, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = RO_COL_MAX
Y = RO_BE_MAX
Z =  KS_MAX
XLABEL = 'Column Steel Reinforcement Ratio [%]'   
YLABEL = 'Beam Steel Reinforcement Ratio [%]' 
ZLABEL = 'Lateral Stiffness in X Dir. [N/mm] from Pushover Analysis'
PLOT_3D(16, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = RO_COL_MAX
Y = RO_BE_MAX
Z =  KS_MAX_d
XLABEL = 'Column Steel Reinforcement Ratio [%]'   
YLABEL = 'Beam Steel Reinforcement Ratio [%]' 
ZLABEL = 'Lateral Stiffness in X Dir. [N/mm] from Dynamic Analysis'
PLOT_3D(17, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = RO_COL_MAX
Y = RO_BE_MAX
Z =  np.divide(KS_MAX_d, KS_MAX)
XLABEL = 'Column Steel Reinforcement Ratio [%]'   
YLABEL = 'Beam Steel Reinforcement Ratio [%]' 
ZLABEL = 'Lateral Stiffness Ratio in X Dir.'
PLOT_3D(18, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = RO_COL_MAX
Y = RO_BE_MAX
Z =  KA_MAX
XLABEL = 'Column Steel Reinforcement Ratio [%]'   
YLABEL = 'Beam Steel Reinforcement Ratio [%]' 
ZLABEL = 'Lateral Stiffness in Y Dir. [N/mm] from Pushover Analysis'
PLOT_3D(19, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = RO_COL_MAX
Y = RO_BE_MAX
Z =  KA_MAX_d
XLABEL = 'Column Steel Reinforcement Ratio [%]'   
YLABEL = 'Beam Steel Reinforcement Ratio [%]' 
ZLABEL = 'Lateral Stiffness in Y Dir. [N/mm] from Dynamic Analysis'
PLOT_3D(20, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = RO_COL_MAX
Y = RO_BE_MAX
Z =  np.divide(KA_MAX_d, KA_MAX)
XLABEL = 'Column Steel Reinforcement Ratio [%]'   
YLABEL = 'Beam Steel Reinforcement Ratio [%]' 
ZLABEL = 'Lateral Stiffness Ratio in Y Dir.'
PLOT_3D(21, X, Y, Z, XLABEL, YLABEL, ZLABEL)


X = RO_COL_MAX
Y = RO_BE_MAX
Z =  Elastic_ST_MAX
XLABEL = 'Rebar Diameter'   
YLABEL = 'Confinement Enhancement Ratio' 
ZLABEL = 'Structure Elastic Stiffness [N/mm]'
PLOT_3D(22, X, Y, Z, XLABEL, YLABEL, ZLABEL)    

X = RO_COL_MAX
Y = RO_BE_MAX
Z =  Plastic_ST_MAX
XLABEL = 'Column Steel Reinforcement Ratio [%]'   
YLABEL = 'Beam Steel Reinforcement Ratio [%]'  
ZLABEL = 'Structure Plastic Stiffness [N/mm]'
PLOT_3D(23, X, Y, Z, XLABEL, YLABEL, ZLABEL)    

X = RO_COL_MAX
Y = RO_BE_MAX
Z =  Ductility_Rito_MAX
XLABEL = 'Column Steel Reinforcement Ratio [%]'   
YLABEL = 'Beam Steel Reinforcement Ratio [%]' 
ZLABEL = 'Ductility Rito  [mm/mm]'
PLOT_3D(24, X, Y, Z, XLABEL, YLABEL, ZLABEL)   

X = RO_COL_MAX
Y = RO_BE_MAX
Z =  Over_Strength_Factor_MAX
XLABEL = 'Column Steel Reinforcement Ratio [%]'   
YLABEL = 'Beam Steel Reinforcement Ratio [%]' 
ZLABEL = 'Over Strength Factor [N/N]'
PLOT_3D(25, X, Y, Z, XLABEL, YLABEL, ZLABEL)   

X = RO_COL_MAX
Y = RO_BE_MAX
Z = R_MAX
XLABEL = 'Column Steel Reinforcement Ratio [%]'   
YLABEL = 'Beam Steel Reinforcement Ratio [%]' 
ZLABEL = 'Structural Behavior Coefficient'
PLOT_3D(26, X, Y, Z, XLABEL, YLABEL, ZLABEL)  

X = RO_COL_MAX
Y = RO_BE_MAX
Z = PERIOD_01_MAX
XLABEL = 'Column Steel Reinforcement Ratio [%]'   
YLABEL = 'Beam Steel Reinforcement Ratio [%]' 
ZLABEL = 'Structural Period 01 [s]' 
PLOT_3D(27, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = RO_COL_MAX
Y = RO_BE_MAX
Z = PERIOD_02_MAX
XLABEL = 'Column Steel Reinforcement Ratio [%]'   
YLABEL = 'Beam Steel Reinforcement Ratio [%]' 
ZLABEL = 'Structural Period 02 [s]' 
PLOT_3D(28, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = RO_COL_MAX
Y = RO_BE_MAX
Z = DAMPING_RATIO_MAX
XLABEL = 'Column Steel Reinforcement Ratio [%]'   
YLABEL = 'Beam Steel Reinforcement Ratio [%]' 
ZLABEL = 'Structural Damping Ratio [%]' 
PLOT_3D(29, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = RO_COL_MAX
Y = RO_BE_MAX
Z = velocity_X_MAX_d
XLABEL = 'Column Steel Reinforcement Ratio [%]'   
YLABEL = 'Beam Steel Reinforcement Ratio [%]' 
ZLABEL = 'Velocity in X Dir. [mm/s]'
PLOT_3D(30, X, Y, Z, XLABEL, YLABEL, ZLABEL)   

X = RO_COL_MAX
Y = RO_BE_MAX
Z = velocity_Y_MAX_d
XLABEL = 'Column Steel Reinforcement Ratio [%]'   
YLABEL = 'Beam Steel Reinforcement Ratio [%]' 
ZLABEL = 'Velocity in Y Dir. [mm/s]'
PLOT_3D(31, X, Y, Z, XLABEL, YLABEL, ZLABEL)   

X = RO_COL_MAX
Y = RO_BE_MAX
Z = acceleration_X_MAX_d
XLABEL = 'Column Steel Reinforcement Ratio [%]'   
YLABEL = 'Beam Steel Reinforcement Ratio [%]' 
ZLABEL = 'Acceleration in X Dir. [mm/s^2]'
PLOT_3D(32, X, Y, Z, XLABEL, YLABEL, ZLABEL)   

X = RO_COL_MAX
Y = RO_BE_MAX
Z = acceleration_Y_MAX_d
XLABEL = 'Column Steel Reinforcement Ratio [%]'   
YLABEL = 'Beam Steel Reinforcement Ratio [%]' 
ZLABEL = 'Acceleration in Y Dir. [mm/s^2]'
PLOT_3D(33, X, Y, Z, XLABEL, YLABEL, ZLABEL)  

X = RO_COL_MAX
Y = RO_BE_MAX
Z = DRIFT_MAX
XLABEL = 'Column Steel Reinforcement Ratio [%]'   
YLABEL = 'Beam Steel Reinforcement Ratio [%]' 
ZLABEL = 'Structural Lateral Drift [%] from Pushover Analysis'
PLOT_3D(34, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = RO_COL_MAX
Y = RO_BE_MAX
Z = DRIFT_MAX_d
XLABEL = 'Column Steel Reinforcement Ratio [%]'   
YLABEL = 'Beam Steel Reinforcement Ratio [%]' 
ZLABEL = 'Structural Lateral Drift [%] from Dynamic Analysis'
PLOT_3D(35, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = RO_COL_MAX
Y = RO_BE_MAX
Z = np.divide(DRIFT_MAX_d, DRIFT_MAX)
XLABEL = 'Column Steel Reinforcement Ratio [%]'   
YLABEL = 'Beam Steel Reinforcement Ratio [%]' 
ZLABEL = 'Structural Lateral Drift [%] Ratio in X Dir.'
PLOT_3D(36, X, Y, Z, XLABEL, YLABEL, ZLABEL)
#%%------------------------------------------------------------------------------
def PLOT_LINE(title, xlabel, ylabel, x, y, color='blue', fig_num=1, logx=False, logy=False):
    import matplotlib.pyplot as plt
    plt.figure(fig_num, figsize=(12, 8))
    plt.plot(x, y, color=color, linewidth=2)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if logx:
        plt.semilogx()
    if logy:
        plt.semilogy()
    plt.grid()
    plt.show()

def PLOT_SCATTER(title, xlabel, ylabel, x, y, color='black', fig_num=1, logx=False, logy=False):
    import matplotlib.pyplot as plt
    plt.figure(fig_num, figsize=(12, 8))
    plt.scatter(x, y, color=color, linewidth=2)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if logx:
        plt.semilogx()
    if logy:
        plt.semilogy()
    plt.grid()
    plt.show()

#%%------------------------------------------------------------------------------  
# %% Plot 2D Frame Shapes
S04.PLOT_2D_FRAME(deformed_scale=1)  # Adjust scale factor as needed
#%%------------------------------------------------------------------------------   
# EXPORT DATA TO EXCEL
import pandas as pd
DATA_TOTAL = { 
    'DIAc': DIAc_MAX,
    'Kfc': Kfc_MAX,
    'RO_COL': RO_COL_MAX,
    'RO_BE': RO_BE_MAX,
    'DISP_X': DISP_X_MAX, 
    'DISP_Y': DISP_Y_MAX,
    'ROTATION': ROT_MAX,
    'AXIAL_FORCE': FORCE_A_MAX,
    'SHEAR_FORCE': FORCE_S_MAX,
    'MOMENT_WO': MOMENT_MAX,
    'AXIAL_RIGIDITY': np.abs(FORCE_A_MAX),
    'ROTATIONAL_ST': KI_MAX,
    'LATERAL_ST_Y': KA_MAX,
    'LATERAL_ST_X': KS_MAX,
    'Elastic_ST_MAX': Elastic_ST_MAX,
    'Plastic_ST_MAX': Plastic_ST_MAX,
    'Tangent_ST_MAX': Tangent_ST_MAX,
    'Ductility_Rito_MAX': Ductility_Rito_MAX,
    'Over_Strength_Factor_MAX': Over_Strength_Factor_MAX,
    'Structural_Behavior_Coefficient_MAX': R_MAX,
}
# Convert to DataFrame
results_df = pd.DataFrame(DATA_TOTAL)
# Export the DataFrame to an Excel file
results_df.to_excel('SENSITIVITY_CONFINEMENT_ENHANCEMENT_RATIO_&_REBAR_&_Cdepth_FREE_VIBRATION_RESULTS.xlsx', index=False) 
#%%------------------------------------------------------------------------------
def PLOT_HEATMAP(df):
    import matplotlib.pyplot as plt
    import numpy as np

    # Calculate the correlation matrix
    corr_matrix = df.corr()

    # Create the figure and axis
    fig, ax = plt.subplots(figsize=(12, 12))

    # Create the heatmap
    cax = ax.matshow(corr_matrix, cmap='viridis')

    # Add the colorbar
    fig.colorbar(cax)

    # Set axis labels
    ax.set_xticks(np.arange(len(corr_matrix.columns)))
    ax.set_yticks(np.arange(len(corr_matrix.index)))
    ax.set_xticklabels(corr_matrix.columns, rotation=90)
    ax.set_yticklabels(corr_matrix.index)

    # Annotate the heatmap with the correlation values
    for i in range(len(corr_matrix.columns)):
        for j in range(len(corr_matrix.index)):
            ax.text(j, i, f'{corr_matrix.iloc[i, j]:.2f}', 
                    ha='center', va='center', color='white')

    # Set title and layout
    plt.title('Correlation Heatmap', pad=20)
    plt.xlabel('Variable')
    plt.ylabel('Variable')
    
    # Display the plot
    plt.tight_layout()
    plt.show()

# PLOT HEATMAP FOR CORRELATION 
PLOT_HEATMAP(results_df)
#%%------------------------------------------------------------------------------
# Print out the state of all nodes
ops.printModel("node",1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
# Print out the state of all elements 
ops.printModel("ele", 1, 2, 3, 4, 5, 6, 7, 8, 9)
# Print the Model
#printModel()
ops.printModel("-JSON", "-file", "SENSITIVITY_CONFINEMENT_ENHANCEMENT_RATIO_&_REBAR_&_Cdepth_FREE_VIBRATION.json")
#%%------------------------------------------------------------------------------

    