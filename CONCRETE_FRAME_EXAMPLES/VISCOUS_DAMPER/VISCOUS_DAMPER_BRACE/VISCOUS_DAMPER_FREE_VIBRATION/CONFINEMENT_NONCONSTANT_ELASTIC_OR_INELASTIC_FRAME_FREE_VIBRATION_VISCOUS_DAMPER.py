######################################################################################################################
#                            >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                        #
#                   COMPARATIVE STUDY OF ELASTIC AND INELASTIC STRUCTURAL BEHAVIOR THROUGH PUSHOVER                  #
#                FREE-VIBRATION ANALYSIS USING OPENSEES ANALYSIS OF CONCRETE FRAME WITH VISCOUS DAMPER.              # 
#--------------------------------------------------------------------------------------------------------------------#
#                              THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                            #
#                                       EMAIL: salar.d.ghashghaei@gmail.com                                          #
######################################################################################################################
"""
1. The script performs a comparative study of elastic and inelastic structural behavior using pushover and dynamic analysis in OpenSeesPy.  
2. It models a 2D reinforced concrete frame with columns (500x500mm) and beams (500x300mm) using nonlinear beam-column elements.  
3. Material definitions include confined/unconfined concrete (fc=25MPa) and bilinear steel reinforcement (with/without hardening).  
4. Two analysis types are implemented: static pushover (displacement-controlled) and dynamic time-history analysis.  
5. Rayleigh damping is calculated based on initial modal periods and target damping ratios (5%).  
6. The pushover analysis applies incremental displacements up to 675mm, recording base reactions and drifts.  
7. Dynamic analysis uses Newmark integration with ground motion inputs scaled by 0.01g (El Centro records in X/Y directions).  
8. Key outputs include force-displacement curves, moment-rotation relationships, and stiffness degradation plots.  
9. Eigenvalue analysis tracks period elongation due to inelasticity during dynamic events.  
10. Damage indices are computed for ductility assessment using bilinear curve fitting.  
11. Overstrength factors (Ω), ductility ratios (μ), and R-factors are quantified for seismic performance evaluation.  
12. Real-time monitoring of base shear, axial forces, and interstory drifts is implemented.  
13. The script includes advanced convergence controls (Newton-Raphson, 1e-6 tolerance).  
14. Confinement effects are modeled with variable enhancement ratios (Kc=1.25 for columns).  
15. Results are exported to Excel, including displacements, forces, stiffness, and period data.  
16. Visualization includes deformed shapes, hysteresis loops, and cumulative response envelopes.  
17. Damping ratios are estimated from free vibration decay in dynamic analyses.  
18. Both geometric nonlinearities (P-Delta/Corotational) and material nonlinearities are considered.  
19. The code supports parametric studies by varying steel models (with/without hardening) and element types.  
20. Comprehensive plotting functions enable side-by-side comparison of elastic vs. inelastic responses for code compliance checks.  

-------------------------------------------------------------------
Viscous Dampers:

Definition:
- Devices to reduce vibrations and absorb energy.
- Resist relative displacement between points; dissipate energy as heat.

Applications:
1. Reduce seismic response in tall buildings.
2. Mitigate wind-induced vibrations in tall/light structures.
3. Control vibrations in bridges and special structures.
4. Retrofit existing buildings without major system changes.

Advantages:
- Increase effective damping with minimal stiffness change.
- Usable in new builds and retrofits.
- Long life, low maintenance.
- Perform over small and large vibration amplitudes.
- Stable performance in varying temperature/environment.

Disadvantages:
- High initial cost.
- Requires careful placement and tuning.
- Possible fluid leakage; periodic inspection needed.
- Displacement-dependent effectiveness.

-------------------------------------------------------------------
Force–Displacement Relationship:

Linear viscous damping:
    F = C * v
    C = damping coefficient (N·s/m), v = relative velocity.
    Produces elliptical hysteresis loop (force ~ velocity).

Nonlinear viscous damping:
    F = C * |v|^(α - 1) * v
    α ≈ 0.3–1.0 for broad performance.

Ideal Properties for Structural Engineers:
- Effective at low/high amplitudes.
- No unwanted stiffness.
- Consistent hysteresis over cycles.
- Temperature/environment stability.
- Symmetric tension/compression response.
- Force within member capacity.

Design Note:
- Ideal force–displacement curve: closed ellipse.
- Area = energy dissipated per cycle.
- Select C, α to match energy needs, limit peak force, and target fundamental period.  
"""
#%%------------------------------------------------------------------------------

# YOUTUBE: What is a Fluid Viscous Damper
'https://www.youtube.com/watch?v=OoWtKmNoz6Q'    
# YOUTUBE: Viscous Fluid Damper (abbreviated as VFD or FVD)
'https://www.youtube.com/watch?v=wLLzeJtpWq0'
# YOUTUBE: Fluid Viscous Damping: Uses, Benefits and Methods of Implementation
'https://www.youtube.com/watch?v=q0dJNhB8AKU'
# YOUTUBE: Superior Steel & Ironworkers installation of Viscous Fluid Dampers
'https://www.youtube.com/watch?v=Yh0Uik_WscE'    
# PAPER: A comprehensive study of viscous damper configurations and vertical damping coefficient distributions for enhanced performance in reinforced concrete structures
'https://link.springer.com/article/10.1007/s42107-023-00831-x'
# PAPER: Seismic Upgrade of Existing Buildings with Fluid Viscous Dampers: Design Methodologies and Case Study
'https://ascelibrary.com/doi/10.1061/%28ASCE%29CF.1943-5509.0000671'    
# PAPER: Seismic damage prediction by deterministic methods: Concepts and procedures
'https://onlinelibrary.wiley.com/doi/10.1002/eqe.4290160507'
# PAPER: A strain-consistent approach for determination of bounds of ductility damage index for different performance levels for seismic design of RC frame members
'https://www.sciencedirect.com/science/article/abs/pii/S0141029611005086'
# PAPER: Energy Index Method: Technique for Identification of Structural Damages
'https://ascelibrary.org/doi/10.1061/%28ASCE%290733-9445%282008%29134%3A6%281061%29'
#%%------------------------------------------------------------------------------
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
import EIGENVALUE_ANALYSIS_FUN as S07
import TRUSS_ELEMENT_FUN as S08
import HARMONIC_IMPACT_LOAD_FUN as S09

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
GMfact = 9810    # [mm/s²]standard acceleration of gravity or standard acceleration
FV = 'X'         # Free-Vibration Direction
DR = 0.05        # Intial Guess for Damping ratio
duration = 40.0  # [s] Total simulation duration
dt = 0.01        # [s] Time step
MASS = 12000     # [kg] Mass on the each column

#%% DEFINE PARAMETERS FOR NONLINEAR STATIC ANALYSIS 
DMAX = 675       # [mm] Maximum Displacement
DINCR = 0.05     # [mm] Incremental Displacement

#%% DEFINE DAMPERPARAMETERS
FY = 3650        # [N/mm^2] Yield Strength Viscous Damper
DRD = 0.15       # Intial Guess for Damping ratio for Viscous Damper
area = 35**2     # [mm^2] Damper Secton Area

#%% DEFINE ANALYSIS PROPERTIES
MAX_ITERATIONS = 5000    # Convergence iteration for test
MAX_TOLERANCE = 1.0e-6   # Convergence tolerance for test
#STEEL_KIND: 1 -> WITHOUT HARDENING AND ULTIMATE STRAIN
#STEEL_KIND: 2 -> WITH HARDENING AND ULTIMATE STRAIN
#%%------------------------------------------------------------------------------
# IMPACT ANALYSIS FUNCTION
def PD_ANALYSIS(u0, v0, a0, IA, IU, IV, STEEL_KIND, ANA_KIND, ELE_KIND):
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
    Kc = 1.25        # Confinement Enhancement Ratio
    S03.CONFINED_CONCRETE_SECTION(secTagC01, Hc, Bc, coverC, AsC, STEEL_KIND, fc, Kc, COL=True)
    # COLUMN SECTION FOR ELEMENT 02, 05
    S03.CONFINED_CONCRETE_SECTION(secTagC02, Hc, Bc, coverC, AsC, STEEL_KIND, fc, Kc-0.15, COL=True)
    # BEAM SECTION FOR ELEMENT 07, 09
    S03.CONFINED_CONCRETE_SECTION(secTagB01, Hb, Bb, coverB, AsB, STEEL_KIND, fc, Kc-0.1, COL=False)
    # BEAM SECTION FOR ELEMENT 08
    S03.CONFINED_CONCRETE_SECTION(secTagB02, Hb, Bb, coverB, AsB, STEEL_KIND, fc, Kc-0.2, COL=False)
    # VISCOUS DAMPER
    ele_tag = 500
    ELE_KIND = 'ELASTIC'
    NODE_I = 1
    NODE_J = 4
    # Tension Force-Displacement Relationship (Positive values)
    FY1, DY1 = FY, 0.01             # First yield point
    FY2, DY2 = 1.121*FY, 0.02       # Peak force
    FY3, DY3 = 0.800*FY, 0.04       # Post-peak softening
    FY4, DY4 = 0.660*FY, 0.06       # Plateau
    FY5, DY5 = 0.430*FY, 0.28       # Further softening
    FY6, DY6 = 0.210*FY, 0.51       # Near-zero stiffness
    FY7, DY7 = 0.100*FY, 0.72       # Zero force (failure)
    
    KP = np.array([FY1, DY1, FY2, DY2, FY3, DY3, FY4, DY4, FY5, DY5, FY6, DY6, FY7, DY7])

    # Compression Force-Displacement Relationship (Negative values)
    FY1n, DY1n = -1.121*FY, -0.01    # First yield in compression
    FY2n, DY2n = -0.800*FY, -0.05    # Peak compressive force
    FY3n, DY3n = -0.430*FY, -0.20    # Post-peak softening

    KN = np.array([FY1n, DY1n, FY2n, DY2n, FY3n, DY3n])
    
    # VISCOUS DAMPER ELEMENT
    S08.TRUSS_ELEMENT(ele_tag, NODE_I, NODE_J, KP, KN, area, MASS, ELE_KIND, DRD, PLOT=True)
    
    # Define geometric transformation (corotational for large displacements)
    transfTag = 1
    #ops.geomTransf('Linear', transfTag)
    #ops.geomTransf('PDelta', transfTag)
    ops.geomTransf('Corotational', transfTag)
    
    if ELE_KIND == 'INELASTIC':
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
        
    if ELE_KIND == 'ELASTIC':
        # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/elasticBeamColumn.html
        # LEFT COLUMN
        ops.element('elasticBeamColumn', 1, 1, 5, secTagC01, transfTag, '-mass', Bc*Hc*0.000025) # COLUMN - ELEMENT 01
        ops.element('elasticBeamColumn', 2, 5, 6, secTagC02, transfTag, '-mass', Bc*Hc*0.000025) # COLUMN - ELEMENT 02
        ops.element('elasticBeamColumn', 3, 6, 3, secTagC01, transfTag, '-mass', Bc*Hc*0.000025) # COLUMN - ELEMENT 03
        # RIGHT COLUMN
        ops.element('elasticBeamColumn', 4, 2, 7, secTagC01, transfTag, '-mass', Bc*Hc*0.000025) # COLUMN - ELEMENT 04
        ops.element('elasticBeamColumn', 5, 7, 8, secTagC02, transfTag, '-mass', Bc*Hc*0.000025) # COLUMN - ELEMENT 05
        ops.element('elasticBeamColumn', 6, 8, 4, secTagC01, transfTag, '-mass', Bc*Hc*0.000025) # COLUMN - ELEMENT 06
        # BEAM
        ops.element('elasticBeamColumn', 7, 3, 9, secTagB01, transfTag, '-mass', Bb*Hb*0.000025)  # BEAM - ELEMENT 07
        ops.element('elasticBeamColumn', 8, 9, 10, secTagB02, transfTag, '-mass', Bb*Hb*0.000025) # BEAM - ELEMENT 08
        ops.element('elasticBeamColumn', 9, 10, 4, secTagB01, transfTag, '-mass', Bb*Hb*0.000025) # BEAM - ELEMENT 09

    if ANA_KIND == 'PUSHOVER':
        WEIGHT = MASS * GMfact
        # Data storage
        FORCE_S, FORCE_A, MOMENT = [], [], []
        DISP_X, DISP_Y, ROT = [], [], []
        KA, KS, KI, STEP = [], [], [], []
        DRIFT_X, DRIFT_Y = [], []
    
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
            KS.append(np.abs(S)/np.abs(disp_X))   # LATERAL STIFFNESS IN X
            KA.append(np.abs(A)/np.abs(disp_Y))   # LATERAL STIFFNESS IN Y
            KI.append(np.abs(M)/np.abs(rot))      # ROTATIONAL STIFFNESS IN Z
            STEP.append(step)
            DRIFT_X.append(100*disp_X/LENGTH_COL) # LATERAL STRUCTURE DRIFT IN X
            DRIFT_Y.append(100*disp_Y/LENGTH_BM)  # LATERAL STRUCTURE DRIFT IN Y
            print(step+1, disp_X, S)
            
        return FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y, ROT, KA, KS, KI, DRIFT_X, DRIFT_Y, STEP        
    
    if ANA_KIND == 'DYNAMIC':
        # Define mass
        ops.mass(3, MASS, MASS, 0.0)
        ops.mass(4, MASS, MASS, 0.0)
        
        # Define time series for input motion (Acceleration time history)
        time_series_tag = 1
        pattern_tag = 1
        
        # Apply Loads
        ops.timeSeries('Linear', time_series_tag)
        ops.pattern('Plain', pattern_tag, time_series_tag)
        print('Impact Load Defined Done.')
        if FV == 'X':
            ops.load(3, 1.0, 0.0, 0.0)
            ops.load(4, 1.0, 0.0, 0.0)
        if FV == 'Y': 
            ops.load(3, 0.0, -1.0, 0.0)
            ops.load(4, 0.0, -1.0, 0.0)
        if FV == 'XY':
            ops.load(3, 1.0, -1.0, 0.0)
            ops.load(4, 1.0, -1.0, 0.0)
            
        if IU == True:
            # Define initial displacment
            ops.setNodeDisp(3, 1, u0, '-commit') # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/setNodeDisp.html
        if IV == True:
            # Define initial velocity
            ops.setNodeVel(3, 1, v0, '-commit')  # INFO LINK: https://openseespydoc.readthedocs.io/en/stable/src/setNodeVel.html
        if IA == True:
            # Define initial  acceleration
            ops.setNodeAccel(3, 1, a0, '-commit') # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/setNodeAccel.html
        #print('Free-vibration Defined.')            
        
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
        DRIFT_X, DRIFT_Y = [], []
        PERIOD_MIN, PERIOD_MAX = [], []
            
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
            KS.append(np.abs(S/disp_X))           # LATERAL STIFFNESS IN X
            KA.append(np.abs(A/disp_Y))           # LATERAL STIFFNESS IN Y
            KI.append(np.abs(M/rot))              # ROTATIONAL STIFFNESS IN Z
            DRIFT_X.append(100*disp_X/LENGTH_COL) # LATERAL STRUCTURE DRIFT IN X
            DRIFT_Y.append(100*disp_Y/LENGTH_BM)  # LATERAL STRUCTURE DRIFT IN Y
            # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
            Period_MIN, Period_MAX = S07.EIGENVALUE_ANALYSIS(4, PLOT=True)
            PERIOD_MIN.append(Period_MIN) 
            PERIOD_MAX.append(Period_MAX)
            print(current_time, disp_X, S)
        
        # Calculate Structure Damping Ratio Based on Lateral Displacement
        damping_ratio = S05.DAMPING_RATIO(DISP_X)    
        
        #ops.wipe()  
        return FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y, ROT, KA, KS, KI, time, velocity_X, velocity_Y, acceleration_X, acceleration_Y, PERIOD_01, PERIOD_02, DRIFT_X, DRIFT_Y, damping_ratio, PERIOD_MIN, PERIOD_MAX

#%%------------------------------------------------------------------------------
# FREE-VIBRATION ANALYSIS PARAMETERS:
IU = True         # Free Vibration with Initial Displacement
IV = True         # Free Vibration with Initial Velocity
IA = True         # Free Vibration with Initial Acceleration
u0 = 105.35       # [mm] Initial displacement
v0 = 0.015        # [mm/s] Initial velocity
a0 = 0.0065       # [mm/s^2] Initial acceleration
#%%------------------------------------------------------------------------------
# Analysis Durations for Linear or Nonlinear Static Analysis:
starttime = TI.process_time()

# WITHOUT HARDENING AND ULTIMATE STRAIN -> STEEL_KIND=1
# WITH HARDENING AND ULTIMATE STRAIN -> STEEL_KIND=2

# RUN ELASTIC STATIC ANALYSIS
STEEL_KIND = 2
ANA_KIND = 'PUSHOVER' # 'PUSHOVER' OR 'DYNAMIC' ANALYSIS
ELE_KIND = 'ELASTIC'  # 'ELASTIC' OR 'INELASTIC' ELEMENTS
DATA = PD_ANALYSIS(u0, v0, a0, IA, IU, IV, STEEL_KIND, ANA_KIND, ELE_KIND)
FORCE_Spe, FORCE_Ape, MOMENTpe, DISP_Xpe, DISP_Ype, ROTpe, KApe, KSpe, KIpe, DRIFT_Xpe, DRIFT_Ype, STEPpe = DATA

# Plot 2D Frame Shapes for Elastic Static Analysis
S04.PLOT_2D_FRAME(deformed_scale=1)  # Adjust scale factor as needed
#%%------------------------------------------------------------------------------
# RUN INELASTIC STATIC ANALYSIS
STEEL_KIND = 2
ANA_KIND = 'PUSHOVER' # 'PUSHOVER' OR 'DYNAMIC' ANALYSIS
ELE_KIND = 'INELASTIC'  # 'ELASTIC' OR 'INELASTIC' ELEMENTS
DATA = PD_ANALYSIS(u0, v0, a0, IA, IU, IV, STEEL_KIND, ANA_KIND, ELE_KIND)
FORCE_Spp, FORCE_App, MOMENTpp, DISP_Xpp, DISP_Ypp, ROTpp, KApp, KSpp, KIpp, DRIFT_Xpp, DRIFT_Ypp, STEPpp = DATA

# Plot 2D Frame Shapes for Inelastic Static Analysis
S04.PLOT_2D_FRAME(deformed_scale=1)  # Adjust scale factor as needed
#%%------------------------------------------------------------------------------
totaltime = TI.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n') 
#%%------------------------------------------------------------------------------
# Analysis Durations for Linear or Nonlinear Dynamic Analysis:
starttime = TI.process_time()

# RUN ELASTIC DYNAMIC ANALYSIS
STEEL_KIND = 2
ANA_KIND = 'DYNAMIC'  # 'PUSHOVER' OR 'DYNAMIC' ANALYSIS
ELE_KIND = 'ELASTIC'  # 'ELASTIC' OR 'INELASTIC' ELEMENTS
DATA = PD_ANALYSIS(u0, v0, a0, IA, IU, IV, STEEL_KIND, ANA_KIND, ELE_KIND)
FORCE_Sde, FORCE_Ade, MOMENTde, DISP_Xde, DISP_Yde, ROTde, KAde, KSde, KIde, timede, velocity_Xde, velocity_Yde, acceleration_Xde, acceleration_Yde, PERIOD_01de, PERIOD_02de, DRIFT_Xde, DRIFT_Yde, damping_ratioe, PERIOD_MINe, PERIOD_MAXe = DATA

print(f"\n ELASTIC Period 01: {PERIOD_01de:.4e} (s) - Period 02: {PERIOD_02de:.4e} (s)")

# Plot 2D Frame Shapes for Elastic Dynamic Analysis
S04.PLOT_2D_FRAME(deformed_scale=10)  # Adjust scale factor as needed
#%%-----------------------------------------------------------------------------
# RUN INELASTIC DYNAMIC ANALYSIS
STEEL_KIND = 2
ANA_KIND = 'DYNAMIC'  # 'PUSHOVER' OR 'DYNAMIC' ANALYSIS
ELE_KIND = 'INELASTIC'  # 'ELASTIC' OR 'INELASTIC' ELEMENTS
DATA = PD_ANALYSIS(u0, v0, a0, IA, IU, IV, STEEL_KIND, ANA_KIND, ELE_KIND)
FORCE_Sdp, FORCE_Adp, MOMENTdp, DISP_Xdp, DISP_Ydp, ROTdp, KAdp, KSdp, KIdp, timedp, velocity_Xdp, velocity_Ydp, acceleration_Xdp, acceleration_Ydp, PERIOD_01dp, PERIOD_02dp, DRIFT_Xdp, DRIFT_Ydp, damping_ratiop, PERIOD_MINp, PERIOD_MAXp = DATA

print(f"\n INELASTIC Period 01: {PERIOD_01dp:.4e} (s) - Period 02: {PERIOD_02dp:.4e} (s)")

totaltime = TI.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n') 

# Plot 2D Frame Shapes for Inelastic Dynamic Analysis
S04.PLOT_2D_FRAME(deformed_scale=10)  # Adjust scale factor as needed
#%%-----------------------------------------------------------------------------
plt.figure(0, figsize=(12, 8))
plt.plot(timede, PERIOD_MINe)
plt.plot(timede, PERIOD_MAXe)
plt.plot(timedp, PERIOD_MINp)
plt.plot(timedp, PERIOD_MAXp)
plt.title('Period of Structure')
plt.ylabel('Structural Period [s]')
plt.xlabel('Time [s]')
#plt.semilogy()
plt.grid()
plt.legend([f'ELASTIC PERIOD - MIN VALUES: Min: {np.min(PERIOD_MINe):.3f} (s) - Mean: {np.mean(PERIOD_MINe):.3f} (s) - Max: {np.max(PERIOD_MINe):.3f} (s)', 
            f'ELASTIC PERIOD - MAX VALUES:  Min: {np.min(PERIOD_MAXp):.3f} (s) - Mean: {np.mean(PERIOD_MAXp):.3f} (s) - Max: {np.max(PERIOD_MAXp):.3f} (s)',
            f'INELASTIC PERIOD - MIN VALUES: Min: {np.min(PERIOD_MINp):.3f} (s) - Mean: {np.mean(PERIOD_MINp):.3f} (s) - Max: {np.max(PERIOD_MINp):.3f} (s)', 
            f'INELASTIC PERIOD - MAX VALUES:  Min: {np.min(PERIOD_MAXp):.3f} (s) - Mean: {np.mean(PERIOD_MAXp):.3f} (s) - Max: {np.max(PERIOD_MAXp):.3f} (s)',
            ])
plt.show()

plt.figure(1, figsize=(12, 8))
plt.plot(MOMENTpe, FORCE_Ape)
plt.plot(MOMENTpp, FORCE_App)
plt.plot(MOMENTde, FORCE_Ade)
plt.plot(MOMENTdp, FORCE_Adp)
plt.title('P-M Interaction')
plt.ylabel('Axial Force [N]')
plt.xlabel('Bending Moment [N.mm]')
plt.legend(['ELASTIC-PUSHOVER', 'INELASTIC-PUSHOVER', 'ELASTIC-DYNAMIC', 'INELASTIC-DYNAMIC'])
plt.grid()
plt.show()

plt.figure(2, figsize=(12, 8))
plt.plot(DISP_Xpe, FORCE_Spe,  linewidth=2)
plt.plot(DISP_Xpp, FORCE_Spp,  linewidth=2)
plt.plot(DISP_Xde, FORCE_Sde,  linewidth=2)
plt.plot(DISP_Xdp, FORCE_Sdp,  linewidth=2)

plt.title('SHEAR FORCE-DISPLACEMENT DIAGRAM')
plt.ylabel('Shear Force [N]')
plt.xlabel('Displacement  in X [mm]')
plt.legend(['ELASTIC-PUSHOVER', 'INELASTIC-PUSHOVER', 'ELASTIC-DYNAMIC', 'INELASTIC-DYNAMIC'])
plt.grid()
plt.show()

plt.figure(3, figsize=(12, 8))
plt.plot(DISP_Ype, FORCE_Ape, linewidth=2)
plt.plot(DISP_Ypp, FORCE_App, linewidth=2)
plt.plot(DISP_Yde, FORCE_Ade, linewidth=2)
plt.plot(DISP_Ydp, FORCE_Adp, linewidth=2)
plt.title('AXIAL FORCE-DISPLACEMENT DIAGRAM')
plt.ylabel('Axial Force [N]')
plt.xlabel('Displacement in Y [mm]')
plt.legend(['ELASTIC-PUSHOVER', 'INELASTIC-PUSHOVER', 'ELASTIC-DYNAMIC', 'INELASTIC-DYNAMIC'])
plt.grid()
plt.show()

plt.figure(4, figsize=(12, 8))
plt.plot(ROTpe, MOMENTpe, linewidth=2)
plt.plot(ROTpp, MOMENTpp, linewidth=2)
plt.plot(ROTde, MOMENTde, linewidth=2)
plt.plot(ROTdp, MOMENTdp, linewidth=2)
plt.title('MOMENT-ROTATION DIAGRAM')
plt.ylabel('Moment [kN.mm]')
plt.xlabel('Rotation [rad]')
plt.legend(['ELASTIC-PUSHOVER', 'INELASTIC-PUSHOVER', 'ELASTIC-DYNAMIC', 'INELASTIC-DYNAMIC'])
plt.grid()
plt.show()

plt.figure(5, figsize=(12, 8))
plt.scatter(KIpe, KSpe, linewidth=2)
plt.scatter(KIpp, KSpp, linewidth=2)
plt.scatter(KIde, KSde,  linewidth=2)
plt.scatter(KIdp, KSdp,  linewidth=2)
plt.title('ROTATIONAL STIFFNESS-LATERAL STIFFNESS DIAGRAM')
plt.ylabel('Rotational Stiffness [N.mm/Rad]')
plt.xlabel('Lateral Stiffness in X Dir. [N/mm]')
plt.semilogx()
plt.semilogy()
plt.legend(['ELASTIC-PUSHOVER', 'INELASTIC-PUSHOVER', 'ELASTIC-DYNAMIC', 'INELASTIC-DYNAMIC'])
plt.grid()
plt.show()

plt.figure(6, figsize=(12, 8))
plt.scatter(KIpe, KApe, linewidth=2)
plt.scatter(KIpp, KApp, linewidth=2)
plt.scatter(KIde, KAde, linewidth=2)
plt.scatter(KIdp, KAdp, linewidth=2)
plt.title('ROTATIONAL STIFFNESS-LATERAL STIFFNESS DIAGRAM')
plt.ylabel('Rotational Stiffness [N.mm/Rad]')
plt.xlabel('Lateral Stiffness in Y Dir. [N/mm]')
plt.semilogx()
plt.semilogy()
plt.legend(['ELASTIC-PUSHOVER', 'INELASTIC-PUSHOVER', 'ELASTIC-DYNAMIC', 'INELASTIC-DYNAMIC'])
plt.grid()
plt.show()

plt.figure(7, figsize=(12, 8))
plt.plot(timede, FORCE_Ade, linewidth=2)
plt.plot(timedp, FORCE_Adp, linewidth=2)
plt.title('Axial Force During the Analysis')
plt.ylabel('Axial Force [N]')
plt.xlabel('Times')
plt.legend(['ELASTIC', 'INELASTIC'])
plt.grid()
plt.show()

plt.figure(8, figsize=(12, 8))
plt.plot(timede, FORCE_Sde, linewidth=2)
plt.plot(timedp, FORCE_Sdp, linewidth=2)
plt.title('Shear Force During the Analysis')
plt.ylabel('Shear Force [N]')
plt.xlabel('Times')
plt.legend(['ELASTIC', 'INELASTIC'])
plt.grid()
plt.show()

plt.figure(9, figsize=(12, 8))
plt.plot(timede, MOMENTde, linewidth=2)
plt.plot(timedp, MOMENTdp, linewidth=2)
plt.title('Moment During the Analysis')
plt.ylabel('Moment [N.mm]')
plt.xlabel('Times')
plt.legend(['ELASTIC', 'INELASTIC'])
plt.grid()
plt.show()

plt.figure(10, figsize=(12, 8))
plt.plot(timede, DISP_Xde, linewidth=2)
plt.plot(timedp, DISP_Xdp, linewidth=2)
plt.title('Displacement During the Analysis')
plt.ylabel('Displacement - X [mm]')
plt.xlabel('Times')
plt.legend(['ELASTIC', 'INELASTIC'])
plt.grid()
plt.show()

plt.figure(11, figsize=(12, 8))
plt.plot(timede, DISP_Yde, linewidth=2)
plt.plot(timedp, DISP_Ydp, linewidth=2)
plt.title('Displacement During the Analysis')
plt.ylabel('Displacement - Y [mm]')
plt.xlabel('Times')
plt.legend(['ELASTIC', 'INELASTIC'])
plt.grid()
plt.show()

plt.figure(12, figsize=(12, 8))
plt.plot(timede, ROTde, linewidth=2)
plt.plot(timedp, ROTdp, linewidth=2)
plt.title('Rotation During the Analysis')
plt.ylabel('Rotation [rad]')
plt.xlabel('Times')
plt.legend(['ELASTIC', 'INELASTIC'])
plt.grid()
plt.show()

plt.figure(13, figsize=(12, 8))
plt.plot(STEPpe, DRIFT_Xpe, linewidth=2)
plt.plot(STEPpp, DRIFT_Xpp, linewidth=2)
plt.plot(STEPpe, DRIFT_Ype, linewidth=2)
plt.plot(STEPpp, DRIFT_Ypp, linewidth=2)
plt.title('Structure Lateral Drift During the Pushover Analysis')
plt.ylabel('Lateral Drift [%]')
plt.xlabel('Steps')
plt.legend([f'PUSHOVER IN X DIR.: {np.max(np.abs(DRIFT_Xpe)):.3e}',
            f'PUSHOVER IN X DIR.: {np.max(np.abs(DRIFT_Xpp)):.3e}',
            f'PUSHOVER IN Y DIR.: {np.max(np.abs(DRIFT_Ype)):.3e}',
            f'PUSHOVER IN Y DIR.: {np.max(np.abs(DRIFT_Ypp)):.3e}'])
plt.grid()
plt.show()

plt.figure(14, figsize=(12, 8))
plt.plot(timede, DRIFT_Xde, linewidth=2)
plt.plot(timedp, DRIFT_Xdp, linewidth=2)
plt.plot(timede, DRIFT_Yde, linewidth=2)
plt.plot(timedp, DRIFT_Ydp, linewidth=2)
plt.title('Structure Lateral Drift During the Dynamic Analysis')
plt.ylabel('Lateral Drift [%]')
plt.xlabel('Times')
plt.legend([f'DYN. IN X DIR.: {np.max(np.abs(DRIFT_Xde)):.3e}',
            f'DYN. IN X DIR.: {np.max(np.abs(DRIFT_Xdp)):.3e}',
            f'DYN. IN Y DIR.: {np.max(np.abs(DRIFT_Yde)):.3e}',
            f'DYN. IN Y DIR.: {np.max(np.abs(DRIFT_Ydp)):.3e}'])
plt.grid()
plt.show()
#%%------------------------------------------------------------------------------  
# Compute the Cumulative Maximum Absolute Value of Last Analysis Data
def MAX_ABS(X):
    import numpy as np
    X = np.asarray(X)  # Convert input to a numpy array for faster operations
    X_MAX = np.zeros_like(X)  # Initialize an array to store cumulative max values
    X_MAX[0] = np.abs(X[0])  # Set the first value

    # Compute cumulative maximum absolute values
    for i in range(1, len(X)):
        X_MAX[i] = max(X_MAX[i-1], np.abs(X[i]))
    
    return X_MAX  

DISP_ZXe = MAX_ABS(DISP_Xde)  
DISP_ZXp = MAX_ABS(DISP_Xdp)  
DISP_ZYe = MAX_ABS(DISP_Yde) 
DISP_ZYp = MAX_ABS(DISP_Ydp) 
VELO_Ze = MAX_ABS(velocity_Xde) 
VELO_Zp = MAX_ABS(velocity_Xdp) 
ACCE_Ze = MAX_ABS(acceleration_Xde) 
ACCE_Zp = MAX_ABS(acceleration_Xdp) 
BASE_Ze = MAX_ABS(FORCE_Sde) 
BASE_Zp = MAX_ABS(FORCE_Sdp)
BASE_ZYe = MAX_ABS(FORCE_Ade)
BASE_ZYp = MAX_ABS(FORCE_Adp)

plt.figure(1, figsize=(8, 6))
plt.plot(timede, DISP_Xde, color='blue', linewidth=2)
plt.plot(timede, DISP_ZXe, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Displacement in X [mm]')
plt.title(f'Time vs Displacement - MAX. ABS: {DISP_ZXe[-1]} | ξ (Calculated): {damping_ratioe:.5e} %')
plt.grid()
plt.show()

plt.figure(11, figsize=(8, 6))
plt.plot(timedp, DISP_Xdp, color='blue', linewidth=2)
plt.plot(timedp, DISP_ZXp, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Displacement in X [mm]')
plt.title(f'Time vs Displacement - MAX. ABS: {DISP_ZXp[-1]} | ξ (Calculated): {damping_ratiop:.5e} %')
plt.grid()
plt.show()

plt.figure(2, figsize=(8, 6))
plt.plot(timede, DISP_Yde, color='blue', linewidth=2)
plt.plot(timede, DISP_ZYe, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Displacement in Y [mm]')
plt.title(f'Time vs Displacement - MAX. ABS: {DISP_ZYe[-1]}')
plt.grid()
plt.show()

plt.figure(22, figsize=(8, 6))
plt.plot(timedp, DISP_Ydp, color='blue', linewidth=2)
plt.plot(timedp, DISP_ZYp, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Displacement in Y [mm]')
plt.title(f'Time vs Displacement - MAX. ABS: {DISP_ZYp[-1]}')
plt.grid()
plt.show()

plt.figure(3, figsize=(8, 6))
plt.plot(timede, velocity_Xde, color='blue', linewidth=2)
plt.plot(timede, VELO_Ze, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Velocity in X [mm/s]')
plt.title(f'Time vs Velocity - MAX. ABS: {VELO_Ze[-1]}')
plt.grid()
plt.show()

plt.figure(33, figsize=(8, 6))
plt.plot(timedp, velocity_Xdp, color='blue', linewidth=2)
plt.plot(timedp, VELO_Zp, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Velocity in X [mm/s]')
plt.title(f'Time vs Velocity - MAX. ABS: {VELO_Zp[-1]}')
plt.grid()
plt.show()

plt.figure(4, figsize=(8, 6))
plt.plot(timede, acceleration_Xde, color='blue', linewidth=2)
plt.plot(timede, ACCE_Ze, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Acceleration in X [mm/s^2]')
plt.title(f'Time vs Acceleration - MAX. ABS: {ACCE_Ze[-1]}')
plt.grid()
plt.show()

plt.figure(44, figsize=(8, 6))
plt.plot(timedp, acceleration_Xdp, color='blue', linewidth=2)
plt.plot(timedp, ACCE_Zp, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Acceleration in X [mm/s^2]')
plt.title(f'Time vs Acceleration - MAX. ABS: {ACCE_Zp[-1]}')
plt.grid()
plt.show()

plt.figure(5, figsize=(8, 6))
plt.plot(timede, FORCE_Sde, color='blue', linewidth=2)
plt.plot(timede, BASE_Ze, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Shear Base-reaction [N]')
plt.title(f'Time vs Shear Base-reaction from Inelastic Dynamic Analysis- MAX. ABS: {BASE_Ze[-1]}')
plt.grid()
plt.show()

plt.figure(55, figsize=(8, 6))
plt.plot(timedp, FORCE_Sdp, color='blue', linewidth=2)
plt.plot(timedp, BASE_Zp, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Shear Base-reaction [N]')
plt.title(f'Time vs Shear Base-reaction from Inelastic Dynamic Analysis- MAX. ABS: {BASE_Zp[-1]}')
plt.grid()
plt.show()

plt.figure(6, figsize=(8, 6))
plt.plot(timede, FORCE_Ade, color='blue', linewidth=2)
plt.plot(timede, BASE_ZYe, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Axial Base-reaction [N]')
plt.title(f'Time vs Axial Base-reaction - MAX. ABS: {BASE_ZYe[-1]}')
plt.grid()
plt.show()

plt.figure(66, figsize=(8, 6))
plt.plot(timedp, FORCE_Adp, color='blue', linewidth=2)
plt.plot(timedp, BASE_ZYp, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Axial Base-reaction [N]')
plt.title(f'Time vs Axial Base-reaction - MAX. ABS: {BASE_ZYp[-1]}')
plt.grid()
plt.show()

plt.figure(7, figsize=(8, 6))
plt.plot(DISP_Xde, FORCE_Sde, color='black', linewidth=2)
plt.plot(DISP_Xdp, FORCE_Sdp, color='purple', linewidth=2)
plt.xlabel('Displacement in X [mm]')
plt.ylabel('Shear Base-reaction [N]')
plt.title('Displacement vs Shear Base-reaction')
plt.legend(['ELASTIC', 'INELASTIC'])
plt.grid()
plt.show()

plt.figure(8, figsize=(8, 6))
plt.plot(DISP_Yde, FORCE_Ade, color='brown', linewidth=2)
plt.plot(DISP_Ydp, FORCE_Adp, color='green', linewidth=2)
plt.xlabel('Displacement in Y [mm]')
plt.ylabel('Axial Base-reaction [N]')
plt.title('Displacement vs Axial Base-reaction')
plt.legend(['ELASTIC', 'INELASTIC'])
plt.grid()
plt.show()

plt.figure(9, figsize=(8, 6))
plt.plot(ROTde, MOMENTde, color='red', linewidth=2)
plt.plot(ROTdp, MOMENTdp, color='black', linewidth=2)
plt.xlabel('Rotation [rad]')
plt.ylabel('Moment Base-reaction [N.mm/Rad]')
plt.title('Rotation vs Moment Base-reaction')
plt.legend(['ELASTIC', 'INELASTIC'])
plt.grid()
plt.show()

#%%------------------------------------------------------------------------------

# --------------------------------------
#  Plot BaseShear-Displacement Analysis 
# --------------------------------------
XX = np.abs(DISP_Xpp); YY = np.abs(FORCE_Spp); # ABSOLUTE VALUE
SLOPE_NODE = 10

DATA = BC.BILNEAR_CURVE(XX, YY, SLOPE_NODE)
X, Y, Elastic_ST, Plastic_ST, Tangent_ST, Ductility_Rito, Over_Strength_Factor = DATA

XLABEL = 'Displacement in X [mm]'
YLABEL = 'Base-Shear Reaction [N]'
LEGEND01 = 'Curve'
LEGEND02 = 'Bilinear Fitted'
LEGEND03 = 'Undefined'
TITLE = f'Last Data of BaseShear-Displacement Analysis - Ductility Ratio: {X[2]/X[1]:.4f} - Over Strength Factor: {Y[2]/Y[1]:.4f}'
COLOR = 'black'
BC.PLOT_2D(np.abs(DISP_Xpp), np.abs(FORCE_Spp), X, Y, X, Y, XLABEL, YLABEL, TITLE, LEGEND01, LEGEND02, LEGEND03, COLOR='black', Z=2) 
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
Dd = np.max(np.abs(DISP_Xdp))
DIx = (Dd - X[1]) /(X[2] - X[1])
print(f'Structural Ductility Damage Index in X Direction: {DIx:.4f}')
# ---------------------------------------
#  Plot BaseAxial-Displacement Analysis
# ---------------------------------------
XX = np.abs(DISP_Ypp); YY = np.abs(FORCE_App); # ABSOLUTE VALUE
SLOPE_NODE = 10

DATA = BC.BILNEAR_CURVE(XX, YY, SLOPE_NODE)
X, Y, Elastic_ST, Plastic_ST, Tangent_ST, Ductility_Rito, Over_Strength_Factor = DATA

XLABEL = 'Displacement in Y [mm]'
YLABEL = 'Base-Axial Reaction [N]'
LEGEND01 = 'Curve'
LEGEND02 = 'Bilinear Fitted'
LEGEND03 = 'Undefined'
TITLE = f'Last Data of BaseAxial-Displacement Analysis - Ductility Ratio: {X[2]/X[1]:.4f} - Over Strength Factor: {Y[2]/Y[1]:.4f}'
COLOR = 'black'
BC.PLOT_2D(np.abs(DISP_Ypp), np.abs(FORCE_App), X, Y, X, Y, XLABEL, YLABEL, TITLE, LEGEND01, LEGEND02, LEGEND03, COLOR='black', Z=2) 
#print(f'\t\t Ductility Ratio: {YY[2]/YY[1]:.4f}')

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
Dd = np.max(np.abs(DISP_Ydp))
DIy = (Dd - X[1]) /(X[2] - X[1])
print(f'Structural Ductility Damage Index in Y Direction: {DIy:.4f}')
#%%------------------------------------------------------------------------------   
# EXPORT DATA TO EXCEL
import pandas as pd
DATA_TOTAL = { 
    'DISP_X_E': DISP_Xde,
    'DISP_Y_E': DISP_Yde,
    'ROTATION_E': ROTde,
    'VELO_E': velocity_Xde,
    'ACCEL_E': acceleration_Xde,
    'AXIAL_FORCE_E': FORCE_Ade,
    'SHEAR_FORCE_E': FORCE_Sde,
    'MOMENT_E': MOMENTde,
    'AXIAL_RIGIDITY_E': np.abs(FORCE_Ade),
    'ROTATIONAL_ST_E': KIde,
    'LATERAL_ST_Y_E': KAde,
    'LATERAL_ST_X_E': KSde,
    'PERIOD_MIN_E': PERIOD_MINe,
    'PERIOD_MAX_E': PERIOD_MAXe,
    'DISP_X_INE': DISP_Xdp,
    'DISP_Y_INE': DISP_Ydp,
    'ROTATION_INE': ROTdp,
    'VELO_INE': velocity_Xdp,
    'ACCEL_INE': acceleration_Xdp,
    'AXIAL_FORCE_INE': FORCE_Adp,
    'SHEAR_FORCE_INE': FORCE_Sdp,
    'MOMENT_INE': MOMENTdp,
    'AXIAL_RIGIDITY_INE': np.abs(FORCE_Adp),
    'ROTATIONAL_ST_INE': KIdp,
    'LATERAL_ST_Y_INE': KAdp,
    'LATERAL_ST_X_INE': KSdp,
    'PERIOD_MIN_INE': PERIOD_MINp,
    'PERIOD_MAX_INE': PERIOD_MAXp,
}
# Convert to DataFrame
results_df = pd.DataFrame(DATA_TOTAL)
# Export the DataFrame to an Excel file
results_df.to_excel('CONFINEMENT_NONCONSTANT_ELASTIC_OR_INELASTIC_FRAME_FREE_VIBRATION_VISCOUS_DAMPER_RESULTS.xlsx', index=False) 
#%%------------------------------------------------------------------------------
# Print out the state of nodes 3 and 4
ops.printModel("node",3, 4)
# Print out the state of element 1 , 2 and 3
ops.printModel("ele", 1, 2 , 3)
# Print the Model
#printModel()
ops.printModel("-JSON", "-file", "CONFINEMENT_NONCONSTANT_ELASTIC_OR_INELASTIC_FRAME_FREE_VIBRATION_VISCOUS_DAMPER.json")
#%%-------------------------------------------------------------------------------
    
