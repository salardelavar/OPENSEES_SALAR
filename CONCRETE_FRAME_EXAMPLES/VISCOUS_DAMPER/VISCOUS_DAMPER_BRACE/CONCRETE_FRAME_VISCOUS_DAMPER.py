######################################################################################################################
#                            >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                        #
#                           ASSESSMENTS OF THE STRUCTURAL DUCTILITY DAMAGE INDEX WITH DIFFERENT                      #
#                   CONFINEMENT ENHANCEMENT RATIO OF CONCRETE FRAME WITH VISCOUS DAMPER USING OPENSEES               #
#--------------------------------------------------------------------------------------------------------------------#
#                              THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                            #
#                                       EMAIL: salar.d.ghashghaei@gmail.com                                          #
######################################################################################################################
"""
[1] Nonlinear Frame Modeling
- 2D RC frame with distributed plasticity using `nonlinearBeamColumn` elements.
- Fiber sections for beams/columns: confined core and unconfined cover concrete.

[2] Material Laws
- Concrete: `Concrete01` for confined/unconfined zones.
- Steel: `Hysteretic` model with pinching, hardening, and cyclic degradation.

[3] Seismic Loads
- Pushover: displacement-controlled lateral loading to failure.
- Dynamic: uniform excitation with user-defined ground motions (X/Y).

[4] Damping
- Rayleigh damping (a0, a1) calibrated via eigenvalue analysis (modes 1–2).

[5] Performance Metrics
- Ductility ratio μ: from bilinearized pushover curve.
- Overstrength Ω₀: yield vs. ultimate capacity.
- Damage Index (DI): normalized displacement demand/capacity.

[6] Advanced Solver
- HHT-α integrator (unconditionally stable) with Newton-Raphson iterations.

[7] Outputs
- Hysteretic responses: P-M, V-Δ, M-θ.
- Time-history plots: displacement, base shear.
- Stiffness degradation tracking.

[8] Validation
- Logarithmic decrement method to verify damping ratio.

[9] Ductility Damage Index (DDI)
- After bilinear fit (Δy = yield disp, Δu = ultimate disp):
    Dd = max dynamic displacement demand.
    DI = (Dd - Dy) / (Du - Dy)
    DI ≈ 0: elastic.
    DI ≥ 1: collapse.

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
import TRUSS_ELEMENT_FUN as S07

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
SSF_X = 1.0      # Seismic Acceleration Scale Factor in X Direction
SSF_Y = 1.0      # Seismic Acceleration Scale Factor in Y Direction
iv0_X = 0.0005   # [mm/s] Initial velocity applied to the node  in X Direction
iv0_Y = 0.0005   # [mm/s] Initial velocity applied to the node  in Y Direction
st_iv0 = 0.0     # [s] Initial velocity applied starting time
SEI = 'X'        # Seismic Direction
DR = 0.05        # Intial Guess for Damping ratio for Structure
duration = 100.0 # [s] Total simulation duration
dt = 0.01        # [s] Time step
MASS = 12000     # [kg] Mass on the each column

FY = 3650        # [N/mm^2] Yield Strength Viscous Damper
DRD = 0.15       # Intial Guess for Damping ratio for Viscous Damper
area = 35**2     # [mm^2] Damper Secton Area

#%% DEFINE PARAMETERS FOR NONLINEAR STATIC ANALYSIS 
DMAX = 220       # [mm] Maximum Displacement
DINCR = 0.05     # [mm] Incremental Displacement

#%% DEFINE ANALYSIS PROPERTIES
MAX_ITERATIONS = 20000    # Convergence iteration for test
MAX_TOLERANCE = 1.0e-6   # Convergence tolerance for test
#STEEL_KIND: 1 -> WITHOUT HARDENING AND ULTIMATE STRAIN
#STEEL_KIND: 2 -> WITH HARDENING AND ULTIMATE STRAIN
#%%------------------------------------------------------------------------------
def PD_ANALYSIS(area, FY, MASS, STEEL_KIND, ANA_KIND, DAMPER_KIND, DR, DRD):
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
    S07.TRUSS_ELEMENT(ele_tag, NODE_I, NODE_J, KP, KN, area, MASS, ELE_KIND, DRD, PLOT=True)
    
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
        
        # Define Loads
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
        ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS)
        ops.integrator('Newmark', 0.5, 0.25) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/newmark.html
        #alpha = 1;gamma=1.5-alpha; beta=((2-alpha)**2)/4;
        #ops.integrator('HHT', alpha, gamma, beta) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/hht.html
        ops.algorithm('Newton') # INFO LINK: https://openseespydoc.readthedocs.io/en/stable/src/algorithm.html
        ops.analysis('Transient')
        
        # Calculate Rayleigh damping factors
        PERIOD_01, PERIOD_02 = S06.RAYLEIGH_DAMPING(2, DR, 0.6*DR, 0, 1)
        
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
            
        # Data storage
        FORCE_S, FORCE_A, MOMENT = [], [], []
        DISP_X, DISP_Y, ROT = [], [], []
        KA, KS, KI, STEP = [], [], [], []
        time = []
        velocity_X, velocity_Y = [], []
        acceleration_X, acceleration_Y = [], []
        DRIFT_X, DRIFT_Y = [], []
            
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
            print(current_time, disp_X, S)
        
        # Calculate Structure Damping Ratio Based on Lateral Displacement
        damping_ratio = S05.DAMPING_RATIO(DISP_X)    
        
        #ops.wipe()  
        return FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y, ROT, KA, KS, KI, time, velocity_X, velocity_Y, acceleration_X, acceleration_Y, PERIOD_01, PERIOD_02, DRIFT_X, DRIFT_Y, damping_ratio


#%%------------------------------------------------------------------------------
# Analysis Durations for Nonlinear Static Analysis:
starttime = TI.process_time()

# WITHOUT HARDENING AND ULTIMATE STRAIN -> STEEL_KIND=1
# WITH HARDENING AND ULTIMATE STRAIN -> STEEL_KIND=2

# RUN NONLINEAR STATIC ANALYSIS
STEEL_KIND = 2
ANA_KIND = 'PUSHOVER'
DAMPER_KIND = 'INELASIC'
DATA = PD_ANALYSIS(area, FY, MASS, STEEL_KIND, ANA_KIND, DAMPER_KIND, DR, DRD)
FORCE_Sp, FORCE_Ap, MOMENTp, DISP_Xp, DISP_Yp, ROTp, KAp, KSp, KIp, DRIFT_Xp, DRIFT_Yp, STEPp = DATA

totaltime = TI.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n') 

# %% Plot 2D Frame Shapes for Nonlinear Static Analysis
S04.PLOT_2D_FRAME(deformed_scale=1)  # Adjust scale factor as needed
#%%------------------------------------------------------------------------------
# Analysis Durations for Nonlinear Dynamic Analysis:
starttime = TI.process_time()

# RUN NONLINEAR DYNAMIC ANALYSIS
STEEL_KIND = 2
ANA_KIND = 'DYNAMIC'
DAMPER_KIND = 'INELASIC'
DATA = PD_ANALYSIS(area, FY, MASS, STEEL_KIND, ANA_KIND, DAMPER_KIND, DR, DRD)
FORCE_Sd, FORCE_Ad, MOMENTd, DISP_Xd, DISP_Yd, ROTd, KAd, KSd, KId, timed, velocity_Xd, velocity_Yd, acceleration_Xd, acceleration_Yd, PERIOD_01d, PERIOD_02d, DRIFT_Xd, DRIFT_Yd, damping_ratio = DATA

print(f"\n Period 01: {PERIOD_01d:.4e} (s) - Period 02: {PERIOD_02d:.4e} (s)")

totaltime = TI.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n') 

# %% Plot 2D Frame Shapes for Nonlinear Dynamic Analysis
S04.PLOT_2D_FRAME(deformed_scale=10000)  # Adjust scale factor as needed
#%%------------------------------------------------------------------------------
plt.figure(1, figsize=(12, 8))
plt.plot(MOMENTp, FORCE_Ap, color='black')
plt.plot(MOMENTd, FORCE_Ad, color='purple')
#plt.scatter(MOMENTp, FORCE_Ap, color='black', linewidth=2)
#plt.scatter(MOMENTd, FORCE_Ad, color='purple', linewidth=2)
plt.title('P-M Interaction')
plt.ylabel('Axial Force [N]')
plt.xlabel('Bending Moment [N.mm]')
plt.legend(['PUSHOVER', 'DYNAMIC'])
plt.grid()
plt.show()

plt.figure(2, figsize=(12, 8))
plt.plot(DISP_Xp, FORCE_Sp, color='green', linewidth=2)
plt.plot(DISP_Xd, FORCE_Sd, color='lime', linewidth=2)
#plt.scatter(DISP_Xp, FORCE_Sp, color='green', linewidth=2)
#plt.scatter(DISP_Xd, FORCE_Sd, color='lime', linewidth=2)
plt.title('SHEAR FORCE-DISPLACEMENT DIAGRAM')
plt.ylabel('Shear Force [N]')
plt.xlabel('Displacement  in X [mm]')
plt.legend(['PUSHOVER', 'DYNAMIC'])
plt.grid()
plt.show()

plt.figure(3, figsize=(12, 8))
plt.plot(DISP_Yp, FORCE_Ap, color='purple', linewidth=2)
plt.plot(DISP_Yd, FORCE_Ad, color='green', linewidth=2)
#plt.scatter(DISP_Yp, FORCE_Ap, color='purple', linewidth=2)
#plt.scatter(DISP_Yd, FORCE_Ad, color='green', linewidth=2)
plt.title('AXIAL FORCE-DISPLACEMENT DIAGRAM')
plt.ylabel('Axial Force [N]')
plt.xlabel('Displacement in Y [mm]')
plt.legend(['PUSHOVER', 'DYNAMIC'])
plt.grid()
plt.show()

plt.figure(4, figsize=(12, 8))
plt.plot(ROTp, MOMENTp, color='red', linewidth=2)
plt.plot(ROTd, MOMENTd, color='pink', linewidth=2)
#plt.scatter(ROTp, MOMENTp, color='red', linewidth=2)
#plt.scatter(ROTd, MOMENTd, color='pink', linewidth=2)
plt.title('MOMENT-ROTATION DIAGRAM')
plt.ylabel('Moment [kN.mm]')
plt.xlabel('Rotation [rad]')
plt.legend(['PUSHOVER', 'DYNAMIC'])

plt.grid()
plt.show()

plt.figure(5, figsize=(12, 8))
#plt.plot(KIp, KSp, color='black', linewidth=2)
#plt.plot(KId, KSd, color='grey', linewidth=2)
plt.scatter(KIp, KSp, color='black', linewidth=2)
plt.scatter(KId, KSd, color='grey', linewidth=2)
plt.title('ROTATIONAL STIFFNESS-LATERAL STIFFNESS DIAGRAM')
plt.ylabel('Rotational Stiffness [N.mm/Rad]')
plt.xlabel('Lateral Stiffness in X Dir. [N/mm]')
plt.semilogx()
plt.semilogy()
plt.legend(['PUSHOVER', 'DYNAMIC'])
plt.grid()
plt.show()

plt.figure(6, figsize=(12, 8))
#plt.plot(KIp, KAp, color='black', linewidth=2)
#plt.plot(KId, KAd, color='grey', linewidth=2)
plt.scatter(KIp, KAp, color='black', linewidth=2)
plt.scatter(KId, KAd, color='grey', linewidth=2)
plt.title('ROTATIONAL STIFFNESS-LATERAL STIFFNESS DIAGRAM')
plt.ylabel('Rotational Stiffness [N.mm/Rad]')
plt.xlabel('Lateral Stiffness in Y Dir. [N/mm]')
plt.semilogx()
plt.semilogy()
plt.legend(['PUSHOVER', 'DYNAMIC'])
plt.grid()
plt.show()

plt.figure(7, figsize=(12, 8))
plt.plot(timed, FORCE_Ad, color='brown', linewidth=2)
#plt.scatter(timed, FORCE_Ad, color='brown', linewidth=2)
plt.title('Axial Force During the Analysis')
plt.ylabel('Axial Force [N]')
plt.xlabel('Times')
plt.grid()
plt.show()

plt.figure(8, figsize=(12, 8))
plt.plot(timed, FORCE_Sd, color='purple', linewidth=2)
#plt.scatter(timed, FORCE_Sd, color='purple', linewidth=2)
plt.title('Shear Force During the Analysis')
plt.ylabel('Shear Force [N]')
plt.xlabel('Times')
plt.grid()
plt.show()

plt.figure(9, figsize=(12, 8))
plt.plot(timed, MOMENTd, color='green', linewidth=2)
#plt.scatter(timed, MOMENTd, color='green', linewidth=2)
plt.title('Moment During the Analysis')
plt.ylabel('Moment [N.mm]')
plt.xlabel('Times')
plt.grid()
plt.show()

plt.figure(10, figsize=(12, 8))
plt.plot(timed, DISP_Xd, color='brown', linewidth=2)
#plt.scatter(timed, DISP_Xd, color='brown', linewidth=2)
plt.title('Displacement During the Analysis')
plt.ylabel('Displacement - X [mm]')
plt.xlabel('Times')
plt.grid()
plt.show()

plt.figure(11, figsize=(12, 8))
plt.plot(timed, DISP_Yd, color='blue', linewidth=2)
#plt.scatter(timed, DISP_Yd, color='brown', linewidth=2)
plt.title('Displacement During the Analysis')
plt.ylabel('Displacement - Y [mm]')
plt.xlabel('Times')
plt.grid()
plt.show()

plt.figure(12, figsize=(12, 8))
plt.plot(timed, ROTd, color='black', linewidth=2)
#plt.scatter(timed, ROTd, color='black', linewidth=2)
plt.title('Rotation During the Analysis')
plt.ylabel('Rotation [rad]')
plt.xlabel('Times')
plt.grid()
plt.show()

plt.figure(13, figsize=(12, 8))
plt.plot(STEPp, DRIFT_Xp, color='green', linewidth=2)
plt.plot(STEPp, DRIFT_Yp, color='purple', linewidth=2)
plt.title('Structure Lateral Drift During the Pushover Analysis')
plt.ylabel('Lateral Drift [%]')
plt.xlabel('Steps')
plt.legend([f'IN X DIR.: {np.max(np.abs(DRIFT_Xp)):.3e}', f'IN Y DIR.: {np.max(np.abs(DRIFT_Yp)):.3e}'])
plt.grid()
plt.show()

plt.figure(14, figsize=(12, 8))
plt.plot(timed, DRIFT_Xd, color='red', linewidth=2)
plt.plot(timed, DRIFT_Yd, color='blue', linewidth=2)
plt.title('Structure Lateral Drift During the Dynamic Analysis')
plt.ylabel('Lateral Drift [%]')
plt.xlabel('Times')
plt.legend([f'IN X DIR.: {np.max(np.abs(DRIFT_Xd)):.3e}', f'IN Y DIR.: {np.max(np.abs(DRIFT_Yd)):.3e}'])
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

DISP_ZX = MAX_ABS(DISP_Xd)  
DISP_ZY = MAX_ABS(DISP_Yd) 
VELO_Z = MAX_ABS(velocity_Xd) 
ACCE_Z = MAX_ABS(acceleration_Xd) 
BASE_Z = MAX_ABS(FORCE_Sd) 
BASE_ZY = MAX_ABS(FORCE_Ad)

plt.figure(1, figsize=(8, 6))
plt.plot(timed, DISP_Xd, color='blue', linewidth=2)
plt.plot(timed, DISP_ZX, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Displacement in X [mm]')
plt.title(f'Time vs Displacement - MAX. ABS: {DISP_ZX[-1]} | ξ (Calculated): {damping_ratio:.5e} %')
plt.grid()
plt.show()

plt.figure(2, figsize=(8, 6))
plt.plot(timed, DISP_Yd, color='blue', linewidth=2)
plt.plot(timed, DISP_ZY, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Displacement in Y [mm]')
plt.title(f'Time vs Displacement - MAX. ABS: {DISP_ZY[-1]}')
plt.grid()
plt.show()

plt.figure(3, figsize=(8, 6))
plt.plot(timed, velocity_Xd, color='blue', linewidth=2)
plt.plot(timed, VELO_Z, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Velocity in X [mm/s]')
plt.title(f'Time vs Velocity - MAX. ABS: {VELO_Z[-1]}')
plt.grid()
plt.show()

plt.figure(4, figsize=(8, 6))
plt.plot(timed, acceleration_Xd, color='blue', linewidth=2)
plt.plot(timed, ACCE_Z, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Acceleration in X [mm/s^2]')
plt.title(f'Time vs Acceleration - MAX. ABS: {ACCE_Z[-1]}')
plt.grid()
plt.show()

plt.figure(5, figsize=(8, 6))
plt.plot(timed, FORCE_Sd, color='blue', linewidth=2)
plt.plot(timed, BASE_Z, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Base-reaction [N]')
plt.title(f'Time vs Base-reaction - MAX. ABS: {BASE_Z[-1]}')
plt.grid()
plt.show()

plt.figure(6, figsize=(8, 6))
plt.plot(timed, FORCE_Ad, color='blue', linewidth=2)
plt.plot(timed, BASE_ZY, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Axial Base-reaction [N]')
plt.title(f'Time vs Axial Base-reaction - MAX. ABS: {BASE_ZY[-1]}')
plt.grid()
plt.show()

plt.figure(7, figsize=(8, 6))
plt.plot(DISP_Xd, FORCE_Sd, color='black', linewidth=2)
plt.xlabel('Displacement in X [mm]')
plt.ylabel('Shear Base-reaction [N]')
plt.title('Displacement vs Shear Base-reaction')
plt.grid()
plt.show()

plt.figure(8, figsize=(8, 6))
plt.plot(DISP_Yd, FORCE_Ad, color='brown', linewidth=2)
plt.xlabel('Displacement in Y [mm]')
plt.ylabel('Axial Base-reaction [N]')
plt.title('Displacement vs Axial Base-reaction')
plt.grid()
plt.show()

plt.figure(9, figsize=(8, 6))
plt.plot(ROTd, MOMENTd, color='red', linewidth=2)
plt.xlabel('Rotation [rad]')
plt.ylabel('Moment Base-reaction [N.mm/Rad]')
plt.title('Rotation vs Moment Base-reaction')
plt.grid()
plt.show()

#%%------------------------------------------------------------------------------

# --------------------------------------
#  Plot BaseShear-Displacement Analysis 
# --------------------------------------
XX = np.abs(DISP_Xp); YY = np.abs(FORCE_Sp); # ABSOLUTE VALUE
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
BC.PLOT_2D(np.abs(DISP_Xp), np.abs(FORCE_Sp), X, Y, X, Y, XLABEL, YLABEL, TITLE, LEGEND01, LEGEND02, LEGEND03, COLOR='black', Z=2) 
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
Dd = np.max(np.abs(DISP_Xd))
DIx = (Dd - X[1]) /(X[2] - X[1])
print(f'Structural Ductility Damage Index in X Direction: {DIx:.4f}')
# ---------------------------------------
#  Plot BaseAxial-Displacement Analysis
# ---------------------------------------
XX = np.abs(DISP_Yp); YY = np.abs(FORCE_Ap); # ABSOLUTE VALUE
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
BC.PLOT_2D(np.abs(DISP_Yp), np.abs(FORCE_Ap), X, Y, X, Y, XLABEL, YLABEL, TITLE, LEGEND01, LEGEND02, LEGEND03, COLOR='black', Z=2) 
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
Dd = np.max(np.abs(DISP_Yd))
DIy = (Dd - X[1]) /(X[2] - X[1])
print(f'Structural Ductility Damage Index in Y Direction: {DIy:.4f}')
#%%------------------------------------------------------------------------------   
# EXPORT DATA TO EXCEL
import pandas as pd
DATA_TOTAL = { 
    'DISP_X_WO': DISP_Xd,
    'DISP_Y_WO': DISP_Yd,
    'ROTATION_WO': ROTd,
    'VELO_WO': velocity_Xd,
    'ACCEL_WO': acceleration_Xd,
    'AXIAL_FORCE_WO': FORCE_Ad,
    'SHEAR_FORCE_WO': FORCE_Sd,
    'MOMENT_WO': MOMENTd,
    'AXIAL_RIGIDITY_WO': np.abs(FORCE_Ad),
    'ROTATIONAL_ST_WO': KId,
    'LATERAL_ST_Y_WO': KAd,
    'LATERAL_ST_X_WO': KSd,
}
# Convert to DataFrame
results_df = pd.DataFrame(DATA_TOTAL)
# Export the DataFrame to an Excel file
results_df.to_excel('CONCRETE_FRAME_VISCOUS_DAMPER_RESULTS.xlsx', index=False) 
#%%------------------------------------------------------------------------------
# Print out the state of nodes 3 and 4
ops.printModel("node",3, 4)
# Print out the state of element 1 , 2 and 3
ops.printModel("ele", 1, 2 , 3)
# Print the Model
#printModel()
ops.printModel("-JSON", "-file", "CONCRETE_FRAME_VISCOUS_DAMPER.json")
#%%-------------------------------------------------------------------------------
    
