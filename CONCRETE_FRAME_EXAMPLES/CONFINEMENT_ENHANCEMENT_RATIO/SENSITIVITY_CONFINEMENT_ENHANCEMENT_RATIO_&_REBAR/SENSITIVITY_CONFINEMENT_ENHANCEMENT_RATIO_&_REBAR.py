################################################################################################################
#                          >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                    #
#  SENSITIVITY ANALYSIS OF CONCRETE FRAME BY CHANGING COLUMN REBAR DIAMETER AND CONFINEMENT ENHANCEMENT RATIO  # 
#                        USING OPENSEES FOR STRUCTURAL BEHAVIOR COEFFICIENT CALCULATION                        #
#--------------------------------------------------------------------------------------------------------------#
#                             THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                       #
#                                        EMAIL: salar.d.ghashghaei@gmail.com                                   #
################################################################################################################
"""
1. Objective: The code performs a sensitivity analysis on a 2D reinforced concrete frame by varying column rebar
 diameter and confinement enhancement ratios to evaluate structural behavior coefficients using OpenSees.  

2. Model Setup: A nonlinear 2D frame model is created with columns, beams, and distributed plasticity elements,
 incorporating geometric transformations (Corotational) for large displacements.  

3. Material Modeling: Confined and unconfined concrete behaviors are modeled using modified Kent-Scott-Park
 formulations, while steel reinforcement follows bilinear or hardening models.  

4. Analysis Types: Both pushover (static) and dynamic analyses are supported, with Rayleigh damping calibrated to modal properties for dynamic cases.  

5. Key Outputs: The code extracts base shear, displacement, stiffness, ductility ratios, overstrength
 factors, and structural behavior coefficients (R).  

6. Sensitivity Parameters: Rebar diameters (20–32 mm) and confinement enhancement ratios (1.15–1.35)
 are systematically varied to assess their impact on performance.  

7. Bilinear Fitting: Pushover curves are post-processed to derive elastic/plastic stiffness, ductility,
 and R factors using a bilinear approximation algorithm.  

8. Visualization: 3D contour plots and 2D graphs illustrate relationships between rebar, confinement,
 and structural responses (e.g., stiffness, R-factors).  

9. Validation: Eigenvalue analysis ensures realistic dynamic properties (periods/damping), while
 convergence checks enhance numerical robustness.  

10. Applications: The tool aids in performance-based design, code compliance (e.g., FEMA, Eurocode),
 and optimizing reinforcement strategies for seismic resilience.  

This approach combines advanced nonlinear modeling, parametric studies, and post-processing to quantify
 how material and geometric variables influence seismic performance.
"""
import openseespy.opensees as ops
import matplotlib.pyplot as plt
import numpy as np
import time as TI
import ANALYSIS_FUNCTION as S02
import CONCRETE_SECTION_FUN_TWO as S03
import PLOT_2D as S04
import BILINEAR_CURVE as BC

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
SSF_X = 0.01     # Seismic Acceleration Scale Factor in X Direction
SSF_Y = 0.01     # Seismic Acceleration Scale Factor in Y Direction
iv0_X = 0.0005   # [mm/s] Initial velocity applied to the node  in X Direction
iv0_Y = 0.0005   # [mm/s] Initial velocity applied to the node  in Y Direction
st_iv0 = 0.0     # [s] Initial velocity applied starting time
SEI = 'X'        # Seismic Direction
DR = 0.05        # Intial Guess for Damping ratio
duration = 50.0  # [s] Total simulation duration
dt = 0.01        # [s] Time step
MASS = 12000     # [kg] Mass on the each column

#%% DEFINE PARAMETERS FOR NONLINEAR STATIC ANALYSIS 
DMAX = 675       # [mm] Maximum Displacement
DINCR = 0.05     # [mm] Incremental Displacement

#%% DEFINE ANALYSIS PROPERTIES
MAX_ITERATIONS = 5000    # Convergence iteration for test
MAX_TOLERANCE = 1.0e-6   # Convergence tolerance for test
#STEEL_KIND: 1 -> WITHOUT HARDENING AND ULTIMATE STRAIN
#STEEL_KIND: 2 -> WITH HARDENING AND ULTIMATE STRAIN
#%%------------------------------------------------------------------------------
def PD_ANALYSIS(DIAc, Kfc, STEEL_KIND, ANA_KIND):
    
    #DIAc = 25                # [mm] # Rebar Size Diameter
    AsC = np.pi*(DIAc**2)/4   # [mm²] Area of Rebar
    RO_COL = (8*AsC)/(Hc*Bc)  # Column Rebar Ratio
    
    DIAb = DIAc-4             # [mm] # Rebar Size Diameter
    AsB = np.pi*(DIAb**2)/4   # [mm²] Area of Rebar
    RO_BE = (6*AsB)/(Hb*Bb)   # Beam Rebar Ratio
    
    # Reset model
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)

    # CORNE NODES
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
        
        return RO_COL, RO_BE, FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y, ROT, KA, KS, KI, Elastic_ST, Plastic_ST, Tangent_ST, Ductility_Rito, Over_Strength_Factor, R, STEP        
    
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
        ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS)
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
        print('Structure First Period:  ', PERIOD_01)
        print('Structure Second Period: ', PERIOD_02) 
        
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
        DISP_X, DISP_Y, ROT = [], [], []
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
        
        # Calculating Damping Ratio Using Logarithmic Decrement Analysis 
        displacement = np.array(DISP_X)
        peaks = np.array([displacement[i] for i in range(1, len(displacement)-1) if displacement[i] > displacement[i-1] and displacement[i] > displacement[i+1]])
        # Natural logarithm
        delta = np.log(peaks[:-1] / peaks[1:])    
        
        #ops.wipe()  
        return RO_COL, RO_BE, FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y, ROT, KA, KS, KI, time, velocity_X, velocity_Y, acceleration_X, acceleration_Y, PERIOD_01, PERIOD_02, delta


#%%------------------------------------------------------------------------------
# SENSTIVITY ANALYSIS BY CHANGING REBAR DIAMETER AND CONFINEMENT ENHANCEMENT RATIO 
# Analysis Durations:
current_time = TI.strftime("%H:%M:%S", TI.localtime())
print(f"Current time (HH:MM:SS): {current_time}\n\n")

RO_COL_MAX, RO_BE_MAX =  [], []
FORCE_S_MAX, FORCE_A_MAX, MOMENT_MAX = [], [], []
DISP_X_MAX, DISP_Y_MAX, ROT_MAX = [], [], []
KA_MAX, KS_MAX, KI_MAX, R_MAX = [], [], [], []
DIAc_MAX, Kfc_MAX = [], []
Elastic_ST_MAX, Plastic_ST_MAX, Tangent_ST_MAX, Ductility_Rito_MAX, Over_Strength_Factor_MAX = [], [], [], [], []

# REBAR DIAMETER FOR COLUMN SECTION
#DIAc = [8, 10, 12, 14, 16, 18, 20, 22, 25, 28, 30, 32]
DIAc = [20, 22, 28, 30, 32]
# CONFINEMENT ENHANCEMENT FACTOR (K)
Kfc = [1.15, 1.20, 1.25, 1.30, 1.35]
II = 0
for dia in DIAc:
    for kf in Kfc:
        II = II + 1
        print(f'\n STEP: {II} - REBAR DIAMETER: {dia} - K: {kf} \n')
        DIAc_MAX.append(dia)
        Kfc_MAX.append(kf)
        DATA = PD_ANALYSIS(dia, kf, STEEL_KIND=2, ANA_KIND='PUSHOVER')
        (RO_COL, RO_BE, FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y,
         ROT, KA, KS, KI, Elastic_ST, Plastic_ST,
         Tangent_ST, Ductility_Rito, Over_Strength_Factor, R, STEP) = DATA
        
        RO_COL_MAX.append(RO_COL)
        RO_BE_MAX.append(RO_BE)
        FORCE_S_MAX.append(np.max(np.abs(FORCE_S)))
        FORCE_A_MAX.append(np.max(np.abs(FORCE_A)))
        MOMENT_MAX.append(np.max(np.abs(MOMENT)))
        DISP_X_MAX.append(np.max(np.abs(DISP_X)))
        DISP_Y_MAX.append(np.max(np.abs(DISP_Y)))
        ROT_MAX.append(np.max(np.abs(ROT)))
        KA_MAX.append(np.max(np.abs(KA)))
        KS_MAX.append(np.max(np.abs(KS)))
        KI_MAX.append(np.max(np.abs(KI)))
        Elastic_ST_MAX.append(Elastic_ST)
        Plastic_ST_MAX.append(Plastic_ST)
        Tangent_ST_MAX.append(Tangent_ST)
        Ductility_Rito_MAX.append(Ductility_Rito)
        Over_Strength_Factor_MAX.append(Over_Strength_Factor)
        R_MAX.append(R)
        
current_time = TI.strftime("%H:%M:%S", TI.localtime())
print(f"Current time (HH:MM:SS): {current_time}\n\n")
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
X = DIAc_MAX
Y = Kfc_MAX
Z = FORCE_S_MAX
XLABEL = 'Rebar Diameter'   
YLABEL = 'Confinement Enhancement Ratio' 
ZLABEL = 'Structure Shear Reaction in X Dir. [N]'
PLOT_3D(1, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = DIAc_MAX
Y = Kfc_MAX
Z =  FORCE_A_MAX
XLABEL = 'Rebar Diameter'   
YLABEL = 'Confinement Enhancement Ratio' 
ZLABEL = 'Structure Axial Reaction in Y Dir. [N]'
PLOT_3D(2, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = DIAc_MAX
Y = Kfc_MAX
Z =  MOMENT_MAX
XLABEL = 'Rebar Diameter'   
YLABEL = 'Confinement Enhancement Ratio' 
ZLABEL = 'Structure Moment Reaction in Y Dir. [N.mm]'
PLOT_3D(3, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = DIAc_MAX
Y = Kfc_MAX
Z =  DISP_X_MAX
XLABEL = 'Rebar Diameter'   
YLABEL = 'Confinement Enhancement Ratio' 
ZLABEL = 'Structure Displacement in X Dir. [mm]'
PLOT_3D(4, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = DIAc_MAX
Y = Kfc_MAX
Z =  DISP_Y_MAX
XLABEL = 'Rebar Diameter'   
YLABEL = 'Strength Enhancement Factor' 
ZLABEL = 'Structure Displacement in Y Dir. [mm]'
PLOT_3D(5, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = DIAc_MAX
Y = Kfc_MAX
Z =  Elastic_ST_MAX
XLABEL = 'Rebar Diameter'   
YLABEL = 'Confinement Enhancement Ratio' 
ZLABEL = 'Structure Elastic Stiffness [N/mm]'
PLOT_3D(6, X, Y, Z, XLABEL, YLABEL, ZLABEL)    

X = DIAc_MAX
Y = Kfc_MAX
Z =  Plastic_ST_MAX
XLABEL = 'Rebar Diameter'   
YLABEL = 'Confinement Enhancement Ratio' 
ZLABEL = 'Structure Plastic Stiffness [N/mm]'
PLOT_3D(7, X, Y, Z, XLABEL, YLABEL, ZLABEL)    

X = DIAc_MAX
Y = Kfc_MAX
Z =  Ductility_Rito_MAX
XLABEL = 'Rebar Diameter'   
YLABEL = 'Confinement Enhancement Ratio' 
ZLABEL = 'Ductility Rito  [mm/mm]'
PLOT_3D(8, X, Y, Z, XLABEL, YLABEL, ZLABEL)   

X = DIAc_MAX
Y = Kfc_MAX
Z =  Over_Strength_Factor_MAX
XLABEL = 'Rebar Diameter'   
YLABEL = 'Confinement Enhancement Ratio' 
ZLABEL = 'Over Strength Factor [N/N]'
PLOT_3D(9, X, Y, Z, XLABEL, YLABEL, ZLABEL)   

X = DIAc_MAX
Y = Kfc_MAX
Z = R_MAX
XLABEL = 'Rebar Diameter'   
YLABEL = 'Strength Enhancement Factor' 
ZLABEL = 'Structural Behavior Coefficient'
PLOT_3D(10, X, Y, Z, XLABEL, YLABEL, ZLABEL)  
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

# Call plots using the functions
PLOT_LINE('P-M Interaction', 'Bending Moment [N.mm]', 'Axial Force [N]', MOMENT, FORCE_A, color='black', fig_num=1)
PLOT_LINE('SHEAR FORCE-DISPLACEMENT DIAGRAM', 'Displacement in X [mm]', 'Shear Force [N]', DISP_X, FORCE_S, color='green', fig_num=2)
PLOT_LINE('AXIAL FORCE-DISPLACEMENT DIAGRAM', 'Displacement in Y [mm]', 'Axial Force [N]', DISP_Y, FORCE_A, color='purple', fig_num=3)
PLOT_LINE('MOMENT-ROTATION DIAGRAM', 'Rotation [rad]', 'Moment [kN.mm]', ROT, MOMENT, color='red', fig_num=4)

PLOT_SCATTER('ROTATIONAL STIFFNESS-LATERAL STIFFNESS DIAGRAM (X Dir)', 'Lateral Stiffness in X Dir. [N/mm]', 'Rotational Stiffness [N.mm/Rad]', KI, KS, color='black', fig_num=5, logx=True, logy=True)
PLOT_SCATTER('ROTATIONAL STIFFNESS-LATERAL STIFFNESS DIAGRAM (Y Dir)', 'Lateral Stiffness in Y Dir. [N/mm]', 'Rotational Stiffness [N.mm/Rad]', KI, KA, color='black', fig_num=6, logx=True, logy=True)

PLOT_LINE('Axial Force During the Analysis', 'Steps', 'Axial Force [N]', STEP, FORCE_A, color='brown', fig_num=7)
PLOT_LINE('Shear Force During the Analysis', 'Steps', 'Shear Force [N]', STEP, FORCE_S, color='purple', fig_num=8)
PLOT_LINE('Moment During the Analysis', 'Steps', 'Moment [kN.mm]', STEP, MOMENT, color='green', fig_num=9)
PLOT_LINE('Displacement During the Analysis (X)', 'Steps', 'Displacement - X [mm]', STEP, DISP_X, color='brown', fig_num=10)
PLOT_LINE('Displacement During the Analysis (Y)', 'Steps', 'Displacement - Y [mm]', STEP, DISP_Y, color='blue', fig_num=11)
PLOT_LINE('Rotation During the Analysis', 'Steps', 'Rotation [rad]', STEP, ROT, color='black', fig_num=12)

#%%------------------------------------------------------------------------------
#import BILINEAR_CURVE as BC

# ---------------------------------------
#  Plot BaseShear-Displacement Analysis
# ---------------------------------------
XX = np.abs(DISP_X); YY = np.abs(FORCE_S); # ABSOLUTE VALUE
SLOPE_NODE = 10

DATA = BC.BILNEAR_CURVE(XX, YY, SLOPE_NODE)
X, Y, Elastic_ST, Plastic_ST, Tangent_ST, Ductility_Rito, Over_Strength_Factor = DATA

XLABEL = 'Displacement in X [mm]'
YLABEL = 'Base-Shear [N]'
LEGEND01 = 'Curve'
LEGEND02 = 'Bilinear Fitted'
LEGEND03 = 'Undefined'
TITLE = f'Last Data of BaseShear-Displacement Analysis - Ductility Ratio: {X[2]/X[1]:.4f} - Over Strength Factor: {Y[2]/Y[1]:.4f}'
COLOR = 'black'
BC.PLOT_2D(np.abs(DISP_X), np.abs(FORCE_S), X, Y, X, Y, XLABEL, YLABEL, TITLE, LEGEND01, LEGEND02, LEGEND03, COLOR='black', Z=2) 
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

# ---------------------------------------
#  Plot BaseAxial-Displacement Analysis
# ---------------------------------------
XX = np.abs(DISP_Y); YY = np.abs(FORCE_A); # ABSOLUTE VALUE
SLOPE_NODE = 10

DATA = BC.BILNEAR_CURVE(XX, YY, SLOPE_NODE)
X, Y, Elastic_ST, Plastic_ST, Tangent_ST, Ductility_Rito, Over_Strength_Factor = DATA

XLABEL = 'Displacement in Y [mm]'
YLABEL = 'Base-Axial [N]'
LEGEND01 = 'Curve'
LEGEND02 = 'Bilinear Fitted'
LEGEND03 = 'Undefined'
TITLE = f'Last Data of BaseAxial-Displacement Analysis - Ductility Ratio: {X[2]/X[1]:.4f} - Over Strength Factor: {Y[2]/Y[1]:.4f}'
COLOR = 'black'
BC.PLOT_2D(np.abs(DISP_Y), np.abs(FORCE_A), X, Y, X, Y, XLABEL, YLABEL, TITLE, LEGEND01, LEGEND02, LEGEND03, COLOR='black', Z=2) 
print(f'\t\t Ductility Ratio: {YY[2]/YY[1]:.4f}')

# Calculate Over Strength Coefficient (Ω0)
Omega_0 = Y[2] / Y[1]
# Calculate Displacement Ductility Ratio (μ)
mu = X[2] / X[1]
# Calculate Ductility Coefficient (Rμ)
R_mu = (2 * mu - 1) ** 0.5 / mu ** 0.5
# Calculate Structural Behavior Coefficient (R)
R = Omega_0 * R_mu
print(f'Over Strength Coefficient (Ω0):      {Omega_0:.4f}')
print(f'Displacement Ductility Ratio (μ):    {mu:.4f}')
print(f'Ductility Coefficient (Rμ):          {R_mu:.4f}')
print(f'Structural Behavior Coefficient (R): {R:.4f}')

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
results_df.to_excel('SENSITIVITY_CONFINEMENT_ENHANCEMENT_RATIO_&_REBAR_RESULTS.xlsx', index=False) 
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
ops.printModel("-JSON", "-file", "SENSITIVITY_CONFINEMENT_ENHANCEMENT_RATIO_&_REBAR.json")
#%%------------------------------------------------------------------------------

    