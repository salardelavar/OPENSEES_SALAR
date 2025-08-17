################################################################################################################
#                          >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                    #
#     SENSITIVITY ANALYSIS OF ELASTIC CONCRETE FRAME BEHAVIOR: INVESTIGATING THE IMPACT OF COLUMN HEIGHT AND   #
#           BEAM LENGTH AND MASS ON STRUCTURAL PERIOD AND OTHER KEY PARAMETERS USING OPENSEES AND PYTHON       #
#--------------------------------------------------------------------------------------------------------------#
#                             THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                       #
#                                        EMAIL: salar.d.ghashghaei@gmail.com                                   #
################################################################################################################
"""
# Linear Dynamic and Sensitivity Analysis of a Concrete Frame Using OpenSees  
---------------------------------------------------------------
This study performs a comprehensive linear dynamic analysis and sensitivity assessment
 of a reinforced concrete frame structure using OpenSees. 
 The research focuses on evaluating the structural response by varying two key parameters:  
1. Coumn Height – Examining how different column lengths influence dynamic behavior  
2. Beam length – Examining how different span lengths influence dynamic behavior  
3. Structural mass – Investigating the effect of mass variation on seismic performance  
---------------------------------------------------------------
## Methodology  
1. Model Development  
   - Create a linear finite element model of a reinforced concrete moment-resisting frame  
   - Implement fiber sections with appropriate material models (Concrete02, Steel02)  
   - Include geometric linearities (P-Delta effects)  

2. Parameter Variation  
   - Beam lengths: ±20% variation from baseline design  
   - Mass modifications: ±30% variation to represent different loading conditions  

3. Analysis Procedures  
   - linear Dynamic Analysis:  
     - Apply earthquake ground motions (e.g., El Centro, Kobe records)  
     - Evaluate displacement demands, story drifts, and damage progression  
   - Sensitivity Analysis:  
     - Perform parametric studies by systematically varying beam length and mass  
     - Quantify influence on key response parameters (peak displacements, base shear)  

4. Implementation in OpenSees  
   - Utilize python scripting for automated parametric analysis  
   - Employ Python post-processing for results visualization  
---------------------------------------------------------------
## Expected Outcomes  
- Identification of critical parameters affecting seismic performance  
- Quantitative assessment of beam length and mass sensitivity  
- Practical insights for performance-based design of concrete frames  

This analysis framework provides a robust approach for evaluating parameter sensitivity in
 structural systems using OpenSees' advanced linear analysis capabilities.
"""
#%%------------------------------------------------------------------------------

import openseespy.opensees as ops
import matplotlib.pyplot as plt
import numpy as np
import time as TI
import ANALYSIS_FUNCTION as S02
import CONCRETE_SECTION_FUN as S03
import PLOT_2D as S04


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
fy = 4000         # [N/mm²] Steel Rebar Yield Strength   
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
DIAc = 25                # [mm] # Rebar Size Diameter
AsC = np.pi*(DIAc**2)/4  # [mm²] Area of Rebar

# Beam Section
Bb = 500                 # [mm] Depth of the Section 
Hb = 300                 # [mm] Height of the Section  
coverB = 50              # [mm] Concrete Section Cover
DIAb = 18                # [mm] # Rebar Size Diameter
AsB = np.pi*(DIAb**2)/4  # [mm²] Area of Rebar

#%% DEFINE THE ELEMENTS LENGTH
LENGTH_COL = 3000        # [mm] Column Length 
LENGTH_BM = 7000         # [mm] Beam Length 

#%% DEFINE PARAMEETRS FOR NONLINEAR DYNAMIC ANALYSIS
GMfact = 9810    # [mm/s²]standard acceleration of gravity or standard acceleration
SSF_X = 1.0      # Seismic Acceleration Scale Factor in X Direction
SSF_Y = 1.0      # Seismic Acceleration Scale Factor in Y Direction
iv0_X = 0.0005   # [mm/s] Initial velocity applied to the node  in X Direction
iv0_Y = 0.0005   # [mm/s] Initial velocity applied to the node  in Y Direction
st_iv0 = 0.0     # [s] Initial velocity applied starting time
SEI = 'X'        # Seismic Direction
DR = 0.05        # Intial Guess for Damping ratio
duration = 15.0  # [s] Total simulation duration
dt = 0.01        # [s] Time step
MASS = 12000     # [kg] Mass on the each column

#%% DEFINE PARAMETERS FOR NONLINEAR STATIC ANALYSIS 
DMAX = 675       # [mm] Maximum Displacement
DINCR = 0.05     # [mm] Incremental Displacement

#%% DEFINE ANALYSIS PROPERTIES
MAX_ITERATIONS = 20000   # Convergence iteration for test
MAX_TOLERANCE = 1.0e-4   # Convergence tolerance for test
#STEEL_KIND: 1 -> WITHOUT HARDENING AND ULTIMATE STRAIN
#STEEL_KIND: 2 -> WITH HARDENING AND ULTIMATE STRAIN
#%%------------------------------------------------------------------------------
def PD_ANALYSIS(LENGTH_COL, LENGTH_BM, MASS, STEEL_KIND, ANA_KIND):
    # Reset model
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)

    # Nodes
    ops.node(1, 0.0, 0.0)
    ops.node(2, LENGTH_BM, 0.0)
    ops.node(3, 0.0, LENGTH_COL)
    ops.node(4, LENGTH_BM, LENGTH_COL)
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
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/elasticBeamColumn.html
    ops.element('elasticBeamColumn', 1, 1, 3, secTagC, transfTag, '-mass', Bc*Hc*0.000025) # COLUMN 01
    ops.element('elasticBeamColumn', 2, 2, 4, secTagC, transfTag, '-mass', Bc*Hc*0.000025) # COLUMN 02
    ops.element('elasticBeamColumn', 3, 3, 4, secTagB, transfTag, '-mass', Bb*Hb*0.000025) # BEAM 01
    
    if ANA_KIND == 'PUSHOVER':
        WEIGHT = MASS * GMfact
        # Data storage
        FORCE_S, FORCE_A, MOMENT = [], [], []
        DISP_X, DISP_Y, ROT = [], [], []
        KA, KS, KI, STEP, DRIFT = [], [], [], [], []
    
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
            DRIFT.append(100*disp_X/LENGTH_COL)
            print(step+1, disp_X, S)
        return FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y, ROT, KA, KS, KI, STEP, DRIFT        
    
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
            DRIFT.append(100*disp_X/LENGTH_COL)
            #print(current_time, disp_X, S)
        
        # Calculating Damping Ratio Using Logarithmic Decrement Analysis 
        displacement = np.array(DISP_X)
        peaks = np.array([displacement[i] for i in range(1, len(displacement)-1) if displacement[i] > displacement[i-1] and displacement[i] > displacement[i+1]])
        # Natural logarithm
        delta = np.log(peaks[:-1] / peaks[1:]) 
        
        #ops.wipe()  
        return FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y, ROT, KA, KS, KI, time, velocity_X, velocity_Y, acceleration_X, acceleration_Y, PERIOD_01, PERIOD_02, delta, DRIFT


#%%------------------------------------------------------------------------------
# SENSTIVITY ANALYSIS BY CHANGING COLUMN HEIGHT AND BEAM LENGTH AND MASS
# Analysis Durations:
current_time = TI.strftime("%H:%M:%S", TI.localtime())
print(f"Current time (HH:MM:SS): {current_time}\n\n")

RO_COL_MAX, RO_BE_MAX =  [], []
FORCE_S_MAX, FORCE_A_MAX, MOMENT_MAX = [], [], []
DISP_X_MAX, DISP_Y_MAX, ROT_MAX = [], [], []
KA_MAX, KS_MAX, KI_MAX = [], [], []
HL_MAX, BL_MAX, MASS_MAX = [], [], []
velocity_X_MAX, velocity_Y_MAX, acceleration_X_MAX, acceleration_Y_MAX, PERIOD_01_MAX, PERIOD_02_MAX = [], [], [], [], [], []
DRIFT_MAX = []
# COLUMN HEIGHT
HL = [ 2000, 2500, 3000, 3500, 4000, 5000, 6000, 7000]
# BEAM LENGTH
BL = [ 2000, 2500, 3000, 3500, 4000, 5000, 6000, 7000]
# MASS ON THE TOP OF EACH COLUMN
mass = [4000, 6000, 8000, 10000, 12000, 14000]
      
II = 0
for hl in HL:
    for bl in BL:
        for m in mass:
            II = II + 1
            print(f'\n STEP: {II} - COL HEIGHT: {hl} - BEAM LENGTH: {bl} - MASS: {m} \n')
            HL_MAX.append(hl)
            BL_MAX.append(bl)
            MASS_MAX.append(m)
            DATA = PD_ANALYSIS(hl, bl, m, STEEL_KIND=2, ANA_KIND='DYNAMIC')
            FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y, ROT, KA, KS, KI, time, velocity_X, velocity_Y, acceleration_X, acceleration_Y, PERIOD_01, PERIOD_02, delta, DRIFT = DATA
            FORCE_S_MAX.append(np.max(np.abs(FORCE_S)))
            FORCE_A_MAX.append(np.max(np.abs(FORCE_A)))
            MOMENT_MAX.append(np.max(np.abs(MOMENT)))
            DISP_X_MAX.append(np.max(np.abs(DISP_X)))
            DISP_Y_MAX.append(np.max(np.abs(DISP_Y)))
            ROT_MAX.append(np.max(np.abs(ROT)))
            KA_MAX.append(np.max(np.abs(KA)))
            KS_MAX.append(np.max(np.abs(KS)))
            KI_MAX.append(np.max(np.abs(KI)))
            
            velocity_X_MAX.append(np.max(np.abs(velocity_X)))
            velocity_Y_MAX.append(np.max(np.abs(velocity_Y))) 
            acceleration_X_MAX.append(np.max(np.abs(acceleration_X)))
            acceleration_Y_MAX.append(np.max(np.abs(acceleration_Y)))
            
            PERIOD_01_MAX.append(np.max(np.abs(PERIOD_01)))
            PERIOD_02_MAX.append(np.max(np.abs(PERIOD_02)))
    
            DRIFT_MAX.append(np.max(np.abs(DRIFT)))
        
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

X = BL_MAX
Y = MASS_MAX
Z = DISP_X_MAX
XLABEL = 'Beam Length [mm]'   
YLABEL = 'Mass [kg]' 
ZLABEL = 'Structure Displacement in X Dir. [mm]'
PLOT_3D(1, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = BL_MAX
Y = MASS_MAX
Z =  DISP_Y_MAX
XLABEL = 'Beam Length [mm]'   
YLABEL = 'Mass [kg]' 
ZLABEL = 'Structure Displacement in Y Dir. [mm]'
PLOT_3D(2, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = BL_MAX
Y = MASS_MAX
Z =  KS_MAX
XLABEL = 'Beam Length [mm]'   
YLABEL = 'Mass [kg]' 
ZLABEL = 'Structure Lateral Stiffness in X Dir. [N/mm]'
PLOT_3D(3, X, Y, Z, XLABEL, YLABEL, ZLABEL)    

X = BL_MAX
Y = MASS_MAX
Z =  velocity_X_MAX
XLABEL = 'Beam Length [mm]'   
YLABEL = 'Mass [kg]' 
ZLABEL = 'Velocity in X Dir. [mm/s]'
PLOT_3D(4, X, Y, Z, XLABEL, YLABEL, ZLABEL)   

X = BL_MAX
Y = MASS_MAX
Z =  velocity_Y_MAX
XLABEL = 'Beam Length [mm]'   
YLABEL = 'Mass [kg]' 
ZLABEL = 'Velocity in Y Dir. [mm/s]'
PLOT_3D(5, X, Y, Z, XLABEL, YLABEL, ZLABEL)   

X = BL_MAX
Y = MASS_MAX
Z =  acceleration_X_MAX
XLABEL = 'Beam Length [mm]'   
YLABEL = 'Mass [kg]' 
ZLABEL = 'Acceleration in X Dir. [mm/s^2]'
PLOT_3D(6, X, Y, Z, XLABEL, YLABEL, ZLABEL)   

X = BL_MAX
Y = MASS_MAX
Z =  acceleration_Y_MAX
XLABEL = 'Beam Length [mm]'   
YLABEL = 'Mass [kg]' 
ZLABEL = 'Acceleration in Y Dir. [mm/s^2]'
PLOT_3D(7, X, Y, Z, XLABEL, YLABEL, ZLABEL)   

X = BL_MAX
Y = MASS_MAX
Z =  PERIOD_01_MAX
XLABEL = 'Beam Length [mm]'   
YLABEL = 'Mass [kg]' 
ZLABEL = 'Min. Period [s]'
PLOT_3D(8, X, Y, Z, XLABEL, YLABEL, ZLABEL)   

X = BL_MAX
Y = MASS_MAX
Z =  PERIOD_02_MAX
XLABEL = 'Beam Length [mm]'   
YLABEL = 'Mass [kg]' 
ZLABEL = 'Max. Period [s]'
PLOT_3D(9, X, Y, Z, XLABEL, YLABEL, ZLABEL)   

X = BL_MAX
Y = MASS_MAX
Z = DRIFT_MAX
XLABEL = 'Beam Length [mm]'   
YLABEL = 'Mass [kg]'
ZLABEL = 'Structural Lateral Drift [%]'
PLOT_3D(10, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = BL_MAX
Y = HL_MAX
Z =  PERIOD_01_MAX
XLABEL = 'Beam Length [mm]'   
YLABEL = 'Column Height [mm]' 
ZLABEL = 'Min. Period [s]'
PLOT_3D(11, X, Y, Z, XLABEL, YLABEL, ZLABEL)   

X = BL_MAX
Y = HL_MAX
Z =  PERIOD_02_MAX
XLABEL = 'Beam Length [mm]'   
YLABEL = 'Column Height [mm]' 
ZLABEL = 'Max. Period [s]'
PLOT_3D(12, X, Y, Z, XLABEL, YLABEL, ZLABEL) 

X = PERIOD_01_MAX
Y = PERIOD_02_MAX
Z = DRIFT_MAX
XLABEL = 'Min. Period [s]'  
YLABEL = 'Max. Period [s]'  
ZLABEL = 'Structural Lateral Drift [%]'
PLOT_3D(13, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = PERIOD_01_MAX
Y = PERIOD_02_MAX
Z = KS_MAX
XLABEL = 'Min. Period [s]'  
YLABEL = 'Max. Period [s]'  
ZLABEL = 'Structural Lateral Stiffness in X Dir. [N/mm]'
PLOT_3D(14, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = PERIOD_01_MAX
Y = PERIOD_02_MAX
Z = KA_MAX
XLABEL = 'Min. Period [s]'  
YLABEL = 'Max. Period [s]'  
ZLABEL = 'Structural Lateral Stiffness in Y Dir. [N/mm]'
PLOT_3D(14, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = PERIOD_01_MAX
Y = PERIOD_02_MAX
Z = KI_MAX
XLABEL = 'Min. Period [s]'  
YLABEL = 'Max. Period [s]'  
ZLABEL = 'Structural Rotational Stiffness in Z Dir. [N/mm]'
PLOT_3D(15, X, Y, Z, XLABEL, YLABEL, ZLABEL)
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

PLOT_LINE('Axial Force During the Analysis', 'Times', 'Axial Force [N]', time, FORCE_A, color='brown', fig_num=7)
PLOT_LINE('Shear Force During the Analysis', 'Times', 'Shear Force [N]', time, FORCE_S, color='purple', fig_num=8)
PLOT_LINE('Moment During the Analysis', 'Times', 'Moment [N.mm]', time, MOMENT, color='green', fig_num=9)
PLOT_LINE('Displacement During the Analysis (X)', 'Times', 'Displacement - X [mm]', time, DISP_X, color='brown', fig_num=10)
PLOT_LINE('Displacement During the Analysis (Y)', 'Times', 'Displacement - Y [mm]', time, DISP_Y, color='blue', fig_num=11)
PLOT_LINE('Rotation During the Analysis', 'Steps', 'Rotation [rad]', time, ROT, color='black', fig_num=12)
#%%------------------------------------------------------------------------------      
# %% Plot 2D Frame Shapes for Nonlinear Dynamic Analysis
S04.PLOT_2D_FRAME(deformed_scale=10)  # Adjust scale factor as needed
#%%------------------------------------------------------------------------------
plt.figure(1, figsize=(12, 8))
plt.plot(MOMENT, FORCE_A, color='purple')
#plt.scatter(MOMENT, FORCE_A, color='purple', linewidth=2)
plt.title('P-M Interaction')
plt.ylabel('Axial Force [N]')
plt.xlabel('Bending Moment [N.mm]')
plt.grid()
plt.show()

plt.figure(2, figsize=(12, 8))
plt.plot(DISP_X, FORCE_S, color='lime', linewidth=2)
#plt.scatter(DISP_X, FORCE_S, color='green', linewidth=2)
plt.title('SHEAR FORCE-DISPLACEMENT DIAGRAM')
plt.ylabel('Shear Force [N]')
plt.xlabel('Displacement  in X [mm]')
plt.grid()
plt.show()

plt.figure(3, figsize=(12, 8))
plt.plot(DISP_Y, FORCE_A, color='green', linewidth=2)
#plt.scatter(DISP_Y, FORCE_A, color='purple', linewidth=2)
plt.title('AXIAL FORCE-DISPLACEMENT DIAGRAM')
plt.ylabel('Axial Force [N]')
plt.xlabel('Displacement in Y [mm]')
plt.grid()
plt.show()

plt.figure(4, figsize=(12, 8))
plt.plot(ROT, MOMENT, color='pink', linewidth=2)
#plt.scatter(ROT, MOMENT, color='pink', linewidth=2)
plt.title('MOMENT-ROTATION DIAGRAM')
plt.ylabel('Moment [N.mm]')
plt.xlabel('Rotation [rad]')

plt.grid()
plt.show()

plt.figure(5, figsize=(12, 8))
#plt.plot(KI, KS, color='grey', linewidth=2)
plt.scatter(KI, KS, color='grey', linewidth=2)
plt.title('ROTATIONAL STIFFNESS-LATERAL STIFFNESS DIAGRAM')
plt.ylabel('Rotational Stiffness [N.mm/Rad]')
plt.xlabel('Lateral Stiffness in X Dir. [N/mm]')
plt.semilogx()
plt.semilogy()
plt.grid()
plt.show()

plt.figure(6, figsize=(12, 8))
#plt.plot(KI, KA, color='grey', linewidth=2)
plt.scatter(KI, KA, color='grey', linewidth=2)
plt.title('ROTATIONAL STIFFNESS-LATERAL STIFFNESS DIAGRAM')
plt.ylabel('Rotational Stiffness [N.mm/Rad]')
plt.xlabel('Lateral Stiffness in Y Dir. [N/mm]')
plt.semilogx()
plt.semilogy()
plt.grid()
plt.show()

plt.figure(7, figsize=(12, 8))
plt.plot(time, FORCE_A, color='brown', linewidth=2)
#plt.scatter(time, FORCE_A, color='brown', linewidth=2)
plt.title('Axial Force During the Analysis')
plt.ylabel('Axial Force [N]')
plt.xlabel('Times')
plt.grid()
plt.show()

plt.figure(8, figsize=(12, 8))
plt.plot(time, FORCE_S, color='purple', linewidth=2)
#plt.scatter(time, FORCE_S, color='purple', linewidth=2)
plt.title('Shear Force During the Analysis')
plt.ylabel('Shear Force [N]')
plt.xlabel('Times')
plt.grid()
plt.show()

plt.figure(9, figsize=(12, 8))
plt.plot(time, MOMENT, color='green', linewidth=2)
#plt.scatter(time, MOMENT, color='green', linewidth=2)
plt.title('Moment During the Analysis')
plt.ylabel('Moment [N.mm]')
plt.xlabel('Times')
plt.grid()
plt.show()

plt.figure(10, figsize=(12, 8))
plt.plot(time, DISP_X, color='brown', linewidth=2)
#plt.scatter(time, DISP_X, color='brown', linewidth=2)
plt.title('Displacement During the Analysis')
plt.ylabel('Displacement - X [mm]')
plt.xlabel('Times')
plt.grid()
plt.show()

plt.figure(11, figsize=(12, 8))
plt.plot(time, DISP_Y, color='blue', linewidth=2)
#plt.scatter(time, DISP_Y, color='brown', linewidth=2)
plt.title('Displacement During the Analysis')
plt.ylabel('Displacement - Y [mm]')
plt.xlabel('Times')
plt.grid()
plt.show()

plt.figure(12, figsize=(12, 8))
plt.plot(time, ROT, color='black', linewidth=2)
#plt.scatter(time, ROT, color='black', linewidth=2)
plt.title('Rotation During the Analysis')
plt.ylabel('Rotation [rad]')
plt.xlabel('Times')
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

DISP_ZX = MAX_ABS(DISP_X)  
DISP_ZY = MAX_ABS(DISP_Y) 
VELO_Z = MAX_ABS(velocity_X) 
ACCE_Z = MAX_ABS(acceleration_X) 
BASE_Z = MAX_ABS(FORCE_S) 

plt.figure(1, figsize=(8, 6))
plt.plot(time, DISP_X, color='blue', linewidth=2)
plt.plot(time, DISP_ZX, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Displacement in X [mm]')
plt.title(f'Time vs Displacement - MAX. ABS: {DISP_ZX[-1]}')
plt.grid()
plt.show()

plt.figure(2, figsize=(8, 6))
plt.plot(time, DISP_Y, color='blue', linewidth=2)
plt.plot(time, DISP_ZY, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Displacement in Y [mm]')
plt.title(f'Time vs Displacement - MAX. ABS: {DISP_ZY[-1]}')
plt.grid()
plt.show()

plt.figure(3, figsize=(8, 6))
plt.plot(time, velocity_X, color='blue', linewidth=2)
plt.plot(time, VELO_Z, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Velocity in X [mm/s]')
plt.title(f'Time vs Velocity - MAX. ABS: {VELO_Z[-1]}')
plt.grid()
plt.show()

plt.figure(4, figsize=(8, 6))
plt.plot(time, acceleration_X, color='blue', linewidth=2)
plt.plot(time, ACCE_Z, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Acceleration in X [mm/s^2]')
plt.title(f'Time vs Acceleration - MAX. ABS: {ACCE_Z[-1]}')
plt.grid()
plt.show()

plt.figure(5, figsize=(8, 6))
plt.plot(time, FORCE_S, color='blue', linewidth=2)
plt.plot(time, BASE_Z, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Base-reaction [N]')
plt.title(f'Time vs Base-reaction - MAX. ABS: {BASE_Z[-1]}')
plt.grid()
plt.show()

plt.figure(6, figsize=(8, 6))
plt.plot(DISP_X, FORCE_S, color='black', linewidth=2)
plt.xlabel('Displacement in X [mm]')
plt.ylabel('Shear Base-reaction [N]')
plt.title(f'Displacement vs Shear Base-reaction')
plt.grid()
plt.show()

plt.figure(7, figsize=(8, 6))
plt.plot(DISP_Y, FORCE_A, color='brown', linewidth=2)
plt.xlabel('Displacement in Y [mm]')
plt.ylabel('Axial Base-reaction [N]')
plt.title(f'Displacement vs Axial Base-reaction')
plt.grid()
plt.show()

#%%------------------------------------------------------------------------------   
# EXPORT DATA TO EXCEL
import pandas as pd
DATA_TOTAL = { 
    'COLUMN_HEIGHT': HL_MAX,
    'BEAM_LENGTH': BL_MAX,
    'MASS': MASS_MAX,
    'DISP_X': DISP_X_MAX, 
    'DISP_Y': DISP_Y_MAX,
    'ROTATION': ROT_MAX,
    'AXIAL_FORCE': FORCE_A_MAX,
    'SHEAR_FORCE': FORCE_S_MAX,
    'MOMENT_WO': MOMENT_MAX,
    'ROTATIONAL_ST': KI_MAX,
    'LATERAL_ST_Y': KA_MAX,
    'LATERAL_ST_X': KS_MAX,
    
    'velocity_X': velocity_X_MAX,
    'velocity_Y': velocity_Y_MAX,
    'acceleration_X': acceleration_X_MAX,
    'acceleration_Y': acceleration_Y_MAX,
    'PERIOD_MIN': PERIOD_01_MAX,
    'PERIOD_MAX': PERIOD_02_MAX,
}
# Convert to DataFrame
results_df = pd.DataFrame(DATA_TOTAL)
# Export the DataFrame to an Excel file
results_df.to_excel('ELASTIC_CONCRETE_FRAME_SENSITIVITY_COLUMN_HEIGHT_&_BEAM_LENGTH_&_MASS_RESULTS.xlsx', index=False) 
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
ops.printModel("node",1, 2, 3, 4)
# Print out the state of element 1 , 2 and 3
ops.printModel("ele", 1, 2 , 3)
# Print the Model
#printModel()
ops.printModel("-JSON", "-file", "ELASTIC_CONCRETE_FRAME_SENSITIVITY_COLUMN_HEIGHT_&_BEAM_LENGTH_&_MASS.json")
#%%------------------------------------------------------------------------------  