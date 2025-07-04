################################################################################################################
#                                                  IN THE NAME OF ALLAH                                        #
#                      SEQUENTIAL THERMAL AND SEISMIC ANALYSIS OF A CONCRETE FRAME USING OPENSEES              #
#                                         APPLY THERMAL LOAD ON ALL ELEMENTS                                   #
#--------------------------------------------------------------------------------------------------------------#
#                             THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                       #
#                                        EMAIL: salar.d.ghashghaei@gmail.com                                   #
################################################################################################################

"""
Models and analyzes a 2D concrete frame subjected to both thermal gradients
 and distributed mechanical loads using OpenSees.

Key Features:
[1] Model Definition:
   - A two-dimensional frame is defined with user-specified nodal coordinates for multiple stories and bays.
   - Base supports are fixed.
   - Concrete material properties, including thermal expansion effects, are defined.
   - Beam and column cross-sections are modeled using fiber sections with rectangular concrete patches to capture nonlinear behavior.

[2] Element and Load Assignment:
   - Beam-column elements are created with corotational geometric transformations to accommodate large displacements.
   - Lobatto integration is used for accurate numerical integration of the fiber sections.
   - Uniformly distributed loads are applied to beam elements.
   - A thermal gradient is imposed on the beams of the first story to simulate temperature-induced effects.

[3] Analysis Configuration:
   - Static analysis is conducted using load control with incremental thermal loading.
   - A Newton-Raphson solution algorithm is employed for iterative convergence.
   - Analysis tolerances and maximum iteration limits are set to ensure stability and accuracy.

[4] Output and Post-processing:
   - Node displacements, support reactions, and element deformations are recorded throughout the analysis.
   - Output data is post-processed to extract and plot base reactions (axial force, shear force, and bending moment)
   and nodal displacements as functions of temperature and applied load.

[5] Visualization:
   - Both the undeformed and deformed configurations of the frame are plotted.
   - Key results such as displacement-temperature relationships and base reaction trends are visualized to aid
   interpretation of structural response under thermal and mechanical loading.

"""
#%% --------------------------------------------------------------------------
# PAPER: Implementation of new elements and material models in OpenSees software to account for post-earthquacke fire damage
'https://www.sciencedirect.com/science/article/abs/pii/S2352012420304124'
# PAPER: A framework for performance-based assessment in post-earthquake fire: Methodology and case study
'https://www.sciencedirect.com/science/article/abs/pii/S0141029623011811'
#%% --------------------------------------------------------------------------
import openseespy.opensees as ops
import matplotlib.pyplot as plt
import numpy as np
import time as TI
import ANALYSIS_FUNCTION as S02
import CONCRETE_THERMAL_SECTION_FUN as S03
#import CONCRETE_FIBERTHERMAL_SECTION as S03
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
lamdaC = 0.1        # Ratio between unloading slope at epsU and initial slope
ftC = 0.005         # [N/mm²] Tensile strength
EtsC = ftC/0.002    # [N/mm²] Tension softening stiffness 


# Cover concrete (unconfined)
fcU = -18           # [N/mm²] Concrete Compressive Strength
ec0U = -0.0025      # [mm/mm] Concrete Compressive Strain
fcUU = -2           # [N/mm²] Concrete Compressive Ultimate Strength
ecuU = -0.008       # [mm/mm] Concrete Compressive Ultimate Strain
lamdaU = 0.1        # Ratio between unloading slope at epsU and initial slope
ftU = 0.005         # [N/mm²] Tensile strength
EtsU = ftU/0.002    # [N/mm²] Tension softening stiffness 

 
# STEEL
# Reinforcing steel
fy = 400            # [N/mm²] Steel Rebar Yield Strength   
Es = 2e5            # [N/mm²] Modulus of Elasticity
ey = fy/Es          # [mm/mm] Steel Rebar Yield Strain
fu = 1.1818*fy      # [N/mm²] Steel Rebar Ultimate Strength
esu = ey*75.2       # [mm/mm] Steel Rebar Ultimate Strain
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

# Define Thermal and Distributed load
Max_Thermal = 900.0      # [°C] Temperature
dl = -0.5                # [N/mm] Distributed Load

LENGTH_COL = 3000        # [mm] Column Length 
LENGTH_BM = 7000         # [mm] Beam Length 
# Define Analysis Properties
MAX_ITERATIONS = 10000   # Convergence iteration for test
MAX_TOLERANCE = 1.0e-10  # Convergence tolerance for test

Nstep = 800              # Number of incremental steps
Incr_Temp = 1/Nstep      # Incremental temperature step

GMfact = 9810    # standard acceleration of gravity or standard acceleration
SSF_X = 0.00001  # Seismic Acceleration Scale Factor in X Direction
SSF_Y = 0.00001  # Seismic Acceleration Scale Factor in Y Direction
iv0_X = 0.00005  # [mm/s] Initial velocity applied to the node  in X Direction
iv0_Y = 0.00005  # [mm/s] Initial velocity applied to the node  in Y Direction
st_iv0 = 0.0     # [s] Initial velocity applied starting time
SEI = 'X'        # Seismic Direction
DR = 0.05        # Damping ratio
duration = 15.0  # [s] Total simulation duration
dt = 0.01        # [s] Time step
MASS = np.abs(dl) * LENGTH_BM    # [kg] Mass of the structure
#%%------------------------------------------------------------------------------
def ALL_ANALYSIS():
    # Reset model
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)

    # CORNER NODES
    ops.node(1, 0.0, 0.0)
    ops.node(2, LENGTH_BM, 0.0)
    ops.node(3, 0.0, LENGTH_COL)
    ops.node(4, LENGTH_BM, LENGTH_COL)
    # BEAM NODES
    ops.node(5, 0.25*LENGTH_BM, LENGTH_COL)
    ops.node(6, 0.50*LENGTH_BM, LENGTH_COL)
    ops.node(7, 0.75*LENGTH_BM, LENGTH_COL)
    # RIGHT COLUMN NODES
    ops.node(8, 0.0, 0.25*LENGTH_COL)
    ops.node(9, 0.0, 0.50*LENGTH_COL)
    ops.node(10, 0.0, 0.75*LENGTH_COL)
    # LEFT COLUMN NODES
    ops.node(11, LENGTH_BM, 0.25*LENGTH_COL)
    ops.node(12, LENGTH_BM, 0.50*LENGTH_COL)
    ops.node(13, LENGTH_BM, 0.75*LENGTH_COL)
    print(' Structure Nodes Done.')    
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 1, 1, 1)

    secTagC = 10
    secTagB = 20
    coreTag = 1
    coverTag = 2
    steelTag = 3
    
    # STEEL REBAR MATERIAL PROPERTIES
    ops.uniaxialMaterial('Steel01Thermal', steelTag, fy, Es, Bs) 
    
    # CONCRETE  MATERIAL PROPERTIES
    """
    Concrete02Thermal:
    Concrete02Thermal is created for modelling concrete, which is derived from
    the standard "Concrete02" material and incorprates with temperature dependent
    properties sugggested in Eurocode 2 (EN1992-1-2)
    """
    ops.uniaxialMaterial('Concrete02Thermal', coreTag, fcC, ec0C, fcUC, ecuC, lamdaC, ftC, EtsC)  # Core concrete (confined)
    ops.uniaxialMaterial('Concrete02Thermal', coverTag, fcU, ec0U, fcUU, ecuU, lamdaU, ftU, EtsU) # Cover concrete (unconfined)
    print(' Thermal Material Done.')
    
    # COLUMN THERMAL SECTION
    S03.CONFINED_CONCRETE_SECTION_THERMAL(secTagC, Hc, Bc, coverC, AsC, coreTag, coverTag, steelTag, COL=True)
    # BEAM THERMAL SECTION
    S03.CONFINED_CONCRETE_SECTION_THERMAL(secTagB, Hb, Bb, coverB, AsB, coreTag, coverTag, steelTag, COL=False)
    
    print(' Thermal Section Done.')
    
    # Define geometric transformation (corotational for large displacements)
    transfTag = 1
    #ops.geomTransf('Linear', transfTag)
    #ops.geomTransf('PDelta', transfTag)
    ops.geomTransf('Corotational', transfTag)
    # Define beam integration (Lobatto integration)
    numIntegrationPoints = 5
    biTagB = 1
    ops.beamIntegration('Lobatto', biTagB, secTagB, numIntegrationPoints)
 
    # Define column integration (Lobatto integration)
    numIntegrationPoints = 5
    biTagC = 2
    ops.beamIntegration('Lobatto', biTagC, secTagC, numIntegrationPoints)
    #--------------------------------------------------------------------
    ops.element('dispBeamColumnThermal', 1, 1, 8, transfTag, biTagC)  # COLUMN 01
    ops.element('dispBeamColumnThermal', 2, 8, 9, transfTag, biTagC)  # COLUMN 02
    ops.element('dispBeamColumnThermal', 3, 9, 10, transfTag, biTagC) # COLUMN 03
    ops.element('dispBeamColumnThermal', 4, 10, 3, transfTag, biTagC) # COLUMN 04
    
    ops.element('dispBeamColumnThermal', 5, 2, 11, transfTag, biTagC)  # COLUMN 05
    ops.element('dispBeamColumnThermal', 6, 11, 12, transfTag, biTagC) # COLUMN 06
    ops.element('dispBeamColumnThermal', 7, 12, 13, transfTag, biTagC) # COLUMN 07
    ops.element('dispBeamColumnThermal', 8, 13, 4, transfTag, biTagC)  # COLUMN 08
    
    ops.element('dispBeamColumnThermal', 9, 3, 5, transfTag, biTagB)  # BEAM 01
    ops.element('dispBeamColumnThermal', 10, 5, 6, transfTag, biTagB) # BEAM 02
    ops.element('dispBeamColumnThermal', 11, 6, 7, transfTag, biTagB) # BEAM 03
    ops.element('dispBeamColumnThermal', 12, 7, 4, transfTag, biTagB) # BEAM 04
    print(' Thermal Elements Done.')
    
    # Data storage
    FORCE_S, FORCE_A, MOMENT = [], [], []
    DISP_X, DISP_Y, ROT = [], [], []
    KA, KS, KI, THERMAL_STEP, STEP = [], [], [], [], []
    FORCE_S01, FORCE_A01, MOMENT01 = [], [], []
    DISP_X01, DISP_Y01, ROT01 = [], [], []

    # Define time series and load pattern
    Time_Tag01 = 1
    ops.timeSeries('Linear', Time_Tag01)
    patternTag01 = 1
    ops.pattern('Plain', patternTag01 , Time_Tag01)
    #ops.load(3, 0.0, 0.0, 0.0)
    #ops.load(4, 0.0, 0.0, 0.0)
    # Apply distributed load to all beams in the structure
    # eleLoad -ele $eleTag1 <$eleTag2 ....> -type -beamUniform $Wz <$Wx>
    ops.eleLoad('-ele', 9, '-type', '-beamUniform', dl) 
    ops.eleLoad('-ele', 10, '-type', '-beamUniform', dl)
    ops.eleLoad('-ele', 11, '-type', '-beamUniform', dl)
    ops.eleLoad('-ele', 12, '-type', '-beamUniform', dl)
    print(' Distributed load Done.')
    
    # A linear thermal gradient is applied to the elements in the first story.
    # The temperature varies from -Db to Db across the section height.
    # INFO LINK: https://openseesforfire.github.io/Subpages/ThermalActionCmds.html
    # Apply thermal load only to the identified beam elements
    Db = 0.5 * Hb
    Dc = 0.5 * Hc
    # eleLoad -ele $eleTag -type -beamThermal $T1 $y1 $T2 $Y2
    ops.eleLoad('-ele', 1, '-type', '-beamThermal', Max_Thermal, -Dc, Max_Thermal, Dc)  
    ops.eleLoad('-ele', 2, '-type', '-beamThermal', Max_Thermal, -Dc, Max_Thermal, Dc)  
    ops.eleLoad('-ele', 3, '-type', '-beamThermal', Max_Thermal, -Dc, Max_Thermal, Dc)  
    ops.eleLoad('-ele', 4, '-type', '-beamThermal', Max_Thermal, -Dc, Max_Thermal, Dc)  
    ops.eleLoad('-ele', 5, '-type', '-beamThermal', Max_Thermal, -Dc, Max_Thermal, Dc)  
    ops.eleLoad('-ele', 6, '-type', '-beamThermal', Max_Thermal, -Dc, Max_Thermal, Dc)  
    ops.eleLoad('-ele', 7, '-type', '-beamThermal', Max_Thermal, -Dc, Max_Thermal, Dc)  
    ops.eleLoad('-ele', 8, '-type', '-beamThermal', Max_Thermal, -Dc, Max_Thermal, Dc)  
    ops.eleLoad('-ele', 9, '-type', '-beamThermal', Max_Thermal, -Db, Max_Thermal, Db)  
    ops.eleLoad('-ele', 10, '-type', '-beamThermal', Max_Thermal, -Db, Max_Thermal, Db)  
    ops.eleLoad('-ele', 11, '-type', '-beamThermal', Max_Thermal, -Db, Max_Thermal, Db)  
    ops.eleLoad('-ele', 12, '-type', '-beamThermal', Max_Thermal, -Db, Max_Thermal, Db)  
    print(' Thermal Load Done.')

    #%% RUN THERMAL ANALYSIS
    print('THERMAL ANALYSIS STARTED BY OPENSEES')
    ops.system('BandGeneral')
    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/normDispIncr.html
    ops.algorithm('Newton')
    ops.integrator('LoadControl', Incr_Temp)
    ops.analysis('Static')
    print(' Thermal Analysis Properties Done.')
    
    ops.recorder('Node', '-file', "BTH_PUSH_01.txt",'-time', '-node', 1, '-dof', 1,2,3, 'reaction') # Base Reaction Time History Node 1
    ops.recorder('Node', '-file', "BTH_PUSH_02.txt",'-time', '-node', 2, '-dof', 1,2,3, 'reaction') # Base Reaction Time History Node 2
    ops.recorder('Node', '-file', "DTH_PUSH.txt",'-time', '-node', 3, '-dof', 1,2,3, 'disp')        # Displacement Time History Node 3
    #ops.recorder('Element', '-file', 'STRESS_STRAIN_BOT_CONCRETE.txt', '-time', '-ele', 9, 'section', secTagB, 'fiber', -Db, 0, 'stressStrain')
    #ops.recorder('Element', '-file', 'STRESS_STRAIN_TOP_CONCRETE.txt', '-time', '-ele', 9, 'section', secTagB, 'fiber', +Db, 0, 'stressStrain')
    #ops.recorder('Element', '-file', 'STRESS_STRAIN_BOT_REBAR.txt', '-time', '-ele', 9, 'section', secTagB, 'fiber', -Db+coverB, -0.5*Bb+coverB, 'stressStrain')
    #ops.recorder('Element', '-file', 'STRESS_STRAIN_TOP_REBAR.txt', '-time', '-ele', 9, 'section', secTagB, 'fiber', +Db-coverB, -0.5*Bb+coverB, 'stressStrain')
    
    for step in range(Nstep):
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
        # OUTPUT JUST THERMAL ANALYSIS
        FORCE_S01.append(S)
        FORCE_A01.append(A)
        MOMENT01.append(M)
        DISP_X01.append(disp_X)
        DISP_Y01.append(disp_Y)
        ROT01.append(rot)
        KS.append(np.abs(S)/np.abs(disp_X)) # LATERAL STIFFNESS IN X
        KA.append(np.abs(A)/np.abs(disp_Y)) # LATERAL STIFFNESS IN Y
        KI.append(np.abs(M)/np.abs(rot))    # ROTATIONAL STIFFNESS IN Z
        THERMAL_STEP.append(ops.getLoadFactor(patternTag01) * Max_Thermal) # THERMAL LOAD 
        STEP.append(step+1) # STEPS FOR THERMAL ANALYSIS
        #print(step+1, disp_X, S)
        print(step+1, disp_X, disp_Y)
    
    # Reset analysis (clear analysis info, not the model)
    ops.loadConst('-time', 0.0)
    ops.remove('loadPattern', patternTag01)  # remove pushover loads
    #ops.wipeAnalysis()
    #%% RUN SEISMIC ANALYSIS
    print('SEISMIC ANALYSIS STARTED BY OPENSEES')
    # Define mass
    ops.mass(3, MASS, MASS, 0.0)
    ops.mass(4, MASS, MASS, 0.0)
    # Static 
    Time_Tag02 = 2
    patternTag02 = 2
    ops.timeSeries('Linear', Time_Tag02)
    ops.pattern('Plain', patternTag02 , Time_Tag02)
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
    PERIOD_01 = np.pi / Omega01 # Structure First Period
    PERIOD_02 = np.pi / Omega02 # Structure Second Period
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
    print('Seismic Defined Done.')
    
    time02 = []
    #displacement = []
    velocity_X, velocity_Y = [], []
    acceleration_X, acceleration_Y = [], []
    DISP_X02, DISP_Y02 = [], []
    velocity_X02, velocity_Y02 = [], []
    acceleration_X02, acceleration_Y02 = [], []
    FORCE_S02, FORCE_A02, MOMENT02 = [], [], []
        
    stable = 0
    current_time = 0.0
    
    while stable == 0 and current_time < duration:
        step = step + 1
        STEP.append(step) # STEPS FOR SEISMIC ANALYSIS
        stable = ops.analyze(1, dt)
        S02.ANALYSIS(stable, 1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS
        current_time = ops.getTime()
        time02.append(current_time)
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
        # OUTPUT JUST SEISMIC ANALYSIS
        velocity_X02.append(ops.nodeVel(3, 1))       # LATERAL VELOCITY IN X FOR NODE 3
        acceleration_X02.append(ops.nodeAccel(3, 1)) # LATERAL ACCELERATION IN X FOR NODE 3
        velocity_Y02.append(ops.nodeVel(3, 2))       # LATERAL VELOCITY IN Y FOR NODE 3
        acceleration_Y02.append(ops.nodeAccel(3, 2)) # LATERAL ACCELERATION IN Y FOR NODE 3
        FORCE_S02.append(S)
        FORCE_A02.append(A)
        MOMENT02.append(M)
        DISP_X02.append(disp_X)
        DISP_Y02.append(disp_Y)
        
        ROT.append(rot)
        KS.append(np.abs(S)/np.abs(disp_X)) # LATERAL STIFFNESS IN X
        KA.append(np.abs(A)/np.abs(disp_Y)) # LATERAL STIFFNESS IN Y
        KI.append(np.abs(M)/np.abs(rot))    # ROTATIONAL STIFFNESS IN Z
        print(current_time, disp_X, S)
    
    # Calculating Damping Ratio Using Logarithmic Decrement Analysis 
    displacement = np.array(DISP_X02)
    peaks = np.array([displacement[i] for i in range(1, len(displacement)-1) if displacement[i] > displacement[i-1] and displacement[i] > displacement[i+1]])
    # Natural logarithm
    delta = np.log(peaks[:-1] / peaks[1:])    
    #ops.wipe()    

    return FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y, ROT, KA, KS, KI, THERMAL_STEP, time02, velocity_X, velocity_Y, acceleration_X, acceleration_Y, PERIOD_01, PERIOD_02, delta, STEP, FORCE_S01, FORCE_A01, MOMENT01,DISP_X01, DISP_Y01, ROT01, FORCE_S02, FORCE_A02, MOMENT02, DISP_X02, DISP_Y02, velocity_X02, velocity_Y02, acceleration_X02, acceleration_Y02
    

#%%------------------------------------------------------------------------------
# Analysis Durations:
starttime = TI.process_time()

# WITHOUT HARDENING AND ULTIMATE STRAIN
DATA = ALL_ANALYSIS()
FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y, ROT, KA, KS, KI, THERMAL_STEP, time02, velocity_X, velocity_Y, acceleration_X, acceleration_Y, PERIOD_01, PERIOD_02, delta, STEP, FORCE_S01, FORCE_A01, MOMENT01,DISP_X01, DISP_Y01, ROT01, FORCE_S02, FORCE_A02, MOMENT02, DISP_X02, DISP_Y02, velocity_X02, velocity_Y02, acceleration_X02, acceleration_Y02 = DATA

totaltime = TI.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n') 
#%%------------------------------------------------------------------------------
"""
def OUTPUT_SECOND_COLUMN(X, COLUMN):
    import numpy as np
    # Time History
    filename = f"{X}.txt"
    data_collected = np.loadtxt(filename)
    X = data_collected[:, COLUMN]   
    return X 

disp_X = OUTPUT_SECOND_COLUMN('DTH_PUSH', 1) # Reading Disp from Text file - X Direaction
disp_Y = OUTPUT_SECOND_COLUMN('DTH_PUSH', 2) # Reading Disp from Text file - Y Direaction
disp_Z = OUTPUT_SECOND_COLUMN('DTH_PUSH', 3) # Reading Disp from Text file - Z Direaction (Rotation)
# AXIAL
base01_Y = OUTPUT_SECOND_COLUMN('BTH_PUSH_01', 1) # Reading base reaction from Text file - NODE 1
base02_Y = OUTPUT_SECOND_COLUMN('BTH_PUSH_02', 1) # Reading base reaction from Text file - NODE 2
BASES_AXIAL = base01_Y + base02_Y
# SHEAR
base01_X = OUTPUT_SECOND_COLUMN('BTH_PUSH_01', 2) # Reading base reaction from Text file - NODE 1
base02_X = OUTPUT_SECOND_COLUMN('BTH_PUSH_02', 2) # Reading base reaction from Text file - NODE 2
BASES_SHEAR = base01_X + base02_X
# MOMENT
base01_Z = OUTPUT_SECOND_COLUMN('BTH_PUSH_01', 3) # Reading base reaction from Text file - NODE 1
base02_Z = OUTPUT_SECOND_COLUMN('BTH_PUSH_02', 3) # Reading base reaction from Text file - NODE 2
BASES_MOMENT = base01_Z + base02_Z
# STRESS AND STRAIN OF ELEMENT 6
# CONCRETE
strain_B_C = OUTPUT_SECOND_COLUMN('STRESS_STRAIN_BOT_CONCRETE', 2) # Reading bottom concrete strain from Text file
stress_B_C = OUTPUT_SECOND_COLUMN('STRESS_STRAIN_BOT_CONCRETE', 1) # Reading bottom concrete stress from Text file
strain_T_C = OUTPUT_SECOND_COLUMN('STRESS_STRAIN_TOP_CONCRETE', 2) # Reading top concrete strain from Text file
stress_T_C = OUTPUT_SECOND_COLUMN('STRESS_STRAIN_TOP_CONCRETE', 1) # Reading top concrete stress from Text file
# STEEL REBAR
strain_B_R = OUTPUT_SECOND_COLUMN('STRESS_STRAIN_BOT_REBAR', 2) # Reading bottom steel rebar strain from Text file
stress_B_R = OUTPUT_SECOND_COLUMN('STRESS_STRAIN_BOT_REBAR', 1) # Reading bottom steel rebar stress from Text file
strain_T_R = OUTPUT_SECOND_COLUMN('STRESS_STRAIN_TOP_REBAR', 2) # Reading top steel rebar strain from Text file
stress_T_R = OUTPUT_SECOND_COLUMN('STRESS_STRAIN_TOP_REBAR', 1) # Reading top steel rebar stress from Text file
#%%------------------------------------------------------------------------------

# Plot results
plt.figure(1, figsize=(8, 6))
plt.plot(disp_Y, BASES_AXIAL, color='blue', linewidth=2)
plt.xlabel('Node 6 Displacement Y (mm)')
plt.ylabel('Base Axial Reaction (kN)')
plt.title('Base Axial Reaction vs Displacement Y')
plt.grid()
plt.show()

plt.figure(2, figsize=(8, 6))
plt.plot(disp_X, BASES_SHEAR, color='red', linewidth=2)
plt.xlabel('Node 6 Displacement X (mm)')
plt.ylabel('Base Shear Reaction (kN)')
plt.title('Base Shear Reaction vs Displacement X')
plt.grid()
plt.show()

plt.figure(3, figsize=(8, 6))
plt.plot(disp_Z, BASES_MOMENT, color='green', linewidth=2)
plt.xlabel('Node 6 Rotation (rad)')
plt.ylabel('Base Moment Reaction (kN.mm)')
plt.title('Base Moment Reaction vs Rotation')
plt.grid()
plt.show()

plt.figure(4, figsize=(8, 6))
plt.plot(THERMAL_STEP, DISP_X, color='purple', linewidth=2)
plt.xlabel('Temperature (°C)')
plt.ylabel(f'Node {3} Displacement X (mm)')
plt.title('Temperature vs Displacement X')
plt.grid()
plt.show()

plt.figure(5, figsize=(8, 6))
plt.plot(THERMAL_STEP, DISP_Y, color='black', linewidth=2)
plt.xlabel('Temperature (°C)')
plt.ylabel(f'Node {3} Displacement Y (mm)')
plt.title('Temperature vs Displacement Y')
plt.grid()
plt.show()

plt.figure(6, figsize=(8, 6))
plt.plot(strain_B_C, stress_B_C, color='blue', label='Bottom Fiber', linewidth=2)
plt.plot(strain_T_C, stress_T_C, color='red', label='Top Fiber', linewidth=2)
plt.xlabel('Strain (mm/mm)')
plt.ylabel('Stress (kN/mm^2)')
plt.title(f'Stress-Strain Relation of Element {6} Concrete Top & Bottom Fibers')
plt.grid()
plt.legend()
plt.show()

plt.figure(7, figsize=(8, 6))
plt.plot(strain_B_R, stress_B_R, color='blue', label='Bottom Fiber', linewidth=2)
plt.plot(strain_T_R, stress_T_R, color='red', label='Top Fiber', linewidth=2)
plt.xlabel('Strain (mm/mm)')
plt.ylabel('Stress (kN/mm^2)')
plt.title(f'Stress-Strain Relation of Element {6} Steel Rebar Top & Bottom Fibers')
plt.grid()
plt.legend()
plt.show()
"""
#%%------------------------------------------------------------------------------

plt.figure(1, figsize=(12, 8))
plt.plot(MOMENT, FORCE_A, color='black')
#plt.scatter(MOMENT, FORCE_A, color='blue', linewidth=2)
plt.title('P-M Interaction')
plt.ylabel('Axial Force [N]')
plt.xlabel('Bending Moment [N.mm]')
plt.grid()
plt.show()

plt.figure(2, figsize=(12, 8))
plt.plot(DISP_X, FORCE_S, color='green', linewidth=2)
#plt.scatter(DISP_X, FORCE_S, color='purple', linewidth=2)
plt.title('SHEAR FORCE-DISPLACEMENT DIAGRAM')
plt.ylabel('Shear Force [N]')
plt.xlabel('Displacement  in X [mm]')
plt.grid()
plt.show()

plt.figure(3, figsize=(12, 8))
plt.plot(DISP_Y, FORCE_A, color='purple', linewidth=2)
#plt.scatter(DISP_Y, FORCE_A, color='purple', linewidth=2)
plt.title('AXIAL FORCE-DISPLACEMENT DIAGRAM')
plt.ylabel('Axial Force [N]')
plt.xlabel('Displacement in Y [mm]')
plt.grid()
plt.show()

plt.figure(4, figsize=(12, 8))
plt.plot(ROT, MOMENT, color='red', linewidth=2)
#plt.scatter(ROT, MOMENT, color='red', linewidth=2)
plt.title('MOMENT-ROTATION DIAGRAM')
plt.ylabel('Moment [kN.mm]')
plt.xlabel('Rotation [rad]')
plt.grid()
plt.show()

plt.figure(5, figsize=(12, 8))
#plt.plot(KI, KS, color='black', linewidth=2)
plt.scatter(KI, KS, color='black', linewidth=2)
plt.title('ROTATIONAL STIFFNESS-LATERAL STIFFNESS DIAGRAM')
plt.ylabel('Rotational Stiffness [N.mm/Rad]')
plt.xlabel('Lateral Stiffness in X Dir. [N/mm]')
plt.semilogx()
plt.semilogy()
plt.grid()
plt.show()

plt.figure(6, figsize=(12, 8))
#plt.plot(KI, KA, color='black', linewidth=2)
plt.scatter(KI, KA, color='black', linewidth=2)
plt.title('ROTATIONAL STIFFNESS-LATERAL STIFFNESS DIAGRAM')
plt.ylabel('Rotational Stiffness [N.mm/Rad]')
plt.xlabel('Lateral Stiffness in Y Dir. [N/mm]')
plt.semilogx()
plt.semilogy()
plt.grid()
plt.show()

plt.figure(7, figsize=(12, 8))
plt.plot(STEP, FORCE_A, color='brown', linewidth=2)
#plt.scatter(STEP, FORCE_A, color='brown', linewidth=2)
plt.title('Axial Force During the All Analysises')
plt.ylabel('Axial Force [N]')
plt.xlabel('Load Steps')
plt.grid()
plt.show()

plt.figure(8, figsize=(12, 8))
plt.plot(STEP, FORCE_S, color='purple', linewidth=2)
#plt.scatter(STEP, FORCE_S, color='purple', linewidth=2)
plt.title('Shear Force During the All Analysises')
plt.ylabel('Shear Force [N]')
plt.xlabel('Load Steps')
plt.grid()
plt.show()

plt.figure(9, figsize=(12, 8))
plt.plot(STEP, MOMENT, color='green', linewidth=2)
#plt.scatter(STEP, MOMENT, color='green', linewidth=2)
plt.title('Moment During the All Analysises')
plt.ylabel('Moment [kN.mm]')
plt.xlabel('Load Steps')
plt.grid()
plt.show()

plt.figure(10, figsize=(12, 8))
plt.plot(STEP, DISP_X, color='brown', linewidth=2)
#plt.scatter(STEP, DISP_X, color='brown', linewidth=2)
plt.title('Displacement During the All Analysises')
plt.ylabel('Displacement - X [mm]')
plt.xlabel('Load Steps')
plt.grid()
plt.show()

plt.figure(11, figsize=(12, 8))
plt.plot(STEP, DISP_Y, color='blue', linewidth=2)
#plt.scatter(STEP, DISP_X, color='brown', linewidth=2)
plt.title('Displacement During the All Analysises')
plt.ylabel('Displacement - Y [mm]')
plt.xlabel('Load Steps')
plt.grid()
plt.show()

plt.figure(12, figsize=(12, 8))
plt.plot(STEP, ROT, color='black', linewidth=2)
#plt.scatter(STEP, ROT, color='black', linewidth=2)
plt.title('Rotation During the All Analysises')
plt.ylabel('Rotation [rad]')
plt.xlabel('Load Steps')
plt.grid()
plt.show()

#%% PLOT JUST THERMAL ANALYSIS RESULTS
plt.figure(13, figsize=(12, 8))
plt.plot(THERMAL_STEP, FORCE_A01, color='brown', linewidth=2)
#plt.scatter(THERMAL_STEP, FORCE_A01, color='brown', linewidth=2)
plt.title('Axial Force During the Thermal Analysis')
plt.ylabel('Axial Force [N]')
plt.xlabel('Thermal Load Temperature [°C]')
plt.grid()
plt.show()

plt.figure(14, figsize=(12, 8))
plt.plot(THERMAL_STEP, FORCE_S01, color='purple', linewidth=2)
#plt.scatter(THERMAL_STEP, FORCE_S01, color='purple', linewidth=2)
plt.title('Shear Force During the Thermal Analysis')
plt.ylabel('Shear Force [N]')
plt.xlabel('Thermal Load Temperature [°C]')
plt.grid()
plt.show()

plt.figure(15, figsize=(12, 8))
plt.plot(THERMAL_STEP, MOMENT01, color='green', linewidth=2)
#plt.scatter(THERMAL_STEP, MOMENT01, color='green', linewidth=2)
plt.title('Moment During the Thermal Analysis')
plt.ylabel('Moment [kN.mm]')
plt.xlabel('Thermal Load Temperature [°C]')
plt.grid()
plt.show()

plt.figure(16, figsize=(12, 8))
plt.plot(THERMAL_STEP, DISP_X01, color='brown', linewidth=2)
#plt.scatter(THERMAL_STEP, DISP_X01, color='brown', linewidth=2)
plt.title('Displacement During the Thermal Analysis')
plt.ylabel('Displacement - X [mm]')
plt.xlabel('Thermal Load Temperature [°C]')
plt.grid()
plt.show()

plt.figure(17, figsize=(12, 8))
plt.plot(THERMAL_STEP, DISP_Y01, color='blue', linewidth=2)
#plt.scatter(THERMAL_STEP, DISP_X01, color='brown', linewidth=2)
plt.title('Displacement During the Thermal Analysis')
plt.ylabel('Displacement - Y [mm]')
plt.xlabel('Thermal Load Temperature [°C]')
plt.grid()
plt.show()

plt.figure(18, figsize=(12, 8))
plt.plot(THERMAL_STEP, ROT01, color='black', linewidth=2)
#plt.scatter(THERMAL_STEP, ROT01, color='black', linewidth=2)
plt.title('Rotation During the Thermal Analysis')
plt.ylabel('Rotation [rad]')
plt.xlabel('Thermal Load Temperature [°C]')
plt.grid()
plt.show()

#%%------------------------------------------------------------------------------  
# EXACT SOLUTION OF DAMPING RATIO DURING SEISMIC ANALYSIS:
from scipy.optimize import fsolve

# Define the equation for natural logarithm of this ratio, called the logarithmic decrement, we denote by δ
def EQUATION(x, delta):
    if np.any(x == 0):  # Avoid division by zero
        return np.inf  
        
    # Calculate the value of the equation
    A = x**2 - 1 + ((2 * np.pi * x) / np.mean(delta)) ** 2
    #print(f"x: {x}, A: {A}")  # Debugging output
    # Return the difference (for root finding)
    return A
      
# Initial guess for root(s)
x0 = 1  # Intial Guess for Damping Ratio
# Solve for x
solution = fsolve(EQUATION, x0, args=(delta))
print(f"Exact Damping Ratio: {solution[0]:.8e}")
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

DISP_ZX = MAX_ABS(DISP_X02)  
DISP_ZY = MAX_ABS(DISP_Y02) 
VELO_Z = MAX_ABS(velocity_X02) 
ACCE_Z = MAX_ABS(acceleration_X02) 
BASE_Z = MAX_ABS(FORCE_S02) 

plt.figure(1, figsize=(8, 6))
plt.plot(time02, DISP_X02, color='blue', linewidth=2)
#plt.plot(time02, DISP_ZX, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Displacement in X [mm]')
plt.title(f'Time vs Displacement - MAX. ABS: {DISP_ZX[-1]} | ξ (Calculated): {100*solution[0]:.5e} %')
plt.grid()
plt.show()

plt.figure(2, figsize=(8, 6))
plt.plot(time02, DISP_Y02, color='blue', linewidth=2)
#plt.plot(time02, DISP_ZY, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Displacement in Y [mm]')
plt.title(f'Time vs Displacement - MAX. ABS: {DISP_ZY[-1]}')
plt.grid()
plt.show()

plt.figure(3, figsize=(8, 6))
plt.plot(time02, velocity_X02, color='blue', linewidth=2)
plt.plot(time02, VELO_Z, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Velocity in X [mm/s]')
plt.title(f'Time vs Velocity - MAX. ABS: {VELO_Z[-1]}')
plt.grid()
plt.show()

plt.figure(4, figsize=(8, 6))
plt.plot(time02, acceleration_X02, color='blue', linewidth=2)
plt.plot(time02, ACCE_Z, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Acceleration in X [mm/s^2]')
plt.title(f'Time vs Acceleration - MAX. ABS: {ACCE_Z[-1]}')
plt.grid()
plt.show()

plt.figure(5, figsize=(8, 6))
plt.plot(time02, FORCE_S02, color='blue', linewidth=2)
plt.plot(time02, BASE_Z, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Base-reaction [N]')
plt.title(f'Time vs Base-reaction - MAX. ABS: {BASE_Z[-1]}')
plt.grid()
plt.show()

plt.figure(6, figsize=(8, 6))
plt.plot(DISP_X02, FORCE_S02, color='black', linewidth=5)
plt.xlabel('Displacement in X [mm]')
plt.ylabel('Base Shear Reaction [N]')
plt.title(f'Base Shear Reaction vs Displacement - MAX. ABS: {np.max(np.abs(DISP_X02)): .3f}')
plt.grid()
plt.show()

plt.figure(7, figsize=(8, 6))
plt.plot(DISP_Y02, FORCE_A02, color='black', linewidth=5)
plt.xlabel('Displacement in Y [mm]')
plt.ylabel('Base Axial Reaction [N]')
plt.title(f'Base Axial Reaction vs Displacement - MAX. ABS: {np.max(np.abs(DISP_Y02)): .3f}')
plt.grid()
plt.show()

plt.figure(8, figsize=(8, 6))
plt.plot(DISP_X, FORCE_S, color='black', linewidth=5)
plt.xlabel('Displacement in X [mm]')
plt.ylabel('Base Shear Reaction [N]')
plt.title(f'Base Shear Reaction vs Displacement - MAX. ABS: {np.max(np.abs(DISP_X02)): .3f}')
plt.grid()
plt.show()

plt.figure(9, figsize=(8, 6))
plt.plot(DISP_Y, FORCE_A, color='black', linewidth=5)
plt.xlabel('Displacement in Y [mm]')
plt.ylabel('Base Axial Reaction [N]')
plt.title(f'Base Axial Reaction vs Displacement - MAX. ABS: {np.max(np.abs(DISP_Y02)): .3f}')
plt.grid()
plt.show()
#%%------------------------------------------------------------------------------    
# %% Plot 2D Frame Shapes
S04.PLOT_2D_FRAME(deformed_scale=10)  # Adjust scale factor as needed
#%%------------------------------------------------------------------------------   
# EXPORT DATA TO EXCEL
import pandas as pd
DATA_TOTAL = { 
    'DISP_X_WO': DISP_X, 
    'DISP_Y_WO': DISP_Y,
    'AXIAL_FORCE_WO': FORCE_A,
    'SHEAR_FORCE_WO': FORCE_S,
    'MOMENT_WO': MOMENT,
    'AXIAL_RIGIDITY_WO': np.abs(FORCE_A),
    'ROTATIONAL_ST_WO': KI,
    'LATERAL_ST_Y_WO': KA,
    'LATERAL_ST_X_WO': KS,
}
# Convert to DataFrame
results_df = pd.DataFrame(DATA_TOTAL)
# Export the DataFrame to an Excel file
results_df.to_excel('CONCRETE_FRAME_THERMAL_LOAD-SEISMIC_RESULTS.xlsx', index=False) 
#%%------------------------------------------------------------------------------
# Print out the state of nodes 3 and 4
ops.printModel("node",3, 4)
# Print out the state of element 1 , 2 and 3
ops.printModel("ele", 1, 2 , 3)
# Print the Model
#printModel()
ops.printModel("-JSON", "-file", "CONCRETE_FRAME_THERMAL_LOAD-SEISMSIC.json")
#%%------------------------------------------------------------------------------   
    
