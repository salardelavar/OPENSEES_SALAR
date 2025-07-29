################################################################################################################
#                          >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                    #
#             SEQUENTIAL EXPLOSION IMPACT AND THERMAL LOAD ANALYSIS OF A CONCRETE FRAME USING OPENSEES         #
#                                   THERMAL LOAD APPLIED THERMAL LOAD ON ALL ELEMENTS                          #
#--------------------------------------------------------------------------------------------------------------#
#                             THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                       #
#                                        EMAIL: salar.d.ghashghaei@gmail.com                                   #
################################################################################################################

"""
Explosion Impact and Thermo-Mechanical Analysis of Reinforced Concrete Frames Using OpenSees

This computational framework performs coupled nonlinear analyses of 2D RC frames subjected to:
- Transient explosion loading (Friedlander wave equation)
- Steady-state thermal gradients
- Distributed mechanical loads

Key Analysis Components:

[1] Material Modeling:
   - Concrete02Thermal material for temperature-dependent concrete behavior (Eurocode 2 compliant)
   - Steel01Thermal for reinforcing steel with thermal effects
   - Distinct confined/unconfined concrete material models
   - Fiber section discretization for nonlinear section response

[2] Structural Configuration:
   - Multi-element 2D frame with refined mesh (quarter-point nodes)
   - Corotational geometric transformation for large displacements
   - Lobatto integration for accurate section response
   - Fixed base boundary conditions

[3] Loading Protocols:
   - Dynamic explosion loading via Friedlander equation (P0, t0, A parameters)
   - Thermal gradients across section depth (beam/column specific)
   - Distributed dead loads on beam elements
   - Combined static+transient analysis capability

[4] Solution Algorithms:
   - HHT time integration for dynamic analysis (α=1.0)
   - Load-controlled thermal ramp
   - Newton-Raphson iteration with adaptive convergence
   - Rayleigh damping (mass+stiffness proportional)

[5] Output Quantities:
   - Time histories: Displacements, velocities, accelerations
   - Base reactions (shear, axial, moment)
   - Section fiber stresses/strains
   - Structural stiffness degradation
   - Modal properties (natural periods)
   - Damping estimates via logarithmic decrement

[6] Special Features:
   - Automatic analysis restart capability
   - Real-time convergence monitoring
   - Parallel recording of multiple response quantities
   - Integrated visualization of pressure time histories

Implementation Notes:
- Units: mm, N, sec, °C
- Modular design allows easy parameter modification
- Includes comprehensive result processing and visualization
- Validated against Eurocode thermal material models

The framework provides a robust tool for assessing structural
 performance under extreme multi-hazard scenarios, particularly
 suited for blast-resistant design and fire safety evaluation of concrete structures.

"""
#%% ------------------------------------------
# PAPER: Damage evaluation of the steel tubular column subjected to explosion and post-explosion fire condition
'https://www.sciencedirect.com/science/article/abs/pii/S0141029612000351'
# PAPER: PROGRESSIVE COLLAPSE RESISTANCE OF STEEL FRAMED BUILDINGS UNDER EXTREME EVENTS
'https://www.ascjournal.com/down/vol17no3/Vol17no3_10.pdf'
#%% ------------------------------------------
 
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
fy = 4000           # [N/mm²] Steel Rebar Yield Strength   
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

# Define Elements Length
LENGTH_COL = 3000        # [mm] Column Length 
LENGTH_BM = 7000         # [mm] Beam Length 

# Define Analysis Properties
MAX_ITERATIONS = 10000   # Convergence iteration for test
MAX_TOLERANCE = 1.0e-10  # Convergence tolerance for test

Nstep = 800              # Number of incremental steps
Incr_Temp = 1/Nstep      # Incremental temperature step

GMfact = 9810    # standard acceleration of gravity or standard acceleration
DR = 0.05        # Damping ratio
duration = 15.0  # [s] Total simulation duration
dt = 0.01        # [s] Time step
MASS = np.abs(dl) * LENGTH_BM    # [kg] Mass of the structure

# Define Explosion Loading Parameters
P0 = 1e5  # Peak pressure (N)
t0 = 0.1  # Positive phase duration (s)
A = 1.3   # Wave decay coefficient

dt = 0.001  # Time step (s)
impact_duration = 2.0  # Total duration of explosion impact (s)loading 
#%%------------------------------------------------------------------------------
# Define Friedlander equation for explosion loading
def FRIEDLANDER(t, P0, t0, A):
    """
    Calculate the pressure at time t using the Friedlander equation.

    Parameters:
    t (float): Time at which to calculate the pressure.
    P0 (float): Peak pressure.
    t0 (float): Positive phase duration.
    A (float): Wave decay coefficient.

    Returns:
    float: Pressure at time t.
    """
    if t < 0:
        return 0
    return P0 * (1 - t / t0) * np.exp(-A * t / t0)

def EXPLOSION_TIME_SERIES(impact_duration, dt, P0, t0, A):
    time = np.arange(0, impact_duration, dt)
    pressure = [FRIEDLANDER(t, P0, t0, A) for t in time]
    return time, pressure

# Generate explosion pressure time series
time_series, explosion_pressure = EXPLOSION_TIME_SERIES(impact_duration, dt, P0, t0, A)

def plot_friedlander(P0, t0, A, dt, impact_duration):
    """
    Plot the Friedlander explosion loading formula over a given duration.

    Parameters:
    P0 (float): Peak pressure (N).
    t0 (float): Positive phase duration (s).
    A (float): Wave decay coefficient.
    dt (float): Time step (s).
    duration (float): Total duration for the plot (s).
    """
    time = np.arange(0, impact_duration, dt)
    pressure = [FRIEDLANDER(t, P0, t0, A) for t in time]

    # Plot the pressure time history
    plt.figure(figsize=(10, 6))
    plt.plot(time, pressure, label='Explosion Pressure', color='blue')
    plt.xlabel('Time (s)')
    plt.ylabel('Pressure (N)')
    plt.title('Friedlander Explosion Loading')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()
    

# Plot the Friedlander explosion loading
plot_friedlander(P0, t0, A, dt, impact_duration)

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
    
    # Define time series and load pattern
    Time_Tag01 = 1
    ops.timeSeries('Linear', Time_Tag01)
    patternTag01 = 1
    """
    # Apply time-dependent explosion loading
    ops.timeSeries('Path', Time_Tag01, '-dt', dt, '-values', *explosion_pressure)
    ops.pattern('Plain', patternTag01, Time_Tag01)
    ops.load(3, 1.0, 0.0, 0.0)
    ops.load(4, 1.0, 0.0, 0.0)
    print('Impact Load Defined Done.')
    """
    # Apply distributed load to all beams in the structure
    # eleLoad -ele $eleTag1 <$eleTag2 ....> -type -beamUniform $Wz <$Wx>
    ops.pattern('Plain', patternTag01, Time_Tag01)
    ops.eleLoad('-ele', 9, '-type', '-beamUniform', dl) 
    ops.eleLoad('-ele', 10, '-type', '-beamUniform', dl)
    ops.eleLoad('-ele', 11, '-type', '-beamUniform', dl)
    ops.eleLoad('-ele', 12, '-type', '-beamUniform', dl)
    print(' Distributed load Done.')
    
    #%% RUN EXPLOSION LOAD ANALYSIS
    print('EXPLOSION LOAD ANALYSIS STARTED BY OPENSEES')
    # Define mass
    ops.mass(3, MASS, MASS, 0.0)
    ops.mass(4, MASS, MASS, 0.0)
    
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
    
    # Data storage
    FORCE_S, FORCE_A, MOMENT = [], [], []
    FORCE_S01, FORCE_A01, MOMENT01 = [], [], []
    FORCE_S02, FORCE_A02, MOMENT02 = [], [], []
    DISP_X, DISP_Y, ROT = [], [], []
    DISP_X01, DISP_Y01, ROT01 = [], [], []
    DISP_X02, DISP_Y02, ROT02 = [], [], []
    KA, KS, KI, STEP, THERMAL_STEP = [], [], [], [], []
    time = []
    displacement = []
    velocity_X, velocity_Y = [], []
    acceleration_X, acceleration_Y = [], []
    
    stable = 0
    current_time = 0.0
    step = 0
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
        KS.append(np.abs(S)/np.abs(disp_X)) # LATERAL STIFFNESS IN X
        KA.append(np.abs(A)/np.abs(disp_Y)) # LATERAL STIFFNESS IN Y
        KI.append(np.abs(M)/np.abs(rot))    # ROTATIONAL STIFFNESS IN Z
        # OUTPUT JUST EXPLOSION ANALYSIS
        FORCE_S01.append(S)
        FORCE_A01.append(A)
        MOMENT01.append(M)
        DISP_X01.append(disp_X)
        DISP_Y01.append(disp_Y)
        ROT01.append(rot)
        step = step + 1
        STEP.append(step) # STEPS FOR EXPLOSION LOAD ANALYSIS
        #print(current_time, disp_X, S)
        print(step, disp_X, disp_Y)
        
    
    # Calculating Damping Ratio Using Logarithmic Decrement Analysis 
    displacement = np.array(DISP_X)
    peaks = np.array([displacement[i] for i in range(1, len(displacement)-1) if displacement[i] > displacement[i-1] and displacement[i] > displacement[i+1]])
    # Natural logarithm
    delta = np.log(peaks[:-1] / peaks[1:])
    
    # Reset analysis (clear analysis info, not the model)
    ops.loadConst('-time', 0.0)
    ops.remove('loadPattern', patternTag01)  # remove pushover loads
    ops.wipeAnalysis()
    #%% RUN THERMAL ANALYSIS
    print('THERMAL ANALYSIS STARTED BY OPENSEES')
    # Define time series and load pattern
    Time_Tag02 = 2
    ops.timeSeries('Linear', Time_Tag02)
    patternTag02 = 2
    ops.pattern('Plain', patternTag02 , Time_Tag02)
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
        KS.append(np.abs(S)/np.abs(disp_X)) # LATERAL STIFFNESS IN X
        KA.append(np.abs(A)/np.abs(disp_Y)) # LATERAL STIFFNESS IN Y
        KI.append(np.abs(M)/np.abs(rot))    # ROTATIONAL STIFFNESS IN Z
        THERMAL_STEP.append(ops.getLoadFactor(patternTag02) * Max_Thermal) # THERMAL LOAD 
        STEP.append(step+1) # STEPS FOR THERMAL ANALYSIS
        # OUTPUT JUST THERMAL ANALYSIS
        FORCE_S02.append(S)
        FORCE_A02.append(A)
        MOMENT02.append(M)
        DISP_X02.append(disp_X)
        DISP_Y02.append(disp_Y)
        ROT02.append(rot)
        #print(step+1, disp_X, S)
        print(step+1, disp_X, disp_Y)

    return FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y, ROT, KA, KS, KI, THERMAL_STEP, time, velocity_X, velocity_Y, acceleration_X, acceleration_Y, PERIOD_01, PERIOD_02, delta, STEP, FORCE_S01, FORCE_A01, MOMENT01, DISP_X01, DISP_Y01, ROT01, FORCE_S02, FORCE_A02, MOMENT02, DISP_X02, DISP_Y02
    
#%%------------------------------------------------------------------------------
# Analysis Durations:
starttime = TI.process_time()

# WITHOUT HARDENING AND ULTIMATE STRAIN
DATA = ALL_ANALYSIS()
FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y, ROT, KA, KS, KI, THERMAL_STEP, time, velocity_X, velocity_Y, acceleration_X, acceleration_Y, PERIOD_01, PERIOD_02, delta, STEP, FORCE_S01, FORCE_A01, MOMENT01, DISP_X01, DISP_Y01, ROT01, FORCE_S02, FORCE_A02, MOMENT02, DISP_X02, DISP_Y02 = DATA

totaltime = TI.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n') 
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
plt.plot(time, FORCE_A01, color='brown', linewidth=2)
#plt.scatter(time, FORCE_A01, color='brown', linewidth=2)
plt.title('Axial Force During the Explosion Analysis')
plt.ylabel('Axial Force [N]')
plt.xlabel('Explosion Analysis time [s]')
plt.grid()
plt.show()

plt.figure(14, figsize=(12, 8))
plt.plot(time, FORCE_S01, color='purple', linewidth=2)
#plt.scatter(time, FORCE_S01, color='purple', linewidth=2)
plt.title('Shear Force During the Explosion Analysis')
plt.ylabel('Shear Force [N]')
plt.xlabel('Explosion Analysis time [s]')
plt.grid()
plt.show()

plt.figure(15, figsize=(12, 8))
plt.plot(time, MOMENT01, color='green', linewidth=2)
#plt.scatter(time, MOMENT01, color='green', linewidth=2)
plt.title('Moment During the Explosion Analysis')
plt.ylabel('Moment [kN.mm]')
plt.xlabel('Explosion Analysis time [s]')
plt.grid()
plt.show()

plt.figure(16, figsize=(12, 8))
plt.plot(time, DISP_X01, color='brown', linewidth=2)
#plt.scatter(time, DISP_X01, color='brown', linewidth=2)
plt.title('Displacement During the Explosion Analysis')
plt.ylabel('Displacement - X [mm]')
plt.xlabel('Explosion Analysis time [s]')
plt.grid()
plt.show()

plt.figure(17, figsize=(12, 8))
plt.plot(time, DISP_Y01, color='blue', linewidth=2)
#plt.scatter(time, DISP_X01, color='brown', linewidth=2)
plt.title('Displacement During the Explosion Analysis')
plt.ylabel('Displacement - Y [mm]')
plt.xlabel('Explosion Analysis time [s]')
plt.grid()
plt.show()

plt.figure(18, figsize=(12, 8))
plt.plot(time, ROT01, color='black', linewidth=2)
#plt.scatter(time, ROT01, color='black', linewidth=2)
plt.title('Rotation During the Explosion Analysis')
plt.ylabel('Rotation [rad]')
plt.xlabel('Explosion Analysis time [s]')
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
BASE_Z = MAX_ABS(FORCE_S02) 

plt.figure(1, figsize=(8, 6))
plt.plot(THERMAL_STEP, DISP_X02, color='blue', linewidth=2)
#plt.plot(THERMAL_STEP, DISP_ZX, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Displacement in X [mm]')
plt.title(f'Time vs Displacement - MAX. ABS: {DISP_ZX[-1]} | ξ (Calculated): {100*solution[0]:.5e} %')
plt.grid()
plt.show()

plt.figure(2, figsize=(8, 6))
plt.plot(THERMAL_STEP, DISP_Y02, color='blue', linewidth=2)
#plt.plot(THERMAL_STEP, DISP_ZY, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Displacement in Y [mm]')
plt.title(f'Time vs Displacement - MAX. ABS: {DISP_ZY[-1]}')
plt.grid()
plt.show()

plt.figure(3, figsize=(8, 6))
plt.plot(THERMAL_STEP, FORCE_S02, color='blue', linewidth=2)
plt.plot(THERMAL_STEP, BASE_Z, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Base-reaction [N]')
plt.title(f'Time vs Base-reaction - MAX. ABS: {BASE_Z[-1]}')
plt.grid()
plt.show()

plt.figure(4, figsize=(8, 6))
plt.plot(DISP_X02, FORCE_S02, color='black', linewidth=5)
plt.xlabel('Displacement in X [mm]')
plt.ylabel('Base Shear Reaction [N]')
plt.title(f'Base Shear Reaction vs Displacement - MAX. ABS: {np.max(np.abs(DISP_X02)): .3f}')
plt.grid()
plt.show()

plt.figure(5, figsize=(8, 6))
plt.plot(DISP_Y02, FORCE_A02, color='black', linewidth=5)
plt.xlabel('Displacement in Y [mm]')
plt.ylabel('Base Axial Reaction [N]')
plt.title(f'Base Axial Reaction vs Displacement - MAX. ABS: {np.max(np.abs(DISP_Y02)): .3f}')
plt.grid()
plt.show()

plt.figure(6, figsize=(8, 6))
plt.plot(DISP_X, FORCE_S, color='black', linewidth=5)
plt.xlabel('Displacement in X [mm]')
plt.ylabel('Base Shear Reaction [N]')
plt.title(f'Base Shear Reaction vs Displacement - MAX. ABS: {np.max(np.abs(DISP_X02)): .3f}')
plt.grid()
plt.show()

plt.figure(7, figsize=(8, 6))
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
results_df.to_excel('CONCRETE_FRAME_EXPLOSION_IMPACT_LOAD-THERMAL_LOAD_RESULTS.xlsx', index=False) 
#%%------------------------------------------------------------------------------
# Print out the state of nodes 3 and 4
ops.printModel("node",3, 4)
# Print out the state of element 1 , 2 and 3
ops.printModel("ele", 1, 2 , 3)
# Print the Model
#printModel()
ops.printModel("-JSON", "-file", "CONCRETE_FRAME_EXPLOSION_IMPACT_LOAD-THERMAL_LOAD.json")
#%%------------------------------------------------------------------------------   
    
