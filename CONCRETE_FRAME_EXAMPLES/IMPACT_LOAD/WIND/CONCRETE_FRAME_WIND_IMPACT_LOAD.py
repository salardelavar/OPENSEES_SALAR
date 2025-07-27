######################################################################################
#                                   IN THE NAME OF ALLAH                             #
#                     WIND IMPACT LOADING ANALYSIS OF CONCRETE FRAME                 # 
#      EVALUATING STRAIN HARDENING AND ULTIMATE STRAIN CRITERIA USING OPENSEES       #
#------------------------------------------------------------------------------------#
#            THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)              #
#                       EMAIL: salar.d.ghashghaei@gmail.com                          #
######################################################################################
"""
--------------------------------------------------------
 Analysis of Wind Impact Loading on Concrete Frame
--------------------------------------------------------

This program performs a nonlinear dynamic analysis of a concrete frame subjected to Wind impact loading using OpenSees.
 The analysis evaluates two different steel material models:
1. Steel01 (without hardening and ultimate strain criteria)
2. Hysteretic (with hardening and ultimate strain criteria)

Key Features:

[1] Material Definitions
Concrete: 
  - Confined core concrete (higher strength and ductility)
  - Unconfined cover concrete
Steel Reinforcement:
  - Two material models compared
  - Includes yield strength, hardening, and ultimate strain parameters

[2] Structural Modeling
- 2D frame with columns and beam
- Nonlinear beam-column elements
- Corotational geometric transformation for large displacements

[3] Loading
- Generate a combined time series of wind pressure including gusts and turbulence for Wind pressure-time history
- Dynamic analysis with Rayleigh damping

[4] Analysis Features
- Transient analysis with HHT integrator
- Convergence criteria and iteration controls
- Extensive result recording (forces, displacements, stiffness)

[5] Post-Processing
- Multiple comparison plots between the two material models
- P-M interaction diagrams
- Force-displacement relationships
- Time-history responses
- Frame deformation visualization

[5] Advanced Calculations
- Natural period estimation
- Damping ratio calculation using logarithmic decrement
- Cumulative maximum response values

The analysis provides valuable insights into how strain hardening and ultimate strain criteria affect
 the frame's response to wind loading, which is crucial for wind-resistant design.

The results are exported to Excel for further processing and include comprehensive visualizations of all key response parameters.
"""
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

LENGTH_COL = 3000        # [mm] Column Length 
LENGTH_BM = 7000         # [mm] Beam Length 

GMfact = 9810    # standard acceleration of gravity or standard acceleration
DR = 0.05        # Intial Guess for Damping ratio
duration = 80.0  # [s] Total simulation duration
impact_duration = 5.0  # [s] Total duration of wind impact loading 
dt = 0.01        # [s] Time step
MASS = 12000     # [kg] Mass of the structure

# Define Analysis Properties
MAX_ITERATIONS = 10000     # Convergence iteration for test
MAX_TOLERANCE = 1.0e-8     # Convergence tolerance for test
#STEEL_KIND: 1 -> WITHOUT HARDENING AND ULTIMATE STRAIN
#STEEL_KIND: 2 -> WITH HARDENING AND ULTIMATE STRAIN
#%%------------------------------------------------------------------------------
def generate_wind_pressure_time_series(
    impact_duration, dt, wave_type='sin+cos', frequency_factor=1.0, amplitude=500.0, turbulence_intensity=100.0):
    """
    Generate a combined time series of wind pressure including gusts and turbulence.

    Parameters:
    - impact_duration: Duration of the impact event (seconds).
    - dt: Time step for the time series (seconds).
    - wave_type: Type of waveform ('sine', 'cos', 'sin+cos'). Default is 'sin+cos'.
    - frequency_factor: Factor to modify the frequency of the waveform. Default is 1.0.
    - amplitude: Amplitude of the waveform. Default is 500.0.
    - turbulence_intensity: Intensity of the random turbulence. Default is 100.0.

    Returns:
    - time: Array of time values.
    - pressure: Array of composite pressure values.
    """
    time = np.arange(0, impact_duration, dt)
    frequency = frequency_factor / impact_duration  # Adjust frequency based on the impact duration

    # Base waveform
    if wave_type == 'sine':
        base_pressure = np.sin(2 * np.pi * frequency * time) * amplitude
    elif wave_type == 'cos':
        base_pressure = np.cos(2 * np.pi * frequency * time) * amplitude
    elif wave_type == 'sin+cos':
        base_pressure = (np.sin(2 * np.pi * frequency * time) + np.cos(1.5 * np.pi * frequency * time)) * amplitude
    else:
        raise ValueError("Unsupported wave type. Choose from 'sine', 'cos', or 'sin+cos'.")

    # Gust-induced load
    gust_pressure = 0.8 * amplitude * np.sin(0.3 * np.pi * time / impact_duration)

    # Random turbulent fluctuations
    turbulent_pressure = np.random.normal(0, turbulence_intensity, size=time.shape)

    # Better simulation of wind pressure with exponential decay factor for gusts
    decay_factor = np.exp(-0.5 * time / impact_duration)
    gust_pressure_with_decay = gust_pressure * decay_factor

    # Combine all loads
    total_pressure = base_pressure + gust_pressure_with_decay + turbulent_pressure

    return time, total_pressure

# Generate composite loading
time_series, wind_pressure = generate_wind_pressure_time_series(impact_duration, dt, wave_type='sin+cos', frequency_factor=0.5, amplitude=500.0, turbulence_intensity=200.0)

# Plot the composite loading
def plot_wind_pressure_time_series(time, pressure):
    # Plot the pressure time history
    plt.figure(figsize=(10, 6))
    plt.plot(time, pressure, label='Wind Pressure', color='purple')
    plt.xlabel('Time (s)')
    plt.ylabel('Pressure Force (N)')
    plt.title('Wind Pressure Loading')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

plot_wind_pressure_time_series(time_series, wind_pressure)

# IMPACT ANALYSIS FUNCTION
def IMPACT_ANALYSIS(LENGTH_COL, LENGTH_BM, STEEL_KIND):
    # Initialize OpenSees model
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)

    # Define nodes
    ops.node(1, 0.0, 0.0)
    ops.node(2, LENGTH_BM, 0.0)
    ops.node(3, 0.0, LENGTH_COL)
    ops.node(4, LENGTH_BM, LENGTH_COL)
    # Define boundary conditions
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 1, 1, 1)
    # Define mass
    ops.mass(3, MASS, MASS, 0.0)
    ops.mass(4, MASS, MASS, 0.0)

    secTagC = 10
    secTagB = 20
    coreTag = 1
    coverTag = 2
    steelTag = 3
    steelPlateTag = 4
    numBarsTop, barAreaTop = 5, np.pi *(18**2)/4
    numBarsBot, barAreaBot = 5, np.pi *(20**2)/4
    numBarsIntTot, barAreaInt = 4, np.pi *(5**2)/4
    
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
    ops.element('nonlinearBeamColumn', 1, 1, 3, numIntgrPts, secTagC, transfTag) # COLUMN 01
    ops.element('nonlinearBeamColumn', 2, 2, 4, numIntgrPts, secTagC, transfTag) # COLUMN 02
    ops.element('nonlinearBeamColumn', 3, 3, 4, numIntgrPts, secTagB, transfTag) # BEAM 01

    # Data storage
    FORCE_S, FORCE_A, MOMENT = [], [], []
    DISP_X, DISP_Y, ROT = [], [], []
    KA, KS, KI, STEP = [], [], [], []
    
    # Static analysis
    time_series_tag = 1
    pattern_tag = 1
    # Apply time-dependent wind loading
    ops.timeSeries('Path', time_series_tag, '-dt', dt, '-values', *wind_pressure)
    ops.pattern('Plain', pattern_tag, time_series_tag)
    ops.load(3, 1.0, 0.0, 0.0)
    #ops.load(4, 1.0, 0.0, 0.0)
    print('Impact Load Defined Done.')

    
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
        KS.append(np.abs(S)/np.abs(disp_X)) # LATERAL STIFFNESS IN X
        KA.append(np.abs(A)/np.abs(disp_Y)) # LATERAL STIFFNESS IN Y
        KI.append(np.abs(M)/np.abs(rot))    # ROTATIONAL STIFFNESS IN Z
        print(current_time, disp_X, S)
    
    # Calculating Damping Ratio Using Logarithmic Decrement Analysis 
    displacement = np.array(DISP_X)
    peaks = np.array([displacement[i] for i in range(1, len(displacement)-1) if displacement[i] > displacement[i-1] and displacement[i] > displacement[i+1]])
    # Natural logarithm
    delta = np.log(peaks[:-1] / peaks[1:])
        
    #ops.wipe()    

    return FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y, ROT, KA, KS, KI, time, velocity_X, velocity_Y, acceleration_X, acceleration_Y, PERIOD_01, PERIOD_02, delta

#%%------------------------------------------------------------------------------
# Analysis Durations:
starttime = TI.process_time()

# WITHOUT HARDENING AND ULTIMATE STRAIN
DATA = IMPACT_ANALYSIS(LENGTH_COL, LENGTH_BM, STEEL_KIND=1)
FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y, ROT, KA, KS, KI, time, velocity_X, velocity_Y, acceleration_X, acceleration_Y, PERIOD_01, PERIOD_02, delta = DATA
# WITH HARDENING AND ULTIMATE STRAIN
DATA02 = IMPACT_ANALYSIS(LENGTH_COL, LENGTH_BM, STEEL_KIND=2)
FORCE_S02, FORCE_A02, MOMENT02, DISP_X02, DISP_Y02, ROT02, KA02, KS02, KI02, time02, velocity_X02, velocity_Y02, acceleration_X02, acceleration_Y02, PERIOD_012, PERIOD_022, delta02 = DATA02

print(f"WITHOUT HARDENING AND ULTIMATE STRAIN: \n Period 01: {PERIOD_01:.4e}  - Period 02: {PERIOD_02:.4e}")
print(f"WITH HARDENING AND ULTIMATE STRAIN: \n Period 01: {PERIOD_012:.4e}  - Period 02: {PERIOD_022:.4e}")

totaltime = TI.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n') 
#%%------------------------------------------------------------------------------
plt.figure(1, figsize=(12, 8))
plt.plot(MOMENT, FORCE_A, color='black')
plt.plot(MOMENT02, FORCE_A02, color='cyan', linestyle='--', linewidth=2)
#plt.scatter(MOMENT, FORCE_A, color='blue', linewidth=2)
#plt.scatter(MOMENT02, FORCE_A02, color='cyan', linestyle='--', linewidth=2)
plt.title('P-M Interaction')
plt.ylabel('Axial Force [N]')
plt.xlabel('Bending Moment [N.mm]')
plt.legend(['Steel01: WITHOUT HARDENING AND ULTIMATE STRAIN', 'Hysteretic: WITH HARDENING AND ULTIMATE STRAIN'])
plt.grid()
plt.show()

plt.figure(2, figsize=(12, 8))
plt.plot(DISP_X, FORCE_S, color='green', linewidth=2)
plt.plot(DISP_X02, FORCE_S02, color='lime', linestyle='--', linewidth=2)
#plt.scatter(DISP_X, FORCE_S, color='purple', linewidth=2)
#plt.scatter(DISP_X02, FORCE_S02, color='magenta', linestyle='--', linewidth=2)
plt.title('SHEAR FORCE-DISPLACEMENT DIAGRAM')
plt.ylabel('Shear Force [N]')
plt.xlabel('Displacement  in X [mm]')
plt.legend(['Steel01: WITHOUT HARDENING AND ULTIMATE STRAIN', 'Hysteretic: WITH HARDENING AND ULTIMATE STRAIN'])
plt.grid()
plt.show()

plt.figure(3, figsize=(12, 8))
plt.plot(DISP_Y, FORCE_A, color='purple', linewidth=2)
plt.plot(DISP_Y02, FORCE_A02, color='magenta', linestyle='--', linewidth=2)
#plt.scatter(DISP_Y, FORCE_A, color='purple', linewidth=2)
#plt.scatter(DISP_Y02, FORCE_A02, color='magenta', linestyle='--', linewidth=2)
plt.title('AXIAL FORCE-DISPLACEMENT DIAGRAM')
plt.ylabel('Axial Force [N]')
plt.xlabel('Displacement in Y [mm]')
plt.legend(['Steel01: WITHOUT HARDENING AND ULTIMATE STRAIN', 'Hysteretic: WITH HARDENING AND ULTIMATE STRAIN'])
plt.grid()
plt.show()

plt.figure(4, figsize=(12, 8))
plt.plot(ROT, MOMENT, color='red', linewidth=2)
plt.plot(ROT02, MOMENT02, color='orange', linestyle='--', linewidth=2)
#plt.scatter(ROT, MOMENT, color='red', linewidth=2)
#plt.scatter(ROT02, MOMENT02, color='orange', linestyle='--', linewidth=2)
plt.title('MOMENT-ROTATION DIAGRAM')
plt.ylabel('Moment [kN.mm]')
plt.xlabel('Rotation [rad]')
plt.legend(['Steel01: WITHOUT HARDENING AND ULTIMATE STRAIN', 'Hysteretic: WITH HARDENING AND ULTIMATE STRAIN'])
plt.grid()
plt.show()

plt.figure(5, figsize=(12, 8))
#plt.plot(KI, KS, color='black', linewidth=2)
#plt.plot(KI02, KS02, color='grey', linestyle='--', linewidth=2)
plt.scatter(KI, KS, color='black', linewidth=2)
plt.scatter(KI02, KS02, color='grey', linestyle='--', linewidth=2)
plt.title('ROTATIONAL STIFFNESS-LATERAL STIFFNESS DIAGRAM')
plt.ylabel('Rotational Stiffness [N.mm/Rad]')
plt.xlabel('Lateral Stiffness in X Dir. [N/mm]')
plt.legend(['Steel01: WITHOUT HARDENING AND ULTIMATE STRAIN', 'Hysteretic: WITH HARDENING AND ULTIMATE STRAIN'])
plt.semilogx()
plt.semilogy()
plt.grid()
plt.show()

plt.figure(6, figsize=(12, 8))
#plt.plot(KI, KA, color='black', linewidth=2)
#plt.plot(KI02, KA02, color='grey', linestyle='--', linewidth=2)
plt.scatter(KI, KA, color='black', linewidth=2)
plt.scatter(KI02, KA02, color='grey', linestyle='--', linewidth=2)
plt.title('ROTATIONAL STIFFNESS-LATERAL STIFFNESS DIAGRAM')
plt.ylabel('Rotational Stiffness [N.mm/Rad]')
plt.xlabel('Lateral Stiffness in Y Dir. [N/mm]')
plt.legend(['Steel01: WITHOUT HARDENING AND ULTIMATE STRAIN', 'Hysteretic: WITH HARDENING AND ULTIMATE STRAIN'])
plt.semilogx()
plt.semilogy()
plt.grid()
plt.show()

plt.figure(7, figsize=(12, 8))
plt.plot(time, FORCE_A, color='brown', linewidth=2)
plt.plot(time02, FORCE_A02, color='gold', linestyle='--', linewidth=2)
#plt.scatter(time, FORCE_A, color='brown', linewidth=2)
#plt.scatter(time02, FORCE_A02, color='gold', linestyle='--', linewidth=2)
plt.title('Axial Force During the Analysis')
plt.ylabel('Axial Force [N]')
plt.xlabel('Times')
plt.legend(['Steel01: WITHOUT HARDENING AND ULTIMATE STRAIN', 'Hysteretic: WITH HARDENING AND ULTIMATE STRAIN'])
plt.grid()
plt.show()

plt.figure(8, figsize=(12, 8))
plt.plot(time, FORCE_S, color='purple', linewidth=2)
plt.plot(time02, FORCE_S02, color='#BF77F6', linestyle='--', linewidth=2)
#plt.scatter(time, FORCE_S, color='purple', linewidth=2)
#plt.scatter(time02, FORCE_S02, color='#BF77F6', linestyle='--', linewidth=2)
plt.title('Shear Force During the Analysis')
plt.ylabel('Shear Force [N]')
plt.xlabel('Times')
plt.legend(['Steel01: WITHOUT HARDENING AND ULTIMATE STRAIN', 'Hysteretic: WITH HARDENING AND ULTIMATE STRAIN'])
plt.grid()
plt.show()

plt.figure(9, figsize=(12, 8))
plt.plot(time, MOMENT, color='green', linewidth=2)
plt.plot(time02, MOMENT02, color='lime', linestyle='--', linewidth=2)
#plt.scatter(time, MOMENT, color='green', linewidth=2)
#plt.scatter(time02, MOMENT02, color='lime', linestyle='--', linewidth=2)
plt.title('Moment During the Analysis')
plt.ylabel('Moment [kN.mm]')
plt.xlabel('Times')
plt.legend(['Steel01: WITHOUT HARDENING AND ULTIMATE STRAIN', 'Hysteretic: WITH HARDENING AND ULTIMATE STRAIN'])
plt.grid()
plt.show()

plt.figure(10, figsize=(12, 8))
plt.plot(time, DISP_X, color='brown', linewidth=2)
plt.plot(time02, DISP_X02, color='gold', linestyle='--', linewidth=2)
#plt.scatter(time, DISP_X, color='brown', linewidth=2)
#plt.scatter(time02, DISP_X02, color='gold', linestyle='--', linewidth=2)
plt.title('Displacement During the Analysis')
plt.ylabel('Displacement - X [mm]')
plt.xlabel('Times')
plt.legend(['Steel01: WITHOUT HARDENING AND ULTIMATE STRAIN', 'Hysteretic: WITH HARDENING AND ULTIMATE STRAIN'])
plt.grid()
plt.show()

plt.figure(11, figsize=(12, 8))
plt.plot(time, DISP_Y, color='blue', linewidth=2)
plt.plot(time02, DISP_Y02, color='cyan', linestyle='--', linewidth=2)
#plt.scatter(time, DISP_X, color='brown', linewidth=2)
#plt.scatter(time02, DISP_X02, color='gold', linestyle='--', linewidth=2)
plt.title('Displacement During the Analysis')
plt.ylabel('Displacement - Y [mm]')
plt.xlabel('Times')
plt.legend(['Steel01: WITHOUT HARDENING AND ULTIMATE STRAIN', 'Hysteretic: WITH HARDENING AND ULTIMATE STRAIN'])
plt.grid()
plt.show()

plt.figure(12, figsize=(12, 8))
plt.plot(time, ROT, color='black', linewidth=2)
plt.plot(time02, ROT02, color='grey', linestyle='--', linewidth=2)
#plt.scatter(time, ROT, color='black', linewidth=2)
#plt.scatter(time02, ROT02, color='grey', linestyle='--', linewidth=2)
plt.title('Rotation During the Analysis')
plt.ylabel('Rotation [rad]')
plt.xlabel('Times')
plt.legend(['Steel01: WITHOUT HARDENING AND ULTIMATE STRAIN', 'Hysteretic: WITH HARDENING AND ULTIMATE STRAIN'])
plt.grid()
plt.show()
#%%------------------------------------------------------------------------------  
# EXACT SOLUTION:
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
plt.title(f'Time vs Displacement - MAX. ABS: {DISP_ZX[-1]} | ξ (Calculated): {100*solution[0]:.5e} %')
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

#%%------------------------------------------------------------------------------  
# %% Plot 2D Frame Shapes
S04.PLOT_2D_FRAME(deformed_scale=100000)  # Adjust scale factor as needed
#%%------------------------------------------------------------------------------   
# EXPORT DATA TO EXCEL
import pandas as pd
DATA_TOTAL = { 
    'DISP_X_WO': DISP_X,# WITHOUT HARDENING AND ULTIMATE STRAIN
    'DISP_X_W': DISP_X02,# WITH HARDENING AND ULTIMATE STRAIN
    'DISP_Y_WO': DISP_Y,
    'DISP_Y_W': DISP_Y02,
    'ROTATION_WO': ROT,
    'ROTATION_W': ROT02,
    'VELO_WO': velocity_X,
    'VELO_W': velocity_X02,
    'ACCEL_WO': acceleration_X,
    'ACCEL_W': acceleration_X02,
    'AXIAL_FORCE_WO': FORCE_A,
    'AXIAL_FORCE_W': FORCE_A02,
    'SHEAR_FORCE_WO': FORCE_S,
    'SHEAR_FORCE_W': FORCE_S02,
    'MOMENT_WO': MOMENT,
    'MOMENT_W': MOMENT02,
    'AXIAL_RIGIDITY_WO': np.abs(FORCE_A),
    'AXIAL_RIGIDITY_W': np.abs(FORCE_A02),
    'ROTATIONAL_ST_WO': KI,
    'ROTATIONAL_ST_W': KI02,
    'LATERAL_ST_Y_WO': KA,
    'LATERAL_ST_Y_W': KA02,
    'LATERAL_ST_X_WO': KS,
    'LATERAL_ST_X_W': KS02,
}
# Convert to DataFrame
results_df = pd.DataFrame(DATA_TOTAL)
# Export the DataFrame to an Excel file
results_df.to_excel('CONCRETE_FRAME_IMPACT_LOAD_WIND_RESULTS.xlsx', index=False) 
#%%------------------------------------------------------------------------------

    
