######################################################################################################################
#                          >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                          #
#                                  MODELING OF PENDULUM MDOF STRUCTURE USING OPENSEES                                #
#                                                   P(t) = -m * ug**(t)                                              #
#--------------------------------------------------------------------------------------------------------------------#
#                   EVALUATION OF DAMPING FORCE (fD), SPRING FORCE (fS) AND INERTIA FORCE (fI)                       #
#--------------------------------------------------------------------------------------------------------------------#
#                              THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                            #
#                                       EMAIL: salar.d.ghashghaei@gmail.com                                          #
######################################################################################################################
"""
Performs time-dependent loading analysis of a Multi Degree of Freedom (MDOF)
 structure using OpenSeesPy, comparing elastic and inelastic spring behavior. 
 Key features include:

1. Implements both elastic (linear) and hysteretic (nonlinear) material models for
 structural springs.
2. Supports initial conditions for displacement, velocity, and acceleration.
3. Uses Newmark's method for time integration with Newton-Raphson iteration.
4. Calculates damping ratios using logarithmic decrement from response peaks.
5. Generates force-displacement backbone curves for inelastic material.
6. Tracks and plots time-history responses (displacement, velocity, acceleration, reactions).
7. Compares elastic vs inelastic system performance.
8. Includes convergence checks and analysis stability monitoring.
9. Outputs model data in JSON format for post-processing.
10. Provides theoretical validation through natural frequency calculations.

Particularly useful for earthquake engineering applications, 
allowing evaluation of structural response under time-dependent loading
 with different material nonlinearities and damping characteristics.
 The hysteretic material model captures energy dissipation 
 inelastic deformation, while the elastic case serves as a reference for linear behavior.
----------------------------------------------------------------------------
This code performs a comparative nonlinear dynamic analysis of a multi-degree-of-freedom (MDOF) system
 using OpenSeesPy, representing a simplified structure (e.g., a pendulum or column) with
 a concentrated mass and a zero-length spring at the base. The structural behavior is
 defined through a bilinear moment-rotation relationship with post‑yield hardening and
 degradation (defined by `pos_disp`/`pos_Moment` and negative counterparts), allowing
 both elastic (`Elastic` material) and inelastic (`HystereticSM` material) spring models.
 The system is subjected to three types of time‑dependent loads:
     (1) a sinusoidal excitation,
     (2) a damped sinusoidal excitation, and 
     (3) an artificial ground acceleration record generated to match a target spectrum (max 0.2g).
Damping is introduced via a viscous material (coefficient `C` computed from the elastic period and
updated each time step using the current instantaneous period). 
Transient analysis is performed using the Newmark‑Beta integrator with Newton‑Raphson
 iterations. For both elastic and inelastic cases, the code records displacements,
 velocities, accelerations, base reactions, and instantaneous stiffness/period.
 It also computes the damping ratio from the displacement decay using logarithmic 
 decrement and performs eigenvalue analyses at each step to track period evolution.
 The results are visualized extensively: time histories, force‑displacement hysteresis,
 damping‑force vs. velocity, inertia‑force vs. acceleration, and the evolution of stiffness and period. 
 Finally, outputs are saved to an Excel file for further processing. 
 The comparison highlights how inelastic behavior (yielding and stiffness degradation)
 alters the dynamic response, period elongation, and energy dissipation compared to a 
 purely elastic system.
"""
#%%-----------------------------------------------------------------------------
# YOUUBE: Simple Pendulum
'https://www.youtube.com/watch?v=fnvGVsxPuLs'
# YOUTUBE: Everything You Need To Know About Pendulums: Physics Help Room
'https://www.youtube.com/watch?v=0q0L7Fj4dk8'
# YOUTUBE: How a Giant Pendulum Made Taipei101 Possible
'https://www.youtube.com/watch?v=mGe9zjwK2gQ'
# BOOK: Differential Equations for Engineers - Wei-Chau Xie - CAMBRIDGE
'https://www.cambridge.org/core/books/differential-equations-for-engineers/1B8F1A62BF6F98EB972CFE9114FA8B84'
#%%-----------------------------------------------------------------------------
import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
import time as TI
import ANALYSIS_FUNCTION as S01
import PERIOD_FUN as S02
import DAMPING_RATIO_FUN as S05
import EIGENVALUE_ANALYSIS_FUN as S03
import GENERATE_ARTIFICIAL_ACCEL_FUN as S04
import SALAR_MATH as S06
#%%------------------------------------------------------------------------------------------------
# Define parameters (units: m, N)
FY = 85000.0                      # [N] Yield Moment of Structure
FU = 1.1 * FY                     # [N] Ultimate Moment of Structure
Ke = 4500000.0                    # [N/m] Spring Elastic Stiffness
DY = FY / Ke                      # [m] Yield Rotation
DSU = 0.36                        # [m] Ultimate Rotation
Ksh = (FU - FY) / (DSU - DY)      # [N/m] Rotation Hardening Modulus
Kp = FU / DSU                     # [N/m] Spring Plastic Stiffness
b = Ksh / Ke                      # Rotation Hardening Ratio

M = 5000.0                        # [kg] Mass
zi = 0.05                         # Damping ratio
duration = 20.0                   # [s] Analysis duration
dt = 0.001                        # [s] Time step
#%%------------------------------------------------------------------------------------------------
# DEFINE ANALYSIS PROPERTIES
MAX_ITERATIONS = 20000    # Convergence iteration for test
MAX_TOLERANCE = 1.0e-6   # Convergence tolerance for test
#SPRING_KIND: 1 -> 'ELASTIC'
#SPRING_KIND: 2 -> 'INELASTIC'

#%%------------------------------------------------------------------------------------------------
# Positive branch points
pos_disp = [0, DY, DSU, 1.1*DSU, 1.25*DSU]
pos_Moment = [0, FY, FU, 0.2*FU, 0.1*FU]
KP = np.array([FY, DY, FU, DSU, 0.2*FU, 1.1*DSU, 0.1*FU, 1.25*DSU])

# Negative branch points
neg_disp = [0, -DY, -DSU, -1.1*DSU, -1.25*DSU]
neg_Moment = [0, -FY, -FU, -0.2*FU, -0.1*FU]
KN = np.array([-FY, -DY, -FU, -DSU, -0.2*FU, -1.1*DSU, -0.1*FU, -1.25*DSU])
#%%------------------------------------------------------------------------------------------------
# Plot
plt.plot(pos_disp, pos_Moment, marker='o', color='red')
plt.plot(neg_disp, neg_Moment, marker='o', color='black')

plt.xlabel("Rotation [m]")
plt.ylabel("Moment [N]")
plt.title("Moment–Rotation Curve")
plt.grid(True)
plt.axhline(0, linewidth=0.5)
plt.axvline(0, linewidth=0.5)
plt.show()
# ELASIC PERIOD:
ELAS_PERIOD = 2*np.pi * np.sqrt(M/Ke)
print(f'ELASIC PERIOD: {ELAS_PERIOD:.3f} (s)')     
# PLASIC PERIOD:
PLAS_PERIOD = 2*np.pi * np.sqrt(M/Kp)  
print(f'PLASIC PERIOD: {PLAS_PERIOD:.3f} (s)') 

# INITIAL MASS FOR RESPONSE SPECTRUM ANALYSIS
mi = (PLAS_PERIOD/2*np.pi)**2 * Kp 
print(mi)
#%%------------------------------------------------------------------------------------------------   
omega = np.sqrt(Ke/M)    # [N/m] Elastic Stiffness
C = 2 * zi * omega * M   # Damping Coefficient
#%% DEFINE EXTERNAL TIME-DEPENDENT LOADING PROPERTIES
DT = dt                   # [s] Time step
DT_time = 5.0             # [s] Total external Load Analysis Durations [*******]
omega_DT = 5.0715         # [rad/s] Natural angular frequency

def EXTERNAL_TIME_DEPENDENT(m, omega_DT, DT, DT_time): # P(t) = -m * ug**(t)
    import numpy as np
    import matplotlib.pyplot as plt
    # External Load Durations
    num_steps = int(DT_time / DT)
    load_time = np.linspace(0, DT_time, num_steps) 
    target_frequency = 1.0 * omega_DT  # Target excitation frequency
    DT_load = m * np.sin(target_frequency * load_time)
    # Plot External Time-dependent Loading
    plt.figure(figsize=(10, 6))
    plt.plot(load_time, DT_load, label=f'External Loading - Max: {np.max(DT_load):.3f}', linewidth=5)
    plt.title('External Time-dependent Loading Over Time [P(t) = -m * ug**(t)]')
    plt.xlabel('Time (s)')
    plt.ylabel('Force (N)')
    plt.grid(True)
    plt.legend()
    plt.show()
    return DT_load

def EXTERNAL_TIME_DEPENDENT_02(m, omega_DT, DT, DT_time): # P(t) = -m * ug**(t)
    import numpy as np
    import matplotlib.pyplot as plt
    # External Load Durations
    num_steps = int(DT_time / DT)
    load_time = np.linspace(0, DT_time, num_steps) 
    target_frequency = 1.0 * omega_DT  # Target excitation frequency
    DT_load = -m * np.exp(-0.05*target_frequency * load_time) * np.sin(target_frequency * load_time)
    # Plot External Time-dependent Loading
    plt.figure(figsize=(10, 6))
    plt.plot(load_time, DT_load, label=f'External Loading - Max: {np.max(DT_load):.3f}', linewidth=5)
    plt.title('External Time-dependent Loading Over Time [P(t) = -m * ug**(t)]')
    plt.xlabel('Time (s)')
    plt.ylabel('Force (N)')
    plt.grid(True)
    plt.legend()
    plt.show()
    return DT_load

#%%-------------------------------------------------------------------------------
# Generate Artificial Acceleration Record for Response Spectrum
print("Generating artificial acceleration record...")
t_acc, acc_record = S04.GENERATE_ARTIFICIAL_ACCEL(duration=DT_time, dt=DT, max_accel=0.2*9.81)

plt.figure(figsize=(10, 4))
plt.plot(t_acc, acc_record, color='black', linewidth=1)
plt.xlabel('Time [s]')
plt.ylabel('Earthquake Ground Acceleration [m/s²]')
plt.title('Artificial Acceleration Record for Response Spectrum Analysis')
plt.grid(True)
plt.tight_layout()
plt.show()

# Plot Histogram
X = acc_record
HISTO_COLOR = 'cyan' 
LABEL = 'Generate Artificial Acceleration Record Histogram'
S06.HISROGRAM_BOXPLOT(X, HISTO_COLOR, LABEL)

# Plot Histogram of First derivative (changes)
changes = np.diff(X)
HISTO_COLOR = 'lime' 
LABEL = 'First derivative Histogram of Generate Artificial Acceleration Record'
S06.HISROGRAM_BOXPLOT(changes, HISTO_COLOR, LABEL)

# Calculate cumulative sum of absolute acceleration and normalize it
cum_sum_acc = np.cumsum(np.abs(acc_record))
#cum_sum_acc = np.cumsum(np.abs(np.sort(acc_record)))
sum_acc = np.sum(np.abs(acc_record))
Normalized = cum_sum_acc / sum_acc

plt.figure(figsize=(12, 8))
plt.plot(t_acc, Normalized, color='red', linewidth=4.5)
plt.xlabel('Time [s]')
plt.ylabel('Normalized Cumulative Value')
plt.title('Normalized Cumulative Absolute Acceleration')
plt.grid(True, alpha=0.3)
plt.axhline(y=0.1, color='blue', linestyle='--', alpha=0.7, label='10% of total energy')
plt.axhline(y=0.9, color='gray', linestyle='--', alpha=0.7, label='90% of total energy')
plt.legend()
plt.tight_layout()
plt.show()
"""
NOTICE: 
This plot shows the Normalized Cumulative Absolute Acceleration,
 which represents the accumulation of earthquake energy over time.
"""

def EXTERNAL_TIME_DEPENDENT_03(m, acc_record, DT, DT_time): # P(t) = -m * ug**(t)
    import numpy as np
    import matplotlib.pyplot as plt
    # External Load Durations
    num_steps = int(DT_time / DT)
    load_time = np.linspace(0, DT_time, num_steps) 
    DT_load = -m * acc_record
    # Plot External Time-dependent Loading
    plt.figure(figsize=(10, 6))
    plt.plot(load_time, DT_load, label=f'External Loading - Max: {np.max(DT_load):.3f}', linewidth=1)
    plt.title('External Time-dependent Loading Over Time [P(t) = -m * ug**(t)]')
    plt.xlabel('Time (s)')
    plt.ylabel('Force (N)')
    plt.grid(True)
    plt.legend()
    plt.show()
    return DT_load
#%% ---------------------------------------------
def EXTERNAL_TIME_DEPENDENT_LOAD_MDOF(SPRING_TYPE, duration, dt):
    # Create model
    ops.wipe()
    # Create a 2D model with 3 DOF per node
    ops.model('Basic', '-ndm', 2, '-ndf', 3)
    
    L = 1.0 # [m] Element Length
    
    # Add nodes
    ops.node(1, 0.0, 0.0)
    ops.node(2, 0.0, L)
    
    ops.node(11, 0.0, 0.0)#-1.0e-8
    
    # Define mass to node 2
    ops.mass(2, M, M, M)
    
    # Fix node 2
    ops.fix(2, 1, 1, 0)
    
    # Fix spring node
    ops.fix(11, 1, 1, 1)
        
    # Elastic Element
    B = 0.05                # [m] Element Section Width
    H = 0.05                # [m] Element Section Height
    AREA = B*H              # [m^2] Element Section Area
    Iz = B*(H**3) / 12
    Es = 200000.0/10**-6    # [N/m^2] Element Young's Modulus
    
    # Geometric transformation
    transfTag = 1
    #ops.geomTransf('Linear', transfTag)
    #ops.geomTransf('PDelta', transfTag)
    ops.geomTransf('Corotational', transfTag)
    ops.element('elasticBeamColumn', 100, 1, 2, AREA, Es, Iz, transfTag, '-mass', 0.0001)
    #INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/elasticBeamColumn.html
    
    # Define material and element
    if SPRING_TYPE == 'ELASTIC': 
        ops.uniaxialMaterial('Elastic', 1, Ke, C)
        ops.element('zeroLength', 10, 11, 1,'-mat', 1,'-dir', 3) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/ZeroLength.html
    if SPRING_TYPE == 'INELASTIC':
        ops.uniaxialMaterial('HystereticSM', 1, '-posEnv', *KP.flatten(), '-negEnv', *KN.flatten(), '-pinch', 1, 1)
        #INFO LINK: https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/HystereticSM.html
        ops.uniaxialMaterial('Viscous', 2, C, 1.0)  # Material for C (alpha=1.0 for linear)
        ops.element('zeroLength', 10, 11, 1,'-mat', 1, 2,'-dir', 3, 3)
        
    # Define analysis
    ops.constraints('Plain')# 'Transformation'
    ops.numberer('Plain') # 'RCM'
    ops.system('BandGeneral')#  'Umfpack'
    ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/normDispIncr.html
    #ops.integrator('CentralDifference')  # JUST FOR LINEAR ANALYSIS - INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/centralDifference.html
    alpha=0.5; beta=0.25;
    ops.integrator('Newmark', alpha, beta) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/newmark.html
    #alpha=2/3;gamma=1.5-alpha; gamma=1.5-alpha;beta=(2-alpha)**2/4;
    #ops.integrator('HHT', alpha, gamma, beta) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/hht.html
    ops.algorithm('Newton')  # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/algorithm.html
    ops.analysis('Transient') # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/analysis.html
    
    # P(t) = -m * ug**(t)
    #DT_load = EXTERNAL_TIME_DEPENDENT(M, omega_DT, DT, DT_time)
    #DT_load = EXTERNAL_TIME_DEPENDENT_02(M, omega_DT, DT, DT_time)
    DT_load = EXTERNAL_TIME_DEPENDENT_03(M, acc_record, DT, DT_time)
    
    # Static analysis
    time_series_tag = 1
    pattern_tag = 1
    # Apply time-dependent explosion loading
    ops.timeSeries('Path', time_series_tag, '-dt', dt, '-values', *DT_load)
    ops.pattern('Plain', pattern_tag, time_series_tag)
    ops.load(2, 1.0, -1.0, -1.0)

    # Perform analysis
    time = []
    dispX, dispY, dispZ = [], [], []
    velX, velY, velZ = [], [], []
    accelX, accelY, accelZ = [], [], []
    reactionX, reactionY, reactionZ = [], [], []
    stiffnessX, stiffnessY, stiffnessZ = [], [], []
    OMEGA, PERIOD = [], []
    PERIOD_MIN, PERIOD_MAX = [], []
    FDx, FSx, FIx = [], [], []
    FDy, FSy, FIy = [], [], []
    FDz, FSz, FIz = [], [], []
    
    stable = 0
    current_time = 0.0
        
    while stable == 0 and current_time < duration:
        ops.analyze(1, dt)
        S01.ANALYSIS(stable, 1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS
        current_time = ops.getTime()
        time.append(current_time)
        ops.reactions()
        reactionX.append(ops.nodeReaction(1, 1)) # BASE REACTION IN X DIR.
        reactionY.append(ops.nodeReaction(1, 2)) # BASE REACTION IN Y DIR.
        reactionZ.append(ops.nodeReaction(1, 3)) # BASE REACTION IN Z DIR.
        dispX.append(ops.nodeDisp(2, 1))         # DISPLACEMENT IN X DIR.  
        dispY.append(ops.nodeDisp(2, 2))         # DISPLACEMENT IN Y DIR.
        dispZ.append(ops.nodeDisp(2, 3))         # ROTATION IN Z DIR.
        velX.append(ops.nodeVel(2, 1))           # VELOCITY IN X DIR.
        velY.append(ops.nodeVel(2, 2))           # VELOCITY IN Y DIR.
        velZ.append(ops.nodeVel(2, 3))           # VELOCITY IN Z DIR.
        accelX.append(ops.nodeAccel(2, 1))       # ACCELERATION IN X DIR.
        accelY.append(ops.nodeAccel(2, 2))       # ACCELERATION IN Y DIR.
        accelZ.append(ops.nodeAccel(2, 3))       # ACCELERATION IN Z DIR.
        stiffnessX.append(np.abs(reactionX[-1]) / np.abs(dispX[-1]))
        stiffnessY.append(np.abs(reactionY[-1]) / np.abs(dispY[-1]))
        stiffnessZ.append(np.abs(reactionZ[-1]) / np.abs(dispZ[-1]))
        OMEGA.append(np.sqrt(stiffnessZ[-1]/M))
        PERIOD.append((np.pi * 2) / OMEGA[-1])
        # IN EACH STEP STRUCTURAL PERIOD WILL BE CALCULATED
        PERIODmin, PERIODmax = S03.EIGENVALUE_ANALYSIS(1, PLOT=True)
        PERIOD_MIN.append(PERIODmin)
        PERIOD_MAX.append(PERIODmax)
        Ci = 2 * zi * (2*np.pi/PERIODmax) * M         # Damping COEFFICIENT - UPDATED AND CHANGES IN EACH STEP
        FDx.append(Ci * velX[-1])                     # DAMPING FORCE IN X DIR.
        FSx.append(stiffnessX[-1] * dispX[-1])        # SPRING FORCE IN X DIR.
        FIx.append(M * accelX[-1])                    # INERTIA FORCE IN X DIR.
        FDy.append(Ci * velY[-1])                     # DAMPING FORCE IN Y DIR.
        FSy.append(stiffnessY[-1] * dispX[-1])        # SPRING FORCE IN Y DIR.
        FIy.append(M * accelY[-1])                    # INERTIA FORCE IN Y DIR.
        FDz.append(Ci * velZ[-1])                     # DAMPING FORCE IN Z DIR.
        FSz.append(stiffnessZ[-1] * dispZ[-1])        # SPRING FORCE IN Z DIR.
        FIz.append(M * accelZ[-1])                    # INERTIA FORCE IN Z DIR.
        print(time[-1], dispZ[-1], velZ[-1])

    # Compute modal properties
    ops.modalProperties("-print", "-file", "SALAR_ModalReport.txt", "-unorm")
        
    # Calculate Damping Ratio in X Dir.
    displacementX = np.array(dispX)
    damping_ratioX = S05.DAMPING_RATIO(displacementX)
    
    # Calculate Damping Ratio in Y Dir.
    displacementY = np.array(dispY)
    damping_ratioY = S05.DAMPING_RATIO(displacementY)
    
    # Calculate Damping Ratio in Z Dir.
    displacementZ = np.array(dispY)
    damping_ratioZ = S05.DAMPING_RATIO(displacementZ)
    
    # OUTPUTED DATA
    DATA = (time, reactionX, dispX, velX, accelX,
            reactionY, dispY, velY, accelY,
            reactionZ, dispZ, velZ, accelZ,
            stiffnessX, stiffnessY, stiffnessZ,
            PERIOD, 
            damping_ratioX, damping_ratioY, damping_ratioZ,
            FDx, FSx, FIx,
            FDy, FSy, FIy,
            FDz, FSz, FIz,
            )
    
    return DATA

#%% ---------------------------------------------
# Analysis Durations for Dynamic Analysis:
starttime = TI.process_time()

SPRING_TYPE = 'ELASTIC'
DATA = EXTERNAL_TIME_DEPENDENT_LOAD_MDOF(SPRING_TYPE, duration, dt)
(time, reactionXE, dispXE, velXE, accelXE,
 reactionYE, dispYE, velYE, accelYE,
 reactionZE, dispZE, velZE, accelZE,
 stiffnessXE, stiffnessYE, stiffnessZE,
 periodE,
 E_damping_ratioXE, E_damping_ratioYE, E_damping_ratioZE,
 FDXe, FSXe, FIXe,
 FDYe, FSYe, FIYe,
 FDZe, FSZe, FIZe,
 ) = DATA

S02.PERIOD_FUN(dispZE, dt)

SPRING_TYPE = 'INELASTIC'
DATA = EXTERNAL_TIME_DEPENDENT_LOAD_MDOF(SPRING_TYPE, duration, dt)
(time, reactionXI, dispXI, velXI, accelXI,
 reactionYI, dispYI, velYI, accelYI,
 reactionZI, dispZI, velZI, accelZI,
 stiffnessXI, stiffnessYI, stiffnessZI,
 periodI,
 E_damping_ratioXI, E_damping_ratioYI, E_damping_ratioZI,
 FDXi, FSXi, FIXi,
 FDYi, FSYi, FIYi,
 FDZi, FSZi, FIZi,
 ) = DATA

S02.PERIOD_FUN(dispZI, dt)

totaltime = TI.process_time() - starttime
print(f'\nTotal Analysis Durations (s): {totaltime:.4f} \n\n')
#%% ---------------------------------------------
# IN X DIR.
# Plot Results
plt.figure(2, figsize=(18, 14))

# Define colors
elastic_color = 'black'
inelastic_color = 'red'

# Reaction plot
plt.subplot(6, 1, 1)
plt.plot(time, reactionXE, color=elastic_color, linewidth=1.5)
plt.plot(time, reactionXI, color=inelastic_color, linewidth=1.5)
plt.title('Reaction Forces vs Time in X Dir. ', fontsize=12, pad=10)
plt.ylabel('Reaction (N)', fontsize=10)
plt.grid(alpha=0.3)


# Displacement plot
plt.subplot(6, 1, 2)
plt.plot(time, dispXE, color=elastic_color, linewidth=1.5, label=f'Elastic Damping Ratio: {E_damping_ratioXE:.3e} %')
plt.plot(time, dispXI, color=inelastic_color, linewidth=1.5, label=f'Inelastic Damping Ratio: {E_damping_ratioXI:.3e} %')
plt.title('Displacement vs Time in X Dir. ', fontsize=12, pad=10)
plt.ylabel('Displacement (m)', fontsize=10)
plt.legend(loc='upper right', framealpha=1)
plt.grid(alpha=0.3)

# Velocity plot
plt.subplot(6, 1, 3)
plt.plot(time, velXE, color=elastic_color, linewidth=1.5)
plt.plot(time, velXI, color=inelastic_color, linewidth=1.5)
plt.title('Velocity vs Time in X Dir. ', fontsize=12, pad=10)
plt.ylabel('Velocity (m/s)', fontsize=10)
plt.grid(alpha=0.3)

# Acceleration plot
plt.subplot(6, 1, 4)
plt.plot(time, accelXE, color=elastic_color, linewidth=1.5)
plt.plot(time, accelXI, color=inelastic_color, linewidth=1.5)
plt.title('Acceleration vs Time in X Dir. ', fontsize=12, pad=10)
plt.ylabel('Acceleration (m/s²)', fontsize=10)
plt.grid(alpha=0.3)

# Stiffness plot
plt.subplot(6, 1, 5)
plt.plot(time, stiffnessXE, color=elastic_color, linewidth=1.5)
plt.plot(time, stiffnessXI, color=inelastic_color, linewidth=1.5)
plt.title('Stiffness vs Time in X Dir. ', fontsize=12, pad=10)
plt.ylabel('Stiffness (N/m)', fontsize=10)
plt.semilogy()
plt.grid(alpha=0.3)

# Period plot
plt.subplot(6, 1, 6)
plt.plot(time, periodE, color=elastic_color, linewidth=1.5, label='Elastic Period')
plt.plot(time, periodI, color=inelastic_color, linewidth=1.5, label='Inelastic Period')
plt.title('Period vs Time in X Dir. ', fontsize=12, pad=10)
plt.ylabel('Period  (s)', fontsize=10)
plt.xlabel('Time (s)', fontsize=10)
plt.legend(loc='upper right', framealpha=1)
plt.semilogy()
plt.grid(alpha=0.3)

plt.tight_layout()
plt.show()
#%%--------------------------------------------------------
# IN Y DIR.
plt.figure(2, figsize=(18, 14))

# Define colors
elastic_color = 'black'
inelastic_color = 'red'

# Reaction plot
plt.subplot(6, 1, 1)
plt.plot(time, reactionYE, color=elastic_color, linewidth=1.5)
plt.plot(time, reactionYI, color=inelastic_color, linewidth=1.5)
plt.title('Reaction Forces vs Time in Y Dir. ', fontsize=12, pad=10)
plt.ylabel('Reaction (N)', fontsize=10)
plt.grid(alpha=0.3)


# Displacement plot
plt.subplot(6, 1, 2)
plt.plot(time, dispYE, color=elastic_color, linewidth=1.5, label=f'Elastic Damping Ratio: {E_damping_ratioYE:.3e} %')
plt.plot(time, dispYI, color=inelastic_color, linewidth=1.5, label=f'Inelastic Damping Ratio: {E_damping_ratioYI:.3e} %')
plt.title('Displacement vs Time in Y Dir. ', fontsize=12, pad=10)
plt.ylabel('Displacement (m)', fontsize=10)
plt.legend(loc='upper right', framealpha=1)
plt.grid(alpha=0.3)

# Velocity plot
plt.subplot(6, 1, 3)
plt.plot(time, velYE, color=elastic_color, linewidth=1.5)
plt.plot(time, velYI, color=inelastic_color, linewidth=1.5)
plt.title('Velocity vs Time in Y Dir. ', fontsize=12, pad=10)
plt.ylabel('Velocity (m/s)', fontsize=10)
plt.grid(alpha=0.3)

# Acceleration plot
plt.subplot(6, 1, 4)
plt.plot(time, accelYE, color=elastic_color, linewidth=1.5)
plt.plot(time, accelYI, color=inelastic_color, linewidth=1.5)
plt.title('Acceleration vs Time in Y Dir. ', fontsize=12, pad=10)
plt.ylabel('Acceleration (m/s²)', fontsize=10)
plt.grid(alpha=0.3)

# Stiffness plot
plt.subplot(6, 1, 5)
plt.plot(time, stiffnessYE, color=elastic_color, linewidth=1.5)
plt.plot(time, stiffnessYI, color=inelastic_color, linewidth=1.5)
plt.title('Stiffness vs Time in Y Dir. ', fontsize=12, pad=10)
plt.ylabel('Stiffness (N/m)', fontsize=10)
plt.semilogy()
plt.grid(alpha=0.3)

# Period plot
plt.subplot(6, 1, 6)
plt.plot(time, periodE, color=elastic_color, linewidth=1.5, label='Elastic Period')
plt.plot(time, periodI, color=inelastic_color, linewidth=1.5, label='Inelastic Period')
plt.title('Period vs Time in Y Dir. ', fontsize=12, pad=10)
plt.ylabel('Period  (s)', fontsize=10)
plt.xlabel('Time (s)', fontsize=10)
plt.legend(loc='upper right', framealpha=1)
plt.semilogy()
plt.grid(alpha=0.3)

plt.tight_layout()
plt.show()
#%%--------------------------------------------------------
# IN Z DIR.
plt.figure(2, figsize=(18, 14))

# Define colors
elastic_color = 'black'
inelastic_color = 'red'

# Reaction plot
plt.subplot(6, 1, 1)
plt.plot(time, reactionZE, color=elastic_color, linewidth=1.5)
plt.plot(time, reactionZI, color=inelastic_color, linewidth=1.5)
plt.title('Reaction Forces vs Time in Z Dir. ', fontsize=12, pad=10)
plt.ylabel('Reaction (N.m)', fontsize=10)
plt.grid(alpha=0.3)


# Displacement plot
plt.subplot(6, 1, 2)
plt.plot(time, dispZE, color=elastic_color, linewidth=1.5, label=f'Elastic Damping Ratio: {E_damping_ratioZE:.3e} %')
plt.plot(time, dispZI, color=inelastic_color, linewidth=1.5, label=f'Inelastic Damping Ratio: {E_damping_ratioZI:.3e} %')
plt.title('Displacement vs Time in Y Dir. ', fontsize=12, pad=10)
plt.ylabel('Displacement (rad)', fontsize=10)
plt.legend(loc='upper right', framealpha=1)
plt.grid(alpha=0.3)

# Velocity plot
plt.subplot(6, 1, 3)
plt.plot(time, velZE, color=elastic_color, linewidth=1.5)
plt.plot(time, velZI, color=inelastic_color, linewidth=1.5)
plt.title('Velocity vs Time in Y Dir. ', fontsize=12, pad=10)
plt.ylabel('Velocity (rad/s)', fontsize=10)
plt.grid(alpha=0.3)

# Acceleration plot
plt.subplot(6, 1, 4)
plt.plot(time, accelZE, color=elastic_color, linewidth=1.5)
plt.plot(time, accelZI, color=inelastic_color, linewidth=1.5)
plt.title('Acceleration vs Time in Y Dir. ', fontsize=12, pad=10)
plt.ylabel('Acceleration (rad/s²)', fontsize=10)
plt.grid(alpha=0.3)

# Stiffness plot
plt.subplot(6, 1, 5)
plt.plot(time, stiffnessZE, color=elastic_color, linewidth=1.5)
plt.plot(time, stiffnessZI, color=inelastic_color, linewidth=1.5)
plt.title('Stiffness vs Time in Y Dir. ', fontsize=12, pad=10)
plt.ylabel('Stiffness (N/m/rad)', fontsize=10)
plt.semilogy()
plt.grid(alpha=0.3)

# Period plot
plt.subplot(6, 1, 6)
plt.plot(time, periodE, color=elastic_color, linewidth=1.5, label='Elastic Period')
plt.plot(time, periodI, color=inelastic_color, linewidth=1.5, label='Inelastic Period')
plt.title('Period vs Time in Z Dir. ', fontsize=12, pad=10)
plt.ylabel('Period  (s)', fontsize=10)
plt.xlabel('Time (s)', fontsize=10)
plt.legend(loc='upper right', framealpha=1)
plt.semilogy()
plt.grid(alpha=0.3)

plt.tight_layout()
plt.show()
#%%--------------------------------------------------------
plt.figure(3, figsize=(8, 6))
plt.plot(dispXE, reactionXE, color='black', linewidth=2)
plt.plot(dispXI, reactionXI, color='purple', linewidth=2)
plt.plot(dispXE, reactionXE, color='green', linewidth=2)
plt.plot(dispXI, reactionXI, color='red', linewidth=2)
plt.xlabel('Displacement [m]', fontsize=10)
plt.ylabel('Base-reaction [N]', fontsize=10)
plt.title('Displacement vs Base-reaction', fontsize=10)
plt.legend(['ELASTIC IN X DIR.', 'INELASTIC IN X DIR.', 'ELASTIC IN Y DIR.', 'INELASTIC IN Y DIR.'])
plt.grid()
plt.show()

plt.figure(33, figsize=(8, 6))
plt.plot(dispZE, reactionZE, color='black', linewidth=2)
plt.plot(dispZI, reactionZI, color='purple', linewidth=2)
plt.xlabel('Rotation [rad]', fontsize=10)
plt.ylabel('Base-reaction [N.m]', fontsize=10)
plt.title('Rotation vs Base-reaction', fontsize=10)
plt.legend(['ELASTIC IN Z DIR.', 'INELASTIC IN Z DIR.'])
plt.grid()
plt.show()

plt.figure(4, figsize=(8, 6))
plt.plot(velXE, FDXe, color='black', linewidth=2)
plt.plot(velXI, FDXi, color='purple', linewidth=2)
plt.plot(velYE, FDYe, color='green', linewidth=2)
plt.plot(velYI, FDYi, color='red', linewidth=2)
plt.xlabel('Velocity (m/s)', fontsize=10)
plt.ylabel('Damping Force (fD) [N]', fontsize=10)
plt.title('Damping Force (fD) vs Velocity Curve', fontsize=10)
plt.legend(['ELASTIC IN X DIR.', 'INELASTIC IN X DIR.', 'ELASTIC IN Y DIR.', 'INELASTIC IN Y DIR.'])
plt.grid()
plt.show()

plt.figure(44, figsize=(8, 6))
plt.plot(velZE, FDZe, color='black', linewidth=2)
plt.plot(velZI, FDZi, color='purple', linewidth=2)
plt.xlabel('Velocity (rad/s)', fontsize=10)
plt.ylabel('Damping Force (fD) [N.m]', fontsize=10)
plt.title('Damping Force (fD) vs Velocity Curve', fontsize=10)
plt.legend(['ELASTIC IN Z DIR.', 'INELASTIC IN Z DIR.'])
plt.grid()
plt.show()

plt.figure(5, figsize=(8, 6))
plt.plot(accelXE, FIXe, color='black', linewidth=2)
plt.plot(accelXI, FIXi, color='purple', linewidth=2)
plt.plot(accelYE, FIYe, color='green', linewidth=2)
plt.plot(accelYI, FIYi, color='red', linewidth=2)
plt.xlabel('Acceleration [m/s²]', fontsize=10)
plt.ylabel('Inertia Force (fI) [N]', fontsize=10)
plt.title('Inertia Force (fI) vs Acceleration Curve', fontsize=10)
plt.legend(['ELASTIC IN X DIR.', 'INELASTIC IN X DIR.', 'ELASTIC IN Y DIR.', 'INELASTIC IN Y DIR.'])
plt.grid()
plt.show()

plt.figure(55, figsize=(8, 6))
plt.plot(accelZE, FIZe, color='black', linewidth=2)
plt.plot(accelZI, FIZi, color='purple', linewidth=2)
plt.xlabel('Acceleration [rad/s²]', fontsize=10)
plt.ylabel('Inertia Force (fI) [N.m]', fontsize=10)
plt.title('Inertia Force (fI) vs Acceleration Curve', fontsize=10)
plt.legend(['ELASTIC IN Z DIR.', 'INELASTIC IN Z DIR.'])
plt.grid()
plt.show()

plt.figure(6, figsize=(8, 6))
plt.plot(dispXE, FSXe, color='black', linewidth=2)
plt.plot(dispXI, FSXi, color='purple', linewidth=2)
plt.plot(dispYE, FSYe, color='green', linewidth=2)
plt.plot(dispYI, FSYi, color='red', linewidth=2)
plt.xlabel('Displacement [m]', fontsize=10)
plt.ylabel('Spring Force (fS) [N]', fontsize=10)
plt.title('Spring Force (fS) vs Displacement', fontsize=10)
plt.legend(['ELASTIC IN X DIR.', 'INELASTIC IN X DIR.', 'ELASTIC IN Y DIR.', 'INELASTIC IN Y DIR.'])
plt.grid()
plt.show()

plt.figure(66, figsize=(8, 6))
plt.plot(dispZE, FSZe, color='black', linewidth=2)
plt.plot(dispZI, FSZi, color='purple', linewidth=2)
plt.xlabel('Rotation [rad]', fontsize=10)
plt.ylabel('Spring Force (fS) [N.m]', fontsize=10)
plt.title('Spring Force (fS) vs Displacement', fontsize=10)
plt.legend(['ELASTIC IN Z DIR.', 'INELASTIC IN Z DIR.'])
plt.grid()
plt.show()

plt.figure(7, figsize=(8, 6))
plt.plot(dispXE, FDXe, color='black', linewidth=2)
plt.plot(dispXI, FDXi, color='purple', linewidth=2)
plt.plot(dispYE, FDYe, color='green', linewidth=2)
plt.plot(dispYI, FDYi, color='red', linewidth=2)
plt.xlabel('Displacement [m]', fontsize=10)
plt.ylabel('Damping Force (fD) [N]', fontsize=10)
plt.title('Damping Force (fD) vs Displacement Curve')
plt.legend(['ELASTIC IN X DIR.', 'INELASTIC IN X DIR.', 'ELASTIC IN Y DIR.', 'INELASTIC IN Y DIR.'])
plt.grid()
plt.show()

plt.figure(77, figsize=(8, 6))
plt.plot(dispZE, FDZe, color='black', linewidth=2)
plt.plot(dispZI, FDZi, color='purple', linewidth=2)
plt.xlabel('Rotation [rad]', fontsize=10)
plt.ylabel('Damping Force (fD) [N.m]', fontsize=10)
plt.title('Damping Force (fD) vs Displacement Curve')
plt.legend(['ELASTIC IN Z DIR.', 'INELASTIC IN Z DIR.'])
plt.grid()
plt.show()
#%% ---------------------------------------------
# Print out the state of nodes 1 and 2
ops.printModel("node",1, 2)
# Print out the state of element 1
ops.printModel("ele", 1)
# Print the Model
#printModel()
ops.printModel("-JSON", "-file", "MDF_PENDULUM_TWO.json")
#%%-------------------------------------------------------------------------------
# EXCEL OUTPUT
import pandas as pd

# Create DataFrame function
def create_df(dispXE, dispXI, velXE, velXI, accelXE, accelXI, reactionXE, reactionXI, FDXe, FDXi, FIXe, FIXi, FSXe, FSXi):
    df = pd.DataFrame({
        "DISPLACEMENT [m] - ELASTIC": dispXE,
        "VELOCITY [m/s] - ELASTIC": velXE,
        "ACCELERATION [m/s^2] - ELASTIC": accelXE,
        "BASE-REACTION [N] - ELASTIC": reactionXE,
        "DAMPING-FORCE [N] - ELASTIC": FDXe,
        "INERTIA-FORCE [N] - ELASTIC": FIXe,
        "SPRING-FORCE [N] - ELASTIC": FSXe,
        "DISPLACEMENT [m] - INELASTIC": dispXI,
        "VELOCITY [m/s] - INELASTIC": velXI,
        "ACCELERATION [m/s^2] - INELASTIC": accelXI,
        "BASE-REACTION [N] - INELASTIC": reactionXI,
        "DAMPING-FORCE [N] - INELASTIC": FDXi,
        "INERTIA-FORCE [N] - INELASTIC": FIXi,
        "SPRING-FORCE [N] - INELASTIC": FSXi,   
    })
    return df


# Save to Excel
with pd.ExcelWriter("MDF_PENDULUM_TWO_OUTPUT.xlsx", engine='openpyxl') as writer:
    
    df1 = create_df(dispXE, dispXI, velXE, velXI, accelXE, accelXI, reactionXE, reactionXI, FDXe, FDXi, FIXe, FIXi, FSXe, FSXi)
    df1.to_excel(writer, sheet_name="OUTPUT", index=False)    
#%%------------------------------------------------------------------------------------------------