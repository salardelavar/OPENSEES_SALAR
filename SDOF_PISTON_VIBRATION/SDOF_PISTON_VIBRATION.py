######################################################################################################################
#                          >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                          #
#                              MODELING OF PISTON VIBRATION SDOF STRUCTURE USING OPENSEES                            #
#                                                P(t) = P0 cos(wt)                                                   #
#                                         P(t) = P0 exp(-0.05wt) cos(wt)                                             #
#--------------------------------------------------------------------------------------------------------------------#
#                   EVALUATION OF DAMPING FORCE (fD), SPRING FORCE (fS) AND INERTIA FORCE (fI)                       #
#--------------------------------------------------------------------------------------------------------------------#
#                              THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                            #
#                                       EMAIL: salar.d.ghashghaei@gmail.com                                          #
######################################################################################################################
"""
1. This code models a single-degree-of-freedom (SDOF) piston–spring–damper system using nonlinear dynamic analysis.
2. The system is excited by harmonic and exponentially decaying harmonic external forces ( P(t) ).
3. Both elastic and inelastic (hysteretic) material behaviors are considered for the spring.
4. Inelastic behavior captures yielding, post-yield hardening, stiffness degradation, and energy dissipation.
5. The model is implemented in OpenSees using a zeroLength element to isolate material response.
6. Time integration is performed with the Newmark-β method, ensuring numerical stability for nonlinear dynamics.
7. Displacement, velocity, and acceleration responses are recorded at every time step.
8. The total resisting force is decomposed into inertia, damping, and spring forces to verify dynamic equilibrium.
9. Instantaneous stiffness and natural period are continuously updated, showing period elongation in the inelastic range.
10. The results allow advanced evaluation of nonlinear vibration behavior, damping mechanisms, and energy dissipation.
"""
#%%-----------------------------------------------------------------------------
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
FY = 850.0                        # [N] Yield Force of Structure
FU = 1.1 * FY                     # [N] Ultimate Force of Structure
Ke = 45000.0                      # [N/m] Spring Elastic Stiffness
DY = FY / Ke                      # [m] Yield Displacement
DSU = 0.36                        # [m] Ultimate Displacement
Ksh = (FU - FY) / (DSU - DY)      # [N/m] Displacement Hardening Modulus
Kp = FU / DSU                     # [N/m] Spring Plastic Stiffness
b = Ksh / Ke                      # Displacement Hardening Ratio

M = 1.0                           # [kg] Mass
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
pos_force = [0, FY, FU, 0.2*FU, 0.1*FU]
KP = np.array([FY, DY, FU, DSU, 0.2*FU, 1.1*DSU, 0.1*FU, 1.25*DSU])

# Negative branch points
neg_disp = [0, -DY, -DSU, -1.1*DSU, -1.25*DSU]
neg_force = [0, -FY, -FU, -0.2*FU, -0.1*FU]
KN = np.array([-FY, -DY, -FU, -DSU, -0.2*FU, -1.1*DSU, -0.1*FU, -1.25*DSU])
#%%------------------------------------------------------------------------------------------------
# Plot
plt.plot(pos_disp, pos_force, marker='o', color='red')
plt.plot(neg_disp, neg_force, marker='o', color='black')

plt.xlabel("Displacement [m]")
plt.ylabel("Force [N]")
plt.title("Force–Displacement Curve")
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
force_amplitude = 99700.0 # [N] Harmonic Load
omega_DT = 5.0715         # [rad/s] Natural angular frequency

# Check function
def CHECK_FUN(DT_time ,duration):
    if DT_time > duration:
        print('\n\nAnalysis Duration Must be greater than External Load Duration\n\n')        
        exit()
    return -1

CHECK_FUN(DT_time ,duration)

def EXTERNAL_TIME_DEPENDENT(force_amplitude, omega_DT, DT, DT_time):
    # P(t) = P0 cos(wt)
    import numpy as np
    import matplotlib.pyplot as plt
    # External Load Durations
    num_steps = int(DT_time / DT)
    load_time = np.linspace(0, DT_time, num_steps) 
    target_frequency = 1.0 * omega_DT  # Target excitation frequency
    DT_load = force_amplitude * np.cos(target_frequency * load_time)
    # Plot External Time-dependent Loading
    plt.figure(figsize=(10, 6))
    plt.plot(load_time, DT_load, label=f'External Loading - Max: {np.max(DT_load):.3f}', linewidth=5)
    plt.title('External Time-dependent Loading Over Time')
    plt.xlabel('Time (s)')
    plt.ylabel('Force (N)')
    plt.grid(True)
    plt.legend()
    plt.show()
    return DT_load

#DT_load = EXTERNAL_TIME_DEPENDENT(force_amplitude, omega_DT, DT, DT_time)

def EXTERNAL_TIME_DEPENDENT_02(force_amplitude, omega_DT, DT, DT_time):
    # P(t) = P0 exp(-0.05wt) cos(wt)
    import numpy as np
    import matplotlib.pyplot as plt
    # External Load Durations
    num_steps = int(DT_time / DT)
    load_time = np.linspace(0, DT_time, num_steps) 
    target_frequency = 1.0 * omega_DT  # Target excitation frequency
    DT_load = force_amplitude * np.exp(-0.05*target_frequency * load_time) * np.cos(target_frequency * load_time)
    # Plot External Time-dependent Loading
    plt.figure(figsize=(10, 6))
    plt.plot(load_time, DT_load, label=f'External Loading - Max: {np.max(DT_load):.3f}', linewidth=5)
    plt.title('External Time-dependent Loading Over Time')
    plt.xlabel('Time (s)')
    plt.ylabel('Force (N)')
    plt.grid(True)
    plt.legend()
    plt.show()
    print(load_time, DT_load)
    return DT_load

#DT_load02 = EXTERNAL_TIME_DEPENDENT_02(force_amplitude, omega_DT, DT, DT_time)

#%% ---------------------------------------------
def EXTERNAL_TIME_DEPENDENT_LOAD_SDOF(SPRING_TYPE, duration, dt):
    # Create model
    ops.wipe()
    # Create a 2D model with 3 DOF per node
    ops.model('Basic', '-ndm', 1, '-ndf', 1)
    
    # Add nodes
    ops.node(1, 0.0)
    ops.node(2, 0.0)
    
    
    # Fix node 1
    ops.fix(1, 1)
        
    # Define material and element
    if SPRING_TYPE == 'ELASTIC': 
        ops.uniaxialMaterial('Elastic', 1, Ke, C)
        ops.element('zeroLength', 10, 1, 2,'-mat', 1,'-dir', 1) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/ZeroLength.html
    if SPRING_TYPE == 'INELASTIC':
        ops.uniaxialMaterial('HystereticSM', 1, '-posEnv', *KP.flatten(), '-negEnv', *KN.flatten(), '-pinch', 1, 1)
        #INFO LINK: https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/HystereticSM.html
        ops.uniaxialMaterial('Viscous', 2, C, 1.0)  # Material for C (alpha=1.0 for linear)
        ops.element('zeroLength', 10, 1, 2,'-mat', 1, 2,'-dir', 1, 1)
    
    # Define mass to node 2
    ops.mass(2, M)
    
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
    
    # P(t) = P0 cos(wt)  
    #DT_load = EXTERNAL_TIME_DEPENDENT(M, omega_DT, DT, DT_time)
    # P(t) = P0 exp(-0.05wt) cos(wt)
    DT_load = EXTERNAL_TIME_DEPENDENT_02(M, omega_DT, DT, DT_time)
    
    # Static analysis
    time_series_tag = 1
    pattern_tag = 1
    # Apply time-dependent explosion loading
    ops.timeSeries('Path', time_series_tag, '-dt', dt, '-values', *DT_load)
    ops.pattern('Plain', pattern_tag, time_series_tag)
    ops.load(2, 1.0)

    # Perform analysis
    time = []
    disp = []
    vel = []
    accel = []
    reaction = []
    stiffness = []
    OMEGA, PERIOD = [], []
    PERIOD_MIN, PERIOD_MAX = [], []
    FD, FS, FI = [], [], []
    
    stable = 0
    current_time = 0.0
        
    while stable == 0 and current_time < duration:
        ops.analyze(1, dt)
        S01.ANALYSIS(stable, 1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS
        current_time = ops.getTime()
        time.append(current_time)
        ops.reactions()
        reaction.append(ops.nodeReaction(1, 1)) # BASE REACTION
        disp.append(ops.nodeDisp(2, 1))         # DISPLACEMENT  
        vel.append(ops.nodeVel(2, 1))           # VELOCITY
        accel.append(ops.nodeAccel(2, 1))       # ACCELERATION
        stiffness.append(np.abs(reaction[-1]) / np.abs(disp[-1]))
        OMEGA.append(np.sqrt(stiffness[-1]/M))
        PERIOD.append((np.pi * 2) / OMEGA[-1])
        # IN EACH STEP STRUCTURAL PERIOD WILL BE CALCULATED
        PERIODmin, PERIODmax = S03.EIGENVALUE_ANALYSIS(1, PLOT=True)
        PERIOD_MIN.append(PERIODmin)
        PERIOD_MAX.append(PERIODmax)
        Ci = 2 * zi * (2*np.pi/PERIODmin) * M    # Damping COEFFICIENT - UPDATED AND CHANGES IN EACH STEP
        FD.append(Ci * vel[-1])                  # DAMPING FORCE
        FS.append(stiffness[-1] * disp[-1])      # SPRING FORCE
        FI.append(M * accel[-1])                 # INERTIA FORCE
        print(time[-1], disp[-1], vel[-1])

    # Compute modal properties
    ops.modalProperties("-print", "-file", "SALAR_ModalReport.txt", "-unorm")
        
    # Calculate Damping Ratio
    displacement = np.array(disp)
    damping_ratio = S05.DAMPING_RATIO(displacement)
    
    # OUTPUTED DATA
    DATA = (time, reaction, disp, vel, accel,
            stiffness, PERIOD, damping_ratio,
            FD, FS, FI)
    
    return DATA

#%% ---------------------------------------------
# Analysis Durations for Dynamic Analysis:
starttime = TI.process_time()

SPRING_TYPE = 'ELASTIC'
DATA = EXTERNAL_TIME_DEPENDENT_LOAD_SDOF(SPRING_TYPE, duration, dt)
(time, reactionE, dispE, velE, accelE,
 stiffnessE, periodE, E_damping_ratioE,
 FDe, FSe, FIe) = DATA

S02.PERIOD_FUN(dispE, dt)

SPRING_TYPE = 'INELASTIC'
DATA = EXTERNAL_TIME_DEPENDENT_LOAD_SDOF(SPRING_TYPE, duration, dt)
(time, reactionI, dispI, velI, accelI,
 stiffnessI, periodI, E_damping_ratioI,
 FDi, FSi, FIi) = DATA

S02.PERIOD_FUN(dispI, dt)

totaltime = TI.process_time() - starttime
print(f'\nTotal Analysis Durations (s): {totaltime:.4f} \n\n')
#%% ---------------------------------------------
# Plot Results
plt.figure(2, figsize=(18, 14))

# Define colors
elastic_color = 'black'
inelastic_color = 'red'

# Reaction plot
plt.subplot(6, 1, 1)
plt.plot(time, reactionE, color=elastic_color, linewidth=1.5)
plt.plot(time, reactionI, color=inelastic_color, linewidth=1.5)
plt.title('Reaction Forces vs Time', fontsize=12, pad=10)
plt.ylabel('Reaction (N)', fontsize=10)
plt.grid(alpha=0.3)


# Displacement plot
plt.subplot(6, 1, 2)
plt.plot(time, dispE, color=elastic_color, linewidth=1.5, label=f'Elastic Damping Ratio: {E_damping_ratioE:.3e} %')
plt.plot(time, dispI, color=inelastic_color, linewidth=1.5, label=f'Inelastic Damping Ratio: {E_damping_ratioI:.3e} %')
plt.title('Displacement vs Time', fontsize=12, pad=10)
plt.ylabel('Displacement (m)', fontsize=10)
plt.legend(loc='upper right', framealpha=1)
plt.grid(alpha=0.3)

# Velocity plot
plt.subplot(6, 1, 3)
plt.plot(time, velE, color=elastic_color, linewidth=1.5)
plt.plot(time, velI, color=inelastic_color, linewidth=1.5)
plt.title('Velocity vs Time', fontsize=12, pad=10)
plt.ylabel('Velocity (m/s)', fontsize=10)
plt.grid(alpha=0.3)

# Acceleration plot
plt.subplot(6, 1, 4)
plt.plot(time, accelE, color=elastic_color, linewidth=1.5)
plt.plot(time, accelI, color=inelastic_color, linewidth=1.5)
plt.title('Acceleration vs Time', fontsize=12, pad=10)
plt.ylabel('Acceleration (m/s²)', fontsize=10)
plt.grid(alpha=0.3)

# Stiffness plot
plt.subplot(6, 1, 5)
plt.plot(time, stiffnessE, color=elastic_color, linewidth=1.5)
plt.plot(time, stiffnessI, color=inelastic_color, linewidth=1.5)
plt.title('Stiffness vs Time', fontsize=12, pad=10)
plt.ylabel('Stiffness (N/m)', fontsize=10)
plt.semilogy()
plt.grid(alpha=0.3)

# Period plot
plt.subplot(6, 1, 6)
plt.plot(time, periodE, color=elastic_color, linewidth=1.5, label='Elastic Period')
plt.plot(time, periodI, color=inelastic_color, linewidth=1.5, label='Inelastic Period')
plt.title('Period vs Time', fontsize=12, pad=10)
plt.ylabel('Period  (s)', fontsize=10)
plt.xlabel('Time (s)', fontsize=10)
plt.legend(loc='upper right', framealpha=1)
plt.semilogy()
plt.grid(alpha=0.3)

plt.tight_layout()
plt.show()

plt.figure(3, figsize=(8, 6))
plt.plot(dispE, reactionE, color='black', linewidth=2)
plt.plot(dispI, reactionI, color='purple', linewidth=2)
plt.xlabel('Displacement [m]', fontsize=10)
plt.ylabel('Base-reaction [N]', fontsize=10)
plt.title('Displacement vs Base-reaction', fontsize=10)
plt.legend(['ELASTIC', 'INELASTIC'])
plt.grid()
plt.show()

plt.figure(4, figsize=(8, 6))
plt.plot(velE, FDe, color='black', linewidth=2)
plt.plot(velI, FDi, color='purple', linewidth=2)
plt.xlabel('Velocity (m/s)', fontsize=10)
plt.ylabel('Damping Force (fD) [N]', fontsize=10)
plt.title('Damping Force (fD) vs Velocity Curve', fontsize=10)
plt.legend(['ELASTIC', 'INELASTIC'])
plt.grid()
plt.show()

plt.figure(5, figsize=(8, 6))
plt.plot(accelE, FIe, color='black', linewidth=2)
plt.plot(accelI, FIi, color='purple', linewidth=2)
plt.xlabel('Acceleration [m/s²]', fontsize=10)
plt.ylabel('Inertia Force (fI) [N]', fontsize=10)
plt.title('Inertia Force (fI) vs Acceleration Curve', fontsize=10)
plt.legend(['ELASTIC', 'INELASTIC'])
plt.grid()
plt.show()

plt.figure(6, figsize=(8, 6))
plt.plot(dispE, FSe, color='black', linewidth=2)
plt.plot(dispI, FSi, color='purple', linewidth=2)
plt.xlabel('Displacement [m]', fontsize=10)
plt.ylabel('Spring Force (fS) [N]', fontsize=10)
plt.title('Spring Force (fS) vs Displacement', fontsize=10)
plt.legend(['ELASTIC', 'INELASTIC'])
plt.grid()
plt.show()

plt.figure(7, figsize=(8, 6))
plt.plot(dispE, FDe, color='black', linewidth=2)
plt.plot(dispI, FDi, color='purple', linewidth=2)
plt.xlabel('Displacement [m]', fontsize=10)
plt.ylabel('Damping Force (fD) [N]', fontsize=10)
plt.title('Damping Force (fD) vs Displacement Curve')
plt.legend(['ELASTIC', 'INELASTIC'])
plt.grid()
plt.show()
#%% ---------------------------------------------
# Print out the state of nodes 1 and 2
ops.printModel("node",1, 2)
# Print out the state of element 1
ops.printModel("ele", 1)
# Print the Model
#printModel()
ops.printModel("-JSON", "-file", "SDOF_PISTON_VIBRATION.json")
#%%-------------------------------------------------------------------------------
# EXCEL OUTPUT
import pandas as pd

# Create DataFrame function
def create_df(dispE, dispI, velE, velI, accelE, accelI, reactionE, reactionI, FDe, FDi, FIe, FIi, FSe, FSi):
    df = pd.DataFrame({
        "DISPLACEMENT [m] - ELASTIC": dispE,
        "VELOCITY [m/s] - ELASTIC": velE,
        "ACCELERATION [m/s^2] - ELASTIC": accelE,
        "BASE-REACTION [N] - ELASTIC": reactionE,
        "DAMPING-FORCE [N] - ELASTIC": FDe,
        "INERTIA-FORCE [N] - ELASTIC": FIe,
        "SPRING-FORCE [N] - ELASTIC": FSe,
        "DISPLACEMENT [m] - INELASTIC": dispI,
        "VELOCITY [m/s] - INELASTIC": velI,
        "ACCELERATION [m/s^2] - INELASTIC": accelI,
        "BASE-REACTION [N] - INELASTIC": reactionI,
        "DAMPING-FORCE [N] - INELASTIC": FDi,
        "INERTIA-FORCE [N] - INELASTIC": FIi,
        "SPRING-FORCE [N] - INELASTIC": FSi,   
    })
    return df


# Save to Excel
with pd.ExcelWriter("SDOF_PISTON_VIBRATION_OUTPUT.xlsx", engine='openpyxl') as writer:
    
    df1 = create_df(dispE, dispI, velE, velI, accelE, accelI, reactionE, reactionI, FDe, FDi, FIe, FIi, FSe, FSi)
    df1.to_excel(writer, sheet_name="OUTPUT", index=False)
#%%------------------------------------------------------------------------------------------------