######################################################################################################################
#                          >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                          #
#        COMPARATIVE FREE-VIBRATION ANALYSIS OF A MDOF STRUCTURE: ELASTIC VS INELASTIC RESPONSE USING OPENSEES       #
#--------------------------------------------------------------------------------------------------------------------#
#                              THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                            #
#                                       EMAIL: salar.d.ghashghaei@gmail.com                                          #
######################################################################################################################
"""
Performs free-vibration analysis of a Multi Degree of Freedom (MDOF)
 structure using OpenSeesPy, comparing elastic and inelastic spring behavior. 
---------------------------------------------------------------------------- 
Key features include:
1. Implements both elastic (linear) and hysteretic (nonlinear) material models for structural springs.
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
allowing evaluation of structural response under free vibration
 with different material nonlinearities and damping characteristics.
 The hysteretic material model captures energy dissipation 
 inelastic deformation, while the elastic case serves as a reference for linear behavior.
"""
import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
import time as TI
import ANALYSIS_FUNCTION as S01
import PERIOD_FUN as S02
import DAMPING_RATIO_FUN as S03
import EIGENVALUE_ANALYSIS_FUN as S04
#%%-------------------------------------------------------------------------------
# Define parameters
M = [150.0, 100.0, 250.0, 350.0]   # [kg] Mass
zi = [0.05, 0.01, 0.02, 0.03]      # Damping ratio
u0 = -0.035                        # [m] Initial displacement
v0 = 0.015                         # [m/s] Initial velocity
a0 = 0.0065                        # [m/s^2] Initial acceleration
duration = 10.0                    # [s] Analysis duration
dt = 0.01                          # [s] Time step

#%%-------------------------------------------------------------------------------
#%% Inelastic Force-Displacement Parameters for 4 Different Inelastic Springs

# Element 1 - Primary Spring (Original parameters)
FY1_1, DY1_1 = 2772.0, 0.01     # First yield point
FY2_1, DY2_1 = 3104.6, 0.05     # Peak force
FY3_1, DY3_1 = 1663.2, 0.08     # Post-peak softening
FY4_1, DY4_1 = 1663.2, 0.17     # Plateau
FY5_1, DY5_1 = 277.2, 0.28      # Further softening
FY6_1, DY6_1 = 200.0, 0.41      # Near-zero stiffness
FY7_1, DY7_1 = 0.0, 0.52        # Zero force (failure)

KP_1 = np.array([FY1_1, DY1_1, FY2_1, DY2_1, FY3_1, DY3_1, FY4_1, DY4_1, 
                FY5_1, DY5_1, FY6_1, DY6_1, FY7_1, DY7_1])

# Compression for Element 1
FY1n_1, DY1n_1 = -2772.0, -0.01    # First yield in compression
FY2n_1, DY2n_1 = -3104.6, -0.20    # Peak compressive force
FY3n_1, DY3n_1 = -1663.2, -0.76    # Post-peak softening
KN_1 = np.array([FY1n_1, DY1n_1, FY2n_1, DY2n_1, FY3n_1, DY3n_1])

# Element 2 - Stiffer Spring (50% stiffer)
stiffness_factor_2 = 1.5
FY1_2, DY1_2 = 2772.0 * stiffness_factor_2, 0.01
FY2_2, DY2_2 = 3104.6 * stiffness_factor_2, 0.05
FY3_2, DY3_2 = 1663.2 * stiffness_factor_2, 0.08
FY4_2, DY4_2 = 1663.2 * stiffness_factor_2, 0.17
FY5_2, DY5_2 = 277.2 * stiffness_factor_2, 0.28
FY6_2, DY6_2 = 200.0 * stiffness_factor_2, 0.41
FY7_2, DY7_2 = 0.0, 0.52

KP_2 = np.array([FY1_2, DY1_2, FY2_2, DY2_2, FY3_2, DY3_2, FY4_2, DY4_2, 
                FY5_2, DY5_2, FY6_2, DY6_2, FY7_2, DY7_2])

# Compression for Element 2
FY1n_2, DY1n_2 = -2772.0 * stiffness_factor_2, -0.01
FY2n_2, DY2n_2 = -3104.6 * stiffness_factor_2, -0.20
FY3n_2, DY3n_2 = -1663.2 * stiffness_factor_2, -0.76
KN_2 = np.array([FY1n_2, DY1n_2, FY2n_2, DY2n_2, FY3n_2, DY3n_2])

# Element 3 - Softer Spring (30% softer)
stiffness_factor_3 = 0.7
FY1_3, DY1_3 = 2772.0 * stiffness_factor_3, 0.01
FY2_3, DY2_3 = 3104.6 * stiffness_factor_3, 0.05
FY3_3, DY3_3 = 1663.2 * stiffness_factor_3, 0.08
FY4_3, DY4_3 = 1663.2 * stiffness_factor_3, 0.17
FY5_3, DY5_3 = 277.2 * stiffness_factor_3, 0.28
FY6_3, DY6_3 = 200.0 * stiffness_factor_3, 0.41
FY7_3, DY7_3 = 0.0, 0.52

KP_3 = np.array([FY1_3, DY1_3, FY2_3, DY2_3, FY3_3, DY3_3, FY4_3, DY4_3, 
                FY5_3, DY5_3, FY6_3, DY6_3, FY7_3, DY7_3])

# Compression for Element 3
FY1n_3, DY1n_3 = -2772.0 * stiffness_factor_3, -0.01
FY2n_3, DY2n_3 = -3104.6 * stiffness_factor_3, -0.20
FY3n_3, DY3n_3 = -1663.2 * stiffness_factor_3, -0.76
KN_3 = np.array([FY1n_3, DY1n_3, FY2n_3, DY2n_3, FY3n_3, DY3n_3])

# Element 4 - Different behavior (modified yield points)
FY1_4, DY1_4 = 2000.0, 0.008
FY2_4, DY2_4 = 3500.0, 0.04
FY3_4, DY3_4 = 1800.0, 0.10
FY4_4, DY4_4 = 1800.0, 0.20
FY5_4, DY5_4 = 300.0, 0.30
FY6_4, DY6_4 = 150.0, 0.45
FY7_4, DY7_4 = 0.0, 0.55

KP_4 = np.array([FY1_4, DY1_4, FY2_4, DY2_4, FY3_4, DY3_4, FY4_4, DY4_4, 
                FY5_4, DY5_4, FY6_4, DY6_4, FY7_4, DY7_4])

# Compression for Element 4
FY1n_4, DY1n_4 = -2000.0, -0.008
FY2n_4, DY2n_4 = -3500.0, -0.18
FY3n_4, DY3n_4 = -1800.0, -0.70
KN_4 = np.array([FY1n_4, DY1n_4, FY2n_4, DY2n_4, FY3n_4, DY3n_4])

# Store all elements in lists for easy access
KP_list = [KP_1, KP_2, KP_3, KP_4]
KN_list = [KN_1, KN_2, KN_3, KN_4]

#%%% Plotting Force-Displacement Diagrams for All 4 Elements
plt.figure(1, figsize=(12, 8))
colors = ['red', 'blue', 'green', 'orange']
labels = ['Element 1 (Original)', 'Element 2 (Stiffer)', 'Element 3 (Softer)', 'Element 4 (Modified)']

for i, (KP, KN) in enumerate(zip(KP_list, KN_list)):
    # Separate into Force and Displacement
    force_p = KP[0::2]
    disp_p = KP[1::2]
    force_n = KN[0::2]
    disp_n = KN[1::2]
    
    plt.plot(disp_p, force_p, '-o', color=colors[i], label=f'{labels[i]} - Tension')
    plt.plot(disp_n, force_n, '--o', color=colors[i], label=f'{labels[i]} - Compression')

plt.axhline(0, color='gray', linewidth=0.5)
plt.axvline(0, color='gray', linewidth=0.5)
plt.xlabel('Displacement [m]')
plt.ylabel('Force [N]')
plt.title('Force-Displacement Diagrams for 4 Different Inelastic Springs')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True)
plt.tight_layout()
plt.show()

#%% Calculate Stiffness and Damping for Each Element
K_elastic = []
C = []
print("Element Properties:")
print("=" * 50)
for i, KP in enumerate(KP_list):
    # Calculate elastic stiffness (initial slope)
    K_elastic.append(KP[0] / KP[1])  # [N/m] Elastic Stiffness
    
    # Calculate natural frequency and damping coefficient
    omega = np.sqrt(K_elastic[i] / M[i])
    C.append(2 * zi[i] * omega * M[i])  # Damping Coefficient
    
    print(f"Element {i+1}:")
    print(f"  Elastic Stiffness (K): {K_elastic[-1]:.2f} N/m")
    print(f"  Natural Frequency (ω): {omega:.4f} rad/s")
    print(f"  Damping Coefficient (C): {C[-1]:.2f} N·s/m")
    print("-" * 30)

#%% Individual plots for each element
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
axes = axes.flatten()

for i, (KP, KN) in enumerate(zip(KP_list, KN_list)):
    force_p = KP[0::2]
    disp_p = KP[1::2]
    force_n = KN[0::2]
    disp_n = KN[1::2]
    
    axes[i].plot(disp_p, force_p, 'r-o', label='Tension')
    axes[i].plot(disp_n, force_n, 'b-o', label='Compression')
    axes[i].axhline(0, color='gray', linewidth=0.5)
    axes[i].axvline(0, color='gray', linewidth=0.5)
    axes[i].set_xlabel('Displacement [m]')
    axes[i].set_ylabel('Force [N]')
    axes[i].set_title(f'Element {i+1}: {labels[i]}')
    axes[i].legend()
    axes[i].grid(True)

plt.tight_layout()
plt.show()
#%%-------------------------------------------------------------------------------
#%% DEFINE ANALYSIS PROPERTIES
MAX_ITERATIONS = 20000    # Convergence iteration for test
MAX_TOLERANCE = 1.0e-10   # Convergence tolerance for test
#SPRING_KIND: 1 -> 'ELASTIC'
#SPRING_KIND: 2 -> 'INELASTIC'

#%%-------------------------------------------------------------------------------
def FREE_VIBRATION_SDOF(IA, IU, IV, SPRING_KIND, M, C, K_elastic, KP, KN, u0, v0, a0, duration, dt):
    # Create model
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)
    
    # Add nodes
    ops.node(1, 0.0)
    ops.node(2, 0.0)
    ops.node(3, 0.0)
    ops.node(4, 0.0)
    ops.node(5, 0.0)
    
    # Fix node 1
    ops.fix(1, 1)
    
    # Define material and element
    if SPRING_KIND == 'ELASTIC': 
        ops.uniaxialMaterial('Elastic', 11, K_elastic[0], C[0])
        ops.uniaxialMaterial('Elastic', 22, K_elastic[1], C[1])
        ops.uniaxialMaterial('Elastic', 33, K_elastic[2], C[2])
        ops.uniaxialMaterial('Elastic', 44, K_elastic[3], C[3])
        # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/ZeroLength.html
        ops.element('zeroLength', 1, 1, 2,'-mat', 11,'-dir', 1) # ELEMENT 01
        ops.element('zeroLength', 2, 2, 3,'-mat', 22,'-dir', 1) # ELEMENT 02
        ops.element('zeroLength', 3, 3, 4,'-mat', 33,'-dir', 1) # ELEMENT 03
        ops.element('zeroLength', 4, 4, 5,'-mat', 44,'-dir', 1) # ELEMENT 04
    if SPRING_KIND == 'INELASTIC': 
        ops.uniaxialMaterial('HystereticSM', 11, '-posEnv', *KP_1.flatten(), '-negEnv', *KN_1.flatten(), '-pinch', 1, 1)
        ops.uniaxialMaterial('HystereticSM', 22, '-posEnv', *KP_2.flatten(), '-negEnv', *KN_2.flatten(), '-pinch', 1, 1)
        ops.uniaxialMaterial('HystereticSM', 33, '-posEnv', *KP_3.flatten(), '-negEnv', *KN_3.flatten(), '-pinch', 1, 1)
        ops.uniaxialMaterial('HystereticSM', 44, '-posEnv', *KP_4.flatten(), '-negEnv', *KN_4.flatten(), '-pinch', 1, 1)
        #INFO LINK: https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/HystereticSM.html
        ops.uniaxialMaterial('Viscous', 311, C[0], 1.0)  # Material for C (alpha=1.0 for linear)
        ops.uniaxialMaterial('Viscous', 322, C[1], 1.0)  # Material for C (alpha=1.0 for linear)
        ops.uniaxialMaterial('Viscous', 333, C[2], 1.0)  # Material for C (alpha=1.0 for linear)
        ops.uniaxialMaterial('Viscous', 344, C[3], 1.0)  # Material for C (alpha=1.0 for linear)
        ops.element('zeroLength', 1, 1, 2,'-mat', 11, 311,'-dir', 1, 1) # ELEMENT 01
        ops.element('zeroLength', 2, 2, 3,'-mat', 22, 322,'-dir', 1, 1) # ELEMENT 02
        ops.element('zeroLength', 3, 3, 4,'-mat', 33, 333,'-dir', 1, 1) # ELEMENT 03
        ops.element('zeroLength', 4, 4, 5,'-mat', 44, 344,'-dir', 1, 1) # ELEMENT 04
    
    
    # Define masses
    ops.mass(2, M[0])
    ops.mass(3, M[1])
    ops.mass(4, M[2])
    ops.mass(5, M[3])
    
    # Define analysis
    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.system('BandGeneral')
    ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/normDispIncr.html
    #ops.integrator('CentralDifference')  # JUST FOR LINEAR ANALYSIS - INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/centralDifference.html
    alpha=0.5; beta=0.5;
    ops.integrator('Newmark', alpha, beta) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/newmark.html
    #alpha=2/3;gamma=1.5-alpha; gamma=1.5-alpha;beta=(2-alpha)**2/4;
    #ops.integrator('HHT', alpha, gamma, beta) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/hht.html
    ops.algorithm('Newton')  # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/algorithm.html
    ops.analysis('Transient') # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/analysis.html
    
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(5, 1.0)

    if IU == True:
        # Define initial displacment
        ops.setNodeDisp(5, 1, u0, '-commit') # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/setNodeDisp.html
    if IV == True:
        # Define initial velocity
        ops.setNodeVel(5, 1, v0, '-commit')  # INFO LINK: https://openseespydoc.readthedocs.io/en/stable/src/setNodeVel.html
    if IA == True:
        # Define initial  acceleration
        ops.setNodeAccel(5, 1, a0, '-commit') # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/setNodeAccel.html
    
    # Perform analysis
    time = []
    disp = []
    vel = []
    accel = []
    reaction = []
    stiffness = []
    OMEGA, PERIOD = [], []
    PERIOD_MIN, PERIOD_MAX = [], []
    # Initialize lists for each node's displacement
    node_displacements = {
        2: [],  # DISP02
        3: [],  # DISP03
        4: [],  # DISP04
        5: []   # DISP05
    }
    
    stable = 0
    current_time = 0.0
        
    while stable == 0 and current_time < duration:
        ops.analyze(1, dt)
        S01.ANALYSIS(stable, 1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS
        current_time = ops.getTime()
        time.append(current_time)
        ops.reactions()
        reaction.append(ops.nodeReaction(1, 1)) # SHEAR BASE REACTION
        disp.append(ops.nodeDisp(5, 1))         # DISPLACEMENT NODE 05 
        vel.append(ops.nodeVel(5, 1))           # VELOCITY NODE 05
        accel.append(ops.nodeAccel(5, 1))       # ACCELERATION NODE 05
        stiffness.append(np.abs(reaction[-1] / disp[-1]))
        OMEGA.append(np.sqrt(stiffness[-1]/np.sum(M)))
        PERIOD.append((np.pi * 2) / OMEGA[-1])
        # Store displacements
        for node_id in node_displacements.keys():
            node_displacements[node_id].append(ops.nodeDisp(node_id, 1))
        # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
        PERIODmin, PERIODmax = S04.EIGENVALUE_ANALYSIS(4, PLOT=True)
        PERIOD_MIN.append(PERIODmin)
        PERIOD_MAX.append(PERIODmax)
        #print(time[-1], disp[-1], vel[-1])
    
    
    # CALCULATING MEAN PERIOD CALCULATED FROM SIMPLE FORMULA
    deltaT_mean = np.mean(PERIOD)
    #deltaT_median = np.median(PERIOD)
    print(f"Min. Period:        {np.min(PERIOD):.8e} [s]") 
    print(f"Mean Period:        {deltaT_mean:.8e} [s]")   
    print(f"Max. Period:        {np.max(PERIOD):.8e} [s]")   

    # Calculating Damping Ratio and Period Using Logarithmic Decrement Analysis 
    damping_ratio = S03.DAMPING_RATIO(disp)   
    
    return time, reaction, disp, vel, accel, stiffness, PERIOD, damping_ratio, deltaT_mean, np.array(PERIOD_MIN), np.array(PERIOD_MAX), node_displacements

#%%-------------------------------------------------------------------------------
# Analysis Durations for Dynamic Analysis:
starttime = TI.process_time()

IU = True        # Free Vibration with Initial Displacement
IV = True        # Free Vibration with Initial Velocity
IA = True        # Free Vibration with Initial Acceleration
SPRING_KIND = 'ELASTIC'
DATA = FREE_VIBRATION_SDOF(IA, IU, IV, SPRING_KIND, M, C, K_elastic, KP, KN, u0, v0, a0, duration, dt)
time, reactionE, dispE, velE, accelE, stiffnessE, periodE, E_damping_ratioE, E_periodE, period_minE, period_maxE, node_displacementsE = DATA
S02.PERIOD_FUN(dispE, dt)

SPRING_KIND = 'INELASTIC'
DATA = FREE_VIBRATION_SDOF(IA, IU, IV, SPRING_KIND, M, C, K_elastic, KP, KN, u0, v0, a0, duration, dt)
time, reactionI, dispI, velI, accelI, stiffnessI, periodI, E_damping_ratioI, E_periodI, period_minP, period_maxP, node_displacementsP = DATA
S02.PERIOD_FUN(dispI, dt)

totaltime = TI.process_time() - starttime
print(f'\nTotal Analysis Durations (s): {totaltime:.4f} \n\n')
#%%-------------------------------------------------------------------------------
# PLOT STRUCTURAL PERIOD DURING THE ANALYSIS
plt.figure(0, figsize=(12, 8))
plt.plot(time, period_minE)
plt.plot(time, period_maxE)
plt.plot(time, period_minP)
plt.plot(time, period_maxP)
plt.title('Period of Structure')
plt.ylabel('Structural Period [s]')
plt.xlabel('Time [s]')
#plt.semilogy()
plt.grid()
plt.legend([f'ELASTIC PERIOD - MIN VALUES: Min: {np.min(period_minE):.3f} (s) - Mean: {np.mean(period_minE):.3f} (s) - Max: {np.max(period_minE):.3f} (s)', 
            f'ELASTIC PERIOD - MAX VALUES:  Min: {np.min(period_maxE):.3f} (s) - Mean: {np.mean(period_maxE):.3f} (s) - Max: {np.max(period_maxE):.3f} (s)',
            f'INELASTIC PERIOD - MIN VALUES: Min: {np.min(period_minP):.3f} (s) - Mean: {np.mean(period_minP):.3f} (s) - Max: {np.max(period_minP):.3f} (s)', 
            f'INELASTIC PERIOD - MAX VALUES:  Min: {np.min(period_maxP):.3f} (s) - Mean: {np.mean(period_maxP):.3f} (s) - Max: {np.max(period_maxP):.3f} (s)',
            ])
plt.show()
#%%-------------------------------------------------------------------------------
# Plot Results
plt.figure(2, figsize=(18, 14))

# Define colors
elastic_color = 'black'
inelastic_color = 'red'

# Reaction plot
plt.subplot(6, 1, 1)
plt.plot(time, reactionE, color=elastic_color, linewidth=1.5, label=f'Elastic Damping Ratio: {100*E_damping_ratioE:.3e} %')
plt.plot(time, reactionI, color=inelastic_color, linewidth=1.5, label=f'Inelastic Damping Ratio: {100*E_damping_ratioI:.3e} %')
plt.title('Reaction Forces vs Time', fontsize=12, pad=10)
plt.ylabel('Reaction (N)', fontsize=10)
plt.grid(alpha=0.3)
plt.legend(loc='upper right', framealpha=1)

# Displacement plot
plt.subplot(6, 1, 2)
plt.plot(time, dispE, color=elastic_color, linewidth=1.5)
plt.plot(time, dispI, color=inelastic_color, linewidth=1.5)
plt.title('Displacement vs Time', fontsize=12, pad=10)
plt.ylabel('Displacement (m)', fontsize=10)
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
plt.plot(time, periodE, color=elastic_color, linewidth=1.5, label=f'Elastic Period: {E_periodE:.3e}')
plt.plot(time, periodI, color=inelastic_color, linewidth=1.5, label=f'Inelastic Period: {E_periodI:.3e}')
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
plt.xlabel('Displacement [m]')
plt.ylabel('Base-reaction [N]')
plt.title(f'Displacement vs Base-reaction')
plt.legend(['ELASTIC', 'INELASTIC'])
plt.grid()
plt.show()

#%%-------------------------------------------------------------------------------
# Plotting Nodal Displacements
def PLOT_DISPLAEMENTS(time_steps, displacements_dict, TITLE):
    plt.figure(figsize=(10, 6))
    
    for node_id, disp_values in displacements_dict.items():
        plt.plot(time_steps, disp_values, label=f'Node {node_id} - MAX. ABS. : {np.max(np.abs(disp_values)): 0.4e}', linewidth=2)
    
    plt.xlabel('Time [s]')
    plt.ylabel('Displacement [m]')
    plt.title(TITLE)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

# Usage
PLOT_DISPLAEMENTS(time, node_displacementsE, TITLE = "Node Displacements vs Time for Elastic Structure") # ELASTIC STRUCTURE
PLOT_DISPLAEMENTS(time, node_displacementsP, TITLE = "Node Displacements vs Time for Inelastic Structure") # INELASTIC STRUCTURE
#%%-------------------------------------------------------------------------------
# Print out the state of all nodes
ops.printModel("node",1, 2, 3, 4, 5)
# Print out the state of all elements
ops.printModel("ele", 1, 2, 3, 4, 5)
# Print the Model
#printModel()
ops.printModel("-JSON", "-file", "FREE-VIBRATION_U0_VO_MDF.json")
#%%-------------------------------------------------------------------------------
