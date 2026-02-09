######################################################################################################################
#                          >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                          #
# COMPARATIVE FREE-VIBRATION AND PUSHOVER ANALYSIS OF AN SDOF SYSTEM: ELASTIC VERSUS INELASTIC RESPONSE WITH PARALLEL# 
#                         SPRINGS FOR STRUCTURAL AND NON-STRUCTURAL ELEMENTS USING OPENSEES                          #
#--------------------------------------------------------------------------------------------------------------------#
#                              THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                            #
#                                       EMAIL: salar.d.ghashghaei@gmail.com                                          #
######################################################################################################################
"""
Performs free-vibration and pushover analysis of a Single Degree of Freedom (SDOF)
 structure using OpenSeesPy, comparing elastic and inelastic spring behavior. 
---------------------------------------------------------------------------- 
- Structural element stiffness refers to how beams, columns, slabs,
 and shear walls resist deformation under loads and directly control
 the building’s strength, stability, and natural period.
- It is intentionally designed and calculated, strongly influencing seismic
 forces, drift, and load redistribution.
- Non-structural element stiffness comes from components like infill walls,
 partitions, façades, and cladding that are not part of the main load-resisting system.
- Although often neglected in design, these elements can significantly increase
 initial lateral stiffness and alter dynamic behavior.
- Damage or failure of non-structural elements may change stiffness suddenly,
 affecting seismic response and causing unexpected performance issues.
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
M = [15000.0, 10000.0] # [kg] Mass of each node - STRUCTURAL AND NON-STRUCTURAL ELEMENT
zi = [0.05, 0.005]     # Damping ratio - STRUCTURAL AND NON-STRUCTURAL ELEMENT
MASS = np.sum(M)       # [kg] Total Mass of Structure
#%% DEFINE PARAMETٍERS FOR DYNAMIC ANALYSIS 
u0 = 0.035                         # [m] Initial displacement
v0 = 0.015                         # [m/s] Initial velocity
a0 = 0.0065                        # [m/s^2] Initial acceleration
duration = 20.0                    # [s] Analysis duration
dt = 0.001                         # [s] Time step

#%%-------------------------------------------------------------------------------
#%% DEFINE PARAMETٍERS FOR STATIC ANALYSIS 
DMAX = 0.045            # [m] Maximum Displacement
DINCR = 0.0001          # [m] Incremental Displacement

#%%-------------------------------------------------------------------------------
#%% DEFINE ANALYSIS PROPERTIES
MAX_ITERATIONS = 20000    # Convergence iteration for test
MAX_TOLERANCE = 1.0e-10   # Convergence tolerance for test
#SPRING_TYPE: 1 -> 'ELASTIC'
#SPRING_TYPE: 2 -> 'INELASTIC'

#%% Inelastic Force-Displacement Parameters for 4 Different Inelastic Springs

# Element 1 - Primary Spring (Original parameters)
stiffness_factor_1 = 1.0
FY1_1, DY1_1 = 2772.0 * stiffness_factor_1, 0.01    # First yield in compression
FY2_1, DY2_1 = 3104.6 * stiffness_factor_1, 0.20    # Peak compressive force
FY3_1, DY3_1 = 1663.2 * stiffness_factor_1, 0.56    # Post-peak softening
KP_1 = np.array([FY1_1, DY1_1, FY2_1, DY2_1, FY3_1, DY3_1])


# Compression for Element 1
FY1n_1, DY1n_1 = -2772.0 * stiffness_factor_1, -0.01     # First yield point
FY2n_1, DY2n_1 = -3104.6 * stiffness_factor_1, -0.05     # Peak force
FY3n_1, DY3n_1 = -1663.2 * stiffness_factor_1, -0.08     # Post-peak softening
FY4n_1, DY4n_1 = -1663.2 * stiffness_factor_1, -0.17     # Plateau
FY5n_1, DY5n_1 = -277.2 * stiffness_factor_1, -0.28      # Further softening
FY6n_1, DY6n_1 = -200.0 * stiffness_factor_1, -0.41      # Near-zero stiffness
FY7n_1, DY7n_1 = -0.0, -0.52                             # Zero force (failure)
KN_1 = np.array([FY1n_1, DY1n_1, FY2n_1, DY2n_1, FY3n_1, DY3n_1, FY4n_1, DY4n_1, 
                FY5n_1, DY5n_1, FY6n_1, DY6n_1, FY7n_1, DY7n_1])

# Element 2 - Stiffer Spring (50% stiffer)
stiffness_factor_2 = 1.5
FY2_1, DY2_1 = 2772.0 * stiffness_factor_2, 0.01
FY2_2, DY2_2 = 3104.6 * stiffness_factor_2, 0.05
FY2_3, DY2_3 = 2663.2 * stiffness_factor_2, 0.06

KP_2 = np.array([FY2_1, DY2_1, FY2_2, DY2_2, FY2_3, DY2_3])

# Compression for Element 2
FY1n_2, DY1n_2 = -2772.0 * stiffness_factor_2, -0.01
FY2n_2, DY2n_2 = -3104.6 * stiffness_factor_2, -0.05
FY3n_2, DY3n_2 = -2663.2 * stiffness_factor_2, -0.12
KN_2 = np.array([FY1n_2, DY1n_2, FY2n_2, DY2n_2, FY3n_2, DY3n_2])


# Store all elements in lists for easy access
KP_list = [KP_1, KP_2]
KN_list = [KN_1, KN_2]

#%%% Plotting Force-Displacement Diagrams for All 2 Elements
plt.figure(1, figsize=(12, 8))
colors = ['red', 'blue']
labels = ['Element 1 (Original) - STRUCTURAL ELEMENT', 'Element 2 (Stiffer) - NON-STRUCTURAL ELEMENT']

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
plt.title('Force-Displacement Diagrams for 2 Different Inelastic Springs')
#plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.legend()
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
fig, axes = plt.subplots(1, 2, figsize=(12, 10))
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
def ANALYSIS_SDOF(ANA_KIND, DMAX, DINCR, IA, IU, IV, SPRING_TYPE, M, C, K_elastic, KP, KN, u0, v0, a0, duration, dt):
    # Create model
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)
    
    # Add nodes
    ops.node(1, 0.0)
    ops.node(2, 0.0)
    
    # Fix node 1
    ops.fix(1, 1)
    
    # Define material and element
    if SPRING_TYPE == 'ELASTIC': 
        ops.uniaxialMaterial('Elastic', 11, K_elastic[0], C[0]) # STRUCTURAL ELEMENT
        ops.uniaxialMaterial('Elastic', 22, K_elastic[1], C[1]) # NON-STRUCTURAL ELEMENT
        # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/ZeroLength.html
        ops.element('zeroLength', 1, 1, 2,'-mat', 11, 22, '-dir', 1, 1) # ELEMENT 01
        
    if SPRING_TYPE == 'INELASTIC': 
        ops.uniaxialMaterial('HystereticSM', 11, '-posEnv', *KP_1.flatten(), '-negEnv', *KN_1.flatten(), '-pinch', 1, 1) # STRUCTURAL ELEMENT
        ops.uniaxialMaterial('HystereticSM', 22, '-posEnv', *KP_2.flatten(), '-negEnv', *KN_2.flatten(), '-pinch', 1, 1) # NON-STRUCTURAL ELEMENT
        #INFO LINK: https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/HystereticSM.html
        ops.uniaxialMaterial('Viscous', 311, C[0], 1.0)  # Material for C (alpha=1.0 for linear)
        ops.uniaxialMaterial('Viscous', 322, C[1], 1.0)  # Material for C (alpha=1.0 for linear)
        ops.element('zeroLength', 1, 1, 2,'-mat', 11, 311, 22, 322,'-dir', 1, 1, 1, 1) # ELEMENT 01
    
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(2, 1.0)
    
    # Define masses
    ops.mass(2, MASS)
    
    if ANA_KIND == 'PUSHOVER':
        # Total steps
        steps = int(np.abs(DMAX)/np.abs(DINCR))

        ops.constraints('Plain')
        ops.numberer('Plain')
        ops.system('BandGeneral')
        ops.algorithm('Newton')
        ops.analysis('Static')
        # Perform analysis
        STEP = []
        disp = []
        reaction = []
        stiffness = []
        PERIOD_MIN, PERIOD_MAX = [], []
        # Initialize lists for each node's displacement
        node_displacements = {
            2: [],  # DISP02
        }
        
        for step in range(steps):
            ops.integrator('DisplacementControl', 2, 1, DINCR) 
            OK = ops.analyze(1)
            S01.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS
            ops.reactions()
            reaction.append(ops.nodeReaction(1, 1)) # SHEAR BASE REACTION
            disp.append(ops.nodeDisp(2, 1))         # DISPLACEMENT NODE 05
            stiffness.append(np.abs(reaction[-1])/np.abs(disp[-1])) # LATERAL STIFFNESS IN X
            # Store displacements
            for node_id in node_displacements.keys():
                node_displacements[node_id].append(ops.nodeDisp(node_id, 1))
            # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
            PERIODmin, PERIODmax = S04.EIGENVALUE_ANALYSIS(1, PLOT=True)
            PERIOD_MIN.append(PERIODmin)
            PERIOD_MAX.append(PERIODmax)
            STEP.append(step)
            print(step+1, disp[-1], reaction[-1])  
        else:
            print('Analysis is Finish Successfully \n\n')    
            
        DATA = (STEP, reaction, disp, stiffness,
                np.array(PERIOD_MIN), np.array(PERIOD_MAX),
                node_displacements)
        
        return DATA
        
    if ANA_KIND == 'FREE_VIBRATION':         
        
        # Define analysis
        ops.constraints('Plain')
        ops.numberer('Plain')
        ops.system('BandGeneral')
        ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/normDispIncr.html
        #ops.integrator('CentralDifference')  # JUST FOR LINEAR ANALYSIS - INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/centralDifference.html
        alpha=0.5; beta=0.25;
        ops.integrator('Newmark', alpha, beta) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/newmark.html
        #alpha=2/3;gamma=1.5-alpha; gamma=1.5-alpha;beta=(2-alpha)**2/4;
        #ops.integrator('HHT', alpha, gamma, beta) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/hht.html
        ops.algorithm('Newton')  # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/algorithm.html
        ops.analysis('Transient') # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/analysis.html
        
    
    
        if IU == True:
            # Define initial displacment
            ops.setNodeDisp(2, 1, u0, '-commit') # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/setNodeDisp.html
        if IV == True:
            # Define initial velocity
            ops.setNodeVel(2, 1, v0, '-commit')  # INFO LINK: https://openseespydoc.readthedocs.io/en/stable/src/setNodeVel.html
        if IA == True:
            # Define initial  acceleration
            ops.setNodeAccel(2, 1, a0, '-commit') # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/setNodeAccel.html
        
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
            disp.append(ops.nodeDisp(2, 1))         # DISPLACEMENT NODE 02 
            vel.append(ops.nodeVel(2, 1))           # VELOCITY NODE 02
            accel.append(ops.nodeAccel(2, 1))       # ACCELERATION NODE 02
            stiffness.append(np.abs(reaction[-1] / disp[-1]))
            OMEGA.append(np.sqrt(stiffness[-1]/MASS))
            PERIOD.append((np.pi * 2) / OMEGA[-1])
            # Store displacements
            for node_id in node_displacements.keys():
                node_displacements[node_id].append(ops.nodeDisp(node_id, 1))
            # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
            PERIODmin, PERIODmax = S04.EIGENVALUE_ANALYSIS(1, PLOT=True)
            PERIOD_MIN.append(PERIODmin)
            PERIOD_MAX.append(PERIODmax)
            print(time[-1], disp[-1], vel[-1])
        else:
            print('Analysis is Finish Successfully \n\n')
        
        # CALCULATING MEAN PERIOD CALCULATED FROM SIMPLE FORMULA
        deltaT_mean = np.mean(PERIOD)
        #deltaT_median = np.median(PERIOD)
        print(f"Min. Period:        {np.min(PERIOD):.8e} [s]") 
        print(f"Mean Period:        {deltaT_mean:.8e} [s]")   
        print(f"Max. Period:        {np.max(PERIOD):.8e} [s]")   
    
        # Calculating Damping Ratio and Period Using Logarithmic Decrement Analysis 
        damping_ratio = S03.DAMPING_RATIO(disp)   
        
        # Compute modal properties
        ops.modalProperties("-print", "-file", "SALAR_ModalReport.txt", "-unorm") 
        
        DATA = (time, reaction, disp, vel, accel,
                stiffness, PERIOD, damping_ratio, deltaT_mean,
                np.array(PERIOD_MIN), np.array(PERIOD_MAX),
                node_displacements)
        
        return DATA

#%%-------------------------------------------------------------------------------
# Analysis Durations for Dynamic Analysis:
starttime = TI.process_time()
# FREE-VIBRATION
IU = True        # Free Vibration with Initial Displacement
IV = True        # Free Vibration with Initial Velocity
IA = True        # Free Vibration with Initial Acceleration
SPRING_TYPE = 'ELASTIC'
ANA_KIND = 'FREE_VIBRATION'
DATA = ANALYSIS_SDOF(ANA_KIND, DMAX, DINCR, IA, IU, IV, SPRING_TYPE, M, C, K_elastic, KP, KN, u0, v0, a0, duration, dt)
time, reactionE, dispE, velE, accelE, stiffnessE, periodE, E_damping_ratioE, E_periodE, period_minE, period_maxE, node_displacementsE = DATA
S02.PERIOD_FUN(dispE, dt)

SPRING_TYPE = 'INELASTIC'
ANA_KIND = 'FREE_VIBRATION'
DATA = ANALYSIS_SDOF(ANA_KIND, DMAX, DINCR, IA, IU, IV, SPRING_TYPE, M, C, K_elastic, KP, KN, u0, v0, a0, duration, dt)
time, reactionI, dispI, velI, accelI, stiffnessI, periodI, E_damping_ratioI, E_periodI, period_minP, period_maxP, node_displacementsP = DATA
S02.PERIOD_FUN(dispI, dt)
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
ops.printModel("-JSON", "-file", "SDOF_FREE-VIBRATION_&_PUSHOVER_STRUCTURAL&NON-STRUCTURAL.json")
#%%-------------------------------------------------------------------------------
# PUSHOVER
SPRING_TYPE = 'ELASTIC'
ANA_KIND = 'PUSHOVER'
DATA = ANALYSIS_SDOF(ANA_KIND, DMAX, DINCR, IA, IU, IV, SPRING_TYPE, M, C, K_elastic, KP, KN, u0, v0, a0, duration, dt)
step, reactionpE, disppE, stiffnessE, period_minE, period_maxE, node_displacementsE = DATA

SPRING_TYPE = 'INELASTIC'
ANA_KIND = 'PUSHOVER'
DATA = ANALYSIS_SDOF(ANA_KIND, DMAX, DINCR, IA, IU, IV, SPRING_TYPE, M, C, K_elastic, KP, KN, u0, v0, a0, duration, dt)
step, reactionpI, disppI, stiffnessI, period_minI, period_maxI, node_displacementsI = DATA

totaltime = TI.process_time() - starttime
print(f'\nTotal Analysis Durations (s): {totaltime:.4f} \n\n')

#%%-------------------------------------------------------------------------------
plt.figure(0, figsize=(12, 8))
plt.plot(dispE, reactionE, color='blue', linewidth=5.5)
plt.plot(dispI, reactionI, color='purple', linewidth=5.5)
plt.plot(disppE, reactionpE, color=elastic_color, linewidth=5.5)
plt.plot(disppI, reactionpI, color=inelastic_color, linewidth=5.5)
plt.title('Free-vibration and Pushover of Structure')
plt.ylabel('Base-raction [N]')
plt.xlabel('Displacement [m]')
plt.legend(['FREE-VIBRATION ELASTIC', 'FREE-VIBRATION INELASTIC','PUSHOVER ELASTIC', 'PUSHOVER INELASTIC'])
plt.grid()
plt.show()
#%%-------------------------------------------------------------------------------
# PLOT STRUCTURAL PERIOD DURING THE ANALYSIS
plt.figure(1, figsize=(12, 8))
plt.plot(step, period_minE)
plt.plot(step, period_maxE)
plt.plot(step, period_minI)
plt.plot(step, period_maxI)
plt.title('Period of Structure')
plt.ylabel('Structural Period [s]')
plt.xlabel('Step')
#plt.semilogy()
plt.grid()
plt.legend([f'ELASTIC PERIOD - MIN VALUES: Min: {np.min(period_minE):.3f} (s) - Mean: {np.mean(period_minE):.3f} (s) - Max: {np.max(period_minE):.3f} (s)', 
            f'ELASTIC PERIOD - MAX VALUES:  Min: {np.min(period_maxE):.3f} (s) - Mean: {np.mean(period_maxE):.3f} (s) - Max: {np.max(period_maxE):.3f} (s)',
            f'INELASTIC PERIOD - MIN VALUES: Min: {np.min(period_minP):.3f} (s) - Mean: {np.mean(period_minP):.3f} (s) - Max: {np.max(period_minP):.3f} (s)', 
            f'INELASTIC PERIOD - MAX VALUES:  Min: {np.min(period_maxP):.3f} (s) - Mean: {np.mean(period_maxP):.3f} (s) - Max: {np.max(period_maxP):.3f} (s)',
            ])
plt.show()
#%%-------------------------------------------------------------------------------
# Usage
PLOT_DISPLAEMENTS(step, node_displacementsE, TITLE = "Node Displacements vs Time for Elastic Structure") # ELASTIC STRUCTURE
PLOT_DISPLAEMENTS(step, node_displacementsI, TITLE = "Node Displacements vs Time for Inelastic Structure") # INELASTIC STRUCTURE
#%%-------------------------------------------------------------------------------
def DISSIPATED_ENERGY_FUN(disp, force):
    import numpy as np
    disp = np.array(disp)
    force = np.array(force)
    # Incremental energy
    energy = 0.0
    for i in range(1, len(disp)):
        d_disp = disp[i] - disp[i-1]
        avg_force = 0.5 * (force[i] + force[i-1])
        energy += abs(avg_force * d_disp)
    
    print(f"Dissipated Energy = {energy:.3f}")
    return energy

print('ELASTIC STRUCTURE FREE-VIBRATION DISSIPATED ENERGY: ')
energyE = DISSIPATED_ENERGY_FUN(dispE, reactionE)
print('INELASTIC STRUCTURE FREE-VIBRATION DISSIPATED ENERGY: ')
energyI = DISSIPATED_ENERGY_FUN(dispI, reactionI)
print('ELASTIC STRUCTURE PUSHOVER DISSIPATED ENERGY: ')
energypE = DISSIPATED_ENERGY_FUN(disppE, reactionpE)
print('INELASTIC STRUCTURE PUSHOVER DISSIPATED ENERGY: ')
energypI = DISSIPATED_ENERGY_FUN(disppI, reactionpI)

print('DISSIPATED ENERGY CAPACITY INDEX:')
print(f'{100*energyI/energypI: .4f} [%]')
#%%-------------------------------------------------------------------------------
def CUMULATIVE_DISSIPATED_ENERGY_FUN(disp, force, TITLE, COLOR):
    import numpy as np
    import matplotlib.pyplot as plt
    disp = np.array(disp)
    force = np.array(force)
    cum_energy = np.zeros(len(disp))
    for i in range(1, len(disp)):
        d_disp = disp[i] - disp[i-1]
        avg_force = 0.5 * (force[i] + force[i-1])
        cum_energy[i] = cum_energy[i-1] + abs(avg_force * d_disp)
    
    plt.figure()
    plt.plot(disp, cum_energy, color=COLOR)
    plt.xlabel('Displacement')
    plt.ylabel(TITLE)
    plt.grid(True)
    plt.show()

CUMULATIVE_DISSIPATED_ENERGY_FUN(dispE, reactionE, 'FREE-VIBRATION - Cumulative Dissipated Energy', 'red')
CUMULATIVE_DISSIPATED_ENERGY_FUN(dispI, reactionI, 'FREE-VIBRATION - Cumulative Dissipated Energy', 'black')

CUMULATIVE_DISSIPATED_ENERGY_FUN(disppE, reactionpE, 'PUSHOVER - Cumulative Dissipated Energy', 'brown')
CUMULATIVE_DISSIPATED_ENERGY_FUN(disppI, reactionpI, 'PUSHOVER - Cumulative Dissipated Energy', 'orange')
#%%-------------------------------------------------------------------------------
