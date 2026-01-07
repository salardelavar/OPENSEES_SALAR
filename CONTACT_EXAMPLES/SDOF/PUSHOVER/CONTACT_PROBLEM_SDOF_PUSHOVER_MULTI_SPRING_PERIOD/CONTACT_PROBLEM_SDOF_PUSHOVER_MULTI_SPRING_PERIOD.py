######################################################################################################################
#                          >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                          #
#         CONTACT-DRIVEN PUSHOVER ANALYSIS OF INELASTIC SDOF SYSTEMS: MONITORING PERIOD SHIFTS DURING                #
#                                        SECONDARY SPRING ACTIVATION IN OPENSEES                                     #
#--------------------------------------------------------------------------------------------------------------------#
#                              THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                            #
#                                       EMAIL: salar.d.ghashghaei@gmail.com                                          #
######################################################################################################################
"""
This script simulates the nonlinear pushover response of a single-degree-of-freedom system with a
 contact/gap mechanism. The structure has a primary spring (elastic or hysteretic) that activates
 immediately, while a secondary parallel spring engages only when displacement exceeds a specified 
 gap distance. This models structural components that come into contact only after certain
 deformation thresholds, such as gap-opening in masonry infills, pounding between adjacent structures,
 or secondary bracing systems activating during strong seismic events.

The analysis tracks force-displacement response, stiffness degradation, and period elongation
 as damage accumulates. The eigenvalue analysis at each step captures how the natural period
 increases with structural softening, a critical indicator of seismic vulnerability during
 progressive damage. Contact activation causes a sudden stiffness increase when the gap closes,
 followed by further period evolution as the system yields.
"""
import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
import time as TI
import ANALYSIS_FUNCTION as S01
import PERIOD_FUN as S02
import EIGENVALUE_ANALYSIS_FUN as S04
import BILINEAR_CURVE as S05
#%%-------------------------------------------------------------------------------
# Define parameters (units: m, N)
# Force–Displacement Relationship for Each Column
FY = 85000.0                                     # [N] Yield Force of Structure
FU = 1.18 * FY                                   # [N] Ultimate Force of Structure
Ke = 4500000.0                                   # [N/m] Spring Elastic Stiffness
DY = FY / Ke                                     # [m] Yield Displacement
DSU = 0.36                                       # [m] Ultimate Displacement
Ksh = (FU - FY) / (DSU - DY)                     # [N/m] Displacement Hardening Modulus
Kp = FU / DSU                                    # [N/m] Spring Plastic Stiffness
b = Ksh / Ke                                     # Displacement Hardening Ratio 

XC = [0.0, 3.4, 8.1, 11.3]                     # [m] Distance of each column
M = [150000.0, 100000.0, 250000.0, 350000.0]   # [kg] Mass
zi = [0.05, 0.01, 0.02, 0.03]                  # Damping ratio
Xcm = (XC[0]*M[0] + XC[1]*M[1] + XC[2]*M[2] + XC[3]*M[3])/np.sum(M)
print(f'Center of Mass: {Xcm:.3f} [m]')
#exit

duration = 10.0                    # [s] Analysis duration
dt = 0.01                          # [s] Time step
#%%-------------------------------------------------------------------------------
#%% DEFINE ANALYSIS PROPERTIES
MAX_ITERATIONS = 20000    # Convergence iteration for test
MAX_TOLERANCE = 1.0e-10   # Convergence tolerance for test
#SPRING_KIND: 1 -> 'ELASTIC'
#SPRING_KIND: 2 -> 'INELASTIC'

#%%-------------------------------------------------------------------------------
#%% Inelastic Force-Displacement Parameters for 4 Different Inelastic Springs 

# Element 1 - Primary Spring (Original parameters)
stiffness_factor_1 = 1.0
FY1_1, DY1_1 = FY * stiffness_factor_1, DY               # First yield point
FY2_1, DY2_1 = FU * stiffness_factor_1, DSU              # Peak force
FY3_1, DY3_1 = 0.8*FU * stiffness_factor_1, 1.10*DSU     # Post-peak softening
FY4_1, DY4_1 = 0.8*FU * stiffness_factor_1, 1.40*DSU     # Plateau
FY5_1, DY5_1 = 0.5*FU * stiffness_factor_1, 1.45*DSU     # Further softening
FY6_1, DY6_1 = 0.3*FU * stiffness_factor_1, 1.70*DSU     # Near-zero stiffness
FY7_1, DY7_1 = 0.0*FU * stiffness_factor_1, 1.95*DSU     # Zero force (failure)

KP_1 = np.array([FY1_1, DY1_1, FY2_1, DY2_1, FY3_1, DY3_1, FY4_1, DY4_1, 
                FY5_1, DY5_1, FY6_1, DY6_1, FY7_1, DY7_1])
KE01 = FY1_1/DY1_1
# Compression for Element 1
FY1n_1, DY1n_1 = -FY * stiffness_factor_1, -DY              # First yield in compression
FY2n_1, DY2n_1 = -FU * stiffness_factor_1, -DSU             # Peak compressive force
FY3n_1, DY3n_1 = -0.3*FU * stiffness_factor_1, -1.35*DSU    # Post-peak softening
KN_1 = np.array([FY1n_1, DY1n_1, FY2n_1, DY2n_1, FY3n_1, DY3n_1])

# Element 2 - Stiffer Spring (50% stiffer)
stiffness_factor_2 = 1.5
FY1_2, DY1_2 = FY * stiffness_factor_2, DY               # First yield point
FY2_2, DY2_2 = FU * stiffness_factor_2, DSU              # Peak force
FY3_2, DY3_2 = 0.8*FU * stiffness_factor_2, 1.10*DSU     # Post-peak softening
FY4_2, DY4_2 = 0.8*FU * stiffness_factor_2, 1.40*DSU     # Plateau
FY5_2, DY5_2 = 0.5*FU * stiffness_factor_2, 1.45*DSU     # Further softening
FY6_2, DY6_2 = 0.3*FU * stiffness_factor_2, 1.70*DSU     # Near-zero stiffness
FY7_2, DY7_2 = 0.0*FU * stiffness_factor_2, 1.95*DSU     # Zero force (failure)

KP_2 = np.array([FY1_2, DY1_2, FY2_2, DY2_2, FY3_2, DY3_2, FY4_2, DY4_2, 
                FY5_2, DY5_2, FY6_2, DY6_2, FY7_2, DY7_2])

# Compression for Element 2
FY1n_2, DY1n_2 = -FY * stiffness_factor_2, -DY              # First yield in compression
FY2n_2, DY2n_2 = -FU * stiffness_factor_2, -DSU             # Peak compressive force
FY3n_2, DY3n_2 = -0.3*FU * stiffness_factor_2, -1.35*DSU    # Post-peak softening
KN_2 = np.array([FY1n_2, DY1n_2, FY2n_2, DY2n_2, FY3n_2, DY3n_2])

# Element 3 - Soft Spring (30% softer)
stiffness_factor_3 = 0.7
FY1_3, DY1_3 = FY * stiffness_factor_3, DY               # First yield point
FY2_3, DY2_3 = FU * stiffness_factor_3, DSU              # Peak force
FY3_3, DY3_3 = 0.8*FU * stiffness_factor_3, 1.10*DSU     # Post-peak softening
FY4_3, DY4_3 = 0.8*FU * stiffness_factor_3, 1.40*DSU     # Plateau
FY5_3, DY5_3 = 0.5*FU * stiffness_factor_3, 1.45*DSU     # Further softening
FY6_3, DY6_3 = 0.3*FU * stiffness_factor_3, 1.70*DSU     # Near-zero stiffness
FY7_3, DY7_3 = 0.0*FU * stiffness_factor_3, 1.95*DSU     # Zero force (failure)

KP_3 = np.array([FY1_3, DY1_3, FY2_3, DY2_3, FY3_3, DY3_3, FY4_3, DY4_3, 
                FY5_3, DY5_3, FY6_3, DY6_3, FY7_3, DY7_3])

# Compression for Element 3
FY1n_3, DY1n_3 = -FY * stiffness_factor_3, -DY              # First yield in compression
FY2n_3, DY2n_3 = -FU * stiffness_factor_3, -DSU             # Peak compressive force
FY3n_3, DY3n_3 = -0.3*FU * stiffness_factor_3, -1.35*DSU    # Post-peak softening
KN_3 = np.array([FY1n_3, DY1n_3, FY2n_3, DY2n_3, FY3n_3, DY3n_3])

# Element 4 - Softer Spring (60% softer)
stiffness_factor_4 = 0.4
FY1_4, DY1_4 = FY * stiffness_factor_4, DY               # First yield point
FY2_4, DY2_4 = FU * stiffness_factor_4, DSU              # Peak force
FY3_4, DY3_4 = 0.8*FU * stiffness_factor_4, 1.10*DSU     # Post-peak softening
FY4_4, DY4_4 = 0.8*FU * stiffness_factor_4, 1.40*DSU     # Plateau
FY5_4, DY5_4 = 0.5*FU * stiffness_factor_4, 1.45*DSU     # Further softening
FY6_4, DY6_4 = 0.3*FU * stiffness_factor_4, 1.70*DSU     # Near-zero stiffness
FY7_4, DY7_4 = 0.0*FU * stiffness_factor_4, 1.95*DSU     # Zero force (failure)

KP_4 = np.array([FY1_4, DY1_4, FY2_4, DY2_4, FY3_4, DY3_4, FY4_4, DY4_4, 
                FY5_4, DY5_4, FY6_4, DY6_4, FY7_4, DY7_4])

# Compression for Element 4
FY1n_4, DY1n_4 = -FY * stiffness_factor_4, -DY              # First yield in compression
FY2n_4, DY2n_4 = -FU * stiffness_factor_4, -DSU             # Peak compressive force
FY3n_4, DY3n_4 = -0.3*FU * stiffness_factor_4, -1.35*DSU    # Post-peak softening
KN_4 = np.array([FY1n_4, DY1n_4, FY2n_4, DY2n_4, FY3n_4, DY3n_4])

# Store all elements in lists for easy access
KP_list = np.multiply([KP_1, KP_2, KP_3, KP_4], 1.0)
KN_list = np.multiply([KN_1, KN_2, KN_3, KN_4], 1.0)

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
def PUSHOVER_SDOF(DINCR, DMAX, SPRING_KIND, M, C, K_elastic, KP, KN, duration, dt, GAP_DIST): 
    G = 9.81  # [m/s^2] Acceleration due to Gravity
    # Create model
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)
    
    # Add nodes
    ops.node(1, 0.0)
    ops.node(2, 0.0)
    
    # Fix node 1
    ops.fix(1, 1)
    
    # Add Contact Nodes
    ops.node(100, 0.0)
    
    # Fix Contact Node
    ops.fix(100, 1)
    
    # Define material and element
    if SPRING_KIND == 'ELASTIC': 
        ops.uniaxialMaterial('Elastic', 11, K_elastic[0], C[0])
        ops.uniaxialMaterial('Elastic', 22, K_elastic[1], C[1])
        ops.uniaxialMaterial('Elastic', 33, K_elastic[2], C[2])
        ops.uniaxialMaterial('Elastic', 44, K_elastic[3], C[3])
        # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/ElasticUni.html
        ops.element('zeroLength', 1, 1, 2,'-mat', 11, 22, 33, 44,'-dir', 1, 1, 1, 1) # ELEMENT 01 TO 04
    
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
        ops.element('zeroLength', 1, 1, 2,
                    '-mat', 11, 311, 22, 322, 33, 333, 44, 344,
                    '-dir', 1, 1, 1, 1, 1, 1, 1, 1) # ELEMENT 01 TO 04

    # Define masses
    ops.mass(2, M[0] + M[1] + M[2] + M[3])
    #INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/mass.html
        
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(2, (M[0] + M[1] + M[2] + M[3])*G)
    
    # Define analysis
    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.system('BandGeneral')
    ops.test('EnergyIncr', MAX_TOLERANCE, MAX_ITERATIONS) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/normDispIncr.html
    ops.algorithm('Newton')  # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/algorithm.html
    IDctrlDOF = 1   ## INCREMENTAL DISPLACEMENT DOF IN X DIRECTION
    IDctrlNODE = 2  ## INCREMENTAL DISPLACEMENT NODE IN X DIRECTION
    ops.integrator('DisplacementControl', IDctrlNODE, IDctrlDOF, DINCR)
    ops.analysis('Static')
    
    # Perform analysis
    time = []
    disp = []
    reaction = []
    stiffness = []
    OMEGA, PERIOD = [], []
    force1 = []; force2 = [];
    PERIOD_MIN, PERIOD_MAX = [], []
        
    # Initialize lists for each node's displacement
    node_displacements = {
        2: [],  # DISP02
    }
    
    node_reactions = {
        2: [],  # REACTION 02
    }
    
    current_time = 0.0
    Nsteps =  int(np.abs(DMAX/ DINCR))  
    
    CONTACT_ADDED = False
    
    for step in range(Nsteps):
        OK = ops.analyze(1, dt)
        S01.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS
        current_time = ops.getTime()
        time.append(current_time)
        #ops.reactions()
        #reaction.append(ops.nodeReaction(1, 1)) # BASE REACTION
        disp.append(ops.nodeDisp(2, 1))         # DISPLACEMENT NODE 02 
        
        if disp[-1] >= GAP_DIST and not CONTACT_ADDED: # CHECK FOR CONTACT
            print(f'IN STEP {step} CONTACT DONE!')
            if SPRING_KIND == 'ELASTIC':
                # CONTACT FORCE-DISPALCEMENT RELATION FOR ELASTIC SPRING
                ops.uniaxialMaterial('Elastic', 110, K_elastic[0], C[0])
                ops.uniaxialMaterial('Elastic', 220, K_elastic[1], C[1])
                ops.uniaxialMaterial('Elastic', 330, K_elastic[2], C[2])
                ops.uniaxialMaterial('Elastic', 440, K_elastic[3], C[3])
                #ops.domainChange()
                # CONTACT ELEMENT FOR ELASTIC SPRING
                ops.element('zeroLength', 2, 100, 2,
                            '-mat', 110, 220, 330, 440,
                            '-dir', 1, 1, 1, 1) # ELEMENT 01 TO 04
                ops.domainChange() # update the solver with new element
                CONTACT_ADDED = True
            if SPRING_KIND == 'INELASTIC':
                # CONTACT FORCE-DISPALCEMENT RELATION FOR INELASTIC SPRING
                ops.uniaxialMaterial('HystereticSM', 110, '-posEnv', *KP_1.flatten(), '-negEnv', *KN_1.flatten(), '-pinch', 1, 1)
                ops.uniaxialMaterial('HystereticSM', 220, '-posEnv', *KP_2.flatten(), '-negEnv', *KN_2.flatten(), '-pinch', 1, 1)
                ops.uniaxialMaterial('HystereticSM', 330, '-posEnv', *KP_3.flatten(), '-negEnv', *KN_3.flatten(), '-pinch', 1, 1)
                ops.uniaxialMaterial('HystereticSM', 440, '-posEnv', *KP_4.flatten(), '-negEnv', *KN_4.flatten(), '-pinch', 1, 1)
                ops.uniaxialMaterial('Viscous', 3110, C[0], 1.0)  # Material for C (alpha=1.0 for linear)
                ops.uniaxialMaterial('Viscous', 3220, C[1], 1.0)  # Material for C (alpha=1.0 for linear)
                ops.uniaxialMaterial('Viscous', 3330, C[2], 1.0)  # Material for C (alpha=1.0 for linear)
                ops.uniaxialMaterial('Viscous', 3440, C[3], 1.0)  # Material for C (alpha=1.0 for linear)
                #ops.domainChange()
                # CONTACT ELEMENT FOR INELASTIC SPRING
                ops.element('zeroLength', 2, 100, 2,
                            '-mat', 110, 3110, 220, 3220, 330, 3330, 440, 3440,
                            '-dir', 1, 1, 1, 1, 1, 1, 1, 1) # ELEMENT 01 TO 04
                ops.domainChange() # update the solver with new element
                CONTACT_ADDED = True   
            CONTACT_ADDED = True 
        
        f1 = ops.eleResponse(1, 'force')[1]    # queries 'force' of zeroLength
        force1.append(f1) # Internal Force for element 1
        #%% If Contact is active, record its force
        if CONTACT_ADDED:
            f2 = ops.eleResponse(2, 'force')[0]
        else:
            f2 = 0.0
        force2.append(f2) # Internal Force for element 2
            
        reaction.append(np.abs((f1+f2))) # Total Reaction Force    
        # Store displacements and reacions
        for node_id in node_displacements.keys():
            node_displacements[node_id].append(ops.nodeDisp(node_id, 1))
            node_reactions[node_id].append(reaction[-1])  # Reaction force
        
        stiffness.append(np.abs(reaction[-1]) / np.abs(disp[-1])) # Structural Stiffness During Analysis  
        OMEGA.append(np.sqrt(stiffness[-1]/np.sum(M)))
        PERIOD.append((np.pi * 2) / OMEGA[-1])
        # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
        PERIODmin, PERIODmax = S04.EIGENVALUE_ANALYSIS(1, PLOT=True)
        PERIOD_MIN.append(PERIODmin)
        PERIOD_MAX.append(PERIODmax)
        
        print(step+1, disp[-1], reaction[-1])
    else:
        print('Analysis is Finish Successfully \n\n')
    
    # CALCULATING MEAN PERIOD CALCULATED FROM SIMPLE FORMULA
    deltaT_mean = np.mean(PERIOD)  
    
    # Compute modal properties
    ops.modalProperties("-print", "-file", "SALAR_ModalReport.txt", "-unorm")
    
    DATA = (time, reaction, disp, stiffness, PERIOD,
            deltaT_mean, np.array(PERIOD_MIN), np.array(PERIOD_MAX),
            node_displacements, node_reactions, force1, force2)
    
    return DATA 

#%%-------------------------------------------------------------------------------
# Analysis Durations for Static Analysis:
starttime = TI.process_time()

DMAX = 0.35            # [m] Max. Pushover Incremental Displacement
DINCR = 0.001          # [m] Pushover Increment
GAP_DIST = 0.25        # [m] Gap Distance (Conact Distance)

SPRING_KIND = 'ELASTIC'
DATA = PUSHOVER_SDOF(DINCR, DMAX, SPRING_KIND, M, C, K_elastic, KP, KN, duration, dt, GAP_DIST)
time, reactionE, dispE, stiffnessE, periodE, E_periodE, period_minE, period_maxE, node_displacementsE, node_reactionsE, force1E, force2E = DATA
S02.PERIOD_FUN(dispE, dt)

SPRING_KIND = 'INELASTIC'
DATA = PUSHOVER_SDOF(DINCR, DMAX, SPRING_KIND, M, C, K_elastic, KP, KN, duration, dt, GAP_DIST)
time, reactionI, dispI, stiffnessI, periodI, E_periodI, period_minP, period_maxP, node_displacementsP, node_reactionsP, force1P, force2P = DATA
S02.PERIOD_FUN(dispI, dt)

totaltime = TI.process_time() - starttime
print(f'\nTotal Analysis Durations (s): {totaltime:.4f} \n\n')
#%%-------------------------------------------------------------------------------
# PLOT STRUCTURAL PERIOD DURING THE ANALYSIS
plt.figure(0, figsize=(12, 8))
plt.plot(100*np.array(dispE)/DSU, period_minE)
plt.plot(100*np.array(dispE)/DSU, period_maxE)
plt.plot(100*np.array(dispI)/DSU, period_minP)
plt.plot(100*np.array(dispI)/DSU, period_maxP)
plt.title('Period of Structure vs Damage level During Pushover Analysis')
plt.ylabel('Structural Period [s]')
plt.xlabel('Damage level [%]')
#plt.xlabel('Displacement [m]')
#plt.semilogy()
plt.grid()
plt.legend([f'ELASTIC PERIOD - MIN VALUES: Min: {np.min(period_minE):.3f} (s) - Mean: {np.mean(period_minE):.3f} (s) - Max: {np.max(period_minE):.3f} (s)', 
            f'ELASTIC PERIOD - MAX VALUES:  Min: {np.min(period_maxE):.3f} (s) - Mean: {np.mean(period_maxE):.3f} (s) - Max: {np.max(period_maxE):.3f} (s)',
            f'INELASTIC PERIOD - MIN VALUES: Min: {np.min(period_minP):.3f} (s) - Mean: {np.mean(period_minP):.3f} (s) - Max: {np.max(period_minP):.3f} (s)', 
            f'INELASTIC PERIOD - MAX VALUES:  Min: {np.min(period_maxP):.3f} (s) - Mean: {np.mean(period_maxP):.3f} (s) - Max: {np.max(period_maxP):.3f} (s)',
            ])
plt.show()

plt.figure(1, figsize=(12, 8))
plt.plot(stiffnessE, period_minE)
plt.plot(stiffnessE, period_maxE)
plt.plot(stiffnessI, period_minP)
plt.plot(stiffnessI, period_maxP)
plt.title('Period of Structure vs Structural Stiffness During Pushover Analysis')
plt.ylabel('Structural Period [s]')
plt.xlabel('Structural Stiffness [N/m]')
#plt.semilogy()
plt.grid()
plt.legend([f'ELASTIC PERIOD - MIN VALUES: Min: {np.min(period_minE):.3f} (s) - Mean: {np.mean(period_minE):.3f} (s) - Max: {np.max(period_minE):.3f} (s)', 
            f'ELASTIC PERIOD - MAX VALUES:  Min: {np.min(period_maxE):.3f} (s) - Mean: {np.mean(period_maxE):.3f} (s) - Max: {np.max(period_maxE):.3f} (s)',
            f'INELASTIC PERIOD - MIN VALUES: Min: {np.min(period_minP):.3f} (s) - Mean: {np.mean(period_minP):.3f} (s) - Max: {np.max(period_minP):.3f} (s)', 
            f'INELASTIC PERIOD - MAX VALUES:  Min: {np.min(period_maxP):.3f} (s) - Mean: {np.mean(period_maxP):.3f} (s) - Max: {np.max(period_maxP):.3f} (s)',
            ])
plt.show()
#%%-------------------------------------------------------------------------------
plt.figure(1, figsize=(12, 8))
plt.plot(dispE, force1E, label='Spring 1', color='black', linewidth=4)
plt.plot(dispE, force2E, label='Spring 2', color='red', linewidth=4)
plt.axvline(GAP_DIST, color='purple', label='Gap Distance', linestyle='--', linewidth=4)
plt.xlabel('Displacement (Node 2) [m]')
plt.ylabel('Spring Force [N]')
plt.title('Elements Force vs Displacement During Elastic Pushover Analysis')
plt.legend(); plt.grid(True); plt.show()
#%%-------------------------------------------------------------------------------
plt.figure(2, figsize=(12, 8))
plt.plot(dispI, force1P, label='Spring 1', color='black', linewidth=4)
plt.plot(dispI, force2P, label='Spring 2', color='red', linewidth=4)
plt.axvline(GAP_DIST, color='purple', label='Gap Distance', linestyle='--', linewidth=4)
plt.xlabel('Displacement (Node 2) [m]')
plt.ylabel('Spring Force [N]')
plt.title('Elements Force vs Displacement During Inelastic Pushover Analysis')
plt.legend(); plt.grid(True); plt.show()
#%%-------------------------------------------------------------------------------
plt.figure(3, figsize=(8, 6))
plt.plot(np.abs(dispE), np.abs(reactionE), color='black', linewidth=2)
plt.plot(np.abs(dispI), np.abs(reactionI), color='purple', linewidth=2)
plt.xlabel('Displacement [m]')
plt.ylabel('Base-reaction [N]')
plt.title(f'Displacement vs Base-reaction')
plt.legend(['ELASTIC', 'INELASTIC'])
plt.grid()
plt.show()

#%%-------------------------------------------------------------------------------
# Plotting Nodal Base-Reaction and Displacements
def PLOT_BASEREACTION_DISPLACEMENT(displacements_dict, reactions_dict, TITLE):
    plt.figure(figsize=(10, 6))
    
    for node_id, disp_values in displacements_dict.items():
        # Check if node exists in reactions dictionary and both arrays have same length
        if node_id in reactions_dict and len(disp_values) == len(reactions_dict[node_id]):
            plt.plot(disp_values, reactions_dict[node_id], 
                    label=f'Node {node_id} - MAX. ABS. : {np.max(np.abs(disp_values)): 0.4e}', 
                    linewidth=2)
        else:
            print(f"Warning: Node {node_id} not found in reactions dict or array length mismatch")
    
    plt.xlabel('Displacement [m]')
    plt.ylabel('Base Reaction [N]')
    plt.title(TITLE)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()


PLOT_BASEREACTION_DISPLACEMENT(node_displacementsE, node_reactionsE, TITLE="Node Displacements vs Base Reactions for Elastic Structure") # ELASTIC STRUCTURE
PLOT_BASEREACTION_DISPLACEMENT(node_displacementsP, node_reactionsP, TITLE="Node Displacements vs Base Reactions for Inelastic Structure") # INELASTIC STRUCTURE
#%%------------------------------------------------------------------------------

# --------------------------------------
#  Plot BaseShear-Displacement Analysis 
# --------------------------------------
XX = np.abs(dispI); YY = np.abs(reactionI); # ABSOLUTE VALUE
SLOPE_NODE = 10

DATA = S05.BILNEAR_CURVE(XX, YY, SLOPE_NODE)
X, Y, Elastic_ST, Plastic_ST, Tangent_ST, Ductility_Rito, Over_Strength_Factor = DATA

XLABEL = 'Displacement in X [m]'
YLABEL = 'Base-Shear Reaction [N]'
LEGEND01 = 'Curve'
LEGEND02 = 'Bilinear Fitted'
LEGEND03 = 'Undefined'
TITLE = f'Last Data of BaseShear-Displacement Analysis - Ductility Ratio: {X[2]/X[1]:.4f} - Over Strength Factor: {Y[2]/Y[1]:.4f}'
COLOR = 'black'
S05.PLOT_2D(np.abs(dispI), np.abs(reactionI), X, Y, X, Y, XLABEL, YLABEL, TITLE, LEGEND01, LEGEND02, LEGEND03, COLOR='black', Z=2) 
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
#%%-------------------------------------------------------------------------------
# Print out the state of all nodes
ops.printModel("node",1, 2)
# Print out the state of all elements
ops.printModel("ele", 1)
# Print the Model
#printModel()
ops.printModel("-JSON", "-file", "CONTACT_PROBLEM_SDOF_PUSHOVER_MULTI_SPRING_PERIOD.json")
#%%-------------------------------------------------------------------------------
