######################################################################################################################
#                          >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                          #
#    SIMULATE THE FREE-VIBRATION RESPONSE OF A SINGLE-DEGREE-OF-FREEDOM (SDOF) STRUCTURE INCORPORATING A PARALLEL    #
#                    CONTACT/GAP MECHANISM TO MODEL STAGED STIFFNESS ACTIVATION USING OPENSEES                       #
#--------------------------------------------------------------------------------------------------------------------#
#                              THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                            #
#                                       EMAIL: salar.d.ghashghaei@gmail.com                                          #
######################################################################################################################
"""
This script simulates the nonlinear dynamic response (Free-vibration Analysis) of a single-degree-of-freedom system with a
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
Performs free-vibration analysis of a Single Degree of Freedom (SDOF)
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
import DAMPING_RATIO_FUN as S05
import EIGENVALUE_ANALYSIS_FUN as S06
#%%------------------------------------------------------------------------------------------------
# Define parameters (units: m, N)
FY = 85000.0                                     # [N] Yield Force of Structure
FU = 1.5 * FY                                    # [N] Ultimate Force of Structure
Ke = 4500000.0                                   # [N/m] Spring Elastic Stiffness
DY = FY / Ke                                     # [m] Yield Displacement
DSU = 0.36                                       # [m] Ultimate Displacement
Ksh = (FU - FY) / (DSU - DY)                     # [N/m] Displacement Hardening Modulus
Kp = FU / DSU                                    # [N/m] Spring Plastic Stiffness
b = Ksh / Ke                                     # Displacement Hardening Ratio

M = 50000.0       # [kg] Mass
zi = 0.05         # Damping ratio
duration = 10.0   # [s] Analysis duration
dt = 0.005         # [s] Time step
#%%------------------------------------------------------------------------------------------------
# DEFINE ANALYSIS PROPERTIES
MAX_ITERATIONS = 20000    # Convergence iteration for test
MAX_TOLERANCE = 1.0e-10   # Convergence tolerance for test
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
# Calculate Over Strength Coefficient (Ω0)
Omega_0 = FU / FY
# Calculate Displacement Ductility Ratio (μ)
mu = DSU / DY
# Calculate Ductility Coefficient (Rμ)
R_mu = (2 * mu - 1) ** 0.5 / mu ** 0.5
# Calculate Structural Behavior Coefficient (R)
R = Omega_0 * R_mu
print(f'Over Strength Coefficient (Ω0):      {Omega_0:.2f}')
print(f'Displacement Ductility Ratio (μ):    {mu:.2f}')
print(f'Ductility Coefficient (Rμ):          {R_mu:.2f}')
print(f'Structural Behavior Coefficient (R): {R:.2f}')
#%%    
omega = np.sqrt(Ke/M)    # [N/m] Elastic Stiffness
C = 2 * zi * omega * M   # Damping Coefficient
#%% ---------------------------------------------
def FREE_VIBRATION_SDOF(IA, IU, IV, SPRING, u0, v0, a0, duration, dt, GAP_DIST):
    #%% Create model
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)
    
    #%% Add nodes
    ops.node(1, 0.0)
    ops.node(2, 0.0)
    ops.node(3, 0.0) # Contact Node
    
    #%% Fix node 1
    ops.fix(1, 1)
    ops.fix(3, 1) # Contact Node
    
    #%% Define material and element
    MatTag01, MatTag02 = 1, 2
    if SPRING == 'ELASTIC': 
        ops.uniaxialMaterial('Elastic', MatTag01, Ke, C)
        ops.element('zeroLength', 1, 1, 2,'-mat', MatTag01,'-dir', 1) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/ZeroLength.html
    if SPRING == 'INELASTIC':
        ops.uniaxialMaterial('HystereticSM', MatTag01, '-posEnv', *KP.flatten(), '-negEnv', *KN.flatten(), '-pinch', 1, 1)
        #INFO LINK: https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/HystereticSM.html
        ops.uniaxialMaterial('Viscous', MatTag02, C, 1.0)  # Material for C (alpha=1.0 for linear)
        ops.element('zeroLength', 1, 1, 2,'-mat', MatTag01, MatTag02,'-dir', 1, 1)
    
    
    # Define mass to node 2
    ops.mass(2, M)
    
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
    ops.load(2, 1.0)

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
    force1 = []; force2 = [];
    PERIOD_MIN, PERIOD_MAX = [], []
    
    stable = 0
    current_time = 0.0
    CONTACT_ADDED = False  
    
    while stable == 0 and current_time < duration:
        ops.analyze(1, dt)
        S01.ANALYSIS(stable, 1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS
        current_time = ops.getTime()
        time.append(current_time)
        #ops.reactions()
        #reaction.append(ops.nodeReaction(1, 1)) # BASE REACTION
        disp.append(ops.nodeDisp(2, 1))   # DISPLACEMENT  
        vel.append(ops.nodeVel(2, 1))     # VELOCITY
        accel.append(ops.nodeAccel(2, 1)) # ACCELERATION
        # Get force in Spring 1 (element 1)
        f1 = ops.eleResponse(1, 'force')[1]    # queries 'force' of zeroLength
        force1.append(f1) # Internal Force for element 1
        
        #%% Check activation condition
        if (disp[-1] >= GAP_DIST) and not CONTACT_ADDED:
            # Create Spring 2 in parallel (element tag 2)
            if SPRING == 'ELASTIC':
                ops.element('zeroLength', 2, 3, 2, '-mat', MatTag01, '-dir', 1)
                ops.domainChange()    # update the solver with new element
            if SPRING == 'INELASTIC':    
                ops.element('zeroLength', 2, 3, 2, '-mat', MatTag01, MatTag02, '-dir', 1, 1)
                ops.domainChange()    # update the solver with new element
            # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/domainChange.html
            CONTACT_ADDED = True
            
        #%% If Spring 2 is active, record its force
        if CONTACT_ADDED:
            f2 = ops.eleResponse(2, 'force')[0]
        else:
            f2 = 0.0
            
        force2.append(f2) # Internal Force for element 2
        reaction.append(np.abs((f1+f2))) # Total Reaction Force
        stiffness.append(np.abs(reaction[-1] / disp[-1]))
        OMEGA.append(np.sqrt(stiffness[-1]/M))
        PERIOD.append((np.pi * 2) / OMEGA[-1])
        print(time[-1], disp[-1], vel[-1])
        #%% IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
        PERIODmin, PERIODmax = S06.EIGENVALUE_ANALYSIS(1, PLOT=True)
        PERIOD_MIN.append(PERIODmin)
        PERIOD_MAX.append(PERIODmax)
        
    # Mean Period 
    deltaT_mean = np.mean(PERIOD)
    #deltaT_median = np.median(PERIOD)
    print(f"Min. Period:        {np.min(PERIOD):.8e} [s]") 
    print(f"Mean Period:        {deltaT_mean:.8e} [s]")   
    print(f"Max. Period:        {np.max(PERIOD):.8e} [s]")   

    # Calculate Damping Ratio of SDOF the system
    damping_ratio = S05.DAMPING_RATIO(disp) 
    
    # Compute modal properties
    ops.modalProperties("-print", "-file", "SALAR_ModalReport.txt", "-unorm") 
    
    DATA = (time, reaction, disp, vel, accel, stiffness,
            PERIOD, deltaT_mean, damping_ratio,
            force1, force2,
            np.array(PERIOD_MIN), np.array(PERIOD_MAX))
    
    return DATA

#%% ---------------------------------------------
# Analysis Durations for Dynamic Analysis:
starttime = TI.process_time()

u0 = -0.031       # [m] Initial displacement
v0 = 0.015        # [m/s] Initial velocity
a0 = 0.0065       # [m/s^2] Initial acceleration
IU = True         # Free Vibration with Initial Displacement
IV = True         # Free Vibration with Initial Velocity
IA = True         # Free Vibration with Initial Acceleration
GAP_DIST = 0.025  # [m] Gap Distance (Conact Distance)

SPRING_KIND = 'ELASTIC'
DATA = FREE_VIBRATION_SDOF(IA, IU, IV, SPRING_KIND, u0, v0, a0, duration, dt, GAP_DIST)
(time, reactionE, dispE, velE, accelE, stiffnessE,
 periodE, E_periodE, E_damping_ratioE,
 force1E, force2E,
 period_minE, period_maxE) = DATA

#S02.PERIOD_FUN(dispE, dt)

SPRING_KIND = 'INELASTIC'
DATA = FREE_VIBRATION_SDOF(IA, IU, IV, SPRING_KIND, u0, v0, a0, duration, dt, GAP_DIST)
(time, reactionI, dispI, velI, accelI, stiffnessI,
 periodI, E_periodI, E_damping_ratioI,
 force1I, force2I,
 period_minI, period_maxI) = DATA

#S02.PERIOD_FUN(dispI, dt)

totaltime = TI.process_time() - starttime
print(f'\nTotal Analysis Durations (s): {totaltime:.4f} \n\n')
#%% ---------------------------------------------
plt.figure(1, figsize=(12, 8))
plt.plot(dispE, force1E, label='Spring 1', color='black', linewidth=4)
plt.plot(dispE, force2E, label='Spring 2', color='red', linewidth=4)
plt.axvline(GAP_DIST, color='purple', label='Gap Distance', linestyle='--', linewidth=4)
plt.xlabel('Displacement (Node 2) [m]')
plt.ylabel('Spring Force [N]')
plt.title('Elastic Elements Force vs Displacement During Free-Vibration Analysis')
plt.legend(); plt.grid(True); plt.show()
#%% ---------------------------------------------
plt.figure(2, figsize=(12, 8))
plt.plot(dispI, force1I, label='Spring 1', color='black', linewidth=4)
plt.plot(dispI, force2I, label='Spring 2', color='red', linewidth=4)
plt.axvline(GAP_DIST, color='purple', label='Gap Distance', linestyle='--', linewidth=4)
plt.xlabel('Displacement (Node 2) [m]')
plt.ylabel('Spring Force [N]')
plt.title('Inelastic Elements Force vs Displacement During Free-Vibration Analysis')
plt.legend(); plt.grid(True); plt.show()
#%% ---------------------------------------------
# PLOT STRUCTURAL PERIOD DURING THE ANALYSIS
plt.figure(0, figsize=(12, 8))
plt.plot(dispE, period_minE, color='black', linewidth=4)
plt.plot(dispE, period_maxE, color='red', linewidth=4)
plt.plot(dispI, period_minI, color='green', linewidth=4)
plt.plot(dispI, period_maxI, color='purple', linewidth=4)
plt.title('Period of Structure vs Displacement During Free-Vibration Analysis')
plt.ylabel('Structural Period [s]')
plt.xlabel('Displacement [m]')
#plt.semilogy()
plt.grid()
plt.legend([f'ELASTIC PERIOD - MIN VALUES: Min: {np.min(period_minE):.3f} (s) - Mean: {np.mean(period_minE):.3f} (s) - Max: {np.max(period_minE):.3f} (s)', 
            f'ELASTIC PERIOD - MAX VALUES:  Min: {np.min(period_maxE):.3f} (s) - Mean: {np.mean(period_maxE):.3f} (s) - Max: {np.max(period_maxE):.3f} (s)',
            f'INELASTIC PERIOD - MIN VALUES: Min: {np.min(period_minI):.3f} (s) - Mean: {np.mean(period_minI):.3f} (s) - Max: {np.max(period_minI):.3f} (s)', 
            f'INELASTIC PERIOD - MAX VALUES:  Min: {np.min(period_maxI):.3f} (s) - Mean: {np.mean(period_maxI):.3f} (s) - Max: {np.max(period_maxI):.3f} (s)',
            ])
plt.show()
#%% ---------------------------------------------
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
#%% ---------------------------------------------
# Print out the state of nodes 1 and 2
ops.printModel("node",1, 2)
# Print out the state of element 1
ops.printModel("ele", 1)
# Print the Model
#printModel()
ops.printModel("-JSON", "-file", "CONTACT_PROBLEM_SDOF_FREE-VIBRATION_U0_VO_A0.json")
#%%-------------------------------------------------------------------------------
