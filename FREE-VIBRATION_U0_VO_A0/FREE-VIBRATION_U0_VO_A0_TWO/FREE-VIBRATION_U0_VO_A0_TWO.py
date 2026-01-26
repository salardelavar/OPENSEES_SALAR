######################################################################################################################
#                          >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                          #
#                               FREE-VIBRATION ANALYSIS OF SDOF STRUCTURE USING OPENSEES                             #
#--------------------------------------------------------------------------------------------------------------------#
#                              THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                            #
#                                       EMAIL: salar.d.ghashghaei@gmail.com                                          #
######################################################################################################################
"""
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
#%%------------------------------------------------------------------------------------------------
# Define parameters (units: m, N)
FY = 85000.0                         # [N] Yield Force of Structure
FU = 1.5 * FY                        # [N] Ultimate Force of Structure
Ke = 4500000.0                       # [N/m] Spring Elastic Stiffness
DY = FY / Ke                         # [m] Yield Displacement
DSU = 0.36                           # [m] Ultimate Displacement
Ksh = (FU - FY) / (DSU - DY)         # [N/m] Displacement Hardening Modulus
Kp = FU / DSU                        # [N/m] Spring Plastic Stiffness
b = Ksh / Ke                         # Displacement Hardening Ratio

M = 50000.0                          # [kg] Mass
zi = 0.05                            # Damping ratio
u0 = -0.044                          # [m] Initial displacement
v0 = 0.015                           # [m/s] Initial velocity
a0 = 0.0065                          # [m/s^2] Initial acceleration
duration = 10.0                      # [s] Analysis duration
dt = 0.01                            # [s] Time step
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
def FREE_VIBRATION_SDOF(IA, IU, IV, SPRING_KIND, u0, v0, a0, duration, dt):
    # Create model
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)
    
    # Add nodes
    ops.node(1, 0.0)
    ops.node(2, 0.0)
    
    # Fix node 1
    ops.fix(1, 1)
    
    # Define material and element
    if SPRING_KIND == 'ELASTIC': 
        ops.uniaxialMaterial('Elastic', 1, Ke, C)
        ops.element('zeroLength', 1, 1, 2,'-mat', 1,'-dir', 1) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/ZeroLength.html
    if SPRING_KIND == 'INELASTIC':
        ops.uniaxialMaterial('HystereticSM', 1, '-posEnv', *KP.flatten(), '-negEnv', *KN.flatten(), '-pinch', 1, 1)
        #INFO LINK: https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/HystereticSM.html
        ops.uniaxialMaterial('Viscous', 2, C, 1.0)  # Material for C (alpha=1.0 for linear)
        ops.element('zeroLength', 1, 1, 2,'-mat', 1, 2,'-dir', 1, 1)
    
    
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
    DI = []
    stable = 0
    current_time = 0.0
        
    while stable == 0 and current_time < duration:
        ops.analyze(1, dt)
        S01.ANALYSIS(stable, 1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS
        current_time = ops.getTime()
        time.append(current_time)
        ops.reactions()
        reaction.append(ops.nodeReaction(1, 1)) # BASE REACTION
        disp.append(ops.nodeDisp(2, 1))   # DISPLACEMENT  
        vel.append(ops.nodeVel(2, 1))     # VELOCITY
        accel.append(ops.nodeAccel(2, 1)) # ACCELERATION
        stiffness.append(np.abs(reaction[-1]) / np.abs(disp[-1]))
        OMEGA.append(np.sqrt(stiffness[-1]/M))
        PERIOD.append((np.pi * 2) / OMEGA[-1])
        DI.append(100*(np.abs(disp[-1])-DY)/(DSU-DY)) # DAMAGE INDEX
        #print(time[-1], disp[-1], vel[-1])
    
    # Mean Period 
    deltaT_mean = np.mean(PERIOD)
    #deltaT_median = np.median(PERIOD)
    print(f"Min. Period:        {np.min(PERIOD):.8e} [s]") 
    print(f"Mean Period:        {deltaT_mean:.8e} [s]")   
    print(f"Max. Period:        {np.max(PERIOD):.8e} [s]")   

    # Calculating Damping Ratio and Period Using Logarithmic Decrement Analysis
    displacement = np.array(disp)
    damping_ratio = S05.DAMPING_RATIO(displacement)   
    
    return time, reaction, disp, vel, accel, stiffness, PERIOD, damping_ratio, deltaT_mean, DI

#%% ---------------------------------------------
# Analysis Durations for Dynamic Analysis:
starttime = TI.process_time()

IU = True        # Free Vibration with Initial Displacement
IV = True        # Free Vibration with Initial Velocity
IA = True        # Free Vibration with Initial Acceleration
SPRING_KIND = 'ELASTIC'
time, reactionE, dispE, velE, accelE, stiffnessE, periodE, E_damping_ratioE, E_periodE, diE = FREE_VIBRATION_SDOF(IA, IU, IV, SPRING_KIND, u0, v0, a0, duration, dt)
S02.PERIOD_FUN(dispE, dt)

SPRING_KIND = 'INELASTIC'
time, reactionI, dispI, velI, accelI, stiffnessI, periodI, E_damping_ratioI, E_periodI, diI = FREE_VIBRATION_SDOF(IA, IU, IV, SPRING_KIND, u0, v0, a0, duration, dt)
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
plt.plot(time, reactionE, color=elastic_color, linewidth=1.5, label=f'Elastic Damping Ratio: {E_damping_ratioE:.3e} %')
plt.plot(time, reactionI, color=inelastic_color, linewidth=1.5, label=f'Inelastic Damping Ratio: {E_damping_ratioI:.3e} %')
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

plt.figure(4, figsize=(8, 6))
plt.plot(dispE, diE, color='black', linewidth=2)
plt.plot(dispI, diI, color='purple', linewidth=2)
plt.xlabel('Displacement [m]')
plt.ylabel('Structural Damage Index [%]')
plt.title(f'Displacement vs Structural Damage Index')
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
ops.printModel("-JSON", "-file", "FREE-VIBRATION_U0_VO.json")
#%%-------------------------------------------------------------------------------
