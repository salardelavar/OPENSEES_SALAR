######################################################################################################################
#                          >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                          #
#        OPTIMIZATION OF A SINGLE-DEGREE-OF-FREEDOM STRUCTURE THROUGH FREE VIBRATION ANALYSIS USING OPENSEES         #
#--------------------------------------------------------------------------------------------------------------------#
#                                  OPTIMIZATION ALGORITHM: NEWTON-RAPHSON METHOD                                     #
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
import DAMPING_RATIO_FUN as S05

#%% DEFINE ANALYSIS PROPERTIES
MAX_ITERATIONS = 20000    # Convergence iteration for test
MAX_TOLERANCE = 1.0e-10   # Convergence tolerance for test
#SPRING_KIND: 1 -> 'ELASTIC'
#SPRING_KIND: 2 -> 'INELASTIC'

#%% ---------------------------------------------
def FREE_VIBRATION_SDOF(area, IA, IU, IV, SPRING_KIND, M, u0, v0, duration, dt):
    # Tension Force-Displacement Relationship (Positive values)
    FY1, DY1 = 2.7720*area, 0.01    # First yield point
    FY2, DY2 = 3.1046*area, 0.02    # Peak force
    FY3, DY3 = 1.6632*area, 0.04    # Post-peak softening
    FY4, DY4 = 1.6632*area, 0.06    # Plateau
    FY5, DY5 = 0.2772*area, 0.28    # Further softening
    FY6, DY6 = 0.2000*area, 0.41    # Near-zero stiffness
    FY7, DY7 = 0.0*area, 0.52       # Zero force (failure)
    
    KP = np.array([FY1, DY1, FY2, DY2, FY3, DY3, FY4, DY4, FY5, DY5, FY6, DY6, FY7, DY7])

    # Compression Force-Displacement Relationship (Negative values)
    FY1n, DY1n = -2.7720*area, -0.01    # First yield in compression
    FY2n, DY2n = -3.1046*area, -0.02    # Peak compressive force
    FY3n, DY3n = -1.6632*area, -0.04    # Post-peak softening

    KN = np.array([FY1n, DY1n, FY2n, DY2n, FY3n, DY3n])
    
    #%%
    K = KP[0] / KP[1]        # [N/m] Elastic Stiffness

    omega = np.sqrt(K/M)
    C = 2 * zi * omega * M   # Damping Coefficient

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
        ops.uniaxialMaterial('Elastic', 1, K, C)
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
    
    stable = 0
    current_time = 0.0
        
    while stable == 0 and current_time < duration:
        ops.analyze(1, dt)
        S01.ANALYSIS(stable, 1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS
        current_time = ops.getTime()
        time.append(current_time)
        ops.reactions()
        reaction.append(ops.nodeReaction(1, 1)) # SHEAR BASE REACTION
        disp.append(ops.nodeDisp(2, 1))   # DISPLACEMENT  
        vel.append(ops.nodeVel(2, 1))     # VELOCITY
        accel.append(ops.nodeAccel(2, 1)) # ACCELERATION
        stiffness.append(np.abs(reaction[-1] / disp[-1]))
        OMEGA.append(np.sqrt(stiffness[-1]/M))
        PERIOD.append((np.pi * 2) / OMEGA[-1])
        #print(time[-1], reaction[-1], disp[-1])
    
    # Calculating Damping Ratio and Period Using Logarithmic Decrement Analysis 
    displacement = np.array(disp)
    TIME = np.array(time)
    peaks = np.array([displacement[i] for i in range(1, len(displacement)-1) if displacement[i] > displacement[i-1] and displacement[i] > displacement[i+1]])
    # Natural logarithm
    delta = np.log(peaks[:-1] / peaks[1:])   # Damping Ratio
    # Mean Period 
    deltaT_mean = np.mean(PERIOD)
    #deltaT_median = np.median(PERIOD)
    print(f"Min. Period:        {np.min(PERIOD):.8e} [s]") 
    print(f"Mean Period:        {deltaT_mean:.8e} [s]")   
    print(f"Max. Period:        {np.max(PERIOD):.8e} [s]")   

    damping_ratio = S05.DAMPING_RATIO(disp)     
    
    return time, reaction, disp, vel, accel, stiffness, PERIOD, damping_ratio, deltaT_mean


#%%------------------------------------------------------------------------------
# FIND BEST SPRING AREA WITH STRUCTURAL BEHAVIOR COEFFICIENT OPTIMIZATION:
# Define parameters
M = 50000.0       # [kg] Mass
zi = 0.05         # Damping ratio
u0 = -0.044       # [m] Initial displacement
v0 = 0.015        # [m/s] Initial velocity
a0 = 0.0065       # [m/s²] Initial acceleration
duration = 100.0  # [s] Analysis duration
dt = 0.01         # [s] Time step
area = 100        # [m²] Cross-sectional area    

IU = True         # Free Vibration with Initial Displacement
IV = True         # Free Vibration with Initial Velocity
IA = True         # Free Vibration with Initial Acceleration

X = area           # [m] Intial Guess Spring Area   
ESP = 1e-5         # Finite difference derivative Convergence Tolerance
TOLERANCE = 1e-6   # Convergence Tolerance
RESIDUAL = 100     # Convergence Residual 
IT = 0             # Intial Iteration
ITMAX = 100000     # Max. Iteration
TARGET_D = 0.03    # [m] Target Spring Displacement

# Target Spring Displacement must less than Initial Displacement
if np.abs(TARGET_D) >= np.abs(u0): 
    print(' \n\nTarget Spring Displacement must less than Initial Displacement.')
    exit()



SPRING_KIND = 'INELASTIC'
# Analysis Durations:
starttime = TI.process_time()

### FIND THE OPTIMUM VALUE 
while (RESIDUAL > TOLERANCE):
    # X -------------------------------------------------------
    DATA  = FREE_VIBRATION_SDOF(X, IA, IU, IV, SPRING_KIND, M, u0, v0, duration, dt)
    time, reaction, disp, vel, accel, stiffness, PERIOD, damping_ratio, period = DATA
    SUPPLY = np.max(np.abs(disp))
    print(f' SUPPLY: {SUPPLY:.8f}')
    F = SUPPLY - TARGET_D
    print('F: ', F)
    # Xmin -------------------------------------------------------
    XMIN = X - ESP  
    DATA  = FREE_VIBRATION_SDOF(XMIN, IA, IU, IV, SPRING_KIND, M, u0, v0, duration, dt)
    time, reaction, disp, vel, accel, stiffness, PERIOD, damping_ratio, period = DATA
    SUPPLYmin = np.max(np.abs(disp))
    Fmin = SUPPLYmin - TARGET_D
    print('Fmin: ', Fmin)
    # Xmax -------------------------------------------------------
    XMAX = X + ESP 
    DATA  = FREE_VIBRATION_SDOF(XMAX, IA, IU, IV, SPRING_KIND, M, u0, v0, duration, dt)
    time, reaction, disp, vel, accel, stiffness, PERIOD, damping_ratio, period = DATA
    SUPPLYmax = np.max(np.abs(disp))
    Fmax = SUPPLYmax - TARGET_D
    print('Fmax: ', Fmax)
    # DF -------------------------------------------------------
    DF = (Fmax - Fmin)/(2 * ESP);# Calculate the Finite difference derivative of F
    print('DF: ', DF)
    # DX -------------------------------------------------------
    DX = F / DF;        # Calculate dx
    print('DX: ', DX)
    # RESIDUAL -------------------------------------------------
    RESIDUAL = abs(DX); # Calculate residual
    X -= DX;            # update X
    IT += 1;            # update iteration
    print('IT: ', IT,' - RESIDUAL: ', RESIDUAL,' - X: ', X,'\n')
                
    if IT == ITMAX:
        print('\t\t Iteration reached to Max. Iteration')
        print('\t\t Change ESP and TOLERANCE for better Convergence')
        X = -1
        break;
    if RESIDUAL < TOLERANCE:
        print(f'\t\t Optimum Spring Area :                     {X:.6f}')
        print(f'\t\t Iteration Counts:                         {IT}')
        print(f'\t\t Convergence Residual:                     {RESIDUAL:.10e}')

    

totaltime = TI.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n')


#%% ---------------------------------------------
# Tension Force-Displacement Relationship (Positive values)
FY1, DY1 = 2.7720*area, 0.0001    # First yield point
FY2, DY2 = 3.1046*area, 0.0002    # Peak force
FY3, DY3 = 1.6632*area, 0.0004    # Post-peak softening
FY4, DY4 = 1.6632*area, 0.0006    # Plateau
FY5, DY5 = 0.2772*area, 0.0028    # Further softening
FY6, DY6 = 0.2000*area, 0.0041    # Near-zero stiffness
FY7, DY7 = 0.0*area, 0.0052       # Zero force (failure)
    
KP = np.array([FY1, DY1, FY2, DY2, FY3, DY3, FY4, DY4, FY5, DY5, FY6, DY6, FY7, DY7])

# Compression Force-Displacement Relationship (Negative values)
FY1n, DY1n = -2.7720*area, -0.0001    # First yield in compression
FY2n, DY2n = -3.1046*area, -0.0002    # Peak compressive force
FY3n, DY3n = -1.6632*area, -0.0004    # Post-peak softening

KN = np.array([FY1n, DY1n, FY2n, DY2n, FY3n, DY3n])

# Separate into Force and Displacement
force_p = KP[0::2]
disp_p = KP[1::2]

force_n = KN[0::2]
disp_n = KN[1::2]

#%%% Plotting Force-Displacement Diagram for Inelastic Spring
plt.figure(1, figsize=(8, 6))
plt.plot(disp_p, force_p, 'r-o', label='Tension')
plt.plot(disp_n, force_n, 'b-o', label='Compression')
plt.axhline(0, color='gray', linewidth=0.5)
plt.axvline(0, color='gray', linewidth=0.5)
plt.xlabel('Displacement')
plt.ylabel('Force')
plt.title('Force-Displacement Diagram for Inelastic Spring')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
#%% ---------------------------------------------
# Plot Results
plt.figure(2, figsize=(18, 14))

# Define colors
elastic_color = 'black'
inelastic_color = 'red'

# Reaction plot
plt.subplot(6, 1, 1)
plt.plot(time, reaction, color=inelastic_color, linewidth=1.5, label=f'Damping Ratio: {100*damping_ratio:.3e} %')
plt.title('Reaction Forces vs Time', fontsize=12, pad=10)
plt.ylabel('Reaction (N)', fontsize=10)
plt.grid(alpha=0.3)
plt.legend(loc='upper right', framealpha=1)

# Displacement plot
plt.subplot(6, 1, 2)
plt.plot(time, disp, color=inelastic_color, linewidth=1.5)
plt.title('Displacement vs Time', fontsize=12, pad=10)
plt.ylabel('Displacement (m)', fontsize=10)
plt.grid(alpha=0.3)

# Velocity plot
plt.subplot(6, 1, 3)
plt.plot(time, vel, color=inelastic_color, linewidth=1.5)
plt.title('Velocity vs Time', fontsize=12, pad=10)
plt.ylabel('Velocity (m/s)', fontsize=10)
plt.grid(alpha=0.3)

# Acceleration plot
plt.subplot(6, 1, 4)
plt.plot(time, accel, color=inelastic_color, linewidth=1.5)
plt.title('Acceleration vs Time', fontsize=12, pad=10)
plt.ylabel('Acceleration (m/s²)', fontsize=10)
plt.xlabel('Time (s)', fontsize=10)
plt.grid(alpha=0.3)

# Stiffness plot
plt.subplot(6, 1, 5)
plt.plot(time, stiffness, color=inelastic_color, linewidth=1.5)
plt.title('Stiffness vs Time', fontsize=12, pad=10)
plt.ylabel('Stiffness (N/m)', fontsize=10)
plt.semilogy()
plt.grid(alpha=0.3)

# Period plot
plt.subplot(6, 1, 6)
plt.plot(time, PERIOD, color=inelastic_color, linewidth=1.5, label=f'Structure Period: {period:.3e}')
plt.title('Period vs Time', fontsize=12, pad=10)
plt.ylabel('Period  (s)', fontsize=10)
plt.xlabel('Time (s)', fontsize=10)
plt.legend(loc='upper right', framealpha=1)
plt.semilogy()
plt.grid(alpha=0.3)

plt.tight_layout()
plt.show()

plt.figure(3, figsize=(8, 6))
plt.plot(disp, reaction, color='purple', linewidth=2)
plt.xlabel('Displacement [m]')
plt.ylabel('Base-reaction [N]')
plt.title(f'Displacement vs Base-reaction')
plt.grid()
plt.show()
#%% ---------------------------------------------
# Print out the state of nodes 1 and 2
ops.printModel("node",1, 2)
# Print out the state of element 1
ops.printModel("ele", 1)
# Print the Model
#printModel()
ops.printModel("-JSON", "-file", "FREE-VIBRATION_U0_V0_A0_OPTIMIZATION.json")
#%%-------------------------------------------------------------------------------
