###########################################################################################################
#                                         IN THE NAME OF ALLAH                                            #
#   DYNAMIC RESPONSE ANALYSIS OF A SINGLE-DEGREE-OF-FREEDOM (SDOF) SYSTEM WITH ACTIVE MASS DAMPER (AMD)   #
#                                    UNDER FREE VIBRATION CONDITIONS                                      #
#---------------------------------------------------------------------------------------------------------#
#                         THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                      #
#                                      EMAIL: salar.d.ghashghaei@gmail.com                                #
###########################################################################################################
""" 
Dynamic response analysis of a single-degree-of-freedom (SDOF) system equipped with an
 active mass damper (AMD) under free vibration conditions.
 The system includes an elastic structure with inherent viscous damping, enhanced by an AMD
 that actively controls vibrations using velocity feedback. The analysis leverages computational
 tools like OpenSeesPy to simulate the system's response, offering insights into vibration mitigation.
Methods to Increase Damping in Structural Systems
To increase the damping ratio of an SDOF structure in dynamic analysis, various techniques
 can be employed based on design, materials, and project needs. 
 Below are some practical approaches:

Active Control Systems:
Devices like active mass dampers (AMDs) use sensors and actuators to apply real-time counteracting forces,
 significantly reducing vibrations through velocity feedback control.


Passive Dampers:
Installing viscous dampers, friction dampers, or tuned mass dampers (TMDs) passively
 dissipates kinetic energy without requiring external power.


Increasing Internal Friction of Materials:
Selecting materials with higher damping properties (e.g., composites or concrete with additives)
 can enhance the damping ratio.


Using Damping Layers:
Adding layers of damping materials (e.g., polymers or rubber) between structural components reduces vibrations and boosts damping.


Modifying Structural Design:
Incorporating energy-dissipating connections or frictional interfaces can improve damping behavior.

In this analysis, the focus is on active control using an AMD, which adapts to dynamic conditions and provides superior vibration mitigation, especially for seismic applications.
Use of Active Mass Dampers in Seismic Design and Retrofitting of Bridges
This section explores the application of active mass dampers (AMDs) in the seismic design and retrofitting of urban bridges, emphasizing their critical role in ensuring functionality
 during and after earthquakes. 
 
Advantages of Active Mass Dampers:
Enhanced Energy Dissipation: AMDs adapt to varying dynamic loads, outperforming passive systems in vibration control.
Real-Time Control: Sensors enable AMDs to adjust forces instantly, optimizing performance across different seismic conditions.
Minimal Structural Interference: AMDs can be retrofitted into existing bridges with little modification, ideal for upgrading older structures.
High Reliability: Modern AMDs feature robust control algorithms and fail-safe mechanisms, ensuring performance under extreme events.

AMDs are typically installed between the deck and piers, working alongside other systems like sliding bearings to mimic seismic isolation behavior.
 By reducing displacement demands and enhancing energy dissipation, they enable more resilient and cost-effective bridge designs. Globally, AMDs have
 been integrated into numerous bridges, including retrofits following events like the 1989 Loma Prieta earthquake, underscoring their effectiveness in seismic protection.
Conclusion
This program demonstrates the efficacy of an active mass damper in mitigating vibrations in an SDOF system under free vibration conditions.
 By employing velocity feedback, the AMD enhances damping, reducing displacement, velocity, and acceleration responses. This approach is particularly
 valuable for the seismic protection of critical infrastructure like bridges, ensuring their functionality during and after earthquakes. 
""" 
#------------------------------------------------------------------------------------------------
import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
import time as ti
import SALAR_MATH as S01
import ANALYSIS_FUNCTION as S02

# Define parameters (units: mm, N)
Ms = 1000.0  # [kg] Mass of the main structure
Mb = 500.0   # [kg] Mass of the active mass damper
ks = 1e6     # [N/mm] Stiffness between ground and structure
Cs = 1e4     # [N-s/mm] Damping coefficient between ground and structure
kb = 15e5    # [N/mm] Stiffness between structure and damper mass
Cb = 5e3     # [N-s/mm] Damping coefficient between structure and damper mass
DR = 0.02    # Inherent Damping Ratio
u0 = 1.7     # [mm] Initial displacement applied to node 3
duration = 10.0  # [s] Total simulation duration
dt = 0.01        # [s] Time step
K_control = 1e4  # [N-s/mm] Control gain for AMD (tunable parameter)

# FOR NONLINEAR STRUCTURES
NONLINER = 'False' # 'True' active Nonlinear Behaviour for main Strucure Spring
DY = 0.001         # [mm] Yield Displacement
FY = 110.0         # [N] Yield Force
DU = 3.200         # [mm] Ultimate Displacement
FU = 423.0         # [N] Ultimate Force

if NONLINER == 'True' and np.abs(u0) > DY:
    print('\n\nStructure behaves Nonlinear.\n\n')

# Define Analysis Properties
MAX_ITERATIONS = 1000000   # Convergence iteration limit
MAX_TOLERANCE = 1.0e-10    # Convergence tolerance
#------------------------------------------------------------------------------------------------
# %% Perform dynamic response analysis of an SDOF system with an active mass damper
def ANALYSIS_SDOF_AMD():
    # Initialize OpenSees model
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)

    # Define nodes
    ops.node(1, 0.0)  # Ground node
    ops.node(2, 0.0)  # Main structure node (Ms)
    ops.node(3, 0.0)  # AMD node (Mb)

    # Fix ground node
    ops.fix(1, 1)

    # Assign masses
    ops.mass(2, Ms)
    ops.mass(3, Mb)

    if NONLINER == 'False':
        # Define elastic materials - Spring between ground and structure
        ops.uniaxialMaterial('Elastic', 1, ks)
    if NONLINER == 'True':    
        # Define hysteretic material - Spring between ground and structure
        pinchX = 0.8           # Pinching factor in X direction
        pinchY = 0.5           # Pinching factor in Y direction
        damage1 = 0.0          # Damage due to ductility
        damage2 = 0.0          # Damage due to energy
        beta = 0.1             # Stiffness degradation parameter
        ops.uniaxialMaterial('HystereticSM', 1, FY , DY, FU , DU, 0.2*FU, 1.1*DU, 0.1*FU , 1.2*DU,
                         -FY, -DY, -FU , -DU, -0.2*FU, -1.05*DU, -0.0*FU, -1.1*DU,
                         pinchX, pinchY, damage1, damage2, beta)# HystereticSM
        
    ops.uniaxialMaterial('Elastic', 3, kb)       # Spring between structure and AMD
    ops.uniaxialMaterial('Viscous', 2, Cs, 1.0)  # Damper between ground and structure
    ops.uniaxialMaterial('Viscous', 4, Cb, 1.0)  # Damper between structure and AMD
    #ops.uniaxialMaterial('BilinearOilDamper', 4, kb, Cb, 10.0, 0.)  # Damper between structure and AMD

    # Define elements
    ops.element('zeroLength', 1, 1, 2, '-mat', 1, '-dir', 1)  # Spring ks
    ops.element('zeroLength', 2, 1, 2, '-mat', 2, '-dir', 1)  # Damper Cs
    ops.element('zeroLength', 3, 2, 3, '-mat', 3, '-dir', 1)  # Spring kb
    ops.element('zeroLength', 4, 2, 3, '-mat', 4, '-dir', 1)  # Damper Cb

    # Static analysis to apply initial displacement
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(3, 1.0)
    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.system('BandGeneral')
    ops.algorithm('Newton')
    ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS, 1)
    ops.integrator('DisplacementControl', 2, 1, u0)
    ops.analysis('Static')
    ops.analyze(1)
    ops.setTime(0.0)
    ops.wipeAnalysis()
    ops.remove('loadPattern', 1)

    # Dynamic analysis setup
    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.system('BandGeneral')
    ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS)
    ops.integrator('Newmark', 0.5, 0.25)
    ops.algorithm('Newton')
    ops.analysis('Transient')

    # Rayleigh damping
    Lambda01 = ops.eigen('-fullGenLapack', 1)
    Omega01 = np.power(max(Lambda01), 0.5)
    a0 = 2 * DR * Omega01  # Mass-proportional damping
    a1 = 2 * DR / Omega01  # Stiffness-proportional damping
    ops.rayleigh(a0, a1, 0, 0)
    PERIOD = 2 * np.pi / Omega01

    # Define load pattern for AMD control force
    ops.timeSeries('Constant', 2)
    ops.pattern('Plain', 2, 2)

    # Initialize response lists
    time = []
    displacement = []
    velocity = []
    acceleration = []
    base_reaction = []
    control_force = []
    current_time = 0.0
    stable = 0

    # Dynamic analysis loop with AMD
    while stable == 0 and current_time < duration:
        # Compute control force based on structure velocity
        vel_node2 = ops.nodeVel(2, 1)
        F_control = -K_control * vel_node2

        # Apply control forces
        ops.load(2, -F_control)  # Negative force on structure
        ops.load(3, F_control)   # Positive force on AMD mass

        # Perform one time step analysis
        stable = ops.analyze(1, dt)
        S02.ANALYSIS(stable, 1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS


        current_time = ops.getTime()
        time.append(current_time)
        displacement.append(ops.nodeDisp(2, 1))
        velocity.append(ops.nodeVel(2, 1))
        acceleration.append(ops.nodeAccel(2, 1))
        # Total base reaction from ground elements
        force_spring = ops.eleResponse(1, 'force')[0]
        force_damper = ops.eleResponse(2, 'force')[0]
        total_reaction = -(force_spring + force_damper)
        base_reaction.append(total_reaction)
        control_force.append(F_control)

    ops.wipe()
    return time, displacement, velocity, acceleration, base_reaction, control_force, PERIOD
#------------------------------------------------------------------------------------------------
# Run analysis
starttime = ti.process_time()
time, displacement, velocity, acceleration, base_reaction, control_force, PERIOD = ANALYSIS_SDOF_AMD()
totaltime = ti.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n')
#------------------------------------------------------------------------------------------------
# Compute cumulative maximum absolute values
def MAX_ABS(X):
    X = np.asarray(X)
    X_MAX = np.zeros_like(X)
    X_MAX[0] = np.abs(X[0])
    for i in range(1, len(X)):
        X_MAX[i] = max(X_MAX[i-1], np.abs(X[i]))
    return X_MAX

DISP_Z = MAX_ABS(displacement)
VELO_Z = MAX_ABS(velocity)
ACCE_Z = MAX_ABS(acceleration)
BASE_Z = MAX_ABS(base_reaction)
CTRL_Z = MAX_ABS(control_force)

# Plotting function
def PLOT_DATA(x, y, y_max, xlabel, ylabel, title):
    plt.figure(figsize=(8, 6))
    plt.plot(x, y, color='blue', linewidth=2)# , label='Response'
    plt.plot(x, y_max, color='red', linewidth=2) #, label='Max Abs'
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    #plt.legend()
    plt.grid()
    plt.show()

# Plot results
PLOT_DATA(time, displacement, DISP_Z, 'Time [s]', 'Displacement [mm]', 
          f'Time vs Displacement (Structure) - MAX. ABS: {DISP_Z[-1]:.3f}')
PLOT_DATA(time, velocity, VELO_Z, 'Time [s]', 'Velocity [mm/s]', 
          f'Time vs Velocity (Structure) - MAX. ABS: {VELO_Z[-1]:.3f}')
PLOT_DATA(time, acceleration, ACCE_Z, 'Time [s]', 'Acceleration [mm/s^2]', 
          f'Time vs Acceleration (Structure) - MAX. ABS: {ACCE_Z[-1]:.3f}')
PLOT_DATA(time, base_reaction, BASE_Z, 'Time [s]', 'Base Reaction [N]', 
          f'Time vs Base Reaction - MAX. ABS: {BASE_Z[-1]:.3f}')
PLOT_DATA(time, control_force, CTRL_Z, 'Time [s]', 'Control Force [N]', 
          f'Time vs Control Force - MAX. ABS: {CTRL_Z[-1]:.3f}')
#------------------------------------------------------------------------------------------------
# Plot Base Reaction vs Displacement
plt.figure(figsize=(8, 6))
plt.plot(displacement, base_reaction, color='black', linewidth=2)
plt.xlabel('Displacement [mm]')
plt.ylabel('Base Reaction [N]')
plt.title(f'Base Reaction vs Displacement - MAX. ABS Disp: {np.max(np.abs(displacement)):.3f}')
plt.grid()
plt.show()
#------------------------------------------------------------------------------------------------
# PLOT THE TIME-HISTORY
S01.PLOT_TIME_HISTORY(time, displacement, velocity, acceleration, base_reaction)
#------------------------------------------------------------------------------------------------
# %% Run the file loading Damping Ratio Function
#exec(open("DAMPING_RATIO_FUN.py").read())
import DAMPING_RATIO_FUN as DRF
  
DISP = displacement    
solution = DRF.DAMPING_RATIO(DISP)    
print(f"Exact Damping Ratio: {100*solution[0]:.8e} (%)")
#------------------------------------------------------------------------------------------------
