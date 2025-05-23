###########################################################################################################
#                                         IN THE NAME OF ALLAH                                            #
#  DYNAMIC RESPONSE ANALYSIS OF A SINGLE-DEGREE-OF-FREEDOM (SDOF) SYSTEM WITH VISCOUS DAMPING AND ELASTIC #
#                               STRUCTURE UNDER FREE VIBRATION CONDITIONS                                 #
#---------------------------------------------------------------------------------------------------------#
#                          THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                     #
#                                   EMAIL: salar.d.ghashghaei@gmail.com                                   #
###########################################################################################################
#
"""
To increase the damping ratio of a single degree of freedom (SDOF) structure in the dynamic analysis of structures,
 various methods and techniques can be employed, depending on the structural design, materials used, and specific project requirements.
 Below are some general and practical approaches:

1. Using Dampers:
   - Installing mechanical dampers such as viscous dampers, friction dampers, or tuned mass dampers (TMD) can significantly
   increase the damping ratio. These devices absorb and dissipate the kinetic energy of the structure.

2. Increasing Internal Friction of Materials:
   - By selecting materials with inherently higher damping properties (e.g., concrete with specific additives or composites),
   the damping ratio can be improved.

3. Modifying Structural Design:
   - Increasing the number of frictional connections or using energy-dissipating connections (e.g., semi-rigid connections) can enhance damping.

4. Using Damping Layers:
   - Adding layers of damping materials (such as polymers or rubber) between structural components can help reduce vibrations
   and increase the damping ratio.

5. Adjusting Mass and Stiffness:
   - Altering the distribution of mass and stiffness in the structure (optimizing dynamic ratios) can influence damping behavior,
   though this method tends to work more indirectly.

In practice, the damping ratio depends on the type of structure and the dynamic loading (e.g., earthquake or wind) and is typically
 determined and optimized through experimental tests or numerical simulations (e.g., modal analysis).

Use of Viscous Dampers in Seismic Design and Retrofitting of Bridges
The article explores the application of viscous dampers in the seismic design and retrofitting
 of urban bridges, emphasizing their critical role in maintaining functionality during and after earthquakes.
 Bridges are vital for transportation and emergency response, as demonstrated by past earthquakes like Loma Prieta (1989),
 Northridge (1994), and Kobe (1995), where bridge failures severely disrupted crisis management. Older bridges, designed under outdated codes,
 often lack the capacity to dissipate seismic energy, making retrofitting essential in earthquake-prone regions. Viscous dampers, often used in
 combination with sliding bearings, offer an effective solution by dissipating seismic energy, reducing input energy to the structure, and ensuring
 the bridge's primary components (deck and piers) remain elastic or near-elastic post-earthquake.

Adding viscous damping is a highly effective passive control method for bridges, providing benefits such as minimal interference with thermal
 expansion, temperature-independent performance in modern designs, and high reliability due to comprehensive testing. Viscous dampers are also durable under service loads, easy to maintain, and cost-effective over their lifecycle. Typically, nonlinear viscous dampers are used, with a force-velocity relationship defined as 
(usually 0.1 to 0.3) ensures energy dissipation even at low velocities, mimicking yielding dampers while maintaining efficiency.
 This setup reduces maximum forces in piers and foundations, lowers displacement demands on bearings, and enables more economical designs.

Viscous dampers are commonly installed between the deck and piers, often as hinged connections to prevent bending deformations and
 ensure sealing integrity. Their integration can mimic seismic isolation behavior, especially in bridges with sliding bearings, by
 reducing displacement demands and enhancing energy dissipation. Globally, numerous bridges, including the Golden Gate Bridge after
 the 1989 Loma Prieta earthquake, have been retrofitted with viscous dampers, highlighting their growing acceptance and effectiveness
 in seismic protection.
"""
#------------------------------------------------------------------------------------------------
import openseespy.opensees as ops
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time as ti
import SALAR_MATH as S01
import ANALYSIS_FUNCTION as S02

#------------------------------------------------------------------------------------------------
# Define parameters (units: mm, N)
# Define  Steel Material Properties for Brace Element
Ms = 1000.0  # [kg] Mass Ms
Mb = 500.0   # [kg] Mass Mb
ks = 1e6     # [N/mm] Stiffness ks
Cs = 1e4     # [N-s/mm] Damping coefficient Cs
kb = 5e5     # [N/m] Stiffness kb
Cb = 5e3     # [N-s/mm] Damping coefficient Cb
Cd = 15e3    # [N-s/mm] Damping coefficient Cd

DR = 0.02    # Damping Ratio
u0 = 0.3     # [mm] Initial displacement applied to the node 3

duration = 10.0  # [s] Total simulation duration
dt = 0.01        # [s] Time step

#------------------------------------------------------------------------------------------------
# Define Analysis Properties
MAX_ITERATIONS = 1000000   # Convergence iteration for test
MAX_TOLERANCE = 1.0e-10    # Convergence tolerance for test
#------------------------------------------------------------------------------------------------
### OPENSEES FUNCTION
def ANALYSIS_SDOF():
    # Initialize OpenSees model
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)
    #---------------------------------------------------------
    # Define nodes
    ops.node(1, 0.0)  # Node for the ground
    ops.node(2, 0.0)  # Node for Ms
    ops.node(3, 0.0)  # Node for Mb
    #---------------------------------------------------------
    # Fix node 1 (ground)
    ops.fix(1, 1)  # Fixed in the vertical direction (1D)
    #---------------------------------------------------------
    # Assign masses
    ops.mass(2, Ms)  # Mass Ms at node 2
    ops.mass(3, Mb)  # Mass Mb at node 3
    #---------------------------------------------------------
    # Define materials
    # Materials for springs
    
    # Bridge Column Substucture Force-Displacement Nonlinear Relation
    KP = np.array([2772.0, 0.01, 3104.6, 0.02, 1663.2, 0.04, 1663.2, 0.06, 277.2, 0.28, 200.0, 0.41, 0, 0.52])         # TENSION STRESS-STRAIN RELATION
    KN = np.array([-2772.0, -0.01, -3104.6, -0.02, -1663.2, -0.04])                                                    # COMPRESSION STRESS-STRAIN RELATION

    #ops.uniaxialMaterial('HystereticSM', 1, '-posEnv', *KP.flatten(), '-negEnv', *KN.flatten(), '-pinch', 1, 1)
    #INFO LINK: https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/HystereticSM.html
    
    # Bridge Column Substucture Force-Displacement Linear Relation
    ops.uniaxialMaterial('Elastic', 1, ks)  # Material for ks
    
    ops.uniaxialMaterial('Elastic', 3, kb)  # Material for kb
    # Materials for dampers (linear viscous damping)
    ops.uniaxialMaterial('Viscous', 2, Cs, 1.0)  # Material for Cs (alpha=1.0 for linear)
    ops.uniaxialMaterial('Viscous', 4, Cb, 1.0)  # Material for Cb (alpha=1.0 for linear)
    ops.uniaxialMaterial('Viscous', 5, Cd, 1.0)  # Material for Cd (alpha=1.0 for linear)
    #---------------------------------------------------------
    # Define elements
    # Between node 1 and 2 (ground to Ms)
    ops.element('zeroLength', 1, 1, 2, '-mat', 1, '-dir', 1)  # Spring ks
    ops.element('zeroLength', 2, 1, 2, '-mat', 2, '-dir', 1)  # Damper Cs
    # Between node 2 and 3 (Ms to Mb)
    ops.element('zeroLength', 3, 2, 3, '-mat', 3, '-dir', 1)  # Spring kb
    ops.element('zeroLength', 4, 2, 3, '-mat', 4, '-dir', 1)  # Damper Cb
    ops.element('zeroLength', 5, 2, 3, '-mat', 5, '-dir', 1)  # Damper Cd
    #---------------------------------------------------------        
    # Static analysis to apply initial displacement
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(3, 1.0)

    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.system('BandGeneral')
    ops.algorithm('Newton')
    ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS)
    ops.integrator('DisplacementControl', 3, 1, u0) # Initial displacement applied to the node 3
    ops.analysis('Static')
    ops.analyze(1)
    
    ops.setTime(0.0)
        
    # Output data
    #ops.recorder('Node', '-file', f"DTH_DYN_{i}.txt",'-time', '-node', 2, '-dof', 1, 'disp')     # Displacement Time History Node 2
    #ops.recorder('Node', '-file', f"VTH_DYN_{i}.txt",'-time', '-node', 2, '-dof', 1, 'vel')      # Velocity Time History Node 2
    #ops.recorder('Node', '-file', f"ATH_DYN_{i}.txt",'-time', '-node', 2, '-dof', 1, 'accel')    # Acceleration Time History Node 2
    #ops.recorder('Node', '-file', f"BTH_DYN_{i}.txt",'-time', '-node', 1, '-dof', 1, 'reaction') # Base Reaction Time History Node 1
    
    # Wipe analysis and reset time
    ops.wipeAnalysis()
    ops.remove('loadPattern', 1)
    #---------------------------------------------------------
    # Dynamic analysis
    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.system('BandGeneral')
    ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS)
    ops.integrator('Newmark', 0.5, 0.25)
    ops.algorithm('Newton')
    ops.analysis('Transient')
    #---------------------------------------------------------
    # Calculate Rayleigh damping factors
    Lambda01 = ops.eigen('-fullGenLapack', 1) # eigenvalue mode 1
    #Lambda01 = ops.eigen('-genBandArpack', 1) # eigenvalue mode 1
    Omega01 = np.power(max(Lambda01), 0.5)
    a0 = (2 * Omega01 * DR) / Omega01 # c = a0 * m : Mass-proportional damping
    a1 = (DR * 2) / Omega01 # c = a1 * k : Stiffness-proportional damping
    # Apply Rayleigh damping
    ops.rayleigh(a0, a1, 0, 0)# INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/reyleigh.html
    #ops.rayleigh(0, 0, 2 * DR * Omega01, 0) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/reyleigh.html
    PERIOD = np.pi / Omega01   # Structure Period   
    #---------------------------------------------------------
    # Run dynamic analysis with new loading
    time = []
    displacement = []
    velocity = []
    acceleration = []
    base_reactionW, base_reaction  = [], []
    STIFF = []
        
    stable = 0
    current_time = 0.0
    while stable == 0 and current_time < duration:
        stable = ops.analyze(1, dt)
        S02.ANALYSIS(stable, 1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS
        current_time = ops.getTime()
        time.append(current_time)
        displacement.append(ops.nodeDisp(3, 1))
        velocity.append(ops.nodeVel(3, 1))
        acceleration.append(ops.nodeAccel(3, 1))                # Structure Acceleration
        base_reactionW.append(-ops.eleResponse(1, 'force')[0])  # Reaction force for weld
        base_reaction.append(-ops.eleResponse(2, 'force')[0])   # Reaction force for brce
        STIFF.append(base_reactionW[-1] / displacement[-1])     # Siffness of Brace element
        
    ops.wipe()
    return time, displacement, velocity, acceleration, base_reaction, PERIOD, STIFF

#------------------------------------------------------------------------------------------------
# Analysis Durations:
starttime = ti.process_time()

time, displacement, velocity, acceleration, base_reaction, PERIOD, STIFF = ANALYSIS_SDOF()  

totaltime = ti.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n') 
#------------------------------------------------------------------------------------------------
# Tension force-displacement relation
KP = np.array([2772.0, 0.01, 3104.6, 0.02, 1663.2, 0.04, 1663.2, 0.06, 277.2, 0.28, 200.0, 0.41, 0, 0.52])         # TENSION STRESS-STRAIN RELATION
KPT = KP.reshape(-1, 2)  # Reshape into pairs of (stress, strain)

# Compression force-displacement relation
KN = np.array([-2772.0, -0.01, -3104.6, -0.02, -1663.2, -0.04])
KNC = KN.reshape(-1, 2)  # Reshape into pairs of (stress, strain)

# Extract force and displacement values
strain_tension = KPT[:, 1]
stress_tension = KPT[:, 0]
strain_compression = KNC[:, 1]
stress_compression = KNC[:, 0]

# Plot force-displacement relation
plt.figure(0, figsize=(8, 6))
plt.plot(strain_tension, stress_tension, 'bo-', label='Tension')
plt.plot(strain_compression, stress_compression, 'ro-', label='Compression')
plt.axhline(0, color='black', linewidth=1)
plt.axvline(0, color='black', linewidth=1)
plt.xlabel("Displacement")
plt.ylabel("Force")
plt.title("Bridge Column Substucture Force-Displacement Nonlinear Relation")
plt.legend()
plt.grid(True)
plt.show()
#------------------------------------------------------------------------------------------------
# Compute the Cumulative Maximum Absolute Value of Last Analysis Data
def MAX_ABS(X):
    X = np.asarray(X)  # Convert input to a numpy array for faster operations
    X_MAX = np.zeros_like(X)  # Initialize an array to store cumulative max values
    X_MAX[0] = np.abs(X[0])  # Set the first value

    # Compute cumulative maximum absolute values
    for i in range(1, len(X)):
        X_MAX[i] = max(X_MAX[i-1], np.abs(X[i]))
    
    return X_MAX  
#------------------------------------------------------------------------------------------------
# Plotting Function
def PLOT_DATA(x, y, y_max, xlabel, ylabel, title):
    plt.figure(figsize=(8, 6))
    plt.plot(x, y, color='blue', linewidth=2)
    plt.plot(x, y_max, color='red', linewidth=2)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid()
    plt.show()
#------------------------------------------------------------------------------------------------
# Plot Base Reaction vs Displacement
plt.figure(figsize=(8, 6))
plt.plot(displacement, base_reaction, color='black', linewidth=2)
plt.xlabel('Displacement [mm]')
plt.ylabel('Base Reaction [N]')
plt.title(f'Result for Base Reaction vs Displacement - MAX. ABS: {np.max(np.abs(displacement)):.3f}')
plt.grid()
plt.show()
#------------------------------------------------------------------------------------------------
# Calculate Cumulative Maximum Absolute Values
DISP_Z = MAX_ABS(displacement)  
VELO_Z = MAX_ABS(velocity) 
ACCE_Z = MAX_ABS(acceleration) 
BASE_Z = MAX_ABS(base_reaction) 
#------------------------------------------------------------------------------------------------
# Plot Time History Data
PLOT_DATA(time, displacement, DISP_Z, 'Time [s]', 'Displacement [mm]', f'Time vs Displacement - MAX. ABS: {DISP_Z[-1]}')
PLOT_DATA(time, velocity, VELO_Z, 'Time [s]', 'Velocity [mm/s]', f'Time vs Velocity - MAX. ABS: {VELO_Z[-1]}')
PLOT_DATA(time, acceleration, ACCE_Z, 'Time [s]', 'Acceleration [mm/s^2]', f'Time vs Acceleration - MAX. ABS: {ACCE_Z[-1]}')
PLOT_DATA(time, base_reaction, BASE_Z, 'Time [s]', 'Base-reaction [N]', f'Time vs Base-reaction - MAX. ABS: {BASE_Z[-1]}')
#------------------------------------------------------------------------------------------------
# PLOT THE TIME-HISTORY
S01.PLOT_TIME_HISTORY(time, displacement, velocity, acceleration, base_reaction)
#------------------------------------------------------------------------------------------------


