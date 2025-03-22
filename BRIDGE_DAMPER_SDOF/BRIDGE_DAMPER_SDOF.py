###########################################################################################################
#                                         IN THE NAME OF ALLAH                                            #
#                    DYNAMIC RESPONSE ANALYSIS OF A SINGLE-DEGREE-OF-FREEDOM (SDOF)                       #
#                                ELASTIC STRUCTURE UNDER FREE-VIBRATION                                   #
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
import ANALYSIS_FUNCTION as S02

#------------------------------------------------------------------------------------------------
# Define parameters (units: mm, N)
# Define  Steel Material Properties for Brace Element
Ms = 1000.0  # [kg] Mass Ms
Mb = 500.0   # [kg] Mass Mb
ks = 1e6     # [N/mm] Stiffness ks
Cs = 1e4     # [N-s/mm] Damping coefficient Cs
kb = 5e5     # Stiffness kb (N/m)
Cb = 5e3     # [N-s/mm] Damping coefficient Cb
Cd = 15e3    # [N-s/mm] Damping coefficient Cd

DR = 0.02    # Damping Ratio
u0 = 0.1     # [mm] Initial displacement applied to the node 3

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
    
    # Define nodes
    ops.node(1, 0.0)  # Node for the ground
    ops.node(2, 0.0)  # Node for Ms
    ops.node(3, 0.0)  # Node for Mb
    
    # Fix node 1 (ground)
    ops.fix(1, 1)  # Fixed in the vertical direction (1D)
    
    # Assign masses
    ops.mass(2, Ms)  # Mass Ms at node 2
    ops.mass(3, Mb)  # Mass Mb at node 3
    
    # Define materials
    # Materials for springs
    ops.uniaxialMaterial('Elastic', 1, ks)  # Material for ks
    ops.uniaxialMaterial('Elastic', 3, kb)  # Material for kb
    # Materials for dampers (linear viscous damping)
    ops.uniaxialMaterial('Viscous', 2, Cs, 1.0)  # Material for Cs (alpha=1.0 for linear)
    ops.uniaxialMaterial('Viscous', 4, Cb, 1.0)  # Material for Cb (alpha=1.0 for linear)
    ops.uniaxialMaterial('Viscous', 5, Cd, 1.0)  # Material for Cd (alpha=1.0 for linear)
    
    # Define elements
    # Between node 1 and 2 (ground to Ms)
    ops.element('zeroLength', 1, 1, 2, '-mat', 1, '-dir', 1)  # Spring ks
    ops.element('zeroLength', 2, 1, 2, '-mat', 2, '-dir', 1)  # Damper Cs
    # Between node 2 and 3 (Ms to Mb)
    ops.element('zeroLength', 3, 2, 3, '-mat', 3, '-dir', 1)  # Spring kb
    ops.element('zeroLength', 4, 2, 3, '-mat', 4, '-dir', 1)  # Damper Cb
    ops.element('zeroLength', 5, 2, 3, '-mat', 5, '-dir', 1)  # Damper Cd
            
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
    
    # Dynamic analysis
    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.system('BandGeneral')
    ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS)
    ops.integrator('Newmark', 0.5, 0.25)
    ops.algorithm('Newton')
    ops.analysis('Transient')
    
    # Calculate Rayleigh damping factors
    Lambda01 = ops.eigen('-fullGenLapack', 1) # eigenvalue mode 1
    #Lambda01 = ops.eigen('-genBandArpack', 1) # eigenvalue mode 1
    Omega01 = np.power(max(Lambda01), 0.5)
    Omega02 = 2 * Omega01
    a0 = DR * (2 * Omega01 * Omega02) / (Omega01 + Omega02) # c = a0 * m : Mass-proportional damping
    a1 = DR * 2 / (Omega01 + Omega02) # c = a1 * k : Stiffness-proportional damping
    # Apply Rayleigh damping
    #ops.rayleigh(a0, a1, 0, 0)# INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/reyleigh.html
    ops.rayleigh(0, 0, 2 * DR * Omega01, 0) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/reyleigh.html
    PERIOD = np.pi / Omega01   # Structure Period   
    
    # Re-run dynamic analysis with new loading
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
# Compute the Cumulative Maximum Absolute Value of Last Analysis Data
def MAX_ABS(X):
    import numpy as np
    X = np.asarray(X)  # Convert input to a numpy array for faster operations
    X_MAX = np.zeros_like(X)  # Initialize an array to store cumulative max values
    X_MAX[0] = np.abs(X[0])  # Set the first value

    # Compute cumulative maximum absolute values
    for i in range(1, len(X)):
        X_MAX[i] = max(X_MAX[i-1], np.abs(X[i]))
    
    return X_MAX  

DISP_Z = MAX_ABS(displacement)  
VELO_Z = MAX_ABS(velocity) 
ACCE_Z = MAX_ABS(acceleration) 
BASE_Z = MAX_ABS(base_reaction) 

plt.figure(1, figsize=(8, 6))
plt.plot(time, displacement, color='blue', linewidth=2)
plt.plot(time, DISP_Z, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Displacement [mm]')
plt.title(f'Time vs Displacement - MAX. ABS: {DISP_Z[-1]}')
plt.grid()
plt.show()

plt.figure(2, figsize=(8, 6))
plt.plot(time, velocity, color='blue', linewidth=2)
plt.plot(time, VELO_Z, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Velocity [mm/s]')
plt.title(f'Time vs Velocity - MAX. ABS: {VELO_Z[-1]}')
plt.grid()
plt.show()

plt.figure(3, figsize=(8, 6))
plt.plot(time, acceleration, color='blue', linewidth=2)
plt.plot(time, ACCE_Z, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Acceleration [mm/s^2]')
plt.title(f'Time vs Acceleration - MAX. ABS: {ACCE_Z[-1]}')
plt.grid()
plt.show()

plt.figure(4, figsize=(8, 6))
plt.plot(time, base_reaction, color='blue', linewidth=2)
plt.plot(time, BASE_Z, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Base-reaction [N]')
plt.title(f'Time vs Base-reaction - MAX. ABS: {BASE_Z[-1]}')
plt.grid()
plt.show()
#------------------------------------------------------------------------------------------------
# PLOT THE TIME-HISTORY
S01.PLOT_TIME_HISTORY(time, displacement, velocity, acceleration, base_reaction)
#------------------------------------------------------------------------------------------------


