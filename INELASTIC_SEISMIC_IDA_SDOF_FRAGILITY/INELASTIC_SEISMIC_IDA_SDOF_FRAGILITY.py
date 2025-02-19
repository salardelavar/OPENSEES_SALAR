###########################################################################################################
#                                          IN THE NAME OF ALLAH                                           #
#                     FRAGILITY CURVES BASED ON INCREMENTAL DYNAMIC ANALYSIS (IDA)                        #
#                               OF A SINGLE-DEGREE-OF-FREEDOM (SDOF) SYSTEM                               #
#---------------------------------------------------------------------------------------------------------#
# This program performs Incremental Dynamic Analysis (IDA) on a Single-Degree-of-Freedom (SDOF) system    #
# subjected to seismic ground motions. The analysis evaluates the structural response under varying       #
# levels of seismic intensity.                                                                            #
# The framework is designed to support researchers and engineers in assessing the probabilistic seismic   #
# performance of structures, with a focus on understanding the impact of uncertainty on structural        #
# response and design.                                                                                    #
#---------------------------------------------------------------------------------------------------------#
# Key Features:                                                                                           #
# - Simulation of SDOF system using OpenSees.                                                             #
# - Incremental scaling of ground motions for IDA.                                                        #
# - Probabilistic fragility assessment based on predefined damage states.                                 #
# - Visualization of structural response and fragility curves.                                            #
# - Export of results for further analysis.                                                               #
#---------------------------------------------------------------------------------------------------------#
#                          THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                     #
#                                   EMAIL: salar.d.ghashghaei@gmail.com                                   #
########################################################################################################### 

import openseespy.opensees as ops
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time as ti
import SALAR_MATH as S01
import Analysis_Function as S02
from scipy.stats import lognorm

#------------------------------------------------------------------------------------------------
# Define parameters (units: m, N)
J_MAX = 50     # Incremental Dynamic Analysis number for simulation
AREA = 0.0055  # [m^2] Spring Section Area
LENGTH = 10.1  # [m] Spring Length
#------------------------------------------------------------------------------------------------
# Define  Steel Material Properties (Steel01)
fy = 0.41        # [N] Yield force of structure
fu = 1.5 * fy    # [N] Ultimate force of structure
Es = 2.1e2       # [N/m] Spring stiffness
ey = fy / Es     # [m] Yield displacement
esu = 0.36       # [m] Ultimate displacement
Esh = (fu - fy) / (esu - ey)                    # [N/m] Displacement hardening modulus
b = Esh / Es                                    # Displacement hardening ratio
K = (Es * AREA) / LENGTH                        # [N/m] Stiffness of the structure
M = 3600          # [kg] Mass of the structure
DR = 0.02         # Damping ratio

duration = 15.0  # [s] Total simulation duration
dt = 0.01        # [s] Time step
#------------------------------------------------------------------------------------------------
# Define Analysis Properties
MAX_ITERATIONS = 1000000   # Convergence iteration for test
MAX_TOLERANCE = 1.0e-10    # Convergence tolerance for test
#------------------------------------------------------------------------------------------------
### OPENSEES FUNCTION
def ANALYSIS_IDA_SDOF(j, J_MAX):
    # Initialize OpenSees model
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)
    GMfact = 9.81* ((j+1) / J_MAX) # [m/s^2] standard acceleration of gravity or standard acceleration 
        
    # Define nodes
    ops.node(1, 0.0)  # Fixed base
    ops.node(2, 0.0)  # Mass node
        
    # Define boundary conditions
    ops.fix(1, 1)
    
    # Define mass
    ops.mass(2, M)
        
    # Define material properties
    ops.uniaxialMaterial('Steel01', 1, fy, Es, b) # Steel with bilinear kinematic hardening Material
        
    # Define element
    ops.element('zeroLength', 1, 1, 2, '-mat', 1, '-dir', 1)  # DOF[1] LATERAL SPRING
    
    # Apply seismic accelerations    
    # Define time series for input motion (Acceleration time history)
    gm_accels = np.loadtxt('Ground_Acceleration_1.txt')  # Assumes acceleration in m/s²
    ops.timeSeries('Path', 1, '-dt', dt, '-values', *gm_accels.tolist(), '-factor', GMfact) # SEISMIC-X
    #ops.timeSeries('Path', 1, '-dt', dt, '-filePath', f'Ground_Acceleration_1.txt', '-factor', GMfact) # SEISMIC-X
        
    # Define load patterns
    # pattern UniformExcitation $patternTag $dof -accel $tsTag <-vel0 $vel0> <-fact $cFact>
    ops.pattern('UniformExcitation', 200, 1, '-accel', 1) # SEISMIC-X
    
    # Output data
    #ops.recorder('Node', '-file', f"DTH_DYN_{i}.txt",'-time', '-node', 2, '-dof', 1, 'disp')     # Displacement Time History Node 2
    #ops.recorder('Node', '-file', f"VTH_DYN_{i}.txt",'-time', '-node', 2, '-dof', 1, 'vel')      # Velocity Time History Node 2
    #ops.recorder('Node', '-file', f"ATH_DYN_{i}.txt",'-time', '-node', 2, '-dof', 1, 'accel')    # Acceleration Time History Node 2
    #ops.recorder('Node', '-file', f"BTH_DYN_{i}.txt",'-time', '-node', 1, '-dof', 1, 'reaction') # Base Reaction Time History Node 1
        
    # Set analysis parameters
    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.system('BandGeneral')
    ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS)
    ops.algorithm('Newton')
    ops.integrator('Newmark', 0.5, 0.25)
    ops.analysis('Transient')
    
    # Calculate Rayleigh damping factors
    Lambda01 = ops.eigen('-fullGenLapack', 1) # eigenvalue mode 1
    Omega = np.power(max(Lambda01), 0.5)
    #Lambda02 = ops.eigen('-genBandArpack', 1) # eigenvalue mode 1
    #Omega = np.power(max(min(Lambda01), min(Lambda02)), 0.5)
    betaKcomm = 2 * (DR/Omega)
    # Apply Rayleigh damping
    ops.rayleigh(0, 0, 0, betaKcomm) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/reyleigh.html
        
    # Re-run dynamic analysis with new loading
    time = []
    displacement = []
    velocity = []
    acceleration = []
    base_reaction = []
    DI = []
        
    stable = 0
    current_time = 0.0
    step = 0
    while stable == 0 and current_time < duration:
        stable = ops.analyze(1, dt)
        S02.ANALYSIS(stable, 1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS
        current_time = ops.getTime()
        time.append(current_time)
        displacement.append(ops.nodeDisp(2, 1))
        velocity.append(ops.nodeVel(2, 1))
        if step <= len(gm_accels)-1:
            acceleration.append(ops.nodeAccel(2, 1) + gm_accels[step]) # Structure Acceleration and Seismic Acceleration
        else:
            acceleration.append(ops.nodeAccel(2, 1)) # Structure Acceleration
        base_reaction.append(-ops.eleResponse(1, 'force')[0])  # Reaction force
        DI.append((displacement[-1] - ey) / (esu - ey))        # Structural Ductility Damage Index 
        step += 1
    ops.wipe()
    return time, displacement, velocity, acceleration, base_reaction, DI

#------------------------------------------------------------------------------------------------
# Analysis Durations:
starttime = ti.process_time()

# Initialize lists to store max values
max_time = []
max_displacement = []
max_velocity = []
max_acceleration = []
max_base_reaction = []
max_DI = []

# IDA ANALYSIS
for j in range(J_MAX):
    time, displacement, velocity, acceleration, base_reaction, DI = ANALYSIS_IDA_SDOF(j, J_MAX)
    # Calculate and store the max absolute values
    max_time.append(np.max(np.abs(time)))
    max_displacement.append(np.max(np.abs(displacement)))
    max_velocity.append(np.max(np.abs(velocity)))
    max_acceleration.append(np.max(np.abs(acceleration)))
    max_base_reaction.append(np.max(np.abs(base_reaction)))
    max_DI.append(np.max(np.abs(DI)))
    print(f'STEP {j + 1} DONE')              

totaltime = ti.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n') 
#------------------------------------------------------------------------------------------------
# Print the last results
print("Maximum Absolute Values Across Simulations:")
print("Time:", max_time[-1])
print("Displacement:", max_displacement[-1])
print("Velocity:", max_velocity[-1])
print("Acceleration:", max_acceleration[-1])
print("Base Reaction:", max_base_reaction[-1])
print("Ductility Damage Index:", max_DI[-1])
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
plt.ylabel('Displacement [m]')
plt.title(f'Time vs Displacement - MAX. ABS: {DISP_Z[-1]}')
plt.grid()
plt.show()

plt.figure(2, figsize=(8, 6))
plt.plot(time, velocity, color='blue', linewidth=2)
plt.plot(time, VELO_Z, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Velocity [m/s]')
plt.title(f'Time vs Velocity - MAX. ABS: {VELO_Z[-1]}')
plt.grid()
plt.show()

plt.figure(3, figsize=(8, 6))
plt.plot(time, acceleration, color='blue', linewidth=2)
plt.plot(time, ACCE_Z, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Acceleration [m/s^2]')
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
# EXPORT DATA TO EXCEL
DATA_TOTAL = {
    'Max_displacement': max_displacement,
    'Max_velocity': max_velocity,
    'Max_acceleration': max_acceleration,
    'Max_Base_Reaction': max_base_reaction,
    'Ductility_Damage_Index': max_DI
}
# Convert to DataFrame
results_df = pd.DataFrame(DATA_TOTAL)
# Export the DataFrame to an Excel file
results_df.to_excel('INELASTIC_SEISMIC_IDA_SDOF_FRAGILITY_RESULTS.xlsx', index=False)
#------------------------------------------------------------------------------------------------  
XLABEL = 'Displacement'
YLABEL = 'Base Reaction'
TITLE = f'{YLABEL} and {XLABEL} scatter chart'
COLOR = 'orange'
X = max_displacement
Y = max_base_reaction
S01.PLOT_SCATTER(X, Y , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

# CLUSTER DATA
S01.CLUSTER_DATA(X, Y, XLABEL, YLABEL, MAX_CLUSTERS=3)
#------------------------------------------------------------------------------------------------
XLABEL = 'Velocity'
YLABEL = 'Base Reaction'
TITLE = f'{YLABEL} and {XLABEL} scatter chart'
COLOR = 'cyan'
X = max_velocity
Y = max_base_reaction
S01.PLOT_SCATTER(X, Y , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

# CLUSTER DATA
S01.CLUSTER_DATA(X, Y, XLABEL, YLABEL, MAX_CLUSTERS=3)
#------------------------------------------------------------------------------------------------
XLABEL = 'Acceleration'
YLABEL = 'Base Reaction'
TITLE = f'{YLABEL} and {XLABEL} scatter chart'
COLOR = 'lime'
X = max_acceleration
Y = max_base_reaction
S01.PLOT_SCATTER(X, Y , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

# CLUSTER DATA
S01.CLUSTER_DATA(X, Y, XLABEL, YLABEL, MAX_CLUSTERS=3)
#------------------------------------------------------------------------------------------------
XLABEL = 'Displacement'
YLABEL = 'Structural Ductility Damage Index'
TITLE = f'{YLABEL} and {XLABEL} scatter chart'
COLOR = 'purple'
X = max_displacement
Y = max_DI
S01.PLOT_SCATTER(X, Y , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

# CLUSTER DATA
S01.CLUSTER_DATA(X, Y, XLABEL, YLABEL, MAX_CLUSTERS=3)
#------------------------------------------------------------------------------------------------
# PLOT THE TIME-HISTORY
S01.PLOT_TIME_HISTORY(time, displacement, velocity, acceleration, base_reaction)
#------------------------------------------------------------------------------------------------
# FRAGILITY ANALYSIS FRAMEWORK
def FRAGILITY_PROBABILITY(damage_states, im_values):
    prob_dict = {}
    for ds, (median, beta) in damage_states.items():
        prob = lognorm(s=beta, scale=median).cdf(im_values)
        prob_dict[ds] = prob
    return prob_dict   
# ----------------------------
# Fragility Assessment
# ----------------------------
# Define damage states per FEMA P-58
# INFO LINK: https://www.fema.gov/sites/default/files/documents/fema_p-58-2-se_volume2_implementation.pdf
damage_params = {
'DS1_Slight': (0.15, 0.4),    # Median PGA=0.15g, β=0.4
'DS2_Moderate': (0.30, 0.5),
'DS3_Extensive': (0.60, 0.6),
'DS4_Complete': (1.00, 0.7)
}
"""
im_values = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
probabilities = {
    'DS1': [0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99],
    'DS2': [0.0, 0.01, 0.05, 0.1, 0.3, 0.5, 0.7, 0.85, 0.95, 0.99],
    'DS3': [0.0, 0.0, 0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.95]
}
"""  
im_values = max_acceleration
probabilities = FRAGILITY_PROBABILITY(damage_params, im_values) 
# ----------------------------
# Visualization
# ----------------------------
plt.figure(figsize=(12, 6))
    
# Response plot
plt.subplot(1, 2, 1)
plt.plot(time, acceleration, lw=1, color='black')
plt.xlabel('Time (s)')
plt.ylabel('Acceleration (g)')
plt.title('Structural Response\nGround Motion')
plt.grid(True)
    
# Fragility curves
plt.subplot(1, 2, 2)
for ds, prob in probabilities.items():
    plt.plot(im_values, prob, lw=2, label=ds)
plt.xlabel('Peak Ground Acceleration (g)')
plt.ylabel('Probability of Exceedance')
plt.title('Fragility Curves')
plt.legend()
plt.grid(True)
    
plt.tight_layout()
plt.show()    
#------------------------------------------------------------------------------------------------
    
