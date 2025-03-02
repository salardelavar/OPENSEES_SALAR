###########################################################################################################
#                                          IN THE NAME OF ALLAH                                           #
#           FRAGILITY ANALYSIS BASED ON ACCELERATION AND STRUCTURAL DUCTILITY DAMAGE INDEX WITH           #
#         RESPONSE SPECTRUM BASED ON DYNAMIC ANALYSIS OF A SINGLE-DEGREE-OF-FREEDOM (SDOF) SYSTEM         #
#---------------------------------------------------------------------------------------------------------#
# This program performs Dynamic Analysis on a Single-Degree-of-Freedom (SDOF) system                      #
# subjected to seismic ground motions. The analysis evaluates the structural response under varying       #
# levels of seismic intensity.                                                                            #
# The framework is designed to support researchers and engineers in assessing the probabilistic seismic   #
# performance of structures, with a focus on understanding the impact of uncertainty on structural        #
# response and design.                                                                                    #
#---------------------------------------------------------------------------------------------------------#
# Key Features:                                                                                           #
# - Simulation of SDOF system using OpenSees.                                                             #
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
import time as TI
import SALAR_MATH as S01
import Analysis_Function as S02
from scipy.stats import norm

#------------------------------------------------------------------------------------------------
# Define parameters (units: m, N)
J_MAX = 400      # Incremental Dynamic Analysis steps for the simulation
#------------------------------------------------------------------------------------------------
# Define  Steel Material Properties (Steel01)
fy = 0.41        # [N] Yield force of structure
fu = 1.5 * fy    # [N] Ultimate force of structure
KE = 2.1e2       # [N/m] Spring elastic stiffness
ey = fy / KE     # [m] Yield displacement
esu = 0.36       # [m] Ultimate displacement
Esh = (fu - fy) / (esu - ey)                    # [N/m] Displacement hardening modulus
b = Esh / KE                                    # Displacement hardening ratio
M = 36           # [kg] Mass of the structure
DR = 0.02         # Damping ratio
KP = fu/esu

T_ELASTIC = 2 * np.pi * np.sqrt(M/KE)  # Period of Elastic Structure
T_PLASTIC = 2 * np.pi * np.sqrt(M/KP)  # Period of Plastic Structure

print('Period of Elastic Structure: ', T_ELASTIC)
print('Period of Plastic Structure: ', T_PLASTIC)

duration = 15.0  # [s] Total simulation duration
dt = 0.01        # [s] Time step
#------------------------------------------------------------------------------------------------
# Define Analysis Properties
MAX_ITERATIONS = 1000000   # Convergence iteration for test
MAX_TOLERANCE = 1.0e-10    # Convergence tolerance for test
#------------------------------------------------------------------------------------------------
### OPENSEES FUNCTION
def ANALYSIS_SPECTRUM_SDOF(j, J_MAX):
    # Initialize OpenSees model
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)
    GMfact = 9.81 # [m/s^2] standard acceleration of gravity or standard acceleration 
    m = ((j+1) / J_MAX) * (T_PLASTIC/2*np.pi)**2 * KP  
    # Define nodes
    ops.node(1, 0.0)  # Fixed base
    ops.node(2, 0.0)  # Mass node
        
    # Define boundary conditions
    ops.fix(1, 1)
    
    # Define mass
    ops.mass(2, m)
        
    # Define material properties
    ops.uniaxialMaterial('Steel01', 1, fy, KE, b) # Steel with bilinear kinematic hardening Material
        
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
    PERIOD = (2*np.pi)/Omega
    
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
    return time, displacement, velocity, acceleration, base_reaction, DI, PERIOD

#------------------------------------------------------------------------------------------------
# Analysis Durations:
starttime = TI.process_time()

# Initialize lists to store max values
max_time = []
max_displacement = []
max_velocity = []
max_acceleration = []
max_base_reaction = []
max_DI = []
max_T = []

# IDA ANALYSIS
for j in range(J_MAX):
    time, displacement, velocity, acceleration, base_reaction, DI, T = ANALYSIS_SPECTRUM_SDOF(j, J_MAX)
    # Calculate and store the max absolute values
    max_time.append(np.max(np.abs(time)))
    max_displacement.append(np.max(np.abs(displacement)))
    max_velocity.append(np.max(np.abs(velocity)))
    max_acceleration.append(np.max(np.abs(acceleration)))
    max_base_reaction.append(np.max(np.abs(base_reaction)))
    #max_DI.append(np.max(np.abs(DI)))
    max_DI.append(DI[-1])
    max_T.append(T)
    print(f'STEP {j + 1} DONE') 
else:
    print('Analysis completed successfully')             

totaltime = TI.process_time() - starttime
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
print("Period of structure:", max_T[-1])
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
    'Ductility_Damage_Index': max_DI,
    'Structure_Peiod': max_T,
}
# Convert to DataFrame
results_df = pd.DataFrame(DATA_TOTAL)
# Export the DataFrame to an Excel file
results_df.to_excel('INELASTIC_SEISMIC_RESPONSE_SDOF_FRAGILITY_RESULTS.xlsx', index=False)
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
S01.PLOT_SCATTER(X, Y , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 7)

# CLUSTER DATA
S01.CLUSTER_DATA(X, Y, XLABEL, YLABEL, MAX_CLUSTERS=3)
#------------------------------------------------------------------------------------------------
XLABEL = 'Acceleration'
YLABEL = 'Base Reaction'
TITLE = f'{YLABEL} and {XLABEL} scatter chart'
COLOR = 'lime'
X = max_acceleration
Y = max_base_reaction
S01.PLOT_SCATTER(X, Y , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 7)

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
XLABEL = 'Displacement'
YLABEL = 'Acceleration'
TITLE = f'{YLABEL} and {XLABEL} scatter chart'
COLOR = 'pink'
X = max_displacement
Y = max_acceleration
S01.PLOT_SCATTER(X, Y , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 7)

# CLUSTER DATA
S01.CLUSTER_DATA(X, Y, XLABEL, YLABEL, MAX_CLUSTERS=3)
#------------------------------------------------------------------------------------------------
XLABEL = 'Period'
YLABEL = 'displacement'
TITLE = f'{YLABEL} and {XLABEL} scatter chart'
COLOR = 'yellow'
X = max_T
Y = max_displacement
S01.PLOT_SCATTER(X, Y , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 7)

# CLUSTER DATA
S01.CLUSTER_DATA(X, Y, XLABEL, YLABEL, MAX_CLUSTERS=3)
#------------------------------------------------------------------------------------------------
XLABEL = 'Period'
YLABEL = 'Velocity'
TITLE = f'{YLABEL} and {XLABEL} scatter chart'
COLOR = 'brown'
X = max_T
Y = max_velocity

S01.PLOT_SCATTER(X, Y , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 7)

# CLUSTER DATA
S01.CLUSTER_DATA(X, Y, XLABEL, YLABEL, MAX_CLUSTERS=3)
#------------------------------------------------------------------------------------------------
XLABEL = 'Period'
YLABEL = 'Acceleration'
TITLE = f'{YLABEL} and {XLABEL} scatter chart'
COLOR = 'pink'
X = max_T
Y = max_acceleration

S01.PLOT_SCATTER(X, Y , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 7)

# CLUSTER DATA
S01.CLUSTER_DATA(X, Y, XLABEL, YLABEL, MAX_CLUSTERS=3)
#------------------------------------------------------------------------------------------------
# PLOT THE TIME-HISTORY
S01.PLOT_TIME_HISTORY(time, displacement, velocity, acceleration, base_reaction)
#------------------------------------------------------------------------------------------------
####  FRAGILITY ANALYSIS BASED ON ACCELERATION AND STRUCTURAL DUCTILITY DAMAGE INDEX
 
####  FRAGILITY ANALYSIS BASED ON ACCELERATION : 
# ----------------------------
# Fragility Assessment
# ----------------------------
# Define damage states per FEMA P-58
# INFO LINK: https://www.fema.gov/sites/default/files/documents/fema_p-58-2-se_volume2_implementation.pdf
damage_states = {
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
# --------------
# Visualization
# --------------
plt.figure(1, figsize=(10, 6))
# Response plot
plt.plot(time, acceleration, lw=1, color='black')
plt.xlabel('Time (s)')
plt.ylabel('Acceleration (g)')
plt.title(f'Last Analysis Structural Response + Ground Motion ::: MAX. ABS. : {np.max(np.abs(acceleration)):.4f}')
plt.grid(True)
plt.show()    

# Fragility curves
plt.figure(2, figsize=(10, 6))
# Calculate and plot fragility curves for each damage state
for damage_state, (median, beta) in damage_states.items():
    # Calculate log-normal probabilities
    ln_im = np.log(im_values)
    ln_median = np.log(median)
    probabilities = norm.cdf((ln_im - ln_median) / beta)
    plt.scatter(im_values, probabilities, marker='o', label=f'{damage_state} (η={median}, β={beta}')
    #plt.plot(im_values, probabilities, lw=2, label=f'{damage_state} (η={median}, β={beta})')
plt.xlabel('Peak Ground Acceleration (g)  [IM]')
plt.ylabel('Probability of Exceedance')
plt.title('Fragility Curves')
plt.legend()
plt.semilogy()
plt.ylim(0, 1.0)
plt.grid(True)
plt.tight_layout()
plt.show()    

#===========================================================
####  FRAGILITY ANALYSIS BASED ON STRUCTURAL DUCTILITY DAMAGE INDEX:

# Define damage state parameters: {Damage State: (median_IM, beta)}
damage_states = {
    'Minor Damage Level': (0.2, 0.4),# Median DI=0.2, β=0.4
    'Moderate Damage Level': (0.4, 0.4),
    'Severe Damage Level': (0.6, 0.5),
    'Failure Level': (1.0, 0.5)
}

# Generate intensity measure (IM) values from 0.0 to 1.0
im_values = max_DI # Structural Ductility Damage Index
# --------------
# Visualization
# --------------
# Create plot
plt.figure(figsize=(10, 6))
# Calculate and plot fragility curves for each damage state
for damage_state, (median, beta) in damage_states.items():
    # Calculate log-normal probabilities
    ln_im = np.log(im_values)
    ln_median = np.log(median)
    probabilities = norm.cdf((ln_im - ln_median) / beta)
    plt.scatter(im_values, probabilities, marker='o', label=f'{damage_state} (η={median}, β={beta}')
    #plt.plot(im_values, probabilities, lw=2, label=f'{damage_state} (η={median}, β={beta})')

# Format plot
plt.xlabel('Structural Ductility Damage Index [IM]', fontsize=12)
plt.ylabel('Probability of Exceedance', fontsize=12)
plt.title('Fragility Curves', fontsize=14)
plt.legend(loc='lower right', fontsize=10)
plt.grid(True)
plt.semilogy()
plt.ylim(0, 1.0)
plt.tight_layout()
plt.show() 

#------------------------------------------------------------------------------------------------   

