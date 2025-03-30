###########################################################################################################
#                                         IN THE NAME OF ALLAH                                            #
#                    DYNAMIC RESPONSE ANALYSIS OF A SINGLE-DEGREE-OF-FREEDOM (SDOF)                       #
#         INELASTIC STRUCTURE UNDER SEISMIC LOADING WITH UNCERTAINTY USING MONTE CARLO SIMULATION:        #
#                 INCORPORATING BETA PROBABILITY DISTRIBUTION FOR STOCHASTIC PARAMETERS                   #
#---------------------------------------------------------------------------------------------------------#
# This program models and analyzes the dynamic response of a single-degree-of-freedom (SDOF) structural   #
# system subjected to seismic accelerations while incorporating uncertainties in structural properties.   #
# The framework supports researchers and engineers in assessing the probabilistic performance of          #
# structures under seismic excitation, emphasizing the role of uncertainty in seismic response and        #
# design.                                                                                                 #
#---------------------------------------------------------------------------------------------------------#
#                          THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                     #
#                                   EMAIL: salar.d.ghashghaei@gmail.com                                   #
###########################################################################################################
#
# Objectives:
# 1. Stochastic Parameter Modeling:
#    - Use Beta probability distribution functions to generate random values for:
#        Stiffness (k)
#        Mass (m)
#        Damping ratio (ζ)
#        Seismic accelerations
#    - Provide a statistical representation of uncertainties in structural properties and loading conditions.
#
# 2. Structural Model Development:
#    - Build an SDOF system in OpenSeesPy with stochastic parameters:
#        * Variable stiffness, mass, and damping ratio
#    - Apply seismic accelerations as dynamic input to the system.
#
# 3. Monte Carlo Simulation:
#    - Perform transient dynamic analysis iteratively over numerous realizations of stochastic parameters.
#    - Capture variability in system behavior and response.
#
# 4. Response Analysis:
#    - Monitor and record critical system responses:
#        Displacement, velocity, and acceleration
#        Base reaction forces over time
#    - Aggregate simulation results to compute probabilistic metrics of structural response.
#
# 5. Visualization and Statistical Assessment:
#    - Generate graphical representations for enhanced insights:
#        Histograms and boxplots for distributions of maximum displacement, velocity, acceleration, and base reaction forces.
#        Time-history plots illustrating system response dynamics for representative cases.
#
# 6. Seismic Safety Assessment:
#    - Evaluate seismic safety using simulation data by training a Random Forest Classifier to
#        predict system safety under varying conditions.
#
# 7. Correlation Analysis:
#    - Create a heatmap to visualize correlations between key parameters and responses.
#
# 8. Multiple Regression Modeling:
#    - Develop a regression model to estimate system responses based on stochastic parameters.
#
# 9. Machine Learning: Long Short-Term Memory (LSTM):
#    - Implement LSTM networks to predict dynamic responses over time for advanced probabilistic modeling.
#
# 10. Reliability Analysis:
#    - Perform reliability assessments of base reactions and element capacities to quantify structural safety margins.
#
# 11. Markov Chain Model:
#    - Structural damage analysis by evaluating displacement
#
#---------------------------------------------------------------------------------------------------------
# This framework integrates stochastic modeling, dynamic simulation, and machine learning to provide a   
# robust tool for evaluating seismic performance and ensuring structural reliability under uncertainty.   


#------------------------------------------------------------------------------------------------
import openseespy.opensees as ops
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time as ti
import SALAR_MATH as S01
import Analysis_Function as S02
import MARKOV_CHAIN as S03
from scipy.stats import norm

#------------------------------------------------------------------------------------------------
# Define parameters (units: m, N)
NUM_SIM = 6000                                   # Total number for simulation
#------------------------------------------------------------------------------------------------
# Define  Steel Material Properties (Steel01)
fy = S01.BETA_PDF(0.39, 0.41, 1, 2, NUM_SIM)     # [N] Yield Force of Structure
fu = 1.5 * fy                                    # [N] Ultimate Force of Structure
Es = S01.BETA_PDF(2.0e2, 2.1e2, 1, 2, NUM_SIM)   # [N/m] Spring Stiffness
ey = fy / Es                                     # [m] Yield Displacement
esu = S01.BETA_PDF(0.33, 0.36, 1, 2, NUM_SIM)    # [m] Ultimate Displacement
Esh = (fu - fy) / (esu - ey)                     # [N/m] Displacement Hardening Modulus
b = Esh / Es                                     # Displacement Hardening Ratio
M = S01.BETA_PDF(1500.0, 1700.0, 2, 1, NUM_SIM)  # [kg] Mass of the Structure
DR = S01.BETA_PDF(0.01, 0.03, 1, 1, NUM_SIM)     # Damping Ratio

duration = 15.0  # [s] Total simulation duration
dt = 0.01        # [s] Time step
#------------------------------------------------------------------------------------------------
# Define Analysis Properties
MAX_ITERATIONS = 1000000   # Convergence iteration for test
MAX_TOLERANCE = 1.0e-10    # Convergence tolerance for test
#------------------------------------------------------------------------------------------------
# SPRING STIFFNESS PROPERTIES
S01.HISROGRAM_BOXPLOT(Es, HISTO_COLOR='pink', LABEL='Spring Stiffness [N/m]')
S01.HISROGRAM_BOXPLOT(fy, HISTO_COLOR='blue', LABEL='Yield Force of Structure [N]')
S01.HISROGRAM_BOXPLOT(fu, HISTO_COLOR='gold', LABEL='Ultimae Force of Structure [N]')
S01.HISROGRAM_BOXPLOT(ey, HISTO_COLOR='purple', LABEL='Yield Displacement [m]')
S01.HISROGRAM_BOXPLOT(esu, HISTO_COLOR='green', LABEL='Ultimate Displacement [m]')
S01.HISROGRAM_BOXPLOT(b, HISTO_COLOR='blue', LABEL='Displacement Hardening Ratio')
# SPRING MASS AND DAMPING RATIO PROPERTIES
S01.HISROGRAM_BOXPLOT(M, HISTO_COLOR='yellow', LABEL='Mass')
S01.HISROGRAM_BOXPLOT(DR, HISTO_COLOR='cyan', LABEL='Damping Ratio')
#------------------------------------------------------------------------------------------------
### OPENSEES FUNCTION
def ANALYSIS_SDOF(i):
    # Initialize OpenSees model
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)
    GMfact = 9.81 # [m/s^2] standard acceleration of gravity or standard acceleration 
        
    # Define nodes
    ops.node(1, 0.0)  # Fixed base
    ops.node(2, 0.0)  # Mass node
        
    # Define boundary conditions
    ops.fix(1, 1)
    
    # Define mass
    ops.mass(2, M[i])
        
    # Define material properties
    """
    MatTag = 1
    ops.uniaxialMaterial('Steel01', MatTag, fy[i], Es[i], b[i]) # Steel with bilinear kinematic hardening Material
    """ 
    """
    MatTag = 1
    KP = np.array([[fy[i] , ey[i]], [fu[i], esu[i]], [0.2*fu[i], 1.1*esu[i]], [0.1*fu[i] , 1.2*esu[i]]])           # TENSION STRESS-STRAIN RELATION
    KN = np.array([[-fy[i] , -ey[i]], [-fu[i], -esu[i]], [-0.2*fu[i], -1.05*esu[i]], [-0.0*fu[i] , -1.1*esu[i]]])  # COMPRESSION STRESS-STRAIN RELATION
    ops.uniaxialMaterial('HystereticSM', MatTag, '-posEnv', *KP.flatten(), '-negEnv', *KN.flatten(), '-pinch', 1, 1, '-damage', 0.1, 0.01, '-beta', 0,'-defoLimitStates', ey[i], -ey[i], esu[i], -esu[i], '-forceLimitStates', fy[i], -fy[i], fu[i], -fu[i])# ,'printInput'
    #INFO LINK: https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/HystereticSM.html
    """
    MatTag = 1
    pinchX = 0.8           # Pinching factor in X direction
    pinchY = 0.5           # Pinching factor in Y direction
    damage1 = 0.0          # Damage due to ductility
    damage2 = 0.0          # Damage due to energy
    beta = 0.1             # Stiffness degradation parameter
    ops.uniaxialMaterial('Hysteretic', MatTag, fy[i], ey[i], fu[i], esu[i], 0.2*fu[i], 1.1*esu[i], -fy[i], -ey[i], -fu[i], -esu[i], -0.2*fu[i], -1.1*esu[i], pinchX, pinchY, damage1, damage2, beta)
    # INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Hysteretic_Material
    
    # Define element
    ops.element('zeroLength', 1, 1, 2, '-mat', MatTag, '-dir', 1)  # DOF[1] LATERAL SPRING
    
    # Apply seismic accelerations    
    # Define time series for input motion (Acceleration time history)
    gm_accels = np.loadtxt(f'Ground_Acceleration_{i+1}.txt')  # Assumes acceleration in m/s²
    ops.timeSeries('Path', 1, '-dt', dt, '-values', *gm_accels.tolist(), '-factor', GMfact) # SEISMIC-X
    #ops.timeSeries('Path', 1, '-dt', dt, '-filePath', f'Ground_Acceleration{i+1}.txt', '-factor', GMfact) # SEISMIC-X
        
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
    #Lambda01 = ops.eigen('-genBandArpack', 1) # eigenvalue mode 1
    Omega01 = np.power(max(Lambda01), 0.5)
    a0 = (2 * Omega01 * DR[i]) / Omega01 # c = a0 * m : Mass-proportional damping
    a1 = (DR[i] * 2) / Omega01 # c = a1 * k : Stiffness-proportional damping
    # Apply Rayleigh damping
    ops.rayleigh(a0, a1, 0, 0)# INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/reyleigh.html
    #ops.rayleigh(0, 0, 2 * DR * Omega01, 0) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/reyleigh.html
    PERIOD = np.pi / Omega01   # Structure Period  
    
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
        DI.append((displacement[-1] - ey[i]) / (esu[i] - ey[i]))        # Structural Ductility Damage Index 
        step += 1
        
    ops.wipe()
    return time, displacement, velocity, acceleration, base_reaction, DI, PERIOD

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
max_T = []
# NUM_SIM is the number of simulations
for i in range(NUM_SIM):
    time, displacement, velocity, acceleration, base_reaction, DI, PERIOD = ANALYSIS_SDOF(i)
    # Calculate and store the max absolute values
    max_time.append(np.max(np.abs(time)))
    max_displacement.append(np.max(np.abs(displacement)))
    max_velocity.append(np.max(np.abs(velocity)))
    max_acceleration.append(np.max(np.abs(acceleration)))
    max_base_reaction.append(np.max(np.abs(base_reaction)))
    #max_DI.append(np.max(np.abs(DI)))
    max_DI.append(DI[-1])
    max_T.append(PERIOD)
    print(f'STEP {i + 1} DONE') 
else:
    print('Analysis completed successfully')    

totaltime = ti.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n') 
#------------------------------------------------------------------------------------------------
# Print the last results
print("Maximum Absolute Values Across Simulations:")
print("Period:", max_time[-1])
print("Displacement:", max_displacement[-1])
print("Velocity:", max_velocity[-1])
print("Acceleration:", max_acceleration[-1])
print("Base Reaction:", max_base_reaction[-1])
print("Ductility Damage Index:", max_DI[-1])
print("Period :", max_T)
#------------------------------------------------------------------------------------------------
S01.HISROGRAM_BOXPLOT(max_displacement, HISTO_COLOR='blue', LABEL='Displacement')
S01.HISROGRAM_BOXPLOT(max_velocity, HISTO_COLOR='purple', LABEL='Velocity')
S01.HISROGRAM_BOXPLOT(max_acceleration, HISTO_COLOR='green', LABEL='Acceleration')
S01.HISROGRAM_BOXPLOT(max_base_reaction, HISTO_COLOR='gold', LABEL='Base Reaction')
S01.HISROGRAM_BOXPLOT(max_DI, HISTO_COLOR='pink', LABEL='Ductility Damage Index')
S01.HISROGRAM_BOXPLOT(max_T, HISTO_COLOR='lime', LABEL='Sructure Period')
#------------------------------------------------------------------------------------------------
####  FRAGILITY ANALYSIS
  
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
    'Spring_Stiffness': Es,
    'Yield_Force_of_Structure': fy,
    'Ultimate_Force_of_Structure': fu,
    'Yield_Displacement': ey,
    'Ultimate_Displacement': esu,
    'Displacement_Hardening_Ratio': b,
    'Mass': M,
    'Damping_Ratio': DR,
    'Max_displacement': max_displacement,
    'Max_velocity': max_velocity,
    'Max_acceleration': max_acceleration,
    'Max_Base_Reaction': max_base_reaction,
    'Ductility_Damage_Index': max_DI,
    'Period': max_T,
}
# Convert to DataFrame
results_df = pd.DataFrame(DATA_TOTAL)
# Export the DataFrame to an Excel file
results_df.to_excel('INELASTIC_UNCERTAINTY_SEISMIC_SDOF_RESULTS.xlsx', index=False)
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

# RANDOM FOREST ANALYSIS
"""
This code predicts the seismic safety of a structure using simulation data by training a Random Forest Classifier to
 classify whether the system is "safe" or "unsafe" based on features like maximum displacement, velocity, acceleration,
 and base reaction. A regression model is also trained to estimate safety likelihood. It evaluates model performance using
 metrics like classification accuracy, mean squared error, and R² score. Additionally, it identifies key features influencing
 safety through feature importance analysis. The tool aids in seismic risk assessment, structural optimization, and understanding
 critical safety parameters.
"""

data = {
    'Max_displacement': max_displacement,
    'Max_velocity': max_velocity,
    'Max_acceleration': max_acceleration,
    'Max_Base_Reaction': max_base_reaction,
    'Damping_Ratio': DR,
    'Ductility_Damage_Index': max_DI,
    'Period': max_T,
}


# Convert to DataFrame
df = pd.DataFrame(data)
#print(df)
S01.RANDOM_FOREST(df)
#------------------------------------------------------------------------------------------------
# PLOT HEATMAP FOR CORRELATION 
S01.PLOT_HEATMAP(df)
#------------------------------------------------------------------------------------------------
# MULTIPLE REGRESSION MODEL
S01.MULTIPLE_REGRESSION(df) 
#------------------------------------------------------------------------------------------------
# MACHINE LEARNING: LONG SHORT-TREM MEMERY (LSTM) METHOD
x = max_displacement 
y = max_acceleration 
Demand_X = x[-1]
look_back = 500#int(NUM_SIM * 0.5)
ITERATION = 200
XLABEL = 'Max Displacement'
YLABEL = 'Max Acceleration'
#S01.PREDICT_LSTM(x, y, Demand_X, look_back, ITERATION, XLABEL, YLABEL)
#------------------------------------------------------------------------------------------------
# PERFORM RELIABILITY ANALYSIS FOR BASE REACTION AND ELEMENT CAPACITY
mean_capacity = np.mean(fu)    # Mean Element Ultimate Capacity
std_dev_capacity = np.std(fu)  # Std Element Ultimate Capacity
num_sim = NUM_SIM
S01.RELIABILITY_ANALYSIS(max_base_reaction, num_sim, mean_capacity, std_dev_capacity)
#------------------------------------------------------------------------------------------------
# NEURAL NETWORK FOR FAILURE PROBABILIYY ESTIMATION
X1 = mean_capacity
X2 = max_base_reaction
S01.NEURAL_NETWORK_FAILURE_PROBABILIYY_ESTIMATION(max_base_reaction, X2, NUM_SIM)
#------------------------------------------------------------------------------------------------
# MARKOV CHAIN MODEl (structural damage analysis by evaluating Structural Ductility Damage Index)
FILE_TF = False         # Indicate whether to read data from a file or use provided data
file_path = None        # Not used when 'file_tf' is False
DATA = max_DI # If not using a file, replace None with a NumPy array of data

S03.MARKOV_CHAIN(FILE_TF, file_path, DATA)
#------------------------------------------------------------------------------------------------