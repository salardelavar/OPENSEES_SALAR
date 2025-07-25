###########################################################################################################
#                                          IN THE NAME OF ALLAH                                           #
# DAMPING RATIO OPTIMIZATION BASED ON FREE-VIBRATION ANALYSIS OF A SINGLE-DEGREE-OF-FREEDOM (SDOF) SYSTEM #
#                                       DECISION TREE THEORY METHOD                                       #
#---------------------------------------------------------------------------------------------------------#
#                                                                                                         #
# FACTORS INFLUENCING DAMPING IN STRUCTURES:                                                              #
#                                                                                                         #
# DAMPING PLAYS A CRITICAL ROLE IN CONTROLLING STRUCTURAL RESPONSE DURING EARTHQUAKES.                    #
# SEVERAL KEY FACTORS INFLUENCE THE DAMPING CHARACTERISTICS OF A STRUCTURE:                               #
#                                                                                                         #
# 1. MATERIAL PROPERTIES:                                                                                 #
# - DIFFERENT MATERIALS EXHIBIT VARYING LEVELS OF INHERENT DAMPING.                                       #
# - CONCRETE: HIGHER DAMPING DUE TO MICRO-CRACKING AND ENERGY DISSIPATION.                                #
# - STEEL: LOWER INHERENT DAMPING BUT IMPROVES WITH FRICTION AT CONNECTIONS.                              #
#                                                                                                         #
# 2. STRUCTURAL CONFIGURATION:                                                                            #
# - THE GEOMETRY AND STIFFNESS DISTRIBUTION AFFECT ENERGY DISSIPATION.                                    # 
# - SLENDER STRUCTURES: GENERALLY LOWER DAMPING.                                                          #
# - MASSIVE STRUCTURES: MORE INTERNAL ENERGY LOSS.                                                        #
#                                                                                                         #
# 3. CONNECTIONS AND JOINTS:                                                                              #
# - BOLTED AND RIVETED CONNECTIONS INTRODUCE ADDITIONAL FRICTIONAL DAMPING.                               #
# - WELDED JOINTS MAY REDUCE DAMPING DUE TO RIGID BEHAVIOR.                                               #
#                                                                                                         #
# 4. NON-STRUCTURAL COMPONENTS:                                                                           #
# - PARTITIONS, FACADES, AND INTERNAL FINISHES CAN INCREASE DAMPING.                                      #
# - DAMAGE TO NON-STRUCTURAL ELEMENTS DURING EARTHQUAKES ENHANCES ENERGY DISSIPATION.                     # 
#                                                                                                         #
# 5. LOADING AND DYNAMIC EFFECTS:                                                                         #
# - HIGHER AMPLITUDE VIBRATIONS TEND TO INCREASE EFFECTIVE DAMPING DUE TO MATERIAL                        #
# YIELDING AND FRICTIONAL EFFECTS.                                                                        #
# - REPEATED LOADING CYCLES CAN ALTER DAMPING CHARACTERISTICS OVER TIME.                                  #
#                                                                                                         #
# UNDERSTANDING AND OPTIMIZING DAMPING MECHANISMS IS CRUCIAL FOR SEISMIC DESIGN,                          #
# HELPING ENGINEERS IMPROVE STRUCTURAL RESILIENCE AND EARTHQUAKE PERFORMANCE.                             #
#---------------------------------------------------------------------------------------------------------#
#                          THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                     #
#                                   EMAIL: salar.d.ghashghaei@gmail.com                                   #
########################################################################################################### 

import openseespy.opensees as ops
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time as TI
#import SALAR_MATH as S01
import ANALYSIS_FUNCTION as S02

#------------------------------------------------------------------------------------------------
# Define parameters (units: m, N)
# Define  Steel Material Properties (Steel01)
FY = 41                       # [N] Yield force of structure
FU = 1.5 * FY                 # [N] Ultimate force of structure
KE = 2.1e4                    # [N/m] Spring elastic stiffness
DY = FY / KE                  # [m] Yield displacement
DU = 3.6                      # [m] Ultimate displacement
Esh = (FU - FY) / (DU - DY)   # [N/m] Displacement hardening modulus
b = Esh / KE                  # Displacement hardening ratio
M = 80                        # [kg] Mass of the structure
KP = FU/DU                    # [N/m] Spring plastic stiffness

u0 = 0.005                     # [m] Initial displacement applied to the node
DR = 0.01                      # Intial Guess for Damping ratio
#C = 50.0                      # Intial Guess Damping coefficient (N/(m/s)^alpha)
#alpha = 1.0                   # Velocity exponent (linear damper)

duration = 100.0              # [s] Total simulation duration
dt = 0.01                     # [s] Time step

T_ELASTIC = 2 * np.pi * np.sqrt(M/KE)  # Period of Elastic Structure
T_PLASTIC = 2 * np.pi * np.sqrt(M/KP)  # Period of Plastic Structure

print('Period of Elastic Structure: ', T_ELASTIC)
print('Period of Plastic Structure: ', T_PLASTIC)
#------------------------------------------------------------------------------------------------
# Calculate Over Strength Coefficient (Ω0)
Omega_0 = FU / FY
# Calculate Displacement Ductility Ratio (μ)
mu = DU / DY
# Calculate Ductility Coefficient (Rμ)
#R_mu = (2 * mu - 1) ** 0.5
#R_mu = 1
R_mu = mu
# Calculate Structural Behavior Coefficient (R)
R = Omega_0 * R_mu
print(f'Over Strength Coefficient (Ω0):      {Omega_0:.2f}')
print(f'Displacement Ductility Ratio (μ):    {mu:.2f}')
print(f'Ductility Coefficient (Rμ):          {R_mu:.2f}')
print(f'Structural Behavior Coefficient (R): {R:.2f}')
#------------------------------------------------------------------------------------------------
# Define Analysis Properties
MAX_ITERATIONS = 1000000   # Convergence iteration for test
MAX_TOLERANCE = 1.0e-10    # Convergence tolerance for test
#------------------------------------------------------------------------------------------------
### OPENSEES FUNCTION
def ANALYSIS_DYN_SDOF(Damping_Ratio):
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
    ops.mass(2, M)
            
    # Define material properties
    MatTag = 1
    #ops.uniaxialMaterial('Steel01', MatTag, FY, KE, b)           # Steel with bilinear kinematic hardening Material
    # Define hysteretic material
    pinchX = 0.8           # Pinching factor in X direction
    pinchY = 0.5           # Pinching factor in Y direction
    damage1 = 0.0          # Damage due to ductility
    damage2 = 0.0          # Damage due to energy
    beta = 0.1             # Stiffness degradation parameter
    ops.uniaxialMaterial('HystereticSM', MatTag, FY , DY, FU , DU, 0.2*FU, 1.1*DU, 0.1*FU , 1.2*DU,
                         -FY, -DY, -FU , -DU, -0.2*FU, -1.05*DU, -0.0*FU, -1.1*DU,
                         pinchX, pinchY, damage1, damage2, beta)
    
    #ops.uniaxialMaterial('ViscousDamper', 2, KE, C, alpha)  # Viscous damper
    # INFO LINK: https://opensees.berkeley.edu/wiki/index.php?title=ViscousDamper_Material

    # Define element combining spring and damper
    #ops.element('zeroLength', 1, 1, 2, '-mat', 1, 2, '-dir', 1, 1)  # Both materials act in X-direction
 
    # Define element
    ops.element('zeroLength', 1, 1, 2, '-mat', 1, '-dir', 1)  # DOF[1] LATERAL SPRING
    
    # Static analysis to apply initial displacement
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(2, 1.0)

    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.system('BandGeneral')
    ops.algorithm('Newton')
    ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS)
    ops.integrator('DisplacementControl', 2, 1, u0) # Initial displacement applied to the node 2
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
    a0 = (2 * Omega01 * DR) / Omega01 # c = a0 * m : Mass-proportional damping
    a1 = (DR * 2) / Omega01 # c = a1 * k : Stiffness-proportional damping
    # Apply Rayleigh damping
    ops.rayleigh(a0, a1, 0, 0)# INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/reyleigh.html
    #ops.rayleigh(0, 0, 2 * DR * Omega01, 0) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/reyleigh.html
    PERIOD = np.pi / Omega01   # Structure Period  
    
    # Dynamic analysis
    time = []
    displacement = []
    velocity = []
    acceleration = []
    base_reaction = []
    DI = []
        
    stable = 0
    current_time = 0.0
    while stable == 0 and current_time < duration:
        stable = ops.analyze(1, dt)
        S02.ANALYSIS(stable, 1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS
        current_time = ops.getTime()
        time.append(current_time)
        displacement.append(ops.nodeDisp(2, 1))
        velocity.append(ops.nodeVel(2, 1))
        acceleration.append(ops.nodeAccel(2, 1)) # Structure Acceleration
        base_reaction.append(-ops.eleResponse(1, 'force')[0])  # Reaction force
        DI.append((displacement[-1] - DY) / (DU - DY))        # Structural Ductility Damage Index 
        
    # Calculating Damping Ratio Using Logarithmic Decrement Analysis    
    peaks = np.array([displacement[i] for i in range(1, len(displacement)-1) if displacement[i] > displacement[i-1] and displacement[i] > displacement[i+1]])
    # Natural logarithm
    delta = np.log(peaks[:-1] / peaks[1:])
    #----------------------------------------------------------
    # APPROXIMATE SOLUTION:
    # Compute average logarithmic decrement - approximate equation
    xi_calculated = np.mean(delta) / (2 * np.pi)
    C_xi_calculated = 2 * xi_calculated * Omega01 * M  # [N/(m/s)] Damping coefficient 
    print(f'C: {C_xi_calculated:.8e}')
    ops.wipe()
    return time, displacement, velocity, acceleration, base_reaction, DI, PERIOD, xi_calculated

#------------------------------------------------------------------------------------------------
# DAMPING RATIO OPTIMIZATION DECISION TREE METHOD:

X_bounds = (0.0, 0.01)  # Bounds for Intial Guess Damping Ratio 
NUM_SAMPLES = 20        # Number of samples for data generation
TOLERANCE = 1e-10       # Convergence Tolerance
RESIDUAL = 100          # Convergence Residual 
IT = 0                  # Intial Iteration
ITMAX = 100000          # Max. Iteration


# -----------------------------------------------------------------------------
# FIND THE OPTIMUM VALUE (DECISION TREE REGRESSION FOR OPTIMAL DAMPING RAIO)
# -----------------------------------------------------------------------------

def DECISION_TREE_OPTIMIZER():
    from sklearn.tree import DecisionTreeRegressor
    # Generate sample data
    X_samples = np.linspace(X_bounds[0], X_bounds[1], NUM_SAMPLES)
    XI_samples = []
    it = 0;
    for X in X_samples:
        # Run structural analysis for current diameter
        time, displacement, velocity, acceleration, base_reaction, DI, T, XI = ANALYSIS_DYN_SDOF(X)
        F = XI - X
        print(f'\n Iteration {it+1} - Damping Ratio: {X} - F: {F}\n\n')
        it += 1
        XI_samples.append(X)
    
    X_samples = X_samples.reshape(-1, 1)
    XI_samples = np.array(XI_samples)
    
    # Train Decision Tree Regressor
    regressor = DecisionTreeRegressor()
    regressor.fit(X_samples, XI_samples)
    
    # Fine grid search within bounds
    fine_damping_ratios = np.linspace(X_bounds[0], X_bounds[1], 1000).reshape(-1,1)
    predicted_damping_ratios = regressor.predict(fine_damping_ratios)
    residuals = np.abs(predicted_damping_ratios - XI)
    
    # Find diameter with minimum residual
    min_index = np.argmin(residuals)
    best_damping_ratio = fine_damping_ratios[min_index][0]
    final_residual = residuals[min_index]
    
    return time, displacement, velocity, acceleration, base_reaction, DI, T, best_damping_ratio, final_residual

# Analysis Durations:
starttime = TI.process_time()

# Run Decision Tree optimization
time, displacement, velocity, acceleration, base_reaction, DI, T, best_damping_ratio, final_residual = DECISION_TREE_OPTIMIZER()

# Results output
print(f'\n\t\t Optimum Damping Ratio: {best_damping_ratio:.6e}')
print(f'\t\t Final Residual:         {final_residual:.6e}')
if final_residual <= TOLERANCE:
    print("\t\t CONVERGENCE ACHIEVED")
else:
    print("\t\t OPTIMAL SOLUTION FOUND WITHIN TOLERANCE LIMITS")

totaltime = TI.process_time() - starttime
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
plt.ylabel('Displacement [m]')
plt.title(f'Time vs Displacement - MAX. ABS: {DISP_Z[-1]} | ξ (Calculated): {100*best_damping_ratio:.5e} %')
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

plt.figure(5, figsize=(8, 6))
plt.plot(displacement, base_reaction, color='black', linewidth=2)
plt.xlabel('Displacement [m]')
plt.ylabel('Base-reaction [N]')
plt.title(f'Displacement vs Base-reaction Digram')
plt.grid()
plt.show()

plt.figure(6, figsize=(8, 6))
plt.plot(displacement, DI, color='brown', linewidth=2)
plt.xlabel('Displacement [m]')
plt.ylabel('Ductility Damage Index')
plt.title(f'Displacement vs Ductility Damage Index Digram')
plt.grid()
plt.show()
#------------------------------------------------------------------------------------------------  
# EXPORT DATA TO EXCEL
DATA_TOTAL = {
    'Max_displacement': displacement,
    'Max_velocity': velocity,
    'Max_acceleration': acceleration,
    'Max_Base_Reaction': base_reaction,
    'Ductility_Damage_Index': DI,
    'Structure_Peiod': T,
    'Damping_Ratio': best_damping_ratio,
}
# Convert to DataFrame
results_df = pd.DataFrame(DATA_TOTAL)
# Export the DataFrame to an Excel file
results_df.to_excel('INELASTIC_FREE_VIBRATION_DAMPING_RATIO_SDOF_OPTIMIZATION_RESULTS.xlsx', index=False)
#------------------------------------------------------------------------------------------------

