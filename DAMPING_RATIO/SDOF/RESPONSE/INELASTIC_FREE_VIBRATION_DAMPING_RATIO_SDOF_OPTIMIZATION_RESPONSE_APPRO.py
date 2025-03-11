###########################################################################################################
#                                          IN THE NAME OF ALLAH                                           #
# DAMPING RATIO OPTIMIZATION BASED ON FREE-VIBRATION ANALYSIS OF A SINGLE-DEGREE-OF-FREEDOM (SDOF) SYSTEM #
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
import SALAR_MATH as S01
import ANALYSIS_FUNCTION as S02
from scipy.stats import norm

#------------------------------------------------------------------------------------------------
# Define parameters (units: m, N)
# Define  Steel Material Properties (Steel01)
FY = 41                       # [N] Yield force of structure
FU = 1.5 * FY                 # [N] Ultimate force of structure
KE = 2.1e4                    # [N/m] Spring elastic stiffness
EY = FY / KE                  # [m] Yield displacement
EU = 3.6                      # [m] Ultimate displacement
Esh = (FU - FY) / (EU - EY)   # [N/m] Displacement hardening modulus
b = Esh / KE                  # Displacement hardening ratio
MASS = 80                     # [kg] Mass of the structure
KP = FU/EU                    # [N/m] Spring plastic stiffness

u0 = 0.5                      # [m] Initial displacement applied to the node
DR = 0.01                     # Intial Guess for Damping ratio
C = 50.0                      # Intial Guess Damping coefficient (N/(m/s)^alpha)
alpha = 1.0                   # Velocity exponent (linear damper)

duration = 100.0              # [s] Total simulation duration
dt = 0.01                     # [s] Time step

T_ELASTIC = 2 * np.pi * np.sqrt(MASS/KE)  # Period of Elastic Structure
T_PLASTIC = 2 * np.pi * np.sqrt(MASS/KP)  # Period of Plastic Structure

print('Period of Elastic Structure: ', T_ELASTIC)
print('Period of Plastic Structure: ', T_PLASTIC)

#------------------------------------------------------------------------------------------------
# Define Analysis Properties
MAX_ITERATIONS = 1000000   # Convergence iteration for test
MAX_TOLERANCE = 1.0e-10    # Convergence tolerance for test
#------------------------------------------------------------------------------------------------
### OPENSEES FUNCTION
def ANALYSIS_DYN_SDOF(Damping_Ratio, u0, M):
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
    ops.uniaxialMaterial('Steel01', 1, FY, KE, b)           # Steel with bilinear kinematic hardening Material
    # INFO LINK: https://opensees.berkeley.edu/wiki/index.php?title=ViscousDamper_Material
    #ops.uniaxialMaterial('ViscousDamper', 2, KE, C, alpha)  # Viscous damper

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
    Omega01 = np.power(max(Lambda01), 0.5)
    Omega02 = 2 * Omega01
    a0 = Damping_Ratio * (2 * Omega01 * Omega02) / (Omega01 + Omega02) # c = a0 * m : Mass-proportional damping
    a1 = Damping_Ratio * 2 / (Omega01 + Omega02) # c = a1 * k : Stiffness-proportional damping
    # Apply Rayleigh damping
    #ops.rayleigh(a0, a1, 0, 0)# INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/reyleigh.html
    ops.rayleigh(0, 0, 2 * Damping_Ratio * Omega01, 0) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/reyleigh.html
    PERIOD = (2 * np.pi) / Omega01 
    
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
        DI.append((displacement[-1] - EY) / (EU - EY))        # Structural Ductility Damage Index 
        
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
    return time, displacement, velocity, acceleration, base_reaction, DI, PERIOD, xi_calculated, delta

#------------------------------------------------------------------------------------------------
# DAMPING RATIO OPTIMIZATION:
    
X = DR             # Intial Guess Damping Ratio 
ESP = 1e-5         # Finite difference derivative Convergence Tolerance
TOLERANCE = 1e-10  # Convergence Tolerance
ITMAX = 100000     # Max. Iteration
ITERATION = 1      # Steps Count
NI = 20            # Steps for Intial Displacement 
NJ = 10            # Steps for Structure Mass

UI, MI = [], []
DISP, VELO, ACCEL, REACTION, DII, TII, XXI = [], [], [], [], [], [], []

#starttime = TI.process_time()
current_time = TI.strftime("%H:%M:%S", TI.localtime())
print(f"Current time (HH:MM:SS): {current_time}\n\n")
    
for I in range(NI):
    U = EU * (I+1) / NI  # Intial displacement of each increment                     
    for J in range(NJ):
        print('\n')
        print('------------------------')
        print('  ITERATION       ', ITERATION)
        print('------------------------')
        M =  MASS * (J+1) / NJ  # Mass of each increment                    
        ### FIND THE OPTIMUM VALUE 
        RESIDUAL = 100     # Convergence Residual 
        IT = 0             # Intial Iteration
        X = DR             # Intial Guess Damping Ratio
        while (RESIDUAL > TOLERANCE):
            # X -------------------------------------------------------
            time, displacement, velocity, acceleration, base_reaction, DI, T, XI = ANALYSIS_DYN_SDOF(X, U, M)
            print(f'XI: {XI:.8e}')
            F = XI - X
            print('F: ', F)
            # Xmin -------------------------------------------------------
            XMIN = X - ESP  
            time, displacement, velocity, acceleration, base_reaction, DI, T, XImin = ANALYSIS_DYN_SDOF(XMIN, U, M)
            Fmin = XImin - XMIN
            print('Fmin: ', Fmin)
            # Xmax -------------------------------------------------------
            XMAX = X + ESP  
            time, displacement, velocity, acceleration, base_reaction, DI, T, XImax = ANALYSIS_DYN_SDOF(XMAX, U, M)
            Fmax = XImax - XMAX
            print('Fmax: ', Fmax)
            # DF -------------------------------------------------------
            DF = (Fmax - Fmin)/(2 * ESP);# Calculate the Finite difference derivative of F
            print('DF: ', DF)
            # DX -------------------------------------------------------
            DX = F / DF;        # Calculate dx
            print('DX: ', DX)
            # RESIDUAL -------------------------------------------------
            RESIDUAL = abs(DX); # Calculate residual
            
            print('IT: ', IT+1,' - RESIDUAL: ', RESIDUAL,' - DAMPING RATIO: ', X,'\n')
            X -= DX;            # update X
            IT += 1;            # update iteration
            if IT == ITMAX:
                print('\t\t Iteration reached to Max. Iteration')
                print('\t\t Change ESP and TOLERANCE for better Convergence')
                X = -1
                break;
            if RESIDUAL < TOLERANCE:
                print(f'\t\t Optimum Damping Ratio:      {X:.6e}')
                print(f'\t\t Iteration Counts:           {IT}')
                print(f'\t\t Convergence Residual:       {RESIDUAL:.10e}')
                
        ITERATION += 1     
        if X >= 0.0: # WE ARE GOING TO OUTPUT JUST POSITIVE DAMPING RATIOS DATA 
            MI.append(M)
            UI.append(U)
            XXI.append(X)
            TII.append(T)
            DISP.append(np.max(np.abs(displacement)))
            VELO.append(np.max(np.abs(velocity)))
            ACCEL.append(np.max(np.abs(acceleration)))
            REACTION.append(np.max(np.abs(base_reaction)))
            DII.append(np.max(np.abs(DI)))
            

#totaltime = TI.process_time() - starttime
#print(f'\nTotal time (s): {totaltime:.4f}')
current_time = TI.strftime("%H:%M:%S", TI.localtime())
print(f"Current time (HH:MM:SS): {current_time}\n\n")
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
plt.title(f'Time vs Displacement - MAX. ABS: {DISP_Z[-1]} | ξ (Calculated): {100*XI:.5e} %')
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
    'Initial_displacement': UI,
    'Structure_Mass': MI,
    'Max_displacement': DISP,
    'Max_velocity': VELO,
    'Max_acceleration': ACCEL,
    'Max_Base_Reaction': REACTION,
    'Ductility_Damage_Index': DII,
    'Structure_Period': TII,
    'Damping_Ratio': XXI,
}
# Convert to DataFrame
results_df = pd.DataFrame(DATA_TOTAL)
# Export the DataFrame to an Excel file
results_df.to_excel('INELASTIC_FREE_VIBRATION_DAMPING_RATIO_SDOF_OPTIMIZATION_RESPONSE_RESULTS.xlsx', index=False)
#------------------------------------------------------------------------------------------------
# PLOT HEATMAP FOR CORRELATION 
S01.PLOT_HEATMAP(results_df)
#------------------------------------------------------------------------------------------------  
XLABEL = 'Displacement'
YLABEL = 'Base Reaction'
TITLE = f'{YLABEL} and {XLABEL} scatter chart'
COLOR = 'orange'
X = DISP
Y = REACTION
S01.PLOT_SCATTER(X, Y , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

# CLUSTER DATA
S01.CLUSTER_DATA(X, Y, XLABEL, YLABEL, MAX_CLUSTERS=3)
#------------------------------------------------------------------------------------------------
XLABEL = 'Velocity'
YLABEL = 'Base Reaction'
TITLE = f'{YLABEL} and {XLABEL} scatter chart'
COLOR = 'cyan'
X = VELO
Y = REACTION
S01.PLOT_SCATTER(X, Y , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

# CLUSTER DATA
S01.CLUSTER_DATA(X, Y, XLABEL, YLABEL, MAX_CLUSTERS=3)
#------------------------------------------------------------------------------------------------
XLABEL = 'Acceleration'
YLABEL = 'Base Reaction'
TITLE = f'{YLABEL} and {XLABEL} scatter chart'
COLOR = 'lime'
X = ACCEL
Y = REACTION
S01.PLOT_SCATTER(X, Y , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

# CLUSTER DATA
S01.CLUSTER_DATA(X, Y, XLABEL, YLABEL, MAX_CLUSTERS=3)
#------------------------------------------------------------------------------------------------
XLABEL = 'Displacement'
YLABEL = 'Structural Ductility Damage Index'
TITLE = f'{YLABEL} and {XLABEL} scatter chart'
COLOR = 'purple'
X = DISP
Y = DII
S01.PLOT_SCATTER(X, Y , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

# CLUSTER DATA
S01.CLUSTER_DATA(X, Y, XLABEL, YLABEL, MAX_CLUSTERS=3)
#------------------------------------------------------------------------------------------------
XLABEL = 'Period'
YLABEL = 'Initial Displacement'
TITLE = f'{YLABEL} and {XLABEL} scatter chart'
COLOR = 'pink'
X = TII
Y = UI
S01.PLOT_SCATTER(X, Y , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

# CLUSTER DATA
S01.CLUSTER_DATA(X, Y, XLABEL, YLABEL, MAX_CLUSTERS=3)
#------------------------------------------------------------------------------------------------
XLABEL = 'Period'
YLABEL = 'Damping Ratio'
TITLE = f'{YLABEL} and {XLABEL} scatter chart'
COLOR = 'brown'
X = TII
Y = XXI
S01.PLOT_SCATTER(X, Y , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

# CLUSTER DATA
S01.CLUSTER_DATA(X, Y, XLABEL, YLABEL, MAX_CLUSTERS=3)
#------------------------------------------------------------------------------------------------
XLABEL = 'Structural Ductility Damage Index'
YLABEL = 'Damping Ratio'
TITLE = f'{YLABEL} and {XLABEL} scatter chart'
COLOR = 'green'
X = DII
Y = XXI
S01.PLOT_SCATTER(X, Y , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

# CLUSTER DATA
S01.CLUSTER_DATA(X, Y, XLABEL, YLABEL, MAX_CLUSTERS=3)
#------------------------------------------------------------------------------------------------
XLABEL = 'Structural Ductility Damage Index'
YLABEL = 'Period'
TITLE = f'{YLABEL} and {XLABEL} scatter chart'
COLOR = 'blue'
X = DII
Y = TII
S01.PLOT_SCATTER(X, Y , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

# CLUSTER DATA
S01.CLUSTER_DATA(X, Y, XLABEL, YLABEL, MAX_CLUSTERS=3)
#------------------------------------------------------------------------------------------------
XLABEL = 'Velocity'
YLABEL = 'Damping Ratio'
TITLE = f'{YLABEL} and {XLABEL} scatter chart'
COLOR = 'red'
X = VELO
Y = XXI
S01.PLOT_SCATTER(X, Y , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

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
im_values = ACCEL
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
#plt.semilogy()
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
im_values = DII # Structural Ductility Damage Index
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
#plt.semilogy()
plt.ylim(0, 1.0)
plt.tight_layout()
plt.show() 

#------------------------------------------------------------------------------------------------   
