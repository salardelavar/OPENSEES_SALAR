###########################################################################################################
#                                          IN THE NAME OF ALLAH                                           #
#  DAMPING RATIO OTIMIZATION BASED ON FREE-IBRATION ANALYSIS OF A SINGLE-DEGREE-OF-FREEDOM (SDOF) SYSTEM  #
#                                            GRAPH THEORY METHOD                                          #
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
import networkx as nx

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
M = 8000                      # [kg] Mass of the structure
KP = FU/EU                    # [N/m] Spring plastic stiffness
u0 = 0.5                      # [m] Initial displacement applied to the node

DR = 0.01                     # Intial Guess for Damping ratio
C = 50.0                      # Intial Guess Damping coefficient (N/(m/s)^alpha)
alpha = 1.0                   # Velocity exponent (linear damper)

duration = 100.0              # [s] Total simulation duration
dt = 0.01                     # [s] Time step

T_ELASTIC = 2 * np.pi * np.sqrt(M/KE)  # Period of Elastic Structure
T_PLASTIC = 2 * np.pi * np.sqrt(M/KP)  # Period of Plastic Structure

print('Period of Elastic Structure: ', T_ELASTIC)
print('Period of Plastic Structure: ', T_PLASTIC)

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
            
    C = 100
    # Define material properties
    ops.uniaxialMaterial('Steel01', 1, FY, KE, b)           # Steel with bilinear kinematic hardening Material
    # INFO LINK: https://opensees.berkeley.edu/wiki/index.php?title=ViscousDamper_Material
    ops.uniaxialMaterial('ViscousDamper', 2, KE, C, alpha)  # Viscous damper

    # Define element combining spring and damper
    ops.element('zeroLength', 1, 1, 2, '-mat', 1, 2, '-dir', 1, 1)  # Both materials act in X-direction
 
    # Define element
    #ops.element('zeroLength', 1, 1, 2, '-mat', 2, '-dir', 1)  # DOF[1] LATERAL SPRING
    
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
    
    C = 2 * Damping_Ratio * Omega01 * M  # [N/(m/s)] Damping coefficient 

    
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
    ops.wipe()
    return time, displacement, velocity, acceleration, base_reaction, DI, PERIOD, xi_calculated

#------------------------------------------------------------------------------------------------
# DAMPING RATIO OPTIMIZATION WITH GRAPH THEORY:
    
X = DR             # Intial Guess Damping Ratio 
ESP = 1e-5         # Finite difference derivative Convergence Tolerance
TOLERANCE = 1e-10  # Convergence Tolerance
RESIDUAL = 100     # Convergence Residual 
IT = 0             # Intial Iteration
ITMAX = 100000     # Max. Iteration


# Analysis Durations:
starttime = TI.process_time()

### FIND THE OPTIMUM VALUE 
# Create a directed graph to record the iterations
G = nx.DiGraph()
G.add_node(IT, X_value=X, residual=RESIDUAL)

# Execute the Newton-Raphson algorithm with graph visualization
while (RESIDUAL > TOLERANCE):
    # X -------------------------------------------------------
    time, displacement, velocity, acceleration, base_reaction, DI, T, XI = ANALYSIS_DYN_SDOF(X)
    print(f'XI: {XI:.8e}')
    F = XI - X
    print('F: ', F)
    # Xmin -------------------------------------------------------
    XMIN = X - ESP  
    time, displacement, velocity, acceleration, base_reaction, DI, T, XImin = ANALYSIS_DYN_SDOF(XMIN)
    Fmin = XImin - XMIN
    print('Fmin: ', Fmin)
    # Xmax -------------------------------------------------------
    XMAX = X + ESP  
    time, displacement, velocity, acceleration, base_reaction, DI, T, XImax = ANALYSIS_DYN_SDOF(XMAX)
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
    
    G.add_node(IT, X_value=X, residual=RESIDUAL)
    G.add_edge(IT - 1, IT, weight=DX)
    
    if IT == ITMAX:
        print('\t\t Iteration reached to Max. Iteration')
        print('\t\t Change ESP and TOLERANCE for better Convergence')
        X = -1
        break;
    if RESIDUAL < TOLERANCE:
        print(f'\t\t Optimum Damping Ratio:      {X:.6e}')
        print(f'\t\t Iteration Counts:           {IT}')
        print(f'\t\t Convergence Residual:       {RESIDUAL:.10e}')
        

# Visualize the iterative process as a directed graph
pos = nx.spring_layout(G)
node_labels = {node: f"X={data['X_value']:.4f}" for node, data in G.nodes(data=True)}
edge_labels = nx.get_edge_attributes(G, 'weight')

plt.figure(figsize=(10, 6))
nx.draw_networkx_nodes(G, pos, node_color='lightblue', node_size=500)
nx.draw_networkx_edges(G, pos, arrowstyle='->', arrowsize=10, edge_color='gray')
nx.draw_networkx_labels(G, pos, labels=node_labels, font_size=10)
nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_color='red', font_size=8)
plt.title("Graph-Theoretic Representation of Damping Ratio Optimization")
plt.axis('off')
plt.show()
           

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
    'Max_displacement': displacement,
    'Max_velocity': velocity,
    'Max_acceleration': acceleration,
    'Max_Base_Reaction': base_reaction,
    'Ductility_Damage_Index': DI,
    'Structure_Peiod': T,
    'Damping_Ratio': X,
}
# Convert to DataFrame
results_df = pd.DataFrame(DATA_TOTAL)
# Export the DataFrame to an Excel file
results_df.to_excel('INELASTIC_FREE_VIBRATION_DAMPING_RATIO_SDOF_OPTIMIZATION_RESULTS.xlsx', index=False)
#------------------------------------------------------------------------------------------------

