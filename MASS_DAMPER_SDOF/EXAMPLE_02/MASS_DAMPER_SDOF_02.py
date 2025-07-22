#####################################################################################
#                                  IN THE NAME OF ALLAH                             #
#       DYNAMIC RESPONSE ANALYSIS OF A SINGLE-DEGREE-OF-FREEDOM (SDOF) SYSTEM       #
#           WITH ACTIVE MASS DAMPER (AMD) UNDER FREE VIBRATION CONDITIONS           #                                    
#           INELASTIC SPRING WITH INTIAL VELOCITY AND TWO SEISMIC LOADING           #
#-----------------------------------------------------------------------------------#
#              THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)           #
#                       EMAIL: salar.d.ghashghaei@gmail.com                         #
#####################################################################################
"""
- A 1D single-degree-of-freedom model with a lumped mass, a nonlinear “MultiLinear” spring (yielding at 240 N),
    and optional 5 % Rayleigh damping.  
- Two real earthquake records are applied in sequence (at 20 s and 50 s), and the Newmark average-acceleration method
    integrates the nonlinear equations of motion.  
- The damped case shows faster decay of vibrations and lower peak displacements/accelerations than the undamped system.  
- Nonlinear yielding leads to stiffness degradation—seen as period elongation—and produces residual displacements, with
    hysteresis loops dissipating energy.  
- These results inform performance-based design and retrofitting strategies, helping engineers limit inelastic deformations
    and improve seismic resilience.
"""
#$$-------------------------------------------------------------------------
# Load the image
def PLOT_IMAGE(image):
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg
    image = mpimg.imread(image_path)

    # Display the image
    plt.figure(figsize=(10, 6))
    plt.imshow(image)
    plt.axis('off')  # Hide axes
    plt.show()
    
image_path = 'MASS_DAMPER_SDOF_02.png'    
PLOT_IMAGE(image_path)
#$$-------------------------------------------------------------------------
"""
When OK equals -1, it generally indicates that the command or operation was not executed
because it was already in progress or had already been completed. This can happen if you
try to run a command that is already running or has been completed in a previous step.

When OK equals -2, it typically indicates that the command or operation was not executed
because it was not recognized or not implemented. This could mean that the command
is either misspelled, not available in the current version of OpenSees, or not applicable to the current context.

When OK equals -3, it typically means that the command or operation failed.
This could be due to various reasons, such as incorrect input parameters,
syntax errors, or issues with the model setup.
"""
def ANALYSIS(OK, INCREMENT, TOLERANCE, MAX_ITERAIONS):
    import openseespy.opensees as op
    test = {1:'NormDispIncr', 2: 'RelativeEnergyIncr', 4: 'RelativeNormUnbalance',5: 'RelativeNormDispIncr', 6: 'NormUnbalance'}
    algorithm = {1:'KrylovNewton', 2: 'SecantNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'BFGS', 7: 'Broyden', 8: 'NewtonLineSearch'}

    for i in test:
        for j in algorithm:
            if OK != 0:
                if j < 4:
                    op.algorithm(algorithm[j], '-initial')

                else:
                    op.algorithm(algorithm[j])

                op.test(test[i], TOLERANCE, MAX_ITERAIONS) 
                OK = op.analyze(INCREMENT)                            
                print(test[i], algorithm[j], OK)             
                if OK == 0:
                    break
            else:
                continue
                
# ---------------------------- 

def PLOT_2D(X1, Y1, X2, Y2, XLABEL, YLABEL, TITLE):
    plt.figure(figsize=(10, 6))
    plt.plot(X1, Y1, label='Damped', color='black')
    plt.plot(X2, Y2, label='Undamped', color='grey', linestyle='--')
    plt.xlabel(XLABEL)
    plt.ylabel(YLABEL)
    plt.title(TITLE)
    plt.grid(True)
    #plt.semilogy()
    plt.legend()
    plt.show()
    
# ----------------------------

def PLOT_SPRING(X1, Y1, XLABEL, YLABEL, TITLE):
    plt.figure(figsize=(8, 6))
    plt.plot(X1, Y1, color='black')
    plt.xlabel(XLABEL)
    plt.ylabel(YLABEL)
    plt.title(TITLE)
    plt.grid(True)
    #plt.semilogy()
    plt.show()  
        
# ----------------------------

def plot_chart(time_damped, base_reaction_damped, displacement_damped, velocity_damped, acceleration_damped, omega_damped,
               time_undamped, base_reaction_undamped, displacement_undamped, velocity_undamped, acceleration_undamped, omega_undamped):
    fig, axs = plt.subplots(5, 1, figsize=(10, 12))

    # Plot base reaction
    axs[0].plot(time_damped, base_reaction_damped, label='Damped', color='blue')
    axs[0].plot(time_undamped, base_reaction_undamped, label='Undamped', color='cyan', linestyle='--')
    axs[0].set_title(f'Base Reaction vs Time - Damped Max Abs: {np.max(np.abs(base_reaction_damped)):.5f}, Undamped Max Abs: {np.max(np.abs(base_reaction_undamped)):.5f}')
    axs[0].set_xlabel('Time [s]')
    axs[0].set_ylabel('Base Reaction [N]')
    axs[0].grid(True)
    axs[0].legend()

    # Plot displacement
    axs[1].plot(time_damped, displacement_damped, label='Damped', color='red')
    axs[1].plot(time_undamped, displacement_undamped, label='Undamped', color='orange', linestyle='--')
    axs[1].set_title(f'Displacement vs Time - Damped Max Abs: {np.max(np.abs(displacement_damped)):.5f}, Undamped Max Abs: {np.max(np.abs(displacement_undamped)):.5f}')
    axs[1].set_xlabel('Time [s]')
    axs[1].set_ylabel('Displacement [m]')
    axs[1].grid(True)
    axs[1].legend()

    # Plot velocity
    axs[2].plot(time_damped, velocity_damped, label='Damped', color='green')
    axs[2].plot(time_undamped, velocity_undamped, label='Undamped', color='lime', linestyle='--')
    axs[2].set_title(f'Velocity vs Time - Damped Max Abs: {np.max(np.abs(velocity_damped)):.5f}, Undamped Max Abs: {np.max(np.abs(velocity_undamped)):.5f}')
    axs[2].set_xlabel('Time [s]')
    axs[2].set_ylabel('Velocity [m/s]')
    axs[2].grid(True)
    axs[2].legend()

    # Plot acceleration
    axs[3].plot(time_damped, acceleration_damped, label='Damped', color='purple')
    axs[3].plot(time_undamped, acceleration_undamped, label='Undamped', color='magenta', linestyle='--')
    axs[3].set_title(f'Acceleration vs Time - Damped Max Abs: {np.max(np.abs(acceleration_damped)):.5f}, Undamped Max Abs: {np.max(np.abs(acceleration_undamped)):.5f}')
    axs[3].set_xlabel('Time [s]')
    axs[3].set_ylabel('Acceleration [m/s²]')
    axs[3].grid(True)
    axs[3].legend()

    # Plot frequncy
    axs[4].plot(time_damped, omega_damped, label='Damped', color='brown')
    axs[4].plot(time_undamped, omega_undamped, label='Undamped', color='grey', linestyle='--')
    axs[4].set_title(f'Natural Fequncy vs Time - Damped Max Abs: {np.max(np.abs(omega_damped)):.5f}, Undamped Max Abs: {np.max(np.abs(omega_undamped)):.5f}')
    axs[4].set_xlabel('Time [s]')
    axs[4].set_ylabel('Natural Frequncy [Hertz]')
    axs[4].grid(True)
    axs[4].semilogy()
    axs[4].legend()

    plt.tight_layout()
    plt.show()    
#$$-------------------------------------------------------------------------
import openseespy.opensees as ops
import matplotlib.pyplot as plt
import numpy as np

def run_analysis(m1, m2, k1, k2, ks, zeta1, zeta2, zetaS, dt, Tfinal, iv0, damped=True):
    # Wipe existing model
    ops.wipe()

    # Define the model with 1 dimension and 1 degree of freedom per node
    ops.model('basic', '-ndm', 1, '-ndf', 1)
    GMfact = 9.81 # standard acceleration of gravity or standard acceleration
    # Natural frequency (rad/s)
    wn1 = (k1 / m2) ** 0.5
    wn2 = (k2 / m2) ** 0.5
    wnS = (ks / m1) ** 0.5
    # Damping coefficient (Ns/m)
    c1 = 2 * wn1 * m2 * zeta1 if damped else 0.0
    c2 = 2 * wn2 * m2 * zeta2 if damped else 0.0
    cs = 2 * wnS * m1 * zetaS if damped else 0.0

    # Define nodes and boundary conditions
    ops.node(1, 0.0)
    ops.fix(1, 1)
    ops.node(2, 0.0)
    ops.mass(2, m1) 
    ops.node(3, 0.0)
    ops.mass(3, m2)

    MAX_ITERATIONS = 5000  # convergence iteration for test
    TOLERANCE = 1.0e-10    # convergence tolerance for test
    
 
    ops.uniaxialMaterial('MultiLinear', 1, *kh1.flatten()) # Horizontal spring 01
    ops.uniaxialMaterial('MultiLinear', 2, *kh2.flatten()) # Horizontal spring 02
    ops.uniaxialMaterial('MultiLinear', 3, *khS.flatten()) # Horizontal spring - Structure
    
    # Define materials for damper
    ops.uniaxialMaterial('Elastic', 4, 0.0, c1)
    ops.uniaxialMaterial('Elastic', 5, 0.0, c2)
    ops.uniaxialMaterial('Elastic', 6, 0.0, cs)
    
    
    # Define elements
    ops.element('zeroLength', 1, 1, 2, '-mat', 3, 6, '-dir', 1, 1)
    ops.element('zeroLength', 2, 2, 3, '-mat', 1, 3, '-dir', 1, 1)
    ops.element('zeroLength', 3, 2, 3, '-mat', 2, 4, '-dir', 1, 1)
    
    # Dynamic analysis setup
    ops.constraints('Transformation')
    ops.numberer('Plain')
    ops.system('UmfPack')
    ops.test('EnergyIncr', TOLERANCE, MAX_ITERATIONS)
    #ops.integrator('CentralDifference')
    #ops.integrator('HHT', 0.9)
    ops.integrator('Newmark', 0.5, 0.25)
    ops.algorithm('ModifiedNewton')
    
    # Rayleigh damping
    Lambda01 = ops.eigen('-fullGenLapack', 1)
    Omega01 = np.power(max(Lambda01), 0.5)
    a0 = 2 * zetaS * Omega01  # Mass-proportional damping
    a1 = 2 * zetaS / Omega01  # Stiffness-proportional damping
    ops.rayleigh(a0, a1, 0, 0)
    PERIOD = 2 * np.pi / Omega01
    print(f'Peroid of Stucture: {PERIOD:.4f}\n')
    
    # Define time series for input motion (Acceleration time history)
    ops.timeSeries('Path', 1, '-dt', 0.01, '-filePath', 'OPENSEES_SPRING_SEISMIC_01.txt', '-factor', GMfact, '-startTime', 20)
    ops.timeSeries('Path', 2, '-dt', 0.01, '-filePath', 'OPENSEES_SPRING_SEISMIC_02.txt', '-factor', GMfact, '-startTime', 50)

    # Define load patterns
    # pattern UniformExcitation $patternTag $dof -accel $tsTag <-vel0 $vel0> <-fact $cFact>
    ops.pattern('UniformExcitation', 1, 1, '-accel', 1, '-vel0', iv0, '-fact', 1.0)
    ops.pattern('UniformExcitation', 2, 1, '-accel', 2)
    
    # Define analysis type
    ops.analysis('Transient')

    # Time Integration Parameters
    Nsteps = int(Tfinal / dt)

    # Lists to Store Results
    time = []
    base_reaction = []
    displacement = []
    velocity = []
    acceleration = []
    fr = []

    # Perform Analysis
    for i in range(Nsteps):
        OK = ops.analyze(1, dt)
        ANALYSIS(OK, 1, TOLERANCE, MAX_ITERATIONS)
        time.append(ops.getTime())
        base_reaction.append(-ops.eleResponse(1, 'force')[0] - ops.eleResponse(2, 'force')[0])  # AXIAL REACTION
        displacement.append(ops.nodeDisp(3, 1))
        velocity.append(ops.nodeVel(3, 1))
        acceleration.append(ops.nodeAccel(3, 1))
        KE = np.abs(base_reaction[-1] / displacement[-1]) # Calculate structure stiffness
        omega =  np.sqrt(KE / m1) # Calculate angular frequency (omega)
        T = (2 * np.pi) / omega  # Calculate period (T)
        f = 1 / T                # Calculate natural frequency
        fr.append(f)

    return time, base_reaction, displacement, velocity, acceleration, fr


# Parameters
m1 = 1000.0  # [kg] Spring Mass 01
m2 = 20.0    # [kg] Spring Mass 02
k1 = 0.1957  # [N/m] Spring Stiffness 01
k2 = 0.31    # [N/m] Spring Stiffness 02
ks = 1.527   # [N/m] Spring Stiffness - Structure
Fy1 = 4      # [N] Yield strength 01
Fy2 = 8.5    # [N] Yield strength 02
Fys = 5.70   # [N] Yield strength - Structure
zeta1 = 0.05   # Damping ratio 01
zeta2 = 0.02   # Damping ratio 02
zetaS = 0.02   # Damping ratio - Structure
dt = 0.01      # Time step in seconds
Tfinal = 75.0  # [s] Total time in seconds
iv0 = 0.0005   # [m/s] Initial velocity 

# Define the MultiLinear material properties for springs
d1 = Fy1 * k1;
kh1 = np.array([[d1, Fy1],
                [1.3*d1, 1.8*Fy1],
                [1.6*d1, 1.5*Fy1],
                [1.6*d1, 1.5*Fy1],
                [2.1*d1, 1.2*Fy1],
                [2.5*d1, 0.7*Fy1]]) # SPRING 01 FORCE-DISPLACEMENT REALTIONS

d2 = Fy2 * k2;
kh2 = np.array([[d2, Fy2],
                [1.3*d2, 1.8*Fy2],
                [1.6*d2, 1.5*Fy2],
                [1.6*d2, 1.5*Fy2],
                [2.1*d2, 1.2*Fy2],
                [2.5*d2, 0.7*Fy2]]) # SPRING 02 FORCE-DISPLACEMENT REALTIONS

ds = Fys * ks;
khS = np.array([[ds, Fys],
                [1.3*ds, 1.8*Fys],
                [1.6*ds, 1.5*Fys],
                [1.6*ds, 1.5*Fys],
                [2.1*ds, 1.2*Fys],
                [2.5*ds, 0.7*Fys]]) # STRUCTURE FORCE-DISPLACEMENT REALTIONS

# PLOT SPRING FORCE-DISPLACEMENT RELATION
displacement_kh, force_kh = kh1[:, 0], kh1[:, 1]
X1 = displacement_kh
Y1 = force_kh
XLABEL = 'Displacement'
YLABEL = 'Force'
TITLE = 'Spring 01 Force and Displacement Relation'
PLOT_SPRING(X1, Y1, XLABEL, YLABEL, TITLE)

# PLOT SPRING FORCE-DISPLACEMENT RELATION
displacement_kh, force_kh = kh2[:, 0], kh2[:, 1]
X1 = displacement_kh
Y1 = force_kh
XLABEL = 'Displacement'
YLABEL = 'Force'
TITLE = 'Spring 02 Force and Displacement Relation'
PLOT_SPRING(X1, Y1, XLABEL, YLABEL, TITLE)

# PLOT SPRING FORCE-DISPLACEMENT RELATION
displacement_kh, force_kh = khS[:, 0], khS[:, 1]
X1 = displacement_kh
Y1 = force_kh
XLABEL = 'Displacement'
YLABEL = 'Force'
TITLE = 'Structure Force and Displacement Relation'
PLOT_SPRING(X1, Y1, XLABEL, YLABEL, TITLE)

# Run analysis for damped and undamped cases
time_damped, base_reaction_damped, displacement_damped, velocity_damped, acceleration_damped, fr_damped = run_analysis(m1, m2, k1, k2, ks, zeta1, zeta2, zetaS, dt, Tfinal, iv0, damped=True)
time_undamped, base_reaction_undamped, displacement_undamped, velocity_undamped, acceleration_undamped, fr_undamped = run_analysis(m1, m2, k1, k2, ks, zeta1, zeta2, zetaS, dt, Tfinal, iv0, damped=False)


# Plot the results
plot_chart(time_damped, base_reaction_damped, displacement_damped, velocity_damped, acceleration_damped, fr_damped,
           time_undamped, base_reaction_undamped, displacement_undamped, velocity_undamped, acceleration_undamped, fr_undamped)

#$$-------------------------------------------------------------------------
### BASE REACTION & DISPALCEMENT:
X1, Y1 = base_reaction_damped, displacement_damped
X2, Y2 = base_reaction_undamped, displacement_undamped
XLABEL = 'Displacement'
YLABEL = 'Base-Reaction'
TITLE = 'Base-Reaction and Displacement'

PLOT_2D(X1, Y1, X2, Y2, XLABEL, YLABEL, TITLE)
#$$-------------------------------------------------------------------------   
# EXPORT DATA TO EXCEL
import pandas as pd
DATA_TOTAL = { 
    'base_reaction_undamped': base_reaction_undamped,
    'displacement_undamped': displacement_undamped,
    'velocity_undamped': velocity_undamped,
    'acceleration_undamped': acceleration_undamped,
    'fr_undamped': fr_damped,
    'base_reaction_damped': base_reaction_damped,
    'displacement_damped': displacement_damped,
    'velocity_damped': velocity_damped,
    'acceleration_damped': acceleration_damped,
    'fr_damped': fr_damped,
}
# Convert to DataFrame
results_df = pd.DataFrame(DATA_TOTAL)
# Export the DataFrame to an Excel file
results_df.to_excel('MASS_DAMPER_SDOF_02_RESULTS.xlsx', index=False) 
#%%------------------------------------------------------------------------------
# Print out the state of nodes 1 , 2 and 3
ops.printModel("node",1, 2, 3)
# Print out the state of element 1 , 2 and 3
ops.printModel("ele", 1, 2 , 3)
# Print the Model
#printModel()
ops.printModel("-JSON", "-file", "MASS_DAMPER_SDOF_02.json")
#%%-------------------------------------------------------------------------------