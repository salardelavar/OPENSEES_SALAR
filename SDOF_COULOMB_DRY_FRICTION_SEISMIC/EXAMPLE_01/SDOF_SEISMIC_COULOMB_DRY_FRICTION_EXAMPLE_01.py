######################################################################################################################
#                          >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                          #
#            SIMULATE THE SEISMIC RESPONSE OF A SINGLE-DEGREE-OF-FREEDOM (SDOF) STRUCTURE WITH COULOMB               #
#                                                DRY FRICTION USING OPENSEES                                         #
#--------------------------------------------------------------------------------------------------------------------#
#                                                         EXAMPLE 01                                                 #
#--------------------------------------------------------------------------------------------------------------------#
#                              THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                            #
#                                       EMAIL: salar.d.ghashghaei@gmail.com                                          #
######################################################################################################################
"""
This Python code simulates a single-degree-of-freedom (SDOF) system with both a linear spring and a Coulomb friction damper.
The system's mass is 4000 tons, and the structural spring stiffness is 10000 kN/m. The Coulomb damper is modeled using two different approaches: a `Steel01` material for an elasto-plastic representation and a dedicated `CoulombDamper` material.
A transient analysis is performed, subjecting the system to a sinusoidal acceleration input of 0.4g at 1.5 Hz for 2000 steps with a time step of 0.01 seconds.
The code records and plots the displacement, damper forces from both models, and hysteresis loops against displacement for comparison. This allows for analysis of how each damping model behaves under seismic excitation..
"""
# WIKIPEDIA
'https://en.wikipedia.org/wiki/Friction#:~:text=Coulomb%20friction%2C%20named%20after%20Charles,to%20the%20net%20applied%20force.'
# BOOK: Differential Equations for Engineers - Wei-Chau Xie - CAMBRIDGE
'https://www.cambridge.org/core/books/differential-equations-for-engineers/1B8F1A62BF6F98EB972CFE9114FA8B84'
# portwooddigital Website by Michael H. Scott
'https://portwooddigital.com/2024/02/10/two-sprung-masses-and-some-friction-force/'
#%%------------------------------------------------------------------------------------------------
import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
import time as TI
import ANALYSIS_FUNCTION as S01
import PERIOD_FUN as S02
import DAMPING_RATIO_FUN as S05
import EIGENVALUE_ANALYSIS_FUN as S06
#%%------------------------------------------------------------------------------------------------
# Define parameters (units: m, N)
FY = 85000.0                                     # [N] Yield Force of Structure
FU = 1.5 * FY                                    # [N] Ultimate Force of Structure
Ke = 4500000.0                                   # [N/m] Spring Elastic Stiffness
DY = FY / Ke                                     # [m] Yield Displacement
DSU = 0.36                                       # [m] Ultimate Displacement
Ksh = (FU - FY) / (DSU - DY)                     # [N/m] Displacement Hardening Modulus
Kp = FU / DSU                                    # [N/m] Spring Plastic Stiffness
b = Ksh / Ke                                     # Displacement Hardening Ratio

M = 50000.0       # [kg] Mass
zi = 0.05         # Damping ratio
duration = 10.0   # [s] Analysis duration
dt = 0.005        # [s] Time step

mu = 0.25         # Coefficient of kinetic friction
#%%------------------------------------------------------------------------------------------------
# DEFINE ANALYSIS PROPERTIES
MAX_ITERATIONS = 20000    # Convergence iteration for test
MAX_TOLERANCE = 1.0e-10   # Convergence tolerance for test
#SPRING_TYPE: 1 -> 'ELASTIC'
#SPRING_TYPE: 2 -> 'INELASTIC'

#%%------------------------------------------------------------------------------------------------
# Positive branch points
pos_disp = [0, DY, DSU, 1.1*DSU, 1.25*DSU]
pos_force = [0, FY, FU, 0.2*FU, 0.1*FU]
KP = np.array([FY, DY, FU, DSU, 0.2*FU, 1.1*DSU, 0.1*FU, 1.25*DSU])

# Negative branch points
neg_disp = [0, -DY, -DSU, -1.1*DSU, -1.25*DSU]
neg_force = [0, -FY, -FU, -0.2*FU, -0.1*FU]
KN = np.array([-FY, -DY, -FU, -DSU, -0.2*FU, -1.1*DSU, -0.1*FU, -1.25*DSU])
#%%------------------------------------------------------------------------------------------------
# Plot
plt.plot(pos_disp, pos_force, marker='o', color='red')
plt.plot(neg_disp, neg_force, marker='o', color='black')

plt.xlabel("Displacement [m]")
plt.ylabel("Force [N]")
plt.title("Force–Displacement Curve")
plt.grid(True)
plt.axhline(0, linewidth=0.5)
plt.axvline(0, linewidth=0.5)
plt.show()
# ELASIC PERIOD:
ELAS_PERIOD = 2*np.pi * np.sqrt(M/Ke)
print(f'ELASIC PERIOD: {ELAS_PERIOD:.3f} (s)')     
# PLASIC PERIOD:
PLAS_PERIOD = 2*np.pi * np.sqrt(M/Kp)  
print(f'PLASIC PERIOD: {PLAS_PERIOD:.3f} (s)') 

# INITIAL MASS FOR RESPONSE SPECTRUM ANALYSIS
mi = (PLAS_PERIOD/2*np.pi)**2 * Kp 
print(mi)
#%%------------------------------------------------------------------------------------------------
# Calculate Over Strength Coefficient (Ω0)
Omega_0 = FU / FY
# Calculate Displacement Ductility Ratio (μ)
mu = DSU / DY
# Calculate Ductility Coefficient (Rμ)
R_mu = (2 * mu - 1) ** 0.5 / mu ** 0.5
# Calculate Structural Behavior Coefficient (R)
R = Omega_0 * R_mu
print(f'Over Strength Coefficient (Ω0):      {Omega_0:.2f}')
print(f'Displacement Ductility Ratio (μ):    {mu:.2f}')
print(f'Ductility Coefficient (Rμ):          {R_mu:.2f}')
print(f'Structural Behavior Coefficient (R): {R:.2f}')
#%%    
#omega = np.sqrt(Ke/M)    # [N/m] Elastic Stiffness
#C = 2 * zi * omega * M   # Damping Coefficient
alpha = 1.0               # For simulation of friction, we add alpha = 0
Fr = mu * M * 9.81        # the Friction Force Opposite to Displacement Rate.
#%% ---------------------------------------------
def SEISMIC_SDOF(SPRING_TYPE, N_STEPS, duration, dt):
    #%% Create model
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)
    
    #%% Add nodes
    ops.node(0, 0.0)
    ops.node(1, 0.0)
    ops.node(2, 0.0)
    
    #%% Fix nodes
    ops.fix(1, 1)
    
    #%% Define material and element
    MatTag01, MatTag02 = 1, 2
    if SPRING_TYPE == 'ELASTIC': 
        ops.uniaxialMaterial('Elastic', MatTag01, Ke)
        ops.element('zeroLength', 1, 1, 2,'-mat', MatTag01,'-dir', 1) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/ZeroLength.html
        # Define Friction (Fr = μ × W)
        ops.uniaxialMaterial('CoulombDamper', 10, 0.0, Fr)
        # INFO LINK: https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/CoulombDamper.html
        ops.element('zeroLength', 20, 0, 2,'-mat', 10, '-dir', 1)
    if SPRING_TYPE == 'INELASTIC':
        ops.uniaxialMaterial('HystereticSM', MatTag01, '-posEnv', *KP.flatten(), '-negEnv', *KN.flatten(), '-pinch', 1, 1)
        #INFO LINK: https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/HystereticSM.html
        #ops.uniaxialMaterial('Viscous', MatTag02, Cd, alpha)  # Material for Cd (alpha=1.0 for linear)
        ops.element('zeroLength', 1, 1, 2,'-mat', MatTag01, '-dir', 1)
        # Define Friction (Fr = μ × W)
        ops.uniaxialMaterial('CoulombDamper', 10, 0.0, Fr)
        # INFO LINK: https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/CoulombDamper.html
        ops.element('zeroLength', 20, 0, 2, '-mat', 10, '-dir', 1)
    
    
    # Define mass to node 2
    ops.mass(2, M)
    
    # Define analysis
    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.system('BandGeneral')
    # ops.system('BandSPD')
    # ops.system('UmfPack')
    # ops.system('ProfileSPD')
    ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/normDispIncr.html
    #ops.integrator('CentralDifference')  # JUST FOR LINEAR ANALYSIS - INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/centralDifference.html
    ALPHA=0.60; BETA=0.3025;
    ops.integrator('Newmark', ALPHA, BETA) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/newmark.html
    #alpha=2/3;gamma=1.5-alpha; gamma=1.5-alpha;beta=(2-alpha)**2/4;
    #ops.integrator('HHT', alpha, gamma, beta) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/hht.html
    ops.algorithm('Newton')  # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/algorithm.html
    ops.analysis('Transient') # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/analysis.html
    
    #ops.timeSeries('Linear', 1)
    #ops.pattern('Plain', 1, 1)
    #ops.load(2, 1.0)

    TIME = np.arange(0, N_STEPS*dt, dt)
    ACCEL = 9.81 * 0.01 * np.sin(2*np.pi*1.5*TIME) # 0.4g, 1.5 Hz
    ops.timeSeries('Path', 1, '-dt', dt, '-values', *ACCEL.tolist())
    ops.pattern('UniformExcitation', 1, 1, '-accel', 1)
    
    # Plot Sesimic Ground Motion
    plt.figure(00, figsize=(12, 8))
    plt.plot(TIME, ACCEL, color='black', linewidth=4)
    plt.title('Ground Motion Seismic Acceleration')
    plt.xlabel('Time [s]')
    plt.ylabel('Acceleration [m/s^2]')
    plt.grid()
    plt.show()
    
    # Perform analysis
    time = []
    disp = []
    vel = []
    accel = []
    reaction = []
    stiffness = []
    OMEGA, PERIOD = [], []
    force1 = []; force2 = [];
    FRICTION_FORCE = []
    PERIOD_MIN, PERIOD_MAX = [], []
    DI = []
    stable = 0
    current_time = 0.0
    
    while stable == 0 and current_time < duration:
        ops.analyze(1, dt)
        S01.ANALYSIS(stable, 1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS
        current_time = ops.getTime()
        time.append(current_time)
        ops.reactions()
        reaction.append(ops.nodeReaction(1, 1))  # BASE REACTION
        disp.append(ops.nodeDisp(2, 1))          # DISPLACEMENT  
        vel.append(ops.nodeVel(2, 1))            # VELOCITY
        accel.append(ops.nodeAccel(2, 1))        # ACCELERATION
        stiffness.append(np.abs(reaction[-1]) / np.abs(disp[-1]))
        OMEGA.append(np.sqrt(stiffness[-1]/M))
        PERIOD.append((np.pi * 2) / OMEGA[-1])
        FRICTION_FORCE.append(ops.eleResponse(20, 'force')[1])
        print(time[-1], disp[-1], vel[-1])
        #%% IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
        PERIODmin, PERIODmax = S06.EIGENVALUE_ANALYSIS(1, PLOT=True)
        PERIOD_MIN.append(PERIODmin)
        PERIOD_MAX.append(PERIODmax)
        DI.append(100*(np.abs(disp[-1])-DY)/(DSU-DY)) # DAMAGE INDEX
        if DI[-1] <= 0.0: # DAMAGE INDEX -> 0.0 <= DI <=100 
            DI[-1] = 0.0
        if DI[-1] >= 100: 
            DI[-1] = 100.0
        
    # Mean Period 
    deltaT_mean = np.mean(PERIOD)
    #deltaT_median = np.median(PERIOD)
    print(f"Min. Period:        {np.min(PERIOD):.8e} [s]") 
    print(f"Mean Period:        {deltaT_mean:.8e} [s]")   
    print(f"Max. Period:        {np.max(PERIOD):.8e} [s]")   

    # Calculate Damping Ratio of SDOF the system
    damping_ratio = S05.DAMPING_RATIO(disp) 
    
    # Compute modal properties
    ops.modalProperties("-print", "-file", "SALAR_ModalReport.txt", "-unorm") 
    
    DATA = (time, reaction, disp, vel, accel, stiffness,
            PERIOD, deltaT_mean, damping_ratio,
            force1, force2,
            np.array(PERIOD_MIN), np.array(PERIOD_MAX), np.array(DI), np.array(FRICTION_FORCE))
    
    return DATA

#%% ---------------------------------------------
# Analysis Durations for Dynamic Analysis:
starttime = TI.process_time()



SPRING_TYPE = 'ELASTIC'
N_STEPS = 2000

DATA = SEISMIC_SDOF(SPRING_TYPE, N_STEPS, duration, dt)
(time, reactionE, dispE, velE, accelE, stiffnessE,
 periodE, E_periodE, E_damping_ratioE,
 force1E, force2E,
 period_minE, period_maxE, diE, FFe) = DATA

#S02.PERIOD_FUN(dispE, dt)

SPRING_TYPE = 'INELASTIC'
N_STEPS = 2000

DATA = SEISMIC_SDOF(SPRING_TYPE, N_STEPS, duration, dt)
(time, reactionI, dispI, velI, accelI, stiffnessI,
 periodI, E_periodI, E_damping_ratioI,
 force1I, force2I,
 period_minI, period_maxI, diI, FFi) = DATA

#S02.PERIOD_FUN(dispI, dt)

totaltime = TI.process_time() - starttime
print(f'\nTotal Analysis Durations (s): {totaltime:.4f} \n\n')
#%% ---------------------------------------------
print("Elastic  max |friction| =", np.max(np.abs(FFe)))
print("Inelastic max |friction| =", np.max(np.abs(FFi)))
print("Max |velocity| elastic   =", np.max(np.abs(velE)))
print("Max |velocity| inelastic =", np.max(np.abs(velI)))
#%% ---------------------------------------------
# PLOT STRUCTURAL PERIOD DURING THE ANALYSIS
plt.figure(0, figsize=(12, 8))
plt.plot(dispE, period_minE, color='black', linewidth=4)
plt.plot(dispE, period_maxE, color='red', linewidth=4)
plt.plot(dispI, period_minI, color='green', linewidth=4)
plt.plot(dispI, period_maxI, color='purple', linewidth=4)
plt.title('Period of Structure vs Displacement During Free-Vibration Analysis')
plt.ylabel('Structural Period [s]')
plt.xlabel('Displacement [m]')
#plt.semilogy()
plt.grid()
plt.legend([f'ELASTIC PERIOD - MIN VALUES: Min: {np.min(period_minE):.3f} (s) - Mean: {np.mean(period_minE):.3f} (s) - Max: {np.max(period_minE):.3f} (s)', 
            f'ELASTIC PERIOD - MAX VALUES:  Min: {np.min(period_maxE):.3f} (s) - Mean: {np.mean(period_maxE):.3f} (s) - Max: {np.max(period_maxE):.3f} (s)',
            f'INELASTIC PERIOD - MIN VALUES: Min: {np.min(period_minI):.3f} (s) - Mean: {np.mean(period_minI):.3f} (s) - Max: {np.max(period_minI):.3f} (s)', 
            f'INELASTIC PERIOD - MAX VALUES:  Min: {np.min(period_maxI):.3f} (s) - Mean: {np.mean(period_maxI):.3f} (s) - Max: {np.max(period_maxI):.3f} (s)',
            ])
plt.show()
#%% ---------------------------------------------
# Plot Results
plt.figure(2, figsize=(18, 14))

# Define colors
elastic_color = 'black'
inelastic_color = 'red'

# Friction Force plot
plt.subplot(7, 1, 1)
plt.plot(time, FFe, color=elastic_color, linewidth=1.5)
plt.plot(time, FFi, color=inelastic_color, linewidth=1.5)
plt.title('Friction Forces vs Time', fontsize=12, pad=10)
plt.ylabel('Friction (N)', fontsize=10)
plt.grid(alpha=0.3)

# Reaction plot
plt.subplot(7, 1, 2)
plt.plot(time, reactionE, color=elastic_color, linewidth=1.5)
plt.plot(time, reactionI, color=inelastic_color, linewidth=1.5)
plt.title('Reaction Forces vs Time', fontsize=12, pad=10)
plt.ylabel('Reaction (N)', fontsize=10)
plt.grid(alpha=0.3)


# Displacement plot
plt.subplot(7, 1, 3)
plt.plot(time, dispE, color=elastic_color, linewidth=1.5, label=f'Elastic Damping Ratio: {E_damping_ratioE:.3e} %')
plt.plot(time, dispI, color=inelastic_color, linewidth=1.5, label=f'Inelastic Damping Ratio: {E_damping_ratioI:.3e} %')
plt.title('Displacement vs Time', fontsize=12, pad=10)
plt.ylabel('Displacement (m)', fontsize=10)
plt.grid(alpha=0.3)
plt.legend(loc='upper right', framealpha=1)

# Velocity plot
plt.subplot(7, 1, 4)
plt.plot(time, velE, color=elastic_color, linewidth=1.5)
plt.plot(time, velI, color=inelastic_color, linewidth=1.5)
plt.title('Velocity vs Time', fontsize=12, pad=10)
plt.ylabel('Velocity (m/s)', fontsize=10)
plt.grid(alpha=0.3)

# Acceleration plot
plt.subplot(7, 1, 5)
plt.plot(time, accelE, color=elastic_color, linewidth=1.5)
plt.plot(time, accelI, color=inelastic_color, linewidth=1.5)
plt.title('Acceleration vs Time', fontsize=12, pad=10)
plt.ylabel('Acceleration (m/s²)', fontsize=10)
plt.grid(alpha=0.3)

# Stiffness plot
plt.subplot(7, 1, 6)
plt.plot(time, stiffnessE, color=elastic_color, linewidth=1.5)
plt.plot(time, stiffnessI, color=inelastic_color, linewidth=1.5)
plt.title('Stiffness vs Time', fontsize=12, pad=10)
plt.ylabel('Stiffness (N/m)', fontsize=10)
plt.semilogy()
plt.grid(alpha=0.3)

# Period plot
plt.subplot(7, 1, 7)
plt.plot(time, periodE, color=elastic_color, linewidth=1.5, label=f'Elastic Period: {E_periodE:.3e}')
plt.plot(time, periodI, color=inelastic_color, linewidth=1.5, label=f'Inelastic Period: {E_periodI:.3e}')
plt.title('Period vs Time', fontsize=12, pad=10)
plt.ylabel('Period  (s)', fontsize=10)
plt.xlabel('Time (s)', fontsize=10)
plt.legend(loc='upper right', framealpha=1)
plt.semilogy()
plt.grid(alpha=0.3)

plt.tight_layout()
plt.show()

plt.figure(2, figsize=(8, 6))
plt.plot(dispE, FFe, color='black', linewidth=2)
plt.plot(dispI, FFi, color='purple', linewidth=2)
plt.xlabel('Displacement [m]')
plt.ylabel('Friction Force [N]')
plt.title(f'Displacement vs Friction Force')
plt.legend(['ELASTIC', 'INELASTIC'])
plt.grid()
plt.show()

plt.figure(3, figsize=(8, 6))
plt.plot(dispE, reactionE, color='black', linewidth=2)
plt.plot(dispI, reactionI, color='purple', linewidth=2)
plt.xlabel('Displacement [m]')
plt.ylabel('Base-reaction [N]')
plt.title(f'Displacement vs Base-reaction')
plt.legend(['ELASTIC', 'INELASTIC'])
plt.grid()
plt.show()

plt.figure(4, figsize=(8, 6))
plt.plot(dispE, diE, color='black', linewidth=2)
plt.plot(dispI, diI, color='purple', linewidth=2)
plt.xlabel('Displacement [m]')
plt.ylabel('Structural Damage Index [%]')
plt.title(f'Displacement vs Structural Damage Index')
plt.legend(['ELASTIC', 'INELASTIC'])
plt.grid()
plt.show()
#%% ---------------------------------------------
# Print out the state of nodes 1 and 2
ops.printModel("node",1, 2)
# Print out the state of element 1
ops.printModel("ele", 1)
# Print the Model
#printModel()
ops.printModel("-JSON", "-file", "SDOF_SEISMIC_COULOMB_DRY_FRICTION.json")
#%%-------------------------------------------------------------------------------
# EXCEL OUTPUT
import pandas as pd

# Create DataFrame function
def create_df(dispE, dispI, velE, velI, accelE, accelI, reactionE, reactionI, diE, diI):
    df = pd.DataFrame({
        "DISPLACEMENT [m] - ELASTIC": dispE,
        "VELOCITY [m/s] - ELASTIC": velE,
        "ACCELERATION [m/s^2] - ELASTIC": accelE,
        "BASE-REACTION [N] - ELASTIC": reactionE,
        "DISPLACEMENT [m] - INELASTIC": dispI,
        "VELOCITY [m/s] - INELASTIC": velI,
        "ACCELERATION [m/s^2] - INELASTIC": accelI,
        "BASE-REACTION [N] - INELASTIC": reactionI,
        "DAMAGE-INDEX [%] - INELASTIC": diI,
    })
    return df


# Save to Excel
with pd.ExcelWriter("SDOF_SEISMIC_COULOMB_DRY_FRICTION_OUTPUT.xlsx", engine='openpyxl') as writer:
    
    df1 = create_df(dispE, dispI, velE, velI, accelE, accelI, reactionE, reactionI, diE, diI)
    df1.to_excel(writer, sheet_name="OUTPUT", index=False)
#%%------------------------------------------------------------------------------------------------
