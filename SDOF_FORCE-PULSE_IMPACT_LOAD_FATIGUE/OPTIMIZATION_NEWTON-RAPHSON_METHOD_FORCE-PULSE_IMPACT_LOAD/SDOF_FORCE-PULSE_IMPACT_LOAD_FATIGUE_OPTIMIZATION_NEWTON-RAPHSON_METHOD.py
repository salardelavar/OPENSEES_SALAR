######################################################################################################################
#                          >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                          #
#   OPTIMIZATION OF A FORCE-PULSE IMPACT LOAD ANALYSIS OF A SINGLE-DEGREE-OF-FREEDOM (SDOF) SYSTEM USING OPENSEES    #
#                                                  WITH FATIGUE MATERIAL                                             #
#--------------------------------------------------------------------------------------------------------------------#
#                                  OPTIMIZATION ALGORITHM: NEWTON-RAPHSON METHOD                                     #
#--------------------------------------------------------------------------------------------------------------------#
#                              THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                            #
#                                       EMAIL: salar.d.ghashghaei@gmail.com                                          #
######################################################################################################################
"""
Nonlinear Dynamic Analysis of SDOF Systems: Hysteretic Behavior, Damping, and Fatigue Effects:
-----------------------------------------------------------------------    
1. OpenSeesPy-Based Simulation: Models a Single-Degree-of-Freedom (SDOF) system under transient force excitation.  
2. Material Nonlinearity: Implements asymmetric hysteresis (tension/compression) via `HystereticSM` material.  
3. Viscous Damping: Uses velocity-dependent damping (linear/nonlinear) with adjustable exponent (α).  
4. Fatigue Damage Tracking: Integrates Fatigue material model to assess cyclic degradation effects.  
5. Dynamic Analysis: Solves equations of motion via Newmark-β (γ=0.5, β=0.25) for stability.  
6. Response Comparison: Contrasts elastic vs. inelastic displacement, velocity, and acceleration time histories.  
7. System Identification: Estimates effective period (T) and damping ratio (ζ) via FFT and log decrement.  
8. Backbone Curves: Defines multi-linear force-displacement envelopes for nonlinear spring behavior.  
9. Support Reactions: Tracks base shear forces to evaluate structural demand.  
10. Visualization: Plots hysteresis loops, time-domain responses, and spectral characteristics.  

Key Insight: The analysis highlights period elongation, energy dissipation, and stiffness
 degradation in inelastic systems, critical for seismic design and performance assessment.  

The equation of motion for an SDOF system is
mü(t) + ců(t) + ku(t) = p(t)

Impulsive Loading Characteristics:
p(t) = P₀ for 0 ≤ t ≤ td
p(t) = 0 for t > td

"""
import openseespy.opensees as ops
import numpy as np
import time as TI
import matplotlib.pyplot as plt
import ANALYSIS_FUNCTION as S01
import ESTIMATE_T_ZETA_FUN as S02

#%% --------------------------
# FORCE PULSE IMPACT LOAD
durationF = 20.0   # [s] Force Pulse Analysis duration
dtF = 0.01         # [s] Time step
timeF = np.arange(0, durationF, dtF)

# Input force (time-dependent pulse)
F_input = np.zeros_like(timeF)
pulse_start = 10  # index for 0.1s
pulse_end = 150   # index for 1.5s

print(f'Pulse Start Time: {pulse_start*dtF} (s)')
print(f'Pulse End Time: {pulse_end*dtF} (s)')


peak_force = 3000.0  # [N] Maximum Impact Value
pulse_duration = pulse_end - pulse_start
#F_input[pulse_start:pulse_end] = peak_force * np.sin(np.pi * np.arange(pulse_duration)/pulse_duration)
F_input[pulse_start:pulse_end] = peak_force * 1.0

#%% DEFINE ANALYSIS PROPERTIES
MAX_ITERATIONS = 20000    # Convergence iteration for test
MAX_TOLERANCE = 1.0e-10   # Convergence tolerance for test
#SPRING_KIND: 1 -> 'ELASTIC'
#SPRING_KIND: 2 -> 'INELASTIC'

def SDOF_FORCE_PULSE(area, SPRING_KIND, duration, dt):
    # Tension Force-Displacement Relationship (Positive values)
    FY1, DY1 = 2.7720*area, 0.01    # First yield point
    FY2, DY2 = 3.1046*area, 0.02    # Peak force
    FY3, DY3 = 1.6632*area, 0.04    # Post-peak softening
    FY4, DY4 = 1.6632*area, 0.06    # Plateau
    FY5, DY5 = 0.2772*area, 0.28    # Further softening
    FY6, DY6 = 0.2000*area, 0.51    # Near-zero stiffness
    FY7, DY7 = 0.0*area, 0.72       # Zero force (failure)
    
    KP = np.array([FY1, DY1, FY2, DY2, FY3, DY3, FY4, DY4, FY5, DY5, FY6, DY6, FY7, DY7])

    # Compression Force-Displacement Relationship (Negative values)
    FY1n, DY1n = -2.7720*area, -0.01    # First yield in compression
    FY2n, DY2n = -3.1046*area, -0.05    # Peak compressive force
    FY3n, DY3n = -1.6632*area, -0.20    # Post-peak softening

    KN = np.array([FY1n, DY1n, FY2n, DY2n, FY3n, DY3n])
    
    #%%
    K = KP[0] / KP[1]         # [N/m] Elastic Stiffness

    omega = np.sqrt(K/M)
    alpha = 1.0               # velocity exponent (usually 0.3–1.0)
    Cd = 2 * zi * omega * M   # [N·s/m] Damping coefficient 
    
    ops.wipe()
    ops.model("Basic", "-ndm", 1, "-ndf", 1)

    # Nodes
    ops.node(1, 0.0)
    ops.node(2, 0.0)
    ops.fix(1, 1)  # fix base

    # Mass
    ops.mass(2, M)

    # Material
    if SPRING_KIND == 'ELASTIC':
        ops.uniaxialMaterial('Elastic', 1, K, Cd)
        ops.uniaxialMaterial('Fatigue', 2, 1, '-E0', 0.4, '-m', -0.458, '-min', -1e16, '-max', 1e16)
        ops.element("zeroLength", 1, 1, 2, '-mat', 1, 2, '-dir', 1, 2)
    elif SPRING_KIND == 'INELASTIC':
        ops.uniaxialMaterial('HystereticSM', 1, '-posEnv', *KP.flatten(), '-negEnv', *KN.flatten(), '-pinch', 1, 1)
        #INFO LINK: https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/HystereticSM.html
        ops.uniaxialMaterial('Viscous', 2, Cd, alpha)  # Material for C (alpha=1.0 for linear)
        ops.uniaxialMaterial('Fatigue', 3, 1, '-E0', 0.4, '-m', -0.458, '-min', -1e16, '-max', 1e16)
        #INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/Fatigue.html
        ops.element('zeroLength', 1, 1, 2,'-mat', 1, 2, 3, '-dir', 1, 1, 1)


    # Load pattern
    ops.timeSeries("Path", 1, "-dt", dtF, "-values", *F_input.tolist())
    ops.pattern("Plain", 1, 1)
    ops.load(2, 1.0)

    # Analysis setup
    ops.integrator("Newmark", 0.5, 0.25)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.algorithm("Newton")
    ops.analysis("Transient")

    time, disp, vel, accel, reaction = [], [], [], [], []  
    stiffness, OMEGA, PERIOD = [], [], []
    
    stable = 0
    current_time = 0.0
            
    while stable == 0 and current_time < duration:
        ops.analyze(1, dt)
        S01.ANALYSIS(stable, 1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS
        current_time = ops.getTime()
        time.append(current_time)
        disp.append(ops.nodeDisp(2, 1))
        vel.append(ops.nodeVel(2, 1))
        accel.append(ops.nodeAccel(2, 1))
        ops.reactions()
        reaction.append(ops.nodeReaction(1, 1)) # BASE REACTION
        stiffness.append(np.abs(reaction[-1]) / np.abs(disp[-1]))
        OMEGA.append(np.sqrt(stiffness[-1]/M))
        PERIOD.append((np.pi * 2) / OMEGA[-1])
        
    T_el, zeta_el = S02.ESTIMATE_T_ZETA(disp, dt)
    print(f"  Period T ≈ {T_el:.3f} s")
    print(f"  Damping ratio ζ ≈ {zeta_el:.4f}")    

    return np.array(time), np.array(disp), np.array(vel), np.array(accel), np.array(reaction), stiffness, PERIOD, zeta_el, T_el

#%%------------------------------------------------------------------------------
# FIND BEST SPRING AREA WITH STRUCTURAL BEHAVIOR COEFFICIENT OPTIMIZATION:
# Define parameters
M = 50000.0       # [kg] Mass
zi = 0.05         # Damping ratio
duration = 50.0   # [s] Analysis duration
dt = 0.01         # [s] Time step
area = 100        # [m²] Cross-sectional area    

X = area           # [m] Intial Guess Spring Area   
ESP = 1e-2         # Finite difference derivative Convergence Tolerance
TOLERANCE = 1e-6   # Convergence Tolerance
RESIDUAL = 100     # Convergence Residual 
IT = 0             # Intial Iteration
ITMAX = 100000     # Max. Iteration
TARGET_D = 0.04    # [m] Target Spring Displacement


SPRING_KIND = 'INELASTIC'
# Analysis Durations:
starttime = TI.process_time()

### FIND THE OPTIMUM VALUE 
while (RESIDUAL > TOLERANCE):
    # X -------------------------------------------------------
    DATA  = SDOF_FORCE_PULSE(X, SPRING_KIND, duration, dt)
    time,  disp, vel, accel, reaction, stiffness, PERIOD, damping_ratio, period = DATA
    SUPPLY = np.max(np.abs(disp))
    print(f' SUPPLY: {SUPPLY:.8f}')
    F = SUPPLY - TARGET_D
    print('F: ', F)
    # Xmin -------------------------------------------------------
    XMIN = X - ESP  
    DATA  = SDOF_FORCE_PULSE(XMIN, SPRING_KIND, duration, dt)
    time,  disp, vel, accel, reaction, stiffness, PERIOD, damping_ratio, period = DATA
    SUPPLYmin = np.max(np.abs(disp))
    Fmin = SUPPLYmin - TARGET_D
    print('Fmin: ', Fmin)
    # Xmax -------------------------------------------------------
    XMAX = X + ESP 
    DATA  = SDOF_FORCE_PULSE(XMAX, SPRING_KIND, duration, dt)
    time, disp, vel, accel, reaction, stiffness, PERIOD, damping_ratio, period = DATA
    SUPPLYmax = np.max(np.abs(disp))
    Fmax = SUPPLYmax - TARGET_D
    print('Fmax: ', Fmax)
    # DF -------------------------------------------------------
    DF = (Fmax - Fmin)/(2 * ESP);# Calculate the Finite difference derivative of F
    print('DF: ', DF)
    # DX -------------------------------------------------------
    DX = F / DF;        # Calculate dx
    print('DX: ', DX)
    # RESIDUAL -------------------------------------------------
    RESIDUAL = abs(DX); # Calculate residual
    X -= DX;            # update X
    IT += 1;            # update iteration
    print('IT: ', IT,' - RESIDUAL: ', RESIDUAL,' - X: ', X,'\n')
                
    if IT == ITMAX:
        print('\t\t Iteration reached to Max. Iteration')
        print('\t\t Change ESP and TOLERANCE for better Convergence')
        X = -1
        break;
    if RESIDUAL < TOLERANCE:
        print(f'\t\t Optimum Spring Area :                     {X:.6f}')
        print(f'\t\t Iteration Counts:                         {IT}')
        print(f'\t\t Convergence Residual:                     {RESIDUAL:.10e}')

    

totaltime = TI.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n')


#%% ---------------------------------------------
# Tension Force-Displacement Relationship (Positive values)
FY1, DY1 = 2.7720*X, 0.01    # First yield point
FY2, DY2 = 3.1046*X, 0.02    # Peak force
FY3, DY3 = 1.6632*X, 0.04    # Post-peak softening
FY4, DY4 = 1.6632*X, 0.06    # Plateau
FY5, DY5 = 0.2772*X, 0.28    # Further softening
FY6, DY6 = 0.2000*X, 0.51    # Near-zero stiffness
FY7, DY7 = 0.0*X, 0.72       # Zero force (failure)
    
KP = np.array([FY1, DY1, FY2, DY2, FY3, DY3, FY4, DY4, FY5, DY5, FY6, DY6, FY7, DY7])

# Compression Force-Displacement Relationship (Negative values)
FY1n, DY1n = -2.7720*X, -0.01    # First yield in compression
FY2n, DY2n = -3.1046*X, -0.05    # Peak compressive force
FY3n, DY3n = -1.6632*X, -0.20    # Post-peak softening

KN = np.array([FY1n, DY1n, FY2n, DY2n, FY3n, DY3n])

# Separate into Force and Displacement
force_p = KP[0::2]
disp_p = KP[1::2]

force_n = KN[0::2]
disp_n = KN[1::2]

#%%% Plotting Force-Displacement Diagram for Inelastic Spring
plt.figure(1, figsize=(8, 6))
plt.plot(disp_p, force_p, 'r-o', label='Tension')
plt.plot(disp_n, force_n, 'b-o', label='Compression')
plt.axhline(0, color='gray', linewidth=0.5)
plt.axvline(0, color='gray', linewidth=0.5)
plt.xlabel('Displacement')
plt.ylabel('Force')
plt.title('Force-Displacement Diagram for Inelastic Spring')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()


#%% --------------------------
# Plot Results
plt.figure(2, figsize=(18, 14))

# Define colors
elastic_color = 'black'
inelastic_color = 'red'

# Force pulse
plt.subplot(7, 1, 1)
plt.plot(timeF, F_input, color=inelastic_color, lw=2)
plt.ylabel("Force (N)")
plt.title("Applied Force Pulse", fontsize=12, fontweight="bold")

# Reaction plot
plt.subplot(7, 1, 2)
plt.plot(time, reaction, color=inelastic_color, linewidth=1.5, label=f'Damping Ratio: {100*damping_ratio:.3e} %')
plt.title('Reaction Forces vs Time', fontsize=12, pad=10)
plt.ylabel('Reaction (N)', fontsize=10)
plt.grid(alpha=0.3)
plt.legend(loc='upper right', framealpha=1)

# Displacement plot
plt.subplot(7, 1, 3)
plt.plot(time, disp, color=inelastic_color, linewidth=1.5)
plt.title('Displacement vs Time', fontsize=12, pad=10)
plt.ylabel('Displacement (m)', fontsize=10)
plt.grid(alpha=0.3)

# Velocity plot
plt.subplot(7, 1, 4)
plt.plot(time, vel, color=inelastic_color, linewidth=1.5)
plt.title('Velocity vs Time', fontsize=12, pad=10)
plt.ylabel('Velocity (m/s)', fontsize=10)
plt.grid(alpha=0.3)

# Acceleration plot
plt.subplot(7, 1, 5)
plt.plot(time, accel, color=inelastic_color, linewidth=1.5)
plt.title('Acceleration vs Time', fontsize=12, pad=10)
plt.ylabel('Acceleration (m/s²)', fontsize=10)
plt.xlabel('Time (s)', fontsize=10)
plt.grid(alpha=0.3)

# Stiffness plot
plt.subplot(7, 1, 6)
plt.plot(time, stiffness, color=inelastic_color, linewidth=1.5)
plt.title('Stiffness vs Time', fontsize=12, pad=10)
plt.ylabel('Stiffness (N/m)', fontsize=10)
plt.semilogy()
plt.grid(alpha=0.3)

# Period plot
plt.subplot(7, 1, 7)
plt.plot(time, PERIOD, color=inelastic_color, linewidth=1.5, label=f'Structure Period: {period:.3e}')
plt.title('Period vs Time', fontsize=12, pad=10)
plt.ylabel('Period  (s)', fontsize=10)
plt.xlabel('Time (s)', fontsize=10)
plt.legend(loc='upper right', framealpha=1)
plt.semilogy()
plt.grid(alpha=0.3)

plt.tight_layout()
plt.show()

#%%-----------------------------------
plt.figure(3, figsize=(8, 6))
plt.plot(disp, reaction, color='purple', linewidth=2)
plt.xlabel('Displacement [m]')
plt.ylabel('Base-reaction [N]')
plt.title(f'Displacement vs Base-reaction')
plt.grid()
plt.show()
#%% ---------------------------------------------
# Print out the state of nodes 1 and 2
ops.printModel("node",1, 2)
# Print out the state of element 1
ops.printModel("ele", 1)
# Print the Model
#printModel()
ops.printModel("-JSON", "-file", "SDOF_FORCE-PULSE_IMPACT_LOAD_FATIGUE_OPTIMIZATION_NEWTON-RAPHSON_METHOD.json")
#%%-------------------------------------------------------------------------------