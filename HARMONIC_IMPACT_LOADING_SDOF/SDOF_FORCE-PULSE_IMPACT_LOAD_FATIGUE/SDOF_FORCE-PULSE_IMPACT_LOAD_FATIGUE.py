######################################################################################################################
#                          >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                          #
#       FORCE-PULSE IMPACT LOAD ANALYSIS OF A SINGLE-DEGREE-OF-FREEDOM (SDOF) SYSTEM USING OPENSEES                  #
#--------------------------------------------------------------------------------------------------------------------#
#                                                  WITH FATIGUE MATERIAL                                             #
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
import matplotlib.pyplot as plt
import ANALYSIS_FUNCTION as S01
import ESTIMATE_T_ZETA_FUN as S02

# --------------------------
# SDOF parameters
M = 1000.0     # [kg] Mass
zi = 0.02      # Damping ratio  
alpha = 1.0    # velocity exponent (usually 0.3–1.0)

duration = 20.0   # [s] Analysis duration
dt = 0.01         # [s] Time step
time = np.arange(0, duration, dt)

# Input force (simple pulse)
F_input = np.zeros_like(time)
F_input[10:150] = 3000.0  # [N] short force pulse

#%% Inelastic Force-Displacement Parameters for Inelastic Spring
# Tension Force-Displacement Relationship (Positive values)
FY1, DY1 = 2772.0, 0.01    # First yield point
FY2, DY2 = 3104.6, 0.02     # Peak force
FY3, DY3 = 1663.2, 0.04     # Post-peak softening
FY4, DY4 = 1663.2, 0.06     # Plateau
FY5, DY5 = 277.2, 0.28      # Further softening
FY6, DY6 = 200.0, 0.41      # Near-zero stiffness
FY7, DY7 = 0.0, 0.52        # Zero force (failure)

KP = np.array([FY1, DY1, FY2, DY2, FY3, DY3, FY4, DY4, FY5, DY5, FY6, DY6, FY7, DY7])

# Compression Force-Displacement Relationship (Negative values)
FY1n, DY1n = -2772.0, -0.01    # First yield in compression
FY2n, DY2n = -3104.6, -0.20    # Peak compressive force
FY3n, DY3n = -1663.2, -0.40    # Post-peak softening

KN = np.array([FY1n, DY1n, FY2n, DY2n, FY3n, DY3n])

# Separate into Force and Displacement
force_p = KP[0::2]
disp_p = KP[1::2]

force_n = KN[0::2]
disp_n = KN[1::2]

#%% Plotting Force-Displacement Diagram for Inelastic Spring
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
plt.show()

#%%
K = KP[0] / KP[1]        # [N/m] Elastic Stiffness
omega = np.sqrt(K/M)
Cd = 2 * zi * omega * M  # [N·s/m] Damping coefficient 

#%% DEFINE ANALYSIS PROPERTIES
MAX_ITERATIONS = 20000    # Convergence iteration for test
MAX_TOLERANCE = 1.0e-10   # Convergence tolerance for test
#SPRING_KIND: 1 -> 'ELASTIC'
#SPRING_KIND: 2 -> 'INELASTIC'
def SDOF_FORCE_PULSE(SPRING_KIND, KP, KN):
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
    ops.timeSeries("Path", 1, "-dt", dt, "-values", *F_input.tolist())
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

    return np.array(time), np.array(disp), np.array(vel), np.array(accel), np.array(reaction)

# Run both cases
time , disp_elastic, vel_elastic, acc_elastic, rec_elastic = SDOF_FORCE_PULSE('ELASTIC', KP, KN)
time , disp_inelastic, vel_inelastic, acc_inelastic, rec_inelastic  = SDOF_FORCE_PULSE('INELASTIC', KP, KN)

#%% --------------------------
# Estimate Structure Period and Damping Ratio
T_el, zeta_el = S02.ESTIMATE_T_ZETA(disp_elastic, dt)
T_bi, zeta_bi = S02.ESTIMATE_T_ZETA(disp_inelastic, dt)

#%% --------------------------
# Display results
print("Elastic model:")
print(f"  Period T ≈ {T_el:.3f} s")
print(f"  Damping ratio ζ ≈ {zeta_el:.4f}")
print("Inelastic model:")
print(f"  Period T ≈ {T_bi:.3f} s")
print(f"  Damping ratio ζ ≈ {zeta_bi:.4f}")

#%% --------------------------
# Plot results
plt.style.use("seaborn-v0_8-darkgrid")
fig, axs = plt.subplots(5, 1, figsize=(20, 15), sharex=True)

colors = {
    "elastic": "black",
    "inelastic": "red",  
    "force": "purple"      
}

# Force pulse
axs[0].plot(time, F_input, color=colors["force"], lw=2)
axs[0].set_ylabel("Force (N)")
axs[0].set_title("Applied Force Pulse", fontsize=12, fontweight="bold")

# Displacement
axs[1].plot(time, disp_elastic, color=colors["elastic"], lw=2, label=f"Elastic (T≈{T_el:.2f}s, ζ≈{zeta_el:.3f})")
axs[1].plot(time, disp_inelastic, color=colors["inelastic"], lw=2, label=f"Inelastic (T≈{T_bi:.2f}s, ζ≈{zeta_bi:.3f})")
axs[1].set_ylabel("Disp (m)")
axs[1].legend(framealpha=0.9)
axs[1].set_title("Displacement Response", fontsize=12, fontweight="bold")

# Velocity
axs[2].plot(time, vel_elastic, color=colors["elastic"], lw=2, label="Elastic")
axs[2].plot(time, vel_inelastic, color=colors["inelastic"], lw=2, label="Inelastic")
axs[2].set_ylabel("Vel (m/s)")
axs[2].legend(framealpha=0.9)
axs[2].set_title("Velocity Response", fontsize=12, fontweight="bold")

# Acceleration
axs[3].plot(time, acc_elastic, color=colors["elastic"], lw=2, label="Elastic")
axs[3].plot(time, acc_inelastic, color=colors["inelastic"], lw=2, label="Inelastic")
axs[3].set_ylabel("Accel (m/s²)")
axs[3].set_xlabel("Time (s)")
axs[3].legend(framealpha=0.9)
axs[3].set_title("Acceleration Response", fontsize=12, fontweight="bold")

# Reaction
axs[4].plot(time, rec_elastic, color=colors["elastic"], lw=2, label="Elastic")
axs[4].plot(time, rec_inelastic, color=colors["inelastic"], lw=2, label="Inelastic")
axs[4].set_ylabel("Reaction (N)")
axs[4].set_xlabel("Time (s)")
axs[4].legend(framealpha=0.9)
axs[4].set_title("Reaction Response", fontsize=12, fontweight="bold")

plt.tight_layout()
plt.show()


plt.figure(3, figsize=(8, 6))
plt.plot(disp_elastic, rec_elastic, color=colors["elastic"], linewidth=2)
plt.plot(disp_inelastic, rec_inelastic, color=colors["inelastic"], linewidth=2)
plt.xlabel('Displacement [m]')
plt.ylabel('Base-reaction [N]')
plt.title('Displacement vs Base-reaction')
plt.legend(['ELASTIC', 'INELASTIC'])
plt.grid()
plt.show()