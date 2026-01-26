###########################################################################################################
#                    >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                     #
#           FRAGILITY ANALYSIS BASED ON ACCELERATION AND STRUCTURAL DUCTILITY DAMAGE INDEX WITH           #
#              INCREMENTAL DYNAMIC ANALYSIS (IDA) OF A SINGLE-DEGREE-OF-FREEDOM (SDOF) SYSTEM             #
#                                UTILIZING 100 GROUND MOTIONS IN OPENSEES                                 #
#---------------------------------------------------------------------------------------------------------#
#                                          PARALLEL COMPUTING VERSION                                     #
#---------------------------------------------------------------------------------------------------------#
# This program performs Incremental Dynamic Analysis (IDA) on a Single-Degree-of-Freedom (SDOF) system    #
# subjected to 100 seismic ground motions. The analysis evaluates the structural response under varying   #
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
# PARALLEL PROCESSING MEANS RUNNING SEVERAL TASKS AT THE SAME TIME INSTEAD OF ONE AFTER ANOTHER.          #
# IN THE CODE, EACH STEP ANALYSIS WAS CALCULATED IN SEQUENCE,                                             #
# SO THE CPU WORKED ON ONLY ONE MODE AT ANY MOMENT. IN THE REWRITTEN VERSION, THE JOBLIB LIBRARY ALLOWS   #
# ALL FOUR MODES TO RUN SIMULTANEOUSLY ON DIFFERENT CPU CORES. EACH CORE PROCESSES ONE MODE INDEPENDENTLY,#
# SO THE TOTAL COMPUTATION TIME BECOMES MUCH SHORTER.                                                     #
#                                                                                                         #
# MODERN COMPUTERS USUALLY HAVE MULTIPLE CORES, FOR EXAMPLE 4, 8, OR EVEN MORE. WHEN WE USE PARALLEL      #
# PROCESSING, WE DIVIDE THE WORKLOAD ACROSS THESE CORES. BECAUSE EACH MODE IS A SEPARATE AND INDEPENDENT  #
# ANALYSIS, THEY ARE PERFECT FOR PARALLEL EXECUTION. INSTEAD OF WAITING FOR MODE 1 TO FINISH BEFORE       #
# STARTING MODE 2, ALL MODES START TOGETHER AND FINISH ALMOST TOGETHER.                                   #
#                                                                                                         #
# IN PRACTICE, THE SPEED IMPROVEMENT DEPENDS ON HOW MANY CORES YOUR CPU HAS. IF YOUR COMPUTER HAS 4 CORES,#
# THE RUNTIME CAN BE UP TO FOUR TIMES FASTER. IN MANY CASES THE SPEEDUP IS AROUND 3–4 TIMES,              #
# BECAUSE THERE IS A SMALL OVERHEAD WHEN STARTING PARALLEL TASKS. THE REWRITTEN CODE USES PARALLEL        #
# AND DELAYED TO AUTOMATICALLY SEND EACH MODE TO A DIFFERENT CORE AND THEN COLLECT ALL RESULTS            #
# IN THE CORRECT ORDER. THIS MAKES THE ANALYSIS MORE EFFICIENT WITHOUT CHANGING THE ENGINEERING RESULTS.  #
#                                                                                                         #
# PARALLEL PROCESSING IS ESPECIALLY HELPFUL IN STRUCTURAL ENGINEERING SIMULATIONS WHERE EACH ANALYSIS     #
# REQUIRES HEAVY NUMERICAL CALCULATION, SUCH AS NONLINEAR POST-BUCKLING. BY USING ALL AVAILABLE CPU POWER,#
# YOU FINISH THE WORK FASTER AND CAN TEST MORE CASES OR MORE MODELS IN THE SAME AMOUNT OF TIME.           #
#---------------------------------------------------------------------------------------------------------#
#                          THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                     #
#                                   EMAIL: salar.d.ghashghaei@gmail.com                                   #
########################################################################################################### 

"""
This code implements a comprehensive nonlinear incremental dynamic analysis framework for
performance-based earthquake engineering assessment of single-degree-of-freedom
(SDOF) systems. The methodology combines traditional nonlinear time-history
analysis with modern probabilistic and machine learning techniques for advanced
structural performance evaluation.

KEY ENGINEERING OBJECTIVES:
1. Comparative assessment of hysteretic models for seismic response prediction
2. Probabilistic seismic demand analysis using multiple ground motions
3. Development of fragility curves for performance-based earthquake engineering
4. Integration of data science methods for structural reliability assessment

ANALYTICAL FEATURES:
- Nonlinear material behavior with pinching and degradation
- Response spectrum analysis across period range
- Real-time structural health monitoring metrics
- Statistical characterization of seismic demands
- Machine learning-based damage prediction
---------------------------------
Model setup:
 - SDOF properties: mass (m), initial stiffness (k), yield displacement (Dy), ultimate displacement (Du), viscous damping (xi).
 - Hysteresis models: HYSTERETICSM (pinching, stiffness degradation, strength decay).
 - Damping: Rayleigh (or equivalent viscous) damping specified by target damping ratio xi for the fundamental mode.

Dynamic response:
 - Natural period T = 2*pi*sqrt(m/k) computed from linearized stiffness.
 - Time-history integration produces displacement, velocity, acceleration and base reaction histories.
 - HYSTERETIC model shows faster decay of amplitude and larger energy dissipation due to pinching and degradation.

Force–displacement behavior:
 - BILINEAR: symmetric hysteresis loops with stable post-yield stiffness; residual displacements are primarily due to plastic offset.
 - HYSTERETIC: pinched loops, reduced unloading/reloading stiffness, strength decay and larger residuals; captures cumulative damage effects.

Stiffness and strength evolution:
 - Effective lateral stiffness reduces during the excitation for both models but degrades faster with HYSTERETIC due to damage mechanisms.
 - Strength deterioration (reduced peak restoring force) in HYSTERETIC leads to reduced re-centering and larger residuals.

Damping estimation:
 - Use logarithmic decrement or energy-based measures from free vibration or post-event cycles.
 - HYSTERETIC typically yields higher equivalent damping (greater energy dissipation) compared with BILINEAR for the same displacement amplitude.

Peak responses:
 - Peak displacement: often lower for HYSTERETIC in early cycles because of softening, but long-term residual displacement may be higher.
 - Peak base shear (reaction): decays faster in HYSTERETIC due to strength loss; BILINEAR sustains higher peak restoring forces for the same drift until hardening or limiting criteria apply.

Visualization:
 - Plot time histories (disp, vel, acc), hysteresis loops (force vs disp), and envelope curves to compare models.
 - Response spectra for displacement, velocity, and acceleration can be constructed from peak responses across parameter sweeps (e.g., varying T or post-yield stiffness).

Implications for seismic assessment:
 - BILINEAR: simple and computationally efficient; may overestimate resilience for severe cyclic demands because it omits degradation.
 - HYSTERETIC: captures important degradation mechanisms (pinching, stiffness/strength loss, ultimate strain) and is recommended for collapse assessment and detailed damage estimation.
 - Model selection should match the performance objective: serviceability checks might use simpler models; collapse and damage-sensitive studies require degraded hysteretic models calibrated with experiment.

Data export and post-processing:
 - Store peak and time-history results (displacement, velocity, acceleration, base reaction) to CSV/Excel for parametric studies.
 - Compute and plot response spectra (disp/vel/acc/reaction) from the stored peak values.

Ductility Damage Index (DDI) — implementation (concept):
 - After identifying yield displacement Dy and ultimate displacement Du from the capacity model:
   Dd = max(|disp_time_history|)  # maximum absolute dynamic displacement demand
   DI = (Dd - Dy) / (Du - Dy)      # Ductility Damage Index in the direction of interest
   Interpretation:
     DI <= 0   : elastic (no damage)
     0 < DI < 1: inelastic damage (serviceability/repairable)
     DI >= 1   : demand reaches or exceeds ultimate capacity (collapse or unacceptable damage)
"""
#%%------------------------------------------------------------------------------------------------
import openseespy.opensees as ops
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time as TI
import SALAR_MATH as S01
import ANALYSIS_FUNCTION as S02
import MARKOV_CHAIN as S03
import RAYLEIGH_DAMPING_FUN as S04
import EIGENVALUE_ANALYSIS_FUN as S05
import DAMPING_RATIO_FUN as S06
import FRAGILITY_CURVE_FUN as S07
from scipy.stats import norm
import pandas as pd
from joblib import Parallel, delayed
#%%------------------------------------------------------------------------------------------------
# Define parameters (units: m, N)
NUM_G = 50                                      # Total number for Standard Acceleration of Gravity in each simulation
NUM_SEISMIC = 100                               # Total number for seismic simulation
#%%------------------------------------------------------------------------------------------------
# Define  Structural Properties
FY = 85000.0                                     # [N] Yield Force of Structure
FU = 1.5 * FY                                    # [N] Ultimate Force of Structure
Ke = 4500000.0                                   # [N/m] Spring Elastic Stiffness
DY = FY / Ke                                     # [m] Yield Displacement
DSU = 0.36                                       # [m] Ultimate Displacement
Ksh = (FU - FY) / (DSU - DY)                     # [N/m] Displacement Hardening Modulus
Kp = FU / DSU                                    # [N/m] Spring Plastic Stiffness
b = Ksh / Ke                                     # Displacement Hardening Ratio
M = 15000.0                                      # [kg] Mass of the Structure
DR = 0.03                                        # Damping Ratio

duration = 15.0  # [s] Total simulation duration
dt = 0.01        # [s] Time step
GMfact = 9.81    # [m/s^2] Standard Acceleration of Gravity
#%%------------------------------------------------------------------------------------------------
# Positive branch points
pos_disp = [0, DY, DSU, 1.1*DSU, 1.25*DSU]
pos_force = [0, FY, FU, 0.2*FU, 0.1*FU]
KP = np.array([FY, DY, FU, DSU, 0.2*FU, 1.1*DSU, 0.1*FU, 1.25*DSU])

# Negative branch points
neg_disp = [0, -DY, -DSU, -1.1*DSU, -1.25*DSU]
neg_force = [0, -FY, -FU, -0.2*FU, -0.1*FU]
KN = np.array([-FY, -DY, -FU, -DSU, -0.2*FU, -1.1*DSU, -0.1*FU, -1.25*DSU])

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

# INITIAL MASS FOR INCREMENTAL RESPONSE SPECTRUM ANALYSIS
mi = 0.2 * (PLAS_PERIOD/2*np.pi)**2 * Kp # Mass of Structure is 20 Percent of Effective Plastic Period
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
#%%------------------------------------------------------------------------------------------------
# Define Analysis Properties
MAX_ITERATIONS = 1000000   # Convergence iteration for test
MAX_TOLERANCE = 1.0e-6    # Convergence tolerance for test
#%%------------------------------------------------------------------------------------------------
### OPENSEES FUNCTION
def ANALYSIS_SDOF(i, m, g):
    # Initialize OpenSees model
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)
    
    GMfact = g  # [m/s^2] standard acceleration of gravity or standard acceleration  
    # Define nodes
    ops.node(1, 0.0)  # Fixed base
    ops.node(2, 0.0)  # Mass node
        
    # Define boundary conditions
    ops.fix(1, 1)

    # Define mass
    ops.mass(2, m)
        
    # Define material properties
    MatTag = 1
    """
    pinchX = 0.8           # Pinching factor in X direction
    pinchY = 0.5           # Pinching factor in Y direction
    damage1 = 0.0          # Damage due to ductility
    damage2 = 0.0          # Damage due to energy
    beta = 0.1             # Stiffness degradation parameter
    ops.uniaxialMaterial('Hysteretic', MatTag, FY, DY, FU, DSU, 0.2*FU, 1.1*DSU, -FY, -DY, -FU, -DSU, -0.2*FU, -1.1*DSU, pinchX, pinchY, damage1, damage2, beta)
    # INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Hysteretic_Material
    """
    ops.uniaxialMaterial('HystereticSM', MatTag, '-posEnv', *KP.flatten(), '-negEnv', *KN.flatten(), '-pinch', 1, 1)
    #INFO LINK: https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/HystereticSM.html
    alpha = 1.0    # velocity exponent (usually 0.3–1.0)
    omega = np.sqrt(Ke/m)
    Cd = 2 * DR * omega * m  # [N·s/m] Damping coefficient 
    ops.uniaxialMaterial('Viscous', 2, Cd, alpha)  # Material for C (alpha=1.0 for linear)
    
    # Define element
    ops.element('zeroLength', 1, 1, 2, '-mat', MatTag, 2, '-dir', 1, 1)  # DOF[1] LATERAL SPRING
    
    # Apply seismic accelerations    
    # Define time series for input motion (Acceleration time history)
    gm_accels = np.loadtxt(f'Ground_Acceleration_{i+1}.txt')  # Assumes acceleration in m/s²
    ops.timeSeries('Path', 1, '-dt', dt, '-values', *gm_accels.tolist(), '-factor', GMfact) # SEISMIC-X
    #ops.timeSeries('Path', 1, '-dt', dt, '-filePath', f'Ground_Acceleration{i+1}.txt', '-factor', GMfact) # SEISMIC-X
        
    # Define load patterns
    # pattern UniformExcitation $patternTag $dof -accel $tsTag <-vel0 $vel0> <-fact $cFact>
    SSF = 10.00 # Seismic Scale Factor
    ops.pattern('UniformExcitation', 200, 1, '-accel', 1, '-fact', SSF) # SEISMIC-X
    
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
    a0 = (2 * Omega01 * DR) / Omega01 # c = a0 * m : Mass-proportional damping
    a1 = (DR * 2) / Omega01 # c = a1 * k : Stiffness-proportional damping
    # Apply Rayleigh damping
    ops.rayleigh(a0, 0.0, 0.0, a1)# INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/reyleigh.html
    #ops.rayleigh(0, 0, 2 * DR * Omega01, 0) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/reyleigh.html
    PERIOD = 2*np.pi / Omega01   # Structure Period 
    print(f'PERIOD: {PERIOD:.6f}')
    
    # Calculate Rayleigh damping factors
    #PERIOD_01, PERIOD_02 = S04.RAYLEIGH_DAMPING(1, 0.5*DR, DR, 0, 1)
    
    # Run Dynamic Analysis
    time = []
    displacement = []
    velocity = []
    acceleration = []
    base_reaction = []
    DI = []
    stiffness = []
    PERIOD_MIN, PERIOD_MAX = [], []

        
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
            acceleration.append(ops.nodeAccel(2, 1))                 # Structure Acceleration
        base_reaction.append(-ops.eleResponse(1, 'force')[0])        # Reaction force
        DI.append(100*(np.abs(displacement[-1]) - DY) / (DSU - DY))  # [%] Structural Ductility Damage Index 
        stiffness.append(np.abs(base_reaction[-1]/displacement[-1])) # [N/m] Stiffness
        # IN EACH STEP STRUCTURAL PERIOD WILL BE CALCULATED
        #PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(1, PLOT=True)
        #PERIOD_MIN.append(PERIODmin)
        #PERIOD_MAX.append(PERIODmax)
        step += 1
        #print(f'{time[-1]} {displacement[-1]:0.4e} {base_reaction[-1]:0.4e}')
        
    # Calculate Structure Damping Ratio Based on Lateral Displacement
    damping_ratio = S06.DAMPING_RATIO(displacement) 
    #ops.wipe()
    return time, displacement, velocity, acceleration, base_reaction, DI, PERIOD, damping_ratio, stiffness

#%%------------------------------------------------------------------------------------------------
# SINGLE PARALLEL ANALYSIS FUNCTION
def run_single_analysis(j, i, mi, GMfact, NUM_G):
    """
    One independent nonlinear SDOF analysis (Safe for parallel computing)
    """

    m  = mi
    gi = 2.0 * GMfact * ((i + 1) / NUM_G)

    data = ANALYSIS_SDOF(j, m, gi)

    (time, displacement, velocity, acceleration,
     base_reaction, DI, PERIOD, damping_ratio, stiffness) = data

    return {
        "disp":  np.max(np.abs(displacement)),
        "vel":   np.max(np.abs(velocity)),
        "acc":   np.max(np.abs(acceleration)),
        "react": np.max(np.abs(base_reaction)),
        "DI":    np.clip(np.max(np.abs(DI)), 0.0, 100.0),
        "K":     np.max(np.abs(stiffness)),
        "G":     gi,
        "damp":  damping_ratio,
        "T":     PERIOD
    }

#%%------------------------------------------------------------------------------------------------
# ----------------------  PARALLEL PROCESSING  ----------------------
# MAIN PARALLEL ANALYSIS LOOP

# Analysis Durations:
current_time = TI.strftime("%H:%M:%S", TI.localtime())
print("Start Time:", current_time)

DATA = {
    1: [],  # DISPLACEMENT
    2: [],  # VELOCITY
    3: [],  # ACCELERATION
    4: [],  # BASE REACTION
    5: [],  # DAMAGE INDEX
    6: [],  # DAMPING RATIO
    7: [],  # STIFFNESS
}

for j in range(NUM_SEISMIC):

    print(f"\n--- START SEISMIC {j + 1} ---")

    results = Parallel(
        n_jobs=-1,            # Use all CPU cores
        backend="loky",       # Process-based (OpenSees safe)
        verbose=5
    )(
        delayed(run_single_analysis)(j, i, mi, GMfact, NUM_G)
        for i in range(NUM_G)
    )


    # COLLECT RESULTS
    max_displacement = [r["disp"]  for r in results]
    max_velocity     = [r["vel"]   for r in results]
    max_acceleration = [r["acc"]   for r in results]
    max_base_reaction= [r["react"] for r in results]
    max_DI           = [r["DI"]    for r in results]
    max_K            = [r["K"]     for r in results]
    max_G            = [r["G"]     for r in results]
    DAMPING_RATIO    = [r["damp"]  for r in results]

    DATA[1].append(max_displacement)
    DATA[2].append(max_velocity)
    DATA[3].append(max_acceleration)
    DATA[4].append(max_base_reaction)
    DATA[5].append(max_DI)
    DATA[6].append(DAMPING_RATIO)
    DATA[7].append(max_K)

    print(f"--- SEISMIC {j + 1} DONE ---")

print("\nAnalysis completed successfully")

current_time = TI.strftime("%H:%M:%S", TI.localtime())
print("Finish Time:", current_time)

#%%------------------------------------------------------------------------------------------------
# EXPORT RESULTS TO EXCEL

labels = {
    1: "DISPLACEMENT",
    2: "VELOCITY",
    3: "ACCELERATION",
    4: "BASE_REACTION",
    5: "DAMAGE_INDEX",
    6: "DAMPING_RATIO",
    7: "STIFFNESS",
}

with pd.ExcelWriter(
    "INELASTIC_SDOF_INCREMENTAL_DYNAMIC_ANALYSIS_SEISMIC_PARALLEL-COMPUTING_RESULTS.xlsx",
    engine="openpyxl"
) as writer:

    for key, value in DATA.items():
        df = pd.DataFrame(value)
        sheet_name = labels[key][:31]
        df.to_excel(writer, sheet_name=sheet_name, index=False)

print("Excel file saved successfully")

#%%------------------------------------------------------------------------------------------------
# PLOTTING FUNCTION

G = np.array(max_G)

def PLOT_2D(COUNT, XLABEL, YLABEL, TITLE):

    plt.figure(COUNT, figsize=(12, 10))

    # All simulations
    for j in range(NUM_SEISMIC):
        plt.plot(G, DATA[COUNT][j], alpha=0.35)

    arr = np.array(DATA[COUNT])

    mean_curve   = np.mean(arr, axis=0)
    median_curve = np.median(arr, axis=0)
    std_curve    = np.std(arr, axis=0)
    q1_curve     = np.percentile(arr, 25, axis=0)
    q3_curve     = np.percentile(arr, 75, axis=0)

    plt.plot(G, mean_curve,   'b-.', linewidth=3, label='Mean')
    plt.plot(G, median_curve, 'r--', linewidth=2, label='Median')
    plt.plot(G, q1_curve,     'g-.', linewidth=2, label='Q1 (25%)')
    plt.plot(G, q3_curve,     'orange', linestyle='-.', linewidth=2, label='Q3 (75%)')

    plt.fill_between(
        G,
        mean_curve - std_curve,
        mean_curve + std_curve,
        color='gray',
        alpha=0.25,
        label='Mean ± Std'
    )

    plt.xlabel(XLABEL)
    plt.ylabel(YLABEL)
    plt.title(TITLE)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

#%%------------------------------------------------------------------------------------------------
# PLOTS

PLOT_2D(1, "Standard Acceleration of Gravity [m/s²]", "Max Displacement (PGD) [m]", "Displacement Response Spectrum")

PLOT_2D(2, "Standard Acceleration of Gravity [m/s²]", "Max Velocity (PGV) [m/s]", "Velocity Response Spectrum")

PLOT_2D(3, "Standard Acceleration of Gravity [m/s²]", "Max Acceleration (PGA) [m/s²]", "Acceleration Response Spectrum")

PLOT_2D(4, "Standard Acceleration of Gravity [m/s²]", "Max Base Reaction [N]", "Base Reaction Response Spectrum")

PLOT_2D(5, "Standard Acceleration of Gravity [m/s²]","Damage Index [%]", "Ductility Damage Index Spectrum")

PLOT_2D(6, "Standard Acceleration of Gravity [m/s²]", "Damping Ratio [%]", "Damping Ratio Spectrum")

PLOT_2D(7, "Standard Acceleration of Gravity [m/s²]", "Structural Stiffness [N/m]", "Structural Stiffness Spectrum")
    

#%%------------------------------------------------------------------------------------------------
