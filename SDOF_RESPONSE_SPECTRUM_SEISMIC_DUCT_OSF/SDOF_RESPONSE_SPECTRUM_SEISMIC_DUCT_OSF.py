###########################################################################################################
#                    >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                     #
#      NONLINEAR DYNAMIC ANALYSIS UNDER A SINGLE GROUND MOTION RECORD COMPUTATION AND VISUALIZATION       #
#      RESPONSE SPECTRA OF ACCELERAION, VELOCITY, DISPLACEMENT DUCTILITY DAMAGE INDEX USING OPENSEES      #
#---------------------------------------------------------------------------------------------------------#
#                          CONSTANT STRUCTURAL DUCTILITY RATIO RESPONSE SPECTRUM                          #
#---------------------------------------------------------------------------------------------------------#
#                          THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                     #
#                                   EMAIL: salar.d.ghashghaei@gmail.com                                   #
###########################################################################################################
"""
This code implements a comprehensive nonlinear dynamic analysis framework for
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
-----------------------------------
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
 - Peak base shear (reaction): decays faster in HYSTERETICSM due to strength loss; BILINEAR sustains higher peak restoring forces for the same drift until hardening or limiting criteria apply.

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

Conclusions:
 - For inelastic SDOF studies, including pinching and degradation in the hysteretic model can change predicted peak responses, residual displacements, and damage indices significantly.
 - Use HYSTERETIC-type models for damage-sensitive or collapse-prone scenarios, and calibrate degradation parameters using test data where possible.
--------------------------------
Constant Structural Ductility Ratio
In seismic engineering, the structural ductility ratio (often denoted as μ) refers to the ability of a structure to undergo inelastic (plastic) deformations without collapsing during an earthquake. It is defined as the ratio of the maximum displacement (Δ_max) to the yield displacement (Δ_yield), i.e., μ = Δ_max / Δ_yield. This measures how much a structure can "stretch" beyond its elastic limit while dissipating energy through yielding, which helps prevent brittle failure. Ductile structures (e.g., those with μ > 5–6) can survive strong ground motions by deforming significantly, reducing the need for overly stiff designs.
A "constant structural ductility ratio" typically relates to constant ductility response spectra. In these spectra, the ductility ratio μ is held constant across different structural periods (T), and the required strength (e.g., pseudoacceleration or force) is plotted against the period to achieve that specific ductility demand. This approach contrasts with elastic response spectra, where no inelastic behavior is assumed. Constant ductility spectra are used to derive ductility reduction factors (Rμ or Rd), which adjust the elastic response spectrum downward to account for energy dissipation through ductility, making designs more economical. They are particularly useful in performance-based seismic design to ensure ductility supply exceeds demand, often calibrated using methods like Newmark-Hall for different period ranges (short: T < 0.2 s; intermediate: 0.2–0.5 s; long: T > 1 s)
"""
#%%------------------------------------------------------------------------------------------------
import openseespy.opensees as ops
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time as ti
import SALAR_MATH as S01
import ANALYSIS_FUNCTION as S02
import MARKOV_CHAIN as S03
import RAYLEIGH_DAMPING_FUN as S04
import EIGENVALUE_ANALYSIS_FUN as S05
import DAMPING_RATIO_FUN as S06
import FRAGILITY_CURVE_FUN as S07
import OPENSEEES_HYSTERETICSM_FORCE_DISP_FUN as S08
from scipy.stats import norm

#%%------------------------------------------------------------------------------------------------
# Define parameters (units: m, N)
NUM_PERIOD = 50                                 # Total number for Period in each simulation
NUM_SIM = 30                                    # Total number for analysis
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
#%%------------------------------------------------------------------------------------------------
"""
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
"""
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
#%%------------------------------------------------------------------------------------------------
# Define Analysis Properties
MAX_ITERATIONS = 1000000   # Convergence iteration for test
MAX_TOLERANCE = 1.0e-6    # Convergence tolerance for test

#%%------------------------------------------------------------------------------------------------
### OPENSEES FUNCTION
def ANALYSIS_SDOF(i, m, duct, osf):
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
    ops.mass(2, m)
        
    # Define material properties
    MatTag_S = 1 # SPRING TAG
    # FORCE-DISPLACEMENT RELATIONSHIP OF LATERAL SPRING AND PLOT 
    DP = [0, 0, 0, 0]
    FP = [0, 0, 0, 0]
    DN = [0, 0, 0, 0]
    FN = [0, 0, 0, 0]
    DSU = DY * duct # IN EACH STEP IT WILL BE CHNAGED
    FU = FY * osf   # IN EACH STEP IT WILL BE CHNAGED
    #print(DSU,"------------" ,FU)
    DP[0], FP[0] = DY, FY
    DP[1], FP[1] = DSU, FU 
    DP[2], FP[2] = 1.1*DSU, 0.20*FU
    DP[3], FP[3] = 1.25*DSU, 0.10*FU
    DN[0], FN[0] = -DY, -FY 
    DN[1], FN[1] = -DSU, -FU
    DN[2], FN[2] = -1.1*DSU, -0.20*FU   
    DN[3], FN[3] = -1.25*DSU, -0.10*FU
    #print(DP, FP)
    #print(DN, FN)
    S08.OPENSEEES_HYSTERETICSM_FORCE_DISP_FUN(MatTag_S, DP, FP, DN, FN, PLOT = False, X_LABEL='Displacement (mm)', Y_LABEL='Force [N]', TITLE='FORCE-DISPLACEMENT CURVE')
    MatTag_C = 2 # SPRING DAMPER
    alpha = 1.0    # velocity exponent (usually 0.3–1.0)
    omega = np.sqrt(Ke/m)
    Cd = 2 * DR * omega * m  # [N·s/m] Damping coefficient 
    ops.uniaxialMaterial('Viscous', 2, Cd, alpha)  # Material for C (alpha=1.0 for linear)
    
    # Define element
    ops.element('zeroLength', 1, 1, 2, '-mat',  MatTag_S, MatTag_C, '-dir', 1, 1)  # DOF[1] LATERAL SPRING
    
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
        DI.append(100*(displacement[-1] - DY) / (DSU - DY))          # [%] Structural Ductility Damage Index 
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
# Analysis Durations:
starttime = ti.process_time()

# Collect Data
DATA = {
    1: [],  # DISP
    2: [],  # VELOCITY
    3: [],  # ACCELERAION
    4: [],  # REACTION
    5: [],  # DI
    6: [],  # DAMPING RATIO
    7: [],  # STRUCTURAL STIFFNESS
    }

# NUM_SIM is the number of simulations
for j in range(NUM_SIM):
    # Initialize lists to store max values
    max_time = []
    max_displacement = []
    max_velocity = []
    max_acceleration = []
    max_base_reaction = []
    max_DI = []
    max_T = []
    max_K = []
    DAMPING_RATIO = []
    duct = np.linspace(2.0, 8.5, NUM_PERIOD).tolist()  # DEFINE MIN. AND MAX. VALUE FOR STRUCURAL DUCTILITY RATIO
    osf = np.linspace(1.05, 1.35, NUM_PERIOD).tolist() # DEFINE MIN. AND MAX. VALUE FOR STRUCURAL OVER STRENGTH FACTOR
    print('Structure Ductility Ratio: ', duct[j], " ------ ",'Sructure Over Strength Factor: ',osf[j])
    for i in range(NUM_PERIOD):
        m = mi * (i+1 / NUM_PERIOD) * 0.04
        data = ANALYSIS_SDOF(0, m, duct[j], osf[j])
        (time, displacement, velocity, acceleration,
        base_reaction, DI, PERIOD, damping_ratio, stiffness) = data
        # Calculate and store the max absolute values
        max_time.append(np.max(np.abs(time)))
        max_displacement.append(np.max(np.abs(displacement)))
        max_velocity.append(np.max(np.abs(velocity)))
        max_acceleration.append(np.max(np.abs(acceleration)))
        max_base_reaction.append(np.max(np.abs(base_reaction)))
        max_DI.append(np.max(np.abs(DI)))
        max_K.append(np.max(np.abs(stiffness)))
        DAMPING_RATIO.append(damping_ratio)

        # CHECK TO DAMAGE INDEX BOUNDARIES
        if max_DI[-1] <= 0.0:
            max_DI[-1] = 0.0
        if max_DI[-1] >= 100.0:
            max_DI[-1] = 100.0    
           
        #max_DI.append(DI[-1])
        max_T.append(PERIOD)
        print(f'STEP {i + 1} DONE')
    # Store Data
    DATA[1].append(max_displacement)
    DATA[2].append(max_velocity)
    DATA[3].append(max_acceleration)
    DATA[4].append(max_base_reaction)
    DATA[5].append(max_DI)   
    DATA[6].append(DAMPING_RATIO)
    DATA[7].append(max_K)
    print(f'\n SIMULATION {j + 1} DONE \n')    
else:
    print('\n\n Analysis completed successfully')    

totaltime = ti.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n') 

#%%------------------------------------------------------------------------------------------------
# OUTPUT DATA FROM RESPONSE SPECTRUM ANALYSIS TO EXCEL FILE
labels = {
    1: "DISPLACEMENT",
    2: "VELOCITY",
    3: "ACCELERATION",
    4: "BASE_REACTION",
    5: "DAMAGE_INDEX",
    6: "DAMPING_RATIO",
    7: "STIFFNESS",
}

with pd.ExcelWriter("SDOF_RESPONSE_SPECTRUM_SEISMIC_DUCT_OSF_RESULTS.xlsx", engine="openpyxl") as writer:
    for key, value in DATA.items():
        df = pd.DataFrame(value)
        sheet_name = labels[key][:31] # SHEEET NAME
        df.to_excel(writer, sheet_name=sheet_name, index=False)

print("Excel file saved: SDOF_RESULTS.xlsx")

#%%------------------------------------------------------------------------------------------------
# PLOT THE RESPONSE SPECTRUMS 
import matplotlib.pyplot as plt
import numpy as np

# Convert to numpy arrays for safety
T = np.array(max_T)

""""
def PLOT_2D(COUNT, X, Y, XLABEL, YLABEL, TITLE):
    plt.figure(COUNT, figsize=(12,10))
    for j in range(NUM_SIM):
        plt.plot(T, DATA[COUNT][j],label = f'Max: {np.max(np.abs(DATA[COUNT][j])):.4e}')
    plt.xlabel(XLABEL)
    plt.ylabel(YLABEL)
    plt.title(TITLE)
    plt.legend()
    plt.grid(True)
    plt.show()
"""    
def PLOT_2D(COUNT, XLABEL, YLABEL, TITLE, SEMILOGY=False):
    plt.figure(COUNT, figsize=(12, 10))

    # Plot all simulations
    for j in range(NUM_SIM):
        #plt.plot(T, DATA[COUNT][j], alpha=0.4, label=f'Sim {j+1} | Max: {np.max(np.abs(DATA[COUNT][j])):.4e}')
        plt.plot(T, DATA[COUNT][j], alpha=0.4)
    
    # Convert to NumPy array for statistical calculations
    arr = np.array(DATA[COUNT])    # Shape: (NUM_SIM, len(T))

    # Compute statistical metrics
    mean_curve   = np.mean(arr, axis=0)
    median_curve = np.median(arr, axis=0)
    std_curve    = np.std(arr, axis=0)
    std_curveDOWN = mean_curve - np.std(arr, axis=0)
    std_curveUP = mean_curve + np.std(arr, axis=0)
    q1_curve = np.percentile(arr, 25, axis=0)
    q3_curve = np.percentile(arr, 75, axis=0)

    # Plot statistical curves
    plt.plot(T, mean_curve, color='navy', linestyle='-.', linewidth=3, label='Mean')
    
    #plt.plot(T, std_curveDOWN, color='steelblue', linestyle='-', linewidth=2, label='Mean - Std')
    #plt.plot(T, std_curveUP,   color='deepskyblue', linestyle='-', linewidth=2, label='Mean + Std')
    
    plt.plot(T, median_curve, color='red', linestyle='--', linewidth=2, label='Median')
    
    plt.plot(T, q1_curve, color='green', linestyle='-.', linewidth=2, label='Q1 (25%)')
    plt.plot(T, q3_curve, color='orange', linestyle='-.', linewidth=2, label='Q3 (75%)')

    # Standard deviation band around the mean
    plt.fill_between(T, mean_curve - std_curve, mean_curve + std_curve,
                     color='gray', alpha=0.25, label='Mean ± Std')

    # Labels and title
    plt.xlabel(XLABEL)
    plt.ylabel(YLABEL)
    plt.title(TITLE)
    
    if SEMILOGY == True:
        plt.semilogy()
        
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()
    return mean_curve, median_curve

    
    
# ----------------------------------------
# Plot 1: Displacement Response Spectrum
# ----------------------------------------
X = T
Y = DATA
XLABEL = "Period T [s]"
YLABEL = "Max Displacement (PGD) [m]"
TITLE =  "Displacement Response Spectrum"
SEMILOGY = True

mean_disp, median_disp = PLOT_2D(1, XLABEL, YLABEL, TITLE, SEMILOGY)

# ----------------------------------------
# Plot 2: Velocity Response Spectrum
# ----------------------------------------
X = T
Y = DATA
XLABEL = "Period T [s]"
YLABEL = "Max Velocity (PGV) [m/s]"
TITLE =  "Velocity Response Spectrum"
SEMILOGY = True

mean_velo, median_velo = PLOT_2D(2, XLABEL, YLABEL, TITLE, SEMILOGY)

# ----------------------------------------
# Plot 3: Acceleration Response Spectrum
# ----------------------------------------
X = T
Y = DATA
XLABEL = "Period T [s]"
YLABEL = "Max Acceleration (PGA) [m/s²]"
TITLE =  "Acceleration Response Spectrum"
SEMILOGY = True

mean_acce, median_acce = PLOT_2D(3, XLABEL, YLABEL, TITLE, SEMILOGY)

# ----------------------------------------
# Plot 4: Base Reaction Spectrum
# ----------------------------------------
X = T
Y = DATA
XLABEL = "Period T [s]"
YLABEL = "Max Base Reaction (N)"
TITLE =  "Base Reaction Response Spectrum"
SEMILOGY = True

mean_reaction, median_reaction = PLOT_2D(4, XLABEL, YLABEL, TITLE, SEMILOGY)

# ----------------------------------------
# Plot 5: Damage Index Spectrum
# ----------------------------------------
X = T
Y = DATA
XLABEL = "Period T [s]"
YLABEL = "Damage Index (DI) [%]"
TITLE =  "Ductility Damage Index Spectrum"
SEMILOGY = False

mean_di, median_di = PLOT_2D(5, XLABEL, YLABEL, TITLE, SEMILOGY)

# ----------------------------------------
# Plot 6: Damping Ratio Spectrum
# ----------------------------------------
X = T
Y = DATA
XLABEL = "Period T [s]"
YLABEL = "Damping Ratio [%]"
TITLE =  "Damping Ratio Spectrum"
SEMILOGY = False

mean_da, median_da = PLOT_2D(6, XLABEL, YLABEL, TITLE, SEMILOGY)

# ----------------------------------------
# Plot 7: Structural Stiffness Spectrum
# ----------------------------------------
X = T
Y = DATA
XLABEL = "Period T [s]"
YLABEL = "Structural Stiffness Spectrum [N/m]"
TITLE =  "Structural Stiffness Spectrum"
SEMILOGY = False

mean_stif, median_stif = PLOT_2D(7, XLABEL, YLABEL, TITLE, SEMILOGY)

#%%------------------------------------------------------------------------------
def PLOT_3D_CONTOUR_XYZ(TAG, X, Y, Z, XLABEL, YLABEL, ZLABEL):
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.interpolate import griddata
    # Convert to NumPy arrays
    X = np.array(X)
    Y = np.array(Y)
    Z = np.array(Z)
    
    # Create grid for contour plot
    xi = np.linspace(min(X), max(X), 100)
    yi = np.linspace(min(Y), max(Y), 100)
    xi, yi = np.meshgrid(xi, yi)
    
    # Interpolate Z values on the grid
    zi = griddata((X, Y), Z, (xi, yi), method='nearest')
    #zi = griddata((X, Y), Z, (xi, yi), method='linear')
    #zi = griddata((X, Y), Z, (xi, yi), method='cubic')
    
    # Plot 3D contour
    fig = plt.figure(TAG, figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')
    contour = ax.plot_surface(xi, yi, zi, cmap='viridis', edgecolor='none')
    
    ax.set_xlabel(XLABEL)
    ax.set_ylabel(YLABEL)
    ax.set_zlabel(ZLABEL,)
    
    fig.colorbar(contour, ax=ax, shrink=0.5, aspect=5)
    plt.title(f'3D Contour Plot of {ZLABEL}')
    plt.show()

#%%------------------------------------------------------------------------------
X = duct
Y = osf
Z = median_disp
XLABEL = 'Structural Ductility Ratio [m/m]'   
YLABEL = 'Structural Over-Strength Factor [N/N]' 
ZLABEL = "Max Displacement (PGD) [m]"
PLOT_3D_CONTOUR_XYZ(0, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = duct
Y = osf
Z = median_velo
XLABEL = 'Structural Ductility Ratio [m/m]'   
YLABEL = 'Structural Over-Strength Factor [N/N]' 
ZLABEL = "Max Velocity (PGV) [m/s]"
PLOT_3D_CONTOUR_XYZ(1, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = duct
Y = osf
Z = median_acce
XLABEL = 'Structural Ductility Ratio [m/m]'   
YLABEL = 'Structural Over-Strength Factor [N/N]' 
ZLABEL = "Max Acceleration (PGA) [m/s²]"
PLOT_3D_CONTOUR_XYZ(2, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = duct
Y = osf
Z = median_reaction
XLABEL = 'Structural Ductility Ratio [m/m]'   
YLABEL = 'Structural Over-Strength Factor [N/N]' 
ZLABEL = "Max Base-rection (PGA) [N]"
PLOT_3D_CONTOUR_XYZ(3, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = duct
Y = osf
Z = median_di
XLABEL = 'Structural Ductility Ratio [m/m]'   
YLABEL = 'Structural Over-Strength Factor [N/N]' 
ZLABEL = "Sructural Ductility Damage Index (DI) [%]"
PLOT_3D_CONTOUR_XYZ(4, X, Y, Z, XLABEL, YLABEL, ZLABEL)
#%%------------------------------------------------------------------------------------------------
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

im_values = median_acce 
X_LABEL = 'Peak Ground Acceleration (g)  [IM]'
S07.FRAGILITY_CURVE(im_values, damage_states, X_LABEL, SEMILOGY=False, PLOT_KIND=False)
#===========================================================
# Define damage state parameters: {Damage State: (median_IM, beta)}
damage_states = {
    'Minor Damage Level': (20.0, 40.0),# Median DI=20, β=40
    'Moderate Damage Level': (40.0, 40.0),
    'Severe Damage Level': (60.0, 50.0),
    'Failure Level': (100.0, 50.0)
}

# Generate intensity measure (IM) values from 0.0 to 100.0
im_values = median_di # Structural Ductility Damage Index
X_LABEL = 'Structural Ductility Damage Index (%)  [IM]'
S07.FRAGILITY_CURVE(im_values, damage_states, X_LABEL, SEMILOGY=False, PLOT_KIND=False)

#%%------------------------------------------------------------------------------------------------
import matplotlib.pyplot as plt

plt.figure(figsize=(12,8))
plt.plot(displacement, base_reaction,color='black')
plt.xlabel("Displacement [m]")
plt.ylabel("Base Reaction [N]")
plt.title("Displacement & Base Reaction Relation From Last Dynamic Analysis")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
#%%------------------------------------------------------------------------------------------------
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
#%%------------------------------------------------------------------------------------------------  
# EXPORT DATA TO EXCEL
DATA_TOTAL = {
    'Spring_Stiffness': Ke,
    'Yield_Force_of_Structure': FY,
    'Ultimate_Force_of_Structure': FU,
    'Yield_Displacement': DY,
    'Ultimate_Displacement': DSU,
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
results_df.to_excel('INELASTIC_RESPONSE_SPECTRUM_SEISMIC_SDOF_LAST_ANALYSIS_DATA_RESULTS.xlsx', index=False)
#%%------------------------------------------------------------------------------------------------  
XLABEL = 'Displacement'
YLABEL = 'Base Reaction'
TITLE = f'{YLABEL} and {XLABEL} scatter chart'
COLOR = 'orange'
X = median_disp
Y = median_reaction
S01.PLOT_SCATTER(X, Y , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

# CLUSTER DATA
S01.CLUSTER_DATA(X, Y, XLABEL, YLABEL, MAX_CLUSTERS=3)
#%%------------------------------------------------------------------------------------------------
XLABEL = 'Velocity'
YLABEL = 'Base Reaction'
TITLE = f'{YLABEL} and {XLABEL} scatter chart'
COLOR = 'cyan'
X = median_velo
Y = median_reaction
S01.PLOT_SCATTER(X, Y , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

# CLUSTER DATA
S01.CLUSTER_DATA(X, Y, XLABEL, YLABEL, MAX_CLUSTERS=3)
#%%------------------------------------------------------------------------------------------------
XLABEL = 'Acceleration'
YLABEL = 'Base Reaction'
TITLE = f'{YLABEL} and {XLABEL} scatter chart'
COLOR = 'lime'
X = median_acce
Y = median_reaction
S01.PLOT_SCATTER(X, Y , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

# CLUSTER DATA
S01.CLUSTER_DATA(X, Y, XLABEL, YLABEL, MAX_CLUSTERS=3)
#%%------------------------------------------------------------------------------------------------
XLABEL = 'Displacement'
YLABEL = 'Structural Ductility Damage Index (%)'
TITLE = f'{YLABEL} and {XLABEL} scatter chart'
COLOR = 'purple'
X = median_disp
Y = median_di
S01.PLOT_SCATTER(X, Y , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

# CLUSTER DATA
S01.CLUSTER_DATA(X, Y, XLABEL, YLABEL, MAX_CLUSTERS=3)
#%%------------------------------------------------------------------------------------------------
# PLOT THE TIME-HISTORY
S01.PLOT_TIME_HISTORY(time, displacement, velocity, acceleration, base_reaction)
#%%------------------------------------------------------------------------------------------------

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
    'Max_displacement': median_disp,
    'Max_velocity': median_velo,
    'Max_acceleration': median_acce,
    'Max_Base_Reaction': median_reaction,
    'Ductility_Damage_Index': median_di,
}


# Convert to DataFrame
df = pd.DataFrame(data)
#print(df)
S01.RANDOM_FOREST(df)
#%%------------------------------------------------------------------------------------------------
# PLOT HEATMAP FOR CORRELATION 
S01.PLOT_HEATMAP(df)
#%%------------------------------------------------------------------------------------------------
# MULTIPLE REGRESSION MODEL
S01.MULTIPLE_REGRESSION(df) 
"""
#%%------------------------------------------------------------------------------------------------
# MACHINE LEARNING: LONG SHORT-TREM MEMERY (LSTM) METHOD
x = max_displacement 
y = max_acceleration 
Demand_X = x[-1]
look_back = 500#int(NUM_SIM * 0.5)
ITERATION = 200
XLABEL = 'Max Displacement'
YLABEL = 'Max Acceleration'
#S01.PREDICT_LSTM(x, y, Demand_X, look_back, ITERATION, XLABEL, YLABEL)
#%%------------------------------------------------------------------------------------------------
# PERFORM RELIABILITY ANALYSIS FOR BASE REACTION AND ELEMENT CAPACITY
mean_capacity = np.mean(FU)    # Mean Element Ultimate Capacity
std_dev_capacity = np.std(FU)  # Std Element Ultimate Capacity
num_sim = NUM_SEISMIC
S01.RELIABILITY_ANALYSIS(max_base_reaction, num_sim, mean_capacity, std_dev_capacity)
#%%------------------------------------------------------------------------------------------------
# NEURAL NETWORK FOR FAILURE PROBABILIYY ESTIMATION
X1 = mean_capacity
X2 = max_base_reaction
S01.NEURAL_NETWORK_FAILURE_PROBABILIYY_ESTIMATION(max_base_reaction, X2, NUM_SEISMIC)
#%%------------------------------------------------------------------------------------------------
# MARKOV CHAIN MODEl (structural damage analysis by evaluating Structural Ductility Damage Index)
FILE_TF = False         # Indicate whether to read data from a file or use provided data
file_path = None        # Not used when 'file_tf' is False
DATA = max_DI # If not using a file, replace None with a NumPy array of data

S03.MARKOV_CHAIN(FILE_TF, file_path, DATA)
#------------------------------------------------------------------------------------------------

"""
