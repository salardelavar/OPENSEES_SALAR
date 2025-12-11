###########################################################################################################
#                    >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                     #
#       NONLINEAR DYNAMIC ANALYSIS OF A TWO-DOF SOIL–STRUCTURE SYSTEM UNDER 20 GROUND MOTIONS WITH        #
#                RESPONSE SPECTRA AND DUCTILITY DAMAGE INDEX CALCULATION IN OPENSEES                      #
#---------------------------------------------------------------------------------------------------------#
#                          THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                     #
#                                   EMAIL: salar.d.ghashghaei@gmail.com                                   #
###########################################################################################################

"""
This code implements a comprehensive nonlinear dynamic analysis framework for
performance-based earthquake engineering assessment of multi-degree-of-freedom
(MDOF) systems. The methodology combines traditional nonlinear time-history
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

Conclusions:
 - For inelastic SDOF studies, including pinching and degradation in the hysteretic model can change predicted peak responses, residual displacements, and damage indices significantly.
 - Use HYSTERETIC-type models for damage-sensitive or collapse-prone scenarios, and calibrate degradation parameters using test data where possible.
Objective:
The study evaluates the dynamic response of a single-degree-of-freedom (SDOF) inelastic system under seismic excitation,
comparing two hysteretic models for the restoring force:
 - HYSTERETICSM (multi-linear / pinching / stiffness degradation including ultimate strain)

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
---------------------------------------------
RECOMMENDATIONS FOR PRACTICE:

1. For serviceability assessment: Elastic analysis may suffice for DI < 20%
2. For life safety evaluation: Nonlinear analysis with degradation essential
3. For collapse prevention: Advanced hysteretic models with ultimate criteria required
4. For design optimization: Response spectrum analysis across period range recommended

LIMITATIONS AND FUTURE ENHANCEMENTS:

1. Material calibration: Laboratory test data needed for parameters
2. Ground motions: Site-specific spectra recommended
3. Validation: Experimental correlation required for confidence
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
from scipy.stats import norm

#%%------------------------------------------------------------------------------------------------
# Define parameters (units: m, N)
NUM_PERIOD = 100                                # Total number for Period in each simulation
NUM_SEISMIC = 20                                # Total number for seismic simulation
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

# Define Soil Properties
Soil_Stiffness = 100000000.0                     # [N/m] Soil Stiffness (No-Tension)
mSOIL = 9400000.0                                # [kg] Mass of the Soil

duration = 15.0  # [s] Total simulation duration
dt = 0.01        # [s] Time step
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
def ANALYSIS_MDOF(i, mSTRUCTURE):
    # Initialize OpenSees model
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)
    GMfact = 9.81 # [m/s^2] standard acceleration of gravity or standard acceleration 
        
    # Define nodes
    ops.node(1, 0.0)  # Fixed base
    ops.node(2, 0.0)  # Soil node
    ops.node(3, 0.0)  # Structure node
        
    # Define boundary conditions
    ops.fix(1, 1)

    # Define mass
    ops.mass(2, mSOIL)      # Soil
    ops.mass(2, mSTRUCTURE) # Structure
        
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
    #%% SOIL
    Soil_Tag = 1
    alpha = 0.25    # velocity exponent (usually 0.3–1.0)
    omega = np.sqrt(Soil_Stiffness/mSOIL)
    C_soil = 2 * DR * omega * mSTRUCTURE  # [N·s/m] Damping coefficient 
    ops.uniaxialMaterial('ENT', Soil_Tag, Soil_Stiffness)
    #INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Elastic-No_Tension_Material
    Soil_Damping_Tag = 2
    ops.uniaxialMaterial('Viscous', Soil_Damping_Tag, C_soil, alpha)  # Material for C (alpha=1.0 for linear)
    #INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Viscous_Material
    #%% STRUCTURE
    Structure_Tag = 3
    ops.uniaxialMaterial('HystereticSM', Structure_Tag, '-posEnv', *KP.flatten(), '-negEnv', *KN.flatten(), '-pinch', 1, 1)
    #INFO LINK: https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/HystereticSM.html
    alpha = 1.0    # velocity exponent (usually 0.3–1.0)
    omega = np.sqrt(Ke/mSTRUCTURE)
    C_stru = 2 * DR * omega * mSTRUCTURE  # [N·s/m] Damping coefficient 
    
    Structure_Damping_Tag = 4
    ops.uniaxialMaterial('Viscous', Structure_Damping_Tag, C_stru, alpha)  # Material for C (alpha=1.0 for linear)
    
    # Define element
    #ops.element('zeroLength', 1, 1, 2, '-mat', Soil_Tag, Soil_Damping_Tag, '-dir', 1, 1)            # DOF[1] SOIL LATERAL SPRING
    #ops.element('zeroLength', 2, 2, 3, '-mat', Structure_Tag, Structure_Damping_Tag, '-dir', 1, 1)  # DOF[2] STRUCTURE LATERAL SPRING
    #%% SOIL
    ops.element('zeroLength', 1, 1, 2, '-mat', Soil_Tag, '-dir', 1)                # DOF[1] SOIL LATERAL SPRING
    ops.element('zeroLength', 2, 1, 2, '-mat', Soil_Damping_Tag, '-dir', 1)        # DOF[1] SOIL LATERAL SPRING
    #%% STRUCTURE
    ops.element('zeroLength', 3, 2, 3, '-mat', Structure_Tag, '-dir', 1)                # DOF[2] STRCUTURE LATERAL SPRING
    ops.element('zeroLength', 4, 2, 3, '-mat', Structure_Damping_Tag, '-dir', 1)        # DOF[2] STRCUTURE LATERAL SPRING
    
    #%% Apply seismic accelerations    
    # Define time series for input motion (Acceleration time history)
    gm_accels = np.loadtxt(f'Ground_Acceleration_{i+1}.txt')  # Assumes acceleration in m/s²
    ops.timeSeries('Path', 1, '-dt', dt, '-values', *gm_accels.tolist(), '-factor', GMfact) # SEISMIC-X
    #ops.timeSeries('Path', 1, '-dt', dt, '-filePath', f'Ground_Acceleration{i+1}.txt', '-factor', GMfact) # SEISMIC-X
        
    #%% Define load patterns
    # pattern UniformExcitation $patternTag $dof -accel $tsTag <-vel0 $vel0> <-fact $cFact>
    SSF = 5.00 # Seismic Scale Factor
    ops.pattern('UniformExcitation', 200, 1, '-accel', 1, '-fact', SSF) # SEISMIC-X
    
    #%% Output data
    #ops.recorder('Node', '-file', f"DTH_DYN_{i}.txt",'-time', '-node', 2, '-dof', 1, 'disp')     # Displacement Time History Node 2
    #ops.recorder('Node', '-file', f"VTH_DYN_{i}.txt",'-time', '-node', 2, '-dof', 1, 'vel')      # Velocity Time History Node 2
    #ops.recorder('Node', '-file', f"ATH_DYN_{i}.txt",'-time', '-node', 2, '-dof', 1, 'accel')    # Acceleration Time History Node 2
    #ops.recorder('Node', '-file', f"BTH_DYN_{i}.txt",'-time', '-node', 1, '-dof', 1, 'reaction') # Base Reaction Time History Node 1
        
    #%% Set analysis parameters
    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.system('BandGeneral')
    ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS)
    ops.algorithm('Newton')
    ops.integrator('Newmark', 0.5, 0.25)
    ops.analysis('Transient')
    
    #%% Calculate Rayleigh damping factors
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
    
    #%% Run Dynamic Analysis
    time = []
    displacement_soil = []
    displacement_stru = []
    velocity_soil = []
    velocity_stru = []
    acceleration_soil = []
    acceleration_stru = []
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
        displacement_soil.append(ops.nodeDisp(2, 1))                      # Soil Displacement
        displacement_stru.append(ops.nodeDisp(3, 1))                      # Structure Displacement
        velocity_soil.append(ops.nodeVel(2, 1))
        velocity_stru.append(ops.nodeVel(3, 1))
        acceleration_soil.append(ops.nodeAccel(2, 1))                     # Soil Acceleration
        if step <= len(gm_accels)-1:
            acceleration_stru.append(ops.nodeAccel(3, 1) + gm_accels[step]) # Structure Acceleration and Seismic Acceleration
        else:
            acceleration_stru.append(ops.nodeAccel(3, 1))                 # Structure Acceleration
        base_reaction.append(-ops.eleResponse(2, 'force')[0])             # Reaction force
        DI.append(100*(displacement_stru[-1] - DY) / (DSU - DY))          # [%] Structural Ductility Damage Index 
        stiffness.append(np.abs(base_reaction[-1]/displacement_stru[-1])) # [N/m] Stiffness
        # IN EACH STEP STRUCTURAL PERIOD WILL BE CALCULATED
        #PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(1, PLOT=True)
        #PERIOD_MIN.append(PERIODmin)
        #PERIOD_MAX.append(PERIODmax)
        step += 1
        #print(f'{time[-1]} {displacement[-1]:0.4e} {base_reaction[-1]:0.4e}')
        
    # Calculate Structure Damping Ratio Based on Lateral Displacement
    damping_ratio = S06.DAMPING_RATIO(displacement_stru) 
    #ops.wipe()
    OUTPUT = (time, displacement_soil, displacement_stru,
            velocity_soil, velocity_stru, acceleration_soil,
            acceleration_stru, base_reaction, DI, PERIOD,
            damping_ratio, stiffness)
    return OUTPUT

#%%------------------------------------------------------------------------------------------------
# Analysis Durations:
starttime = ti.process_time()

# Collect Data
DATA = {
    1: [],  # DISP - SOIL
    2: [],  # VELOCITY - SOIL
    3: [],  # ACCELERAION - SOIL
    4: [],  # DISP - STRUCTURE
    5: [],  # VELOCITY - STRUCTURE
    6: [],  # ACCELERAION - STRUCTURE
    7: [],  # REACTION
    8: [],  # DI
    9: [],  # DAMPING RATIO
    10: [],  # STRUCTURAL STIFFNESS
    }

# NUM_SIM is the number of simulations
for j in range(NUM_SEISMIC):
    # Initialize lists to store max values
    max_time = []
    max_displacement_soil = []
    max_velocity_soil = []
    max_acceleration_soil = []
    max_displacement_stru = []
    max_velocity_stru = []
    max_acceleration_stru = []
    max_base_reaction = []
    max_DI = []
    max_T = []
    max_K = []
    DAMPING_RATIO = []
    for i in range(NUM_PERIOD):
        m = mi * (i+1 / NUM_PERIOD) * 0.2
        data = ANALYSIS_MDOF(j, m)
        (time, displacement_soil, displacement_stru,
                velocity_soil, velocity_stru, acceleration_soil,
                acceleration_stru, base_reaction, DI, PERIOD,
                damping_ratio, stiffness) = data
        # Calculate and store the max absolute values
        # SOIL
        max_time.append(np.max(np.abs(time)))
        max_displacement_soil.append(np.max(np.abs(displacement_soil)))
        max_velocity_soil.append(np.max(np.abs(velocity_soil)))
        max_acceleration_soil.append(np.max(np.abs(acceleration_soil)))
        # STRUCTURE
        max_displacement_stru.append(np.max(np.abs(displacement_stru)))
        max_velocity_stru.append(np.max(np.abs(velocity_stru)))
        max_acceleration_stru.append(np.max(np.abs(acceleration_stru)))
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
    # SOIL
    DATA[1].append(max_displacement_soil)
    DATA[2].append(max_velocity_soil)
    DATA[3].append(max_acceleration_soil)
    # STRUCTURE
    DATA[4].append(max_displacement_stru)
    DATA[5].append(max_velocity_stru)
    DATA[6].append(max_acceleration_stru)
    DATA[7].append(max_base_reaction)
    DATA[8].append(max_DI)   
    DATA[9].append(DAMPING_RATIO)
    DATA[10].append(max_K)
    print(f'\n SEISMIC {j + 1} DONE \n')    
else:
    print('\n\n Analysis completed successfully')    

totaltime = ti.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n') 

#%%------------------------------------------------------------------------------------------------
# OUTPUT DATA FROM RESPONSE SPECTRUM ANALYSIS TO EXCEL FILE
labels = {
    1: "SOIL-DISPLACEMENT",
    2: "SOIL-VELOCITY",
    3: "SOIL-ACCELERATION",
    4: "STRUCTURE-DISPLACEMENT",
    5: "STRUCTURE--VELOCITY",
    6: "STRUCTURE--ACCELERATION",
    7: "BASE_REACTION",
    8: "DAMAGE_INDEX",
    9: "DAMPING_RATIO",
    10: "STIFFNESS",
}

with pd.ExcelWriter("SOIL_INELASTIC_RESPONSE_SPECTRUM_SEISMIC_MDOF_RESULTS.xlsx", engine="openpyxl") as writer:
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
    for j in range(NUM_SEISMIC):
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

    
    
# ----------------------------------------
# Plot 1: Soil Displacement Response Spectrum
# ----------------------------------------
X = T
Y = DATA
XLABEL = "Period T [s]"
YLABEL = "Soil Max Displacement (PGD) [m]"
TITLE =  "Soil Displacement Response Spectrum"
SEMILOGY = False

PLOT_2D(1, XLABEL, YLABEL, TITLE, SEMILOGY)

# ----------------------------------------
# Plot 2: Soil Velocity Response Spectrum
# ----------------------------------------
X = T
Y = DATA
XLABEL = "Period T [s]"
YLABEL = "Soil Max Velocity (PGV) [m/s]"
TITLE =  "Soil Velocity Response Spectrum"
SEMILOGY = False

PLOT_2D(2, XLABEL, YLABEL, TITLE, SEMILOGY)

# ----------------------------------------
# Plot 3: Soil Acceleration Response Spectrum
# ----------------------------------------
X = T
Y = DATA
XLABEL = "Period T [s]"
YLABEL = "Soil Max Acceleration (PGA) [m/s²]"
TITLE =  "Soil Acceleration Response Spectrum"
SEMILOGY = False

PLOT_2D(3, XLABEL, YLABEL, TITLE, SEMILOGY)

# ----------------------------------------
# Plot 4: Structure Displacement Response Spectrum
# ----------------------------------------
X = T
Y = DATA
XLABEL = "Period T [s]"
YLABEL = "Structure Max Displacement (PGD) [m]"
TITLE =  "Structure Displacement Response Spectrum"
SEMILOGY = False

PLOT_2D(4, XLABEL, YLABEL, TITLE, SEMILOGY)

# ----------------------------------------
# Plot 5: Structure Velocity Response Spectrum
# ----------------------------------------
X = T
Y = DATA
XLABEL = "Period T [s]"
YLABEL = "Structure Max Velocity (PGV) [m/s]"
TITLE =  "Structure Velocity Response Spectrum"
SEMILOGY = False

PLOT_2D(5, XLABEL, YLABEL, TITLE, SEMILOGY)

# ----------------------------------------
# Plot 6: Structure Acceleration Response Spectrum
# ----------------------------------------
X = T
Y = DATA
XLABEL = "Period T [s]"
YLABEL = "Structure Max Acceleration (PGA) [m/s²]"
TITLE =  "Structure Acceleration Response Spectrum"
SEMILOGY = False

PLOT_2D(6, XLABEL, YLABEL, TITLE, SEMILOGY)

# ----------------------------------------
# Plot 7: Base Reaction Spectrum
# ----------------------------------------
X = T
Y = DATA
XLABEL = "Period T [s]"
YLABEL = "Max Base Reaction (N)"
TITLE =  "Base Reaction Response Spectrum"
SEMILOGY = True

PLOT_2D(7, XLABEL, YLABEL, TITLE, SEMILOGY)

# ----------------------------------------
# Plot 8: Damage Index Spectrum
# ----------------------------------------
X = T
Y = DATA
XLABEL = "Period T [s]"
YLABEL = "Damage Index (DI) [%]"
TITLE =  "Ductility Damage Index Spectrum"
SEMILOGY = True

PLOT_2D(8, XLABEL, YLABEL, TITLE, SEMILOGY)

# ----------------------------------------
# Plot 9: Damping Ratio Spectrum
# ----------------------------------------
X = T
Y = DATA
XLABEL = "Period T [s]"
YLABEL = "Damping Ratio [%]"
TITLE =  "Damping Ratio Spectrum"
SEMILOGY = True

PLOT_2D(9, XLABEL, YLABEL, TITLE, SEMILOGY)

# ----------------------------------------
# Plot 10: Structural Stiffness Spectrum
# ----------------------------------------
X = T
Y = DATA
XLABEL = "Period T [s]"
YLABEL = "Structural Stiffness Spectrum [N/m]"
TITLE =  "Structural Stiffness Spectrum"
SEMILOGY = True

PLOT_2D(10, XLABEL, YLABEL, TITLE, SEMILOGY)

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
plt.plot(time, acceleration_stru, lw=1, color='black')
plt.xlabel('Time (s)')
plt.ylabel('Acceleration (g)')
plt.title(f'Last Analysis Structural Response + Ground Motion ::: MAX. ABS. : {np.max(np.abs(acceleration_stru)):.4f}')
plt.grid(True)
plt.show()    

im_values = max_acceleration_stru 
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
im_values = max_DI # Structural Ductility Damage Index
X_LABEL = 'Structural Ductility Damage Index (%)  [IM]'
S07.FRAGILITY_CURVE(im_values, damage_states, X_LABEL, SEMILOGY=False, PLOT_KIND=False)

#%%------------------------------------------------------------------------------------------------
import matplotlib.pyplot as plt

plt.figure(figsize=(12,8))
plt.plot(displacement_stru, base_reaction,color='black')
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

DISP_Z = MAX_ABS(displacement_stru)  
VELO_Z = MAX_ABS(velocity_stru) 
ACCE_Z = MAX_ABS(acceleration_stru) 
BASE_Z = MAX_ABS(base_reaction) 

plt.figure(1, figsize=(8, 6))
plt.plot(time, displacement_stru, color='blue', linewidth=2)
plt.plot(time, DISP_Z, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Displacement [m]')
plt.title(f'Time vs Displacement - MAX. ABS: {DISP_Z[-1]}')
plt.grid()
plt.show()

plt.figure(2, figsize=(8, 6))
plt.plot(time, velocity_stru, color='blue', linewidth=2)
plt.plot(time, VELO_Z, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Velocity [m/s]')
plt.title(f'Time vs Velocity - MAX. ABS: {VELO_Z[-1]}')
plt.grid()
plt.show()

plt.figure(3, figsize=(8, 6))
plt.plot(time, acceleration_stru, color='blue', linewidth=2)
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
    'Max_displacement_SOIL': max_displacement_soil,
    'Max_velocity_SOIL': max_velocity_soil,
    'Max_acceleration_SOIL': max_acceleration_soil,
    'Max_displacement_STRUCTURE': max_displacement_stru,
    'Max_velocity_STRUCTURE': max_velocity_stru,
    'Max_acceleration_STRUCTURE': max_acceleration_stru,
    'Max_Base_Reaction': max_base_reaction,
    'Ductility_Damage_Index': max_DI,
    'Period': max_T,
}
# Convert to DataFrame
results_df = pd.DataFrame(DATA_TOTAL)
# Export the DataFrame to an Excel file
results_df.to_excel('SOIL_INELASTIC_RESPONSE_SPECTRUM_SEISMIC_MDOF_LAST_ANALYSIS_DATA_RESULTS.xlsx', index=False)
#%%------------------------------------------------------------------------------------------------  
XLABEL = 'Displacement'
YLABEL = 'Base Reaction'
TITLE = f'{YLABEL} and {XLABEL} scatter chart'
COLOR = 'orange'
X = max_displacement_stru
Y = max_base_reaction
S01.PLOT_SCATTER(X, Y , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

# CLUSTER DATA
S01.CLUSTER_DATA(X, Y, XLABEL, YLABEL, MAX_CLUSTERS=3)
#%%------------------------------------------------------------------------------------------------
XLABEL = 'Velocity'
YLABEL = 'Base Reaction'
TITLE = f'{YLABEL} and {XLABEL} scatter chart'
COLOR = 'cyan'
X = max_velocity_stru
Y = max_base_reaction
S01.PLOT_SCATTER(X, Y , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

# CLUSTER DATA
S01.CLUSTER_DATA(X, Y, XLABEL, YLABEL, MAX_CLUSTERS=3)
#%%------------------------------------------------------------------------------------------------
XLABEL = 'Acceleration'
YLABEL = 'Base Reaction'
TITLE = f'{YLABEL} and {XLABEL} scatter chart'
COLOR = 'lime'
X = max_acceleration_stru
Y = max_base_reaction
S01.PLOT_SCATTER(X, Y , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

# CLUSTER DATA
S01.CLUSTER_DATA(X, Y, XLABEL, YLABEL, MAX_CLUSTERS=3)
#%%------------------------------------------------------------------------------------------------
XLABEL = 'Displacement'
YLABEL = 'Structural Ductility Damage Index (%)'
TITLE = f'{YLABEL} and {XLABEL} scatter chart'
COLOR = 'purple'
X = max_displacement_stru
Y = max_DI
S01.PLOT_SCATTER(X, Y , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

# CLUSTER DATA
S01.CLUSTER_DATA(X, Y, XLABEL, YLABEL, MAX_CLUSTERS=3)
#%%------------------------------------------------------------------------------------------------
# PLOT THE TIME-HISTORY
S01.PLOT_TIME_HISTORY(time, displacement_stru, velocity_stru, acceleration_stru, base_reaction)
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
    'Max_displacement': max_displacement_stru,
    'Max_velocity': max_velocity_stru,
    'Max_acceleration': max_acceleration_stru,
    'Max_Base_Reaction': max_base_reaction,
    'Damping_Ratio': DR,
    'Ductility_Damage_Index': max_DI,
    'Period': max_T,
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

