###########################################################################################################
#                          >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<               #
#               SENSITIVITY ANALYSIS OF SINGLE-DEGREE-FREEOM (SDOF) STRUCTURES USING FREE-VIBRATION       #
#           EFFECTS OF INITIAL DISPLACEMENT, MASS, STRCTURAL DUCTILITY RATIO AND OVER-STRENGTH FACTOR     #
#                ON OUTPUT KEY PARAMETERS FROM NONLINEAR DYNAMIC ANALYSES USING PYTHON AND OPENSEES       #
#---------------------------------------------------------------------------------------------------------#
#                                 FREE VIBRATION ANALYSIS USING INITIAL DISPLACEMENT                      #
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
#                          THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                     #
#                                   EMAIL: salar.d.ghashghaei@gmail.com                                   #
###########################################################################################################
"""
1. This script performs a "sensitivity analysis of a single-degree-of-freedom (SDOF) structure" using "free vibration" in OpenSeesPy.
2. The objective is to study the effects of "initial displacement, mass, ductility ratio, and over-strength factor" on nonlinear dynamic response.
3. Structural properties such as "yield force, ultimate force, elastic stiffness, hardening ratio, mass, and damping ratio" are defined first.
4. Elastic and plastic periods of the system are computed and reported.
5. Key seismic performance parameters are calculated, including "over-strength (Ω₀)", "ductility (μ)", and "behavior factor (R)".
6. A nonlinear "hysteretic spring model" is used to represent structural behavior with strength degradation.
7. A "viscous damper" is added to simulate energy dissipation.
8. The analysis considers "free vibration due to initial displacement" (no external ground motion).
9. The function `ANALYSIS_SDOF` builds the OpenSees model, applies initial conditions, and runs a transient analysis.
10. "Newmark time integration" and "Newton–Raphson iteration" are used for nonlinear solution.
11. "Rayleigh damping coefficients" are computed from the first eigenvalue.
12. Time histories of "displacement, velocity, acceleration, base reaction, stiffness, and damage index" are recorded.
13. A "ductility-based damage index" is calculated at each time step.
14. The effective stiffness degradation is monitored during the response.
15. Multiple simulations are executed by looping over ranges of "mass, initial displacement, ductility, and over-strength".
16. For each simulation, "maximum response values" are extracted.
17. Damage index values are limited between "0% and 100%".
18. All results are stored for post-processing and statistical analysis.
19. The script reports completion of each simulation and total runtime.
20. Overall, the code provides a "parametric nonlinear dynamic assessment" of SDOF systems under free vibration conditions.
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
import OPENSEEES_HYSTERETICSM_FORCE_DISP_FUN as S08
from scipy.stats import norm
import pickle
from joblib import Parallel, delayed
from tqdm import tqdm
import os

#%%------------------------------------------------------------------------------------------------
# Define parameters (units: m, N)
NUM_SIM = 6                                    # Total number for analysis
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

duration = 50.0                                  # [s] Total simulation duration
dt = 0.01                                        # [s] Time step

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
def ANALYSIS_SDOF(m, duct, osf, dy, IU=True):
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
    
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(2, 1.0)
    
    if IU == True:
        # Define initial displacment
        ops.setNodeDisp(2, 1, dy, '-commit') # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/setNodeDisp.html
    """
    if IV == True:
        # Define initial velocity
        ops.setNodeVel(2, 1, v0, '-commit')  # INFO LINK: https://openseespydoc.readthedocs.io/en/stable/src/setNodeVel.html
    if IA == True:
        # Define initial  acceleration
        ops.setNodeAccel(2, 1, a0, '-commit') # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/setNodeAccel.html
    """
    
    # Output data
    #ops.recorder('Node', '-file', f"DTH_DYN_{i}.txt",'-time', '-node', 2, '-dof', 1, 'disp')     # Displacement Time History Node 2
    #ops.recorder('Node', '-file', f"VTH_DYN_{i}.txt",'-time', '-node', 2, '-dof', 1, 'vel')      # Velocity Time History Node 2
    #ops.recorder('Node', '-file', f"ATH_DYN_{i}.txt",'-time', '-node', 2, '-dof', 1, 'accel')    # Acceleration Time History Node 2
    #ops.recorder('Node', '-file', f"BTH_DYN_{i}.txt",'-time', '-node', 1, '-dof', 1, 'reaction') # Base Reaction Time History Node 1
        
    # Set analysis parameters
    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.system('BandGeneral')
    ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/normDispIncr.html
    #ops.integrator('CentralDifference')  # JUST FOR LINEAR ANALYSIS - INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/centralDifference.html
    alpha=0.5; beta=0.5;
    ops.integrator('Newmark', alpha, beta) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/newmark.html
    #alpha=2/3;gamma=1.5-alpha; gamma=1.5-alpha;beta=(2-alpha)**2/4;
    #ops.integrator('HHT', alpha, gamma, beta) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/hht.html
    ops.algorithm('Newton')  # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/algorithm.html
    ops.analysis('Transient') # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/analysis.html
    
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
        acceleration.append(ops.nodeAccel(2, 1))                     # Structure Acceleration
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
# Parameters
CHECKPOINT_FILE = 'SDOF_RESULTS.pkl'

dy   = np.linspace(0.1 * DY, 1.5 * DY, NUM_SIM)
duct = np.linspace(2.0, 8.5, NUM_SIM)
osf  = np.linspace(1.05, 1.35, NUM_SIM)
m    = np.linspace(mi * (1 / NUM_SIM) * 0.04, mi * 0.04, NUM_SIM)

# Generate All Simulation Cases
cases = [
    (m_i, duct_k, osf_l, dy_j)
    for dy_j in dy
    for m_i in m
    for duct_k in duct
    for osf_l in osf
]

TOTAL_SIM = len(cases)
print(f'Total simulations: {TOTAL_SIM}')

# Single Simulation Function (Thread-Safe)
def RUN_SIMULATION(mass, ductility, osf_factor, dy0):

    data = ANALYSIS_SDOF(
        mass,
        ductility,
        osf_factor,
        dy0,
        IU=True
    )

    (time, displacement, velocity, acceleration,
     base_reaction, DI, PERIOD, damping_ratio, stiffness) = data

    DI_max = np.clip(np.max(DI), 0.0, 100.0)

    return {
        'Mass': mass,
        'Dy': dy0,
        'Ductility': ductility,
        'OSF': osf_factor,
        'Max_Time': np.max(np.abs(time)),
        'PGD': np.max(np.abs(displacement)),
        'PGV': np.max(np.abs(velocity)),
        'PGA': np.max(np.abs(acceleration)),
        'Base_Reaction': np.max(np.abs(base_reaction)),
        'Stiffness': np.max(np.abs(stiffness)),
        'Period': PERIOD,
        'Damping_Ratio': damping_ratio,
        'Damage_Index': DI_max
    }

# ----------------------  PARALLEL PROCESSING  ----------------------
# Analysis Durations:
current_time = TI.strftime("%H:%M:%S", TI.localtime())
print("Start Time:", current_time)

# Parallel Execution + Checkpoint
if os.path.exists(CHECKPOINT_FILE):
    print('Loading checkpoint...')
    with open(CHECKPOINT_FILE, 'rb') as f:
        results = pickle.load(f)
else:
    print('Running parallel simulations...')
    results = Parallel(
        n_jobs=-1,
        backend='loky',
        verbose=10
    )(
        delayed(RUN_SIMULATION)(*case)
        for case in tqdm(cases)
    )

    with open(CHECKPOINT_FILE, 'wb') as f:
        pickle.dump(results, f)

df = pd.DataFrame(results)
print(df.head())

current_time = TI.strftime("%H:%M:%S", TI.localtime())
print("Finish Time:", current_time)
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
    #zi = griddata((X, Y), Z, (xi, yi), method='nearest')
    #zi = griddata((X, Y), Z, (xi, yi), method='linear')
    zi = griddata((X, Y), Z, (xi, yi), method='cubic')
    
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
X = df['Ductility']
Y = df['Mass']
Z = df['PGD']
XLABEL = 'Structural Ductility Ratio [m/m]'   
YLABEL = 'Structural Mass [kg]' 
ZLABEL = "Max Displacement (PGD) [m]"
PLOT_3D_CONTOUR_XYZ(0, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = df['Ductility']
Y = df['Mass']
Z = df['PGV']
XLABEL = 'Structural Ductility Ratio [m/m]'   
YLABEL = 'Structural Mass [kg]' 
ZLABEL = "Max Velocity (PGV) [m/s]"
PLOT_3D_CONTOUR_XYZ(1, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = df['Ductility']
Y = df['Mass']
Z = df['PGA']
XLABEL = 'Structural Ductility Ratio [m/m]'   
YLABEL = 'Structural Mass [kg]' 
ZLABEL = "Max Acceleration (PGA) [m/s²]"
PLOT_3D_CONTOUR_XYZ(2, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = df['Ductility']
Y = df['Mass']
Z = df['Base_Reaction']
XLABEL = 'Structural Ductility Ratio [m/m]'   
YLABEL = 'Structural Mass [kg]' 
ZLABEL = "Max Base-rection (PGA) [N]"
PLOT_3D_CONTOUR_XYZ(3, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = df['Ductility']
Y = df['Mass']
Z = df['Damage_Index']
XLABEL = 'Structural Ductility Ratio [m/m]'   
YLABEL = 'Structural Mass [kg]' 
ZLABEL = "Sructural Ductility Damage Index (DI) [%]"
PLOT_3D_CONTOUR_XYZ(4, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = df['Ductility']
Y = df['Mass']
Z = df['Period']
XLABEL = 'Structural Ductility Ratio [m/m]'   
YLABEL = 'Structural Mass [kg]' 
ZLABEL = "Sructural Period [s]"
PLOT_3D_CONTOUR_XYZ(5, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = df['Ductility']
Y = df['Mass']
Z = df['Damping_Ratio']
XLABEL = 'Structural Ductility Ratio [m/m]'   
YLABEL = 'Structural Mass [kg]' 
ZLABEL = "Sructural Damping Ratio [%]"
PLOT_3D_CONTOUR_XYZ(6, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = df['Ductility']
Y = df['Period']
Z = df['Damping_Ratio']
XLABEL = 'Structural Ductility Ratio [m/m]'   
YLABEL = 'Sructural Period [s]' 
ZLABEL = "Sructural Damping Ratio [%]"
PLOT_3D_CONTOUR_XYZ(7, X, Y, Z, XLABEL, YLABEL, ZLABEL)
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

im_values = df['PGA'] 
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
im_values = df['Damage_Index'] # Structural Ductility Damage Index
X_LABEL = 'Structural Ductility Damage Index (%)  [IM]'
S07.FRAGILITY_CURVE(im_values, damage_states, X_LABEL, SEMILOGY=False, PLOT_KIND=False)

#%%------------------------------------------------------------------------------------------------