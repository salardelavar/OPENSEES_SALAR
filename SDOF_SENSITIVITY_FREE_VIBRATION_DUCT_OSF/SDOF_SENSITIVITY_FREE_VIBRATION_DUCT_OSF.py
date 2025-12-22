###########################################################################################################
#                          >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<               #
#               SENSITIVITY ANALYSIS OF SINGLE-DEGREE-FREEOM (SDOF) STRUCTURES USING FREE-VIBRATION       #
#            EFFECTS OF INTIAL DISPLACEMENT, MASS, STRCTURAL DUCTILITY RATIO AND OVER-STRENGTH FACTOR     #
#                ON OUTPUT KEY PARAMETERS FROM NONLINEAR DYNAMIC ANALYSES USING PYTHON AND OPENSEES       #
#---------------------------------------------------------------------------------------------------------#
#                                 FREE VIBRATION ANALYSIS USING INITIAL DISPLACEMENT                      #
#---------------------------------------------------------------------------------------------------------#
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
# Analysis Durations:
starttime = ti.process_time()

# Collect Data
# Initialize lists to store max values
max_time = []
max_displacement = []
max_velocity = []
max_acceleration = []
max_base_reaction = []
max_di = []
max_T = []
max_K = []
max_mass = []
max_duct = []
max_osf = []
max_dy = []
max_DR = []

dy = np.linspace(0.1*DY, 1.5*DY, NUM_SIM).tolist() # DEFINE MIN. AND MAX. VALUE FOR INITIAL DISPLACEMENT
duct = np.linspace(2.0, 8.5, NUM_SIM).tolist()  # DEFINE MIN. AND MAX. VALUE FOR STRUCURAL DUCTILITY RATIO
osf = np.linspace(1.05, 1.35, NUM_SIM).tolist() # DEFINE MIN. AND MAX. VALUE FOR STRUCURAL OVER STRENGTH FACTOR
m = np.linspace(mi*(1 / NUM_SIM)*0.04, mi*0.04, NUM_SIM).tolist() # DEFINE MIN. AND MAX. VALUE FOR STRUCURAL MASS
IT = 0

# NUM_SIM is the number of simulations
for j in range(NUM_SIM):              # INTIAL DISPLACEMENT
    for i in range(NUM_SIM):          # STRCTURAL MASS
        for k in range(NUM_SIM):      # STRCTURAL DUCTILITY RATIO
            for l in range(NUM_SIM):  # STRCTURAL OVER-SRENGTH FACTOR
                data = ANALYSIS_SDOF(m[i], duct[k], osf[l], dy[j], IU=True)
                (time, displacement, velocity, acceleration,
                base_reaction, DI, PERIOD, damping_ratio, stiffness) = data
                # Calculate and store the max absolute values
                max_mass.append(m[i])
                max_dy.append(dy[j])
                max_time.append(np.max(np.abs(time)))
                max_displacement.append(np.max(np.abs(displacement)))
                max_velocity.append(np.max(np.abs(velocity)))
                max_acceleration.append(np.max(np.abs(acceleration)))
                max_base_reaction.append(np.max(np.abs(base_reaction)))
                max_K.append(np.max(np.abs(stiffness)))
                max_DR.append(damping_ratio)
                max_duct.append(duct[k])
                max_osf.append(osf[l])
                max_T.append(PERIOD)
                max_di.append(np.max(DI))
                # CHECK TO DAMAGE INDEX BOUNDARIES
                if max_di[-1] <= 0.0:
                    max_di[-1] = 0.0
                if max_di[-1] >= 100.0:
                    max_di[-1] = 100.0    
                IT += 1;
                print(f'\n SIMULATION {IT} DONE \n')    
else:
    print('\n\n Analysis completed successfully')    

totaltime = ti.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n') 

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
X = max_duct
Y = max_mass
Z = max_displacement
XLABEL = 'Structural Ductility Ratio [m/m]'   
YLABEL = 'Structural Mass [kg]' 
ZLABEL = "Max Displacement (PGD) [m]"
PLOT_3D_CONTOUR_XYZ(0, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = max_duct
Y = max_mass
Z = max_velocity
XLABEL = 'Structural Ductility Ratio [m/m]'   
YLABEL = 'Structural Mass [kg]' 
ZLABEL = "Max Velocity (PGV) [m/s]"
PLOT_3D_CONTOUR_XYZ(1, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = max_duct
Y = max_mass
Z = max_acceleration
XLABEL = 'Structural Ductility Ratio [m/m]'   
YLABEL = 'Structural Mass [kg]' 
ZLABEL = "Max Acceleration (PGA) [m/s²]"
PLOT_3D_CONTOUR_XYZ(2, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = max_duct
Y = max_mass
Z = max_base_reaction
XLABEL = 'Structural Ductility Ratio [m/m]'   
YLABEL = 'Structural Mass [kg]' 
ZLABEL = "Max Base-rection (PGA) [N]"
PLOT_3D_CONTOUR_XYZ(3, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = max_duct
Y = max_mass
Z = max_di
XLABEL = 'Structural Ductility Ratio [m/m]'   
YLABEL = 'Structural Mass [kg]' 
ZLABEL = "Sructural Ductility Damage Index (DI) [%]"
PLOT_3D_CONTOUR_XYZ(4, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = max_duct
Y = max_mass
Z = max_T
XLABEL = 'Structural Ductility Ratio [m/m]'   
YLABEL = 'Structural Mass [kg]' 
ZLABEL = "Sructural Period [s]"
PLOT_3D_CONTOUR_XYZ(5, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = max_duct
Y = max_mass
Z = max_DR
XLABEL = 'Structural Ductility Ratio [m/m]'   
YLABEL = 'Structural Mass [kg]' 
ZLABEL = "Sructural Damping Ratio [%]"
PLOT_3D_CONTOUR_XYZ(6, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = max_duct
Y = max_T
Z = max_DR
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
# --------------
# Visualization
# --------------
plt.figure(1, figsize=(10, 6))
# Response plot
plt.plot(time, acceleration, lw=1, color='black')
plt.xlabel('Time (s)')
plt.ylabel('Acceleration (g)')
plt.title(f'Last Analysis Structural Response ::: MAX. ABS. : {np.max(np.abs(acceleration)):.4f}')
plt.grid(True)
plt.show()    

im_values = max_acceleration 
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
im_values = max_di # Structural Ductility Damage Index
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

    'Damping_Ratio': max_DR,
    'Max_displacement': max_displacement,
    'Max_velocity': max_velocity,
    'Max_acceleration': max_acceleration,
    'Max_Base_Reaction': max_base_reaction,
    'Ductility_Damage_Index': max_di,
    'Period': max_T,
}
# Convert to DataFrame
results_df = pd.DataFrame(DATA_TOTAL)
# Export the DataFrame to an Excel file
results_df.to_excel('SDOF_SENSITIVITY_FREE_VIBRATION_DUCT_OSF_DATA_RESULTS.xlsx', index=False)
#%%------------------------------------------------------------------------------------------------  
XLABEL = 'Displacement'
YLABEL = 'Base Reaction'
TITLE = f'{YLABEL} and {XLABEL} scatter chart'
COLOR = 'orange'
X = max_displacement
Y = max_base_reaction
S01.PLOT_SCATTER(X, Y , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

# CLUSTER DATA
S01.CLUSTER_DATA(X, Y, XLABEL, YLABEL, MAX_CLUSTERS=3)
#%%------------------------------------------------------------------------------------------------
XLABEL = 'Velocity'
YLABEL = 'Base Reaction'
TITLE = f'{YLABEL} and {XLABEL} scatter chart'
COLOR = 'cyan'
X = max_velocity
Y = max_base_reaction
S01.PLOT_SCATTER(X, Y , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

# CLUSTER DATA
S01.CLUSTER_DATA(X, Y, XLABEL, YLABEL, MAX_CLUSTERS=3)
#%%------------------------------------------------------------------------------------------------
XLABEL = 'Acceleration'
YLABEL = 'Base Reaction'
TITLE = f'{YLABEL} and {XLABEL} scatter chart'
COLOR = 'lime'
X = max_acceleration
Y = max_base_reaction
S01.PLOT_SCATTER(X, Y , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

# CLUSTER DATA
S01.CLUSTER_DATA(X, Y, XLABEL, YLABEL, MAX_CLUSTERS=3)
#%%------------------------------------------------------------------------------------------------
XLABEL = 'Displacement'
YLABEL = 'Structural Ductility Damage Index (%)'
TITLE = f'{YLABEL} and {XLABEL} scatter chart'
COLOR = 'purple'
X = max_displacement
Y = max_di
S01.PLOT_SCATTER(X, Y , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

# CLUSTER DATA
S01.CLUSTER_DATA(X, Y, XLABEL, YLABEL, MAX_CLUSTERS=3)
#%%------------------------------------------------------------------------------------------------
# PLOT THE TIME-HISTORY
S01.PLOT_TIME_HISTORY(time, displacement, velocity, acceleration, base_reaction)
#%%------------------------------------------------------------------------------------------------
# RANDOM FOREST ANALYSIS
"""
This code predicts the free-vibration safety of a structure using simulation data by training a Random Forest Classifier to
 classify whether the system is "safe" or "unsafe" based on features like maximum displacement, velocity, acceleration,
 and base reaction. A regression model is also trained to estimate safety likelihood. It evaluates model performance using
 metrics like classification accuracy, mean squared error, and R² score. Additionally, it identifies key features influencing
 safety through feature importance analysis. The tool aids in free-vibration risk assessment, structural optimization, and understanding
 critical safety parameters.
"""

data = {
    'Max_displacement': max_displacement,
    'Max_velocity': max_velocity,
    'Max_acceleration': max_acceleration,
    'Max_Base_Reaction': max_base_reaction,
    'Ductility_Damage_Index': max_di,
    'Damping_Ratio': max_DR,
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
num_sim = NUM_SIM
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