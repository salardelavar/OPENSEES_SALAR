###########################################################################################################
#                    >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                     #
#          CONSTRAINED OPTIMIZATION OF REINFORCED CONCRETE COLUMN DESIGN USING NONLINEAR PUSHOVER         #
#                                             ANALYSIS IN OPENSEES                                        #
#---------------------------------------------------------------------------------------------------------#
#     FIND BEST COLUMN REBAR DIAMETER AND COLUMN SECTION DEPTH WITH TARGET STRUCTURAL DUCTILITY RATIO     #
#---------------------------------------------------------------------------------------------------------#
#                          THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                     #
#                                   EMAIL: salar.d.ghashghaei@gmail.com                                   #
###########################################################################################################
"""
# 1. Import OpenSeesPy, NumPy, SciPy, and custom analysis modules
#    Used for nonlinear FEM analysis, optimization, and post-processing

# 2. Define concrete and steel material properties
#    Includes confined/unconfined concrete and nonlinear reinforcing steel

# 3. Define column and beam geometry (section sizes, cover, rebar diameter)
#    These parameters control stiffness, strength, and ductility

# 4. Define pushover analysis parameters
#    Maximum displacement, increment size, element lengths, and solver settings

# 5. Define PUSHOVER_ANALYSIS() function
#    Builds a 2D RC frame model, applies loads, and performs nonlinear pushover

# 6. Inside PUSHOVER_ANALYSIS():
#    - Create nodes, boundary conditions, sections, and elements
#    - Apply displacement-controlled static analysis
#    - Record forces, displacements, rotations, and stiffness values

# 7. Fit a bilinear curve to base-shear vs displacement
#    Used to compute ductility ratio (μ) and over-strength factor (Ω₀)

# 8. Return structural response data including μ and Ω₀
#    These are the key performance indicators for optimization

# 9. Define target ranges for ductility and over-strength
#    Ensures acceptable seismic performance

# 10. Define design variable bounds
#     Limits rebar diameter and column depth to practical values

# 11. Define a penalty function
#     Penalizes solutions violating ductility or over-strength limits

# 12. Define the objective function
#     Runs pushover analysis and measures deviation from target performance

# 13. Add constraint penalties to the objective
#     Forces the optimizer to stay within performance limits

# 14. Print intermediate optimization results
#     Useful for monitoring convergence and debugging

# 15. Set an initial design guess
#     Starting point for the optimizer

# 16. Use SciPy’s minimize() with L-BFGS-B algorithm
#     A stable, bounded optimizer suitable for expensive nonlinear analyses

# 17. Run the optimization loop automatically
#     Adjusts DIAc and Hc to minimize objective + penalties

# 18. Store optimization results in 'result'
#     Includes optimal variables, iterations, and convergence status

# 19. Print final optimal column rebar diameter and section depth
#     Represents the best seismic-performance design found

# 20. Overall goal:
#     Automatically design an RC column section that satisfies
#     target ductility and over-strength using nonlinear pushover analysis

"""
import openseespy.opensees as ops
import matplotlib.pyplot as plt
import numpy as np
import ANALYSIS_FUNCTION as S02
import CONCRETE_SECTION_FUN as S03
import PLOT_2D as S04
import BILINEAR_CURVE as BC
import time as TI

CONCRETE_KIND = 1
# Define materials for nonlinear columns and beam
# Define parameters (units: mm, N)
# ------------------------------------------
# CONCRETE                  tag   f'c        ec0   f'cu        ecu
if CONCRETE_KIND == 1:
    # Cover concrete (unconfined)
    fcU = -18           # [N/mm²] Concrete Compressive Strength
    ec0U = -0.0025      # [mm/mm] Concrete Compressive Strain
    fcUU = -2           # [N/mm²] Concrete Compressive Ultimate Strength
    ecuU = -0.008       # [mm/mm] Concrete Compressive Ultimate Strain
    # Core concrete (confined)
    Kfc = 1.5;			# ratio of confined to unconfined concrete strength - COLUMN
    fcC = Kfc*fcU;      # [N/mm²] Concrete Compressive Strength
    ec0C = -0.0045      # [mm/mm] Concrete Compressive Strain
    fcUC = -21          # [N/mm²] Concrete Compressive Ultimate Strength
    ecuC = -0.015       # [mm/mm] Concrete Compressive Ultimate Strain
    
if CONCRETE_KIND == 2:
    # Cover concrete (unconfined)
    fcU = -18               # [N/mm²] Concrete Compressive Strength
    fc1U = fcU;			    # [N/mm²] UNCONFINED concrete (todeschini parabolic model), maximum stress
    eps1U = -0.0025;		# [mm/mm] strain at maximum strength of unconfined concrete
    fc2U = 0.2*fc1U;		# [N/mm²] ultimate stress
    eps2U = -0.012;			# [mm/mm] strain at ultimate stress
    Lambda = 0.1;			# ratio between unloading slope at $eps2 and initial slope $Ec
    # Core concrete (confined)
    Kfc = 1.5;			# ratio of confined to unconfined concrete strength - COLUMN
    fc1C = Kfc*fcU;		# [N/mm²] CONFINED concrete (mander model), maximum stress - COLUMN
    Ec = 4700 * np.sqrt(-fcU) # [N/mm^2] Concrete Elastic Modulus
    eps1C = 2*fc1C/Ec;	# [mm/mm] strain at maximum stress 
    fc2C = 0.2*fc1C;	# [N/mm²] ultimate stress
    eps2C = 5*eps1C;	# [mm/mm] strain at ultimate stress 
    
    # tensile-strength properties
    ftC = -0.55*fc1C;		# [N/mm²] tensile strength +tension
    ftU = -0.55*fc1U;		# [N/mm²] tensile strength +tension
    
    # tensile-strength properties
    ftC = -0.55*fc1C;		# [N/mm²] tensile strength +tension
    ftU = -0.55*fc1U;		# [N/mm²] tensile strength +tension
    Ets = ftU/0.002;		# tension softening stiffness
     
# STEEL
# Reinforcing steel
fy = 400          # [N/mm²] Steel Rebar Yield Strength   
Es = 2e5          # [N/mm²] Modulus of Elasticity
ey = fy/Es        # [mm/mm] Steel Rebar Yield Strain
fu = 1.1818*fy    # [N/mm²] Steel Rebar Ultimate Strength
esu = ey*75.2     # [mm/mm] Steel Rebar Ultimate Strain
Esh = (fu - fy)/(esu - ey)
Bs = Esh / Es

# Column Section
Bc = 500                 # [mm] Depth of the Section 
Hc = 500                 # [mm] Height of the Section  
coverC = 50              # [mm] Concrete Section Cover
DIAc = 25                # [mm] # Rebar Size Diameter


# Beam Section
Bb = 500                 # [mm] Depth of the Section 
Hb = 300                 # [mm] Height of the Section  
coverB = 50              # [mm] Concrete Section Cover


DMAX = 175               # [mm] Maximum Lateral Displacement in X Dir.
DINCR = 0.01             # [mm] Increment Displacement


LENGTH_COL = 3000        # [mm] Column Length 
LENGTH_BM = 7000         # [mm] Beam Length 
# Define Analysis Properties
MAX_ITERATIONS = 5000    # Convergence iteration for test
MAX_TOLERANCE = 1.0e-6   # Convergence tolerance for test
#STEEL_KIND: 1 -> WITHOUT HARDENING AND ULTIMATE STRAIN
#STEEL_KIND: 2 -> WITH HARDENING AND ULTIMATE STRAIN
#CONCRETE_KIND: 1 -> Concrete01
#CONCRETE_KIND: 2 -> Concrete02
#%%------------------------------------------------------------------------------
def PUSHOVER_ANALYSIS(DIAc, Hc, LENGTH_COL, LENGTH_BM, DMAX, DINCR, STEEL_KIND, CONCRETE_KIND):
    #DIAc = 25               # [mm] Column Rebar Size Diameter
    AsC = np.pi*(DIAc**2)/4  # [mm²] Column Area of Rebar
    DIAb = DIAc-5            # [mm] Beam Rebar Size Diameter
    AsB = np.pi*(DIAb**2)/4  # [mm²] Beam Area of Rebar
    # Reset model
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)

    # Nodes
    ops.node(1, 0.0, 0.0)
    ops.node(2, LENGTH_BM, 0.0)
    ops.node(3, 0.0, LENGTH_COL)
    ops.node(4, LENGTH_BM, LENGTH_COL)
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 1, 1, 1)

    secTagC = 10
    secTagB = 20
    coreTag = 1
    coverTag = 2
    steelTag = 3
    
    if STEEL_KIND == 1:# WITHOUT HARDENING AND ULTIMATE STRAIN
        ops.uniaxialMaterial('Steel01', steelTag, fy, Es, 0.0) 
    if STEEL_KIND == 2:# WITH HARDENING AND ULTIMATE STRAIN    
        pinchX = 0.8   # Pinching factor in X direction
        pinchY = 0.5   # Pinching factor in Y direction
        damage1 = 0.0  # Damage due to ductility
        damage2 = 0.0  # Damage due to energy
        beta = 0.1 # Stiffness degradation parameter
        ops.uniaxialMaterial('Hysteretic', steelTag, fy, ey, fu, esu, 0.2*fu, 1.1*esu, -fy, -ey, -fu, -esu, -0.2*fu, -1.1*esu, pinchX, pinchY, damage1, damage2, beta)
        # INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Hysteretic_Material

    if CONCRETE_KIND == 1:
        ops.uniaxialMaterial('Concrete01', coreTag, fcC, ec0C, fcUC, ecuC)  # Core concrete (confined)
        ops.uniaxialMaterial('Concrete01', coverTag, fcU, ec0U, fcUU, ecuU) # Cover concrete (unconfined)
    if CONCRETE_KIND == 2:
        ops.uniaxialMaterial('Concrete02', coreTag, fc1C, eps1C, fc2C, eps2C, Lambda, ftC, Ets) # Core concrete (confined)
        ops.uniaxialMaterial('Concrete02', coverTag, fc1U, eps1U, fc2U, eps2U, Lambda, ftU, Ets) # Cover concrete (unconfined)
    
    # COLUMN SECTION
    S03.CONFINED_CONCRETE_SECTION(secTagC, Hc, Bc, coverC, AsC, coreTag, coverTag, steelTag, COL=True)
    # BEAM SECTION
    S03.CONFINED_CONCRETE_SECTION(secTagB, Hb, Bb, coverB, AsB, coreTag, coverTag, steelTag, COL=False)
    
    # Define geometric transformation (corotational for large displacements)
    transfTag = 1
    #ops.geomTransf('Linear', transfTag)
    #ops.geomTransf('PDelta', transfTag)
    ops.geomTransf('Corotational', transfTag)
    numIntgrPts = 5
    ops.element('nonlinearBeamColumn', 1, 1, 3, numIntgrPts, secTagC, transfTag) # COLUMN 01
    ops.element('nonlinearBeamColumn', 2, 2, 4, numIntgrPts, secTagC, transfTag) # COLUMN 02
    ops.element('nonlinearBeamColumn', 3, 3, 4, numIntgrPts, secTagB, transfTag) # BEAM 01

    # Data storage
    FORCE_S, FORCE_A, MOMENT = [], [], []
    DISP_X, DISP_Y, ROT = [], [], []
    KA, KS, KI, STEP = [], [], [], []

    # Define time series and load pattern
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(3, 1.0, -1.0, 0.0)
    ops.load(4, 1.0, -1.0, 0.0)

    # Total steps per half-cycle
    steps = int(np.abs(DMAX)/np.abs(DINCR))

    # Use displacement control on rotational dof (dof 3 at node 2)
    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.system('BandGeneral')
    ops.algorithm('Newton')
    ops.analysis('Static')
    
    for step in range(steps):
        
        ops.integrator('DisplacementControl', 3, 1, DINCR) 
        ops.integrator('DisplacementControl', 4, 1, DINCR) 
        OK = ops.analyze(1)
        S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS
        # Record results
        ops.reactions()
        S = ops.nodeReaction(1, 1) + ops.nodeReaction(2, 1) # SHEAR BASE REACTION
        A = ops.nodeReaction(1, 2) + ops.nodeReaction(2, 2) # AXIAL BASE REACTION
        M = ops.nodeReaction(1, 3) + ops.nodeReaction(2, 3) # MOMENT BASE REACTION
        #print(rot, M)
        disp_X = ops.nodeDisp(3, 1) # LATERAL DISPLACEMENT IN X FOR NODE 3
        disp_Y = ops.nodeDisp(3, 2) # LATERAL DISPLACEMENT IN Y FOR NODE 3
        rot = ops.nodeDisp(3, 3)    # ROTATION IN Z FOR NODE 3
        FORCE_S.append(S)
        FORCE_A.append(A)
        MOMENT.append(M)
        DISP_X.append(disp_X)
        DISP_Y.append(disp_Y)
        ROT.append(rot)
        KS.append(np.abs(S)/np.abs(disp_X)) # LATERAL STIFFNESS IN X
        KA.append(np.abs(A)/np.abs(disp_Y)) # LATERAL STIFFNESS IN Y
        KI.append(np.abs(M)/np.abs(rot))    # ROTATIONAL STIFFNESS IN Z
        STEP.append(step)
        #print(step+1, disp_X, S)
        
    #ops.wipe() 
    # ---------------------------------------
    #  Plot BaseShear-Displacement Analysis
    # ---------------------------------------
    XX = np.abs(DISP_X); YY = np.abs(FORCE_S); # ABSOLUTE VALUE
    SLOPE_NODE = 10

    DATA = BC.BILNEAR_CURVE(XX, YY, SLOPE_NODE)
    X, Y, Elastic_ST, Plastic_ST, Tangent_ST, Ductility_Rito, Over_Strength_Factor = DATA
    """
    XLABEL = 'Displacement in X [mm]'
    YLABEL = 'Base-Shear [N]'
    LEGEND01 = 'Curve'
    LEGEND02 = 'Bilinear Fitted'
    LEGEND03 = 'Undefined'
    TITLE = f'Last Data of BaseShear-Displacement Analysis - Ductility Ratio: {X[2]/X[1]:.4f} - Over Strength Factor: {Y[2]/Y[1]:.4f}'
    COLOR = 'black'
    BC.PLOT_2D(np.abs(DISP_X), np.abs(FORCE_S), X, Y, X, Y, XLABEL, YLABEL, TITLE, LEGEND01, LEGEND02, LEGEND03, COLOR='black', Z=2) 
    print(f'\t\t Ductility Ratio: {Y[2]/Y[1]:.4f}')
    """

    # Calculate Over Strength Coefficient (Ω0)
    Omega_0 = Y[2] / Y[1]
    # Calculate Displacement Ductility Ratio (μ)
    mu = X[2] / X[1]
    # Calculate Ductility Coefficient (Rμ)
    #R_mu = (2 * mu - 1) ** 0.5
    #R_mu = 1
    R_mu = mu
    # Calculate Structural Behavior Coefficient (R)
    R = Omega_0 * R_mu
    """
    print(f'Over Strength Coefficient (Ω0):      {Omega_0:.4f}')
    print(f'Displacement Ductility Ratio (μ):    {mu:.4f}')
    print(f'Ductility Coefficient (Rμ):          {R_mu:.4f}')
    print(f'Structural Behavior Coefficient (R): {R:.4f}') 
    """

    return FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y, ROT, KA, KS, KI, mu, Omega_0, STEP

#%%------------------------------------------------------------------------------
# FIND BEST COLUMN REBAR DIAMETER AND COLUMN SECTION DEPTH WITH TARGET STRUCTURAL DUCTILITY RATIO:
from scipy.optimize import minimize
X = np.zeros([2,1])    
X[0] = DIAc        # [mm] Intial Guess Column Rebar Size Diameter  
X[1] = Hc          # [mm] Intial Guess Column Height
ESP = 1e-5         # Finite difference derivative Convergence Tolerance
TOLERANCE = 1e-6   # Convergence Tolerance
RESIDUAL = 100     # Convergence Residual 
damping = 0.8      # Damping factor for updates
reg_factor = 1e-6   # Regularization to avoid singular Jacobian
IT = 0             # Intial Iteration
ITMAX = 100000     # Max. Iteration

TARGET_DUCT = 7.5  # [mm/mm] Target Structural Behavior Coefficient
TARGET_OSF = 1.1   # [N/N] Over Srength Factor

# TARGETS
DUCT_MIN, DUCT_MAX = 7.0, 7.5
OSF_MIN,  OSF_MAX  = 1.0, 1.1

# Weights (tune if needed)
W_DUCT = 1.0
W_OSF  = 1.0

# ANALYSIS BOUNDS
bounds = [
    (8.0, 35.0),     # [mm] DIAc - Column Rebar Size Diameter 
    (200.0, 700.0)   # [mm] Hc - Column Height  
]

# PENALTY FUNCTION
def penalty(val, vmin, vmax):
    if val < vmin:
        return (vmin - val)**2
    elif val > vmax:
        return (val - vmax)**2
    else:
        return 0.0

# OBJECTIVE FUNCTION
def objective(X):
    DIAc, Hc = X

    try:
        DATA = PUSHOVER_ANALYSIS(
            DIAc, Hc,
            LENGTH_COL, LENGTH_BM,
            DMAX, DINCR,
            STEEL_KIND=2,
            CONCRETE_KIND=1
        )

        (_, _, _, _, _, _, _, _, _,
         DUCT_ANA, OSF_ANA, _) = DATA

    except:
        # Analysis failed → strong penalty
        return 1e6

    # Main objective: match target values
    obj = (
        (DUCT_ANA - TARGET_DUCT)**2 +
        (OSF_ANA  - TARGET_OSF )**2
    )

    # Constraint penalties
    pen = (
        penalty(DUCT_ANA, DUCT_MIN, DUCT_MAX) +
        penalty(OSF_ANA,  OSF_MIN,  OSF_MAX)
    )

    print(f"DIAc={DIAc:.2f} mm | Hc={Hc:.1f} mm | "
          f"DUCT={DUCT_ANA:.5f} | OSF={OSF_ANA:.5f} | Obj={pen:.3e}")
    return obj + 100.0 * pen  # penalty weight

# INITIAL GUESS
X0 = np.array([20.0, 400.0])

# OPTIMIZATION
starttime = TI.process_time()

result = minimize(
    objective,
    X0,                 # Initial guess for the design variables
    method='L-BFGS-B',  # Selects the optimization algorithm
    bounds=bounds,
    options={
        'maxiter': 200, # Maximum number of optimization iterations
        'ftol': 1e-5,   # Function tolerance (convergence criterion)
        'disp': True    # Print optimization progress to the console
    }
)

endtime = TI.process_time()

# RESULTS
print("\n\n\n=== OPTIMIZATION COMPLETED ===")
print(f"Optimal DIAc = {result.x[0]:.3f} mm")
print(f"Optimal Hc   = {result.x[1]:.3f} mm")
print(f"Final Objective = {result.fun:.3e}")
print(f"Iterations = {result.nit}")
print("Success:", result.success)
print("Message:", result.message)
print("Duration (sec):", TI.process_time() - starttime)

#%%------------------------------------------------------------------------------
# Run Final Pushover at Optimal Analysis for Final Evaluation
DIAc_opt, Hc_opt = result.x

DATA_FINAL = PUSHOVER_ANALYSIS(
    DIAc_opt, Hc_opt,
    LENGTH_COL, LENGTH_BM,
    DMAX, DINCR,
    STEEL_KIND=2,
    CONCRETE_KIND=1
)

(FORCE_S, FORCE_A, MOMENT,
 DISP_X, DISP_Y, ROT,
 KA, KS, KI,
 DUCT_ANA, OSF_ANA, STEP) = DATA_FINAL

#%%------------------------------------------------------------------------------
plt.figure(1, figsize=(12, 8))
plt.plot(MOMENT, FORCE_A, color='black')
#plt.scatter(MOMENT, FORCE_A, color='blue', linewidth=2)
plt.title('P-M Interaction')
plt.ylabel('Axial Force [N]')
plt.xlabel('Bending Moment [N.mm]')
#plt.legend(['Steel01: WITHOUT HARDENING AND ULTIMATE STRAIN', 'Hysteretic: WITH HARDENING AND ULTIMATE STRAIN'])
plt.grid()
plt.show()

plt.figure(2, figsize=(12, 8))
plt.plot(DISP_X, FORCE_S, color='green', linewidth=2)
#plt.scatter(DISP_X, FORCE_S, color='purple', linewidth=2)
plt.title('SHEAR FORCE-DISPLACEMENT DIAGRAM')
plt.ylabel('Shear Force [N]')
plt.xlabel('Displacement  in X [mm]')
plt.grid()
plt.show()

plt.figure(3, figsize=(12, 8))
plt.plot(DISP_Y, FORCE_A, color='purple', linewidth=2)
#plt.scatter(DISP_Y, FORCE_A, color='purple', linewidth=2)
plt.title('AXIAL FORCE-DISPLACEMENT DIAGRAM')
plt.ylabel('Axial Force [N]')
plt.xlabel('Displacement in Y [mm]')
plt.grid()
plt.show()

plt.figure(4, figsize=(12, 8))
plt.plot(ROT, MOMENT, color='red', linewidth=2)
#plt.scatter(ROT, MOMENT, color='red', linewidth=2)
plt.title('MOMENT-ROTATION DIAGRAM')
plt.ylabel('Moment [kN.mm]')
plt.xlabel('Rotation [rad]')
plt.grid()
plt.show()

plt.figure(5, figsize=(12, 8))
plt.plot(KI, KS, color='black', linewidth=2)
#plt.scatter(KI, KS, color='black', linewidth=2)
plt.title('ROTATIONAL STIFFNESS-LATERAL STIFFNESS DIAGRAM')
plt.ylabel('Rotational Stiffness [N.mm/Rad]')
plt.xlabel('Lateral Stiffness in X Dir. [N/mm]')
plt.semilogx()
plt.semilogy()
plt.grid()
plt.show()

plt.figure(6, figsize=(12, 8))
plt.plot(KI, KA, color='black', linewidth=2)
#plt.scatter(KI, KA, color='black', linewidth=2)
plt.title('ROTATIONAL STIFFNESS-LATERAL STIFFNESS DIAGRAM')
plt.ylabel('Rotational Stiffness [N.mm/Rad]')
plt.xlabel('Lateral Stiffness in Y Dir. [N/mm]')
plt.semilogx()
plt.semilogy()
plt.grid()
plt.show()

plt.figure(7, figsize=(12, 8))
plt.plot(STEP, FORCE_A, color='brown', linewidth=2)
#plt.scatter(STEP, FORCE_A, color='brown', linewidth=2)
plt.title('Axial Force During the Analysis')
plt.ylabel('Axial Force [N]')
plt.xlabel('Steps')
plt.grid()
plt.show()

plt.figure(8, figsize=(12, 8))
plt.plot(STEP, FORCE_S, color='purple', linewidth=2)
#plt.scatter(STEP, FORCE_S, color='purple', linewidth=2)
plt.title('Shear Force During the Analysis')
plt.ylabel('Shear Force [N]')
plt.xlabel('Steps')
plt.grid()
plt.show()

plt.figure(9, figsize=(12, 8))
plt.plot(STEP, MOMENT, color='green', linewidth=2)
#plt.scatter(STEP, MOMENT, color='green', linewidth=2)
plt.title('Moment During the Analysis')
plt.ylabel('Moment [kN.mm]')
plt.xlabel('Steps')
plt.grid()
plt.show()

plt.figure(10, figsize=(12, 8))
plt.plot(STEP, DISP_X, color='brown', linewidth=2)
#plt.scatter(STEP, DISP_X, color='brown', linewidth=2)
plt.title('Displacement During the Analysis')
plt.ylabel('Displacement - X [mm]')
plt.xlabel('Steps')
plt.grid()
plt.show()

plt.figure(11, figsize=(12, 8))
plt.plot(STEP, DISP_Y, color='blue', linewidth=2)
#plt.scatter(STEP, DISP_X, color='brown', linewidth=2)
plt.title('Displacement During the Analysis')
plt.ylabel('Displacement - Y [mm]')
plt.xlabel('Steps')
plt.grid()
plt.show()

plt.figure(12, figsize=(12, 8))
plt.plot(STEP, ROT, color='black', linewidth=2)
#plt.scatter(STEP, ROT, color='black', linewidth=2)
plt.title('Rotation During the Analysis')
plt.ylabel('Rotation [rad]')
plt.xlabel('Steps')
plt.grid()
plt.show()
#%%------------------------------------------------------------------------------
# ---------------------------------------
#  Plot BaseShear-Displacement Analysis
# ---------------------------------------
XX = np.abs(DISP_X); YY = np.abs(FORCE_S); # ABSOLUTE VALUE
SLOPE_NODE = 10

DATA = BC.BILNEAR_CURVE(XX, YY, SLOPE_NODE)
X, Y, Elastic_ST, Plastic_ST, Tangent_ST, Ductility_Rito, Over_Strength_Factor = DATA

XLABEL = 'Displacement in X [mm]'
YLABEL = 'Base-Shear [N]'
LEGEND01 = 'Curve'
LEGEND02 = 'Bilinear Fitted'
LEGEND03 = 'Undefined'
TITLE = f'Last Data of BaseShear-Displacement Analysis - Ductility Ratio: {X[2]/X[1]:.4f} - Over Strength Factor: {Y[2]/Y[1]:.4f}'
COLOR = 'black'
BC.PLOT_2D(np.abs(DISP_X), np.abs(FORCE_S), X, Y, X, Y, XLABEL, YLABEL, TITLE, LEGEND01, LEGEND02, LEGEND03, COLOR='black', Z=2) 
#print(f'\t\t Ductility Ratio: {Y[2]/Y[1]:.4f}')


# Calculate Over Strength Coefficient (Ω0)
Omega_0 = Y[2] / Y[1]
# Calculate Displacement Ductility Ratio (μ)
mu = X[2] / X[1]
# Calculate Ductility Coefficient (Rμ)
#R_mu = (2 * mu - 1) ** 0.5
#R_mu = 1
R_mu = mu
# Calculate Structural Behavior Coefficient (R)
R = Omega_0 * R_mu
print(f'Over Strength Coefficient (Ω0):  {Omega_0:.4f}')
print(f'Displacement Ductility Ratio (μ):{mu:.4f}')
print(f'Ductility Coefficient (Rμ):  {R_mu:.4f}')
print(f'Structural Behavior Coefficient (R): {R:.4f}') 

#%%------------------------------------------------------------------------------  
# %% Plot 2D Frame Shapes
S04.PLOT_2D_FRAME(deformed_scale=10)  # Adjust scale factor as needed
#%%------------------------------------------------------------------------------   
# EXPORT DATA TO EXCEL
import pandas as pd
DATA_TOTAL = { 
    'DISP_X_WO': DISP_X,
    'DISP_Y_WO': DISP_Y,
    'ROTATION_WO': ROT,
    'AXIAL_FORCE_WO': FORCE_A,
    'SHEAR_FORCE_WO': FORCE_S,
    'MOMENT_WO': MOMENT,
    'AXIAL_RIGIDITY_WO': np.abs(FORCE_A),
    'ROTATIONAL_ST_WO': KI,
}
# Convert to DataFrame
results_df = pd.DataFrame(DATA_TOTAL)
# Export the DataFrame to an Excel file
results_df.to_excel('PUSHOVER_Cdepth_REBAR_DUCT_OPTIMIZATION_RESULTS.xlsx', index=False) 
#%%------------------------------------------------------------------------------

    
