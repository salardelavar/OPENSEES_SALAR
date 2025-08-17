###########################################################################################################
#                   >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                      #
# OPTIMIZATION OF STRUCTURAL BEHAVIOR COEFFICIENT USING PUSHOVER ANALYSIS OF CONCRETE FRAME SECTIONS:     #
# EVALUATING STRAIN HARDENING AND ULTIMATE STRAIN EFFECTS IN OPENSEES. DETERMINING OPTIMAL COLUMN SECTION #
# REBAR DIAMETER FOR A TARGET STRUCTURAL BEHAVIOR COEFFICIENT VIA THE NEWTON-RAPHSON METHOD.              #
#---------------------------------------------------------------------------------------------------------#
#                              OPTIMIZATION ALGORITHM: NEWTON-RAPHSON METHOD                              #
#---------------------------------------------------------------------------------------------------------#
#                          THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                     #
#                                   EMAIL: salar.d.ghashghaei@gmail.com                                   #
###########################################################################################################
"""
1. The script performs pushover analysis on a concrete frame using OpenSees
 to optimize the column rebar diameter for a target ductility ratio.
2. Two steel material models (*Steel01* and *Hysteretic*) and two concrete
 models (*Concrete01* and *Concrete02*) are supported.
3. A frame with beam and column elements is created, and nonlinear beam-column
 elements are used for realistic simulation.
4. Rebar areas are calculated based on input diameters, and sectional properties
 are defined using confined and unconfined concrete.
5. The *PUSHOVER\_ANALYSIS* function incrementally applies lateral displacement
 and records force, displacement, and stiffness data.
6. The response is processed to compute the bilinear approximation and extract
 ductility and strength parameters.
7. A Newton-Raphson root-finding algorithm adjusts the column rebar diameter to
 match the target structural ductility ratio.
8. Finite difference approximation is used to estimate the derivative of the
 ductility function with respect to rebar diameter.
9. Each iteration updates the rebar size until convergence is achieved or the
 maximum number of iterations is reached.
10. Convergence is based on the residual of the diameter update (DX) relative
 to a tolerance threshold.
11. The optimal column and beam rebar diameters are printed upon successful convergence.
12. This method allows automated rebar design optimization based on seismic
 performance criteria.
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
fy = 4000          # [N/mm²] Steel Rebar Yield Strength   
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

#%% DEFINE PARAMETERS FOR NONLINEAR STATIC ANALYSIS 
DMAX = 400.0             # [mm] Maximum Lateral Displacement in X Dir.
DINCR = 0.05             # [mm] Increment Displacement

#%% DEFINE THE ELEMENTS LENGTH
LENGTH_COL = 3000        # [mm] Column Length 
LENGTH_BM = 7000         # [mm] Beam Length 

#%% DEFINE ANALYSIS PROPERTIES
MAX_ITERATIONS = 5000    # Convergence iteration for test
MAX_TOLERANCE = 1.0e-6   # Convergence tolerance for test
#STEEL_KIND: 1 -> WITHOUT HARDENING AND ULTIMATE STRAIN
#STEEL_KIND: 2 -> WITH HARDENING AND ULTIMATE STRAIN
#CONCRETE_KIND: 1 -> Concrete01
#CONCRETE_KIND: 2 -> Concrete02
#%%------------------------------------------------------------------------------
def PUSHOVER_ANALYSIS(DIAc, LENGTH_COL, LENGTH_BM, DMAX, DINCR, STEEL_KIND, CONCRETE_KIND):
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
        KS.append(np.abs(S/disp_X)) # LATERAL STIFFNESS IN X
        KA.append(np.abs(A/disp_Y)) # LATERAL STIFFNESS IN Y
        KI.append(np.abs(M/rot))    # ROTATIONAL STIFFNESS IN Z
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

    return FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y, ROT, KA, KS, KI, R, STEP

#%%------------------------------------------------------------------------------
# WITHOUT HARDENING AND ULTIMATE STRAIN
#FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y, ROT, KA, KS, KI, R, STEP = PUSHOVER_ANALYSIS(DIAc, LENGTH_COL, LENGTH_BM, DMAX, DINCR, STEEL_KIND=1, CONCRETE_KIND=1)
# WITH HARDENING AND ULTIMATE STRAIN
#FORCE_S02, FORCE_A02, MOMENT02, DISP_X02, DISP_Y02, ROT02, KA02, KS02, KI02, R02, STEP02  = PUSHOVER_ANALYSIS(DIAc, LENGTH_COL, LENGTH_BM, DMAX, DINCR, STEEL_KIND=2, CONCRETE_KIND=1)
#%%------------------------------------------------------------------------------
# FIND BEST COLUMN REBAR DIAMETER WITH STRUCTURAL BEHAVIOR COEFFICIENT OPTIMIZATION:
    
X = DIAc           # [mm] Intial Guess Column Rebar Size Diameter  
ESP = 1e-5         # Finite difference derivative Convergence Tolerance
TOLERANCE = 1e-10  # Convergence Tolerance
RESIDUAL = 100     # Convergence Residual 
IT = 0             # Intial Iteration
ITMAX = 100000     # Max. Iteration
TARGET_R = 4.5     # Target Structural Behavior Coefficient

# Analysis Durations:
starttime = TI.process_time()

### FIND THE OPTIMUM VALUE 
while (RESIDUAL > TOLERANCE):
    # X -------------------------------------------------------
    DATA  = PUSHOVER_ANALYSIS(X, LENGTH_COL, LENGTH_BM, DMAX, DINCR, STEEL_KIND=2, CONCRETE_KIND=1)
    FORCE_S02, FORCE_A02, MOMENT02, DISP_X02, DISP_Y02, ROT02, KA02, KS02, KI02, R_ANA, STEP02 = DATA
    print(f'R ANALYSIS: {R_ANA:.8f}')
    F = R_ANA - TARGET_R
    print('F: ', F)
    # Xmin -------------------------------------------------------
    XMIN = X - ESP  
    DATA  = PUSHOVER_ANALYSIS(XMIN, LENGTH_COL, LENGTH_BM, DMAX, DINCR, STEEL_KIND=2, CONCRETE_KIND=1)
    FORCE_S02, FORCE_A02, MOMENT02, DISP_X02, DISP_Y02, ROT02, KA02, KS02, KI02, R_ANA_MIN, STEP02 = DATA
    Fmin = R_ANA_MIN - TARGET_R
    print('Fmin: ', Fmin)
    # Xmax -------------------------------------------------------
    XMAX = X + ESP 
    DATA  = PUSHOVER_ANALYSIS(XMAX, LENGTH_COL, LENGTH_BM, DMAX, DINCR, STEEL_KIND=2, CONCRETE_KIND=1)
    FORCE_S02, FORCE_A02, MOMENT02, DISP_X02, DISP_Y02, ROT02, KA02, KS02, KI02, R_ANA_MAX, STEP02 = DATA
    Fmax = R_ANA_MAX - TARGET_R
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
        print(f'\t\t Optimum Column Rebar Size Diameter :      {X:.6f}')
        print(f'\t\t Optimum Beam Rebar Size Diameter :        {X-5:.6f}')
        print(f'\t\t Iteration Counts:                         {IT}')
        print(f'\t\t Convergence Residual:                     {RESIDUAL:.10e}')

    

totaltime = TI.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n')
#%%------------------------------------------------------------------------------
plt.figure(1, figsize=(12, 8))
#plt.plot(MOMENT, FORCE_A, color='black')
plt.plot(MOMENT02, FORCE_A02, color='cyan', linestyle='--', linewidth=2)
#plt.scatter(MOMENT, FORCE_A, color='blue', linewidth=2)
#plt.scatter(MOMENT02, FORCE_A02, color='cyan', linestyle='--', linewidth=2)
plt.title('P-M Interaction')
plt.ylabel('Axial Force [N]')
plt.xlabel('Bending Moment [N.mm]')
#plt.legend(['Steel01: WITHOUT HARDENING AND ULTIMATE STRAIN', 'Hysteretic: WITH HARDENING AND ULTIMATE STRAIN'])
plt.grid()
plt.show()

plt.figure(2, figsize=(12, 8))
#plt.plot(DISP_X, FORCE_S, color='green', linewidth=2)
plt.plot(DISP_X02, FORCE_S02, color='lime', linestyle='--', linewidth=2)
#plt.scatter(DISP_X, FORCE_S, color='purple', linewidth=2)
#plt.scatter(DISP_X02, FORCE_S02, color='magenta', linestyle='--', linewidth=2)
plt.title('SHEAR FORCE-DISPLACEMENT DIAGRAM')
plt.ylabel('Shear Force [N]')
plt.xlabel('Displacement  in X [mm]')
#plt.legend(['Steel01: WITHOUT HARDENING AND ULTIMATE STRAIN', 'Hysteretic: WITH HARDENING AND ULTIMATE STRAIN'])
plt.grid()
plt.show()

plt.figure(3, figsize=(12, 8))
#plt.plot(DISP_Y, FORCE_A, color='purple', linewidth=2)
plt.plot(DISP_Y02, FORCE_A02, color='magenta', linestyle='--', linewidth=2)
#plt.scatter(DISP_Y, FORCE_A, color='purple', linewidth=2)
#plt.scatter(DISP_Y02, FORCE_A02, color='magenta', linestyle='--', linewidth=2)
plt.title('AXIAL FORCE-DISPLACEMENT DIAGRAM')
plt.ylabel('Axial Force [N]')
plt.xlabel('Displacement in Y [mm]')
#plt.legend(['Steel01: WITHOUT HARDENING AND ULTIMATE STRAIN', 'Hysteretic: WITH HARDENING AND ULTIMATE STRAIN'])
plt.grid()
plt.show()

plt.figure(4, figsize=(12, 8))
#plt.plot(ROT, MOMENT, color='red', linewidth=2)
plt.plot(ROT02, MOMENT02, color='orange', linestyle='--', linewidth=2)
#plt.scatter(ROT, MOMENT, color='red', linewidth=2)
#plt.scatter(ROT02, MOMENT02, color='orange', linestyle='--', linewidth=2)
plt.title('MOMENT-ROTATION DIAGRAM')
plt.ylabel('Moment [kN.mm]')
plt.xlabel('Rotation [rad]')
#plt.legend(['Steel01: WITHOUT HARDENING AND ULTIMATE STRAIN', 'Hysteretic: WITH HARDENING AND ULTIMATE STRAIN'])
plt.grid()
plt.show()

plt.figure(5, figsize=(12, 8))
#plt.plot(KI, KS, color='black', linewidth=2)
#plt.plot(KI02, KS02, color='grey', linestyle='--', linewidth=2)
#plt.scatter(KI, KS, color='black', linewidth=2)
plt.scatter(KI02, KS02, color='grey', linestyle='--', linewidth=2)
plt.title('ROTATIONAL STIFFNESS-LATERAL STIFFNESS DIAGRAM')
plt.ylabel('Rotational Stiffness [N.mm/Rad]')
plt.xlabel('Lateral Stiffness in X Dir. [N/mm]')
#plt.legend(['Steel01: WITHOUT HARDENING AND ULTIMATE STRAIN', 'Hysteretic: WITH HARDENING AND ULTIMATE STRAIN'])
plt.semilogx()
plt.semilogy()
plt.grid()
plt.show()

plt.figure(6, figsize=(12, 8))
#plt.plot(KI, KA, color='black', linewidth=2)
#plt.plot(KI02, KA02, color='grey', linestyle='--', linewidth=2)
#plt.scatter(KI, KA, color='black', linewidth=2)
plt.scatter(KI02, KA02, color='grey', linestyle='--', linewidth=2)
plt.title('ROTATIONAL STIFFNESS-LATERAL STIFFNESS DIAGRAM')
plt.ylabel('Rotational Stiffness [N.mm/Rad]')
plt.xlabel('Lateral Stiffness in Y Dir. [N/mm]')
#plt.legend(['Steel01: WITHOUT HARDENING AND ULTIMATE STRAIN', 'Hysteretic: WITH HARDENING AND ULTIMATE STRAIN'])
plt.semilogx()
plt.semilogy()
plt.grid()
plt.show()

plt.figure(7, figsize=(12, 8))
#plt.plot(STEP, FORCE_A, color='brown', linewidth=2)
plt.plot(STEP02, FORCE_A02, color='gold', linestyle='--', linewidth=2)
#plt.scatter(STEP, FORCE_A, color='brown', linewidth=2)
#plt.scatter(STEP02, FORCE_A02, color='gold', linestyle='--', linewidth=2)
plt.title('Axial Force During the Analysis')
plt.ylabel('Axial Force [N]')
plt.xlabel('Steps')
#plt.legend(['Steel01: WITHOUT HARDENING AND ULTIMATE STRAIN', 'Hysteretic: WITH HARDENING AND ULTIMATE STRAIN'])
plt.grid()
plt.show()

plt.figure(8, figsize=(12, 8))
#plt.plot(STEP, FORCE_S, color='purple', linewidth=2)
plt.plot(STEP02, FORCE_S02, color='#BF77F6', linestyle='--', linewidth=2)
#plt.scatter(STEP, FORCE_S, color='purple', linewidth=2)
#plt.scatter(STEP02, FORCE_S02, color='#BF77F6', linestyle='--', linewidth=2)
plt.title('Shear Force During the Analysis')
plt.ylabel('Shear Force [N]')
plt.xlabel('Steps')
#plt.legend(['Steel01: WITHOUT HARDENING AND ULTIMATE STRAIN', 'Hysteretic: WITH HARDENING AND ULTIMATE STRAIN'])
plt.grid()
plt.show()

plt.figure(9, figsize=(12, 8))
#plt.plot(STEP, MOMENT, color='green', linewidth=2)
plt.plot(STEP02, MOMENT02, color='lime', linestyle='--', linewidth=2)
#plt.scatter(STEP, MOMENT, color='green', linewidth=2)
#plt.scatter(STEP02, MOMENT02, color='lime', linestyle='--', linewidth=2)
plt.title('Moment During the Analysis')
plt.ylabel('Moment [N.mm]')
plt.xlabel('Steps')
#plt.legend(['Steel01: WITHOUT HARDENING AND ULTIMATE STRAIN', 'Hysteretic: WITH HARDENING AND ULTIMATE STRAIN'])
plt.grid()
plt.show()

plt.figure(10, figsize=(12, 8))
#plt.plot(STEP, DISP_X, color='brown', linewidth=2)
plt.plot(STEP02, DISP_X02, color='gold', linestyle='--', linewidth=2)
#plt.scatter(STEP, DISP_X, color='brown', linewidth=2)
#plt.scatter(STEP02, DISP_X02, color='gold', linestyle='--', linewidth=2)
plt.title('Displacement During the Analysis')
plt.ylabel('Displacement - X [mm]')
plt.xlabel('Steps')
#plt.legend(['Steel01: WITHOUT HARDENING AND ULTIMATE STRAIN', 'Hysteretic: WITH HARDENING AND ULTIMATE STRAIN'])
plt.grid()
plt.show()

plt.figure(11, figsize=(12, 8))
#plt.plot(STEP, DISP_Y, color='blue', linewidth=2)
plt.plot(STEP02, DISP_Y02, color='cyan', linestyle='--', linewidth=2)
#plt.scatter(STEP, DISP_X, color='brown', linewidth=2)
#plt.scatter(STEP02, DISP_X02, color='gold', linestyle='--', linewidth=2)
plt.title('Displacement During the Analysis')
plt.ylabel('Displacement - Y [mm]')
plt.xlabel('Steps')
#plt.legend(['Steel01: WITHOUT HARDENING AND ULTIMATE STRAIN', 'Hysteretic: WITH HARDENING AND ULTIMATE STRAIN'])
plt.grid()
plt.show()

plt.figure(12, figsize=(12, 8))
#plt.plot(STEP, ROT, color='black', linewidth=2)
plt.plot(STEP02, ROT02, color='grey', linestyle='--', linewidth=2)
#plt.scatter(STEP, ROT, color='black', linewidth=2)
#plt.scatter(STEP02, ROT02, color='grey', linestyle='--', linewidth=2)
plt.title('Rotation During the Analysis')
plt.ylabel('Rotation [rad]')
plt.xlabel('Steps')
#plt.legend(['Steel01: WITHOUT HARDENING AND ULTIMATE STRAIN', 'Hysteretic: WITH HARDENING AND ULTIMATE STRAIN'])
plt.grid()
plt.show()
#%%------------------------------------------------------------------------------
# ---------------------------------------
#  Plot BaseShear-Displacement Analysis
# ---------------------------------------
XX = np.abs(DISP_X02); YY = np.abs(FORCE_S02); # ABSOLUTE VALUE
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
BC.PLOT_2D(np.abs(DISP_X02), np.abs(FORCE_S02), X, Y, X, Y, XLABEL, YLABEL, TITLE, LEGEND01, LEGEND02, LEGEND03, COLOR='black', Z=2) 
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
S04.PLOT_2D_FRAME(deformed_scale=1)  # Adjust scale factor as needed
#%%------------------------------------------------------------------------------   
# EXPORT DATA TO EXCEL
import pandas as pd
DATA_TOTAL = { 
    #'DISP_X_WO': DISP_X,# WITHOUT HARDENING AND ULTIMATE STRAIN
    'DISP_X_W': DISP_X02,# WITH HARDENING AND ULTIMATE STRAIN
    #'DISP_Y_WO': DISP_Y,
    'DISP_Y_W': DISP_Y02,
    #'ROTATION_WO': ROT,
    'ROTATION_W': ROT02,
    #'AXIAL_FORCE_WO': FORCE_A,
    'AXIAL_FORCE_W': FORCE_A02,
    #'SHEAR_FORCE_WO': FORCE_S,
    'SHEAR_FORCE_W': FORCE_S02,
    #'MOMENT_WO': MOMENT,
    'MOMENT_W': MOMENT02,
    #'AXIAL_RIGIDITY_WO': np.abs(FORCE_A),
    'AXIAL_RIGIDITY_W': np.abs(FORCE_A02),
    #'ROTATIONAL_ST_WO': KI,
    'ROTATIONAL_ST_W': KI02,
}
# Convert to DataFrame
results_df = pd.DataFrame(DATA_TOTAL)
# Export the DataFrame to an Excel file
results_df.to_excel('CONCRETE_FRAME_PUSHOVER_STRUCTURAL_BEHAVIOR_COEFFICIENT_RESULTS.xlsx', index=False) 
#%%------------------------------------------------------------------------------

    
