################################################################################################################
#                                                  IN THE NAME OF ALLAH                                        #
#       SENSITIVITY ANALYSIS OF CONCRETE FRAME BY CHANGING COLUMN REBAR DIAMETER AND COLUMN SECTION DEPTH      # 
#                        USING OPENSEES FOR STRUCTURAL BEHAVIOR COEFFICIENT CALCULATION                        #
#--------------------------------------------------------------------------------------------------------------#
#                             THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                       #
#                                        EMAIL: salar.d.ghashghaei@gmail.com                                   #
################################################################################################################
"""
1. The analysis compares nonlinear rotational behavior of concrete beam-column
 elements under pushover lateral displacements using OpenSees.
2. Two material models—*Steel01* (bilinear without degradation) and *Hysteretic*
 (tri-linear with pinching and strength/stiffness degradation)—are used.
3. Both models are subjected to identical loading protocols to investigate pushover
 response under increasing drift demands.
4. The *Steel01* model exhibits stable hysteresis loops with no degradation, reflecting
 idealized elastic–plastic behavior.
5. In contrast, the *Hysteretic* model shows strength and stiffness degradation, capturing
 post-peak deterioration and pinching effects.
6. Element rotation histories reveal increasing divergence as inelastic demand accumulates
 across cycles.
7. The *Hysteretic* model produces reduced energy dissipation capacity due to pinching and 
cumulative damage.
8. Peak rotation capacity is reduced in the *Hysteretic* model, indicating realistic modeling
 of damage and failure modes.
9. The comparison highlights the limitations of bilinear idealizations in capturing cyclic
 degradation in seismic applications.
10. Advanced modeling with calibrated degradation parameters is essential for accurate
 seismic performance prediction and collapse assessment.
"""
import openseespy.opensees as ops
import matplotlib.pyplot as plt
import numpy as np
import time as TI
import ANALYSIS_FUNCTION as S02
import CONCRETE_SECTION_FUN as S03
import BILINEAR_CURVE as BC
import PLOT_2D as S04


# Define materials for nonlinear columns and beam
# Define parameters (units: mm, N)
# ------------------------------------------
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


# Beam Section
Bb = 500                 # [mm] Depth of the Section 
Hb = 300                 # [mm] Height of the Section  
coverB = 50              # [mm] Concrete Section Cover


DMAX = 200               # [mm] Maximum Displacement
DINCR = 0.001            # [mm] Increment Displacement


LENGTH_COL = 3000        # [mm] Column Length 
LENGTH_BM = 7000         # [mm] Beam Length 
# Define Analysis Properties
MAX_ITERATIONS = 5000    # Convergence iteration for test
MAX_TOLERANCE = 1.0e-10  # Convergence tolerance for test
#STEEL_KIND: 1 -> WITHOUT HARDENING AND ULTIMATE STRAIN
#STEEL_KIND: 2 -> WITH HARDENING AND ULTIMATE STRAIN
#CONCRETE_KIND: 1 -> 'Concrete01'
#CONCRETE_KIND: 2 -> 'Concrete02'
#%%------------------------------------------------------------------------------
def PUSHOVER_ANALYSIS(DIAc, Hc, LENGTH_COL, LENGTH_BM, DMAX, DINCR, STEEL_KIND, CONCRETE_KIND):
    
    #DIAc = 25                # [mm] # Rebar Size Diameter
    AsC = np.pi*(DIAc**2)/4   # [mm²] Area of Rebar
    RO_COL = (8*AsC)/(Hc*Bc)  # Column Rebar Ratio
    
    DIAb = DIAc-7             # [mm] # Rebar Size Diameter
    AsB = np.pi*(DIAb**2)/4   # [mm²] Area of Rebar
    RO_BE = (6*AsB)/(Hb*Bb)   # Beam Rebar Ratio
    
    if CONCRETE_KIND == 1:        
        # Cover concrete (unconfined)
        fcU = -18           # [N/mm²] Concrete Compressive Strength
        ec0U = -0.0025      # [mm/mm] Concrete Compressive Strain
        fcUU = -2           # [N/mm²] Concrete Compressive Ultimate Strength
        ecuU = -0.008       # [mm/mm] Concrete Compressive Ultimate Strain
        # Core concrete (confined)
        Kfc = 1.2;			# ratio of confined to unconfined concrete strength
        fcC = Kfc*fcU       # [N/mm²] Concrete Compressive Strength
        ec0C = -0.0045      # [mm/mm] Concrete Compressive Strain
        fcUC = -21          # [N/mm²] Concrete Compressive Ultimate Strength
        ecuC = -0.015       # [mm/mm] Concrete Compressive Ultimate Strain
        
    if CONCRETE_KIND == 2:
        # Cover concrete (unconfined)
        fc1U = fcU;			    # [N/mm²] UNCONFINED concrete (todeschini parabolic model), maximum stress
        eps1U = -0.0025;		# [mm/mm] strain at maximum strength of unconfined concrete
        fc2U = 0.2*fc1U;		# [N/mm²] ultimate stress
        eps2U = -0.012;			# [mm/mm] strain at ultimate stress
        Lambda = 0.1;			# ratio between unloading slope at $eps2 and initial slope $Ec
        # Core concrete (confined)
        Kfc = 1.2;			# ratio of confined to unconfined concrete strength
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
        # CONCRETE                            tag   f'c        ec0   f'cu        ecu
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

    # Total steps
    steps = int(np.abs(DMAX)/np.abs(DINCR))

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
    
    return RO_COL, RO_BE, FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y, ROT, KA, KS, KI, Elastic_ST, Plastic_ST, Tangent_ST, Ductility_Rito, Over_Strength_Factor, R, STEP

#%%------------------------------------------------------------------------------
# SENSTIVITY ANALYSIS BY CHANGING REBAR DIAMETER AND COLUMN SECTION DEPTH
# Analysis Durations:
current_time = TI.strftime("%H:%M:%S", TI.localtime())
print(f"Current time (HH:MM:SS): {current_time}\n\n")

RO_COL_MAX, RO_BE_MAX =  [], []
FORCE_S_MAX, FORCE_A_MAX, MOMENT_MAX = [], [], []
DISP_X_MAX, DISP_Y_MAX, ROT_MAX = [], [], []
KA_MAX, KS_MAX, KI_MAX, R_MAX = [], [], [], []
DIAc_MAX, Hc_MAX = [], []
Elastic_ST_MAX, Plastic_ST_MAX, Tangent_ST_MAX, Ductility_Rito_MAX, Over_Strength_Factor_MAX = [], [], [], [], []

# REBAR DIAMETER FOR COLUMN SECTION
#DIAc = [8, 10, 12, 14, 16, 18, 20, 22, 25, 28, 30, 32]
DIAc = [ 8, 10, 12, 14, 16, 18]
# DEPTH OF COLUMN SECTION
#Hc = [300, 325, 350, 370, 390, 400, 415, 435, 455, 475, 500, 510, 525, 550]
Hc = [300, 315, 325, 335, 350]

II = 0
for dia in DIAc:
    for hc in Hc:
        II = II + 1
        print(f'\n STEP: {II} - REBAR DIAMETER: {dia} - HC: {hc} \n')
        DIAc_MAX.append(dia)
        Hc_MAX.append(hc)
        DATA = PUSHOVER_ANALYSIS(dia, hc, LENGTH_COL, LENGTH_BM, DMAX, DINCR, STEEL_KIND=2, CONCRETE_KIND=1)
        RO_COL, RO_BE, FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y, ROT, KA, KS, KI, Elastic_ST, Plastic_ST, Tangent_ST, Ductility_Rito, Over_Strength_Factor, R, STEP = DATA
        RO_COL_MAX.append(RO_COL)
        RO_BE_MAX.append(RO_BE)
        FORCE_S_MAX.append(np.max(np.abs(FORCE_S)))
        FORCE_A_MAX.append(np.max(np.abs(FORCE_A)))
        MOMENT_MAX.append(np.max(np.abs(MOMENT)))
        DISP_X_MAX.append(np.max(np.abs(DISP_X)))
        DISP_Y_MAX.append(np.max(np.abs(DISP_Y)))
        ROT_MAX.append(np.max(np.abs(ROT)))
        KA_MAX.append(np.max(np.abs(KA)))
        KS_MAX.append(np.max(np.abs(KS)))
        KI_MAX.append(np.max(np.abs(KI)))
        Elastic_ST_MAX.append(Elastic_ST)
        Plastic_ST_MAX.append(Plastic_ST)
        Tangent_ST_MAX.append(Tangent_ST)
        Ductility_Rito_MAX.append(Ductility_Rito)
        Over_Strength_Factor_MAX.append(Over_Strength_Factor)
        R_MAX.append(R)
        
current_time = TI.strftime("%H:%M:%S", TI.localtime())
print(f"Current time (HH:MM:SS): {current_time}\n\n")
#%%------------------------------------------------------------------------------
def PLOT_3D(TAG, X, Y, Z, XLABEL, YLABEL, ZLABEL):
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
X = DIAc_MAX
Y = Hc_MAX
Z =  DISP_X_MAX
XLABEL = 'Rebar Diameter [mm]'   
YLABEL = 'Column Section Depth [mm]' 
ZLABEL = 'Structure Displacement in X Dir. [mm]'
PLOT_3D(1, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = DIAc_MAX
Y = Hc_MAX
Z =  DISP_Y_MAX
XLABEL = 'Rebar Diameter [mm]'   
YLABEL = 'Column Section Depth [mm]' 
ZLABEL = 'Structure Displacement in Y Dir. [mm]'
PLOT_3D(2, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = DIAc_MAX
Y = Hc_MAX
Z =  Plastic_ST_MAX
XLABEL = 'Rebar Diameter [mm]'   
YLABEL = 'Column Section Depth [mm]' 
ZLABEL = 'Structure Plastic Stiffness [N/mm]'
PLOT_3D(3, X, Y, Z, XLABEL, YLABEL, ZLABEL)    

X = DIAc_MAX
Y = Hc_MAX
Z =  Ductility_Rito_MAX
XLABEL = 'Rebar Diameter [mm]'   
YLABEL = 'Column Section Depth [mm]' 
ZLABEL = 'Ductility Rito  [mm/mm]'
PLOT_3D(4, X, Y, Z, XLABEL, YLABEL, ZLABEL)   

X = DIAc_MAX
Y = Hc_MAX
Z =  Over_Strength_Factor_MAX
XLABEL = 'Rebar Diameter [mm]'   
YLABEL = 'Column Section Depth [mm]' 
ZLABEL = 'Over Strength Factor [N/N]'
PLOT_3D(5, X, Y, Z, XLABEL, YLABEL, ZLABEL)   

X = DIAc_MAX
Y = Hc_MAX
Z = R_MAX
XLABEL = 'Rebar Diameter [mm]'   
YLABEL = 'Column Section Depth [mm]' 
ZLABEL = 'Structural Behavior Coefficient'
PLOT_3D(6, X, Y, Z, XLABEL, YLABEL, ZLABEL) 
#%%------------------------------------------------------------------------------
def PLOT_LINE(title, xlabel, ylabel, x, y, color='blue', fig_num=1, logx=False, logy=False):
    import matplotlib.pyplot as plt
    plt.figure(fig_num, figsize=(12, 8))
    plt.plot(x, y, color=color, linewidth=2)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if logx:
        plt.semilogx()
    if logy:
        plt.semilogy()
    plt.grid()
    plt.show()

def PLOT_SCATTER(title, xlabel, ylabel, x, y, color='black', fig_num=1, logx=False, logy=False):
    import matplotlib.pyplot as plt
    plt.figure(fig_num, figsize=(12, 8))
    plt.scatter(x, y, color=color, linewidth=2)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if logx:
        plt.semilogx()
    if logy:
        plt.semilogy()
    plt.grid()
    plt.show()

# Call plots using the functions
PLOT_LINE('P-M Interaction', 'Bending Moment [N.mm]', 'Axial Force [N]', MOMENT, FORCE_A, color='black', fig_num=1)
PLOT_LINE('SHEAR FORCE-DISPLACEMENT DIAGRAM', 'Displacement in X [mm]', 'Shear Force [N]', DISP_X, FORCE_S, color='green', fig_num=2)
PLOT_LINE('AXIAL FORCE-DISPLACEMENT DIAGRAM', 'Displacement in Y [mm]', 'Axial Force [N]', DISP_Y, FORCE_A, color='purple', fig_num=3)
PLOT_LINE('MOMENT-ROTATION DIAGRAM', 'Rotation [rad]', 'Moment [kN.mm]', ROT, MOMENT, color='red', fig_num=4)

PLOT_SCATTER('ROTATIONAL STIFFNESS-LATERAL STIFFNESS DIAGRAM (X Dir)', 'Lateral Stiffness in X Dir. [N/mm]', 'Rotational Stiffness [N.mm/Rad]', KI, KS, color='black', fig_num=5, logx=True, logy=True)
PLOT_SCATTER('ROTATIONAL STIFFNESS-LATERAL STIFFNESS DIAGRAM (Y Dir)', 'Lateral Stiffness in Y Dir. [N/mm]', 'Rotational Stiffness [N.mm/Rad]', KI, KA, color='black', fig_num=6, logx=True, logy=True)

PLOT_LINE('Axial Force During the Analysis', 'Steps', 'Axial Force [N]', STEP, FORCE_A, color='brown', fig_num=7)
PLOT_LINE('Shear Force During the Analysis', 'Steps', 'Shear Force [N]', STEP, FORCE_S, color='purple', fig_num=8)
PLOT_LINE('Moment During the Analysis', 'Steps', 'Moment [kN.mm]', STEP, MOMENT, color='green', fig_num=9)
PLOT_LINE('Displacement During the Analysis (X)', 'Steps', 'Displacement - X [mm]', STEP, DISP_X, color='brown', fig_num=10)
PLOT_LINE('Displacement During the Analysis (Y)', 'Steps', 'Displacement - Y [mm]', STEP, DISP_Y, color='blue', fig_num=11)
PLOT_LINE('Rotation During the Analysis', 'Steps', 'Rotation [rad]', STEP, ROT, color='black', fig_num=12)

#%%------------------------------------------------------------------------------
#import BILINEAR_CURVE as BC

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
#R_mu = 1
#R_mu = (2 * mu - 1) ** 0.5
R_mu = mu
# Calculate Structural Behavior Coefficient (R)
R = Omega_0 * R_mu
print(f'Over Strength Coefficient (Ω0):      {Omega_0:.4f}')
print(f'Displacement Ductility Ratio (μ):    {mu:.4f}')
print(f'Ductility Coefficient (Rμ):          {R_mu:.4f}')
print(f'Structural Behavior Coefficient (R): {R:.4f}')

# ---------------------------------------
#  Plot BaseAxial-Displacement Analysis
# ---------------------------------------
XX = np.abs(DISP_Y); YY = np.abs(FORCE_A); # ABSOLUTE VALUE
SLOPE_NODE = 10

DATA = BC.BILNEAR_CURVE(XX, YY, SLOPE_NODE)
X, Y, Elastic_ST, Plastic_ST, Tangent_ST, Ductility_Rito, Over_Strength_Factor = DATA

XLABEL = 'Displacement in Y [mm]'
YLABEL = 'Base-Axial [N]'
LEGEND01 = 'Curve'
LEGEND02 = 'Bilinear Fitted'
LEGEND03 = 'Undefined'
TITLE = f'Last Data of BaseAxial-Displacement Analysis - Ductility Ratio: {X[2]/X[1]:.4f} - Over Strength Factor: {Y[2]/Y[1]:.4f}'
COLOR = 'black'
BC.PLOT_2D(np.abs(DISP_Y), np.abs(FORCE_A), X, Y, X, Y, XLABEL, YLABEL, TITLE, LEGEND01, LEGEND02, LEGEND03, COLOR='black', Z=2) 
print(f'\t\t Ductility Ratio: {YY[2]/YY[1]:.4f}')

# Calculate Over Strength Coefficient (Ω0)
Omega_0 = Y[2] / Y[1]
# Calculate Displacement Ductility Ratio (μ)
mu = X[2] / X[1]
# Calculate Ductility Coefficient (Rμ)
R_mu = (2 * mu - 1) ** 0.5 / mu ** 0.5
# Calculate Structural Behavior Coefficient (R)
R = Omega_0 * R_mu
print(f'Over Strength Coefficient (Ω0):      {Omega_0:.4f}')
print(f'Displacement Ductility Ratio (μ):    {mu:.4f}')
print(f'Ductility Coefficient (Rμ):          {R_mu:.4f}')
print(f'Structural Behavior Coefficient (R): {R:.4f}')

#%%------------------------------------------------------------------------------  
# %% Plot 2D Frame Shapes
S04.PLOT_2D_FRAME(deformed_scale=10)  # Adjust scale factor as needed
#%%------------------------------------------------------------------------------   
# EXPORT DATA TO EXCEL
import pandas as pd
DATA_TOTAL = { 
    'RO_COL': RO_COL_MAX,
    'RO_BE': RO_BE_MAX,
    'DISP_X': DISP_X_MAX, 
    'DISP_Y': DISP_Y_MAX,
    'ROTATION': ROT_MAX,
    'AXIAL_FORCE': FORCE_A_MAX,
    'SHEAR_FORCE': FORCE_S_MAX,
    'MOMENT_WO': MOMENT_MAX,
    'AXIAL_RIGIDITY': np.abs(FORCE_A_MAX),
    'ROTATIONAL_ST': KI_MAX,
    'LATERAL_ST_Y': KA_MAX,
    'LATERAL_ST_X': KS_MAX,
    'R_MAX': R_MAX,
    'Elastic_ST_MAX': Elastic_ST_MAX,
    'Plastic_ST_MAX': Plastic_ST_MAX,
    'Tangent_ST_MAX': Tangent_ST_MAX,
    'Ductility_Rito_MAX': Ductility_Rito_MAX,
    'Over_Strength_Factor_MAX': Over_Strength_Factor_MAX,
    'Structural_Behavior_Coefficient_MAX': R_MAX,
}
# Convert to DataFrame
results_df = pd.DataFrame(DATA_TOTAL)
# Export the DataFrame to an Excel file
results_df.to_excel('CONCRETE_FRAME_SENSITIVITY_REBAR_Cdepth_RESULTS.xlsx', index=False) 
#%%------------------------------------------------------------------------------

    
