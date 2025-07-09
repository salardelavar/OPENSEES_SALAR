################################################################################################################
#                                                  IN THE NAME OF ALLAH                                        #
#       SENSITIVITY ANALYSIS OF CONCRETE FRAME BY CHANGING COLUMN REBAR DIAMETER AND COLUMN SECTION DEPTH      #
#              ANALYZING CREEP AND SHRINKAGE OF A CONCRETE FRAME. EVALUATING STRAIN HARDENING                  #
#           AND ULTIMATE STRAIN CRITERIA USING OPENSEES AND CALCULATE STRUCTURAL BEHAVIOR COEFFICIENT          #
#--------------------------------------------------------------------------------------------------------------#
#                             THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                       #
#                                        EMAIL: salar.d.ghashghaei@gmail.com                                   #
################################################################################################################
"""
[1] The analysis compares nonlinear rotational behavior of concrete beam-column
 elements under creep and shrinkage using OpenSees.
[2] Two material models—*Steel01* (bilinear without degradation) and *Hysteretic*
 (tri-linear with pinching and strength/stiffness degradation)—are used.
[3] Both models are subjected to identical loading protocols to investigate pushover
 response under increasing drift demands.
[4] In contrast, the *Hysteretic* model shows strength and stiffness degradation, capturing
 post-peak deterioration and pinching effects.
[5] Element rotation histories reveal increasing divergence as inelastic demand accumulates
 across cycles.
[6] The *Hysteretic* model produces reduced energy dissipation capacity due to pinching and 
cumulative damage.
[7] Peak rotation capacity is reduced in the *Hysteretic* model, indicating realistic modeling
 of damage and failure modes.
[8] The comparison highlights the limitations of bilinear idealizations in capturing cyclic
 degradation in seismic applications.
[9] Advanced modeling with calibrated degradation parameters is essential for accurate
 seismic performance prediction and collapse assessment.
 
BOOK: Creep and Shrinkage, Their Effect on the Behavior of Concrete Structures
'https://link.springer.com/book/10.1007/978-1-4612-5424-9'
WIKOPEDIA: 
'https://en.wikipedia.org/wiki/Creep_and_shrinkage_of_concrete'
PAPER: Experimental investigation on the fundamental behavior of concrete creep
'https://www.sciencedirect.com/science/article/abs/pii/S0950061817313028'
PAPER: Creep and shrinkage of concrete: physical origins and practical measurements
'https://www.sciencedirect.com/science/article/abs/pii/S0029549300003046'
"""
import openseespy.opensees as ops
import matplotlib.pyplot as plt
import numpy as np
import time as TI
import ANALYSIS_FUNCTION as S02
import CONCRETE_SECTION_FUN as S03
import PLOT_2D as S04
import BILINEAR_CURVE as BC
#%%------------------------------------------------------------------------------

# Define materials for nonlinear columns and beam
# Define parameters (units: mm, N)
# ------------------------------------------
# CONCRETE                  tag   f'c        ec0   f'cu        ecu
fc = -25 # [N/mm^2] Nominal concrete compressive strength
Ec = 4700 * np.sqrt(-fc) # [N/mm^2] Concrete Elastic Modulus

# confined concrete - bottom and top section
Kfc = 1.2;			    # ratio of confined to unconfined concrete strength
fc1C = Kfc*fc;		    # CONFINED concrete (mander model), maximum stress
eps1C = 2*fc1C/Ec;	    # strain at maximum stress 
fc2C = 0.2*fc1C;		# ultimate stress
eps2C = 5*eps1C;		# strain at ultimate stress 
EcC = 4700 * np.sqrt(-fc1C) # [N/mm^2] Concrete Elastic Modulus

# unconfined concrete
fc1U = fc;			    # UNCONFINED concrete maximum stress
eps1U = -0.0025;	    # strain at maximum strength of unconfined concrete
fc2U = 0.2*fc1U;		# ultimate stress
eps2U = -0.012;			# strain at ultimate stress
Lambda = 0.1;		    # ratio between unloading slope at $eps2 and initial slope $Ec
EcU = 4700 * np.sqrt(-fc1U) # [N/mm^2] Concrete Elastic Modulus

# tensile-strength properties
ftC = -0.55*fc1C;		# tensile strength +tension
ftU = -0.55*fc1U;		# tensile strength +tension
EtC = ftC/0.002;		# tension softening stiffness
EtU = ftU/0.002;		# tension softening stiffness
 
# STEEL
# Reinforcing steel
fy = 400            # [N/mm²] Steel Rebar Yield Strength   
Es = 2e5            # [N/mm²] Modulus of Elasticity
ey = fy/Es          # [mm/mm] Steel Rebar Yield Strain
fu = 1.1818*fy      # [N/mm²] Steel Rebar Ultimate Strength
esu = ey*75.2       # [mm/mm] Steel Rebar Ultimate Strain
Esh = (fu - fy)/(esu - ey)
Bs = Esh / Es
#%%------------------------------------------------------------------------------
# Column Section
Bc = 500                  # [mm] Depth of the Section 
Hc = 500                  # [mm] Height of the Section  
coverC = 50               # [mm] Concrete Section Cover
DIAc = 25                 # [mm] # Rebar Size Diameter
AsC = np.pi*(DIAc**2)/4   # [mm²] Area of Rebar
#%%------------------------------------------------------------------------------
# Beam Section
Bb = 500                 # [mm] Depth of the Section 
Hb = 300                 # [mm] Height of the Section  
coverB = 50              # [mm] Concrete Section Cover
DIAb = 18                # [mm] # Rebar Size Diameter
AsB = np.pi*(DIAb**2)/4  # [mm²] Area of Rebar
#%%------------------------------------------------------------------------------
N_YEARS = 5              # Number of Years
TMAX = 365 * N_YEARS     # [day] Maximum Days
TINCR = 25               # [day] Increment Day

# TDConcrete
# LINK INFO: https://openseespydoc.readthedocs.io/en/latest/src/TDConcrete.html
BETA = 0.4       # Tension softening parameter (tension softening exponent) (Recommended value)
tD = 14          # Analysis time at initiation of drying (in days)
epsshu = 400e-6  # Ultimate shrinkage strain as per ACI 209R-92 (shrinkage is negative)
psish = 64.174   # Fitting parameter of the shrinkage time evolution function as per ACI 209R-92
Tcr = 28         # Creep model age (in days)

phiu = 3.0       # Ultimate creep coefficient as per ACI 209R-92
psicr1 = 1.0     # Fitting parameter of the creep time evolution function as per ACI 209R-92 (Recommended value) 
psicr2 = 75.4218 # Fitting parameter of the creep time evolution function as per ACI 209R-92 (Based on section dimensions)
tcast = 2        # Analysis time corresponding to concrete casting (in days; minimum value 2.0)

# TDConcreteEXP
# INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/TDConcreteEXP.html
epscru = 0.0002  # Ultimate creep strain (e.g., taken from experimental measurements)
sigCr = -5       # Concrete compressive stress (input as negative) associated with $epscru (e.g., experimentally applied)
psicr1 = 1.0     # Fitting parameter of the creep time evolution function as per ACI 209R-92 (Recommended value) 
psicr2 = 75.4218 # Fitting parameter of the creep time evolution function as per ACI 209R-92 (Based on section dimensions)
tcast = 2        # Analysis time corresponding to concrete casting (in days; minimum value 2.0)

# TDConcreteMC10
# INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/TDConcreteMC10.html

# TDConcreteMC10NL
# INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/TDConcreteMC10NL.html
#%%------------------------------------------------------------------------------
LENGTH_COL = 3000        # [mm] Column Length 
LENGTH_BM = 7000         # [mm] Beam Length
 #%%------------------------------------------------------------------------------
# POINT LOAD IN EACH NODE
#               [N]   [N]     [N.mm]
PX3, PY3, MZ3 = 1.0, -10000, -25e4
PX4, PY4, MZ4 = 1.0, -10000, +25e4
dl = -0.1                # [N/mm] Distributed Load
#%%------------------------------------------------------------------------------
# Define Analysis Properties
MAX_ITERATIONS = 5000    # Convergence iteration for test
MAX_TOLERANCE = 1.0e-10  # Convergence tolerance for test
#STEEL_KIND: 1 -> WITHOUT HARDENING AND ULTIMATE STRAIN
#STEEL_KIND: 2 -> WITH HARDENING AND ULTIMATE STRAIN
#%%------------------------------------------------------------------------------
def CREEP_ANALYSIS(DIAc, Hc, LENGTH_COL, LENGTH_BM, DMAX, DINCR, STEEL_KIND):
    
    #DIAc = 25                # [mm] # Rebar Size Diameter
    AsC = np.pi*(DIAc**2)/4   # [mm²] Area of Rebar
    RO_COL = (8*AsC)/(Hc*Bc)  # Column Rebar Ratio
    
    DIAb = DIAc-4             # [mm] # Rebar Size Diameter
    AsB = np.pi*(DIAb**2)/4   # [mm²] Area of Rebar
    RO_BE = (6*AsB)/(Hb*Bb)   # Beam Rebar Ratio
    # Reset model
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)

    # Nodes
    ops.node(1, 0.0, 0.0)
    ops.node(2, LENGTH_BM, 0.0)
    ops.node(3, 0.0, LENGTH_COL)
    ops.node(4, LENGTH_BM, LENGTH_COL)
    ops.node(5, (1/2)*LENGTH_BM, LENGTH_COL)
    
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

    # Time-dependent concrete material - Core concrete (confined)
    ops.uniaxialMaterial('TDConcrete', coreTag, fc1C, ftC, EcC, BETA, tD, -epsshu, psish, Tcr, phiu, psicr1, psicr2, tcast) 
    # Time-dependent concrete material - Cover concrete (unconfined)
    ops.uniaxialMaterial('TDConcrete', coverTag, fc1U, ftU, EcU, BETA, tD, -epsshu, psish, Tcr, phiu, psicr1, psicr2, tcast) 
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/TDConcrete.html
    
    #ops.uniaxialMaterial('TDConcreteEXP', coverTag, fc1U, ftU, EcU, BETA, tD, -epsshu, psish, Tcr, epscru, sigCr, psicr1, psicr2, tcast)
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/TDConcreteEXP.html
    
    # COLUMN SECTION
    S03.CONFINED_CONCRETE_SECTION(secTagC, Hc, Bc, coverC, AsC, coreTag, coverTag, steelTag, COL=True)
    # BEAM SECTION
    S03.CONFINED_CONCRETE_SECTION(secTagB, Hb, Bb, coverB, AsB, coreTag, coverTag, steelTag, COL=False)
    """
    DD = 0.5 * Hb
    # Output Data
    ops.recorder('Element', '-file', 'STRESS_STRAIN_BOT_CONCRETE.txt', '-time', '-ele', 3, 'section', secTagB, 'fiber', -DD, 0, 'stressStrain')
    ops.recorder('Element', '-file', 'STRESS_STRAIN_TOP_CONCRETE.txt', '-time', '-ele', 3, 'section', secTagB, 'fiber', +DD, 0, 'stressStrain')
    ops.recorder('Element', '-file', 'STRESS_STRAIN_BOT_REBAR.txt', '-time', '-ele', 3, 'section', secTagB, 'fiber', -DD+coverB, -0.5*Bb+coverB, 'stressStrain')
    ops.recorder('Element', '-file', 'STRESS_STRAIN_TOP_REBAR.txt', '-time', '-ele', 3, 'section', secTagB, 'fiber', +DD-coverB, -0.5*Bb+coverB, 'stressStrain')
    """
    # Define geometric transformation (corotational for large displacements)
    transfTag = 1
    #ops.geomTransf('Linear', transfTag)
    #ops.geomTransf('PDelta', transfTag)
    ops.geomTransf('Corotational', transfTag)
    numIntgrPts = 5
    ops.element('nonlinearBeamColumn', 1, 1, 3, numIntgrPts, secTagC, transfTag) # COLUMN 01
    ops.element('nonlinearBeamColumn', 2, 2, 4, numIntgrPts, secTagC, transfTag) # COLUMN 02
    ops.element('nonlinearBeamColumn', 3, 3, 5, numIntgrPts, secTagB, transfTag) # BEAM 01
    ops.element('nonlinearBeamColumn', 4, 5, 4, numIntgrPts, secTagB, transfTag) # BEAM 02

    # Data storage
    FORCE_S, FORCE_A, MOMENT = [], [], []
    DISP_X, DISP_Y, ROT = [], [], []
    KA, KS, KI, STEP, TIME = [], [], [], [], []

    # Define time series and load pattern
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(3, PX3, PY3, MZ3)
    ops.load(4, PX4, PY4, MZ4)
    
    # Apply distributed load to all beams in the structure
    # eleLoad -ele $eleTag1 <$eleTag2 ....> -type -beamUniform $Wz <$Wx>
    #ops.eleLoad('-ele', 3, '-type', '-beamUniform', dl) 

    # Total steps
    steps = int(np.abs(TMAX)/np.abs(TINCR))
    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.system('BandGeneral')
    ops.algorithm('Newton')
    ops.setTime(Tcr)
    ops.integrator('LoadControl', 0)
    ops.analysis('Static')
    ops.analyze(1)
    
    ops.setCreep(1) # Turn creep on    
    ops.integrator('LoadControl', TINCR)
    t = 0
    for step in range(steps):
        OK = ops.analyze(1)
        S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS
        # Record results
        ops.reactions()
        S = ops.nodeReaction(1, 1) + ops.nodeReaction(2, 1) # SHEAR BASE REACTION
        A = ops.nodeReaction(1, 2) + ops.nodeReaction(2, 2) # AXIAL BASE REACTION
        M = ops.nodeReaction(1, 3) + ops.nodeReaction(2, 3) # MOMENT BASE REACTION
        #print(rot, M)
        disp_X = ops.nodeDisp(3, 1) # LATERAL DISPLACEMENT IN X FOR NODE 3
        disp_Y = ops.nodeDisp(5, 2) # LATERAL DISPLACEMENT IN Y FOR NODE 5
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
        t += TINCR
        TIME.append(t)
        #print(step+1, disp_X, disp_Y, S, A)
        
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

    return RO_COL, RO_BE, FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y, ROT, KA, KS, KI, STEP, TIME, Elastic_ST, Plastic_ST, Tangent_ST, Ductility_Rito, Over_Strength_Factor, R

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
DIAc = [18, 20, 22, 28, 30, 32]
# DEPTH OF COLUMN SECTION
#Hc = [300, 325, 350, 370, 390, 400, 415, 435, 455, 475, 500, 510, 525, 550]
Hc = [500, 510, 525, 550, 580, 600]

II = 0
for dia in DIAc:
    for hc in Hc:
        II = II + 1
        print(f'\n STEP: {II} - REBAR DIAMETER: {dia} - HC: {hc} \n')
        DIAc_MAX.append(dia)
        Hc_MAX.append(hc)
        DATA  = CREEP_ANALYSIS(dia, hc, LENGTH_COL, LENGTH_BM, TMAX, TINCR, STEEL_KIND=2)
        
        (RO_COL, RO_BE, FORCE_S, FORCE_A, MOMENT, DISP_X,
         DISP_Y, ROT, KA, KS, KI, STEP, TIME, Elastic_ST,
         Plastic_ST, Tangent_ST, Ductility_Rito, Over_Strength_Factor, R) = DATA
        
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

#%%------------------------------------------------------------------------------
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
Z = FORCE_S_MAX
XLABEL = 'Rebar Diameter [mm]'   
YLABEL = 'Column Section Depth [mm]' 
ZLABEL = 'Structure Shear Reaction in X Dir. [N]'
PLOT_3D(1, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = DIAc_MAX
Y = Hc_MAX
Z =  FORCE_A_MAX
XLABEL = 'Rebar Diameter [mm]'   
YLABEL = 'Column Section Depth [mm]' 
ZLABEL = 'Structure Axial Reaction in Y Dir. [N]'
PLOT_3D(2, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = DIAc_MAX
Y = Hc_MAX
Z =  MOMENT_MAX
XLABEL = 'Rebar Diameter [mm]'   
YLABEL = 'Column Section Depth [mm]' 
ZLABEL = 'Structure Moment Reaction in Y Dir. [N.mm]'
PLOT_3D(3, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = DIAc_MAX
Y = Hc_MAX
Z =  DISP_X_MAX
XLABEL = 'Rebar Diameter [mm]'   
YLABEL = 'Column Section Depth [mm]' 
ZLABEL = 'Structure Displacement in X Dir. [mm]'
PLOT_3D(4, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = DIAc_MAX
Y = Hc_MAX
Z =  DISP_Y_MAX
XLABEL = 'Rebar Diameter [mm]'   
YLABEL = 'Column Section Depth [mm]' 
ZLABEL = 'Structure Displacement in Y Dir. [mm]'
PLOT_3D(5, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = DIAc_MAX
Y = Hc_MAX
Z =  Elastic_ST_MAX
XLABEL = 'Rebar Diameter [mm]'   
YLABEL = 'Column Section Depth [mm]' 
ZLABEL = 'Structure Elastic Stiffness [N/mm]'
PLOT_3D(6, X, Y, Z, XLABEL, YLABEL, ZLABEL)    

X = DIAc_MAX
Y = Hc_MAX
Z =  Plastic_ST_MAX
XLABEL = 'Rebar Diameter [mm]'   
YLABEL = 'Column Section Depth [mm]' 
ZLABEL = 'Structure Plastic Stiffness [N/mm]'
PLOT_3D(7, X, Y, Z, XLABEL, YLABEL, ZLABEL)    

X = DIAc_MAX
Y = Hc_MAX
Z =  Ductility_Rito_MAX
XLABEL = 'Rebar Diameter [mm]'   
YLABEL = 'Column Section Depth [mm]' 
ZLABEL = 'Ductility Rito  [mm/mm]'
PLOT_3D(8, X, Y, Z, XLABEL, YLABEL, ZLABEL)   

X = DIAc_MAX
Y = Hc_MAX
Z =  Over_Strength_Factor_MAX
XLABEL = 'Rebar Diameter [mm]'   
YLABEL = 'Column Section Depth [mm]' 
ZLABEL = 'Over Strength Factor [N/N]'
PLOT_3D(9, X, Y, Z, XLABEL, YLABEL, ZLABEL)   

X = DIAc_MAX
Y = Hc_MAX
Z = R_MAX
XLABEL = 'Rebar Diameter [mm]'   
YLABEL = 'Column Section Depth [mm]' 
ZLABEL = 'Structural Behavior Coefficient'
PLOT_3D(10, X, Y, Z, XLABEL, YLABEL, ZLABEL) 
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
    'DIAc': DIAc_MAX,
    'Hc': Hc_MAX,
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
results_df.to_excel('CONCRETE_FRAME_CREEP_AND_SHRINKAGE_SENSITIVITY_REBAR_Cdepth_RESULTS.xlsx', index=False) 
#%%------------------------------------------------------------------------------
def PLOT_HEATMAP(df):
    import matplotlib.pyplot as plt
    import numpy as np

    # Calculate the correlation matrix
    corr_matrix = df.corr()

    # Create the figure and axis
    fig, ax = plt.subplots(figsize=(12, 12))

    # Create the heatmap
    cax = ax.matshow(corr_matrix, cmap='viridis')

    # Add the colorbar
    fig.colorbar(cax)

    # Set axis labels
    ax.set_xticks(np.arange(len(corr_matrix.columns)))
    ax.set_yticks(np.arange(len(corr_matrix.index)))
    ax.set_xticklabels(corr_matrix.columns, rotation=90)
    ax.set_yticklabels(corr_matrix.index)

    # Annotate the heatmap with the correlation values
    for i in range(len(corr_matrix.columns)):
        for j in range(len(corr_matrix.index)):
            ax.text(j, i, f'{corr_matrix.iloc[i, j]:.2f}', 
                    ha='center', va='center', color='white')

    # Set title and layout
    plt.title('Correlation Heatmap', pad=20)
    plt.xlabel('Variable')
    plt.ylabel('Variable')
    
    # Display the plot
    plt.tight_layout()
    plt.show()

# PLOT HEATMAP FOR CORRELATION 
PLOT_HEATMAP(DATA_TOTAL)
#%%------------------------------------------------------------------------------
# Print out the state of all nodes
ops.printModel("node",1, 2, 3, 4, 5)
# Print out the state of element 1 , 2 and 3
ops.printModel("ele", 1, 2 , 3)
# Print the Model
#printModel()
ops.printModel("-JSON", "-file", "CONCRETE_FRAME_CREEP_AND_SHRINKAGE.json")
#%%------------------------------------------------------------------------------  
    
