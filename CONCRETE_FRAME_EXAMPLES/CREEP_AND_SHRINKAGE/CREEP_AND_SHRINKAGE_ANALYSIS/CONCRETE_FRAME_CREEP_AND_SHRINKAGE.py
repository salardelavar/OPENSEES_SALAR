################################################################################################################
#                                                  IN THE NAME OF ALLAH                                        #
#              ANALYZING CREEP AND SHRINKAGE OF A CONCRETE FRAME. EVALUATING STRAIN HARDENING                  #
#           AND ULTIMATE STRAIN CRITERIA USING OPENSEES AND CALCULATE STRUCTURAL BEHAVIOR COEFFICIENT          #
#--------------------------------------------------------------------------------------------------------------#
#                             THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                       #
#                                        EMAIL: salar.d.ghashghaei@gmail.com                                   #
################################################################################################################
"""
1. The analysis compares nonlinear rotational behavior of concrete beam-column
 elements under creep and shrinkage using OpenSees.
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
def CREEP_ANALYSIS(LENGTH_COL, LENGTH_BM, DMAX, DINCR, STEEL_KIND):
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
        print(step+1, disp_X, disp_Y, S, A)
        
    #ops.wipe()    

    return FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y, ROT, KA, KS, KI, STEP, TIME

#%%------------------------------------------------------------------------------
# Analysis Durations:
starttime = TI.process_time()

# WITHOUT HARDENING AND ULTIMATE STRAIN
FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y, ROT, KA, KS, KI, STEP, TIME = CREEP_ANALYSIS(LENGTH_COL, LENGTH_BM, TMAX, TINCR, STEEL_KIND=1)
# WITH HARDENING AND ULTIMATE STRAIN
FORCE_S02, FORCE_A02, MOMENT02, DISP_X02, DISP_Y02, ROT02, KA02, KS02, KI02, STEP02, TIME02  = CREEP_ANALYSIS(LENGTH_COL, LENGTH_BM, TMAX, TINCR, STEEL_KIND=2)

totaltime = TI.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n') 
#%%------------------------------------------------------------------------------
plt.figure(1, figsize=(12, 8))
plt.plot(MOMENT, FORCE_A, color='black')
plt.plot(MOMENT02, FORCE_A02, color='cyan', linestyle='--', linewidth=2)
#plt.scatter(MOMENT, FORCE_A, color='blue', linewidth=2)
#plt.scatter(MOMENT02, FORCE_A02, color='cyan', linestyle='--', linewidth=2)
plt.title('P-M Interaction')
plt.ylabel('Axial Force [N]')
plt.xlabel('Bending Moment [N.mm]')
plt.legend(['Steel01: WITHOUT HARDENING AND ULTIMATE STRAIN', 'Hysteretic: WITH HARDENING AND ULTIMATE STRAIN'])
plt.grid()
plt.show()

plt.figure(2, figsize=(12, 8))
plt.plot(DISP_X, FORCE_S, color='green', linewidth=2)
plt.plot(DISP_X02, FORCE_S02, color='lime', linestyle='--', linewidth=2)
#plt.scatter(DISP_X, FORCE_S, color='purple', linewidth=2)
#plt.scatter(DISP_X02, FORCE_S02, color='magenta', linestyle='--', linewidth=2)
plt.title('SHEAR FORCE-DISPLACEMENT DIAGRAM')
plt.ylabel('Shear Force [N]')
plt.xlabel('Displacement  in X [mm] NODE 3')
plt.legend(['Steel01: WITHOUT HARDENING AND ULTIMATE STRAIN', 'Hysteretic: WITH HARDENING AND ULTIMATE STRAIN'])
plt.grid()
plt.show()

plt.figure(3, figsize=(12, 8))
plt.plot(DISP_Y, FORCE_A, color='purple', linewidth=2)
plt.plot(DISP_Y02, FORCE_A02, color='magenta', linestyle='--', linewidth=2)
#plt.scatter(DISP_Y, FORCE_A, color='purple', linewidth=2)
#plt.scatter(DISP_Y02, FORCE_A02, color='magenta', linestyle='--', linewidth=2)
plt.title('AXIAL FORCE-DISPLACEMENT DIAGRAM')
plt.ylabel('Axial Force [N]')
plt.xlabel('Displacement in Y [mm] NODE 5')
plt.legend(['Steel01: WITHOUT HARDENING AND ULTIMATE STRAIN', 'Hysteretic: WITH HARDENING AND ULTIMATE STRAIN'])
plt.grid()
plt.show()

plt.figure(4, figsize=(12, 8))
plt.plot(ROT, MOMENT, color='red', linewidth=2)
plt.plot(ROT02, MOMENT02, color='orange', linestyle='--', linewidth=2)
#plt.scatter(ROT, MOMENT, color='red', linewidth=2)
#plt.scatter(ROT02, MOMENT02, color='orange', linestyle='--', linewidth=2)
plt.title('MOMENT-ROTATION DIAGRAM')
plt.ylabel('Moment [kN.mm]')
plt.xlabel('Rotation [rad]')
plt.legend(['Steel01: WITHOUT HARDENING AND ULTIMATE STRAIN', 'Hysteretic: WITH HARDENING AND ULTIMATE STRAIN'])
plt.grid()
plt.show()

plt.figure(5, figsize=(12, 8))
#plt.plot(KI, KS, color='black', linewidth=2)
#plt.plot(KI02, KS02, color='grey', linestyle='--', linewidth=2)
plt.scatter(KI, KS, color='black', linewidth=2)
plt.scatter(KI02, KS02, color='grey', linestyle='--', linewidth=2)
plt.title('ROTATIONAL STIFFNESS-LATERAL STIFFNESS DIAGRAM')
plt.ylabel('Rotational Stiffness [N.mm/Rad]')
plt.xlabel('Lateral Stiffness in X Dir. [N/mm]')
plt.legend(['Steel01: WITHOUT HARDENING AND ULTIMATE STRAIN', 'Hysteretic: WITH HARDENING AND ULTIMATE STRAIN'])
plt.semilogx()
plt.semilogy()
plt.grid()
plt.show()

plt.figure(6, figsize=(12, 8))
#plt.plot(KI, KA, color='black', linewidth=2)
#plt.plot(KI02, KA02, color='grey', linestyle='--', linewidth=2)
plt.scatter(KI, KA, color='black', linewidth=2)
plt.scatter(KI02, KA02, color='grey', linestyle='--', linewidth=2)
plt.title('ROTATIONAL STIFFNESS-LATERAL STIFFNESS DIAGRAM')
plt.ylabel('Rotational Stiffness [N.mm/Rad]')
plt.xlabel('Lateral Stiffness in Y Dir. [N/mm]')
plt.legend(['Steel01: WITHOUT HARDENING AND ULTIMATE STRAIN', 'Hysteretic: WITH HARDENING AND ULTIMATE STRAIN'])
plt.semilogx()
plt.semilogy()
plt.grid()
plt.show()

plt.figure(7, figsize=(12, 8))
plt.plot(STEP, FORCE_A, color='brown', linewidth=2)
plt.plot(STEP02, FORCE_A02, color='gold', linestyle='--', linewidth=2)
#plt.scatter(STEP, FORCE_A, color='brown', linewidth=2)
#plt.scatter(STEP02, FORCE_A02, color='gold', linestyle='--', linewidth=2)
plt.title('Axial Force During the Analysis')
plt.ylabel('Axial Force [N]')
plt.xlabel('Steps')
plt.legend(['Steel01: WITHOUT HARDENING AND ULTIMATE STRAIN', 'Hysteretic: WITH HARDENING AND ULTIMATE STRAIN'])
plt.grid()
plt.show()

plt.figure(8, figsize=(12, 8))
plt.plot(STEP, FORCE_S, color='purple', linewidth=2)
plt.plot(STEP02, FORCE_S02, color='#BF77F6', linestyle='--', linewidth=2)
#plt.scatter(STEP, FORCE_S, color='purple', linewidth=2)
#plt.scatter(STEP02, FORCE_S02, color='#BF77F6', linestyle='--', linewidth=2)
plt.title('Shear Force During the Analysis')
plt.ylabel('Shear Force [N]')
plt.xlabel('Steps')
plt.legend(['Steel01: WITHOUT HARDENING AND ULTIMATE STRAIN', 'Hysteretic: WITH HARDENING AND ULTIMATE STRAIN'])
plt.grid()
plt.show()

plt.figure(9, figsize=(12, 8))
plt.plot(STEP, MOMENT, color='green', linewidth=2)
plt.plot(STEP02, MOMENT02, color='lime', linestyle='--', linewidth=2)
#plt.scatter(STEP, MOMENT, color='green', linewidth=2)
#plt.scatter(STEP02, MOMENT02, color='lime', linestyle='--', linewidth=2)
plt.title('Moment During the Analysis')
plt.ylabel('Moment [kN.mm]')
plt.xlabel('Steps')
plt.legend(['Steel01: WITHOUT HARDENING AND ULTIMATE STRAIN', 'Hysteretic: WITH HARDENING AND ULTIMATE STRAIN'])
plt.grid()
plt.show()

plt.figure(10, figsize=(12, 8))
plt.plot(STEP, DISP_X, color='brown', linewidth=2)
plt.plot(STEP02, DISP_X02, color='gold', linestyle='--', linewidth=2)
#plt.scatter(STEP, DISP_X, color='brown', linewidth=2)
#plt.scatter(STEP02, DISP_X02, color='gold', linestyle='--', linewidth=2)
plt.title('Displacement During the Analysis')
plt.ylabel('Displacement - X [mm] NODE 3')
plt.xlabel('Steps')
plt.legend(['Steel01: WITHOUT HARDENING AND ULTIMATE STRAIN', 'Hysteretic: WITH HARDENING AND ULTIMATE STRAIN'])
plt.grid()
plt.show()

plt.figure(11, figsize=(12, 8))
plt.plot(STEP, DISP_Y, color='blue', linewidth=2)
plt.plot(STEP02, DISP_Y02, color='cyan', linestyle='--', linewidth=2)
#plt.scatter(STEP, DISP_X, color='brown', linewidth=2)
#plt.scatter(STEP02, DISP_X02, color='gold', linestyle='--', linewidth=2)
plt.title('Displacement During the Analysis')
plt.ylabel('Displacement - Y [mm] NODE 5')
plt.xlabel('Steps')
plt.legend(['Steel01: WITHOUT HARDENING AND ULTIMATE STRAIN', 'Hysteretic: WITH HARDENING AND ULTIMATE STRAIN'])
plt.grid()
plt.show()

plt.figure(12, figsize=(12, 8))
plt.plot(STEP, ROT, color='black', linewidth=2)
plt.plot(STEP02, ROT02, color='grey', linestyle='--', linewidth=2)
#plt.scatter(STEP, ROT, color='black', linewidth=2)
#plt.scatter(STEP02, ROT02, color='grey', linestyle='--', linewidth=2)
plt.title('Rotation During the Analysis')
plt.ylabel('Rotation [rad]')
plt.xlabel('Steps')
plt.legend(['Steel01: WITHOUT HARDENING AND ULTIMATE STRAIN', 'Hysteretic: WITH HARDENING AND ULTIMATE STRAIN'])
plt.grid()
plt.show()
#%%------------------------------------------------------------------------------
import BILINEAR_CURVE as BC

# ---------------------------------------
#  Plot BaseShear-Displacement Analysis
# ---------------------------------------
XX = np.abs(DISP_X); YY = np.abs(FORCE_S); # ABSOLUTE VALUE
SLOPE_NODE = 10

DATA = BC.BILNEAR_CURVE(XX, YY, SLOPE_NODE)
X, Y, Elastic_ST, Plastic_ST, Tangent_ST, Ductility_Rito, Over_Strength_Factor = DATA

XLABEL = 'Displacement in X [mm] NODE 3'
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

XLABEL = 'Displacement in Y [mm] NODE 5'
YLABEL = 'Base-Axial [N]'
LEGEND01 = 'Curve'
LEGEND02 = 'Bilinear Fitted'
LEGEND03 = 'Undefined'
TITLE = f'Last Data of BaseAxial-Displacement Analysis - Ductility Ratio: {X[2]/X[1]:.4f} - Over Strength Factor: {Y[2]/Y[1]:.4f}'
COLOR = 'black'
BC.PLOT_2D(np.abs(DISP_Y), np.abs(FORCE_A), X, Y, X, Y, XLABEL, YLABEL, TITLE, LEGEND01, LEGEND02, LEGEND03, COLOR='black', Z=2) 
#print(f'\t\t Ductility Ratio: {YY[2]/YY[1]:.4f}')

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
#%%------------------------------------------------------------------------------
"""
def OUTPUT_SECOND_COLUMN(X, COLUMN):
    import numpy as np
    # Time History
    filename = f"{X}.txt"
    data_collected = np.loadtxt(filename)
    X = data_collected[:, COLUMN]   
    return X 

# STRESS AND STRAIN OF ELEMENT 3
# CONCRETE
strain_B_C = OUTPUT_SECOND_COLUMN('STRESS_STRAIN_BOT_CONCRETE', 2) # Reading bottom concrete strain from Text file
stress_B_C = OUTPUT_SECOND_COLUMN('STRESS_STRAIN_BOT_CONCRETE', 1) # Reading bottom concrete stress from Text file
strain_T_C = OUTPUT_SECOND_COLUMN('STRESS_STRAIN_TOP_CONCRETE', 2) # Reading top concrete strain from Text file
stress_T_C = OUTPUT_SECOND_COLUMN('STRESS_STRAIN_TOP_CONCRETE', 1) # Reading top concrete stress from Text file
# STEEL REBAR
strain_B_R = OUTPUT_SECOND_COLUMN('STRESS_STRAIN_BOT_REBAR', 2) # Reading bottom steel rebar strain from Text file
stress_B_R = OUTPUT_SECOND_COLUMN('STRESS_STRAIN_BOT_REBAR', 1) # Reading bottom steel rebar stress from Text file
strain_T_R = OUTPUT_SECOND_COLUMN('STRESS_STRAIN_TOP_REBAR', 2) # Reading top steel rebar strain from Text file
stress_T_R = OUTPUT_SECOND_COLUMN('STRESS_STRAIN_TOP_REBAR', 1) # Reading top steel rebar stress from Text file

plt.figure(1, figsize=(8, 6))
plt.plot(strain_B_C, stress_B_C, color='blue', label='Bottom Fiber', linewidth=2)
plt.plot(strain_T_C, stress_T_C, color='red', label='Top Fiber', linewidth=2)
plt.xlabel('Strain (mm/mm)')
plt.ylabel('Stress (kN/mm^2)')
plt.title(f'Stress-Strain Relation of Element {3} Concrete Top & Bottom Fibers')
plt.grid()
plt.legend()
plt.show()

plt.figure(2, figsize=(8, 6))
plt.plot(strain_B_R, stress_B_R, color='blue', label='Bottom Fiber', linewidth=2)
plt.plot(strain_T_R, stress_T_R, color='red', label='Top Fiber', linewidth=2)
plt.xlabel('Strain (mm/mm)')
plt.ylabel('Stress (kN/mm^2)')
plt.title(f'Stress-Strain Relation of Element {3} Steel Rebar Top & Bottom Fibers')
plt.grid()
plt.legend()
plt.show()
"""
#%%------------------------------------------------------------------------------  
# %% Plot 2D Frame Shapes
S04.PLOT_2D_FRAME(deformed_scale=10)  # Adjust scale factor as needed
#%%------------------------------------------------------------------------------   
# EXPORT DATA TO EXCEL
import pandas as pd
DATA_TOTAL = { 
    'DISP_X_WO': DISP_X, # WITHOUT HARDENING AND ULTIMATE STRAIN
    'DISP_X_W': DISP_X02,# WITH HARDENING AND ULTIMATE STRAIN
    'DISP_Y_WO': DISP_Y,
    'DISP_Y_W': DISP_Y02,
    'ROTATION_WO': ROT,
    'ROTATION_W': ROT02,
    'AXIAL_FORCE_WO': FORCE_A,
    'AXIAL_FORCE_W': FORCE_A02,
    'SHEAR_FORCE_WO': FORCE_S,
    'SHEAR_FORCE_W': FORCE_S02,
    'MOMENT_WO': MOMENT,
    'MOMENT_W': MOMENT02,
    'AXIAL_RIGIDITY_WO': np.abs(FORCE_A),
    'AXIAL_RIGIDITY_W': np.abs(FORCE_A02),
    'ROTATIONAL_ST_WO': KI,
    'ROTATIONAL_ST_W': KI02,
    'LATERAL_ST_Y_WO': KA,
    'LATERAL_ST_Y_W': KA02,
    'LATERAL_ST_X_WO': KS,
    'LATERAL_ST_X_W': KS02,
}
# Convert to DataFrame
results_df = pd.DataFrame(DATA_TOTAL)
# Export the DataFrame to an Excel file
results_df.to_excel('CONCRETE_FRAME_CREEP_AND_SHRINKAGE_RESULTS.xlsx', index=False) 
#%%------------------------------------------------------------------------------
# Print out the state of all nodes
ops.printModel("node",1, 2, 3, 4, 5)
# Print out the state of element 1 , 2 and 3
ops.printModel("ele", 1, 2 , 3)
# Print the Model
#printModel()
ops.printModel("-JSON", "-file", "CONCRETE_FRAME_CREEP_AND_SHRINKAGE.json")
#%%------------------------------------------------------------------------------  
    
