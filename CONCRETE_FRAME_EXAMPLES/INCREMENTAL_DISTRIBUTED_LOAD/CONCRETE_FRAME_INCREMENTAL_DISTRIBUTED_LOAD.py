################################################################################################################
#                          >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                    #
#       PUSHOVER ANALYSIS WITH INCREMENTAL RECTANGULAR OR TRIANGULAR DISTRIBUTED LOAD OF CONCRETE FRAME.       #
#                  EVALUATING STRAIN HARDENING AND ULTIMATE STRAIN CRITERIA USING OPENSEES                     #
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
import PLOT_2D as S04


# Define materials for nonlinear columns and beam
# Define parameters (units: mm, N)
# ------------------------------------------
# CONCRETE                  tag   f'c        ec0   f'cu        ecu
# Core concrete (confined)
fcC = -27.6         # [N/mm²] Concrete Compressive Strength
ec0C = -0.0045      # [mm/mm] Concrete Compressive Strain
fcUC = -21          # [N/mm²] Concrete Compressive Ultimate Strength
ecuC = -0.015       # [mm/mm] Concrete Compressive Ultimate Strain

# Cover concrete (unconfined)
fcU = -18           # [N/mm²] Concrete Compressive Strength
ec0U = -0.0025      # [mm/mm] Concrete Compressive Strain
fcUU = -2           # [N/mm²] Concrete Compressive Ultimate Strength
ecuU = -0.008       # [mm/mm] Concrete Compressive Ultimate Strain
 
# STEEL
# Reinforcing steel
fy = 4000           # [N/mm²] Steel Rebar Yield Strength   
Es = 2e5            # [N/mm²] Modulus of Elasticity
ey = fy/Es          # [mm/mm] Steel Rebar Yield Strain
fu = 1.1818*fy      # [N/mm²] Steel Rebar Ultimate Strength
esu = ey*75.2       # [mm/mm] Steel Rebar Ultimate Strain
Esh = (fu - fy)/(esu - ey)
Bs = Esh / Es

# Column Section
Bc = 500                 # [mm] Depth of the Section 
Hc = 500                 # [mm] Height of the Section  
coverC = 50              # [mm] Concrete Section Cover
DIAc = 25                # [mm] # Rebar Size Diameter
AsC = np.pi*(DIAc**2)/4   # [mm²] Area of Rebar

# Beam Section
Bb = 500                 # [mm] Depth of the Section 
Hb = 300                 # [mm] Height of the Section  
coverB = 50              # [mm] Concrete Section Cover
DIAb = 18                # [mm] # Rebar Size Diameter
AsB = np.pi*(DIAb**2)/4  # [mm²] Area of Rebar

# Apply the triangular distributed load
F_MAX = -0.375           # [N/mm] Maximum Distributed Load
F_INCR = -0.0001         # [N/mm] Force Increment for Distributed Load
LOAD = 'TRIANGULAR'      # Select type of distributed load

# Define Structure Length
LENGTH_COL = 3000        # [mm] Column Length 
LENGTH_BM = 7000         # [mm] Beam Length 

# Define Analysis Properties
MAX_ITERATIONS = 5000    # Convergence iteration for test
MAX_TOLERANCE = 1.0e-6   # Convergence tolerance for test
#STEEL_KIND: 1 -> WITHOUT HARDENING AND ULTIMATE STRAIN
#STEEL_KIND: 2 -> WITH HARDENING AND ULTIMATE STRAIN
#%%------------------------------------------------------------------------------
def PUSHOVER_ANALYSIS(LENGTH_COL, LENGTH_BM, F_MAX, F_INCR, STEEL_KIND):
    # Reset model
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)

    # CORNER NODES
    ops.node(1, 0.0, 0.0)
    ops.node(2, LENGTH_BM, 0.0)
    ops.node(3, 0.0, LENGTH_COL)
    ops.node(4, LENGTH_BM, LENGTH_COL)
    # BEAM NODES
    ops.node(5, 0.25*LENGTH_BM, LENGTH_COL)
    ops.node(6, 0.50*LENGTH_BM, LENGTH_COL)
    ops.node(7, 0.75*LENGTH_BM, LENGTH_COL)
    # RIGHT COLUMN NODES
    ops.node(8, 0.0, 0.25*LENGTH_COL)
    ops.node(9, 0.0, 0.50*LENGTH_COL)
    ops.node(10, 0.0, 0.75*LENGTH_COL)
    # LEFT COLUMN NODES
    ops.node(11, LENGTH_BM, 0.25*LENGTH_COL)
    ops.node(12, LENGTH_BM, 0.50*LENGTH_COL)
    ops.node(13, LENGTH_BM, 0.75*LENGTH_COL)
        
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


    ops.uniaxialMaterial('Concrete01', coreTag, fcC, ec0C, fcUC, ecuC)  # Core concrete (confined)
    ops.uniaxialMaterial('Concrete01', coverTag, fcU, ec0U, fcUU, ecuU) # Cover concrete (unconfined)
    
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
    ops.element('nonlinearBeamColumn', 1, 1, 8, numIntgrPts, secTagC, transfTag)  # COLUMN 01
    ops.element('nonlinearBeamColumn', 2, 8, 9, numIntgrPts, secTagC, transfTag)  # COLUMN 02
    ops.element('nonlinearBeamColumn', 3, 9, 10, numIntgrPts, secTagC, transfTag) # COLUMN 03
    ops.element('nonlinearBeamColumn', 4, 10, 3, numIntgrPts, secTagC, transfTag) # COLUMN 04
    
    ops.element('nonlinearBeamColumn', 5, 2, 11, numIntgrPts, secTagC, transfTag)  # COLUMN 05
    ops.element('nonlinearBeamColumn', 6, 11, 12, numIntgrPts, secTagC, transfTag) # COLUMN 06
    ops.element('nonlinearBeamColumn', 7, 12, 13, numIntgrPts, secTagC, transfTag) # COLUMN 07
    ops.element('nonlinearBeamColumn', 8, 13, 4, numIntgrPts, secTagC, transfTag)  # COLUMN 08
    
    ops.element('nonlinearBeamColumn', 9, 3, 5, numIntgrPts, secTagB, transfTag)  # BEAM 01
    ops.element('nonlinearBeamColumn', 10, 5, 6, numIntgrPts, secTagB, transfTag) # BEAM 02
    ops.element('nonlinearBeamColumn', 11, 6, 7, numIntgrPts, secTagB, transfTag) # BEAM 03
    ops.element('nonlinearBeamColumn', 12, 7, 4, numIntgrPts, secTagB, transfTag) # BEAM 04

    # Data storage
    FORCE_S, FORCE_A, MOMENT = [], [], []
    DISP_X, DISP_Y, ROT = [], [], []
    KA, KS, KI, STEP = [], [], [], []

    # Define time series and load pattern
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    #ops.load(3, 0.0, 0.0, 0.0)
    #ops.load(4, 0.0, 0.0, 0.0)

        
    # Define uniform distributed load analysis
    ops.system('BandGeneral')
    ops.constraints('Plain')
    ops.numberer('Plain')
    #ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/normDispIncr.html
    ops.algorithm('Newton')
    ops.integrator('LoadControl', F_INCR)
    ops.analysis('Static')
    
    # Total steps
    steps = int(np.abs(F_MAX)/np.abs(F_INCR))
    
    for step in range(steps):
        AA = np.abs(F_INCR) # Absolute Force Increment for Distributed Load
        if LOAD == 'TRIANGULAR':
            # Apply distributed load to all beams in the structure
            # eleLoad -ele $eleTag1 <$eleTag2 ....> -type -beamUniform $Wz <$Wx>
            ops.eleLoad('-ele', 9, '-type', '-beamUniform', step * AA * 0.25) 
            ops.eleLoad('-ele', 10, '-type', '-beamUniform', step * AA * 0.50)
            ops.eleLoad('-ele', 11, '-type', '-beamUniform', step * AA * 0.75)
            ops.eleLoad('-ele', 12, '-type', '-beamUniform', step * AA * 1.00)
        if LOAD == 'RECTANGULAR':
            # Apply distributed load to all beams in the structure
            # eleLoad -ele $eleTag1 <$eleTag2 ....> -type -beamUniform $Wz <$Wx>
            ops.eleLoad('-ele', 9, '-type', '-beamUniform', step*AA) 
            ops.eleLoad('-ele', 10, '-type', '-beamUniform', step*AA)
            ops.eleLoad('-ele', 11, '-type', '-beamUniform', step*AA)
            ops.eleLoad('-ele', 12, '-type', '-beamUniform', step*AA)

        OK = ops.analyze(1)
        S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS
        # Record results
        ops.reactions()
        S = ops.nodeReaction(1, 1) + ops.nodeReaction(2, 1) # SHEAR BASE REACTION
        A = ops.nodeReaction(1, 2) + ops.nodeReaction(2, 2) # AXIAL BASE REACTION
        M = ops.nodeReaction(1, 3) + ops.nodeReaction(2, 3) # MOMENT BASE REACTION
        #print(rot, M)
        disp_X = ops.nodeDisp(6, 1) # LATERAL DISPLACEMENT IN X FOR NODE 6
        disp_Y = ops.nodeDisp(6, 2) # LATERAL DISPLACEMENT IN Y FOR NODE 6
        rot = ops.nodeDisp(6, 3)    # ROTATION IN Z FOR NODE 6
        FORCE_S.append(S)
        FORCE_A.append(A)
        MOMENT.append(M)
        DISP_X.append(disp_X)
        DISP_Y.append(disp_Y)
        ROT.append(rot)
        KS.append(np.abs(S)/np.abs(disp_X)) # LATERAL STIFFNESS IN X
        KA.append(np.abs(A)/np.abs(disp_Y)) # LATERAL STIFFNESS IN Y
        KI.append(np.abs(M)/np.abs(rot))    # ROTATIONAL STIFFNESS IN Z
        STEP.append(step*AA)                # DISTRIBUTED LOAD 
        print(step+1, disp_X, S)
        
    #ops.wipe()    

    return FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y, ROT, KA, KS, KI, STEP

#%%------------------------------------------------------------------------------
# Analysis Durations:
starttime = TI.process_time()

# WITHOUT HARDENING AND ULTIMATE STRAIN
FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y, ROT, KA, KS, KI, STEP = PUSHOVER_ANALYSIS(LENGTH_COL, LENGTH_BM, F_MAX, F_INCR, STEEL_KIND=1)
# WITH HARDENING AND ULTIMATE STRAIN
FORCE_S02, FORCE_A02, MOMENT02, DISP_X02, DISP_Y02, ROT02, KA02, KS02, KI02, STEP02  = PUSHOVER_ANALYSIS(LENGTH_COL, LENGTH_BM, F_MAX, F_INCR, STEEL_KIND=2)

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
plt.xlabel('Displacement  in X [mm] in Middle of Beam')
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
plt.xlabel('Displacement in Y [mm] in Middle of Beam')
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
plt.xlabel('Distributed Load Steps [N/mm]')
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
plt.xlabel('Distributed Load Steps [N/mm]')
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
plt.xlabel('Distributed Load Steps [N/mm]')
plt.legend(['Steel01: WITHOUT HARDENING AND ULTIMATE STRAIN', 'Hysteretic: WITH HARDENING AND ULTIMATE STRAIN'])
plt.grid()
plt.show()

plt.figure(10, figsize=(12, 8))
plt.plot(STEP, DISP_X, color='brown', linewidth=2)
plt.plot(STEP02, DISP_X02, color='gold', linestyle='--', linewidth=2)
#plt.scatter(STEP, DISP_X, color='brown', linewidth=2)
#plt.scatter(STEP02, DISP_X02, color='gold', linestyle='--', linewidth=2)
plt.title('Displacement During the Analysis')
plt.ylabel('Displacement - X [mm] in Middle of Beam')
plt.xlabel('Distributed Load Steps [N/mm]')
plt.legend(['Steel01: WITHOUT HARDENING AND ULTIMATE STRAIN', 'Hysteretic: WITH HARDENING AND ULTIMATE STRAIN'])
plt.grid()
plt.show()

plt.figure(11, figsize=(12, 8))
plt.plot(STEP, DISP_Y, color='blue', linewidth=2)
plt.plot(STEP02, DISP_Y02, color='cyan', linestyle='--', linewidth=2)
#plt.scatter(STEP, DISP_X, color='brown', linewidth=2)
#plt.scatter(STEP02, DISP_X02, color='gold', linestyle='--', linewidth=2)
plt.title('Displacement During the Analysis')
plt.ylabel('Displacement - Y [mm] in Middle of Beam')
plt.xlabel('Distributed Load Steps [N/mm]')
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
plt.xlabel('Distributed Load Steps [N/mm]')
plt.legend(['Steel01: WITHOUT HARDENING AND ULTIMATE STRAIN', 'Hysteretic: WITH HARDENING AND ULTIMATE STRAIN'])
plt.grid()
plt.show()
#%%------------------------------------------------------------------------------
import BILINEAR_CURVE as BC

# ---------------------------------------
#  Plot BaseShear-Displacement Analysis
# ---------------------------------------
XX = np.abs(DISP_X02); YY = np.abs(FORCE_S02); # ABSOLUTE VALUE
SLOPE_NODE = 10

DATA = BC.BILNEAR_CURVE(XX, YY, SLOPE_NODE)
X, Y, Elastic_ST, Plastic_ST, Tangent_ST, Ductility_Rito, Over_Strength_Factor = DATA

XLABEL = 'Displacement in X [mm] in Middle of Beam'
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
XX = np.abs(DISP_Y02); YY = np.abs(FORCE_A02); # ABSOLUTE VALUE
SLOPE_NODE = 10

DATA = BC.BILNEAR_CURVE(XX, YY, SLOPE_NODE)
X, Y, Elastic_ST, Plastic_ST, Tangent_ST, Ductility_Rito, Over_Strength_Factor = DATA

XLABEL = 'Displacement in Y [mm] in Middle of Beam'
YLABEL = 'Base-Axial [N]'
LEGEND01 = 'Curve'
LEGEND02 = 'Bilinear Fitted'
LEGEND03 = 'Undefined'
TITLE = f'Last Data of BaseAxial-Displacement Analysis - Ductility Ratio: {X[2]/X[1]:.4f} - Over Strength Factor: {Y[2]/Y[1]:.4f}'
COLOR = 'black'
BC.PLOT_2D(np.abs(DISP_Y02), np.abs(FORCE_A02), X, Y, X, Y, XLABEL, YLABEL, TITLE, LEGEND01, LEGEND02, LEGEND03, COLOR='black', Z=2) 
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
# %% Plot 2D Frame Shapes
S04.PLOT_2D_FRAME(deformed_scale=1)  # Adjust scale factor as needed
#%%------------------------------------------------------------------------------   
# EXPORT DATA TO EXCEL
import pandas as pd
DATA_TOTAL = { 
    'Distributed_Load_Steps': STEP,# [N/mm] Uniform Distributed Loads
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
results_df.to_excel('CONCRETE_FRAME_INCREMENTAL_DISTRIBUTED_LOAD_RESULTS.xlsx', index=False) 
#%%------------------------------------------------------------------------------
# Print out the state of all nodes
#ops.printModel("node",1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)
# Print out the state of element 1 , 2 and 3
#ops.printModel("ele", 1, 2 , 3)
# Print the Model
ops.printModel()
ops.printModel("-JSON", "-file", "CONCRETE_FRAME_PROGRESSIVE_COLLAPSE.json")
#%%------------------------------------------------------------------------------ 

    
