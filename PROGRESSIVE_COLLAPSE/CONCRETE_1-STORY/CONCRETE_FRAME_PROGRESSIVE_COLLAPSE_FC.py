################################################################################################################
#                                                  IN THE NAME OF ALLAH                                        #
#                                 PROGRESSIVE COLLAPSE ANALYSIS OF CONCRETE FRAME.                             #
#                  EVALUATING STRAIN HARDENING AND ULTIMATE STRAIN CRITERIA USING OPENSEES                     #
#--------------------------------------------------------------------------------------------------------------#
#                                                     LOAD CONTROL                                            #
#--------------------------------------------------------------------------------------------------------------#
#                             THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                       #
#                                        EMAIL: salar.d.ghashghaei@gmail.com                                   #
################################################################################################################
"""
[1] The analysis compares nonlinear rotational behavior of concrete beam-column
 elements under pushover lateral displacements using OpenSees.
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
------------------------------------------------------------------------------------- 
 Progressive collapse of reinforced concrete frames occurs when a local failure—due to accidental
  actions such as impact, explosion or fire—triggers a chain reaction of element removals, leading
  to partial or total structural loss. Advanced assessment hinges on capturing nonlinear material
  behavior, geometric effects, and load‐redistribution mechanisms that dictate whether alternative
  load paths can sustain the imposed demands.

 [1] Modeling Philosophy:
 - Fiber based sections discretize concrete and steel across the cross-section, enabling accurate
  stress–strain representation under combined axial, bending and shear demands. Cover, core concrete,
  and rebar layouts are modeled with uniaxial constitutive laws that include confinement, cracking,
  strain hardening and ultimate strain limits.
 - Nonlinear beam–column elements employ Gauss integration points along member length, paired with
  corotational kinematics to account for large displacements and P-Δ effects in a fully consistent 2D formulation.

 [2] Analysis Strategy:
 - Alternate load‐path method: deliberately remove one or more columns (or beams) after applying
  gravity loads, then trace the static response under incremental displacement control at a critical
  location. The structural response captures bending yielding, shear failure, catenary action and
  eventual loss of load‐bearing capacity.
 - Pushover framework: displacement control at a predefined “attack” node (e.g., mid-height of a key column)
  simulates the increasing drift demands after element removal. Reaction forces at the base yield
  a capacity curve relating force vs. displacement, from which reserve strength and ductility can be assessed.

 [3] Key Response Mechanisms:
 - Flexural yielding and plastic hinge formation in adjacent beams and columns allow moment redistribution.
  Hinge rotation capacity depends on reinforcement ratio, concrete confinement and strain‐hardening characteristics of steel.
 - P-Δ instability magnifies demands when large drifts develop; corotational transforms ensure equilibrium
  accounts for geometric nonlinearity.
 - Catenary action engages once flexural capacity is exhausted and members deform significantly, mobilizing
  tensile forces in reinforcement. Accurate modeling of ultimate tendon strain (eult) is critical to predict post-peak response.
 - Shear failure remains brittle; its prevention through detailing (stirrups, confinement) is vital to allow
  ductile mechanisms to develop.

 [4] Collapse Criteria and Robustness:
 - Vertical and lateral drift limits define collapse thresholds. Exceeding the ultimate drift at a control
  node triggers element deletion, simulating fracture or buckling.
 - Progressive removal tests on different locations probe system robustness, verifying that the structure
  retains sufficient redundancy and alternative load paths.

 [5] Practical Implications:
 - Design against progressive collapse requires enforcing continuity (tie forces), detailing for ductility
  (strong-column–weak-beam hierarchy), and redundancy (multiple load paths).
 - Nonlinear analyses—both static pushover and dynamic removal simulations—inform code provisions
  (e.g., UFC 4-023-03, GSA Guidelines) by quantifying reserve strength margins and post-failure behavior.

 In workflow, the combination of fiber section definitions, corotational beam–column elements,
  displacement-controlled pushover and element removal routines captures the critical phases of
  progressive collapse: initial yielding, redistribution, catenary action and final instability.
  Interpreting capacity curves and post-peak degradation provides insights on member detailing and
  overall frame robustness under accidental collapse scenarios.
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
fy = 400            # [N/mm²] Steel Rebar Yield Strength   
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

# Apply the nodal Point load
PX, PY, MZ = 1000000, -100000, 0.0

# Apply the triangular distributed load
F_MAX = -0.152           # [N/mm] Maximum Distributed Load
F_INCR = -0.0001         # [N/mm] Force Increment for Distributed Load

Collapse_Disp = -10      # [mm] Absolute Value Collapse Vertical Displacement

# Define Structure Length
LENGTH_COL = 3000        # [mm] Column Length 
LENGTH_BM = 7000         # [mm] Beam Length 

# Define Analysis Properties
MAX_ITERATIONS = 5000    # Convergence iteration for test
MAX_TOLERANCE = 1.0e-10  # Convergence tolerance for test

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
    # LEFT COLUMN ELEMENTS
    ops.element('nonlinearBeamColumn', 1, 1, 8, numIntgrPts, secTagC, transfTag)  # COLUMN 01
    ops.element('nonlinearBeamColumn', 2, 8, 9, numIntgrPts, secTagC, transfTag)  # COLUMN 02
    ops.element('nonlinearBeamColumn', 3, 9, 10, numIntgrPts, secTagC, transfTag) # COLUMN 03
    ops.element('nonlinearBeamColumn', 4, 10, 3, numIntgrPts, secTagC, transfTag) # COLUMN 04
    # RIGHT COLUMN ELEMENTS
    ops.element('nonlinearBeamColumn', 5, 2, 11, numIntgrPts, secTagC, transfTag)  # COLUMN 05
    ops.element('nonlinearBeamColumn', 6, 11, 12, numIntgrPts, secTagC, transfTag) # COLUMN 06
    ops.element('nonlinearBeamColumn', 7, 12, 13, numIntgrPts, secTagC, transfTag) # COLUMN 07
    ops.element('nonlinearBeamColumn', 8, 13, 4, numIntgrPts, secTagC, transfTag)  # COLUMN 08
    # BEAM ELEMENTS
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
    #ops.load(3, PX, PY, MZ)
    ops.load(4, PX, PY, MZ)

        
    # Define uniform distributed load analysis
    ops.system('BandGeneral')
    ops.constraints('Plain')
    ops.numberer('Plain')
    #ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS, 2) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/normDispIncr.html
    ops.algorithm('Newton')
    ops.integrator('LoadControl', F_INCR) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/loadControl.html
    ops.analysis('Static')
    
    # Total steps
    steps = int(np.abs(F_MAX/F_INCR))
    delete_element = False
    
    for step in range(steps):
        AA = np.abs(F_INCR) # Absolute Force Increment for Distributed Load
        
        # Apply the triangular distributed load to all beams in the structure
        # eleLoad -ele $eleTag1 <$eleTag2 ....> -type -beamUniform $Wz <$Wx>
        ops.eleLoad('-ele', 9, '-type', '-beamUniform', step * AA * 0.25, 0.0) 
        ops.eleLoad('-ele', 10, '-type', '-beamUniform', step * AA * 0.50, 0.0)
        ops.eleLoad('-ele', 11, '-type', '-beamUniform', step * AA * 0.75, 0.0)
        ops.eleLoad('-ele', 12, '-type', '-beamUniform', step * AA * 1.00, 0.0)
        
        OK = ops.analyze(1)
        S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS
        disp_X = ops.nodeDisp(6, 1) # LATERAL DISPLACEMENT IN X FOR NODE 6
        disp_Y = ops.nodeDisp(6, 2) # LATERAL DISPLACEMENT IN Y FOR NODE 6
        rot = ops.nodeDisp(6, 3)    # ROTATION IN Z FOR NODE 6
        ops.reactions()
        if abs(disp_Y) < abs(Collapse_Disp): # DISPLACEMENT CONTROL FOR COLUMN 3 UNTIL COLUMN IS NOT DELETED
            S = ops.nodeReaction(1, 1) + ops.nodeReaction(2, 1) # SHEAR BASE REACTION
            A = ops.nodeReaction(1, 2) + ops.nodeReaction(2, 2) # AXIAL BASE REACTION
            M = ops.nodeReaction(1, 3) + ops.nodeReaction(2, 3) # MOMENT BASE REACTION
        if abs(disp_Y) == abs(Collapse_Disp) and not delete_element: # DISPLACEMENT CONTROL FOR COLUMN 3 AND DELETE ELEMENT
            print(f"Displacement exceeds {Collapse_Disp} mm. Removing right column at step {step + 1}.") 
            ops.remove('element', 9)  # REMOVE ELEMENT 9
            ops.remove('element', 10) # REMOVE ELEMENT 10
            ops.remove('element', 11) # REMOVE ELEMENT 11
            ops.remove('element', 12) # REMOVE ELEMENT 12
            ops.remove('sp', 5)       # REMOVE NODE 5
            ops.remove('sp', 6)       # REMOVE NODE 6
            ops.remove('sp', 7)       # REMOVE NODE 7
            # LINK INFO: https://openseespydoc.readthedocs.io/en/stable/src/remove.html
            ops.domainChange()    # update the solver without element
            # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/domainChange.html
            ops.integrator('LoadControl', -F_INCR) # SO, AFTER THE COLUMN COLLAPSES, WE APPLY A VERTICAL DISPLACEMENT INCREMENT
            delete_element = True 
            S = ops.nodeReaction(1, 1) + ops.nodeReaction(2, 1) # SHEAR BASE REACTION
            A = ops.nodeReaction(1, 2) + ops.nodeReaction(2, 2) # AXIAL BASE REACTION
            M = ops.nodeReaction(1, 3) + ops.nodeReaction(2, 3) # MOMENT BASE REACTION
        if abs(disp_Y) > abs(Collapse_Disp): # DISPLACEMENT CONTROL FOR RIGHT COLUMN UNTIL COLUMN IS NOT DELETED
            ops.remove('element', 9)  # REMOVE ELEMENT 9
            ops.remove('element', 10) # REMOVE ELEMENT 10
            ops.remove('element', 11) # REMOVE ELEMENT 11
            ops.remove('element', 12) # REMOVE ELEMENT 12
            ops.remove('sp', 5)       # REMOVE NODE 5
            ops.remove('sp', 6)       # REMOVE NODE 6
            ops.remove('sp', 7)       # REMOVE NODE 7
            # LINK INFO: https://openseespydoc.readthedocs.io/en/stable/src/remove.html
            S = ops.nodeReaction(1, 1) + ops.nodeReaction(2, 1) # SHEAR BASE REACTION
            A = ops.nodeReaction(1, 2) + ops.nodeReaction(2, 2) # AXIAL BASE REACTION
            M = ops.nodeReaction(1, 3) + ops.nodeReaction(2, 3) # MOMENT BASE REACTION
        FORCE_S.append(S)
        FORCE_A.append(A)
        MOMENT.append(M)
        DISP_X.append(disp_X)
        DISP_Y.append(disp_Y)
        ROT.append(rot)
        KS.append(np.abs(S/disp_X)) # LATERAL STIFFNESS IN X
        KA.append(np.abs(A/disp_Y)) # LATERAL STIFFNESS IN Y
        KI.append(np.abs(M/rot))    # ROTATIONAL STIFFNESS IN Z
        STEP.append(step*AA)        # DISTRIBUTED LOAD 
        print(f'{step+1} {disp_X:.5f} {S:.5f}')
        
    #ops.wipe()    

    return FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y, ROT, KA, KS, KI, STEP

#%%------------------------------------------------------------------------------
# Analysis Durations:
starttime = TI.process_time()

FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y, ROT, KA, KS, KI, STEP = PUSHOVER_ANALYSIS(LENGTH_COL, LENGTH_BM, F_MAX, F_INCR, STEEL_KIND=2)

totaltime = TI.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n') 
#%%------------------------------------------------------------------------------
plt.figure(1, figsize=(12, 8))
plt.plot(MOMENT, FORCE_A, color='black')
#plt.scatter(MOMENT, FORCE_A, color='blue', linewidth=2)
plt.title('P-M Interaction')
plt.ylabel('Axial Force [N]')
plt.xlabel('Bending Moment [N.mm]')
plt.grid()
plt.show()

plt.figure(2, figsize=(12, 8))
plt.plot(DISP_X, FORCE_S, color='green', linewidth=2)
#plt.scatter(DISP_X, FORCE_S, color='purple', linewidth=2)
plt.title('SHEAR FORCE-DISPLACEMENT DIAGRAM')
plt.ylabel('Shear Force [N]')
plt.xlabel('Displacement  in X [mm] in Middle of Beam')
plt.grid()
plt.show()

plt.figure(3, figsize=(12, 8))
plt.plot(DISP_Y, FORCE_A, color='purple', linewidth=2)
#plt.scatter(DISP_Y, FORCE_A, color='purple', linewidth=2)
plt.title('AXIAL FORCE-DISPLACEMENT DIAGRAM')
plt.ylabel('Axial Force [N]')
plt.xlabel('Displacement in Y [mm] in Middle of Beam')
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
#plt.plot(KI, KS, color='black', linewidth=2)
plt.scatter(KI, KS, color='black', linewidth=2)
plt.title('ROTATIONAL STIFFNESS-LATERAL STIFFNESS DIAGRAM')
plt.ylabel('Rotational Stiffness [N.mm/Rad]')
plt.xlabel('Lateral Stiffness in X Dir. [N/mm]')
plt.semilogx()
plt.semilogy()
plt.grid()
plt.show()

plt.figure(6, figsize=(12, 8))
#plt.plot(KI, KA, color='black', linewidth=2)
plt.scatter(KI, KA, color='black', linewidth=2)
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
plt.xlabel('Distributed Load Steps [N/mm]')
plt.grid()
plt.show()

plt.figure(8, figsize=(12, 8))
plt.plot(STEP, FORCE_S, color='purple', linewidth=2)
#plt.scatter(STEP, FORCE_S, color='purple', linewidth=2)
plt.title('Shear Force During the Analysis')
plt.ylabel('Shear Force [N]')
plt.xlabel('Distributed Load Steps [N/mm]')
plt.grid()
plt.show()

plt.figure(9, figsize=(12, 8))
plt.plot(STEP, MOMENT, color='green', linewidth=2)
#plt.scatter(STEP, MOMENT, color='green', linewidth=2)
plt.title('Moment During the Analysis')
plt.ylabel('Moment [kN.mm]')
plt.xlabel('Distributed Load Steps [N/mm]')
plt.grid()
plt.show()

plt.figure(10, figsize=(12, 8))
plt.plot(STEP, DISP_X, color='brown', linewidth=2)
#plt.scatter(STEP, DISP_X, color='brown', linewidth=2)
plt.title('Displacement During the Analysis')
plt.ylabel('Displacement - X [mm] in Middle of Beam')
plt.xlabel('Distributed Load Steps [N/mm]')
plt.grid()
plt.show()

plt.figure(11, figsize=(12, 8))
plt.plot(STEP, DISP_Y, color='blue', linewidth=2)
#plt.scatter(STEP, DISP_X, color='brown', linewidth=2)
plt.title('Displacement During the Analysis')
plt.ylabel('Displacement - Y [mm] in Middle of Beam')
plt.xlabel('Distributed Load Steps [N/mm]')
plt.grid()
plt.show()

plt.figure(12, figsize=(12, 8))
plt.plot(STEP, ROT, color='black', linewidth=2)
#plt.scatter(STEP, ROT, color='black', linewidth=2)
plt.title('Rotation During the Analysis')
plt.ylabel('Rotation [rad]')
plt.xlabel('Distributed Load Steps [N/mm]')
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
XX = np.abs(DISP_Y); YY = np.abs(FORCE_A); # ABSOLUTE VALUE
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
# %% Plot 2D Frame Shapes
S04.PLOT_2D_FRAME(deformed_scale=10)  # Adjust scale factor as needed
#%%------------------------------------------------------------------------------   
# EXPORT DATA TO EXCEL
import pandas as pd
DATA_TOTAL = { 
    'Distributed_Load_Steps': STEP,# [N/mm] Uniform Distributed Loads
    'DISP_X_WO': DISP_X, 
    'DISP_Y_WO': DISP_Y,
    'ROTATION_WO': ROT,
    'AXIAL_FORCE_WO': FORCE_A,
    'SHEAR_FORCE_WO': FORCE_S,
    'MOMENT_WO': MOMENT,
    'AXIAL_RIGIDITY_WO': np.abs(FORCE_A),
    'ROTATIONAL_ST_WO': KI,
    'LATERAL_ST_Y_WO': KA,
    'LATERAL_ST_X_WO': KS,
}
# Convert to DataFrame
results_df = pd.DataFrame(DATA_TOTAL)
# Export the DataFrame to an Excel file
results_df.to_excel('CONCRETE_FRAME_PROGRESSIVE_COLLAPSE_RESULTS.xlsx', index=False) 
#%%------------------------------------------------------------------------------
# Print out the state of all nodes
#ops.printModel("node",1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)
# Print out the state of element 1 , 2 and 3
#ops.printModel("ele", 1, 2 , 3)
# Print the Model
ops.printModel()
ops.printModel("-JSON", "-file", "CONCRETE_FRAME_PROGRESSIVE_COLLAPSE.json")
#%%------------------------------------------------------------------------------ 

    


