###########################################################################################################
#                                             IN THE NAME OF ALLAH                                        #
#     REVERSED CYCLIC PUSHOVER ANALYSIS OF ECCENTRICALLY BARE BRACED STEEL FRAME WITH FRICTION DAMPER     #
#               EVALUATING STRAIN HARDENING AND ULTIMATE STRAIN CRITERIA USING OPENSEES                   #
#---------------------------------------------------------------------------------------------------------#
#                          THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                     #
#                                   EMAIL: salar.d.ghashghaei@gmail.com                                   #
###########################################################################################################
"""
1. The analysis compares nonlinear rotational behavior of steel beam-column
 elements under cyclic lateral displacements using OpenSees.
2. Two material models—*Steel01* (bilinear without degradation) and *Hysteretic*
 (tri-linear with pinching and strength/stiffness degradation)—are used.
3. Both models are subjected to identical loading protocols to investigate cyclic
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
# PAPER:Investigation of the Hysteresis Performance of Multi-Story Y-Shaped Eccentrically Bare Braced Steel Frame with Block Slit Damper (BSD)
'https://www.mdpi.com/2075-5309/15/3/451' 
# PAPER: A State-of-the-Art Review of Passive Energy Dissipation Systems in Steel Braces
'https://www.mdpi.com/2075-5309/13/4/851'
# PAPER: A modified friction damper for diagonal bracing of structures
'https://www.sciencedirect.com/science/article/pii/S0143974X13001077'
# PAPER: Experimental and Numerical Study on an Innovative Trapezoidal-Shaped Damper to Improve the Behavior of CBF Braces
'https://www.mdpi.com/2075-5309/13/1/140'
   
import openseespy.opensees as ops
import matplotlib.pyplot as plt
import numpy as np
import ANALYSIS_FUNCTION as S02
import STEEL_SECTION_FUN as S03
import PLOT_2D as S04

def CYCLIC_LOADING(num_cycles, samples_per_cycle, MAX_VALUE, EXPO=False, exponent=2):
    import numpy as np
    import matplotlib.pyplot as plt

    t_cycle = np.linspace(0, 2 * np.pi, samples_per_cycle)

    if not EXPO:
        # Linear amplitude growth capped to MAX_VALUE
        increment_linear = 1.0
        max_amplitude = num_cycles * increment_linear
        scale_factor = MAX_VALUE / max_amplitude

        signals_linear = []
        SINGAL = []
        last_end = 0.0
        for i in range(1, num_cycles + 1):
            amplitude = i * scale_factor
            signal = amplitude * np.sin(t_cycle)
            SINGAL.append(signal)
            signal_shifted = signal - last_end  # shift to continue from last value
            signals_linear.append(signal_shifted)
            last_end = signal
        signal_linear = np.concatenate(signals_linear)
        SINGALS = np.concatenate(SINGAL)
        time_linear = np.linspace(0, num_cycles, num_cycles * samples_per_cycle)

        # Plot
        plt.figure()
        plt.plot(time_linear, SINGALS, color='black')
        plt.title('Cyclic Curve with Linear Amplitude Increment')
        plt.xlabel('Cycle Number')
        plt.ylabel('Load or Displacement Amplitude')
        plt.grid(True)
        plt.tight_layout()
        plt.show()

        return time_linear, signal_linear

    else:
        # Polynomial amplitude growth capped to MAX_VALUE
        max_amplitude = (num_cycles) ** exponent
        scale_factor = MAX_VALUE / max_amplitude

        signals_poly = []
        SINGAL = []
        last_end = 0.0
        for i in range(1, num_cycles + 1):
            amplitude = (i ** exponent) * scale_factor
            signal = amplitude * np.sin(t_cycle)
            SINGAL.append(signal)
            signal_shifted = signal - last_end
            signals_poly.append(signal_shifted)
            last_end = signal
        signal_poly = np.concatenate(signals_poly)
        SINGALS = np.concatenate(SINGAL)
        time_poly = np.linspace(0, num_cycles, num_cycles * samples_per_cycle)

        # Plot
        plt.figure()
        plt.plot(time_poly, SINGALS, color='black')
        plt.title(f'Cyclic Curve with Polynomial Amplitude Increment (exponent={exponent})')
        plt.xlabel('Cycle Number')
        plt.ylabel('Load or Displacement Amplitude')
        plt.grid(True)
        plt.tight_layout()
        plt.show()

        return time_poly, signal_poly


    
# Parameters
num_cycles = 200
samples_per_cycle = 1000
MAX_VALUE = 7 # [mm] MAXIMUM DISPLACEMENT VALUE
TIME, AMP = CYCLIC_LOADING(num_cycles, samples_per_cycle, MAX_VALUE, EXPO=True, exponent=2) # Using this for Cyclic Displacement

#%%%------------------------------------------------------------------------------------
# Define materials for nonlinear columns
# Define parameters (units: mm, N)
# ------------------------------------------
# STEEL
# Reinforcing steel
fy = 400          # [N/mm²] Steel Rebar Yield Strength   
Es = 2e5          # [N/mm²] Modulus of Elasticity
ey = fy/Es        # [mm/mm] Steel Yield Strain
fu = 1.1818*fy    # [N/mm²] Steel Ultimate Strength
esu = ey*75.2     # [mm/mm] Steel Ultimate Strain
Esh = (fu - fy)/(esu - ey)
Bs = Esh / Es
# FRACTION DAMPER IN X DIRECTION
FYx = 520e4   # [N] Yield force in X direction
EYx = 1.1     # [mm] Yield displacement in X direction
FUx = 150e4   # [N] Ultimate force in X direction
EUx = 10.5    # [mm] Ultimate displacement in X direction
# FRACTION DAMPER IN Y DIRECTION
FYy = 602e3   # [N] Yield force in Y direction
EYy = 5.5     # [mm] Yield displacement in Y direction
FUy = 762e3   # [N] Ultimate force in Y direction
EUy = 50.3    # [mm] Ultimate displacement in Y direction

LENGTH_COL = 3000           # [mm] Column Length 
LENGTH_BM = 7000            # [mm] Beam Length 
# Define Analysis Properties
MAX_ITERATIONS = 5000      # Convergence iteration for test
MAX_TOLERANCE = 1.0e-8     # Convergence tolerance for test
#STEEL_KIND: 1 -> WITHOUT HARDENING AND ULTIMATE STRAIN
#STEEL_KIND: 2 -> WITH HARDENING AND ULTIMATE STRAIN
#%%------------------------------------------------------------------------------
def CYCLIC_ANALYSIS(LENGTH_COL, LENGTH_BM, CYCL, NUM, STEEL_KIND):
    # Reset model
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)

    # Nodes
    ops.node(1, 0.0, 0.0)
    ops.node(2, LENGTH_BM, 0.0)
    ops.node(3, 0.0, LENGTH_COL)
    ops.node(4, LENGTH_BM, LENGTH_COL)
    ops.node(5, 0.5*LENGTH_BM, LENGTH_COL)
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 1, 1, 1)

    secTagC = 10  # COLUMN
    secTagB = 20  # BEAM
    secTagBR = 30 # BRACE
    FDmatTagX = 40 # FRACTION DAMPER IN X DIRECTION
    FDmatTagY = 50 # FRACTION DAMPER IN Y DIRECTION
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
        
    
    # Pinching parameters in X direection
    pinchXFDx  = 0.1   # 10% pinching in displacement
    pinchYFDx  = 0.2   # 20% strength loss under pinching
    # Damage parameters (no cyclic degradation)
    damage1FDx = 0.0   # Strength deterioration
    damage2FDx = 0.0   # Stiffness deterioration
    # Fractional exponent (smoothness of backbone)
    betaFDx    = 0.2
    ops.uniaxialMaterial('Hysteretic', FDmatTagX, FYx, EYx, FUx, EUx, 0.2*FUx, 1.35*EUx,
                                                -FYx, -EYx, -FUx, -EUx, -0.32*FUx, -1.21*EUx,
                                                pinchXFDx, pinchYFDx, damage1FDx, damage2FDx, betaFDx)
    
    # Pinching parameters in Y direection
    pinchXFDy  = 0.1   # 10% pinching in displacement
    pinchYFDy  = 0.2   # 20% strength loss under pinching
    # Damage parameters (no cyclic degradation)
    damage1FDy = 0.0   # Strength deterioration
    damage2FDy = 0.0   # Stiffness deterioration
    # Fractional exponent (smoothness of backbone)
    betaFDy    = 0.2
    ops.uniaxialMaterial('Hysteretic', FDmatTagY, FYy, EYy, FUy, EUy, 0.2*FUy, 1.12*EUy,
                                                -FYy, -EYy, -FUy, -EUy, -0.32*FUy, -1.21*EUy,
                                                pinchXFDy, pinchYFDy, damage1FDy, damage2FDy, betaFDy)

    # STEEL I SECTION
    d, MASS = S03.STEEL_I_SECTION_QUAD(secTagB, steelTag, PLOT=True, DENSITY=7850/1e9)
    # GREEK CROSS SECTION WITH FLANGES
    #d, MASS = S03.GREEK_CROSS_SECTION_WITH_FLANGES(secTag, steelTag, PLOT=True, DENSITY=7850/1e9)

    # STEEL DOUBLE I SECTIONS WITH PLATE
    d, MASS = S03.DOUBLE_I_SECTION(secTagC, steelTag, PLOT=True, DENSITY=7850/1e9)
    #d, MASS = S03.BOX_SECTION_WITH_ANGLES_PLATES(secTagBR, steelTag, PLOT=True, DENSITY=7850/1e9)
    #d, MASS = S03.C_TUBE_SECTION(secTagBR, steelTag, PLOT=True, DENSITY=7850/1e9)
    d, MASS = S03.R_TUBE_SECTION_QUAD(secTagBR, steelTag, PLOT=True, DENSITY=7850/1e9)
    # 2 UNP SECTION WITH PLATES
    #d, MASS = S03.TWO_UNP_SECTION_WITH_TOP_BOTTOM_PLATES(secTagBR, steelTag, PLOT=True, DENSITY=7850/1e9, include_plates=False)
    
    # Define geometric transformation (corotational for large displacements)
    transfTag = 1
    #ops.geomTransf('Linear', transfTag)
    #ops.geomTransf('PDelta', transfTag)
    ops.geomTransf('Corotational', transfTag)
    numIntgrPts = 5
    ops.element('nonlinearBeamColumn', 1, 1, 3, numIntgrPts, secTagC, transfTag) # COLUMN 01
    ops.element('nonlinearBeamColumn', 2, 2, 4, numIntgrPts, secTagC, transfTag) # COLUMN 02
    ops.element('nonlinearBeamColumn', 3, 3, 5, numIntgrPts, secTagB, transfTag) # BEAM 01
    ops.element('nonlinearBeamColumn', 4, 4, 5, numIntgrPts, secTagB, transfTag) # BEAM 02
    #ops.element('truss', 5, 1, 5, secTagBR)                                      # BRACE 01
    #ops.element('truss', 6, 2, 5, secTagBR)                                      # BRACE 02
    ops.element('trussSection', 5, 1, 5, secTagBR)                               # BRACE 01
    ops.element('trussSection', 6, 2, 5, secTagBR)                               # BRACE 02
    #ops.element('corotTruss', 5, 1, 5, secTagBR)                                  # BRACE 01
    #ops.element('corotTruss', 6, 2, 5, secTagBR)                                  # BRACE 02
    # LINK INFO: https://opensees.berkeley.edu/wiki/index.php/Corotational_Truss_Element
    

    # FRACTIONAL DAMPER 
    # Define element (zeroLength element with hysteretic material)
    ops.element('zeroLength', 7, 5, 5, '-mat', FDmatTagX, FDmatTagY, '-dir', 1, 2)
    #ops.element('zeroLength', 7, 5, 5, '-mat', FDmatTag, '-dir', 1)
    #ops.element('zeroLength', 8, 5, 5, '-mat', FDmatTag, '-dir', 2)
    # LINK INFO: https://opensees.berkeley.edu/wiki/index.php/ZeroLength_Element


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
    steps = NUM

    # Use displacement control on rotational dof (dof 3 at node 2)
    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.system('BandGeneral')
    ops.algorithm('Newton')
    ops.analysis('Static')
    
    for step in range(steps):
        print(step+1, CYCL[step])
        ops.integrator('DisplacementControl', 3, 1, CYCL[step]) 
        ops.integrator('DisplacementControl', 4, 1, CYCL[step]) 
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
         
    #ops.wipe()

    return FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y, ROT, KA, KS, KI, STEP

#%%------------------------------------------------------------------------------
# WITHOUT HARDENING AND ULTIMATE STRAIN
FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y, ROT, KA, KS, KI, STEP = CYCLIC_ANALYSIS(LENGTH_COL, LENGTH_BM, CYCL=AMP, NUM=len(AMP), STEEL_KIND=1)
# WITH HARDENING AND ULTIMATE STRAIN
FORCE_S02, FORCE_A02, MOMENT02, DISP_X02, DISP_Y02, ROT02, KA02, KS02, KI02, STEP02  = CYCLIC_ANALYSIS(LENGTH_COL, LENGTH_BM, CYCL=AMP, NUM=len(AMP), STEEL_KIND=2)
#%%------------------------------------------------------------------------------
# Define piecewise linear envelope for loading and unloading
displacementsX = np.array([-1.21*EUx, -EUx, -EYx, 0, EYx, EUx, 1.12*EUx])
forcesX = np.array([-0.32*FUx, -FUx, -FYx, 0, FYx, FUx, 0.2*FUx])

displacementsY = np.array([-1.21*EUy, -EUy, -EYy, 0, EYy, EUy, 1.12*EUy])
forcesY = np.array([-0.32*FUy, -FUy, -FYy, 0, FYy, FUy, 0.2*FUy])



# Plotting
# Plot the hysteretic curve
plt.figure(0, figsize=(12, 5))
plt.plot(displacementsX, forcesX, 'blue', label="Hysteretic Model in X Dir.", markersize=4)
plt.plot(displacementsY, forcesY, 'cyan', label="Hysteretic Model in Y Dir.", markersize=4)
plt.axhline(0, color='gray', linestyle='--', linewidth=0.8)
plt.axvline(0, color='gray', linestyle='--', linewidth=0.8)
plt.xlabel("Displacement (m)")
plt.ylabel("Force (kN)")
plt.title("Force-Displacement Curve for Stiffness Spring - Hysteretic Uniaxial Material")
plt.legend()
plt.grid(True)
plt.show()
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
plt.xlabel('Displacement  in X [mm]')
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
plt.xlabel('Displacement in Y [mm]')
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
plt.title('ROTATIONAL STIFFNESS-LATERAL STIFFNESS IN X DIAGRAM')
plt.ylabel('Rotational Stiffness [N.mm/Rad]')
plt.xlabel('Lateral Stiffness in X [N/mm]')
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
plt.title('ROTATIONAL STIFFNESS-LATERAL STIFFNESS IN Y DIAGRAM')
plt.ylabel('Rotational Stiffness [N.mm/Rad]')
plt.xlabel('Lateral Stiffness in Y [N/mm]')
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
plt.ylabel('Displacement - X [mm]')
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
plt.ylabel('Displacement - Y [mm]')
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
# %% Plot 2D Frame Shapes
S04.PLOT_2D_FRAME(deformed_scale=100)  # Adjust scale factor as needed
#%%------------------------------------------------------------------------------   
# EXPORT DATA TO EXCEL
import pandas as pd
DATA_TOTAL = { 
    'DISP_X_WO': DISP_X,# WITHOUT HARDENING AND ULTIMATE STRAIN
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
    'LATERAL_S_ST_WO': KS,
    'LATERAL_S_ST_W': KS02,
    'LATERAL_A_ST_WO': KA,
    'LATERAL_A_ST_W': KA02,
    'ROTATIONAL_ST_WO': KI,
    'ROTATIONAL_ST_W': KI02,
}
# Convert to DataFrame
results_df = pd.DataFrame(DATA_TOTAL)
# Export the DataFrame to an Excel file
results_df.to_excel('STEEL_BRACED_FRAME_FRACTIONAL_DAMPER_CYCLIC_RESULTS.xlsx', index=False) 

    
