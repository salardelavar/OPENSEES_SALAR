###########################################################################################################
#                                             IN THE NAME OF ALLAH                                        #
#            REVERSED CYCLIC PUSHOVER ANALYSIS OF COMPOSITE COLUMN REINFORCED CONCRETE SECTION            #
#            WITH PLATES EVALUATING STRAIN HARDENING AND ULTIMATE STRAIN CRITERIA USING OPENSEES          #
#---------------------------------------------------------------------------------------------------------#
#                          THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                     #
#                                   EMAIL: salar.d.ghashghaei@gmail.com                                   #
###########################################################################################################
"""
1. Generates user-defined cyclic displacement or load histories (linear or polynomial growth)
 and plots the loading signal.
 
2. Defines confined (core) and unconfined (cover) concrete via `Concrete01`,
 plus rebar and plate steel (bilinear or hysteretic with pinching/degradation).
 
3. Builds a fiberized rectangular section using `S03.CONCRETE_CONFINED_REC_PLATE_SECTION_QUAD`,
 combining concrete and steel fibers.
 
4. Creates a column element between two nodes in a 2D, 3-DOF model, with 
rotation at the free node as the controlled DOF.

5. Applies displacement‐control on rotation by feeding the cyclic curvature array into
 OpenSees’s `DisplacementControl` integrator.
 
6. At each step, runs a static Newton analysis, checks convergence via `S02.ANALYSIS`,
 and records reactions, displacements, rotations, and instantaneous stiffness.
 
7. Executes two scenarios: (1) bilinear steel without hardening, (2) hysteretic steel
 with strain-hardening and degradation.
 
8. Produces hysteresis and interaction plots: P–M, axial vs. displacement, moment–rotation,
 stiffness vs. rigidity, and time histories of all response quantities.
 
9. Aggregates all results into a pandas DataFrame and exports them to `COMPOSITE_CYCLIC_RESULTS.xlsx`
 for post-processing.
 
10. Offers a modular, extensible framework for advanced cyclic SDOF RC section studies,
 capturing confinement, pinching, and degradation effects.

"""
import openseespy.opensees as ops
import matplotlib.pyplot as plt
import numpy as np
import ANALYSIS_FUNCTION as S02
import COMPOSITE_CONCRETE_SECTION_FUN as S03

def CYCLIC_LOADING(num_cycles, samples_per_cycle, Amplitude_facor=1e-2, EXPO=False, exponent=2):
    import numpy as np
    import matplotlib.pyplot as plt
    
    increment_linear = 1.0  # linear increment per cycle
    #exponent = 2  # polynomial exponent for amplitude growth
    
    # Time vectors for each cycle
    t_cycle = np.linspace(0, 2 * np.pi, samples_per_cycle)
    
    if EXPO == False:
        # Generate cyclic loading with linear increment
        signals_linear = []
        for i in range(1, num_cycles + 1):
            amplitude = i * increment_linear * Amplitude_facor
            signals_linear.append(amplitude * np.sin(t_cycle))
        signal_linear = np.concatenate(signals_linear)
        time_linear = np.linspace(0, num_cycles, num_cycles * samples_per_cycle)
        # Plot: Linear Increment
        plt.figure()
        plt.plot(time_linear, signal_linear, color='black')
        plt.title('Cyclic Loading with Linear Amplitude Increment')
        plt.xlabel('Cycle Number')
        plt.ylabel('Load Amplitude')
        plt.grid(True)
        plt.tight_layout()
        plt.show()
        AAA = time_linear, signal_linear
        
    
    if EXPO == True:
        # Generate cyclic loading with polynomial increment
        signals_poly = []
        for i in range(1, num_cycles + 1):
            amplitude =  Amplitude_facor * i ** exponent
            signals_poly.append(amplitude * np.sin(t_cycle))
        signal_poly = np.concatenate(signals_poly)
        time_poly = np.linspace(0, num_cycles, num_cycles * samples_per_cycle)
        # Plot: Polynomial Increment
        plt.figure()
        plt.plot(time_poly, signal_poly, color='black')
        plt.title(f'Cyclic Loading with Polynomial Amplitude Increment (exponent={exponent})')
        plt.xlabel('Cycle Number')
        plt.ylabel('Load Amplitude')
        plt.grid(True)
        plt.tight_layout()
        plt.show()
        AAA = time_poly, signal_poly
    
    return AAA
    
# Parameters
num_cycles = 60
samples_per_cycle = 500
TIME, AMP = CYCLIC_LOADING(num_cycles, samples_per_cycle) # Using this for Cyclic Displacement

# Define materials for nonlinear columns
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
fy = 400          # [N/mm²] Steel Rebar Yield Strength   
Es = 2e5          # [N/mm²] Modulus of Elasticity
ey = fy/Es        # [mm/mm] Steel Rebar Yield Strain
fu = 1.1818*fy    # [N/mm²] Steel Rebar Ultimate Strength
esu = ey*75.2     # [mm/mm] Steel Rebar Ultimate Strain
Esh = (fu - fy)/(esu - ey)
Bs = Esh / Es

# Steel Plate
fyP = 240          # [N/mm²] Steel Rebar Yield Strength   
EsP = 2e5          # [N/mm²] Modulus of Elasticity
eyP = fyP/EsP      # [mm/mm] Steel Rebar Yield Strain
fuP = 1.1818*fyP   # [N/mm²] Steel Rebar Ultimate Strength
esuP = eyP*41.8    # [mm/mm] Steel Rebar Ultimate Strain
EshP = (fuP - fyP)/(esuP - eyP)
BsP = EshP / EsP

LENGTH = 4000           # [mm] Column Length 
B = 400                 # [mm] Depth of the Section 
H = 500                 # [mm] Height of the Section  
cover = 50              # [mm] Concrete Section Cover
DIA = 25                # [mm] # Rebar Size Diameter
As = np.pi*(DIA**2)/4   # [mm²] Area of Rebar
plateThickness = 10     # [mm] Plate Thickness
# Define Analysis Properties
MAX_ITERATIONS = 5000      # Convergence iteration for test
MAX_TOLERANCE = 1.0e-10    # Convergence tolerance for test
#STEEL_KIND: 1 -> WITHOUT HARDENING AND ULTIMATE STRAIN
#STEEL_KIND: 2 -> WITH HARDENING AND ULTIMATE STRAIN
#%%------------------------------------------------------------------------------
def CYCLIC_ANALYSIS(LENGTH, b, h, cover, Es, fy, As, ECU, CYCL, NUM, STEEL_KIND, num_cycles=3, amp_curv=0.000001):
    # Invert ultimate strain for OpenSees convention
    ECU_op = -ECU

    # Reset model
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)

    # Nodes
    ops.node(1, 0.0, 0.0)
    ops.node(2, 0.0, LENGTH)
    ops.fix(1, 1, 1, 1)

    secTag = 10
    coreTag = 1
    coverTag = 2
    steelTag = 3
    steelPlateTag = 4
    numBarsTop, barAreaTop = 5, np.pi *(18**2)/4
    numBarsBot, barAreaBot = 5, np.pi *(20**2)/4
    numBarsIntTot, barAreaInt = 4, np.pi *(5**2)/4
    
    if STEEL_KIND == 1:# WITHOUT HARDENING AND ULTIMATE STRAIN
        ops.uniaxialMaterial('Steel01', steelTag, fy, Es, 0.0) 
        ops.uniaxialMaterial('Steel01', steelPlateTag, fyP, EsP, 0.0) 
    if STEEL_KIND == 2:# WITH HARDENING AND ULTIMATE STRAIN    
        pinchX = 0.8   # Pinching factor in X direction
        pinchY = 0.5   # Pinching factor in Y direction
        damage1 = 0.0  # Damage due to ductility
        damage2 = 0.0  # Damage due to energy
        beta = 0.1 # Stiffness degradation parameter
        ops.uniaxialMaterial('Hysteretic', steelTag, fy, ey, fu, esu, 0.2*fu, 1.1*esu, -fy, -ey, -fu, -esu, -0.2*fu, -1.1*esu, pinchX, pinchY, damage1, damage2, beta)
        
        ops.uniaxialMaterial('Hysteretic', steelPlateTag , fyP, eyP, fuP, esuP, 0.2*fuP, 1.1*esuP, -fyP, -eyP, -fuP, -esuP, -0.2*fuP, -1.1*esuP, pinchX, pinchY, damage1, damage2, beta)
        # INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Hysteretic_Material


    ops.uniaxialMaterial('Concrete01', coreTag, fcC, ec0C, fcUC, ecuC)  # Core concrete (confined)
    ops.uniaxialMaterial('Concrete01', coverTag, fcU, ec0U, fcUU, ecuU) # Cover concrete (unconfined)

    d, MASS = S03.CONCRETE_CONFINED_REC_PLATE_SECTION_QUAD(secTag, H, B, cover, cover, coreTag, coverTag, steelTag, steelPlateTag,
                                          numBarsTop, barAreaTop, numBarsBot, barAreaBot, numBarsIntTot, barAreaInt,
                                          plateThickness, PLOT=True, CONCRETE_DENSITY=2500/1e9, STEEL_DENSITY=7850/1e9)

    # Define geometric transformation (corotational for large displacements)
    transfTag = 1
    #ops.geomTransf('Linear', transfTag)
    #ops.geomTransf('PDelta', transfTag)
    ops.geomTransf('Corotational', transfTag)
    numIntgrPts = 5
    ops.element('nonlinearBeamColumn', 1, 1, 2, numIntgrPts, secTag, transfTag)

    # Data storage
    FORCE_S, FORCE_A, MOMENT = [], [], []
    DISP_X, DISP_Y, ROT = [], [], []
    KA, KS, KI, STEP = [], [], [], []

    # Define time series and load pattern
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(2, 1.0, -1.0, 0.0)

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
        ops.integrator('DisplacementControl', 2, 1, CYCL[step]) 
        OK = ops.analyze(1)
        S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS
        # Record results
        ops.reactions()
        S = ops.nodeReaction(1, 1) # SHEAR BASE REACTION
        A = ops.nodeReaction(1, 2) # AXIAL BASE REACTION
        M = ops.nodeReaction(1, 3) # MOMENT BASE REACTION
        #print(rot, M)
        disp_X = ops.nodeDisp(2, 1) # LATERAL DISPLACEMENT IN X FOR NODE 2
        disp_Y = ops.nodeDisp(2, 1) # LATERAL DISPLACEMENT IN Y FOR NODE 2
        rot = ops.nodeDisp(2, 3)    # ROTATION IN Z FOR NODE 2
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

    return FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y, ROT, KA, KS, KI, STEP

#%%------------------------------------------------------------------------------
# WITHOUT HARDENING AND ULTIMATE STRAIN
FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y, ROT, KA, KS, KI, STEP = CYCLIC_ANALYSIS(LENGTH, B, H, cover, Es, fy, As, esu, CYCL=AMP, NUM=len(AMP), STEEL_KIND=1)
# WITH HARDENING AND ULTIMATE STRAIN
FORCE_S02, FORCE_A02, MOMENT02, DISP_X02, DISP_Y02, ROT02, KA02, KS02, KI02, STEP02  = CYCLIC_ANALYSIS(LENGTH, B, H, cover, Es, fy, As, esu, CYCL=AMP, NUM=len(AMP), STEEL_KIND=2)

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
plt.title('ROTATIONAL STIFFNESS-LATERAL STIFFNESS DIAGRAM')
plt.ylabel('Rotational Stiffness [N.mm/Rad]')
plt.xlabel('Lateral Stiffness [N/mm]')
plt.legend(['Steel01: WITHOUT HARDENING AND ULTIMATE STRAIN', 'Hysteretic: WITH HARDENING AND ULTIMATE STRAIN'])
plt.semilogx()
plt.semilogy()
plt.grid()
plt.show()

plt.figure(6, figsize=(12, 8))
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

plt.figure(7, figsize=(12, 8))
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

plt.figure(8, figsize=(12, 8))
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

plt.figure(9, figsize=(12, 8))
plt.plot(STEP, DISP_X, color='brown', linewidth=2)
plt.plot(STEP02, DISP_X02, color='gold', linestyle='--', linewidth=2)
#plt.scatter(STEP, DISP_X, color='brown', linewidth=2)
#plt.scatter(STEP02, DISP_X02, color='gold', linestyle='--', linewidth=2)
plt.title('Displacement During the Analysis')
plt.ylabel('Displacement - X [m]')
plt.xlabel('Steps')
plt.legend(['Steel01: WITHOUT HARDENING AND ULTIMATE STRAIN', 'Hysteretic: WITH HARDENING AND ULTIMATE STRAIN'])
plt.grid()
plt.show()

plt.figure(10, figsize=(12, 8))
plt.plot(STEP, DISP_Y, color='blue', linewidth=2)
plt.plot(STEP02, DISP_Y02, color='cyan', linestyle='--', linewidth=2)
#plt.scatter(STEP, DISP_X, color='brown', linewidth=2)
#plt.scatter(STEP02, DISP_X02, color='gold', linestyle='--', linewidth=2)
plt.title('Displacement During the Analysis')
plt.ylabel('Displacement - Y [m]')
plt.xlabel('Steps')
plt.legend(['Steel01: WITHOUT HARDENING AND ULTIMATE STRAIN', 'Hysteretic: WITH HARDENING AND ULTIMATE STRAIN'])
plt.grid()
plt.show()

plt.figure(11, figsize=(12, 8))
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
    'ROTATIONAL_ST_WO': KI,
    'ROTATIONAL_ST_W': KI02,
}
# Convert to DataFrame
results_df = pd.DataFrame(DATA_TOTAL)
# Export the DataFrame to an Excel file
results_df.to_excel('COMPOSITE_COLUMN_CYCLIC_RESULTS.xlsx', index=False) 

    
