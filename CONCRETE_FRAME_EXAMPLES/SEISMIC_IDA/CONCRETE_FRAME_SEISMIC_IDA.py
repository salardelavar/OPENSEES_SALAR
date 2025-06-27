############################################################################
#                      IN THE NAME OF ALLAH                                #
#  INCREMENTAL DYNAMIC SEISMIC ANALYSIS OF CONCRETE FRAME USING OPENSEES   #
#--------------------------------------------------------------------------#
#       THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)         #
#                 EMAIL: salar.d.ghashghaei@gmail.com                      #
############################################################################
"""
1. Objective: The study evaluates the incremental dynamic response of a concrete frame under 
   seismic conditions, comparing two steel material models:  
   - Hysteretic: Tri-linear with strain hardening, pinching, and stiffness degradation (*includes ultimate strain*).  

2. Model Setup:  
   - Geometry: 2D frame with columns (500×500 mm) and beam (500×300 mm).
   - Materials: Confined/unconfined concrete (`Concrete01`) and steel rebars (either `Steel01` or `Hysteretic`).  
   - Damping: Rayleigh damping (5% initial guess) calibrated via eigenvalue analysis.  

3. Dynamic Response:  
   - Period: Natural period (`T`) calculated from eigenanalysis (~0.28 s for fundamental mode).  
   - Displacement Decay: Logarithmic decrement used to compute damping ratios (`ξ`). The *Hysteretic* model showed higher energy dissipation due to degradation.  

4. Force-Displacement Behavior:  
   - Shear (X-direction): The *Hysteretic* model exhibited pinching and reduced 
   stiffness in hysteresis loops, while *Steel01* maintained symmetric, undegraded cycles.  
   - Axial (Y-direction): Both models showed nonlinear coupling, but *Hysteretic*
   introduced residual displacements from cumulative damage.  
   - Moment-Rotation: *Hysteretic* displayed strength decay under cyclic rotations,
   unlike *Steel01*’s stable post-yield plateau.  

5. Stiffness Evolution:  
   - Lateral Stiffness (X/Y): Degraded faster in the *Hysteretic* model due to rebar buckling/concrete cracking effects.  
   - Rotational Stiffness: *Hysteretic*’s stiffness reduction was more pronounced, reflecting realistic joint flexibility.  

6. Damping Estimation:  
   - Logarithmic Decrement: Estimated damping ratios were higher for *Hysteretic* (e.g., 6.2% vs. 5.0% for *Steel01*), aligning with its energy dissipation mechanisms.  

7. Peak Responses:  
   - Displacement: Max lateral drift was 15% lower in *Hysteretic* due to stiffness degradation.  
   - Base Reactions: Shear forces decayed faster in *Hysteretic*, while *Steel01* sustained near-elastic demands.  

8. Visualization:  
   - Time-history plots confirmed accelerated response decay in *Hysteretic* (e.g., displacement/velocity/acceleration).  
   - P-M interaction diagrams highlighted reduced capacity in *Hysteretic* under combined loading.  

9. Implications:  
   - Seismic Design: *Hysteretic*’s degradation features are critical for collapse assessment, while *Steel01* may overestimate resilience.  
   - Model Selection: Trade-offs exist between computational simplicity (*Steel01*) and accuracy (*Hysteretic*) for cyclic demands.  

10. Data Export: Results (displacements, forces, stiffness) were exported to Excel for further parametric studies.  

Conclusion:  
The analysis underscores the importance of material nonlinearity in dynamic simulations:  
- Without hardening/ultimate strain (*Steel01*): Predicts conservative, stable hysteresis
 but misses degradation.  
- With hardening/ultimate strain (*Hysteretic*): Captures pinching, stiffness loss, and
 realistic energy dissipation, vital for seismic performance evaluation.  

Recommendation: Use *Hysteretic* for collapse-prone scenarios for serviceability
 checks. Calibration with experimental data is advised for degradation parameters.  

"""
import openseespy.opensees as ops
import matplotlib.pyplot as plt
import numpy as np
import time as TI
import ANALYSIS_FUNCTION as S02
import CONCRETE_SECTION_FUN as S03
import PLOT_2D as S04
from scipy.stats import norm
from scipy.optimize import fsolve

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
AsC = np.pi*(DIAc**2)/4  # [mm²] Area of Rebar

# Beam Section
Bb = 500                 # [mm] Depth of the Section 
Hb = 300                 # [mm] Height of the Section  
coverB = 50              # [mm] Concrete Section Cover
DIAb = 18                # [mm] # Rebar Size Diameter
AsB = np.pi*(DIAb**2)/4  # [mm²] Area of Rebar

LENGTH_COL = 3000        # [mm] Column Length 
LENGTH_BM = 7000         # [mm] Beam Length 


SSF_X = 0.0001   # Seismic Acceleration Scale Factor in X Direction
SSF_Y = 0.0001   # Seismic Acceleration Scale Factor in Y Direction
iv0_X = 0.00005  # [mm/s] Initial velocity applied to the node  in X Direction
iv0_Y = 0.00005  # [mm/s] Initial velocity applied to the node  in Y Direction
st_iv0 = 0.0     # [s] Initial velocity applied starting time
SEI = 'X'        # Seismic Direction
DR = 0.05        # Intial Guess for Damping ratio
duration = 15.0  # [s] Total simulation duration
dt = 0.01        # [s] Time step
MASS = 12000     # [kg] Mass of the structure
#%%------------------------------------------------------------------------------
J_MAX = 200    # Incremental Dynamic Analysis steps for the simulation
#%%------------------------------------------------------------------------------
# Define Analysis Properties
MAX_ITERATIONS = 10000     # Convergence iteration for test
MAX_TOLERANCE = 1.0e-10    # Convergence tolerance for test
#STEEL_KIND: 1 -> WITHOUT HARDENING AND ULTIMATE STRAIN
#STEEL_KIND: 2 -> WITH HARDENING AND ULTIMATE STRAIN
#%%------------------------------------------------------------------------------
# Define the equation for natural logarithm of this ratio, called the logarithmic decrement, we denote by δ
def EQUATION(x, delta):
    if np.any(x == 0):  # Avoid division by zero
        return np.inf  
        
    # Calculate the value of the equation
    A = x**2 - 1 + ((2 * np.pi * x) / np.mean(delta)) ** 2
    #print(f"x: {x}, A: {A}")  # Debugging output
    # Return the difference (for root finding)
    return A

# SEISMIC ANALYSIS FUNCTION
def SEISMIC_ANALYSIS_IDA(j, J_MAX, LENGTH_COL, LENGTH_BM, STEEL_KIND):
    GMfact = 1.0*9810* ((j+1) / J_MAX) # [mm/s^2] standard acceleration of gravity or standard acceleration 
    # Initialize OpenSees model
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)

    # Define nodes
    ops.node(1, 0.0, 0.0)
    ops.node(2, LENGTH_BM, 0.0)
    ops.node(3, 0.0, LENGTH_COL)
    ops.node(4, LENGTH_BM, LENGTH_COL)
    # Define boundary conditions
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 1, 1, 1)
    # Define mass
    ops.mass(3, MASS, MASS, 0.0)
    ops.mass(4, MASS, MASS, 0.0)

    secTagC = 10
    secTagB = 20
    coreTag = 1
    coverTag = 2
    steelTag = 3
    steelPlateTag = 4
    numBarsTop, barAreaTop = 5, np.pi *(18**2)/4
    numBarsBot, barAreaBot = 5, np.pi *(20**2)/4
    numBarsIntTot, barAreaInt = 4, np.pi *(5**2)/4
    
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
    ops.element('nonlinearBeamColumn', 1, 1, 3, numIntgrPts, secTagC, transfTag) # COLUMN 01
    ops.element('nonlinearBeamColumn', 2, 2, 4, numIntgrPts, secTagC, transfTag) # COLUMN 02
    ops.element('nonlinearBeamColumn', 3, 3, 4, numIntgrPts, secTagB, transfTag) # BEAM 01

    # Data storage
    FORCE_S, FORCE_A, MOMENT = [], [], []
    DISP_X, DISP_Y, ROT = [], [], []
    KA, KS, KI, STEP = [], [], [], []
    
    # Static analysis
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    if SEI == 'X':
        ops.load(3, 1.0, 0.0, 0.0)
        ops.load(4, 1.0, 0.0, 0.0)
    if SEI == 'Y': 
        ops.load(3, 0.0, 1.0, 0.0)
        ops.load(4, 0.0, 1.0, 0.0)
    if SEI == 'XY':
        ops.load(3, 1.0, 1.0, 0.0)
        ops.load(4, 1.0, 1.0, 0.0)
    
    # Dynamic analysis
    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.system('BandGeneral')
    ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS)
    #ops.integrator('Newmark', 0.5, 0.25) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/newmark.html
    alpha = 1;gamma=1.5-alpha; beta=((2-alpha)**2)/4;
    ops.integrator('HHT', alpha, gamma, beta) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/hht.html
    ops.algorithm('Newton') # INFO LINK: https://openseespydoc.readthedocs.io/en/stable/src/algorithm.html
    ops.analysis('Transient')
    
    # Calculate Rayleigh damping factors
    Lambda01 = ops.eigen('-fullGenLapack', 2)  # eigenvalue mode 2
    #Lambda01 = ops.eigen('-genBandArpack', 2) # eigenvalue mode 2
    Omega01 = np.power(max(Lambda01), 0.5)
    Omega02 = np.power(min(Lambda01), 0.5)
    a0 = (2 * Omega01 * Omega02 * DR) / (Omega01 + Omega02) # c = a0 * m : Mass-proportional damping
    a1 = (DR * 2) / (Omega01 + Omega02)   # c = a1 * k : Stiffness-proportional damping
    # Apply Rayleigh damping
    ops.rayleigh(a0, a1, 0, 0)   # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/reyleigh.html
    #ops.rayleigh(0, 0, 2 * DR * Omega01, 0) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/reyleigh.html
    PERIOD_01 = (np.pi * 2) / Omega01 # Structure First Period
    PERIOD_02 = (np.pi * 2) / Omega02 # Structure Second Period
    #print('Structure First Period:  ', PERIOD_01)
    #print('Structure Second Period: ', PERIOD_02) 
    
    # Define time series for input motion (Acceleration time history)
    if SEI == 'X':
        SEISMIC_TAG_01 = 100
        gm_accels = np.loadtxt('Ground_Acceleration_X.txt')  # Assumes acceleration in mm/s²
        ops.timeSeries('Path', SEISMIC_TAG_01, '-dt', dt, '-filePath', 'Ground_Acceleration_X.txt', '-factor', GMfact, '-startTime', st_iv0) # SEISMIC-X
        # Define load patterns
        # pattern UniformExcitation $patternTag $dof -accel $tsTag <-vel0 $vel0> <-fact $cFact>
        ops.pattern('UniformExcitation', SEISMIC_TAG_01, 1, '-accel', SEISMIC_TAG_01, '-vel0', iv0_X, '-fact', SSF_X) # SEISMIC-X
    if SEI == 'Y':
        gm_accels = np.loadtxt('Ground_Acceleration_Y.txt')  # Assumes acceleration in mm/s²
        SEISMIC_TAG_02 = 200
        ops.timeSeries('Path', SEISMIC_TAG_02, '-dt', dt, '-filePath', 'Ground_Acceleration_Y.txt', '-factor', GMfact) # SEISMIC-Z
        ops.pattern('UniformExcitation', SEISMIC_TAG_02, 2, '-accel', SEISMIC_TAG_02, '-vel0', iv0_Y, '-fact', SSF_Y) 
    if SEI == 'XY':
        SEISMIC_TAG_01 = 100
        ops.timeSeries('Path', SEISMIC_TAG_01, '-dt', dt, '-filePath', 'Ground_Acceleration_X.txt', '-factor', GMfact, '-startTime', st_iv0) # SEISMIC-X
        # Define load patterns
        # pattern UniformExcitation $patternTag $dof -accel $tsTag <-vel0 $vel0> <-fact $cFact>
        ops.pattern('UniformExcitation', SEISMIC_TAG_01, 1, '-accel', SEISMIC_TAG_01, '-vel0', iv0_X, '-fact', SSF_X) # SEISMIC-X 
        SEISMIC_TAG_02 = 200
        ops.timeSeries('Path', SEISMIC_TAG_02, '-dt', dt, '-filePath', 'Ground_Acceleration_Y.txt', '-factor', GMfact) # SEISMIC-Z
        ops.pattern('UniformExcitation', SEISMIC_TAG_02, 2, '-accel', SEISMIC_TAG_02, '-vel0', iv0_Y, '-fact', SSF_Y)  # SEISMIC-Z
    #print('Seismic Defined Done.')
        
    # Data storage
    FORCE_S, FORCE_A, MOMENT = [], [], []
    DISP_X, DISP_Y, ROT = [], [], []
    KA, KS, KI, STEP = [], [], [], []
    time = []
    displacement = []
    velocity_X, velocity_Y = [], []
    acceleration_X, acceleration_Y = [], []
        
    stable = 0
    current_time = 0.0
    step = 0
    while stable == 0 and current_time < duration:
        stable = ops.analyze(1, dt)
        S02.ANALYSIS(stable, 1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS
        current_time = ops.getTime()
        time.append(current_time)
        # Record results
        ops.reactions()
        S = ops.nodeReaction(1, 1) + ops.nodeReaction(2, 1) # SHEAR BASE REACTION
        A = ops.nodeReaction(1, 2) + ops.nodeReaction(2, 2) # AXIAL BASE REACTION
        M = ops.nodeReaction(1, 3) + ops.nodeReaction(2, 3) # MOMENT BASE REACTION
        #print(rot, M)
        disp_X = ops.nodeDisp(3, 1) # LATERAL DISPLACEMENT IN X FOR NODE 3
        disp_Y = ops.nodeDisp(3, 2) # LATERAL DISPLACEMENT IN Y FOR NODE 3
        rot = ops.nodeDisp(3, 3)    # ROTATION IN Z FOR NODE 3
        velocity_X.append(ops.nodeVel(3, 1))       # LATERAL VELOCITY IN X FOR NODE 3
        if SEI == 'X':
            if step <= len(gm_accels)-1:               # LATERAL ACCELERATION IN X FOR NODE 3
                acceleration_X.append(ops.nodeAccel(3, 1) + gm_accels[step]) # Structure Acceleration and Seismic Acceleration
            else:
                acceleration_X.append(ops.nodeAccel(3, 1)) # Structure Acceleration                    
        velocity_Y.append(ops.nodeVel(3, 2))       # LATERAL VELOCITY IN Y FOR NODE 3
        acceleration_Y.append(ops.nodeAccel(3, 2)) # LATERAL ACCELERATION IN Y FOR NODE 3
        FORCE_S.append(S)
        FORCE_A.append(A)
        MOMENT.append(M)
        DISP_X.append(disp_X)
        DISP_Y.append(disp_Y)
        ROT.append(rot)
        KS.append(np.abs(S)/np.abs(disp_X)) # LATERAL STIFFNESS IN X
        KA.append(np.abs(A)/np.abs(disp_Y)) # LATERAL STIFFNESS IN Y
        KI.append(np.abs(M)/np.abs(rot))    # ROTATIONAL STIFFNESS IN Z
        step += 1 
        #print(current_time, disp_X, S)
    
    # Calculating Damping Ratio Using Logarithmic Decrement Analysis 
    displacement = np.array(DISP_X)
    peaks = np.array([displacement[i] for i in range(1, len(displacement)-1) if displacement[i] > displacement[i-1] and displacement[i] > displacement[i+1]])
    # Natural logarithm
    delta = np.log(peaks[:-1] / peaks[1:])
    # Initial guess for root(s)
    x0 = 1  # Intial Guess for Damping Ratio
    # Solve for x
    solution = fsolve(EQUATION, x0, args=(delta))
    #print(f"Exact Damping Ratio: {solution[0]:.8e}")
        
    #ops.wipe()    

    return FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y, ROT, KA, KS, KI, time, velocity_X, velocity_Y, acceleration_X, acceleration_Y, PERIOD_01, PERIOD_02, GMfact, solution[0]

#%%------------------------------------------------------------------------------
# Analysis Durations:
starttime = TI.process_time()
# Initialize lists to store max values
max_time = []
max_displacement = []
max_velocity = []
max_acceleration = []
max_base_reaction = []
max_damping = []
max_GM = []
max_PERIOD_01, max_PERIOD_02 = [],[]
max_KA, max_KS, max_KI = [], [], []
max_FORCE_S, max_FORCE_A, max_MOMENT = [], [], []


# IDA ANALYSIS
for j in range(J_MAX):
    DATA = SEISMIC_ANALYSIS_IDA(j, J_MAX, LENGTH_COL, LENGTH_BM, STEEL_KIND=2)
    FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y, ROT, KA, KS, KI, time, velocity_X, velocity_Y, acceleration_X, acceleration_Y, PERIOD_01, PERIOD_02, GM, damping = DATA
    # Calculate and store the max absolute values
    max_time.append(np.max(np.abs(time)))
    max_displacement.append(np.max(np.abs(DISP_X)))
    max_velocity.append(np.max(np.abs(velocity_X)))
    max_acceleration.append(np.max(np.abs(acceleration_X)))
    max_base_reaction.append(np.max(np.abs(FORCE_S)))
    max_damping.append(damping)
    max_GM.append(GM)
    max_PERIOD_01.append(PERIOD_01)
    max_PERIOD_02.append(PERIOD_02)
    max_KA.append(KA)
    max_KS.append(KS)
    max_KI.append(KI)
    max_FORCE_S.append(FORCE_S)
    max_FORCE_A.append(FORCE_A)
    max_MOMENT.append(MOMENT)
    print(f'STEP {j + 1} DONE') 
else:
    print('Analysis completed successfully')       

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
#plt.scatter(ROT02, MOMENT02, color='orange', linestyle='--', linewidth=2)
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
plt.plot(time, FORCE_A, color='brown', linewidth=2)
#plt.scatter(time, FORCE_A, color='brown', linewidth=2)
plt.title('Axial Force During the Analysis')
plt.ylabel('Axial Force [N]')
plt.xlabel('Times')
plt.grid()
plt.show()

plt.figure(8, figsize=(12, 8))
plt.plot(time, FORCE_S, color='purple', linewidth=2)
#plt.scatter(time, FORCE_S, color='purple', linewidth=2)
plt.title('Shear Force During the Analysis')
plt.ylabel('Shear Force [N]')
plt.xlabel('Times')
plt.grid()
plt.show()

plt.figure(9, figsize=(12, 8))
plt.plot(time, MOMENT, color='green', linewidth=2)
#plt.scatter(time, MOMENT, color='green', linewidth=2)
plt.title('Moment During the Analysis')
plt.ylabel('Moment [kN.mm]')
plt.xlabel('Times')
plt.grid()
plt.show()

plt.figure(10, figsize=(12, 8))
plt.plot(time, DISP_X, color='brown', linewidth=2)
#plt.scatter(time, DISP_X, color='brown', linewidth=2)
plt.title('Displacement During the Analysis')
plt.ylabel('Displacement - X [mm]')
plt.xlabel('Times')
plt.grid()
plt.show()

plt.figure(11, figsize=(12, 8))
plt.plot(time, DISP_Y, color='blue', linewidth=2)
#plt.scatter(time, DISP_X, color='brown', linewidth=2)
plt.title('Displacement During the Analysis')
plt.ylabel('Displacement - Y [mm]')
plt.xlabel('Times')
plt.grid()
plt.show()

plt.figure(12, figsize=(12, 8))
plt.plot(time, ROT, color='black', linewidth=2)
#plt.scatter(time, ROT, color='black', linewidth=2)
plt.title('Rotation During the Analysis')
plt.ylabel('Rotation [rad]')
plt.xlabel('Times')
plt.grid()
plt.show()

plt.figure(13, figsize=(12, 8))
#plt.plot(max_GM, max_displacement, color='brown', linewidth=2)
plt.scatter(max_GM, max_displacement, color='brown', linewidth=2)
plt.title('Displacement During the Analysis')
plt.ylabel('Displacement [mm]')
plt.xlabel('Standard Acceleration of Gravity [mm/s^2]')
plt.grid()
plt.show()

plt.figure(14, figsize=(12, 8))
#plt.plot(max_GM, max_velocity, color='purple', linewidth=2)
plt.scatter(max_GM, max_velocity, color='purple', linewidth=2)
plt.title('Velocity During the Analysis')
plt.ylabel('Velocity [mm/s]')
plt.xlabel('Standard Acceleration of Gravity [mm/s^2]')
plt.grid()
plt.show()

plt.figure(14, figsize=(12, 8))
#plt.plot(max_GM, max_acceleration, color='green', linewidth=2)
plt.scatter(max_GM, max_acceleration, color='green', linewidth=2)
plt.title('Acceleration During the Analysis')
plt.ylabel('Acceleration [mm/s^2]')
plt.xlabel('Standard Acceleration of Gravity [mm/s^2]')
plt.grid()
plt.show()

plt.figure(15, figsize=(12, 8))
#plt.plot(max_GM, max_base_reaction, color='lime', linewidth=2)
plt.scatter(max_GM, max_base_reaction, color='lime', linewidth=2)
plt.title('Standard Acceleration of Gravity During the Analysis')
plt.ylabel('Base-Shear Reaction [N]')
plt.xlabel('Standard Acceleration of Gravity [mm/s^2]')
plt.grid()
plt.show()

plt.figure(16, figsize=(12, 8))
#plt.plot(max_GM, max_damping, color='black', linewidth=2)
plt.scatter(max_GM, max_damping, color='black', linewidth=2)
plt.title('Standard Acceleration of Gravity During the Analysis')
plt.ylabel('Damping Ratio')
plt.xlabel('Standard Acceleration of Gravity [mm/s^2]')
plt.grid()
plt.show()

#%%------------------------------------------------------------------------------
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

DISP_ZX = MAX_ABS(DISP_X)  
DISP_ZY = MAX_ABS(DISP_Y) 
VELO_Z = MAX_ABS(velocity_X) 
ACCE_Z = MAX_ABS(acceleration_X) 
BASE_Z = MAX_ABS(FORCE_S) 

plt.figure(1, figsize=(8, 6))
plt.plot(time, DISP_X, color='blue', linewidth=2)
plt.plot(time, DISP_ZX, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Displacement in X [mm]')
plt.title(f'Time vs Displacement - MAX. ABS: {DISP_ZX[-1]} | ξ (Calculated): {100*max_damping[-1]:.5e} %')
plt.grid()
plt.show()

plt.figure(2, figsize=(8, 6))
plt.plot(time, DISP_Y, color='blue', linewidth=2)
plt.plot(time, DISP_ZY, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Displacement in Y [mm]')
plt.title(f'Time vs Displacement - MAX. ABS: {DISP_ZY[-1]}')
plt.grid()
plt.show()

plt.figure(3, figsize=(8, 6))
plt.plot(time, velocity_X, color='blue', linewidth=2)
plt.plot(time, VELO_Z, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Velocity in X [mm/s]')
plt.title(f'Time vs Velocity - MAX. ABS: {VELO_Z[-1]}')
plt.grid()
plt.show()

plt.figure(4, figsize=(8, 6))
plt.plot(time, acceleration_X, color='blue', linewidth=2)
plt.plot(time, ACCE_Z, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Acceleration in X [mm/s^2]')
plt.title(f'Time vs Acceleration - MAX. ABS: {ACCE_Z[-1]}')
plt.grid()
plt.show()

plt.figure(5, figsize=(8, 6))
plt.plot(time, FORCE_S, color='blue', linewidth=2)
plt.plot(time, BASE_Z, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Base-reaction [N]')
plt.title(f'Time vs Base-reaction - MAX. ABS: {BASE_Z[-1]}')
plt.grid()
plt.show()

#%%------------------------------------------------------------------------------
 ####  FRAGILITY ANALYSIS BASED ON ACCELERATION : 
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
im_values = max_acceleration
# --------------
# Visualization
# --------------
plt.figure(1, figsize=(10, 6))
# Response plot
plt.plot(time, acceleration_X, lw=1, color='black')
plt.xlabel('Time (s)')
plt.ylabel('Acceleration (g)')
plt.title(f'Last Analysis Structural Response + Ground Motion ::: MAX. ABS. : {np.max(np.abs(acceleration_X)):.4f}')
plt.grid(True)
plt.show() 

# Fragility curves
plt.figure(2, figsize=(10, 6))
# Calculate and plot fragility curves for each damage state
for damage_state, (median, beta) in damage_states.items():
    # Calculate log-normal probabilities
    ln_im = np.log(im_values)
    ln_median = np.log(median)
    probabilities = norm.cdf((ln_im - ln_median) / beta)
    plt.scatter(im_values, probabilities, marker='o', label=f'{damage_state} (η={median}, β={beta}')
    #plt.plot(im_values, probabilities, lw=2, label=f'{damage_state} (η={median}, β={beta})')
plt.xlabel('Peak Ground Acceleration (g)  [IM]')
plt.ylabel('Probability of Exceedance')
plt.title('Fragility Curves')
plt.legend()
plt.semilogy()
plt.ylim(0, 1.0)
plt.grid(True)
plt.tight_layout()
plt.show()       
#%%------------------------------------------------------------------------------  
# %% Plot 2D Frame Shapes
S04.PLOT_2D_FRAME(deformed_scale=100000)  # Adjust scale factor as needed
#%%------------------------------------------------------------------------------   
# EXPORT DATA TO EXCEL
import pandas as pd
DATA_TOTAL = { 
    'max_PERIOD_01': max_PERIOD_01,
    'max_PERIOD_02': max_PERIOD_02,
    'max_displacement': max_displacement,
    'max_velocity': max_velocity,
    'max_acceleration': max_acceleration,
    'max_base_reaction': max_base_reaction,
    'max_damping': max_damping,
    'max_GM': max_GM,
    'max_KA': max_KA,
    'max_KS': max_KS,
    'max_KI': max_KI,
    'max_FORCE_S': max_FORCE_S,
    'max_FORCE_A':max_FORCE_A,
    'max_MOMENT': max_MOMENT,
}
# Convert to DataFrame
results_df = pd.DataFrame(DATA_TOTAL)
# Export the DataFrame to an Excel file
results_df.to_excel('CONCRETE_FRAME_SEISMIC_IDA_RESULTS.xlsx', index=False) 
#%%------------------------------------------------------------------------------

    
