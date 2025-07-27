######################################################################################################################
#                          >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                          #
#      OPTIMIZATION OF COLUMN SECTION DEPTH BASED ON TARGET STRUCTURAL PERIOD FOR CONCRETE FRAME USING OPENSEES      #
#--------------------------------------------------------------------------------------------------------------------#
#                                  OPTIMIZATION ALGORITHM: NEWTON-RAPHSON METHOD                                     #
#--------------------------------------------------------------------------------------------------------------------#
#                              THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                            #
#                                       EMAIL: salar.d.ghashghaei@gmail.com                                          #
######################################################################################################################
"""
1. This Python script performs nonlinear dynamic and pushover analysis of a 2D reinforced concrete frame using OpenSeesPy.  
2. The code implements a Newton-Raphson optimization algorithm to determine optimal column depth based on target structural period.  
3. Advanced material models include confined/unconfined concrete (Concrete01) and hysteretic steel reinforcement (Hysteretic/Steel01).  
4. The frame consists of two nonlinear beam-column elements (columns) and one beam element with distributed plasticity.  
5. Dynamic analysis features Rayleigh damping calibrated using the first two modal frequencies.  
6. The script can perform both displacement-controlled pushover and time-history analyses under seismic excitation.  
7. Key outputs include force-displacement responses, stiffness degradation, and period-dependent behavior.  
8. The analysis tracks base reactions, displacements, rotations, and section stiffnesses throughout loading.  
9. Advanced features include corotational geometric transformations for large displacements.  
10. The code implements logarithmic decrement analysis to estimate damping ratios from free vibration decay.  
11. Comprehensive plotting functions visualize hysteresis loops, P-M interactions, and time-history responses.  
12. The optimization process automatically adjusts column depth to match a target fundamental period.  
13. Results are exported to Excel for post-processing, including peak responses and cumulative maxima.  
14. The model includes concentrated nodal masses and accounts for element mass contributions.  
15. Both stiffness-proportional and mass-proportional damping components are considered.  
16. The script handles convergence issues with automatic analysis restart capabilities.  
17. Statistical analysis of results includes histograms and box plots of key response parameters.  
18. The implementation demonstrates best practices for nonlinear frame analysis in OpenSees.  
19. The code is well-commented with references to OpenSees documentation and theory manuals.  
20. This provides a robust template for performance-based seismic design optimization of concrete frames.  

The script represents a sophisticated integration of structural dynamics theory, numerical optimization,
 and nonlinear finite element analysis for seismic design applications.
"""
#%%------------------------------------------------------------------------------
import openseespy.opensees as ops
import matplotlib.pyplot as plt
import numpy as np
import time as TI
import ANALYSIS_FUNCTION as S02
import CONCRETE_SECTION_FUN as S03
import PLOT_2D as S04
import SALAR_MATH as S05

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
fy = 4000         # [N/mm²] Steel Rebar Yield Strength   
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

#%% DEFINE THE ELEMENTS LENGTH
LENGTH_COL = 3000        # [mm] Column Length 
LENGTH_BM = 7000         # [mm] Beam Length 

#%% DEFINE PARAMEETRS FOR NONLINEAR DYNAMIC ANALYSIS
GMfact = 9810    # [mm/s²]standard acceleration of gravity or standard acceleration
SSF_X = 1.0      # Seismic Acceleration Scale Factor in X Direction
SSF_Y = 1.0      # Seismic Acceleration Scale Factor in Y Direction
iv0_X = 0.0005   # [mm/s] Initial velocity applied to the node  in X Direction
iv0_Y = 0.0005   # [mm/s] Initial velocity applied to the node  in Y Direction
st_iv0 = 0.0     # [s] Initial velocity applied starting time
SEI = 'X'        # Seismic Direction
DR = 0.05        # Intial Guess for Damping ratio
duration = 15.0  # [s] Total simulation duration
dt = 0.01        # [s] Time step
MASS = 12000     # [kg] Mass on the each column

#%% DEFINE PARAMEETRS FOR NONLINEAR STATIC ANALYSIS 
DMAX = 675       # [mm] Maximum Displacement
DINCR = 0.05     # [mm] Incremental Displacement

#%% DEFINE ANALYSIS PROPERTIES
MAX_ITERATIONS = 5000    # Convergence iteration for test
MAX_TOLERANCE = 1.0e-6   # Convergence tolerance for test
#STEEL_KIND: 1 -> WITHOUT HARDENING AND ULTIMATE STRAIN
#STEEL_KIND: 2 -> WITH HARDENING AND ULTIMATE STRAIN
#%%------------------------------------------------------------------------------
def PD_ANALYSIS(Hc, STEEL_KIND, ANA_KIND):
    Hb = Hc - 100       # [mm] Beam Height of the Section 
    
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
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/nonlinearBeamColumn.html
    ops.element('nonlinearBeamColumn', 1, 1, 3, numIntgrPts, secTagC, transfTag, '-mass', Bc*Hc*0.000025) # COLUMN 01
    ops.element('nonlinearBeamColumn', 2, 2, 4, numIntgrPts, secTagC, transfTag, '-mass', Bc*Hc*0.000025) # COLUMN 02
    ops.element('nonlinearBeamColumn', 3, 3, 4, numIntgrPts, secTagB, transfTag, '-mass', Bb*Hb*0.000025) # BEAM 01
    
    ELE_MASS = 0.000025*((Bc*Hc*2) + (Bb*Hb*0.000025))
    
    if ANA_KIND == 'PUSHOVER':
        WEIGHT = MASS * GMfact
        # Data storage
        FORCE_S, FORCE_A, MOMENT = [], [], []
        DISP_X, DISP_Y, ROT = [], [], []
        KA, KS, KI, STEP = [], [], [], []
    
        # Define time series and load pattern
        ops.timeSeries('Linear', 1)
        ops.pattern('Plain', 1, 1)
        #ops.load(3, 1.0, -WEIGHT, 0.0)
        #ops.load(4, 1.0, -WEIGHT, 0.0)
        ops.load(3, 1.0, -1, 0.0)
        ops.load(4, 1.0, -1, 0.0)
        
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
            print(step+1, disp_X, S)
        return FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y, ROT, KA, KS, KI, STEP        
    
    if ANA_KIND == 'DYNAMIC':
        # Define mass
        ops.mass(3, MASS, MASS, 0.0)
        ops.mass(4, MASS, MASS, 0.0)
        
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
        print('Structure First Period:  ', PERIOD_01)
        print('Structure Second Period: ', PERIOD_02) 
        
        # Define time series for input motion (Acceleration time history)
        if SEI == 'X':
            SEISMIC_TAG_01 = 100
            ops.timeSeries('Path', SEISMIC_TAG_01, '-dt', dt, '-filePath', 'Ground_Acceleration_X.txt', '-factor', GMfact, '-startTime', st_iv0) # SEISMIC-X
            # Define load patterns
            # pattern UniformExcitation $patternTag $dof -accel $tsTag <-vel0 $vel0> <-fact $cFact>
            ops.pattern('UniformExcitation', SEISMIC_TAG_01, 1, '-accel', SEISMIC_TAG_01, '-vel0', iv0_X, '-fact', SSF_X) # SEISMIC-X
        if SEI == 'Y':
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
        OMEGA, PERIOD = [], []
            
        stable = 0
        current_time = 0.0
        
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
            disp_X = ops.nodeDisp(3, 1)                # LATERAL DISPLACEMENT IN X FOR NODE 3
            disp_Y = ops.nodeDisp(3, 2)                # LATERAL DISPLACEMENT IN Y FOR NODE 3
            rot = ops.nodeDisp(3, 3)                   # ROTATION IN Z FOR NODE 3
            velocity_X.append(ops.nodeVel(3, 1))       # LATERAL VELOCITY IN X FOR NODE 3
            acceleration_X.append(ops.nodeAccel(3, 1)) # LATERAL ACCELERATION IN X FOR NODE 3
            velocity_Y.append(ops.nodeVel(3, 2))       # LATERAL VELOCITY IN Y FOR NODE 3
            acceleration_Y.append(ops.nodeAccel(3, 2)) # LATERAL ACCELERATION IN Y FOR NODE 3
            FORCE_S.append(S)
            FORCE_A.append(A)
            MOMENT.append(M)
            DISP_X.append(disp_X)
            DISP_Y.append(disp_Y)
            ROT.append(rot)
            KS.append(np.abs(S/disp_X)) # LATERAL STIFFNESS IN X
            KA.append(np.abs(A/disp_Y)) # LATERAL STIFFNESS IN Y
            KI.append(np.abs(M/rot))    # ROTATIONAL STIFFNESS IN Z
            OMEGA.append(np.sqrt(KS[-1] / (2*MASS+ELE_MASS)))
            PERIOD.append((np.pi * 2) / OMEGA[-1])
            #print(current_time, disp_X, S)
        
        # Calculating Damping Ratio Using Logarithmic Decrement Analysis 
        displacement = np.array(DISP_X)
        peaks = np.array([displacement[i] for i in range(1, len(displacement)-1) if displacement[i] > displacement[i-1] and displacement[i] > displacement[i+1]])
        # Natural logarithm
        delta = np.log(peaks[:-1] / peaks[1:])    
        
        # Mean Period of Structure based of each step lateral structural stiffness
        #Period_Mean = np.mean(PERIOD)
        #Period_Median = np.median(PERIOD) 
        
        #ops.wipe()  
        return FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y, ROT, KA, KS, KI, time, velocity_X, velocity_Y, acceleration_X, acceleration_Y, PERIOD_01, PERIOD_02, delta


#%%------------------------------------------------------------------------------
# FIND BEST COLUMN CONCRETE SECTION DEPTH WITH STRUCTURAL PERIOD:
    
X = Hc             # [mm] Intial Guess Column Section Depth 
ESP = 1e-5         # Finite difference derivative Convergence Tolerance
TOLERANCE = 1e-6   # Convergence Tolerance
RESIDUAL = 100     # Convergence Residual 
IT = 0             # Intial Iteration
ITMAX = 100000     # Max. Iteration
TARGET_P = 7.5     # [s] Target Structural Period

STEEL_KIND, ANA_KIND = 2, 'DYNAMIC'
# Analysis Durations:
starttime = TI.process_time()

### FIND THE OPTIMUM VALUE 
while (RESIDUAL > TOLERANCE):
    # X -------------------------------------------------------
    DATA  = PD_ANALYSIS(X, STEEL_KIND, ANA_KIND)
    FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y, ROT, KA, KS, KI, time, velocity_X, velocity_Y, acceleration_X, acceleration_Y, PERIOD_01, PERIOD_02, delta = DATA
    print(f'PERIOD_01 : {PERIOD_01:.8f} - PERIOD_02 : {PERIOD_02:.8f}')
    P_ANA = max(PERIOD_01, PERIOD_02)
    F = P_ANA - TARGET_P
    print('F: ', F)
    # Xmin -------------------------------------------------------
    XMIN = X - ESP  
    DATA  = PD_ANALYSIS(XMIN, STEEL_KIND, ANA_KIND)
    FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y, ROT, KA, KS, KI, time, velocity_X, velocity_Y, acceleration_X, acceleration_Y, PERIOD_01, PERIOD_02, delta = DATA
    P_ANA = max(PERIOD_01, PERIOD_02)
    Fmin = P_ANA - TARGET_P
    print('Fmin: ', Fmin)
    # Xmax -------------------------------------------------------
    XMAX = X + ESP 
    DATA  = PD_ANALYSIS(XMAX, STEEL_KIND, ANA_KIND)
    FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y, ROT, KA, KS, KI, time, velocity_X, velocity_Y, acceleration_X, acceleration_Y, PERIOD_01, PERIOD_02, delta = DATA
    P_ANA = max(PERIOD_01, PERIOD_02)
    Fmax = P_ANA - TARGET_P
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
        print(f'\t\t Optimum Column Section Depth :            {X:.6f}')
        print(f'\t\t Iteration Counts:                         {IT}')
        print(f'\t\t Convergence Residual:                     {RESIDUAL:.10e}')

    

totaltime = TI.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n')
#%%------------------------------------------------------------------------------
# %% Plot 2D Frame Shapes for Nonlinear Static Analysis
S04.PLOT_2D_FRAME(deformed_scale=10)  # Adjust scale factor as needed
#%%------------------------------------------------------------------------------

plt.figure(1, figsize=(12, 8))
plt.plot(MOMENT, FORCE_A, color='purple')
#plt.scatter(MOMENT, FORCE_A, color='purple', linewidth=2)
plt.title('P-M Interaction')
plt.ylabel('Axial Force [N]')
plt.xlabel('Bending Moment [N.mm]')
plt.grid()
plt.show()

plt.figure(2, figsize=(12, 8))
plt.plot(DISP_X, FORCE_S, color='lime', linewidth=2)
#plt.scatter(DISP_X, FORCE_S, color='green', linewidth=2)
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
plt.ylabel('Moment [N.mm]')
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
#plt.plot(KI, KA, color='grey', linewidth=2)
plt.scatter(KI, KA, color='grey', linewidth=2)
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
#plt.scatter(timed, MOMENTd, color='green', linewidth=2)
plt.title('Moment During the Analysis')
plt.ylabel('Moment [N.mm]')
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
#plt.scatter(time, DISP_Y, color='brown', linewidth=2)
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

#%%------------------------------------------------------------------------------
#%% Plot the Histogram of Data
S05.HISROGRAM_BOXPLOT(DISP_X, HISTO_COLOR='purple', LABEL='Displacement in X [mm]')
S05.HISROGRAM_BOXPLOT(DISP_Y, HISTO_COLOR='green', LABEL='Displacement in Y [mm]')
S05.HISROGRAM_BOXPLOT(ROT, HISTO_COLOR='brown', LABEL='Rotation in Z [rad]')

S05.HISROGRAM_BOXPLOT(FORCE_S, HISTO_COLOR='lime', LABEL='Shear-reaction in X [N]')
S05.HISROGRAM_BOXPLOT(FORCE_A, HISTO_COLOR='cyan', LABEL='Axial-reaction in X [N]')
S05.HISROGRAM_BOXPLOT(MOMENT, HISTO_COLOR='pink', LABEL='Moment-reaction in Z [N.mm]')

S05.HISROGRAM_BOXPLOT(KS, HISTO_COLOR='blue', LABEL='Lateral Stiffness in X Dir. [N/mm]')
S05.HISROGRAM_BOXPLOT(KA, HISTO_COLOR='orange', LABEL='Lateral Stiffness in Y Dir. [N/mm]')
S05.HISROGRAM_BOXPLOT(KI, HISTO_COLOR='grey', LABEL='Rotational Stiffness [N.mm/Rad]')
#%%------------------------------------------------------------------------------  
# EXACT SOLUTION:
from scipy.optimize import fsolve

# Define the equation for natural logarithm of this ratio, called the logarithmic decrement, we denote by δ
def EQUATION(x, delta):
    if np.any(x == 0):  # Avoid division by zero
        return np.inf  
        
    # Calculate the value of the equation
    A = x**2 - 1 + ((2 * np.pi * x) / np.mean(delta)) ** 2
    #print(f"x: {x}, A: {A}")  # Debugging output
    # Return the difference (for root finding)
    return A
      
# Initial guess for root(s)
x0 = 1  # Intial Guess for Damping Ratio
# Solve for x
solution = fsolve(EQUATION, x0, args=(delta))
print(f"Exact Damping Ratio: {solution[0]:.8e}")
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
plt.title(f'Time vs Displacement - MAX. ABS: {DISP_ZX[-1]} | ξ (Calculated): {100*solution[0]:.5e} %')
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

plt.figure(6, figsize=(8, 6))
plt.plot(DISP_X, FORCE_S, color='black', linewidth=2)
plt.xlabel('Displacement in X [mm]')
plt.ylabel('Shear Base-reaction [N]')
plt.title(f'Displacement vs Shear Base-reaction - MAX. ABS: {np.max(np.abs(FORCE_S)):.4f}')
plt.grid()
plt.show()

plt.figure(7, figsize=(8, 6))
plt.plot(DISP_Y, FORCE_A, color='brown', linewidth=2)
plt.xlabel('Displacement in Y [mm]')
plt.ylabel('Axial Base-reaction [N]')
plt.title(f'Displacement vs Axial Base-reaction - MAX. ABS: {np.max(np.abs(FORCE_A)):.4f}')
plt.grid()
plt.show()
#%%------------------------------------------------------------------------------   
# EXPORT DATA TO EXCEL
import pandas as pd
DATA_TOTAL = { 
    'DISP_X': DISP_X,
    'DISP_Y': DISP_Y,
    'ROTATION': ROT,
    'VELO': velocity_X,
    'ACCEL': acceleration_X,
    'AXIAL_FORCE': FORCE_A,
    'SHEAR_FORCE': FORCE_S,
    'MOMENT': MOMENT,
    'ROTATIONAL_ST': KI,
    'LATERAL_ST_Y': KA,
    'LATERAL_ST_X': KS,
}
# Convert to DataFrame
results_df = pd.DataFrame(DATA_TOTAL)
# Export the DataFrame to an Excel file
results_df.to_excel('DYNAMIC_Cdepth_PERIOD_OPTIMIZATION_RESULTS.xlsx', index=False) 
#%%------------------------------------------------------------------------------
# Print out the state of nodes 3 and 4
ops.printModel("node",3, 4)
# Print out the state of element 1 , 2 and 3
ops.printModel("ele", 1, 2 , 3)
# Print the Model
#printModel()
ops.printModel("-JSON", "-file", "DYNAMIC_Cdepth_PERIOD_OPTIMIZATION.json")
#%%-------------------------------------------------------------------------------
    
