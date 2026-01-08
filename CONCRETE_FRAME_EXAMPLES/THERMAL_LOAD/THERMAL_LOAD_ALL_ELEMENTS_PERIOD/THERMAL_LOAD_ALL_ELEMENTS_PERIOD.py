################################################################################################################
#                                                  IN THE NAME OF ALLAH                                        #
#                                  THERMAL ANALYSIS OF CONCRETE FRAME USING OPENSEES                           #
#               EVALUATION OF THE STRUCTURAL PERIOD DURING THERMAL LOAD APPLIED TO ALL ELEMENTS                #
#--------------------------------------------------------------------------------------------------------------#
#                             THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                       #
#                                        EMAIL: salar.d.ghashghaei@gmail.com                                   #
################################################################################################################

"""
Models and Analyzes a 2D Concrete Frame subjected to Thermal and Distributed Loads using OpenSees. 
Key points:  
[1] Model Definition: The 2D frame has specified node coordinates for stories and bays, with fixed supports at the base.
 Material properties for concrete (with thermal effects) and concrete rectangular section geometries are defined using fiber elements.  
[2] Element and Load Setup: Beam-column elements with corotational geometric transformation and Lobatto beam integration
 are created. Distributed loads are applied to beams, and a thermal gradient is applied to the first story beams.  
[3] Analysis Setup: The analysis uses static load control with thermal increments, and the Newton-Raphson algorithm ensures convergence.
 Convergence tolerances and maximum iterations are defined.  
[4] Output and Post-processing: Displacements, reactions, and deformations are recorded during the analysis.
 Data is extracted from output files for plotting base reactions (axial, shear, moment) and node displacements
 against temperature or applied load.  
[5] Visualization: The frame's undeformed and deformed shapes are plotted, and results like temperature-displacement
 relationships and base reactions are visualized.
"""
import openseespy.opensees as ops
import matplotlib.pyplot as plt
import numpy as np
import time as TI
import ANALYSIS_FUNCTION as S02
import CONCRETE_THERMAL_SECTION_FUN as S03
#import CONCRETE_FIBERTHERMAL_SECTION as S03
import PLOT_2D as S04
import EIGENVALUE_ANALYSIS_FUN as S05
import GRAVITY_ANALYSIS_FUN as S06
#%%------------------------------------------------------------------------------

# Define materials for nonlinear columns and beam
# Define parameters (units: mm, N)
# ------------------------------------------
# CONCRETE                  tag   f'c        ec0   f'cu        ecu
# Cover concrete (unconfined)
fc = 18                  # [N/mm²] Concrete Compressive Strength
Kc = 1.35

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

# Define Thermal and Distributed load
Max_Thermal = 900.0      # [°C] Temperature
dl = -0.3                # [N/mm] Distributed Load

LENGTH_COL = 3000        # [mm] Column Length 
LENGTH_BM = 7000         # [mm] Beam Length 
# Define Analysis Properties
MAX_ITERATIONS = 10000   # Convergence iteration for test
MAX_TOLERANCE = 1.0e-10  # Convergence tolerance for test

Nstep = 500              # Number of incremental steps
Incr_Temp = 1/Nstep      # Incremental temperature step

MASS = 10000.0          # [kg] Mass on the Top of each column
#%%------------------------------------------------------------------------------
def THERMAL_ANALYSIS():
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
    print(' Structure Nodes Done.')    
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 1, 1, 1)

    secTagC = 10
    secTagB = 20

    # COLUMN THERMAL SECTION
    S03.CONFINED_CONCRETE_SECTION_THERMAL(secTagC, fc, Kc, Hc, Bc, coverC, AsC, COL=True)
    # BEAM THERMAL SECTION
    S03.CONFINED_CONCRETE_SECTION_THERMAL(secTagB, fc, Kc-0.15, Hb, Bb, coverB, AsB, COL=False)
    
    print(' Thermal Section Done.')
    
    # Define geometric transformation (corotational for large displacements)
    transfTag = 1
    #ops.geomTransf('Linear', transfTag)
    #ops.geomTransf('PDelta', transfTag)
    ops.geomTransf('Corotational', transfTag)
    # Define beam integration (Lobatto integration)
    numIntegrationPoints = 5
    biTagB = 1
    ops.beamIntegration('Lobatto', biTagB, secTagB, numIntegrationPoints)
 
    # Define column integration (Lobatto integration)
    numIntegrationPoints = 5
    biTagC = 2
    ops.beamIntegration('Lobatto', biTagC, secTagC, numIntegrationPoints)
    #--------------------------------------------------------------------
    ops.element('dispBeamColumnThermal', 1, 1, 8, transfTag, biTagC)  # COLUMN 01
    ops.element('dispBeamColumnThermal', 2, 8, 9, transfTag, biTagC)  # COLUMN 02
    ops.element('dispBeamColumnThermal', 3, 9, 10, transfTag, biTagC) # COLUMN 03
    ops.element('dispBeamColumnThermal', 4, 10, 3, transfTag, biTagC) # COLUMN 04
    
    ops.element('dispBeamColumnThermal', 5, 2, 11, transfTag, biTagC)  # COLUMN 05
    ops.element('dispBeamColumnThermal', 6, 11, 12, transfTag, biTagC) # COLUMN 06
    ops.element('dispBeamColumnThermal', 7, 12, 13, transfTag, biTagC) # COLUMN 07
    ops.element('dispBeamColumnThermal', 8, 13, 4, transfTag, biTagC)  # COLUMN 08
    
    ops.element('dispBeamColumnThermal', 9, 3, 5, transfTag, biTagB)  # BEAM 01
    ops.element('dispBeamColumnThermal', 10, 5, 6, transfTag, biTagB) # BEAM 02
    ops.element('dispBeamColumnThermal', 11, 6, 7, transfTag, biTagB) # BEAM 03
    ops.element('dispBeamColumnThermal', 12, 7, 4, transfTag, biTagB) # BEAM 04
    print(' Thermal Elements Done.')
    
    # Data storage
    FORCE_S, FORCE_A, MOMENT = [], [], []
    DISP_X, DISP_Y, ROT = [], [], []
    KA, KS, KI, THERMAL_STEP = [], [], [], []
    PERIOD_MIN, PERIOD_MAX = [], []
    # Define time series and load pattern
    Time_Tag = 1
    ops.timeSeries('Linear', Time_Tag)
    patternTag = 1
    ops.pattern('Plain', patternTag , Time_Tag)
    # Apply distributed load to all beams in the structure
    # eleLoad -ele $eleTag1 <$eleTag2 ....> -type -beamUniform $Wz <$Wx>
    ops.eleLoad('-ele', 9, '-type', '-beamUniform', dl) 
    ops.eleLoad('-ele', 10, '-type', '-beamUniform', dl)
    ops.eleLoad('-ele', 11, '-type', '-beamUniform', dl)
    ops.eleLoad('-ele', 12, '-type', '-beamUniform', dl)
    NstepGravity = 10
    S06.GRAVITY_ANALYSIS_FUN(NstepGravity, MAX_TOLERANCE, MAX_ITERATIONS)
    print(' Distributed load Done.')
    
    # Define masses
    ops.mass(3, MASS, MASS, 0.0)
    ops.mass(4, MASS, MASS, 0.0)
    
    # A linear thermal gradient is applied to the elements in the first story.
    # The temperature varies from -Db to Db across the section height.
    # INFO LINK: https://openseesforfire.github.io/Subpages/ThermalActionCmds.html
    # Apply thermal load only to the identified beam elements
    Db = 0.5 * Hb
    Dc = 0.5 * Hc
    # eleLoad -ele $eleTag -type -beamThermal $T1 $y1 $T2 $Y2
    ops.eleLoad('-ele', 1, '-type', '-beamThermal', Max_Thermal, -Dc, Max_Thermal, Dc)  
    ops.eleLoad('-ele', 2, '-type', '-beamThermal', Max_Thermal, -Dc, Max_Thermal, Dc)  
    ops.eleLoad('-ele', 3, '-type', '-beamThermal', Max_Thermal, -Dc, Max_Thermal, Dc)  
    ops.eleLoad('-ele', 4, '-type', '-beamThermal', Max_Thermal, -Dc, Max_Thermal, Dc)  
    ops.eleLoad('-ele', 5, '-type', '-beamThermal', Max_Thermal, -Dc, Max_Thermal, Dc)  
    ops.eleLoad('-ele', 6, '-type', '-beamThermal', Max_Thermal, -Dc, Max_Thermal, Dc)  
    ops.eleLoad('-ele', 7, '-type', '-beamThermal', Max_Thermal, -Dc, Max_Thermal, Dc)  
    ops.eleLoad('-ele', 8, '-type', '-beamThermal', Max_Thermal, -Dc, Max_Thermal, Dc)  
    ops.eleLoad('-ele', 9, '-type', '-beamThermal', Max_Thermal, -Db, Max_Thermal, Db)  
    ops.eleLoad('-ele', 10, '-type', '-beamThermal', Max_Thermal, -Db, Max_Thermal, Db)  
    ops.eleLoad('-ele', 11, '-type', '-beamThermal', Max_Thermal, -Db, Max_Thermal, Db)  
    ops.eleLoad('-ele', 12, '-type', '-beamThermal', Max_Thermal, -Db, Max_Thermal, Db)  
    print(' Thermal Load Done.')

    # Run Thermal Analysis
    ops.system('BandGeneral')
    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/normDispIncr.html
    ops.algorithm('Newton')
    ops.integrator('LoadControl', Incr_Temp)
    ops.analysis('Static')
    print(' Thermal Analysis Properties Done.')
    
    ops.recorder('Node', '-file', "BTH_PUSH_01.txt",'-time', '-node', 1, '-dof', 1,2,3, 'reaction') # Base Reaction Time History Node 1
    ops.recorder('Node', '-file', "BTH_PUSH_02.txt",'-time', '-node', 2, '-dof', 1,2,3, 'reaction') # Base Reaction Time History Node 2
    ops.recorder('Node', '-file', "DTH_PUSH.txt",'-time', '-node', 3, '-dof', 1,2,3, 'disp')        # Displacement Time History Node 3
    #ops.recorder('Element', '-file', 'STRESS_STRAIN_BOT_CONCRETE.txt', '-time', '-ele', 9, 'section', secTagB, 'fiber', -Db, 0, 'stressStrain')
    #ops.recorder('Element', '-file', 'STRESS_STRAIN_TOP_CONCRETE.txt', '-time', '-ele', 9, 'section', secTagB, 'fiber', +Db, 0, 'stressStrain')
    #ops.recorder('Element', '-file', 'STRESS_STRAIN_BOT_REBAR.txt', '-time', '-ele', 9, 'section', secTagB, 'fiber', -Db+coverB, -0.5*Bb+coverB, 'stressStrain')
    #ops.recorder('Element', '-file', 'STRESS_STRAIN_TOP_REBAR.txt', '-time', '-ele', 9, 'section', secTagB, 'fiber', +Db-coverB, -0.5*Bb+coverB, 'stressStrain')
    
    for step in range(Nstep):
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
        THERMAL_STEP.append(ops.getLoadFactor(patternTag) * Max_Thermal) # THERMAL LOAD 
        # IN EACH STEP, STRUCTURE PERIOD IS CALCULATED
        PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(4, PLOT=True)
        PERIOD_MIN.append(PERIODmin)
        PERIOD_MAX.append(PERIODmax)
        print(step+1, disp_X, disp_Y)
        
    #ops.wipe()    

    return FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y, ROT, KA, KS, KI, THERMAL_STEP, np.array(PERIOD_MIN), np.array(PERIOD_MAX)

#%%------------------------------------------------------------------------------
# Analysis Durations:
starttime = TI.process_time()

# WITHOUT HARDENING AND ULTIMATE STRAIN
FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y, ROT, KA, KS, KI, THERMAL_STEP, PERIOD_MIN, PERIOD_MAX = THERMAL_ANALYSIS()

totaltime = TI.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n') 
#%%------------------------------------------------------------------------------
"""
def OUTPUT_SECOND_COLUMN(X, COLUMN):
    import numpy as np
    # Time History
    filename = f"{X}.txt"
    data_collected = np.loadtxt(filename)
    X = data_collected[:, COLUMN]   
    return X 

disp_X = OUTPUT_SECOND_COLUMN('DTH_PUSH', 1) # Reading Disp from Text file - X Direaction
disp_Y = OUTPUT_SECOND_COLUMN('DTH_PUSH', 2) # Reading Disp from Text file - Y Direaction
disp_Z = OUTPUT_SECOND_COLUMN('DTH_PUSH', 3) # Reading Disp from Text file - Z Direaction (Rotation)
# AXIAL
base01_Y = OUTPUT_SECOND_COLUMN('BTH_PUSH_01', 1) # Reading base reaction from Text file - NODE 1
base02_Y = OUTPUT_SECOND_COLUMN('BTH_PUSH_02', 1) # Reading base reaction from Text file - NODE 2
BASES_AXIAL = base01_Y + base02_Y
# SHEAR
base01_X = OUTPUT_SECOND_COLUMN('BTH_PUSH_01', 2) # Reading base reaction from Text file - NODE 1
base02_X = OUTPUT_SECOND_COLUMN('BTH_PUSH_02', 2) # Reading base reaction from Text file - NODE 2
BASES_SHEAR = base01_X + base02_X
# MOMENT
base01_Z = OUTPUT_SECOND_COLUMN('BTH_PUSH_01', 3) # Reading base reaction from Text file - NODE 1
base02_Z = OUTPUT_SECOND_COLUMN('BTH_PUSH_02', 3) # Reading base reaction from Text file - NODE 2
BASES_MOMENT = base01_Z + base02_Z
# STRESS AND STRAIN OF ELEMENT 6
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
#%%------------------------------------------------------------------------------
# Plot results
plt.figure(1, figsize=(8, 6))
plt.plot(disp_Y, BASES_AXIAL, color='blue', linewidth=2)
plt.xlabel('Node 6 Displacement Y (mm)')
plt.ylabel('Base Axial Reaction (kN)')
plt.title('Base Axial Reaction vs Displacement Y')
plt.grid()
plt.show()

plt.figure(2, figsize=(8, 6))
plt.plot(disp_X, BASES_SHEAR, color='red', linewidth=2)
plt.xlabel('Node 6 Displacement X (mm)')
plt.ylabel('Base Shear Reaction (kN)')
plt.title('Base Shear Reaction vs Displacement X')
plt.grid()
plt.show()

plt.figure(3, figsize=(8, 6))
plt.plot(disp_Z, BASES_MOMENT, color='green', linewidth=2)
plt.xlabel('Node 6 Rotation (rad)')
plt.ylabel('Base Moment Reaction (kN.mm)')
plt.title('Base Moment Reaction vs Rotation')
plt.grid()
plt.show()

plt.figure(4, figsize=(8, 6))
plt.plot(THERMAL_STEP, DISP_X, color='purple', linewidth=2)
plt.xlabel('Temperature (°C)')
plt.ylabel(f'Node {3} Displacement X (mm)')
plt.title('Temperature vs Displacement X')
plt.grid()
plt.show()

plt.figure(5, figsize=(8, 6))
plt.plot(THERMAL_STEP, DISP_Y, color='black', linewidth=2)
plt.xlabel('Temperature (°C)')
plt.ylabel(f'Node {3} Displacement Y (mm)')
plt.title('Temperature vs Displacement Y')
plt.grid()
plt.show()

plt.figure(6, figsize=(8, 6))
plt.plot(strain_B_C, stress_B_C, color='blue', label='Bottom Fiber', linewidth=2)
plt.plot(strain_T_C, stress_T_C, color='red', label='Top Fiber', linewidth=2)
plt.xlabel('Strain (mm/mm)')
plt.ylabel('Stress (kN/mm^2)')
plt.title(f'Stress-Strain Relation of Element {6} Concrete Top & Bottom Fibers')
plt.grid()
plt.legend()
plt.show()

plt.figure(7, figsize=(8, 6))
plt.plot(strain_B_R, stress_B_R, color='blue', label='Bottom Fiber', linewidth=2)
plt.plot(strain_T_R, stress_T_R, color='red', label='Top Fiber', linewidth=2)
plt.xlabel('Strain (mm/mm)')
plt.ylabel('Stress (kN/mm^2)')
plt.title(f'Stress-Strain Relation of Element {6} Steel Rebar Top & Bottom Fibers')
plt.grid()
plt.legend()
plt.show()
"""
#%%------------------------------------------------------------------------------
# PLOT STRUCTURAL PERIOD DURING THE THERMAL ANALYSIS
plt.figure(0, figsize=(12, 8))
plt.plot(DISP_X, PERIOD_MIN, color='black', linewidth=4)
plt.plot(DISP_X, PERIOD_MAX, color='red', linewidth=4)
plt.title('Period of Structure vs Lateral Displacement During Thermal Analysis')
plt.ylabel('Structural Period [s]')
plt.xlabel('Lateral Displacement [mm] (Node 3)')
#plt.semilogy()
plt.grid()
plt.legend([f'PERIOD - MIN VALUES: Min: {np.min(PERIOD_MIN):.3f} (s) - Mean: {np.mean(PERIOD_MIN):.3f} (s) - Max: {np.max(PERIOD_MIN):.3f} (s)', 
            f'PERIOD - MAX VALUES:  Min: {np.min(PERIOD_MAX):.3f} (s) - Mean: {np.mean(PERIOD_MAX):.3f} (s) - Max: {np.max(PERIOD_MAX):.3f} (s)',
            ])
plt.show()

plt.figure(11, figsize=(12, 8))
plt.plot(THERMAL_STEP, PERIOD_MIN, color='black', linewidth=4)
plt.plot(THERMAL_STEP, PERIOD_MAX, color='red', linewidth=4)
plt.title('Period of Structure vs Thermal Load Temperature During Thermal Analysis')
plt.ylabel('Structural Period [s]')
plt.xlabel('Thermal Load Temperature [°C]')
#plt.semilogy()
plt.grid()
plt.legend([f'PERIOD - MIN VALUES: Min: {np.min(PERIOD_MIN):.3f} (s) - Mean: {np.mean(PERIOD_MIN):.3f} (s) - Max: {np.max(PERIOD_MIN):.3f} (s)', 
            f'PERIOD - MAX VALUES:  Min: {np.min(PERIOD_MAX):.3f} (s) - Mean: {np.mean(PERIOD_MAX):.3f} (s) - Max: {np.max(PERIOD_MAX):.3f} (s)',
            ])
plt.show()
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
plt.plot(THERMAL_STEP, FORCE_A, color='brown', linewidth=2)
#plt.scatter(THERMAL_STEP, FORCE_A, color='brown', linewidth=2)
plt.title('Axial Force During the Analysis')
plt.ylabel('Axial Force [N]')
plt.xlabel('Thermal Load Temperature [°C]')
plt.grid()
plt.show()

plt.figure(8, figsize=(12, 8))
plt.plot(THERMAL_STEP, FORCE_S, color='purple', linewidth=2)
#plt.scatter(THERMAL_STEP, FORCE_S, color='purple', linewidth=2)
plt.title('Shear Force During the Analysis')
plt.ylabel('Shear Force [N]')
plt.xlabel('Thermal Load Temperature [°C]')
plt.grid()
plt.show()

plt.figure(9, figsize=(12, 8))
plt.plot(THERMAL_STEP, MOMENT, color='green', linewidth=2)
#plt.scatter(THERMAL_STEP, MOMENT, color='green', linewidth=2)
plt.title('Moment During the Analysis')
plt.ylabel('Moment [kN.mm]')
plt.xlabel('Distributed Load THERMAL_STEPs [N/mm]')
plt.grid()
plt.show()

plt.figure(10, figsize=(12, 8))
plt.plot(THERMAL_STEP, DISP_X, color='brown', linewidth=2)
#plt.scatter(THERMAL_STEP, DISP_X, color='brown', linewidth=2)
plt.title('Displacement During the Analysis')
plt.ylabel('Displacement - X [mm]')
plt.xlabel('Thermal Load Temperature [°C]')
plt.grid()
plt.show()

plt.figure(11, figsize=(12, 8))
plt.plot(THERMAL_STEP, DISP_Y, color='blue', linewidth=2)
#plt.scatter(THERMAL_STEP, DISP_X, color='brown', linewidth=2)
plt.title('Displacement During the Analysis')
plt.ylabel('Displacement - Y [mm]')
plt.xlabel('Thermal Load Temperature [°C]')
plt.grid()
plt.show()

plt.figure(12, figsize=(12, 8))
plt.plot(THERMAL_STEP, ROT, color='black', linewidth=2)
#plt.scatter(THERMAL_STEP, ROT, color='black', linewidth=2)
plt.title('Rotation During the Analysis')
plt.ylabel('Rotation [rad]')
plt.xlabel('Thermal Load Temperature [°C]')
plt.grid()
plt.show()
#%%------------------------------------------------------------------------------  
# %% Plot 2D Frame Shapes
S04.PLOT_2D_FRAME(deformed_scale=10)  # Adjust scale factor as needed
#%%------------------------------------------------------------------------------   
# EXPORT DATA TO EXCEL
import pandas as pd
DATA_TOTAL = { 
    'THERMAL_LOAD_STEPS': THERMAL_STEP,# Thermal Load Temperature [°C]
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
results_df.to_excel('THERMAL_LOAD_ALL_ELEMENTS_PERIOD_RESULTS.xlsx', index=False) 
#%%------------------------------------------------------------------------------

    
