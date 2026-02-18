######################################################################################################################
#                          >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                          #
#                            MODELING OF 2D DOUBLE PENDULUM MDOF STRUCTURE USING OPENSEES                            #
#--------------------------------------------------------------------------------------------------------------------#
#                   EVALUATION OF DAMPING FORCE (fD), SPRING FORCE (fS) AND INERTIA FORCE (fI)                       #
#--------------------------------------------------------------------------------------------------------------------#
#         THIS PROGRAM WRITTEN BY MICHAEL H. SCOTT AND MODIFIED BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                #
#                                       EMAIL: salar.d.ghashghaei@gmail.com                                          #
######################################################################################################################
"""
This code performs nonlinear time history analysis of a 2D double pendulum
 (truss elements) under harmonic base excitation, comparing elastic vs inelastic
 material behavior. It computes dynamic responses (displacement, velocity, acceleration, reactions)
 at each node, extracts time-varying period and stiffness degradation, and calculates 
 damping ratios from response histories. The analysis uses corotational truss elements
 with either Elastic or Hysteretic (steel) material models, Newmark integration, and
 accounts for gravity loads and initial geometric imperfections.
 Key outputs include force-displacement hysteresis, damping force-velocity relationships,
 and period elongation due to inelastic action - essential for understanding seismic
 performance and energy dissipation in nonlinear structures.
 Very helpful Website for better learning:
 https://portwooddigital.com/2025/09/08/double-inverted-pendulum/    
"""
#%%-----------------------------------------------------------------------------
# YOUUBE: Simple Pendulum
'https://www.youtube.com/watch?v=fnvGVsxPuLs'
# YOUTUBE: Everything You Need To Know About Pendulums: Physics Help Room
'https://www.youtube.com/watch?v=0q0L7Fj4dk8'
# YOUTUBE: How a Giant Pendulum Made Taipei101 Possible
'https://www.youtube.com/watch?v=mGe9zjwK2gQ'
# BOOK: Differential Equations for Engineers - Wei-Chau Xie - CAMBRIDGE
'https://www.cambridge.org/core/books/differential-equations-for-engineers/1B8F1A62BF6F98EB972CFE9114FA8B84'
#%%-----------------------------------------------------------------------------
import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
import time as TI
import ANALYSIS_FUNCTION as S01
import PERIOD_FUN as S02
import DAMPING_RATIO_FUN as S05
import EIGENVALUE_ANALYSIS_FUN as S03
import GENERATE_ARTIFICIAL_ACCEL_FUN as S04
import SALAR_MATH as S06
import RAYLEIGH_DAMPING_FUN as S07
import PLOT_2D_TRUSS as S08
#%%----------------------------------------------------------------------------- 
# Define parameters (units: mm, N, kg)
 
L1 = 780                # [mm] Element 01
L2 = 750                # [mm] Element 02
M1 = 2.5                # [kg] Mass of Node 01
M2 = 5.0                # [kg] Mass of Node 02
g = 9810.0              # [mm/s²] standard acceleration of gravity or standard acceleration    
W1 = M1 * g             # [N] Weight of Node 01
W2 = M2 * g             # [N] Weight of Node 02
 
# Elastic Element
R = 5.0                 # [mm] Element Section Radious
AREA = np.pi*R**2       # [mm^2] Element Section Area
Es = 200000.0           # [N/mm^2] Element Young's Modulus
 
f = 12.0                # [Hz] natural frequency
AMPL = 1                # Amplification Factor
#%%----------------------------------------------------------------------------- 
duration = 200.0                  # [s] Analysis duration
dt = 0.001                        # [s] Time step
DR = 0.03                         # Damping Ratio
#%%-----------------------------------------------------------------------------
# DEFINE ANALYSIS PROPERTIES
MAX_ITERATIONS = 20000    # Convergence iteration for test
MAX_TOLERANCE = 1.0e-6   # Convergence tolerance for test
#SPRING_KIND: 1 -> 'ELASTIC'
#SPRING_KIND: 2 -> 'INELASTIC'
#%%-----------------------------------------------------------------------------
def PENDULUM_MDOF(MAT_TYPE, duration, dt, DR):
    # Create model
    ops.wipe()
    
    # Create a 2D model with 2 DOF per node
    ops.model('basic','-ndm',2,'-ndf',2)
    
    # Add nodes 
    ops.node(1, 0.0, 0.0)
    ops.node(2, 0.0, -L1);
    ops.node(3, 0.0, -L1-L2)
    
    # Fix node 1
    ops.fix(1, 1, 1)
    
    # Define mass to node 2 and 3
    ops.mass(2, M1, M1)
    ops.mass(3, M2, M2)
    
    # Material (elastic steel)
    MAT_TAG = 1
    if MAT_TYPE == 'ELASTIC':
        Es = 200000.0           # [N/mm^2] Element Young's Modulus
        #ops.uniaxialMaterial('Elastic', MAT_TAG, Es)              # TESNSION AND COMPRESSION IS SAME VALUES
        ops.uniaxialMaterial('Elastic', MAT_TAG, Es ,0.0, 0.5*Es)  # TESNSION AND COMPRESSION IS NOT SAME VALUES
        # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/ElasticUni.html
    # Material (inelastic steel)    
    if MAT_TYPE == 'INELASTIC':
        Fy = 240.0			    # [N/mm²] Steel yield stress
        Es = 200000.0		    # [N/mm²] Modulus of steel
        ey = Fy/Es			    # [mm/mm] Steel yield strain
        Fu = 1.1818*Fy          # [N/mm²] Steel Ultimate Strength
        esu = 0.12              # [mm/mm] Steel Ultimate Strain
        Esh = (Fu - Fy)/(esu - ey)
        Bs = Esh / Es           # strain-hardening ratio 
        pinchX = 0.8            # Pinching factor in X direction
        pinchY = 0.5            # Pinching factor in Y direction
        damage1 = 0.0           # Damage due to ductility
        damage2 = 0.0           # Damage due to energy
        beta = 0.1              # Stiffness degradation parameter
        ops.uniaxialMaterial('Hysteretic', MAT_TAG, Fy, ey, Fu, esu, 0.2*Fu, 1.1*esu, -Fy, -ey, -Fu, -esu, -0.2*Fu, -1.1*esu, pinchX, pinchY, damage1, damage2, beta)
        # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/Hysteretic.html
     
    # Truss Elements
    ops.element('corotTruss', 1, 1, 2, AREA, MAT_TAG, '-rho', AREA*0.00785 , '-doRayleigh', 1)
    ops.element('corotTruss', 2, 2, 3, AREA, MAT_TAG, '-rho', AREA*0.00785 , '-doRayleigh', 1)
    # INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Corotational_Truss_Element
        
    # Excitation
    if f > 0:
        T = 1/f # [s] Period
        ops.timeSeries('Sine', 1, 0, duration, T,'-factor', AMPL)
        # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/timeSeries.html
        ops.pattern('MultipleSupport', 1)
        # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/multiExcitation.html
        ops.groundMotion(1,'Plain', '-disp', 1)
        # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/groundMotion.html
        ops.imposedMotion(1, 2, 1) # node, dof, gmTag
        # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/imposedMotion.html
     
    # Self-weight
    ops.timeSeries('Constant',2)
    ops.pattern('Plain', 2, 2)
    ops.load(2, 0, -W1)
    ops.load(3, 0 ,-W2)
     
    # Initial out-of-plumb
    ops.setNodeDisp(2, 1, -L1, '-commit')
    ops.setNodeDisp(3, 1, -L1-L2, '-commit')
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/setNodeDisp.html
    
    """
    u0 = -5.35                         # [mm] Initial displacement
    v0 = 0.15                          # [mm/s] Initial velocity
    a0 = 0.065                         # [mm/s^2] Initial acceleration
    IU = True                          # Free Vibration with Initial Displacement
    IV = True                          # Free Vibration with Initial Velocity
    IA = True                          # Free Vibration with Initial Acceleration
    
    if IU == True:
        # Define initial displacment
        ops.setNodeDisp(2, 1, u0, '-commit') # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/setNodeDisp.html
    if IV == True:
        # Define initial velocity
        ops.setNodeVel(2, 1, v0, '-commit')  # INFO LINK: https://openseespydoc.readthedocs.io/en/stable/src/setNodeVel.html
    if IA == True:
        # Define initial  acceleration
        ops.setNodeAccel(2, 1, a0, '-commit') # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/setNodeAccel.html
    """
    
    ops.constraints('Transformation')
    ops.system('UmfPack')
    #ops.system('BandGeneral')#  'Umfpack'
    ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/normDispIncr.html
    #ops.integrator('CentralDifference')  # JUST FOR LINEAR ANALYSIS - INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/centralDifference.html
    alpha=0.5; beta=0.25;
    ops.integrator('Newmark', alpha, beta) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/newmark.html
    #alpha=2/3;gamma=1.5-alpha; gamma=1.5-alpha;beta=(2-alpha)**2/4;
    #ops.integrator('HHT', alpha, gamma, beta) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/hht.html
    ops.algorithm('Newton')  # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/algorithm.html
    ops.analysis('Transient') # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/analysis.html

    # Perform analysis
    time = []
    dispX_01, dispY_01, dispX_02, dispY_02 = [], [], [], []
    velX_01, velY_01, velX_02, velY_02 = [], [], [], []
    accelX_01, accelY_01, accelX_02, accelY_02 = [], [], [], []
    reactionX, reactionY = [], []
    stiffnessX, stiffnessY = [], []
    OMEGA, PERIOD = [], []
    PERIOD_MIN, PERIOD_MAX = [], []
    FDx, FSx, FIx = [], [], []
    FDy, FSy, FIy = [], [], []
        
    stable = 0
    current_time = 0.0
    
    while stable == 0 and current_time < duration:
        ops.analyze(1, dt)
        S01.ANALYSIS(stable, 1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS
        current_time = ops.getTime()
        time.append(current_time)
        ops.reactions()
        reactionX.append(ops.nodeReaction(1, 1)) # BASE REACTION IN X DIR.
        reactionY.append(ops.nodeReaction(1, 2)) # BASE REACTION IN Y DIR.
        dispX_01.append(ops.nodeDisp(2, 1))      # DISPLACEMENT IN X DIR. - NODE 2
        dispY_01.append(ops.nodeDisp(2, 2))      # DISPLACEMENT IN Y DIR. - NODE 2
        velX_01.append(ops.nodeVel(2, 1))        # VELOCITY IN X DIR. - NODE 2
        velY_01.append(ops.nodeVel(2, 2))        # VELOCITY IN Y DIR. - NODE 2
        accelX_01.append(ops.nodeAccel(2, 1))    # ACCELERATION IN X DIR. - NODE 2
        accelY_01.append(ops.nodeAccel(2, 2))    # ACCELERATION IN Y DIR. - NODE 2
        dispX_02.append(ops.nodeDisp(3, 1))      # DISPLACEMENT IN X DIR. - NODE 3
        dispY_02.append(ops.nodeDisp(3, 2))      # DISPLACEMENT IN Y DIR. - NODE 3
        velX_02.append(ops.nodeVel(3, 1))        # VELOCITY IN X DIR. - NODE 3
        velY_02.append(ops.nodeVel(3, 2))        # VELOCITY IN Y DIR. - NODE 3
        accelX_02.append(ops.nodeAccel(3, 1))    # ACCELERATION IN X DIR. - NODE 3
        accelY_02.append(ops.nodeAccel(3, 2))    # ACCELERATION IN Y DIR. - NODE 3
        stiffnessX.append(np.abs(reactionX[-1]) / np.abs(dispX_02[-1]))
        stiffnessY.append(np.abs(reactionY[-1]) / np.abs(dispY_02[-1]))
        OMEGA.append(np.sqrt(stiffnessX[-1]/(M1+M2)))
        PERIOD.append((np.pi * 2) / OMEGA[-1])
        # IN EACH STEP STRUCTURAL PERIOD WILL BE CALCULATED
        PERIODmin, PERIODmax = S03.EIGENVALUE_ANALYSIS(2, PLOT=True)
        #PERIODmin, PERIODmax = S07.RAYLEIGH_DAMPING(2, 0.5*DR, DR, 0, 1)
        PERIOD_MIN.append(PERIODmin)
        PERIOD_MAX.append(PERIODmax)
        Ci = 2 * DR * (2*np.pi/PERIODmax) * (M1+M2) # Damping COEFFICIENT - UPDATED AND CHANGES IN EACH STEP
        FDx.append(Ci * velX_02[-1])                # DAMPING FORCE IN X DIR.
        FSx.append(stiffnessX[-1] * dispX_02[-1])   # SPRING FORCE IN X DIR.
        FIx.append((M1+M2) * accelX_02[-1])         # INERTIA FORCE IN X DIR.
        FDy.append(Ci * velY_02[-1])                # DAMPING FORCE IN Y DIR.
        FSy.append(stiffnessY[-1] * dispX_02[-1])   # SPRING FORCE IN Y DIR.
        FIy.append((M1+M2) * accelY_02[-1])         # INERTIA FORCE IN Y DIR.
        print(time[-1], dispX_02[-1], velX_02[-1])
    
    # Compute modal properties
    ops.modalProperties("-print", "-file", "SALAR_ModalReport.txt", "-unorm")
    
    # Calculate Damping Ratio in X Dir.
    displacementX = np.array(dispX_02)
    damping_ratioX = S05.DAMPING_RATIO(displacementX)
    
    # Calculate Damping Ratio in Y Dir.
    displacementY = np.array(dispY_02)
    damping_ratioY = S05.DAMPING_RATIO(displacementY)
    
    
    # OUTPUTED DATA
    DATA = (time, 
     reactionX, dispX_01, velX_01, accelX_01, dispX_02, velX_02, accelX_02,
    reactionY, dispY_01, velY_01, accelY_01, dispY_02, velY_02, accelY_02,
    stiffnessX, stiffnessY, 
    PERIOD, 
    damping_ratioX, damping_ratioY, 
    FDx, FSx, FIx,
    FDy, FSy, FIy,
    )
    
    return DATA
#%% ---------------------------------------------
# Analysis Durations for Dynamic Analysis:
starttime = TI.process_time()
        
MAT_TYPE = 'ELASTIC'   # 'ELASTIC' OR 'INELASTIC'
DATA = PENDULUM_MDOF(MAT_TYPE, duration, dt, DR)
(time,
 reactionXE, dispX_01E, velX_01E, accelX_01E, dispX_02E, velX_02E, accelX_02E,
 reactionYE, dispY_01E, velY_01E, accelY_01E, dispY_02E, velY_02E, accelY_02E,
 stiffnessXE, stiffnessYE, 
 periodE,
 E_damping_ratioXE, E_damping_ratioYE,
 FDXe, FSXe, FIXe,
 FDYe, FSYe, FIYe,
 ) = DATA

S02.PERIOD_FUN(dispX_02E, dt)

S08.PLOT_2D_FRAME_TRUSS(deformed_scale=1.0)
#%%----------------------------------------------------
MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'
DATA = PENDULUM_MDOF(MAT_TYPE, duration, dt, DR)
(time,
 reactionXI, dispX_01I, velX_01I, accelX_01I, dispX_02I, velX_02I, accelX_02I,
 reactionYI, dispY_01I, velY_01I, accelY_01I, dispY_02I, velY_02I, accelY_02I,
 stiffnessXI, stiffnessYI, 
 periodI,
 E_damping_ratioXI, E_damping_ratioYI,
 FDXi, FSXi, FIXi,
 FDYi, FSYi, FIYi,
 ) = DATA

S02.PERIOD_FUN(dispX_02I, dt)
S08.PLOT_2D_FRAME_TRUSS(deformed_scale=1.0)

totaltime = TI.process_time() - starttime
print(f'\nTotal Analysis Durations (s): {totaltime:.4f} \n\n')
#%% ---------------------------------------------
###   IN X DIR. NODE 02
# Plot Results
plt.figure(2, figsize=(18, 14))

# Define colors
elastic_color = 'black'
inelastic_color = 'red'

# Reaction plot
plt.subplot(6, 1, 1)
plt.plot(time, reactionXE, color=elastic_color, linewidth=1.5)
plt.plot(time, reactionXI, color=inelastic_color, linewidth=1.5)
plt.title('Reaction Forces vs Time in X Dir. ', fontsize=12, pad=10)
plt.ylabel('Reaction (N)', fontsize=10)
plt.grid(alpha=0.3)

# Displacement plot
plt.subplot(6, 1, 2)
plt.plot(time, dispX_01E, color=elastic_color, linewidth=1.5, label=f'Elastic Damping Ratio: {E_damping_ratioXE:.3e} %')
plt.plot(time, dispX_01I, color=inelastic_color, linewidth=1.5, label=f'Inelastic Damping Ratio: {E_damping_ratioXI:.3e} %')
plt.title('Displacement vs Time in X Dir. ', fontsize=12, pad=10)
plt.ylabel('Displacement (mm)', fontsize=10)
plt.legend(loc='upper right', framealpha=1)
plt.grid(alpha=0.3)

# Velocity plot
plt.subplot(6, 1, 3)
plt.plot(time, velX_01E, color=elastic_color, linewidth=1.5)
plt.plot(time, velX_01I, color=inelastic_color, linewidth=1.5)
plt.title('Velocity vs Time in X Dir. ', fontsize=12, pad=10)
plt.ylabel('Velocity (mm/s)', fontsize=10)
plt.grid(alpha=0.3)

# Acceleration plot
plt.subplot(6, 1, 4)
plt.plot(time, accelX_01E, color=elastic_color, linewidth=1.5)
plt.plot(time, accelX_01I, color=inelastic_color, linewidth=1.5)
plt.title('Acceleration vs Time in X Dir. ', fontsize=12, pad=10)
plt.ylabel('Acceleration (mm/s²)', fontsize=10)
plt.grid(alpha=0.3)

# Stiffness plot
plt.subplot(6, 1, 5)
plt.plot(time, stiffnessXE, color=elastic_color, linewidth=1.5)
plt.plot(time, stiffnessXI, color=inelastic_color, linewidth=1.5)
plt.title('Stiffness vs Time in X Dir. ', fontsize=12, pad=10)
plt.ylabel('Stiffness (N/mm)', fontsize=10)
plt.semilogy()
plt.grid(alpha=0.3)

# Period plot
plt.subplot(6, 1, 6)
plt.plot(time, periodE, color=elastic_color, linewidth=1.5, label='Elastic Period')
plt.plot(time, periodI, color=inelastic_color, linewidth=1.5, label='Inelastic Period')
plt.title('Period vs Time in X Dir. ', fontsize=12, pad=10)
plt.ylabel('Period  (s)', fontsize=10)
plt.xlabel('Time (s)', fontsize=10)
plt.legend(loc='upper right', framealpha=1)
plt.semilogy()
plt.grid(alpha=0.3)

plt.tight_layout()
plt.show()
#%%--------------------------------------------------------
###   IN X DIR. NODE 03
# Plot Results
plt.figure(22, figsize=(18, 14))

# Reaction plot
plt.subplot(6, 1, 1)
plt.plot(time, reactionXE, color=elastic_color, linewidth=1.5)
plt.plot(time, reactionXI, color=inelastic_color, linewidth=1.5)
plt.title('Reaction Forces vs Time in X Dir. ', fontsize=12, pad=10)
plt.ylabel('Reaction (N)', fontsize=10)
plt.grid(alpha=0.3)

# Displacement plot
plt.subplot(6, 1, 2)
plt.plot(time, dispX_02E, color=elastic_color, linewidth=1.5, label=f'Elastic Damping Ratio: {E_damping_ratioXE:.3e} %')
plt.plot(time, dispX_02I, color=inelastic_color, linewidth=1.5, label=f'Inelastic Damping Ratio: {E_damping_ratioXI:.3e} %')
plt.title('Displacement vs Time in X Dir. ', fontsize=12, pad=10)
plt.ylabel('Displacement (mm)', fontsize=10)
plt.legend(loc='upper right', framealpha=1)
plt.grid(alpha=0.3)

# Velocity plot
plt.subplot(6, 1, 3)
plt.plot(time, velX_02E, color=elastic_color, linewidth=1.5)
plt.plot(time, velX_02I, color=inelastic_color, linewidth=1.5)
plt.title('Velocity vs Time in X Dir. ', fontsize=12, pad=10)
plt.ylabel('Velocity (mm/s)', fontsize=10)
plt.grid(alpha=0.3)

# Acceleration plot
plt.subplot(6, 1, 4)
plt.plot(time, accelX_02E, color=elastic_color, linewidth=1.5)
plt.plot(time, accelX_02I, color=inelastic_color, linewidth=1.5)
plt.title('Acceleration vs Time in X Dir. ', fontsize=12, pad=10)
plt.ylabel('Acceleration (mm/s²)', fontsize=10)
plt.grid(alpha=0.3)

# Stiffness plot
plt.subplot(6, 1, 5)
plt.plot(time, stiffnessXE, color=elastic_color, linewidth=1.5)
plt.plot(time, stiffnessXI, color=inelastic_color, linewidth=1.5)
plt.title('Stiffness vs Time in X Dir. ', fontsize=12, pad=10)
plt.ylabel('Stiffness (N/mm)', fontsize=10)
plt.semilogy()
plt.grid(alpha=0.3)

# Period plot
plt.subplot(6, 1, 6)
plt.plot(time, periodE, color=elastic_color, linewidth=1.5, label='Elastic Period')
plt.plot(time, periodI, color=inelastic_color, linewidth=1.5, label='Inelastic Period')
plt.title('Period vs Time in X Dir. ', fontsize=12, pad=10)
plt.ylabel('Period  (s)', fontsize=10)
plt.xlabel('Time (s)', fontsize=10)
plt.legend(loc='upper right', framealpha=1)
plt.semilogy()
plt.grid(alpha=0.3)

plt.tight_layout()
plt.show()
#%%--------------------------------------------------------
# IN Y DIR.  NODE 02
plt.figure(3, figsize=(18, 14))

# Define colors
elastic_color = 'black'
inelastic_color = 'red'

# Reaction plot
plt.subplot(6, 1, 1)
plt.plot(time, reactionYE, color=elastic_color, linewidth=1.5)
plt.plot(time, reactionYI, color=inelastic_color, linewidth=1.5)
plt.title('Reaction Forces vs Time in Y Dir. ', fontsize=12, pad=10)
plt.ylabel('Reaction (N)', fontsize=10)
plt.grid(alpha=0.3)


# Displacement plot
plt.subplot(6, 1, 2)
plt.plot(time, dispY_01E, color=elastic_color, linewidth=1.5, label=f'Elastic Damping Ratio: {E_damping_ratioYE:.3e} %')
plt.plot(time, dispY_01I, color=inelastic_color, linewidth=1.5, label=f'Inelastic Damping Ratio: {E_damping_ratioYI:.3e} %')
plt.title('Displacement vs Time in Y Dir. ', fontsize=12, pad=10)
plt.ylabel('Displacement (mm)', fontsize=10)
plt.legend(loc='upper right', framealpha=1)
plt.grid(alpha=0.3)

# Velocity plot
plt.subplot(6, 1, 3)
plt.plot(time, velY_01E, color=elastic_color, linewidth=1.5)
plt.plot(time, velY_01I, color=inelastic_color, linewidth=1.5)
plt.title('Velocity vs Time in Y Dir. ', fontsize=12, pad=10)
plt.ylabel('Velocity (mm/s)', fontsize=10)
plt.grid(alpha=0.3)

# Acceleration plot
plt.subplot(6, 1, 4)
plt.plot(time, accelY_01E, color=elastic_color, linewidth=1.5)
plt.plot(time, accelY_01I, color=inelastic_color, linewidth=1.5)
plt.title('Acceleration vs Time in Y Dir. ', fontsize=12, pad=10)
plt.ylabel('Acceleration (mm/s²)', fontsize=10)
plt.grid(alpha=0.3)

# Stiffness plot
plt.subplot(6, 1, 5)
plt.plot(time, stiffnessYE, color=elastic_color, linewidth=1.5)
plt.plot(time, stiffnessYI, color=inelastic_color, linewidth=1.5)
plt.title('Stiffness vs Time in Y Dir. ', fontsize=12, pad=10)
plt.ylabel('Stiffness (N/mm)', fontsize=10)
plt.semilogy()
plt.grid(alpha=0.3)

# Period plot
plt.subplot(6, 1, 6)
plt.plot(time, periodE, color=elastic_color, linewidth=1.5, label='Elastic Period')
plt.plot(time, periodI, color=inelastic_color, linewidth=1.5, label='Inelastic Period')
plt.title('Period vs Time in Y Dir. ', fontsize=12, pad=10)
plt.ylabel('Period  (s)', fontsize=10)
plt.xlabel('Time (s)', fontsize=10)
plt.legend(loc='upper right', framealpha=1)
plt.semilogy()
plt.grid(alpha=0.3)

plt.tight_layout()
plt.show()
#%%--------------------------------------------------------
# IN Y DIR.  NODE 03
plt.figure(33, figsize=(18, 14))

# Define colors
elastic_color = 'black'
inelastic_color = 'red'

# Reaction plot
plt.subplot(6, 1, 1)
plt.plot(time, reactionYE, color=elastic_color, linewidth=1.5)
plt.plot(time, reactionYI, color=inelastic_color, linewidth=1.5)
plt.title('Reaction Forces vs Time in Y Dir. ', fontsize=12, pad=10)
plt.ylabel('Reaction (N)', fontsize=10)
plt.grid(alpha=0.3)


# Displacement plot
plt.subplot(6, 1, 2)
plt.plot(time, dispY_02E, color=elastic_color, linewidth=1.5, label=f'Elastic Damping Ratio: {E_damping_ratioYE:.3e} %')
plt.plot(time, dispY_02I, color=inelastic_color, linewidth=1.5, label=f'Inelastic Damping Ratio: {E_damping_ratioYI:.3e} %')
plt.title('Displacement vs Time in Y Dir. ', fontsize=12, pad=10)
plt.ylabel('Displacement (mm)', fontsize=10)
plt.legend(loc='upper right', framealpha=1)
plt.grid(alpha=0.3)

# Velocity plot
plt.subplot(6, 1, 3)
plt.plot(time, velY_02E, color=elastic_color, linewidth=1.5)
plt.plot(time, velY_02I, color=inelastic_color, linewidth=1.5)
plt.title('Velocity vs Time in Y Dir. ', fontsize=12, pad=10)
plt.ylabel('Velocity (mm/s)', fontsize=10)
plt.grid(alpha=0.3)

# Acceleration plot
plt.subplot(6, 1, 4)
plt.plot(time, accelY_02E, color=elastic_color, linewidth=1.5)
plt.plot(time, accelY_02I, color=inelastic_color, linewidth=1.5)
plt.title('Acceleration vs Time in Y Dir. ', fontsize=12, pad=10)
plt.ylabel('Acceleration (mm/s²)', fontsize=10)
plt.grid(alpha=0.3)

# Stiffness plot
plt.subplot(6, 1, 5)
plt.plot(time, stiffnessYE, color=elastic_color, linewidth=1.5)
plt.plot(time, stiffnessYI, color=inelastic_color, linewidth=1.5)
plt.title('Stiffness vs Time in Y Dir. ', fontsize=12, pad=10)
plt.ylabel('Stiffness (N/mm)', fontsize=10)
plt.semilogy()
plt.grid(alpha=0.3)

# Period plot
plt.subplot(6, 1, 6)
plt.plot(time, periodE, color=elastic_color, linewidth=1.5, label='Elastic Period')
plt.plot(time, periodI, color=inelastic_color, linewidth=1.5, label='Inelastic Period')
plt.title('Period vs Time in Y Dir. ', fontsize=12, pad=10)
plt.ylabel('Period  (s)', fontsize=10)
plt.xlabel('Time (s)', fontsize=10)
plt.legend(loc='upper right', framealpha=1)
plt.semilogy()
plt.grid(alpha=0.3)

plt.tight_layout()
plt.show()
#%%--------------------------------------------------------
plt.figure(3, figsize=(8, 6))
plt.plot(dispX_01E, reactionXE, color='black', linewidth=2)
plt.plot(dispX_01I, reactionXI, color='purple', linewidth=2)
plt.plot(dispY_01E, reactionYE, color='green', linewidth=2)
plt.plot(dispY_01I, reactionYI, color='red', linewidth=2)
plt.xlabel('Displacement Node 02 [mm]', fontsize=10)
plt.ylabel('Base-reaction [N]', fontsize=10)
plt.title('Displacement Node 02 vs Base-reaction', fontsize=10)
plt.legend(['ELASTIC IN X DIR.', 'INELASTIC IN X DIR.', 'ELASTIC IN Y DIR.', 'INELASTIC IN Y DIR.'])
plt.grid()
plt.show()

plt.figure(4, figsize=(8, 6))
plt.plot(velX_01E, FDXe, color='black', linewidth=2)
plt.plot(velX_01I, FDXi, color='purple', linewidth=2)
plt.plot(velY_01E, FDYe, color='green', linewidth=2)
plt.plot(velY_01I, FDYi, color='red', linewidth=2)
plt.xlabel('Velocity Noe 02 (mm/s)', fontsize=10)
plt.ylabel('Damping Force (fD) [N]', fontsize=10)
plt.title('Damping Force (fD) vs Velocity Curve', fontsize=10)
plt.legend(['ELASTIC IN X DIR.', 'INELASTIC IN X DIR.', 'ELASTIC IN Y DIR.', 'INELASTIC IN Y DIR.'])
plt.grid()
plt.show()

plt.figure(5, figsize=(8, 6))
plt.plot(accelX_01E, FIXe, color='black', linewidth=2)
plt.plot(accelX_01I, FIXi, color='purple', linewidth=2)
plt.plot(accelY_01E, FIYe, color='green', linewidth=2)
plt.plot(accelY_01I, FIYi, color='red', linewidth=2)
plt.xlabel('Acceleration Node 02 [mm/s²]', fontsize=10)
plt.ylabel('Inertia Force (fI) [N]', fontsize=10)
plt.title('Inertia Force (fI) vs Acceleration Curve', fontsize=10)
plt.legend(['ELASTIC IN X DIR.', 'INELASTIC IN X DIR.', 'ELASTIC IN Y DIR.', 'INELASTIC IN Y DIR.'])
plt.grid()
plt.show()

plt.figure(6, figsize=(8, 6))
plt.plot(dispX_01E, FSXe, color='black', linewidth=2)
plt.plot(dispX_01I, FSXi, color='purple', linewidth=2)
plt.plot(dispY_01E, FSYe, color='green', linewidth=2)
plt.plot(dispY_01I, FSYi, color='red', linewidth=2)
plt.xlabel('Displacement Node 02 [mm]', fontsize=10)
plt.ylabel('Spring Force (fS) [N]', fontsize=10)
plt.title('Spring Force (fS) vs Displacement', fontsize=10)
plt.legend(['ELASTIC IN X DIR.', 'INELASTIC IN X DIR.', 'ELASTIC IN Y DIR.', 'INELASTIC IN Y DIR.'])
plt.grid()
plt.show()

plt.figure(7, figsize=(8, 6))
plt.plot(dispX_01E, FDXe, color='black', linewidth=2)
plt.plot(dispX_01I, FDXi, color='purple', linewidth=2)
plt.plot(dispY_01E, FDYe, color='green', linewidth=2)
plt.plot(dispY_01I, FDYi, color='red', linewidth=2)
plt.xlabel('Displacement Node 02 [mm]', fontsize=10)
plt.ylabel('Damping Force (fD) [N]', fontsize=10)
plt.title('Damping Force (fD) vs Displacement Curve')
plt.legend(['ELASTIC IN X DIR.', 'INELASTIC IN X DIR.', 'ELASTIC IN Y DIR.', 'INELASTIC IN Y DIR.'])
plt.grid()
plt.show()
#%% ---------------------------------------------
# Print out the state of nodes 1 and 2
ops.printModel("node",1, 2)
# Print out the state of element 1
ops.printModel("ele", 1)
# Print the Model
#printModel()
ops.printModel("-JSON", "-file", "MDOF_PENDULUM_THREE.json")
#%%-------------------------------------------------------------------------------
# EXCEL OUTPUT
import pandas as pd

# Create DataFrame function
def create_df(dispXE, dispXI, velXE, velXI, accelXE, accelXI, reactionXE, reactionXI, FDXe, FDXi, FIXe, FIXi, FSXe, FSXi):
    df = pd.DataFrame({
        "DISPLACEMENT [mm] - ELASTIC": dispXE,
        "VELOCITY [mm/s] - ELASTIC": velXE,
        "ACCELERATION [mm/s^2] - ELASTIC": accelXE,
        "BASE-REACTION [N] - ELASTIC": reactionXE,
        "DAMPING-FORCE [N] - ELASTIC": FDXe,
        "INERTIA-FORCE [N] - ELASTIC": FIXe,
        "SPRING-FORCE [N] - ELASTIC": FSXe,
        "DISPLACEMENT [mm] - INELASTIC": dispXI,
        "VELOCITY [m/s] - INELASTIC": velXI,
        "ACCELERATION [mm/s^2] - INELASTIC": accelXI,
        "BASE-REACTION [N] - INELASTIC": reactionXI,
        "DAMPING-FORCE [N] - INELASTIC": FDXi,
        "INERTIA-FORCE [N] - INELASTIC": FIXi,
        "SPRING-FORCE [N] - INELASTIC": FSXi,   
    })
    return df


# Save to Excel
with pd.ExcelWriter("MDOF_PENDULUM_THREE_OUTPUT.xlsx", engine='openpyxl') as writer:
    
    df1 = create_df(dispX_02E, dispX_02I, velX_02E, velX_02I, accelX_02E, accelX_02I, reactionXE, reactionXI, FDXe, FDXi, FIXe, FIXi, FSXe, FSXi)
    df1.to_excel(writer, sheet_name="OUTPUT", index=False)    
#%%------------------------------------------------------------------------------------------------