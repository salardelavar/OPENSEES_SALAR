###########################################################################################################
#                    >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                     #
#                  NONLINEAR DYNAMIC ANALYSIS UNDER FREE-VIBRATION AND VISUALIZATION                      #
#     RESPONSE SPECTRA OF ACCELERATION, VELOCITY, DISPLACEMENT DUCTILITY DAMAGE INDEX USING OPENSEES      #
#---------------------------------------------------------------------------------------------------------#
#    INVESTIGATION OF THE EFFECTS OF MDOF STRUCTURAL DUCTILITY ON DAMAGE LEVEL AND STRUCTURAL DAMPING     #
#                                       USING NONLINEAR DYNAMIC ANALYSIS                                  #
#---------------------------------------------------------------------------------------------------------#
#                          CONSTANT STRUCTURAL DUCTILITY RATIO RESPONSE SPECTRUM                          #
#---------------------------------------------------------------------------------------------------------#
# EQUIVALENT SDOF SYSTEM DERIVATION VIA DISPLACEMENT-BASED SEISMIC DESIGN PROCEDURE WITH FREE-VIBRATION   #
#---------------------------------------------------------------------------------------------------------#
#IT IS RECOMMENEDED FOR BETTER SEEING THE STRUCTURAL DUCTILITY AND STRUCTURAL OVER-STRENGTH FACTOR EFFECTS#
# INTIAL DISPLACEMENT HAVE TO GREATER THAN YIELD DISPLACEMENT.                                            #
#---------------------------------------------------------------------------------------------------------#
#                          THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                     #
#                                   EMAIL: salar.d.ghashghaei@gmail.com                                   #
###########################################################################################################

"""
Investigation of the Effects of Structural Ductility on Damage Level and Structural Damping:
This code implements a comprehensive nonlinear dynamic analysis framework for
performance-based earthquake engineering assessment of single-degree-of-freedom
(SDOF) systems. The methodology combines traditional nonlinear time-history
analysis with modern probabilistic and machine learning techniques for advanced
structural performance evaluation with changing structural ductility raio and over strength factor.

KEY ENGINEERING OBJECTIVES:
1. Comparative assessment of hysteretic models for free-vibration response prediction
2. Probabilistic free-vibration demand analysis using multiple ground motions
3. Development of fragility curves for performance-based earthquake engineering
4. Integration of data science methods for structural reliability assessment

ANALYTICAL FEATURES:
- Nonlinear material behavior with pinching and degradation
- Response spectrum analysis across period range
- Real-time structural health monitoring metrics
- Statistical characterization of free-vibration demands
- Machine learning-based damage prediction
-----------------------------------
Model setup:
 - SDOF properties: mass (m), initial stiffness (k), yield displacement (Dy), ultimate displacement (Du), viscous damping (xi).
 - Hysteresis models: HYSTERETICSM (pinching, stiffness degradation, strength decay).
 - Damping: Rayleigh (or equivalent viscous) damping specified by target damping ratio xi for the fundamental mode.

Dynamic response:
 - Natural period T = 2*pi*sqrt(m/k) computed from linearized stiffness.
 - Time-history integration produces displacement, velocity, acceleration and base reaction histories.
 - HYSTERETIC model shows faster decay of amplitude and larger energy dissipation due to pinching and degradation.

Force–displacement behavior:
 - BILINEAR: symmetric hysteresis loops with stable post-yield stiffness; residual displacements are primarily due to plastic offset.
 - HYSTERETIC: pinched loops, reduced unloading/reloading stiffness, strength decay and larger residuals; captures cumulative damage effects.

Stiffness and strength evolution:
 - Effective lateral stiffness reduces during the excitation for both models but degrades faster with HYSTERETIC due to damage mechanisms.
 - Strength deterioration (reduced peak restoring force) in HYSTERETIC leads to reduced re-centering and larger residuals.

Damping estimation:
 - Use logarithmic decrement or energy-based measures from free vibration or post-event cycles.
 - HYSTERETIC typically yields higher equivalent damping (greater energy dissipation) compared with BILINEAR for the same displacement amplitude.

Peak responses:
 - Peak displacement: often lower for HYSTERETIC in early cycles because of softening, but long-term residual displacement may be higher.
 - Peak base shear (reaction): decays faster in HYSTERETICSM due to strength loss; BILINEAR sustains higher peak restoring forces for the same drift until hardening or limiting criteria apply.

Visualization:
 - Plot time histories (disp, vel, acc), hysteresis loops (force vs disp), and envelope curves to compare models.
 - Response spectra for displacement, velocity, and acceleration can be constructed from peak responses across parameter sweeps (e.g., varying T or post-yield stiffness).

Implications for free-vibration assessment:
 - BILINEAR: simple and computationally efficient; may overestimate resilience for severe cyclic demands because it omits degradation.
 - HYSTERETIC: captures important degradation mechanisms (pinching, stiffness/strength loss, ultimate strain) and is recommended for collapse assessment and detailed damage estimation.
 - Model selection should match the performance objective: serviceability checks might use simpler models; collapse and damage-sensitive studies require degraded hysteretic models calibrated with experiment.

Data export and post-processing:
 - Store peak and time-history results (displacement, velocity, acceleration, base reaction) to CSV/Excel for parametric studies.
 - Compute and plot response spectra (disp/vel/acc/reaction) from the stored peak values.

Ductility Damage Index (DDI) — implementation (concept):
 - After identifying yield displacement Dy and ultimate displacement Du from the capacity model:
   Dd = max(|disp_time_history|)  # maximum absolute dynamic displacement demand
   DI = (Dd - Dy) / (Du - Dy)      # Ductility Damage Index in the direction of interest
   Interpretation:
     DI <= 0   : elastic (no damage)
     0 < DI < 1: inelastic damage (serviceability/repairable)
     DI >= 1   : demand reaches or exceeds ultimate capacity (collapse or unacceptable damage)

Conclusions:
 - For inelastic SDOF studies, including pinching and degradation in the hysteretic model can change predicted peak responses, residual displacements, and damage indices significantly.
 - Use HYSTERETIC-type models for damage-sensitive or collapse-prone scenarios, and calibrate degradation parameters using test data where possible.
--------------------------------
Constant Structural Ductility Ratio
In free-vibration engineering, the structural ductility ratio (often denoted as μ) refers to the ability of a structure to undergo inelastic (plastic) deformations without collapsing during an earthquake. It is defined as the ratio of the maximum displacement (Δ_max) to the yield displacement (Δ_yield), i.e., μ = Δ_max / Δ_yield. This measures how much a structure can "stretch" beyond its elastic limit while dissipating energy through yielding, which helps prevent brittle failure. Ductile structures (e.g., those with μ > 5–6) can survive strong ground motions by deforming significantly, reducing the need for overly stiff designs.
A "constant structural ductility ratio" typically relates to constant ductility response spectra. In these spectra, the ductility ratio μ is held constant across different structural periods (T), and the required strength (e.g., pseudoacceleration or force) is plotted against the period to achieve that specific ductility demand. This approach contrasts with elastic response spectra, where no inelastic behavior is assumed. Constant ductility spectra are used to derive ductility reduction factors (Rμ or Rd), which adjust the elastic response spectrum downward to account for energy dissipation through ductility, making designs more economical. They are particularly useful in performance-based free-vibration design to ensure ductility supply exceeds demand, often calibrated using methods like Newmark-Hall for different period ranges (short: T < 0.2 s; intermediate: 0.2–0.5 s; long: T > 1 s)
------------------------------------------------ 
EQUIVALENT SDOF SYSTEM DERIVATION VIA DISPLACEMENT-BASED SEISMIC DESIGN PROCEDURE WITH PUSHOVER ANALYSIS:
Change MDOF to SDOF System with Displacement Based Design Concept

A displacement-based pushover transformation,
 converting a multi-degree-of-freedom (MDOF) system into an equivalent
 single-degree-of-freedom (SDOF) system for seismic assessment.
 It calculates effective modal properties—displacement, mass, and
 stiffness—by weighting element forces and nodal displacements according
 to a presumed deformed shape.
 The derived effective period provides a simplified dynamic characteristic
 for performance-based engineering. The visualizations effectively track the
 evolution of these equivalent parameters throughout the nonlinear analysis steps.
"""
# PAPER: Displacement-based seismic design of buildings theory - M.S. Medhekar, D.J.L. Kennedy - 1999 Elsevier
# YOUTUBE:
'https://www.youtube.com/watch?v=MZUhSHmIUdI'
#%%-------------------------------------------------------------------------------
import openseespy.opensees as ops
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time as ti
import SALAR_MATH as S01
import ANALYSIS_FUNCTION as S02
import MARKOV_CHAIN as S03
import RAYLEIGH_DAMPING_FUN as S04
import EIGENVALUE_ANALYSIS_FUN as S05
import DAMPING_RATIO_FUN as S06
import FRAGILITY_CURVE_FUN as S07
import OPENSEEES_HYSTERETICSM_FORCE_DISP_FUN as S08
import RAYLEIGH_DAMPING_FUN as S09
from scipy.stats import norm
#import COMPUTE_EFFECTIVE_PROPERTIES_FUN as S10

#%%------------------------------------------------------------------------------------------------
# Define parameters (units: m, N)
NUM_SIM = 15                                  # Total number for analysis
NUM_PERIOD = NUM_SIM                          # Total number for Period in each simulation
# NOTICE: NUM_PERIOD AND NUM_SIM MUST BE SAME VALUE
# BASED ON ABOVE DATA, TOTATL NUMBER OF SIMULATION IS (NUM_PERIOD * NUM_SIM)
#%%------------------------------------------------------------------------------------------------
# Define  Structural Properties
FY = 85000.0                                     # [N] Yield Force of Structure
FU = 1.18 * FY                                   # [N] Ultimate Force of Structure
Ke = 4500000.0                                   # [N/m] Spring Elastic Stiffness
DY = FY / Ke                                     # [m] Yield Displacement
DSU = 0.36                                       # [m] Ultimate Displacement
Ksh = (FU - FY) / (DSU - DY)                     # [N/m] Displacement Hardening Modulus
Kp = FU / DSU                                    # [N/m] Spring Plastic Stiffness
b = Ksh / Ke                                     # Displacement Hardening Ratio
M = 55000.0                                      # [kg] Mass of the Structure
DR = 0.03                                        # Damping Ratio


u0 = 0.0144       # [m] Initial displacement
v0 = 0.015        # [m/s] Initial velocity
a0 = 0.0065       # [m/s^2] Initial acceleration

IU = True         # Free Vibration with Initial Displacement
IV = False        # Free Vibration with Initial Velocity
IA = False        # Free Vibration with Initial Acceleration

duration = 2.0   # [s] Total simulation duration
dt = 0.01        # [s] Time step
#%%------------------------------------------------------------------------------------------------
if np.abs(u0) >= DY:
    print('Inelastic Behavior of Spring is going to be Run\n\n')
else:
    print('Elastic Behavior of Spring is going to be Run\n\n')
#%%------------------------------------------------------------------------------------------------
# ELASIC PERIOD:
ELAS_PERIOD = 2*np.pi * np.sqrt(M/Ke)
print(f'ELASIC PERIOD: {ELAS_PERIOD:.3f} (s)')     
# PLASIC PERIOD:
PLAS_PERIOD = 2*np.pi * np.sqrt(M/Kp)  
print(f'PLASIC PERIOD: {PLAS_PERIOD:.3f} (s)') 

# INITIAL MASS FOR RESPONSE SPECTRUM ANALYSIS
mi = (PLAS_PERIOD/2*np.pi)**2 * Kp 
print(mi)
#%%------------------------------------------------------------------------------------------------
# Calculate Over Strength Coefficient (Ω0)
Omega_0 = FU / FY
# Calculate Displacement Ductility Ratio (μ)
mu = DSU / DY
# Calculate Ductility Coefficient (Rμ)
R_mu = (2 * mu - 1) ** 0.5 / mu ** 0.5
# Calculate Structural Behavior Coefficient (R)
R = Omega_0 * R_mu
print(f'Over Strength Coefficient (Ω0):      {Omega_0:.2f}')
print(f'Displacement Ductility Ratio (μ):    {mu:.2f}')
print(f'Ductility Coefficient (Rμ):          {R_mu:.2f}')
print(f'Structural Behavior Coefficient (R): {R:.2f}')
#%%------------------------------------------------------------------------------------------------
# Define Analysis Properties
MAX_ITERATIONS = 1000000   # Convergence iteration for test
MAX_TOLERANCE = 1.0e-6    # Convergence tolerance for test

#%%------------------------------------------------------------------------------------------------
### OPENSEES FUNCTION
def ANALYSIS_MDOF(MAT_TYPE, TOTAL_MASS, duct, osf):
    # Initialize model
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)
    
    MAX_ITERATIONS = 5000   # Maximum number of iterations
    MAX_TOLERANCE = 1.0e-6  # Specified tolerance for convergence
    GMfact = 9.81 # [m/s^2] standard acceleration of gravity or standard acceleration 
    
    # [kg] Mass
    MASS = [0.40 * TOTAL_MASS, 0.20 * TOTAL_MASS, 0.30 * TOTAL_MASS, 0.10 * TOTAL_MASS]  
    # Damping ratio
    DRi = [0.05, 0.01, 0.02, 0.03]      
    
    # Define nodes
    ops.node(1, 0.0)  # Fixed base
    ops.node(2, 0.0)  # Mass node
    ops.node(3, 0.0)  # Mass node
    ops.node(4, 0.0)  # Mass node
    ops.node(5, 0.0)  # Mass node
        
    # Define boundary conditions
    ops.fix(1, 1)

    # Define masses
    for JJ in range(1, 4):
        ops.mass(JJ, MASS[JJ])
    
    #%% FOUR DEGREES OF FREEDOM STRUCTURE, CALCULATE LATERAL STIFFNESS AND DAMPING 
    for II in range(0, 4):
        FY = 85000.0                                     # [N] Yield Force of Structure
        FU = 1.18 * FY                                   # [N] Ultimate Force of Structure
        Ke = 4500000.0 * (4 - II)                        # [N/m] Spring Elastic Stiffness
        #DY = FY / Ke                                     # [m] Yield Displacement
        #DSU = 0.36                                       # [m] Ultimate Displacement
        DSU = DY * duct # IN EACH STEP IT WILL BE CHNAGED
        FU = FY * osf   # IN EACH STEP IT WILL BE CHNAGED
        Ksh = (FU - FY) / (DSU - DY)                     # [N/m] Displacement Hardening Modulus
        Kp = FU / DSU                                    # [N/m] Spring Plastic Stiffness
        b = Ksh / Ke                                     # Displacement Hardening Ratio

        # Define material properties
        MAT_TAG = 1000 + II # SPRING TAG
        # FORCE-DISPLACEMENT RELATIONSHIP OF LATERAL SPRING AND PLOT 
        DP = [0, 0, 0, 0]
        FP = [0, 0, 0, 0]
        DN = [0, 0, 0, 0]
        FN = [0, 0, 0, 0]
        #print(DSU,"------------" ,FU)
        DP[0], FP[0] = DY, FY
        DP[1], FP[1] = DSU, FU 
        DP[2], FP[2] = 1.1*DSU, 0.20*FU
        DP[3], FP[3] = 1.25*DSU, 0.10*FU
        DN[0], FN[0] = -DY, -FY 
        DN[1], FN[1] = -DSU, -FU
        DN[2], FN[2] = -1.1*DSU, -0.20*FU   
        DN[3], FN[3] = -1.25*DSU, -0.10*FU
        #print(DP, FP)
        #print(DN, FN)
        
        if MAT_TYPE == 'INELASTIC':
            S08.OPENSEEES_HYSTERETICSM_FORCE_DISP_FUN(MAT_TAG, DP, FP, DN, FN, PLOT=False, X_LABEL='Displacement (mm)', Y_LABEL='Force [N]', TITLE='FORCE-DISPLACEMENT CURVE')
        if MAT_TYPE == 'ELASTIC':
            #ops.uniaxialMaterial('Elastic', MAT_TAG, Ke)             # TESNSION AND COMPRESSION IS SAME VALUES
            ops.uniaxialMaterial('Elastic', MAT_TAG, Ke ,0.0, 0.5*Ke) # TESNSION AND COMPRESSION IS NOT SAME VALUES
            # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/ElasticUni.html
        
        MatTag_C = 2000 + II # SPRING DAMPER
        alpha = 1.0    # velocity exponent (usually 0.3–1.0)
        omega = np.sqrt(Ke / MASS[II])
        Cd = 2 * DRi[II] * omega * MASS[II]  # [N·s/m] Damping coefficient 
        ops.uniaxialMaterial('Viscous', MatTag_C, Cd, alpha)  # Material for C (alpha=1.0 for linear)
        
        # Define element
        ops.element('zeroLength', II+1, II+1, II+2, '-mat',  MAT_TAG, MatTag_C, '-dir', 1, 1)  # DOF LATERAL SPRING

    center_node = 5
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(2, 1.0)
    ops.load(3, 1.0)
    ops.load(4, 1.0)
    ops.load(center_node, 1.0)
    
    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.system('BandGeneral')
    
    if MAT_TYPE == 'ELASTIC':
        ops.algorithm('Linear')
    if MAT_TYPE == 'INELASTIC':
        ops.algorithm('Newton')
        
    time = []
    disp = []
    velo = []
    acc = []
    reaction = []
    stiffness = []
    OMEGA, PERIOD = [], []
    PERIOD_MIN, PERIOD_MAX = [], []
    DI = []
    ele_force = {
        1: [],  # AXIALFORCE-01
        2: [],  # AXIALFORCE-02
        3: [],  # AXIALFORCE-03
        4: [],  # AXIALFORCE-04
        }
    # Initialize lists for each node's displacement
    node_displacements = {
        1: [],  # DISP01
        2: [],  # DISP02
        3: [],  # DISP03
        4: [],  # DISP04
        5: [],  # DISP05
        }
    
    if IU == True:
        # Define initial displacment
        ops.setNodeDisp(center_node, 1, u0, '-commit') # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/setNodeDisp.html
    if IV == True:
        # Define initial velocity
        ops.setNodeVel(center_node, 1, v0, '-commit')  # INFO LINK: https://openseespydoc.readthedocs.io/en/stable/src/setNodeVel.html
    if IA == True:
        # Define initial  acceleration
        ops.setNodeAccel(center_node, 1, a0, '-commit') # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/setNodeAccel.html
      
    # Set analysis parameters
    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.system('BandGeneral')
    ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/normDispIncr.html
    #ops.integrator('CentralDifference')  # JUST FOR LINEAR ANALYSIS - INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/centralDifference.html
    alpha=0.5; beta=0.5;
    ops.integrator('Newmark', alpha, beta) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/newmark.html
    #alpha=2/3;gamma=1.5-alpha; gamma=1.5-alpha;beta=(2-alpha)**2/4;
    #ops.integrator('HHT', alpha, gamma, beta) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/hht.html
    ops.algorithm('Newton')  # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/algorithm.html
    ops.analysis('Transient') # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/analysis.html
    
    # Calculate Rayleigh damping factors
    Lambda01 = ops.eigen('-fullGenLapack', 1) # eigenvalue mode 1
    #Lambda01 = ops.eigen('-genBandArpack', 1) # eigenvalue mode 1
    Omega01 = np.power(max(Lambda01), 0.5)
    A0 = (2 * Omega01 * DR) / Omega01 # c = a0 * m : Mass-proportional damping
    A1 = (DR * 2) / Omega01 # c = a1 * k : Stiffness-proportional damping
    # Apply Rayleigh damping
    ops.rayleigh(A0, 0.0, 0.0, A1)# INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/reyleigh.html
    #ops.rayleigh(0, 0, 2 * DR * Omega01, 0) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/reyleigh.html
    PERIOD = 2*np.pi / Omega01   # Structure Period 
    print(f'PERIOD: {PERIOD:.6f}')
    
    # Calculate Rayleigh damping factors
    #PERIOD_01, PERIOD_02 = S04.RAYLEIGH_DAMPING(1, 0.5*DR, DR, 0, 1)
    
    # Run Dynamic Analysis
    time = []
    disp_02, disp_03, disp_04, disp_05 = [], [], [], []
    velo_02, velo_03, velo_04, velo_05 = [], [], [], []
    accel_02, accel_03, accel_04, accel_05 = [], [], [], []
    base_reaction = []
    DI = []
    stiffness = []
    PERIOD_MIN, PERIOD_MAX = [], []

        
    stable = 0
    current_time = 0.0
    step = 0
    while stable == 0 and current_time < duration:
        stable = ops.analyze(1, dt)
        S02.ANALYSIS(stable, 1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS
        current_time = ops.getTime()
        time.append(current_time)
        disp_02.append(ops.nodeDisp(2, 1)); velo_02.append(ops.nodeVel(2, 1)); accel_02.append(ops.nodeAccel(2, 1));
        disp_03.append(ops.nodeDisp(3, 1)); velo_03.append(ops.nodeVel(3, 1)); accel_03.append(ops.nodeAccel(3, 1));
        disp_04.append(ops.nodeDisp(4, 1)); velo_04.append(ops.nodeVel(4, 1)); accel_04.append(ops.nodeAccel(4, 1));
        disp_05.append(ops.nodeDisp(5, 1)); velo_05.append(ops.nodeVel(5, 1)); accel_05.append(ops.nodeAccel(5, 1));
        ops.reactions()
        base_reaction.append(ops.nodeReaction(1, 1))                 # Reaction force
        stiffness.append(np.abs(base_reaction[-1]/disp_05[-1]))      # [N/m] Stiffness
        DI.append(100*(disp_05[-1] - DY) / (DSU - DY))               # [%] Structural Ductility Damage Index 
        if DI[-1] <= 0.0: # DAMAGE INDEX -> 0.0 <= DI <=100 
            DI[-1] = 0.0
        if DI[-1] >= 100: 
            DI[-1] = 100.0
        # Store forces and displacements
        for ele_id in ele_force.keys(): 
            ele_force[ele_id].append(ops.eleResponse(ele_id, 'force')[0])        # [N] ELEMENT AXIAL FORCE       # [N] ELEMENT MOMENT FORCE                
        # Store displacements
        for node_id in node_displacements.keys():    
            node_displacements[node_id].append(ops.nodeDisp(node_id, 1))        
        # IN EACH STEP STRUCTURAL PERIOD WILL BE CALCULATED
        #PERIODmin, PERIODmax = S09.RAYLEIGH_DAMPING(2, 0.5*DR, DR, 0, 1)
        PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(2, PLOT=True)
        PERIOD_MIN.append(PERIODmin)
        PERIOD_MAX.append(PERIODmax)
        step += 1
        #print(f'{time[-1]} {displacement[-1]:0.4e} {base_reaction[-1]:0.4e}')
        
    # Calculate Structure Damping Ratio Based on Lateral Displacement
    damping_ratio_02 = S06.DAMPING_RATIO(disp_02) # DAMAPING RATIO FROM DOF 02
    damping_ratio_03 = S06.DAMPING_RATIO(disp_03) # DAMAPING RATIO FROM DOF 03
    damping_ratio_04 = S06.DAMPING_RATIO(disp_04) # DAMAPING RATIO FROM DOF 04
    damping_ratio_05 = S06.DAMPING_RATIO(disp_05) # DAMAPING RATIO FROM DOF 05
    
    # Run the file loading effective properties
    #exec(open("COMPUTE_EFFECTIVE_PROPERTIES_FUN.py").read())
    #S10.COMPUTE_EFFECTIVE_PROPERTIES_FUN
    
    #%% EQUIVALENT SDOF SYSTEM DERIVATION VIA DISPLACEMENT-BASED SEISMIC DESIGN PROCEDURE WITH PUSHOVER ANALYSIS
    # Change MDOF to SDOF System with Displacement Based Design Concept
    # THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
    """
    This script implements a displacement-based pushover transformation,
     converting a multi-degree-of-freedom (MDOF) system into an equivalent
     single-degree-of-freedom (SDOF) system for seismic assessment.
     It calculates effective modal properties—displacement, mass, and
     stiffness—by weighting element forces and nodal displacements according
     to a presumed deformed shape.
     The derived effective period provides a simplified dynamic characteristic
     for performance-based engineering. The visualizations effectively track the
     evolution of these equivalent parameters throughout the nonlinear analysis steps.
    """
    #---------------------------------------------------------------------------

    displacement_X_1 = np.array(list(node_displacements.values())[1])   # DOF 02    
    displacement_X_2 = np.array(list(node_displacements.values())[2])   # DOF 03 
    displacement_X_3 = np.array(list(node_displacements.values())[3])   # DOF 04 
    displacement_X_4 = np.array(list(node_displacements.values())[4])   # DOF 05 

    ele_force_01 = np.array(list(ele_force.values())[0])   # ELEMENT 01 
    ele_force_02 = np.array(list(ele_force.values())[1])   # ELEMENT 02 
    ele_force_03 = np.array(list(ele_force.values())[2])   # ELEMENT 03 
    ele_force_04 = np.array(list(ele_force.values())[3])   # ELEMENT 04 

    STIFF_X_1 = np.abs(ele_force_01 / displacement_X_1)
    STIFF_X_2 = np.abs(ele_force_02 / displacement_X_2)
    STIFF_X_3 = np.abs(ele_force_03 / displacement_X_3)
    STIFF_X_4 = np.abs(ele_force_04 / displacement_X_4)

    MX2 = (MASS[0] * np.square(displacement_X_1) + 
           MASS[1] * np.square(displacement_X_2) + 
           MASS[2] * np.square(displacement_X_3) +
           MASS[3] * np.square(displacement_X_4))

    MX = (MASS[0] * np.array(displacement_X_1) + 
          MASS[1] * np.array(displacement_X_2) + 
          MASS[2] * np.array(displacement_X_3) + 
          MASS[3] * np.array(displacement_X_4))
    EFFECTIVE_DISP_X = MX2 / np.abs(MX)

    EFFECTIVE_MASS_X = np.abs(MX) / EFFECTIVE_DISP_X


    # Effective Stiffness
    KX = (np.array(STIFF_X_1) * np.array(displacement_X_1) + 
          np.array(STIFF_X_2) * np.array(displacement_X_2) + 
          np.array(STIFF_X_3) * np.array(displacement_X_3) + 
          np.array(STIFF_X_4) * np.array(displacement_X_4))

    EFFECTIVE_STIFF_X = np.abs(KX) / EFFECTIVE_DISP_X


    # Effective Period
    EFFECTIVE_PERIOD_X = 2 * np.pi / np.sqrt(EFFECTIVE_STIFF_X/EFFECTIVE_MASS_X)

    print('Median Effective Displacement:           ', np.median(EFFECTIVE_DISP_X))
    print('Median Effective Mass:                   ', np.median(EFFECTIVE_MASS_X))
    print('Median Effective Stiffness:              ', np.median(EFFECTIVE_STIFF_X))
    print('Median Effective Period:                 \n', np.median(EFFECTIVE_PERIOD_X))

    # Effective Damping Ratio
    DRX = (damping_ratio_02 * np.array(displacement_X_1) + 
          damping_ratio_03 * np.array(displacement_X_2) + 
          damping_ratio_04 * np.array(displacement_X_3) + 
          damping_ratio_05 * np.array(displacement_X_4))

    EFFECTIVE_DR_X = np.abs(DRX / EFFECTIVE_DISP_X)

    MED_MASS = np.median(EFFECTIVE_MASS_X)
    MED_STIFF = np.median(EFFECTIVE_STIFF_X)
    MED_PERIOD = np.median(EFFECTIVE_PERIOD_X)
    MED_DR = np.median(EFFECTIVE_DR_X)

    # Create a figure with two subplots (Effective Displacement and Effective Mass)
    fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, 1, figsize=(18, 12))

    ax1.plot(time, EFFECTIVE_DISP_X, color='green', linewidth=3)
    ax1.set_title(f'Effective Displacement - Median: {np.median(EFFECTIVE_DISP_X): .5f}')
    ax1.set_xlabel('Time [s]')
    ax1.set_ylabel('Effective Displacement')
    #ax1.legend(loc='upper right')
    ax1.grid(True)

    ax2.plot(time, EFFECTIVE_MASS_X, color='magenta', linewidth=3)
    ax2.set_title(f'Effective Mass - Median: {np.median(EFFECTIVE_MASS_X): .5f}')
    ax2.set_xlabel('Time [s]')
    ax2.set_ylabel('Effective Mass')
    #ax2.legend(loc='upper right')
    ax2.grid(True)

    ax3.plot(time, EFFECTIVE_STIFF_X, color='cyan', linewidth=3)
    ax3.set_title(f'Effective Stiffness - Median: {np.median(EFFECTIVE_STIFF_X): .5f}')
    ax3.set_xlabel('Time [s]')
    ax3.set_ylabel('Effective Stiffness')
    #ax3.legend(loc='upper right')
    ax3.grid(True)

    ax4.plot(time, EFFECTIVE_PERIOD_X, color='black', linewidth=3)
    ax4.set_title(f'Effective Period - Median: {np.median(EFFECTIVE_PERIOD_X): .5f}')
    ax4.set_xlabel('Time [s]')
    ax4.set_ylabel('Effective Period')
    #ax4.legend(loc='upper right')
    #ax4.semilogy()
    ax4.grid(True)

    ax5.plot(time, EFFECTIVE_DR_X, color='red', linewidth=3)
    ax5.set_title(f'Effective Damping Ratio - Median: {np.median(EFFECTIVE_DR_X): .2f} [%]')
    ax5.set_xlabel('Time [s]')
    ax5.set_ylabel('Effective Damping Ratio')
    #ax5.legend(loc='upper right')
    #ax5.semilogy()
    ax5.grid(True)

    plt.tight_layout()
    plt.show()
    
    #%%--------------------
    
    DATA = (time, 
            disp_02,
            disp_03,
            disp_04,
            disp_05,
            velo_02,
            velo_03,
            velo_04,
            velo_05,
            accel_02,
            accel_03,
            accel_04,
            accel_05,
            base_reaction, DI,
            PERIOD, stiffness,
            damping_ratio_02,
            damping_ratio_03, 
            damping_ratio_04, 
            damping_ratio_05, 
            ele_force, node_displacements, 
            PERIOD_MIN, PERIOD_MAX,
            MED_MASS, MED_STIFF, MED_PERIOD, MED_DR)
    
    return DATA

#%%------------------------------------------------------------------------------------------------
# Analysis Durations:
starttime = ti.process_time()

# Collect Data
DATA = {
    1: [],  # DISP
    2: [],  # VELOCITY
    3: [],  # ACCELERAION
    4: [],  # REACTION
    5: [],  # DI
    6: [],  # DAMPING RATIO
    7: [],  # STRUCTURAL STIFFNESS
    8: [],  # STRUCTURAL DUCTILITY RATIO
    9: [],  # STRUCTURAL OVER-STRENGTH FACTOR
    10:[],  # EFFECTIVE MASS
    11:[],  # EFFECTIVE STIFFNESS
    12:[],  # STRUCTURAL PERIOD
    13:[],  # STRUCTURAL DAMPING RATIO
    }

# NUM_SIM is the number of simulations
MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'
duct = np.linspace(2.0, 8.5, NUM_SIM).tolist()  # DEFINE MIN. AND MAX. VALUE FOR STRUCURAL DUCTILITY RATIO
osf = np.linspace(1.05, 1.35, NUM_SIM).tolist() # DEFINE MIN. AND MAX. VALUE FOR STRUCURAL OVER STRENGTH FACTOR

for j in range(NUM_SIM):
    # Initialize lists to store max values
    max_time = []
    max_disp_02, max_disp_03, max_disp_04, max_disp_05 = [], [], [], []
    max_velo_02, max_velo_03, max_velo_04, max_velo_05 = [], [], [], []
    max_accel_02, max_accel_03, max_accel_04, max_accel_05 = [], [], [], []
    max_dr_02, max_dr_03, max_dr_04, max_dr_05  = [], [], [], []
    max_base_reaction = []
    max_DI = []
    max_T = []
    max_K = []
    max_duct = []
    max_osf = []
    max_MED_MASS, max_MED_STIFF, max_MED_PERIOD, max_MED_DR = [], [], [], []
    print('Structure Ductility Ratio: ', duct[j], " ------ ",'Structure Over-strength Factor: ',osf[j])
    for i in range(NUM_PERIOD):
        m = mi * (i+1 / NUM_PERIOD) * 0.04
        data = ANALYSIS_MDOF(MAT_TYPE, m, duct[j], osf[j])
        (time, 
        disp_02,
        disp_03,
        disp_04,
        disp_05,
        velo_02,
        velo_03,
        velo_04,
        velo_05,
        accel_02,
        accel_03,
        accel_04,
        accel_05,
        base_reaction, DI,
        PERIOD, stiffness,
        damping_ratio_02,
        damping_ratio_03, 
        damping_ratio_04, 
        damping_ratio_05, 
        ele_force, node_displacements, 
        PERIOD_MIN, PERIOD_MAX,
        MED_MASS, MED_STIFF, MED_PERIOD, MED_DR) = data
        # Calculate and store the max absolute values
        max_time.append(np.max(np.abs(time)))
        max_disp_02.append(np.max(np.abs(disp_02))); max_disp_03.append(np.max(np.abs(disp_03))); max_disp_04.append(np.max(np.abs(disp_04))); max_disp_05.append(np.max(np.abs(disp_05)));
        max_velo_02.append(np.max(np.abs(velo_02))); max_velo_03.append(np.max(np.abs(velo_03))); max_velo_04.append(np.max(np.abs(velo_04))); max_velo_05.append(np.max(np.abs(velo_05))); 
        max_accel_02.append(np.max(np.abs(accel_02))); max_accel_03.append(np.max(np.abs(accel_03))); max_accel_04.append(np.max(np.abs(accel_04))); max_accel_05.append(np.max(np.abs(accel_05))); 
        max_base_reaction.append(np.max(base_reaction))
        max_DI.append(np.max(DI))
        max_K.append(np.max(stiffness))
        #max_K.append(0.0)
        max_dr_02.append(damping_ratio_02);
        max_dr_03.append(damping_ratio_03);
        max_dr_04.append(damping_ratio_04);
        max_dr_05.append(damping_ratio_05);
        max_duct.append(duct[j])
        max_osf.append(osf[j])  
        max_T.append(PERIOD)
        max_MED_MASS.append(MED_MASS)
        max_MED_STIFF.append(MED_STIFF)
        max_MED_PERIOD.append(MED_PERIOD)
        max_MED_DR.append(MED_DR)
        print(f'STEP {i + 1} DONE')
        
    # Store Data
    DATA[1].append(max_disp_05)
    DATA[2].append(max_velo_05)
    DATA[3].append(max_accel_05)
    DATA[4].append(max_base_reaction)
    DATA[5].append(max_DI)   
    DATA[6].append(max_dr_05)
    DATA[7].append(max_K)
    DATA[8].append(max_duct)
    DATA[9].append(max_osf)
    DATA[10].append(max_MED_MASS)
    DATA[11].append(max_MED_STIFF)
    DATA[12].append(max_MED_PERIOD)
    DATA[13].append(max_MED_DR)
    print(f'\n SIMULATION {j + 1} DONE \n')    
else:
    print('\n\n Analysis completed successfully')    

totaltime = ti.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n') 

#%%------------------------------------------------------------------------------------------------
# OUTPUT DATA FROM RESPONSE SPECTRUM ANALYSIS TO EXCEL FILE
labels = {
    1: "DISPLACEMENT",
    2: "VELOCITY",
    3: "ACCELERATION",
    4: "BASE_REACTION",
    5: "DAMAGE_INDEX",
    6: "DAMPING_RATIO",
    7: "STIFFNESS",
    8: "DUCTILITY_RATIO",
    9: "OVER_STRENGTH_FACTOR",
    10: "EFFECTIVE_MASS",
    11: "EFFECTIVE_STIFFNESS",
    12: "EFFECTIVE_PERIOD",
    13: "EFFECTIVE_DAMPING_RATIO",
}

with pd.ExcelWriter("MDOF_RESPONSE_SPECTRUM_FREE_VIBRATION_DUCT_OSF_RESULTS.xlsx", engine="openpyxl") as writer:
    for key, value in DATA.items():
        df = pd.DataFrame(value)
        sheet_name = labels[key][:31] # SHEEET NAME
        df.to_excel(writer, sheet_name=sheet_name, index=False)

print("Excel file saved: SDOF_RESULTS.xlsx")

#%%------------------------------------------------------------------------------------------------
# PLOT THE RESPONSE SPECTRUMS 
import matplotlib.pyplot as plt
import numpy as np

# Convert to numpy arrays for safety
T = np.array(max_T)

""""
def PLOT_2D(COUNT, X, Y, XLABEL, YLABEL, TITLE):
    plt.figure(COUNT, figsize=(12,10))
    for j in range(NUM_SIM):
        plt.plot(T, DATA[COUNT][j],label = f'Max: {np.max(np.abs(DATA[COUNT][j])):.4e}')
    plt.xlabel(XLABEL)
    plt.ylabel(YLABEL)
    plt.title(TITLE)
    plt.legend()
    plt.grid(True)
    plt.show()
"""    
def PLOT_2D(COUNT, XLABEL, YLABEL, TITLE, SEMILOGY=False):
    plt.figure(COUNT, figsize=(12, 10))

    # Plot all simulations
    for j in range(NUM_SIM):
        #plt.plot(T, DATA[COUNT][j], alpha=0.4, label=f'Sim {j+1} | Max: {np.max(np.abs(DATA[COUNT][j])):.4e}')
        plt.plot(T, DATA[COUNT][j], alpha=0.4)
    
    # Convert to NumPy array for statistical calculations
    arr = np.array(DATA[COUNT])    # Shape: (NUM_SIM, len(T))

    # Compute statistical metrics
    mean_curve   = np.mean(arr, axis=0)
    median_curve = np.median(arr, axis=0)
    std_curve    = np.std(arr, axis=0)
    std_curveDOWN = mean_curve - np.std(arr, axis=0)
    std_curveUP = mean_curve + np.std(arr, axis=0)
    q1_curve = np.percentile(arr, 25, axis=0)
    q3_curve = np.percentile(arr, 75, axis=0)

    # Plot statistical curves
    plt.plot(T, mean_curve, color='navy', linestyle='-.', linewidth=3, label='Mean')
    
    #plt.plot(T, std_curveDOWN, color='steelblue', linestyle='-', linewidth=2, label='Mean - Std')
    #plt.plot(T, std_curveUP,   color='deepskyblue', linestyle='-', linewidth=2, label='Mean + Std')
    
    plt.plot(T, median_curve, color='red', linestyle='--', linewidth=2, label='Median')
    
    plt.plot(T, q1_curve, color='green', linestyle='-.', linewidth=2, label='Q1 (25%)')
    plt.plot(T, q3_curve, color='orange', linestyle='-.', linewidth=2, label='Q3 (75%)')

    # Standard deviation band around the mean
    plt.fill_between(T, mean_curve - std_curve, mean_curve + std_curve,
                     color='gray', alpha=0.25, label='Mean ± Std')

    # Labels and title
    plt.xlabel(XLABEL)
    plt.ylabel(YLABEL)
    plt.title(TITLE)
    
    if SEMILOGY == True:
        plt.semilogy()
        
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()
    return mean_curve, median_curve

    
    
# ----------------------------------------
# Plot 1: Displacement Response Spectrum
# ----------------------------------------
X = T
Y = DATA
XLABEL = "Period T [s]"
YLABEL = "Max Displacement (PGD) [m]"
TITLE =  "Displacement Response Spectrum"
SEMILOGY = True

mean_disp, median_disp = PLOT_2D(1, XLABEL, YLABEL, TITLE, SEMILOGY)

# ----------------------------------------
# Plot 2: Velocity Response Spectrum
# ----------------------------------------
X = T
Y = DATA
XLABEL = "Period T [s]"
YLABEL = "Max Velocity (PGV) [m/s]"
TITLE =  "Velocity Response Spectrum"
SEMILOGY = True

mean_velo, median_velo = PLOT_2D(2, XLABEL, YLABEL, TITLE, SEMILOGY)

# ----------------------------------------
# Plot 3: Acceleration Response Spectrum
# ----------------------------------------
X = T
Y = DATA
XLABEL = "Period T [s]"
YLABEL = "Max Acceleration (PGA) [m/s²]"
TITLE =  "Acceleration Response Spectrum"
SEMILOGY = True

mean_acce, median_acce = PLOT_2D(3, XLABEL, YLABEL, TITLE, SEMILOGY)

# ----------------------------------------
# Plot 4: Base Reaction Spectrum
# ----------------------------------------
X = T
Y = DATA
XLABEL = "Period T [s]"
YLABEL = "Max Base Reaction (N)"
TITLE =  "Base Reaction Response Spectrum"
SEMILOGY = True

mean_reaction, median_reaction = PLOT_2D(4, XLABEL, YLABEL, TITLE, SEMILOGY)

# ----------------------------------------
# Plot 5: Damage Index Spectrum
# ----------------------------------------
X = T
Y = DATA
XLABEL = "Period T [s]"
YLABEL = "Damage Index (DI) [%]"
TITLE =  "Ductility Damage Index Spectrum"
SEMILOGY = False

mean_di, median_di = PLOT_2D(5, XLABEL, YLABEL, TITLE, SEMILOGY)

# ----------------------------------------
# Plot 6: Damping Ratio Spectrum
# ----------------------------------------
X = T
Y = DATA
XLABEL = "Period T [s]"
YLABEL = "Damping Ratio [%]"
TITLE =  "Damping Ratio Spectrum"
SEMILOGY = False

mean_da, median_da = PLOT_2D(6, XLABEL, YLABEL, TITLE, SEMILOGY)

# ----------------------------------------
# Plot 7: Structural Stiffness Spectrum
# ----------------------------------------
X = T
Y = DATA
XLABEL = "Period T [s]"
YLABEL = "Structural Stiffness Spectrum [N/m]"
TITLE =  "Structural Stiffness Spectrum"
SEMILOGY = False

mean_stif, median_stif = PLOT_2D(7, XLABEL, YLABEL, TITLE, SEMILOGY)

# ----------------------------------------
# Plot 8: Effective Mass Spectrum
# ----------------------------------------
X = T
Y = DATA
XLABEL = "Period T [s]"
YLABEL = "Effective Mass Spectrum [kg]"
TITLE =  "Effective Mass Spectrum"
SEMILOGY = False

mean_efMASS, median_efMASS = PLOT_2D(10, XLABEL, YLABEL, TITLE, SEMILOGY)

# ----------------------------------------
# Plot 9: Effective Stiffness Spectrum
# ----------------------------------------
X = T
Y = DATA
XLABEL = "Period T [s]"
YLABEL = "Effective Stiffness Spectrum [N/m]"
TITLE =  "Effective Stiffness Spectrum"
SEMILOGY = False

mean_efSTIF, median_efSTIF = PLOT_2D(11, XLABEL, YLABEL, TITLE, SEMILOGY)

# ----------------------------------------
# Plot 10: Effective Period Spectrum
# ----------------------------------------
X = T
Y = DATA
XLABEL = "Period T [s]"
YLABEL = "Effective Period Spectrum [s"
TITLE =  "Effective Period Spectrum"
SEMILOGY = False

mean_efPER, median_efPER = PLOT_2D(12, XLABEL, YLABEL, TITLE, SEMILOGY)

# ------------------------------------------
# Plot 11: Effective Damping Ratio Spectrum
# ------------------------------------------
X = T
Y = DATA
XLABEL = "Period T [s]"
YLABEL = "Effective Damping Ratio Spectrum [%]"
TITLE =  "Effective Damping Ratio Spectrum"
SEMILOGY = False

mean_efDR, median_efDR = PLOT_2D(13, XLABEL, YLABEL, TITLE, SEMILOGY)

#%%------------------------------------------------------------------------------
# --------------------------------------------------------------
#   Plot Structural Ductility Ratio and Structural Damage Index
# --------------------------------------------------------------
plt.figure(1, figsize=(12, 8))
plt.plot(duct, median_di, color='purple', linewidth=2)
plt.xlabel('Structural Ductility Ratio [m/m]')
plt.ylabel('Structural Damage Index (DI) [%]')
plt.title('Structural Damage Index vs Structural Ductility Ratio')
plt.legend(['MEDIAN DDI'])
plt.grid()
plt.show()
# -----------------------------------------------------
#   Plot Structural Ductility Ratio and Damping Ratio
# -----------------------------------------------------
plt.figure(2, figsize=(12, 8))
plt.plot(duct, median_da, color='black', linewidth=2)
plt.xlabel('Structural Ductility Ratio [m/m]')
plt.ylabel('Damping Ratio [%]')
plt.title('Damping Ratio vs Structural Ductility Ratio')
plt.legend(['MEDIAN DAMPING RATIO'])
plt.grid()
plt.show()
#%%------------------------------------------------------------------------------
def PLOT_3D_CONTOUR_XYZ(TAG, X, Y, Z, XLABEL, YLABEL, ZLABEL):
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
    zi = griddata((X, Y), Z, (xi, yi), method='nearest')
    #zi = griddata((X, Y), Z, (xi, yi), method='linear')
    #zi = griddata((X, Y), Z, (xi, yi), method='cubic')
    
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
X = duct
Y = osf
Z = median_disp
XLABEL = 'Structural Ductility Ratio [m/m]'   
YLABEL = 'Structural Over-Strength Factor [N/N]' 
ZLABEL = "Max Displacement (PGD) [m]"
PLOT_3D_CONTOUR_XYZ(0, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = duct
Y = osf
Z = median_velo
XLABEL = 'Structural Ductility Ratio [m/m]'   
YLABEL = 'Structural Over-Strength Factor [N/N]' 
ZLABEL = "Max Velocity (PGV) [m/s]"
PLOT_3D_CONTOUR_XYZ(1, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = duct
Y = osf
Z = median_acce
XLABEL = 'Structural Ductility Ratio [m/m]'   
YLABEL = 'Structural Over-Strength Factor [N/N]' 
ZLABEL = "Max Acceleration (PGA) [m/s²]"
PLOT_3D_CONTOUR_XYZ(2, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = duct
Y = osf
Z = median_reaction
XLABEL = 'Structural Ductility Ratio [m/m]'   
YLABEL = 'Structural Over-Strength Factor [N/N]' 
ZLABEL = "Max Base-rection (PGA) [N]"
PLOT_3D_CONTOUR_XYZ(3, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = duct
Y = osf
Z = median_di
XLABEL = 'Structural Ductility Ratio [m/m]'   
YLABEL = 'Structural Over-Strength Factor [N/N]' 
ZLABEL = "Sructural Ductility Damage Index (DI) [%]"
PLOT_3D_CONTOUR_XYZ(4, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = duct
Y = osf
Z = median_efMASS
XLABEL = 'Structural Ductility Ratio [m/m]'   
YLABEL = 'Structural Over-Strength Factor [N/N]' 
ZLABEL = "Effective Mass [kg]"
PLOT_3D_CONTOUR_XYZ(5, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = duct
Y = osf
Z = median_efSTIF
XLABEL = 'Structural Ductility Ratio [m/m]'   
YLABEL = 'Structural Over-Strength Factor [N/N]' 
ZLABEL = "Effective Stiffness [N/m]"
PLOT_3D_CONTOUR_XYZ(6, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = duct
Y = osf
Z = median_efPER
XLABEL = 'Structural Ductility Ratio [m/m]'   
YLABEL = 'Structural Over-Strength Factor [N/N]' 
ZLABEL = "Effective Period [s]"
PLOT_3D_CONTOUR_XYZ(7, X, Y, Z, XLABEL, YLABEL, ZLABEL)

X = duct
Y = osf
Z = median_efDR
XLABEL = 'Structural Ductility Ratio [m/m]'   
YLABEL = 'Structural Over-Strength Factor [N/N]' 
ZLABEL = "Effective Damping Ratio [%]"
PLOT_3D_CONTOUR_XYZ(8, X, Y, Z, XLABEL, YLABEL, ZLABEL)
#%%------------------------------------------------------------------------------------------------
####  FRAGILITY ANALYSIS
  
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
"""
im_values = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
probabilities = {
    'DS1': [0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99],
    'DS2': [0.0, 0.01, 0.05, 0.1, 0.3, 0.5, 0.7, 0.85, 0.95, 0.99],
    'DS3': [0.0, 0.0, 0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.95]
}
"""  
# --------------
# Visualization
# --------------
plt.figure(1, figsize=(10, 6))
# Response plot
plt.plot(time, accel_05, lw=1, color='black')
plt.xlabel('Time (s)')
plt.ylabel('Acceleration (g)')
plt.title(f'Last Analysis Structural Response + Ground Motion ::: MAX. ABS. : {np.max(np.abs(accel_05)):.4f}')
plt.grid(True)
plt.show()    

im_values = median_acce 
X_LABEL = 'Peak Ground Acceleration (g)  [IM]'
S07.FRAGILITY_CURVE(im_values, damage_states, X_LABEL, SEMILOGY=False, PLOT_KIND=True)
#===========================================================
# Define damage state parameters: {Damage State: (median_IM, beta)}
damage_states = {
    'Minor Damage Level': (20.0, 40.0),# Median DI=20, β=40
    'Moderate Damage Level': (40.0, 40.0),
    'Severe Damage Level': (60.0, 50.0),
    'Failure Level': (100.0, 50.0)
}

# Generate intensity measure (IM) values from 0.0 to 100.0
im_values = median_di # Structural Ductility Damage Index
X_LABEL = 'Structural Ductility Damage Index (%)  [IM]'
S07.FRAGILITY_CURVE(im_values, damage_states, X_LABEL, SEMILOGY=False, PLOT_KIND=True)

#%%------------------------------------------------------------------------------------------------
import matplotlib.pyplot as plt

plt.figure(figsize=(12,8))
plt.plot(disp_05, base_reaction,color='black')
plt.xlabel("Displacement [m]")
plt.ylabel("Base Reaction [N]")
plt.title("Displacement & Base Reaction Relation From Last Dynamic Analysis")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
#%%------------------------------------------------------------------------------------------------
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

DISP_Z = MAX_ABS(disp_05)  
VELO_Z = MAX_ABS(velo_05) 
ACCE_Z = MAX_ABS(accel_05) 
BASE_Z = MAX_ABS(base_reaction) 

plt.figure(1, figsize=(8, 6))
plt.plot(time, disp_05, color='blue', linewidth=2)
plt.plot(time, DISP_Z, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Displacement [m]')
plt.title(f'Time vs Displacement - MAX. ABS: {DISP_Z[-1]}')
plt.grid()
plt.show()

plt.figure(2, figsize=(8, 6))
plt.plot(time, velo_05, color='blue', linewidth=2)
plt.plot(time, VELO_Z, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Velocity [m/s]')
plt.title(f'Time vs Velocity - MAX. ABS: {VELO_Z[-1]}')
plt.grid()
plt.show()

plt.figure(3, figsize=(8, 6))
plt.plot(time, accel_05, color='blue', linewidth=2)
plt.plot(time, ACCE_Z, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Acceleration [m/s^2]')
plt.title(f'Time vs Acceleration - MAX. ABS: {ACCE_Z[-1]}')
plt.grid()
plt.show()

plt.figure(4, figsize=(8, 6))
plt.plot(time, base_reaction, color='blue', linewidth=2)
plt.plot(time, BASE_Z, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Base-reaction [N]')
plt.title(f'Time vs Base-reaction - MAX. ABS: {BASE_Z[-1]}')
plt.grid()
plt.show()
#%%------------------------------------------------------------------------------------------------  
# EXPORT DATA TO EXCEL
DATA_TOTAL = {

    'Damping_Ratio_02': max_dr_02,
    'Damping_Ratio_03': max_dr_03,
    'Damping_Ratio_04': max_dr_04,
    'Damping_Ratio_05': max_dr_05,
    'Max_displacement_02': max_disp_02,
    'Max_displacement_03': max_disp_03,
    'Max_displacement_04': max_disp_04,
    'Max_displacement_05': max_disp_05,
    'Max_velocity_02': max_velo_02,
    'Max_velocity_03': max_velo_03,
    'Max_velocity_04': max_velo_04,
    'Max_velocity_05': max_velo_05,
    'Max_acceleration_02': max_accel_02,
    'Max_acceleration_03': max_accel_03,
    'Max_acceleration_04': max_accel_04,
    'Max_acceleration_05': max_accel_05,
    'Max_Base_Reaction': max_base_reaction,
    'Ductility_Damage_Index': max_DI,
    'Period': max_T,
}
# Convert to DataFrame
results_df = pd.DataFrame(DATA_TOTAL)
# Export the DataFrame to an Excel file
results_df.to_excel('MDOF_RESPONSE_SPECTRUM_FREE_VIBRATION_DUCT_OSF_LAST_ANALYSIS_DATA_RESULTS.xlsx', index=False)
#%%------------------------------------------------------------------------------------------------  
XLABEL = 'Displacement'
YLABEL = 'Base Reaction'
TITLE = f'{YLABEL} and {XLABEL} scatter chart'
COLOR = 'orange'
X = median_disp
Y = median_reaction
S01.PLOT_SCATTER(X, Y , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

# CLUSTER DATA
S01.CLUSTER_DATA(X, Y, XLABEL, YLABEL, MAX_CLUSTERS=3)
#%%------------------------------------------------------------------------------------------------
XLABEL = 'Velocity'
YLABEL = 'Base Reaction'
TITLE = f'{YLABEL} and {XLABEL} scatter chart'
COLOR = 'cyan'
X = median_velo
Y = median_reaction
S01.PLOT_SCATTER(X, Y , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

# CLUSTER DATA
S01.CLUSTER_DATA(X, Y, XLABEL, YLABEL, MAX_CLUSTERS=3)
#%%------------------------------------------------------------------------------------------------
XLABEL = 'Acceleration'
YLABEL = 'Base Reaction'
TITLE = f'{YLABEL} and {XLABEL} scatter chart'
COLOR = 'lime'
X = median_acce
Y = median_reaction
S01.PLOT_SCATTER(X, Y , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

# CLUSTER DATA
S01.CLUSTER_DATA(X, Y, XLABEL, YLABEL, MAX_CLUSTERS=3)
#%%------------------------------------------------------------------------------------------------
XLABEL = 'Displacement'
YLABEL = 'Structural Ductility Damage Index (%)'
TITLE = f'{YLABEL} and {XLABEL} scatter chart'
COLOR = 'purple'
X = median_disp
Y = median_di
S01.PLOT_SCATTER(X, Y , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

# CLUSTER DATA
S01.CLUSTER_DATA(X, Y, XLABEL, YLABEL, MAX_CLUSTERS=3)
#%%------------------------------------------------------------------------------------------------
# PLOT THE TIME-HISTORY
S01.PLOT_TIME_HISTORY(time, disp_05, velo_05, accel_05, base_reaction)
#%%------------------------------------------------------------------------------------------------

# RANDOM FOREST ANALYSIS
"""
This code predicts the free-vibration safety of a structure using simulation data by training a Random Forest Classifier to
 classify whether the system is "safe" or "unsafe" based on features like maximum displacement, velocity, acceleration,
 and base reaction. A regression model is also trained to estimate safety likelihood. It evaluates model performance using
 metrics like classification accuracy, mean squared error, and R² score. Additionally, it identifies key features influencing
 safety through feature importance analysis. The tool aids in free-vibration risk assessment, structural optimization, and understanding
 critical safety parameters.
"""

data = {
    'Max_displacement': median_disp,
    'Max_velocity': median_velo,
    'Max_acceleration': median_acce,
    'Max_Base_Reaction': median_reaction,
}


# Convert to DataFrame
df = pd.DataFrame(data)
#print(df)
S01.RANDOM_FOREST(df)
#%%------------------------------------------------------------------------------------------------
# PLOT HEATMAP FOR CORRELATION 
S01.PLOT_HEATMAP(df)
#%%------------------------------------------------------------------------------------------------
# MULTIPLE REGRESSION MODEL
S01.MULTIPLE_REGRESSION(df) 
"""
#%%------------------------------------------------------------------------------------------------
# MACHINE LEARNING: LONG SHORT-TREM MEMERY (LSTM) METHOD
x = max_displacement 
y = max_acceleration 
Demand_X = x[-1]
look_back = 500#int(NUM_SIM * 0.5)
ITERATION = 200
XLABEL = 'Max Displacement'
YLABEL = 'Max Acceleration'
#S01.PREDICT_LSTM(x, y, Demand_X, look_back, ITERATION, XLABEL, YLABEL)
#%%------------------------------------------------------------------------------------------------
# PERFORM RELIABILITY ANALYSIS FOR BASE REACTION AND ELEMENT CAPACITY
mean_capacity = np.mean(FU)    # Mean Element Ultimate Capacity
std_dev_capacity = np.std(FU)  # Std Element Ultimate Capacity
num_sim = NUM_SIM
S01.RELIABILITY_ANALYSIS(max_base_reaction, num_sim, mean_capacity, std_dev_capacity)
#%%------------------------------------------------------------------------------------------------
# NEURAL NETWORK FOR FAILURE PROBABILIYY ESTIMATION
X1 = mean_capacity
X2 = max_base_reaction
S01.NEURAL_NETWORK_FAILURE_PROBABILIYY_ESTIMATION(max_base_reaction, X2, NUM_SEISMIC)
#%%------------------------------------------------------------------------------------------------
# MARKOV CHAIN MODEl (structural damage analysis by evaluating Structural Ductility Damage Index)
FILE_TF = False         # Indicate whether to read data from a file or use provided data
file_path = None        # Not used when 'file_tf' is False
DATA = max_DI # If not using a file, replace None with a NumPy array of data

S03.MARKOV_CHAIN(FILE_TF, file_path, DATA)
#------------------------------------------------------------------------------------------------
"""