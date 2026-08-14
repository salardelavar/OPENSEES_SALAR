###########################################################################################################
#                    >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                     #
#     FRAGILITY ANALYSIS BASED ON ACCELERATION, STRUCTURAL DUCTILITY DAMAGE INDEX, ENERGY DISSIPATION     #
#            CAPACITY INDEX, EQULIVALENT VISCOUS DAMPING RATIO WITH NONLINEAR DYNAMIC ANALYSIS            #
#          OF A SINGLE-DEGREE-OF-FREEDOM (SDOF) SYSTEM UTILIZING 30 GROUND MOTIONS IN OPENSEES            #
#---------------------------------------------------------------------------------------------------------#
#                                   INELASTIC RESPONSE SPECTRUM ANALYSIS                                  #
#---------------------------------------------------------------------------------------------------------#
#                EQUIVALENT VISCOUS DAMPING RATIO: xi_eq = 100 * E_d / (4 * pi * E_s)                     #
#---------------------------------------------------------------------------------------------------------#
#          ENERGY DISSIPATION CAPACITY INDEX = 100 * E_d(earthquake) / E_d(cyclic displacement)           #
#---------------------------------------------------------------------------------------------------------#
# This program performs  Dynamic Analysis on a Single-Degree-of-Freedom (SDOF) system                     #
# subjected to 30 seismic ground motions. The analysis evaluates the structural response under varying    #
# levels of seismic intensity.                                                                            #
# The framework is designed to support researchers and engineers in assessing the probabilistic seismic   #
# performance of structures, with a focus on understanding the impact of uncertainty on structural        #
# response and design.                                                                                    #
#---------------------------------------------------------------------------------------------------------#
# Key Features:                                                                                           #
# - Simulation of SDOF system using OpenSees.                                                             #
# - Incremental scaling of ground motions for Nonlinear Dynamic Analysis.                                 #
# - Probabilistic fragility assessment based on predefined damage states.                                 #
# - Visualization of structural response and fragility curves.                                            #
# - Export of results for further analysis.                                                               #
#---------------------------------------------------------------------------------------------------------#
#                        THIS PYTHON SCRIPT WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                 #
#                                   EMAIL: salar.d.ghashghaei@gmail.com                                   #
########################################################################################################### 

"""
This code implements a comprehensive nonlinear dynamic analysis framework for
performance-based earthquake engineering assessment of single-degree-of-freedom
(SDOF) systems. The methodology combines traditional nonlinear time-history
analysis with modern probabilistic and machine learning techniques for advanced
structural performance evaluation.

KEY ENGINEERING OBJECTIVES:
1. Comparative assessment of hysteretic models for seismic response prediction
2. Probabilistic seismic demand analysis using multiple ground motions
3. Development of fragility curves for performance-based earthquake engineering
4. Integration of data science methods for structural reliability assessment

ANALYTICAL FEATURES:
- Nonlinear material behavior with pinching and degradation
- Response spectrum analysis across period range
- Real-time structural health monitoring metrics
- Statistical characterization of seismic demands
- Machine learning-based damage prediction
---------------------------------
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
 - Peak base shear (reaction): decays faster in HYSTERETIC due to strength loss; BILINEAR sustains higher peak restoring forces for the same drift until hardening or limiting criteria apply.

Visualization:
 - Plot time histories (disp, vel, acc), hysteresis loops (force vs disp), and envelope curves to compare models.
 - Response spectra for displacement, velocity, and acceleration can be constructed from peak responses across parameter sweeps (e.g., varying T or post-yield stiffness).

Implications for seismic assessment:
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

------------------------------------------------ 
Energy Dissipation Capacity Index (EDCI):    
The Energy Dissipation Capacity Index is a quantitative
 measure used in structural engineering to evaluate how
 effectively a structural element (e.g., a beam, column, shear wall, or connection)
 can absorb and dissipate energy during seismic loading compared to its
 performance under controlled cyclic displacement loading.
It compares the actual energy absorbed during an earthquake
 with the maximum energy dissipation capacity that the component
 demonstrates in a laboratory‑style cyclic test.  
 
Why This Index Is Important:
During an earthquake, structures undergo repeated cycles of deformation.
 A system with high energy dissipation capacity can withstand more damage
 without collapsing because it can convert seismic input energy into hysteretic energy, not elastic rebound.

The EDCI helps engineers understand:
[1] Ductility performance
[2] Hysteretic behavior
[3] Damage tolerance
[4] Collapse prevention capability
It is especially used in performance‑based seismic evaluation and retrofit design.   
------------------------------------------------
Equivalent Viscous Damping Ratio (EVDR):
This function computes (EQULIVALENT_VISCOUS_DAMPING_RATIO_FUN)
the equivalent viscous damping ratio xi_eq = 100 * E_d / (4 * pi * E_s) from a hysteresis loop,
where E_s = 0.5 * F_max * u_max uses the secant stiffness at maximum displacement.
Method 2 (Shoelace) gives the exact enclosed area E_d, preferred for real test data; Method 1 (convex hull) is a convex upper bound.
Output xi_eq can be added to inherent damping for capacity-spectrum seismic assessment; plot shows loop, filled area, and energy values.
Ensure the input is a complete, ordered hysteresis cycle, and use Method 2 to capture pinching or degradation accurately.       
"""
#%%------------------------------------------------------
import openseespy.opensees as ops
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time as ti
import SALAR_MATH as S01
import ANALYSIS_FUNCTION as S02
import PERIOD_FUN as S03
import DAMPING_RATIO_FUN as S04
import EIGENVALUE_ANALYSIS_FUN as S05
import RAYLEIGH_DAMPING_FUN as S06
import BILINEAR_CURVE as S07
import OPENSEEES_HYSTERETICSM_FORCE_DISP_FUN as S08
import SALAR_MATH as S09
import MARKOV_CHAIN as S0333
import FRAGILITY_CURVE_FUN as S077
from scipy.stats import norm
import EQULIVALENT_VISCOUS_DAMPING_RATIO_FUN as S055
import DAMAGE_INDEX_FUN as S066
#%%------------------------------------------------------
# Define parameters (units: m, N)
NUM_T = 30                                      # Total number for Structural Period in each simulation
NUM_SEISMIC = 30                                # Total number for seismic simulation


TOTAL_MASS = 1500.0                            # [kg] Total Mass of Structure
GMfact = 9.81                                  # [m/s^2] standard acceleration of gravity
#%%------------------------------------------------------
# Define Analysis Properties
MAX_ITERATIONS = 1000000   # Convergence iteration for test
MAX_TOLERANCE = 1.0e-6    # Convergence tolerance for test
#%%------------------------------------------------------
# EVALUATION OF DISSIPATED ENERGY CAPACITY INDEX
def DISSIPATED_ENERGY_FUN_WITH_PLOT(displacement, base_shear, method, title="Hysteresis Curve"):
    if method == 1:
        """
        Compute dissipated energy using convex hull and plot the hysteresis curve
        with the outer hull area shaded.
    
        Parameters
        ----------
        displacement : array-like
        base_shear  : array-like
        title       : str
    
        Returns
        -------
        float
            Area of convex hull (dissipated energy)
        """
        import numpy as np
        from scipy.spatial import ConvexHull
        import matplotlib.pyplot as plt
        displacement = np.asarray(displacement)
        base_shear  = np.asarray(base_shear)
    
        if displacement.size != base_shear.size:
            raise ValueError("Displacement and base shear arrays must have equal lengths.")
    
        points = np.column_stack((displacement, base_shear))
        hull = ConvexHull(points)
        area = hull.volume   # 2D hull → area
    
        fig, ax = plt.subplots(figsize=(7, 6))
    
        # Plot full hysteresis
        ax.plot(displacement, base_shear, 'k-', linewidth=1, label="Hysteresis Curve")
    
        # Plot convex hull edges
        hull_pts = points[hull.vertices]
        ax.plot(hull_pts[:, 0], hull_pts[:, 1], 'r--', lw=2, label="Convex Hull")
    
        # Shade hull area
        ax.fill(hull_pts[:, 0], hull_pts[:, 1], color='red', alpha=0.25, label="Hull Area (Energy)")
    
        # Labels and style
        ax.set_title(f"{title} - (Convex Hull)")
        ax.set_xlabel("Displacement (m)")
        ax.set_ylabel("Base Shear (N)")
        ax.grid(True, linestyle='--', alpha=0.5)
        ax.legend()
        
    if method == 2:  
        import numpy as np
        import matplotlib.pyplot as plt
        # Data preparation
        disp = np.asarray(displacement, dtype=float)
        shear = np.asarray(base_shear, dtype=float)

        if disp.size != shear.size:
            raise ValueError("Displacement and base shear arrays must have the same length.")
        if disp.size < 3:
            raise ValueError("At least 3 points are required to form a closed loop.")

        # Close the loop if not already closed (important for shoelace)
        if not (disp[0] == disp[-1] and shear[0] == shear[-1]):
            disp = np.append(disp, disp[0])
            shear = np.append(shear, shear[0])

        # Dissipated energy (E_d) via Shoelace formula
        x = disp
        y = shear
        area = 0.5 * np.abs(np.dot(x[:-1], y[1:]) - np.dot(y[:-1], x[1:]))
        
        # Plotting
        fig, ax = plt.subplots(figsize=(7, 6))
    
        # Hysteresis curve
        idx_max = np.argmax(np.abs(disp))
        ax.plot(disp, shear, 'k-', linewidth=1.2, label="Hysteresis Loop")
        ax.scatter(disp[idx_max], shear[idx_max], color='blue', s=80,
                   zorder=5, label=r"$(u_{\rm max}, F_{\rm max})$")
    
        # Fill the enclosed area
        ax.fill(disp, shear, color='red', alpha=0.25, label=f"E$_d$ = {area:.3f} N·m")
    
    
        ax.set_title(title)
        ax.set_xlabel("Displacement (m)")
        ax.set_ylabel("Base Shear (N)")
        ax.grid(True, linestyle='--', alpha=0.5)
        ax.legend(loc='lower right')
        fig.tight_layout()
    
    return area, fig
#%%------------------------------------------------------
### OPENSEES FUNCTION
def SDOF(MAT_TYPE, TOTAL_MASS, ANAL_TYPE, i):
    # Initialize model
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)
    
    MAX_ITERATIONS = 5000   # Maximum number of iterations
    MAX_TOLERANCE = 1.0e-6  # Specified tolerance for convergence
    GMfact = 9.810          # [m/s^2] standard acceleration of gravity or standard acceleration 
        
    # Define nodes
    ops.node(1, 0.0)  # Fixed base
    ops.node(2, 0.0)  # Mass node
        
    # Define boundary conditions
    ops.fix(1, 1)

    # Define mass
    ops.mass(2, TOTAL_MASS)
    
    FY = 85000.0                                     # [N] Yield Force of Structure
    FU = 1.18 * FY                                   # [N] Ultimate Force of Structure
    Ke = 4500000.0                                   # [N/m] Spring Elastic Stiffness
    DY = FY / Ke                                     # [m] Yield Displacement
    DSU = 0.36                                       # [m] Ultimate Displacement
    Ksh = (FU - FY) / (DSU - DY)                     # [N/m] Displacement Hardening Modulus
    Kp = FU / DSU                                    # [N/m] Spring Plastic Stiffness
    b = Ksh / Ke                                     # Displacement Hardening Ratio
    DR = 0.03                                        # Damping Ratio
    # Positive branch points
    pos_disp = [0, DY, DSU, 1.1*DSU, 1.25*DSU]
    pos_force = [0, FY, FU, 0.2*FU, 0.1*FU]
    KP = np.array([FY, DY, FU, DSU, 0.2*FU, 1.1*DSU, 0.1*FU, 1.25*DSU])
    
    # Negative branch points
    neg_disp = [0, -DY, -DSU, -1.1*DSU, -1.25*DSU]
    neg_force = [0, -FY, -FU, -0.2*FU, -0.1*FU]
    KN = np.array([-FY, -DY, -FU, -DSU, -0.2*FU, -1.1*DSU, -0.1*FU, -1.25*DSU])
    
    # Plot
    plt.figure(0, figsize=(12, 10))
    plt.plot(pos_disp, pos_force, marker='o', color='red')
    plt.plot(neg_disp, neg_force, marker='o', color='black')
    
    plt.xlabel("Displacement [m]")
    plt.ylabel("Force [N]")
    plt.title("Force–Displacement Curve")
    plt.grid(True)
    plt.axhline(0, linewidth=0.5)
    plt.axvline(0, linewidth=0.5)
    plt.show()
    
    # Define material properties
    MAT_TAG = 1 # SPRING TAG
    # FORCE-DISPLACEMENT RELATIONSHIP OF LATERAL SPRING AND PLOT 
    DP = [0, 0, 0, 0]
    FP = [0, 0, 0, 0]
    DN = [0, 0, 0, 0]
    FN = [0, 0, 0, 0]
    #DSU = DY * duct # IN EACH STEP IT WILL BE CHNAGED
    #FU = FY * osf   # IN EACH STEP IT WILL BE CHNAGED
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
        S08.OPENSEEES_HYSTERETICSM_FORCE_DISP_FUN(MAT_TAG, DP, FP, DN, FN, PLOT = False, X_LABEL='Displacement (mm)', Y_LABEL='Force [N]', TITLE='FORCE-DISPLACEMENT CURVE')
    if MAT_TYPE == 'ELASTIC':
        #ops.uniaxialMaterial('Elastic', MAT_TAG, Ke)             # TESNSION AND COMPRESSION IS SAME VALUES
        ops.uniaxialMaterial('Elastic', MAT_TAG, Ke ,0.0, 0.5*Ke) # TESNSION AND COMPRESSION IS NOT SAME VALUES
        # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/ElasticUni.html
    
    MatTag_C = 2 # SPRING DAMPER
    alpha = 1.0    # velocity exponent (usually 0.3–1.0)
    omega = np.sqrt(Ke/TOTAL_MASS)
    Cd = 2 * DR * omega * TOTAL_MASS  # [N·s/m] Damping coefficient 
    ops.uniaxialMaterial('Viscous', 2, Cd, alpha)  # Material for C (alpha=1.0 for linear)
    
    # Define element
    ops.element('zeroLength', 1, 1, 2, '-mat',  MAT_TAG, MatTag_C, '-dir', 1, 1)  # DOF[1] LATERAL SPRING
    
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(2, 1.0)
    
    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.system('BandGeneral')
    
    if MAT_TYPE == 'ELASTIC':
        ops.algorithm('Linear')
    if MAT_TYPE == 'INELASTIC':
        ops.algorithm('Newton')
        
    center_node = 2
    time = []
    disp = []
    velo = []
    acc = []
    reaction = []
    stiffness = []
    OMEGA, PERIOD = [], []
    PERIOD_MIN, PERIOD_MAX = [], []
    DI = []
    
    if ANAL_TYPE == 'PERIOD': 
      #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(1, 0.5*DR, DR, 0, 1)
      PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(1, PLOT=True) 
      # Compute modal properties
      ops.modalProperties("-print", "-file", f"SALAR_ModalReport_{ANAL_TYPE}.txt", "-unorm")
      return PERIODmin, PERIODmax
  
    if ANAL_TYPE == 'STATIC': 
        ops.integrator('LoadControl', 1.0)
        ops.analysis('Static')
        OK = ops.analyze(1)
        S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS)
        ops.reactions()
        reaction.append(ops.nodeReaction(1, 1))                         # BASE REACTION
        disp.append(ops.nodeDisp(center_node, 1))                       # DISPLACEMENT NODE 02 IN X DIR                
        # EVALUATION OF DUCTILITY DAMAGE INDEX
        if MAT_TYPE == 'INELASTIC':
            di = S066.DAMAGE_INDEX_FUN(disp[-1], DY, DSU)
            DI.append(di)                   # DAMAGE INDEX
        if MAT_TYPE == 'ELASTIC':
            DI.append(0.0)    
        print('\n\nSTATIC ANALYSIS DONE.\n\n')  
            
        DATA = (reaction, disp, DI)
        
        return  DATA
    
    if ANAL_TYPE == 'PUSHOVER': # STATIC TIME-HISTORY ANALYSIS
        IDctrlDOF = 1   # 1: Horizental Dispalcement - 2: Vertical Dispalcement
        DINCR = -0.001  # [m] Incremental Vertical Displacement
        DMAX = -DSU     # [m] Max. Displacement
        ops.integrator('DisplacementControl', center_node, IDctrlDOF, DINCR)
        ops.analysis('Static')
        Nsteps =  int(np.abs(DMAX/ DINCR)) 
        STEP = 0.0
        for step in range(Nsteps):
            OK = ops.analyze(1)
            S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS)
            ops.reactions()
            reaction.append(ops.nodeReaction(1, 1))                         # BASE REACTION
            disp.append(ops.nodeDisp(center_node, 1))                      # DISPLACEMENT NODE 02 IN X DIR       
            # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
            #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(1, 0.5*DR, DR, 0, 1)
            PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(1, PLOT=True)
            PERIOD_MIN.append(PERIODmin)
            PERIOD_MAX.append(PERIODmax) 
            # EVALUATION OF DUCTILITY DAMAGE INDEX
            if MAT_TYPE == 'INELASTIC':
                di = S066.DAMAGE_INDEX_FUN(disp[-1], DY, DSU)
                DI.append(di)                   # DAMAGE INDEX
            if MAT_TYPE == 'ELASTIC':
                DI.append(0.0)
            STEP += 1
            #print(STEP, disp[-1], reaction[-1])
            print(f"Step: {STEP}, Displacement: {disp[-1]:.4f} m, Reaction: {reaction[-1]:.2f} N")     
        else:
            print('\n\nPUSHOVER ANALYSIS DONE.\n\n')
            
        DATA = (reaction, disp, DI,
                np.array(PERIOD_MIN), np.array(PERIOD_MAX))
    
        return  DATA
    
    if ANAL_TYPE == 'CYCLIC_DISPLACEMENT': # STATIC TIME-HISTORY ANALYSIS
        IDctrlDOF = 1   # 1: Horizental Dispalcement - 2: Vertical Dispalcement
        DMAX = -DSU     # [m] Max. Displacement
        n_points = 10000
        # 4. CYCLIC DISPALCEMENT PROTOCOL
        # Key strain points (same logic as your protocol)
        key_disp = np.array([
             0.0,
             0.1*DMAX,   -0.1*DMAX,
             0.5*DMAX,   -0.5*DMAX,
             0.8*DMAX,   -0.8*DMAX,
             DMAX,       -DMAX,
             0.1*DMAX,   -0.1*DMAX,
             0.5*DMAX,   -0.5*DMAX,
             0.8*DMAX,   -0.8*DMAX,
             DMAX,       -DMAX,
             0.0
        ])
        
        # Generate 1000-point displacement protocol
        disp_protocol = np.interp(
            np.linspace(0, len(key_disp) - 1, n_points),
            np.arange(len(key_disp)),
            key_disp
        )
        
        ops.analysis('Static')
        STEP = 0.0
        for target_disp in disp_protocol:
            current_disp = ops.nodeDisp(center_node, 1) # DISPALCEMENT APPLIED IN MIDDLE NOD IN X DIR.
            dU = target_disp - current_disp
            ops.integrator('DisplacementControl', center_node, IDctrlDOF, dU)
            OK = ops.analyze(1)
            S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS)
            ops.reactions()
            reaction.append(ops.nodeReaction(1, 1))                         # BASE REACTION
            disp.append(ops.nodeDisp(center_node, 1))                       # DISPLACEMENT NODE 02 IN X DIR       
            # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
            #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(1, 0.5*DR, DR, 0, 1)
            PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(1, PLOT=True)
            PERIOD_MIN.append(PERIODmin)
            PERIOD_MAX.append(PERIODmax)
            # EVALUATION OF DUCTILITY DAMAGE INDEX
            if MAT_TYPE == 'INELASTIC':
                di = S066.DAMAGE_INDEX_FUN(disp[-1], DY, DSU)
                DI.append(di)                   # DAMAGE INDEX
            if MAT_TYPE == 'ELASTIC':
                DI.append(0.0)
            STEP += 1
            #print(STEP, disp[-1], reaction[-1])
            print(f"Step: {STEP}, Displacement: {disp[-1]:.4f} m, Reaction: {reaction[-1]:.2f} N")     
        else:
            print('\n\nCYCLIC DISPLAEMENT ANALYSIS DONE.\n\n')
            
        DATA = (reaction, disp, DI,
                np.array(PERIOD_MIN), np.array(PERIOD_MAX))
    
        return  DATA 

    if ANAL_TYPE == 'STATIC_EXTERNAL_TIME-DEPENDENT_LOADING': # STATIC TIME-HISTORY ANALYSIS
        #%% DEFINE EXTERNAL TIME-DEPENDENT LOADING PROPERTIES
        # IN HERE ANALYSIS TIME AND DURATION ARE LOAD STEPS
        duration = 20.0             # [s] Analysis duration
        dt = 0.01                   # [s] Time step
        DT = dt                     # [s] Time step
        DT_time = 5.0               # [s] Total external Load Analysis Durations [*******]
        force_amplitude = 960.0     # [N] Amplitude Force
        omega_DT = 5.0715           # [rad/s] Natural angular frequency
        time_steps = int(duration/dt)
        
        # Check function
        def CHECK_FUN(DT_time ,duration):
            if DT_time > duration:
                print('\n\nAnalysis Duration Must be greater than External Load Duration\n\n')        
                exit()
            return -1

        CHECK_FUN(DT_time ,duration)
        def EXTERNAL_TIME_DEPENDENT(force_amplitude, omega_DT, DT, DT_time): # P(t) = P0 sin(wt)
            import numpy as np
            import matplotlib.pyplot as plt
            # External Load Durations
            num_steps = int(DT_time / DT)
            load_time = np.linspace(0, DT_time, num_steps) 
            target_frequency = 1.0 * omega_DT  # Target excitation frequency
            DT_load = force_amplitude * np.sin(target_frequency * load_time)
            # Plot External Time-dependent Loading
            plt.figure(figsize=(10, 6))
            plt.plot(load_time, DT_load, label=f'External Loading - Max: {np.max(DT_load):.3f}', linewidth=5)
            plt.title('External Time-dependent Loading Over Time')
            plt.xlabel('Time (s)')
            plt.ylabel('Force (N)')
            plt.grid(True)
            plt.legend()
            plt.show()
            return DT_load

        #DT_load = EXTERNAL_TIME_DEPENDENT(force_amplitude, omega_DT, DT, DT_time)

        def EXTERNAL_TIME_DEPENDENT_02(force_amplitude, omega_DT, DT, DT_time): #  P(t) = P0 exp(-0.05wt) sin(wt) 
            import numpy as np
            import matplotlib.pyplot as plt
            # External Load Durations
            num_steps = int(DT_time / DT)
            load_time = np.linspace(0, DT_time, num_steps) 
            target_frequency = 1.0 * omega_DT  # Target excitation frequency
            DT_load = force_amplitude * np.exp(-0.05*target_frequency * load_time) * np.sin(target_frequency * load_time)
            # Plot External Time-dependent Loading
            plt.figure(figsize=(10, 6))
            plt.plot(load_time, DT_load, label=f'External Loading - Max: {np.max(DT_load):.3f}', linewidth=5)
            plt.title('External Time-dependent Loading Over Time')
            plt.xlabel('Time (s)')
            plt.ylabel('Force (N)')
            plt.grid(True)
            plt.legend()
            plt.show()
            #print(load_time, DT_load)
            return DT_load

        DT_load02 = EXTERNAL_TIME_DEPENDENT_02(force_amplitude, omega_DT, DT, DT_time)
        # Static Time-depenent External loading analysis
        TS_TAG = 3
        PATT_TAG = 3
        # Apply time-dependent explosion loading
        ops.timeSeries('Path', TS_TAG, '-dt', dt, '-values', *DT_load02)
        ops.pattern('Plain', TS_TAG, PATT_TAG)
        ops.load(center_node, 1.0)
        
        ops.system('BandGeneral')
        ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/normDispIncr.html
        ops.algorithm('Newton')  # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/algorithm.html
        ops.integrator('LoadControl', dt)
        #ops.integrator('DisplacementControl', center_node, 0.001)
        ops.analysis('Static') # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/analysis.html
        
        STEP = 0.0
        stable = 0
        
        for JJ in range(time_steps):
            stable = ops.analyze(1)
            S02.ANALYSIS(stable, 1, MAX_TOLERANCE, MAX_ITERATIONS)
            ops.reactions()
            reaction.append(ops.nodeReaction(1, 1))                         # BASE REACTION
            disp.append(ops.nodeDisp(center_node, 1))                       # DISPLACEMENT NODE 02 IN X DIR      
            # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
            #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(1, 0.5*DR, DR, 0, 1)
            PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(1, PLOT=True)
            PERIOD_MIN.append(PERIODmin)
            PERIOD_MAX.append(PERIODmax)
            # EVALUATION OF DUCTILITY DAMAGE INDEX
            if MAT_TYPE == 'INELASTIC':
                di = S066.DAMAGE_INDEX_FUN(disp[-1], DY, DSU)
                DI.append(di)                   # DAMAGE INDEX
            if MAT_TYPE == 'ELASTIC':
                DI.append(0.0)
            STEP += 1
            #print(STEP, disp[-1], reaction[-1])
            print(f"Step: {STEP}, Displacement: {disp[-1]:.4f} m, Reaction: {reaction[-1]:.2f} N")     
        else:
            print('\n\nSTATIC EXTERNAL TIME-DEPENDENT LOADING ANALYSIS DONE.\n\n')
            
        DATA = (reaction, disp, DI,
                np.array(PERIOD_MIN), np.array(PERIOD_MAX))
    
        return  DATA
    
    if ANAL_TYPE == 'DYNAMIC_EXTERNAL_TIME-DEPENDENT_LOADING': # DYNAMIC TIME-HISTORY ANALYSIS
        #%% DEFINE EXTERNAL TIME-DEPENDENT LOADING PROPERTIES
        duration = 20.0             # [s] Analysis duration
        dt = 0.01                   # [s] Time step
        DT = dt                     # [s] Time step
        DT_time = 5.0               # [s] Total external Load Analysis Durations [*******]
        force_amplitude = 960.0     # [N] Amplitude Force
        omega_DT = 5.0715           # [rad/s] Natural angular frequency
        DR = 0.03                   # Damping Ratio

        # Check function
        def CHECK_FUN(DT_time ,duration):
            if DT_time > duration:
                print('\n\nAnalysis Duration Must be greater than External Load Duration\n\n')        
                exit()
            return -1

        CHECK_FUN(DT_time ,duration)
        def EXTERNAL_TIME_DEPENDENT(force_amplitude, omega_DT, DT, DT_time): # P(t) = P0 sin(wt)
            import numpy as np
            import matplotlib.pyplot as plt
            # External Load Durations
            num_steps = int(DT_time / DT)
            load_time = np.linspace(0, DT_time, num_steps) 
            target_frequency = 1.0 * omega_DT  # Target excitation frequency
            DT_load = force_amplitude * np.sin(target_frequency * load_time)
            # Plot External Time-dependent Loading
            plt.figure(figsize=(10, 6))
            plt.plot(load_time, DT_load, label=f'External Loading - Max: {np.max(DT_load):.3f}', linewidth=5)
            plt.title('External Time-dependent Loading Over Time')
            plt.xlabel('Time (s)')
            plt.ylabel('Force (N)')
            plt.grid(True)
            plt.legend()
            plt.show()
            return DT_load

        #DT_load = EXTERNAL_TIME_DEPENDENT(force_amplitude, omega_DT, DT, DT_time)

        def EXTERNAL_TIME_DEPENDENT_02(force_amplitude, omega_DT, DT, DT_time): #  P(t) = P0 exp(-0.05wt) sin(wt) 
            import numpy as np
            import matplotlib.pyplot as plt
            # External Load Durations
            num_steps = int(DT_time / DT)
            load_time = np.linspace(0, DT_time, num_steps) 
            target_frequency = 1.0 * omega_DT  # Target excitation frequency
            DT_load = force_amplitude * np.exp(-0.05*target_frequency * load_time) * np.sin(target_frequency * load_time)
            # Plot External Time-dependent Loading
            plt.figure(figsize=(10, 6))
            plt.plot(load_time, DT_load, label=f'External Loading - Max: {np.max(DT_load):.3f}', linewidth=5)
            plt.title('External Time-dependent Loading Over Time')
            plt.xlabel('Time (s)')
            plt.ylabel('Force (N)')
            plt.grid(True)
            plt.legend()
            plt.show()
            #print(load_time, DT_load)
            return DT_load

        DT_load02 = EXTERNAL_TIME_DEPENDENT_02(force_amplitude, omega_DT, DT, DT_time)
        # Static Time-depenent External loading analysis
        TS_TAG = 3
        PATT_TAG = 3
        # Apply time-dependent explosion loading
        ops.timeSeries('Path', TS_TAG, '-dt', dt, '-values', *DT_load02)
        ops.pattern('Plain', TS_TAG, PATT_TAG)
        ops.load(center_node, 1.0)
        
        ops.constraints('Plain')
        ops.numberer('Plain')
        ops.system('BandGeneral')
        ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/normDispIncr.html
        #ops.integrator('CentralDifference')  # JUST FOR LINEAR ANALYSIS - INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/centralDifference.html
        alpha=0.5; beta=0.25;
        ops.integrator('Newmark', alpha, beta) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/newmark.html
        #alpha=2/3;gamma=1.5-alpha; gamma=1.5-alpha;beta=(2-alpha)**2/4;
        #ops.integrator('HHT', alpha, gamma, beta) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/hht.html
        ops.algorithm('Newton')  # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/algorithm.html
        ops.analysis('Transient') # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/analysis.html
        
        stable = 0
        current_time = 0.0
        
        while stable == 0 and current_time < duration:
            ops.analyze(1, dt)
            S02.ANALYSIS(stable, 1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS
            current_time = ops.getTime()
            time.append(current_time)
            ops.reactions()
            reaction.append(ops.nodeReaction(1, 1))               # BASE REACTION
            disp.append(ops.nodeDisp(center_node, 1))             # DISPLACEMENT NODE 02 IN X DIR 
            velo.append(ops.nodeVel(center_node, 1))              # VELOCITY NODE 02
            acc.append(ops.nodeAccel(center_node, 1))             # ACCELERATION NODE 02
            stiffness.append(np.abs(reaction[-1] / disp[-1]))
            OMEGA.append(np.sqrt(stiffness[-1]/TOTAL_MASS))
            PERIOD.append((np.pi * 2) / OMEGA[-1])   
            # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
            #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(1, 0.5*DR, DR, 0, 1)
            PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(1, PLOT=True)
            PERIOD_MIN.append(PERIODmin)
            PERIOD_MAX.append(PERIODmax)
            # EVALUATION OF DUCTILITY DAMAGE INDEX
            if MAT_TYPE == 'INELASTIC':
                di = S066.DAMAGE_INDEX_FUN(disp[-1], DY, DSU)
                DI.append(di)                   # DAMAGE INDEX
            if MAT_TYPE == 'ELASTIC':
                DI.append(0.0)
            #print(time[-1], disp[-1], velo[-1])
            print(f"Time: {time[-1]:.4f}, Displacement: {disp[-1]:.4f} m, Reaction: {reaction[-1]:.2f} N")      
        else:
            print('\n\nDYNAMIC EXTERNAL TIME-DEPENDENT LOADING ANALYSIS DONE.\n\n')  
        # Calculating Damping Ratio and Period Using Logarithmic Decrement Analysis 
        damping_ratio = S04.DAMPING_RATIO(disp)  
        
        # Compute modal properties
        ops.modalProperties("-print", "-file", f"SALAR_ModalReport_{ANAL_TYPE}.txt", "-unorm") 
        
        DATA = (time, reaction, disp, velo, acc, DI,
                stiffness, PERIOD, damping_ratio,
                np.array(PERIOD_MIN), np.array(PERIOD_MAX))
        
        return  DATA
        
    if ANAL_TYPE == 'FREE-VIBRATION': # DYNAMIC TIME-HISTORY ANALYSIS
        #%% DEFINE PARAMETERS FOR FREE-VIBRATION ANALYSIS
        u0 = -0.10                         # [m] Initial displacement
        v0 = 0.15                          # [m/s] Initial velocity
        a0 = 0.065                         # [m/s^2] Initial acceleration
        IU = True                          # Free Vibration with Initial Displacement
        IV = True                          # Free Vibration with Initial Velocity
        IA = True                          # Free Vibration with Initial Acceleration
        duration = 20.0                    # [s] Analysis duration
        dt = 0.001                         # [s] Time step
        DR = 0.03                          # Damping Ratio
        
        if IU == True:
            # Define initial displacment
            ops.setNodeDisp(center_node, 1, u0, '-commit') # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/setNodeDisp.html
        if IV == True:
            # Define initial velocity
            ops.setNodeVel(center_node, 1, v0, '-commit')  # INFO LINK: https://openseespydoc.readthedocs.io/en/stable/src/setNodeVel.html
        if IA == True:
            # Define initial  acceleration
            ops.setNodeAccel(center_node, 1, a0, '-commit') # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/setNodeAccel.html
            
        ops.constraints('Plain')
        ops.numberer('Plain')
        ops.system('BandGeneral')
        ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/normDispIncr.html
        #ops.integrator('CentralDifference')  # JUST FOR LINEAR ANALYSIS - INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/centralDifference.html
        alpha=0.5; beta=0.25;
        ops.integrator('Newmark', alpha, beta) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/newmark.html
        #alpha=2/3;gamma=1.5-alpha; gamma=1.5-alpha;beta=(2-alpha)**2/4;
        #ops.integrator('HHT', alpha, gamma, beta) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/hht.html
        ops.algorithm('Newton')  # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/algorithm.html
        ops.analysis('Transient') # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/analysis.html
        
        stable = 0
        current_time = 0.0
        while stable == 0 and current_time < duration:
            ops.analyze(1, dt)
            S02.ANALYSIS(stable, 1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS
            current_time = ops.getTime()
            time.append(current_time)
            ops.reactions()
            reaction.append(ops.nodeReaction(1, 1))               # BASE REACTION
            disp.append(ops.nodeDisp(center_node, 1))             # DISPLACEMENT NODE 02 IN X DIR 
            velo.append(ops.nodeVel(center_node, 1))              # VELOCITY NODE 02
            acc.append(ops.nodeAccel(center_node, 1))             # ACCELERATION NODE 02
            stiffness.append(np.abs(reaction[-1] / disp[-1]))
            OMEGA.append(np.sqrt(stiffness[-1]/TOTAL_MASS))
            PERIOD.append((np.pi * 2) / OMEGA[-1]) 
            # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
            #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(1, 0.5*DR, DR, 0, 1)
            PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(1, PLOT=True)
            PERIOD_MIN.append(PERIODmin)
            PERIOD_MAX.append(PERIODmax)
            # EVALUATION OF DUCTILITY DAMAGE INDEX
            if MAT_TYPE == 'INELASTIC':
                di = S066.DAMAGE_INDEX_FUN(disp[-1], DY, DSU)
                DI.append(di)                   # DAMAGE INDEX
            if MAT_TYPE == 'ELASTIC':
                DI.append(0.0)
            #print(time[-1], disp[-1], velo[-1])
            print(f"Time: {time[-1]:.4f}, Displacement: {disp[-1]:.4f} m, Reaction: {reaction[-1]:.2f} N")      
        else:
            print('\n\nFREE-VIBRATION ANALYSIS DONE.\n\n')    
        # Calculating Damping Ratio and Period Using Logarithmic Decrement Analysis 
        damping_ratio = S04.DAMPING_RATIO(disp)  
        
        # Compute modal properties
        ops.modalProperties("-print", "-file", f"SALAR_ModalReport_{ANAL_TYPE}.txt", "-unorm") 
        
        DATA = (time, reaction, disp, velo, acc, DI,
                stiffness, PERIOD, damping_ratio,
                np.array(PERIOD_MIN), np.array(PERIOD_MAX))
        
        return DATA
    if ANAL_TYPE == 'SEISMIC': # DYNAMIC TIME-HISTORY ANALYSIS
        #%% DEFINE PARAMETERS FOR SEISMIC ANALYSIS
        duration = 20.0                    # [s] Analysis duration
        dt = 0.01                          # [s] Time step
        GMfact = 9.810                     # [m/s²]standard acceleration of gravity or standard acceleration
        SSF_X = 50.0                       # Seismic Acceleration Scale Factor in X Direction
        SSF_Y = 50.0                       # Seismic Acceleration Scale Factor in Y Direction
        iv0_X = 0.0005                     # [m/s] Initial velocity applied to the node  in X Direction
        iv0_Y = 0.0005                     # [m/s] Initial velocity applied to the node  in Y Direction
        st_iv0 = 0.0                       # [s] Initial velocity applied starting time
        SEI = 'X'                          # Seismic Direction
        DR = 0.03                          # Damping ratio
        
        # Define time series for input motion (Acceleration time history)
        if SEI == 'X':
            SEISMIC_TAG_01 = 100
            gm_accels = np.loadtxt(f'Ground_Acceleration_{i+1}.txt')  # Assumes acceleration in m/s²
            ops.timeSeries('Path', SEISMIC_TAG_01, '-dt', dt, '-values', *gm_accels.tolist(), '-factor', GMfact) # SEISMIC-X
            # Define load patterns
            # pattern UniformExcitation $patternTag $dof -accel $tsTag <-vel0 $vel0> <-fact $cFact>
            ops.pattern('UniformExcitation', SEISMIC_TAG_01, 1, '-accel', SEISMIC_TAG_01, '-vel0', iv0_X, '-fact', SSF_X) # SEISMIC-X
        if SEI == 'Y':
            SEISMIC_TAG_02 = 200
            gm_accels = np.loadtxt(f'Ground_Acceleration_{i+1}.txt')  # Assumes acceleration in m/s²
            ops.timeSeries('Path', SEISMIC_TAG_02, '-dt', dt, '-values', *gm_accels.tolist(), '-factor', GMfact, '-vel0', iv0_Y) # SEISMIC-Y
            # Define load patterns
            # pattern UniformExcitation $patternTag $dof -accel $tsTag <-vel0 $vel0> <-fact $cFact>
            ops.pattern('UniformExcitation', SEISMIC_TAG_02, 2, '-accel', SEISMIC_TAG_02, '-vel0', iv0_Y, '-fact', SSF_Y) # SEISMIC-Y
        if SEI == 'XY':
            SEISMIC_TAG_01 = 100
            ops.timeSeries('Path', SEISMIC_TAG_01, '-dt', dt, '-filePath', 'Ground_Acceleration_X.txt', '-factor', GMfact, '-startTime', st_iv0) # SEISMIC-X
            # Define load patterns
            # pattern UniformExcitation $patternTag $dof -accel $tsTag <-vel0 $vel0> <-fact $cFact>
            ops.pattern('UniformExcitation', SEISMIC_TAG_01, 1, '-accel', SEISMIC_TAG_01, '-vel0', iv0_X, '-fact', SSF_X) # SEISMIC-X 
            SEISMIC_TAG_02 = 200
            ops.timeSeries('Path', SEISMIC_TAG_02, '-dt', dt, '-filePath', 'Ground_Acceleration_Y.txt', '-factor', GMfact) # SEISMIC-Z
            ops.pattern('UniformExcitation', SEISMIC_TAG_02, 2, '-accel', SEISMIC_TAG_02, '-vel0', iv0_Y, '-fact', SSF_Y)  # SEISMIC-Z
        print('Seismic Defined Done.')
        
        ops.constraints('Plain')
        ops.numberer('Plain')
        ops.system('BandGeneral')
        ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/normDispIncr.html
        #ops.integrator('CentralDifference')  # JUST FOR LINEAR ANALYSIS - INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/centralDifference.html
        alpha=0.5; beta=0.25;
        ops.integrator('Newmark', alpha, beta) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/newmark.html
        #alpha=2/3;gamma=1.5-alpha; gamma=1.5-alpha;beta=(2-alpha)**2/4;
        #ops.integrator('HHT', alpha, gamma, beta) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/hht.html
        ops.algorithm('Newton')  # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/algorithm.html
        ops.analysis('Transient') # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/analysis.html
        
        stable = 0
        current_time = 0.0
        while stable == 0 and current_time < duration:
            ops.analyze(1, dt)
            S02.ANALYSIS(stable, 1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS
            current_time = ops.getTime()
            time.append(current_time)
            ops.reactions()
            reaction.append(ops.nodeReaction(1, 1))               # BASE REACTION
            disp.append(ops.nodeDisp(center_node, 1))             # DISPLACEMENT NODE 02 IN X DIR 
            velo.append(ops.nodeVel(center_node, 1))              # VELOCITY NODE 02
            acc.append(ops.nodeAccel(center_node, 1))             # ACCELERATION NODE 02
            stiffness.append(np.abs(reaction[-1] / disp[-1]))
            OMEGA.append(np.sqrt(stiffness[-1]/TOTAL_MASS))
            PERIOD.append((np.pi * 2) / OMEGA[-1])
            # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
            #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(1, 0.5*DR, DR, 0, 1)
            PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(1, PLOT=True)
            PERIOD_MIN.append(PERIODmin)
            PERIOD_MAX.append(PERIODmax)
            # EVALUATION OF DUCTILITY DAMAGE INDEX
            if MAT_TYPE == 'INELASTIC':
                di = S066.DAMAGE_INDEX_FUN(disp[-1], DY, DSU)
                DI.append(di)                   # DAMAGE INDEX
            if MAT_TYPE == 'ELASTIC':
                DI.append(0.0)
            #print(time[-1], disp[-1], velo[-1])
            print(f"Time: {time[-1]:.4f}, Displacement: {disp[-1]:.4f} m, Reaction: {reaction[-1]:.2f} N")

        else:
            print('\n\nSEISMIC ANALYSIS DONE.\n\n')     
        # Calculating Damping Ratio and Period Using Logarithmic Decrement Analysis 
        damping_ratio = S04.DAMPING_RATIO(disp)  
        
        # Compute modal properties
        ops.modalProperties("-print", "-file", f"SALAR_ModalReport_{ANAL_TYPE}.txt", "-unorm") 
        
        DATA = (time, reaction, disp, velo, acc, DI,
                stiffness, PERIOD, damping_ratio,
                np.array(PERIOD_MIN), np.array(PERIOD_MAX))
        
        return DATA    
#%%------------------------------------------------------
# CYCLIC DISPLACEMENT ANALYSIS (STATIC TIME-HISTORY ANALYSIS)
MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'
ANAL_TYPE = 'CYCLIC_DISPLACEMENT'

data = SDOF(MAT_TYPE, TOTAL_MASS, ANAL_TYPE, 0)
(reaction_CP, disp_CP, DI_CP,
 PERIOD_MIN_CP, PERIOD_MAX_CP) = data

Ed_CP, fig_CP = DISSIPATED_ENERGY_FUN_WITH_PLOT(
            disp_CP, reaction_CP, method = 2,
            title="Cyclic Loading – Dissipated Energy"
        )
fig_CP.show()
        
print(f"Dissipated Energy from Cyclic Displacement= {Ed_CP:.2f} N·m")

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
    8: [],  # DISSIPATED ENERGY CAPACITY INDEX
    9: [],  # EQULIVALENT VISCOUS DAMPING RATIO
    }

# NUM_SIM is the number of simulations
for i in range(NUM_SEISMIC):
    # Initialize lists to store max values
    max_time = []
    max_displacement = []
    max_velocity = []
    max_acceleration = []
    max_base_reaction = []
    max_DI = []
    max_T = []
    max_K = []
    max_dr = []
    max_edci = [] # DISSIPATED ENERGY CAPACITY INDEX
    max_zeta = [] # EQULIVALENT VISCOUS DAMPING RATIO
    #g = 2*9.81 * (j+1 / NUM_G) # [m/s^2] standard acceleration of gravity
    for j in range(NUM_T):
        m = TOTAL_MASS * (j+1 / NUM_T)   # [kg] mass of structure
        
        # SEISMIC ANALYSIS (DYNAMIC TIME-HISTORY ANALYSIS)
        MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'
        ANAL_TYPE = 'SEISMIC'

        data = SDOF(MAT_TYPE, m, ANAL_TYPE, i)

        (time_SEI, reaction_SEI, disp_SEI, velo_SEI, acc_SEI, DI_SEI,
         stiffness_SEI, PERIOD_SEI, damping_ratio_SEI,
         PERIOD_MIN_SEI, PERIOD_MAX_SEI) = data
        
        #%% START - Energy Dissipation Capacity Index
        Ed_SEI, fig_SEI = DISSIPATED_ENERGY_FUN_WITH_PLOT(
            disp_SEI, reaction_SEI, method = 2, 
            title="Earthquake Response – Dissipated Energy"
        )
        fig_SEI.show()
        
        print(f"Dissipated Energy from Earthquake= {Ed_SEI:.2f} N·m")
        
        DECI = 100 * Ed_SEI / Ed_CP
        
        if DECI <= 100:
            print(f'\n\tDISSIPATED ENERGY CAPACITY INDEX: {DECI:.3f} [%]\n')
        else:
            print('\n\tFOR EVALUATION OF DISSIPATED ENERGY CAPACITY INDEX:')
            print('\n\tCHECK THE CYCLIC DISPLACEMENT ANALYSIS AND IF IT IS POSSIBLE')
            print('\t\t\tINCREASE THE DISPLACEMENT.\n')
            
        if DECI <= 0:
            print("\n\tZONE 0: NO DAMAGE\n")
        elif DECI > 0 and DECI <= 10:
            print("\n\tZONE 1: VERY MINOR DAMAGE\n")
        elif DECI > 10 and DECI <= 20:
            print("\n\tZONE 2: MINOR DAMAGE\n")
        elif DECI > 20 and DECI <= 30:
            print("\n\tZONE 3: MODERATE–LOW DAMAGE\n")
        elif DECI > 30 and DECI <= 40:
            print("\n\tZONE 4: MODERATE DAMAGE\n")
        elif DECI > 40 and DECI <= 50:
            print("\n\tZONE 5: MODERATE–HIGH DAMAGE\n")
        elif DECI > 50 and DECI <= 60:
            print("\n\tZONE 6: SEVERE–LOW DAMAGE\n")
        elif DECI > 60 and DECI <= 70:
            print("\n\tZONE 7: SEVERE–MEDIUM DAMAGE\n")    
        elif DECI > 70 and DECI <= 80:
            print("\n\tZONE 8: SEVERE–HIGH DAMAGE\n")
        elif DECI > 80 and DECI <= 90:
            print("\n\tZONE 9: VERY SEVERE DAMAGE\n")
        elif DECI > 90 and DECI <= 100:
            print("\n\tZONE 10: FAILURE DAMAGE\n")
       #%% END - Energy Dissipation Capacity Index
       
        #%% START - Energy Dissipation Capacity Index
        XLABEL = "Displacement [m]"
        YLABEL = "Base Reaction [N]"
        TITLE = "Equivalent viscous damping ratio - Incremental Dynamic Analysis Hysteresis"
        method = 1 # 
        zeta = S055.EQULIVALENT_VISCOUS_DAMPING_RATIO_FUN(disp_SEI, reaction_SEI, method, XLABEL, YLABEL, TITLE)
        print(f"Equivalent viscous damping ratio = {zeta:.4f}")
        #%% END - Energy Dissipation Capacity Index
        
        # Calculate and store the max absolute values
        max_time.append(np.max(np.abs(time_SEI)))
        max_displacement.append(np.max(np.abs(disp_SEI)))
        max_velocity.append(np.max(np.abs(velo_SEI)))
        max_acceleration.append(np.max(np.abs(acc_SEI)))
        max_base_reaction.append(np.max(np.abs(reaction_SEI)))
        max_DI.append(np.max(np.abs(DI_SEI)))
        max_K.append(np.max(np.abs(stiffness_SEI)))
        max_T.append(np.max(PERIOD_MAX_SEI)) 
        max_dr.append(damping_ratio_SEI)
        max_edci.append(DECI)
        max_zeta.append(zeta)
        
        print(f'SEISMIC {j + 1} - STEP {i + 1} DONE')
    # Store Data
    DATA[1].append(max_displacement)
    DATA[2].append(max_velocity)
    DATA[3].append(max_acceleration)
    DATA[4].append(max_base_reaction)
    DATA[5].append(max_DI)   
    DATA[6].append(max_dr)
    DATA[7].append(max_K)
    DATA[8].append(max_edci)
    DATA[9].append(max_zeta)
    print(f'\n SEISMIC {j + 1} DONE \n')    
else:
    print('\n\n Analysis completed successfully')    

totaltime = ti.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n') 

#%%------------------------------------------------------
# OUTPUT DATA FROM RESPONSE SPECTRUM ANALYSIS TO EXCEL FILE
labels = {
    1: "DISPLACEMENT",
    2: "VELOCITY",
    3: "ACCELERATION",
    4: "BASE_REACTION",
    5: "DAMAGE_INDEX",
    6: "DAMPING_RATIO",
    7: "STIFFNESS",
    8: "DISSIPATED_ENERGY_CAPACITY_INDEX",
    9: "EQULIVALENT_VISCOUS_DAMPING_RATIO",
}

with pd.ExcelWriter("INELASTIC_RESPONSE_SPECTRUM_SEISMIC_SDOF_RESULTS.xlsx", engine="openpyxl") as writer:
    for key, value in DATA.items():
        df = pd.DataFrame(value)
        sheet_name = labels[key][:31] # SHEEET NAME
        df.to_excel(writer, sheet_name=sheet_name, index=False)

print("Excel file saved: SDOF_RESULTS.xlsx")

#%%------------------------------------------------------
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
def PLOT_2D(COUNT, XLABEL, YLABEL, TITLE):
    plt.figure(COUNT, figsize=(12, 10))

    # Plot all simulations
    for j in range(NUM_SEISMIC):
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

    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

    
    
# ----------------------------------------
# Plot 1: Displacement Response Spectrum
# ----------------------------------------
X = T
Y = DATA
XLABEL = "Structural period [s]"
YLABEL = "Max Displacement (PGD) [m]"
TITLE =  "Displacement Response Spectrum"

PLOT_2D(1, XLABEL, YLABEL, TITLE)

# ----------------------------------------
# Plot 2: Velocity Response Spectrum
# ----------------------------------------
X = T
Y = DATA
XLABEL = "Structural period [s]"
YLABEL = "Max Velocity (PGV) [m/s]"
TITLE =  "Velocity Response Spectrum"

PLOT_2D(2, XLABEL, YLABEL, TITLE)

# ----------------------------------------
# Plot 3: Acceleration Response Spectrum
# ----------------------------------------
X = T
Y = DATA
XLABEL = "Structural period [s]"
YLABEL = "Max Acceleration (PGA) [m/s²]"
TITLE =  "Acceleration Response Spectrum"

PLOT_2D(3, XLABEL, YLABEL, TITLE)

# ----------------------------------------
# Plot 4: Base Reaction Spectrum
# ----------------------------------------
X = T
Y = DATA
XLABEL = "Structural period [s]"
YLABEL = "Max Base Reaction (N)"
TITLE =  "Base Reaction Response Spectrum"

PLOT_2D(4, XLABEL, YLABEL, TITLE)

# ----------------------------------------
# Plot 5: Damage Index Spectrum
# ----------------------------------------
X = T
Y = DATA
XLABEL = "Structural period [s]"
YLABEL = "Damage Index (DI) [%]"
TITLE =  "Ductility Damage Index Spectrum"

PLOT_2D(5, XLABEL, YLABEL, TITLE)

# ----------------------------------------
# Plot 6: Damping Ratio Spectrum
# ----------------------------------------
X = T
Y = DATA
XLABEL = "Structural period [s]"
YLABEL = "Damping Ratio [%]"
TITLE =  "Damping Ratio Spectrum"

PLOT_2D(6, XLABEL, YLABEL, TITLE)

# ----------------------------------------
# Plot 7: Structural Stiffness Spectrum
# ----------------------------------------
X = T
Y = DATA
XLABEL = "Structural period [s]"
YLABEL = "Structural Stiffness Spectrum [N/m]"
TITLE =  "Structural Stiffness Spectrum"

PLOT_2D(7, XLABEL, YLABEL, TITLE)

# ---------------------------------------------------
# Plot 8: Energy Dissipation Capacity Index Spectrum
# ---------------------------------------------------
X = T
Y = DATA
XLABEL = "Structural period [s]"
YLABEL = "Energy Dissipation Capacity Index [%]"
TITLE =  "Energy Dissipation Capacity Index Spectrum"

PLOT_2D(8, XLABEL, YLABEL, TITLE)

# ---------------------------------------------------
# Plot 9: Equivalent viscous damping Ratio Spectrum
# ---------------------------------------------------
X = T
Y = DATA
XLABEL = "Structural period [s]"
YLABEL = "Equivalent viscous damping ratio [%]"
TITLE =  "Equivalent viscous damping ratio Spectrum"

PLOT_2D(9, XLABEL, YLABEL, TITLE)

#%%------------------------------------------------------
# Print the last results
print("Maximum Absolute Values Across Simulations:")
print("Period:", max_time[-1])
print("Displacement:", max_displacement[-1])
print("Velocity:", max_velocity[-1])
print("Acceleration:", max_acceleration[-1])
print("Base Reaction:", max_base_reaction[-1])
print("Ductility Damage Index:", max_DI[-1])
print("Energy Dissipation Capacity Index:", max_edci[-1])
print("Equivalent viscous damping ratio :", max_zeta[-1])
print("Period :", max_T[-1])
#%%------------------------------------------------------
S01.HISROGRAM_BOXPLOT(max_displacement, HISTO_COLOR='blue', LABEL='Displacement')
S01.HISROGRAM_BOXPLOT(max_velocity, HISTO_COLOR='purple', LABEL='Velocity')
S01.HISROGRAM_BOXPLOT(max_acceleration, HISTO_COLOR='green', LABEL='Acceleration')
S01.HISROGRAM_BOXPLOT(max_base_reaction, HISTO_COLOR='gold', LABEL='Base Reaction')
S01.HISROGRAM_BOXPLOT(max_DI, HISTO_COLOR='pink', LABEL='Ductility Damage Index')
S01.HISROGRAM_BOXPLOT(max_T, HISTO_COLOR='lime', LABEL='Structure Period')
S01.HISROGRAM_BOXPLOT(max_edci, HISTO_COLOR='cyan', LABEL='Energy Dissipation Capacity Index')
S01.HISROGRAM_BOXPLOT(max_zeta, HISTO_COLOR='brown', LABEL='Equivalent viscous damping ratio')
#%%------------------------------------------------------
####  FRAGILITY ANALYSIS
  
# ----------------------------
# Fragility Assessment
# ----------------------------
# Define damage states per FEMA P-58
# INFO LINK: https://www.fema.gov/sites/default/files/documents/fema_p-58-2-se_volume2_implementation.pdf
damage_states = {
'DS1_Slight': (0.15, 0.4),    # Median PGA=0.15g, β=0.4
'DS2_Moderate': (0.30, 0.4),
'DS3_Extensive': (0.60, 0.4),
'DS4_Complete': (1.00, 0.4)
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
plt.plot(time_SEI, acc_SEI, lw=1, color='black')
plt.xlabel('Time (s)')
plt.ylabel('Acceleration (g)')
plt.title(f'Last Analysis Structural Response + Ground Motion ::: MAX. ABS. : {np.max(np.abs(acc_SEI)):.4f}')
plt.grid(True)
plt.show()    

im_values = max_acceleration 
X_LABEL = 'Peak Ground Acceleration (g)  [IM]'
S077.FRAGILITY_CURVE(im_values, damage_states, X_LABEL, SEMILOGY=False, PLOT_KIND=True)
#===========================================================
# Define damage state parameters: {Damage State: (median_IM, beta)}
damage_states = {
    'Minor Damage Level': (20.0, 0.4),# Median DI=20, β=0.4
    'Moderate Damage Level': (40.0, 0.4),
    'Severe Damage Level': (60.0, 0.4),
    'Failure Level': (100.0, 0.4)
}

# Generate intensity measure (IM) values from 0.0 to 100.0
im_values = max_edci # Energy Dissipation Capacity Index
X_LABEL = 'Energy Dissipation Capacity Index (%)  [IM]'
S077.FRAGILITY_CURVE(im_values, damage_states, X_LABEL, SEMILOGY=False, PLOT_KIND=True)
#===========================================================
# Define damage state parameters: {Damage State: (median_IM, beta)}
damage_states = {
    'Minor Damage Level': (20.0, 5.0),# Median DI=20, β=5
    'Moderate Damage Level': (40.0, 5.0),
    'Severe Damage Level': (60.0, 5.0),
    'Failure Level': (100.0, 5.0)
}

# Generate intensity measure (IM) values from 0.0 to 100.0
im_values = max_DI # Ductility Damage Index
X_LABEL = 'Ductility Damage Index (%)  [IM]'
S077.FRAGILITY_CURVE(im_values, damage_states, X_LABEL, SEMILOGY=False, PLOT_KIND=True)
#===========================================================
# Define damage state parameters: {Damage State: (median_IM, beta)}
damage_states = {
    'Minor Damage Level': (20.0, 5.0),# Median DI=20, β=5
    'Moderate Damage Level': (40.0, 5.0),
    'Severe Damage Level': (60.0, 5.0),
    'Failure Level': (100.0, 5.0)
}

# Generate intensity measure (IM) values from 0.0 to 100.0
im_values = max_zeta # Energy Dissipation Capacity Index
X_LABEL = 'Equivalent viscous damping ratio (%)  [IM]'
S077.FRAGILITY_CURVE(im_values, damage_states, X_LABEL, SEMILOGY=False, PLOT_KIND=True)
#%%------------------------------------------------------
import matplotlib.pyplot as plt

plt.figure(figsize=(12,8))
plt.plot(disp_SEI, reaction_SEI,color='black')
plt.xlabel("Displacement [m]")
plt.ylabel("Base Reaction [N]")
plt.title("Displacement & Base Reaction Relation From Last Dynamic Analysis")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
#%%------------------------------------------------------
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

DISP_Z = MAX_ABS(disp_SEI)  
VELO_Z = MAX_ABS(velo_SEI) 
ACCE_Z = MAX_ABS(acc_SEI) 
BASE_Z = MAX_ABS(reaction_SEI) 

plt.figure(1, figsize=(8, 6))
plt.plot(time_SEI, disp_SEI, color='blue', linewidth=2)
plt.plot(time_SEI, DISP_Z, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Displacement [m]')
plt.title(f'Time vs Displacement - MAX. ABS: {DISP_Z[-1]}')
plt.grid()
plt.show()

plt.figure(2, figsize=(8, 6))
plt.plot(time_SEI, velo_SEI, color='blue', linewidth=2)
plt.plot(time_SEI, VELO_Z, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Velocity [m/s]')
plt.title(f'Time vs Velocity - MAX. ABS: {VELO_Z[-1]}')
plt.grid()
plt.show()

plt.figure(3, figsize=(8, 6))
plt.plot(time_SEI, acc_SEI, color='blue', linewidth=2)
plt.plot(time_SEI, ACCE_Z, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Acceleration [m/s^2]')
plt.title(f'Time vs Acceleration - MAX. ABS: {ACCE_Z[-1]}')
plt.grid()
plt.show()

plt.figure(4, figsize=(8, 6))
plt.plot(time_SEI, reaction_SEI, color='blue', linewidth=2)
plt.plot(time_SEI, BASE_Z, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Base-reaction [N]')
plt.title(f'Time vs Base-reaction - MAX. ABS: {BASE_Z[-1]}')
plt.grid()
plt.show()
#%%------------------------------------------------------  
# EXPORT DATA TO EXCEL
DATA_TOTAL = {
    'Max_displacement': max_displacement,
    'Max_velocity': max_velocity,
    'Max_acceleration': max_acceleration,
    'Max_Base_Reaction': max_base_reaction,
    'Ductility_Damage_Index': max_DI,
    'Period': max_T,
    'Energy Dissipation Capacity Index': max_edci,
    'Equivalent viscous damping ratio': max_zeta,
}
# Convert to DataFrame
results_df = pd.DataFrame(DATA_TOTAL)
# Export the DataFrame to an Excel file
results_df.to_excel('SDOF_RESPONSE_SPECTRUM_EDCI_EVDR_DATA_RESULTS.xlsx', index=False)
#%%------------------------------------------------------  
XLABEL = 'Displacement'
YLABEL = 'Base Reaction'
TITLE = f'{YLABEL} and {XLABEL} scatter chart'
COLOR = 'orange'
X = max_displacement
Y = max_base_reaction
S01.PLOT_SCATTER(X, Y , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

# CLUSTER DATA
S01.CLUSTER_DATA(X, Y, XLABEL, YLABEL, MAX_CLUSTERS=3)
#%%------------------------------------------------------
XLABEL = 'Velocity'
YLABEL = 'Base Reaction'
TITLE = f'{YLABEL} and {XLABEL} scatter chart'
COLOR = 'cyan'
X = max_velocity
Y = max_base_reaction
S01.PLOT_SCATTER(X, Y , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

# CLUSTER DATA
S01.CLUSTER_DATA(X, Y, XLABEL, YLABEL, MAX_CLUSTERS=3)
#%%------------------------------------------------------
XLABEL = 'Acceleration'
YLABEL = 'Base Reaction'
TITLE = f'{YLABEL} and {XLABEL} scatter chart'
COLOR = 'lime'
X = max_acceleration
Y = max_base_reaction
S01.PLOT_SCATTER(X, Y , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

# CLUSTER DATA
S01.CLUSTER_DATA(X, Y, XLABEL, YLABEL, MAX_CLUSTERS=3)
#%%------------------------------------------------------
XLABEL = 'Displacement'
YLABEL = 'Structural Ductility Damage Index (%)'
TITLE = f'{YLABEL} and {XLABEL} scatter chart'
COLOR = 'purple'
X = max_displacement
Y = max_DI
S01.PLOT_SCATTER(X, Y , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

# CLUSTER DATA
S01.CLUSTER_DATA(X, Y, XLABEL, YLABEL, MAX_CLUSTERS=3)
#%%------------------------------------------------------
# PLOT THE TIME-HISTORY
S01.PLOT_TIME_HISTORY(time_SEI, disp_SEI, velo_SEI, acc_SEI, reaction_SEI)
#%%------------------------------------------------------
# RANDOM FOREST ANALYSIS
"""
This code predicts the seismic safety of a structure using simulation data by training a Random Forest Classifier to
 classify whether the system is "safe" or "unsafe" based on features like maximum displacement, velocity, acceleration,
 and base reaction. A regression model is also trained to estimate safety likelihood. It evaluates model performance using
 metrics like classification accuracy, mean squared error, and R² score. Additionally, it identifies key features influencing
 safety through feature importance analysis. The tool aids in seismic risk assessment, structural optimization, and understanding
 critical safety parameters.
"""

data = {
    'Max_displacement': max_displacement,
    'Max_velocity': max_velocity,
    'Max_acceleration': max_acceleration,
    'Max_Base_Reaction': max_base_reaction,
    'Ductility_Damage_Index': max_DI,
    'Energy_Dissipation_Capacity_Index': max_edci,
    'Equivalent_viscous_damping_ratio': max_zeta,    
}


# Convert to DataFrame
df = pd.DataFrame(data)
#print(df)
S01.RANDOM_FOREST(df)
#%%------------------------------------------------------
# PLOT HEATMAP FOR CORRELATION 
S01.PLOT_HEATMAP(df)
#%%------------------------------------------------------
# MULTIPLE REGRESSION MODEL
S01.MULTIPLE_REGRESSION(df) 
