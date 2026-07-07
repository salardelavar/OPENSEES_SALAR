###########################################################################################################
#                   >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                      #
#     SENSITIVITY ANALYSIS BY CHANGING STRUCTURAL DAMPER STIFFNESS AND DAMPING RATIO USING OPENSEES       #
#---------------------------------------------------------------------------------------------------------#
#                    THIS PYTHON SCRIPT WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                     #
#                                   EMAIL: salar.d.ghashghaei@gmail.com                                   #
###########################################################################################################

"""
Suppose we have an existing steel or concrete structure and wish to
 enhance its seismic performance by adding supplemental dampers.
 The goal is to increase the overall effective damping of the structural 
 system through the contribution of these devices. In this type of study
 a sensitivity analysis is carried out by varying key damper parameters;
 such as the elastic stiffness of the damper and its viscous damping
 coefficient—in order to evaluate the resulting total damping of the
 structure. This total damping combines the inherent damping of the
 original structure and the additional damping provided by the dampers.
 
Nonlinear Seismic Performance Assessment of a SDOF:
An OpenSeesPy Framework for Material and Geometric Nonlinearity Under Static, Cyclic, and Earthquake Loading

This OpenSeesPy script performs rigorous nonlinear static and dynamic analysis of a SDOF
 for performance-based earthquake engineering. The 2D model incorporates both material nonlinearity
 (elastic-perfectly plastic or hysteretic steel with strain hardening)
 and geometric nonlinearity (corotational truss formulation) to capture P-delta effects
 and large displacements—essential for collapse assessment.
 
Eight analysis protocols are implemented:    
(1) [PERIOD] : Structural Period
(2) [STATIC] : Gravity load analysis establishing dead load state
(3) [PUSHOVER] : Displacement-controlled pushover generating full capacity curves
 and plastic mechanism identification
(4) [CYCLIC_DISPLACEMENT] : Symmetric cyclic displacement protocol capturing hysteresis,
 pinching behavior, and energy dissipation degradation
(5) [STATIC_EXTERNAL_TIME-DEPENDENT_LOADING] : Static Analysis of External time-dependent loading P(t) = P0 sin(wt) or P(t) = P0 exp(-0.05wt) sin(wt)  
(6) [DYNAMIC_EXTERNAL_TIME-DEPENDENT_LOADING] : Dynamic Analysis of External time-dependent loading P(t) = P0 sin(wt) or P(t) = P0 exp(-0.05wt) sin(wt)  
(7) [FREE-VIBRATION] : Free-vibration with initial conditions extracting damping ratios
 via logarithmic decrement
(8) [SEISMIC] : Multi-directional seismic excitation with Rayleigh damping (3% ratio)
 and uniform acceleration patterns.

The code continuously records displacements, enabling member-level demand-to-capacity
 ratio tracking. Period monitoring during inelastic response reveals softening
 and potential period shift into resonant ground motion frequency ranges.
 This framework enables vulnerability curve development, seismic fragility assessment,
 and retrofit prioritization—directly supporting next-generation bridge seismic design codes
 and risk-informed asset management decisions.
 
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
"""
"""
OPTIMIZATION PROBLEM FOR ENERGY DISSIPATION CAPACITY INDEX (EDCI)
-----------------------------------------------------------------
1. Problem Formulation – We invert the conventional design process:
   instead of checking a given yield strength Fy, we find the 'minimum'
   Fy that makes the structure's hysteretic energy dissipation under
   a specific earthquake exactly equal to a target fraction (20%) of
   its maximum dissipation capacity measured from a full cyclic test.

2. Energy‑Based Performance Metric – The Energy Dissipation Capacity
   Index (EDCI) = E_seismic / E_cyclic is a scalar that quantifies
   how close the seismic demand comes to saturating the structure's
   plastic energy absorption capability – a direct proxy for cumulative
   damage and collapse margin.

3. Nonlinear Forward Model – Each evaluation of Fy requires solving
   the equation of motion: M·ü + C·u̇ + R(u, Fy) = f_ext(t) with a
   rate‑independent hysteretic material (elastic‑perfectly‑plastic or
   with hardening), where the restoring force R depends sensitively on
   Fy through the yield surface.

4. Root‑Finding via Newton–Raphson – The residual g(Fy) = EDCI(Fy)-0.20
   is driven to zero using a secant‑like Newton iteration, where the
   Jacobian g'(Fy) is approximated by central finite differences
   [g(Fy+ε) - g(Fy-ε)]/(2ε), with ε tuned to balance truncation and
   conditioning errors.

5. Numerical Challenges – The function g is non‑smooth due to abrupt
   yielding, pinching, and stiffness degradation; hence the solver is
   augmented with a line‑search and a fallback bisection step when the
   Newton update overshoots, ensuring robust convergence despite
   discontinuities.

6. Cyclic Capacity Envelope – The denominator E_cyclic is obtained from
   a quasi‑static displacement‑controlled cyclic protocol that
   progressively increases amplitude to capture the full hysteresis
   loop, including the post‑peak softening branch – defining the
   structure's ultimate energy dissipation reservoir.

7. Seismic Demand Evaluation – The numerator E_seismic is computed from
   a full dynamic time‑history analysis with Rayleigh damping (3%) and
   a scaled ground motion, integrating the instantaneous power
   ∫ R·du over the entire shaking duration.

8. Physical Insight – The converged Fy corresponds to the strength
   threshold that allows the structure to exploit its full ductility
   without exceeding the damage index associated with the target EDCI;
   below this value, the system would enter the "failure" zone
   (EDCI > 100%) under the same excitation.

9. Sensitivity and Stability – The finite‑difference step ε is set to
   1e-3·Fy to avoid numerical ill‑conditioning, and the solver terminates
   when the relative change in Fy falls below 1e-6, ensuring that the
   design strength is determined with sub‑percent accuracy.

10. Engineering Decision – This optimization yields a rational,
    risk‑consistent design point that directly links seismic hazard
    (energy input) to structural capacity (energy absorption), enabling
    performance‑based earthquake engineering without empirical
    overstrength factors – a move towards collapse‑prevention
    verification based on explicit energy balance.
"""
"""
The numbers below represent the equivalent viscous damping ratio relative to the critical damping of the overall structure:

- Viscous damper (fluid damper): 15% to 35% (even higher in special projects)
- Friction damper: 15% to 30%
- Tuned mass damper (TMD): Typically adds 1% to 5% (mainly for reducing wind vibrations and resonance)
- Metallic yielding damper (ADAS, TADAS): 10% to 25% (with yielding of steel)
- Viscoelastic damper: 10% to 20%
"""
# PAPER: Solution of Navier-Stokes equations for fluids with magnetorheological compensation used in structures with energy dissipaters
'https://www.researchgate.net/publication/357612213_Solution_of_Navier-Stokes_equations_for_fluids_with_magnetorheological_compensation_used_in_structures_with_energy_dissipaters?_tp=eyJjb250ZXh0Ijp7ImZpcnN0UGFnZSI6Il9kaXJlY3QiLCJwYWdlIjoiX2RpcmVjdCJ9fQ'    
# PAPER: Investigation of the seismic behaviours of three-dimensional high-rise steel frame structures equipped with oil dampers with variable stiffness
'https://www.sciencedirect.com/science/article/abs/pii/S0143974X21000237'    
# PAPER: Experimental Study on Two Full Scale Iranian Viscous Dampers
'https://www.researchgate.net/figure/llustrates-a-patented-viscous-damper-manufactured-by-Behsazan-Larzeh-Davam-Co-As_fig1_330081268'    
# PAPER: A Compact Variable Stiffness and Damping Shock Absorber for Vehicle Suspension
'https://ieeexplore.ieee.org/document/7086317'    
#%%----------------------------------------------------
import numpy as np
import openseespy.opensees as ops
import matplotlib.pyplot as plt 
import PLOT_2D_TRUSS as S01
import ANALYSIS_FUNCTION as S02
import PERIOD_FUN as S03
import DAMPING_RATIO_FUN as S04
import EIGENVALUE_ANALYSIS_FUN as S05
import RAYLEIGH_DAMPING_FUN as S06
import BILINEAR_CURVE as S07
import OPENSEEES_HYSTERETICSM_FORCE_DISP_FUN as S08
import DISSIPATED_ENERGY_WITH_PLOT_FUN as S09
import SALAR_MATH as S10
#%%----------------------------------------------------
def SDOF(Kd, DRx, MAT_TYPE, TOTAL_MASS, ANAL_TYPE):
    # Initialize model
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)
    
    MAX_ITERATIONS = 5000   # Maximum number of iterations
    MAX_TOLERANCE = 1.0e-6  # Specified tolerance for convergence
    GMfact = 9.81 # [m/s^2] standard acceleration of gravity or standard acceleration 
        
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
    DR = 0.02                                        # Damping Ratio
    # Positive branch points
    pos_disp = [0, DY, DSU, 1.1*DSU, 1.25*DSU]
    pos_force = [0, FY, FU, 0.2*FU, 0.1*FU]
    KP = np.array([FY, DY, FU, DSU, 0.2*FU, 1.1*DSU, 0.1*FU, 1.25*DSU])
    
    # Negative branch points
    neg_disp = [0, -DY, -DSU, -1.1*DSU, -1.25*DSU]
    neg_force = [0, -FY, -FU, -0.2*FU, -0.1*FU]
    KN = np.array([-FY, -DY, -FU, -DSU, -0.2*FU, -1.1*DSU, -0.1*FU, -1.25*DSU])
    
    # Plot
    plt.figure(figsize=(10, 6))
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
    
    MatTag_C = 2     # SPRING DAMPER
    alpha = 1.0      # velocity exponent (usually 0.3–1.0)
    #Kd = 100_000_000 # ELASTIC STIFFNESS OF DAMPER
    DRd = DR + DRx    # DAMPER DAMPING RATIO
    omega = np.sqrt(Kd/TOTAL_MASS)
    Cd = 2 * DRd  * omega * TOTAL_MASS  # [N·s/m] Damping coefficient 
    #ops.uniaxialMaterial('Viscous', MatTag_C, Cd, alpha)  # Material for C (alpha=1.0 for linear)
    #ops.uniaxialMaterial('ViscousDamper', MatTag_C, Kd Cd, alpha)
    # INFO LINK: https://openseespydoc.readthedocs.io/en/stable/src/ViscousDamper.html
    ops.uniaxialMaterial('BilinearOilDamper', MatTag_C, Kd, Cd)# , Fr=1.0, p=1.0, LGap=0.0, NM=1, RelTol=1e-6, AbsTol=1e-10, MaxHalf=15
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/BilinearOilDamper.html
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
        DI.append(100*(np.abs(disp[-1])-DY)/(DSU-DY))                   # DAMAGE INDEX
        if DI[-1] <= 0.0: # DAMAGE INDEX -> 0.0 <= DI <=100 
            DI[-1] = 0.0
        if DI[-1] >= 100: 
            DI[-1] = 100.0    
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
            DI.append(100*(np.abs(disp[-1])-DY)/(DSU-DY))                   # DAMAGE INDEX
            if DI[-1] <= 0.0: # DAMAGE INDEX -> 0.0 <= DI <=100 
                DI[-1] = 0.0
            if DI[-1] >= 100: 
                DI[-1] = 100.0
            STEP += 1
            #print(STEP, disp[-1], reaction[-1])
            print(f"Step: {STEP}, Displacement: {disp[-1]:.4f} mm, Reaction: {reaction[-1]:.2f} N")     
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
             0.1*DMAX,     -0.1*DMAX,
             0.5*DMAX,     -0.5*DMAX,
             0.8*DMAX,     -0.8*DMAX,
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
            DI.append(100*(np.abs(disp[-1])-DY)/(DSU-DY))                   # DAMAGE INDEX
            if DI[-1] <= 0.0: # DAMAGE INDEX -> 0.0 <= DI <=100 
                DI[-1] = 0.0
            if DI[-1] >= 100: 
                DI[-1] = 100.0
            STEP += 1
            #print(STEP, disp[-1], reaction[-1])
            print(f"Step: {STEP}, Displacement: {disp[-1]:.4f} mm, Reaction: {reaction[-1]:.2f} N")     
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
            DI.append(100*(np.abs(disp[-1])-DY)/(DSU-DY))                   # DAMAGE INDEX
            if DI[-1] <= 0.0: # DAMAGE INDEX -> 0.0 <= DI <=100 
                DI[-1] = 0.0
            if DI[-1] >= 100: 
                DI[-1] = 100.0
            STEP += 1
            #print(STEP, disp[-1], reaction[-1])
            print(f"Step: {STEP}, Displacement: {disp[-1]:.4f} mm, Reaction: {reaction[-1]:.2f} N")     
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
            DI.append(100*(np.abs(disp[-1])-DY)/(DSU-DY))                   # DAMAGE INDEX
            if DI[-1] <= 0.0: # DAMAGE INDEX -> 0.0 <= DI <=100 
                DI[-1] = 0.0
            if DI[-1] >= 100: 
                DI[-1] = 100.0
            #print(time[-1], disp[-1], velo[-1])
            print(f"Time: {time[-1]:.4f}, Displacement: {disp[-1]:.4f} mm, Reaction: {reaction[-1]:.2f} N")      
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
        u0 = -0.010                      # [m] Initial displacement
        v0 = 0.0                           # [m/s] Initial velocity
        a0 = 0.0                           # [m/s^2] Initial acceleration
        IU = True                          # Free Vibration with Initial Displacement
        IV = True                          # Free Vibration with Initial Velocity
        IA = True                          # Free Vibration with Initial Acceleration
        duration = 50.0                    # [s] Analysis duration
        dt = 0.01                          # [s] Time step
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
            #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(1, 0.5*DRx, DRx, 0, 1)
            PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(1, PLOT=True)
            PERIOD_MIN.append(PERIODmin)
            PERIOD_MAX.append(PERIODmax)
            DI.append(100*(np.abs(disp[-1])-DY)/(DSU-DY))          # DAMAGE INDEX
            if DI[-1] <= 0.0: # DAMAGE INDEX -> 0.0 <= DI <=100 
                DI[-1] = 0.0
            if DI[-1] >= 100: 
                DI[-1] = 100.0
            #print(time[-1], disp[-1], velo[-1])
            print(f"Time: {time[-1]:.4f}, Displacement: {disp[-1]:.4f} mm, Reaction: {reaction[-1]:.2f} N")      
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
        GMfact = 9810                      # [mm/s²]standard acceleration of gravity or standard acceleration
        SSF_X = 0.10                       # Seismic Acceleration Scale Factor in X Direction
        SSF_Y = 0.10                       # Seismic Acceleration Scale Factor in Y Direction
        iv0_X = 0.0005                     # [mm/s] Initial velocity applied to the node  in X Direction
        iv0_Y = 0.0005                     # [mm/s] Initial velocity applied to the node  in Y Direction
        st_iv0 = 0.0                       # [s] Initial velocity applied starting time
        SEI = 'X'                          # Seismic Direction
        DR = 0.03                          # Damping ratio
        
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
            DI.append(100*(np.abs(disp[-1])-DY)/(DSU-DY))         # DAMAGE INDEX
            if DI[-1] <= 0.0: # DAMAGE INDEX -> 0.0 <= DI <=100 
                DI[-1] = 0.0
            if DI[-1] >= 100: 
                DI[-1] = 100.0
            #print(time[-1], disp[-1], velo[-1])
            print(f"Time: {time[-1]:.4f}, Displacement: {disp[-1]:.4f} mm, Reaction: {reaction[-1]:.2f} N")

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

#%%-------------------------------------------------------
def PLOT_TIME_HISTORY(time, reaction, disp, velo, acc):

    import matplotlib.pyplot as plt

    data_dict = {
        "Reaction [N]": reaction,
        "Displacement [m]": disp,
        "Velocity [m/s]": velo,
        "Acceleration [m/s²]": acc,
    }

    fig, axes = plt.subplots(4, 2, figsize=(15, 12))
    axes = axes.flatten()

    for i, (label, data) in enumerate(data_dict.items()):
        axes[i].plot(time, data, color='black', linewidth=2)
        axes[i].set_title(label, fontsize=9)
        axes[i].set_xlabel("Time (s)")
        axes[i].grid(True)

  
    for j in range(len(data_dict), len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout()
    plt.show()


def PLOT(XDATA, YDATA, TITLE, XLABEL, YLABEL, COLOR, SEMILOGY):
    import matplotlib.pyplot as plt
    fig = plt.figure(-1, figsize=(12, 8))
    plt.plot(XDATA, YDATA, color=COLOR, linewidth=2)
    plt.title(TITLE)
    plt.ylabel(YLABEL)
    plt.xlabel(XLABEL)
    if SEMILOGY == True:
        plt.semilogy()
    plt.grid()
    
# Plotting Nodal Displacements
def PLOT_DISPLAEMENTS(time_steps, displacements_dict, TITLE):
    plt.figure(figsize=(10, 6))
    
    for node_id, disp_values in displacements_dict.items():
        plt.plot(time_steps, disp_values, label=f'Node {node_id} - MAX. ABS. : {np.max(np.abs(disp_values)): 0.4e}', linewidth=2)
    
    plt.xlabel('Time [s]')
    plt.ylabel('Displacement [m]')
    plt.title(TITLE)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()  
    
# Plotting Element Forces
def PLOT_FORCES(time_steps, forces_dict, YLABEL, TITLE):
    plt.figure(figsize=(10, 6))
    
    for node_id, force_values in forces_dict.items():
        plt.plot(time_steps, force_values, label=f'Ele. {node_id} - MAX. ABS. : {np.max(np.abs(force_values)): 0.4e}', linewidth=2)
    
    plt.xlabel('Time [s]')
    plt.ylabel(YLABEL)
    plt.title(TITLE)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()  
    
def PLOT_2D(X, Y, Xfit, Yfit, X2, Y2, XLABEL, YLABEL, TITLE, LEGEND01, LEGEND02, LEGEND03, COLOR, Z):
    import matplotlib.pyplot as plt
    plt.figure(figsize=(12, 8))
    if Z == 1:
        # Plot 1 line
        plt.plot(X, Y,color=COLOR)
        plt.xlabel(XLABEL)
        plt.ylabel(YLABEL)
        plt.title(TITLE)
        plt.grid(True)
        plt.show()
    if Z == 2:
        # Plot 2 lines
        plt.plot(X, Y, Xfit, Yfit, 'r--', linewidth=3)
        plt.title(TITLE)
        plt.xlabel(XLABEL)
        plt.ylabel(YLABEL)
        plt.legend([LEGEND01, LEGEND02], loc='lower right')
        plt.grid(True)
        plt.show()
    if Z == 3:
        # Plot 3 lines
        plt.plot(X, Y, Xfit, Yfit, 'r--', X2, Y2, 'g-*', linewidth=3)
        plt.title(TITLE)
        plt.xlabel(XLABEL)
        plt.ylabel(YLABEL)
        plt.legend([LEGEND01, LEGEND02, LEGEND03], loc='lower right')
        plt.grid(True)
        plt.show() 
#%%----------------------------------------------------
import time as TI
import numpy as np
from mpl_toolkits.mplot3d import Axes3D  # for 3D plotting

# Define your parameters (example values – adjust as needed)
MAT_TYPE = 'ELASTIC'

TOTAL_MASS = 5000.0   # [kg] Total Mass of Structure
X1 = 30000            # Intial Damper Stiffness of Structure 
X2 = 0.001            # Intial Damper Damping-ratio

# Number of points in each direction
nK = 5                # Number of Samples for Damper Stiffness of Structure   
nD = 200              # Number of Samples for Damper Damping-ratio   

# Analysis Durations:
starttime = TI.process_time()

# Pre-allocate 2D arrays for WW and ZZ
WW_matrix = np.zeros((nK, nD))
ZZ_matrix = np.zeros((nK, nD))
STEP = 0
# Loop over II (Kd index) and JJ (DRx index)
for i, II in enumerate(range(nK)):
    Kd = X1 * (II+1)   
    for j, JJ in enumerate(range(nD)):
        DRx = X2 * (JJ+1)   
        
        # Call SDOF – assumes it returns the tuple in the order you expect
        (time_FV, reaction_FV, disp_FV, velo_FV, acc_FV, DI_FV,
         stiffness_FV, PERIOD_FV, damping_ratio_FV,
         PERIOD_MIN_FV, PERIOD_MAX_FV) = SDOF(Kd, DRx, MAT_TYPE, TOTAL_MASS, 'FREE-VIBRATION')
        
        # Store the values
        WW_matrix[i, j] = damping_ratio_FV
        ZZ_matrix[i, j] = max(PERIOD_MAX_FV)   # if PERIOD_MAX_FV is a list/array
        STEP += 1
        print('STEP: ', STEP,' DONE.')
        
    # FOR EACH ANALYSIS CYCLE; WE CAN EVALUE AND PLOT THE HISTOGRAM
    X, HISTO_COLOR, LABEL = WW_matrix[:, j], 'cyan', 'STRUCTURAL DAMPING RATIO \n FROM FREE-VIBRATION ANALYSIS'
    S10.HISROGRAM_BOXPLOT(X, HISTO_COLOR, LABEL)
    X, HISTO_COLOR, LABEL = ZZ_matrix[:, j], 'purple', 'STRUCTURAL PERIOD \n FROM FREE-VIBRATION ANALYSIS'
    S10.HISROGRAM_BOXPLOT(X, HISTO_COLOR, LABEL)

# Create meshgrid for plotting
Kd_vals = np.array([X1 * ii for ii in range(nK)])
DRx_vals = np.array([X2 * jj for jj in range(nD)]) 
Kd_grid, DRx_grid = np.meshgrid(Kd_vals, DRx_vals, indexing='ij')

# ---- 3D Plot 1: Structural Damping Ratio vs Kd and DRx ----
fig1 = plt.figure(figsize=(10, 6))
ax1 = fig1.add_subplot(111, projection='3d')
surf1 = ax1.plot_surface(Kd_grid, DRx_grid, WW_matrix, cmap='viridis', edgecolor='none')
ax1.set_xlabel('Damper Elastic Stiffness [N/m] (Kd)')
ax1.set_ylabel('Damper Damping Ratio (DRx)')
ax1.set_zlabel('Structural Damping Ratio')
ax1.set_title('Damping Ratio vs Kd and DRx')
fig1.colorbar(surf1, shrink=0.5, aspect=10)
plt.show()

# ---- 3D Plot 2: Structural Period vs Kd and DRx ----
fig2 = plt.figure(figsize=(10, 6))
ax2 = fig2.add_subplot(111, projection='3d')
surf2 = ax2.plot_surface(Kd_grid, DRx_grid, ZZ_matrix, cmap='plasma', edgecolor='none')
ax1.set_xlabel('Damper Elastic Stiffness [N/m] (Kd)')
ax1.set_ylabel('Damper Damping Ratio (DRx)')
ax2.set_zlabel('Structural Period (s)')
ax2.set_title('Structural Period vs Kd and DRx')
fig2.colorbar(surf2, shrink=0.5, aspect=10)
plt.show()

totaltime = TI.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n')
#%%----------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 5))
cf = ax.contourf(Kd_grid, DRx_grid, WW_matrix, levels=20, cmap='viridis')
cs = ax.contour(Kd_grid, DRx_grid, WW_matrix, levels=10, colors='white', linewidths=0.5)
ax.clabel(cs, inline=True, fontsize=8)
ax.set_xlabel('Damper Elastic Stiffness [N/m] (Kd)')
ax.set_ylabel('Damper Damping Ratio (DRx)')
ax.set_title('Structural Damping Ratio – Topographic Map')
fig.colorbar(cf, ax=ax, label='Damping Ratio')
plt.tight_layout()
plt.show()

fig, ax = plt.subplots(figsize=(7, 5))
cf = ax.contourf(Kd_grid, DRx_grid, ZZ_matrix, levels=20, cmap='plasma')
cs = ax.contour(Kd_grid, DRx_grid, ZZ_matrix, levels=10, colors='black', linewidths=0.5)
ax.clabel(cs, inline=True, fontsize=8)
ax.set_xlabel('Damper Elastic Stiffness [N/m] (Kd)')
ax.set_ylabel('Damper Damping Ratio (DRx)')
ax.set_title('Structural Period – Topographic Map')
fig.colorbar(cf, ax=ax, label='Period (s)')
plt.tight_layout()
plt.show()
#%%----------------------------------------------------

