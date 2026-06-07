###########################################################################################################
#                   >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                      #
#   COMPREHENSIVE NONLINEAR SEISMIC ASSESSMENT OF A PEDESTRIAN BRIDGE STEEL VIERENDEEL GIRDERS            #
#                        AN OPENSEES FRAMEWORK FOR STATIC PUSHOVER, CYCLIC DEGRADATION,                   #
#                        STATIC TIME-HISTORY AND DYNAMIC TIME-HISTORY ANALYSIS                            #
#---------------------------------------------------------------------------------------------------------#
#                                                  P(t) = P0 sin(wt)                                      #
#                                           P(t) = P0 exp(-0.05wt) sin(wt)                                #
#---------------------------------------------------------------------------------------------------------#
#         ASSESSMENT OF DUCTILITY DAMAGE INDICES FOR STRUCTURAL ELEMENTS AND SYSTEMS AND EVALUATION       #
#                                        OF ENERGY DISSIPATION CAPACITY INDEX                             #
#---------------------------------------------------------------------------------------------------------#
#                    THIS PYTHON SCRIPT WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                     #
#                                   EMAIL: salar.d.ghashghaei@gmail.com                                   #
###########################################################################################################
"""
Nonlinear Seismic Performance Assessment of a pedestrain bridge:
An OpenSeesPy Framework for Material and Geometric Nonlinearity Under Static, Cyclic, and Earthquake Loading

This OpenSeesPy script performs rigorous nonlinear static and dynamic analysis of a pedestrain bridge
 for performance-based earthquake engineering. The 2D model incorporates both material nonlinearity
 (elastic-perfectly plastic or hysteretic steel with strain hardening)
 and geometric nonlinearity to capture P-delta effects
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

The code continuously records element axial forces, stresses, strains,
 and full-field nodal displacements, enabling member-level demand-to-capacity
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

---------------------------

- A Vierendeel girder bridge is a unique open-frame bridge with rectangular openings and no diagonal members.
- It relies on rigid connections, so members resist loads primarily through bending (moment), not pure tension or compression.
- Without diagonals, it is structurally more complex and required more material in early designs.
- Notable examples include the Kap Shui Mun Bridge (Hong Kong), Glendale Bridges (USA), and the Bridge of the Waterhoek (Belgium).
- The design allows for unobstructed openings, useful where diagonals would be an eyesore or obstacle.
- Modern computer analysis has made designing these frames more practical.
- It is often used for pedestrian footbridges for its sleek, modern aesthetic.
- The open web also permits vehicles to pass through bridge towers, saving space. 

The Vierendeel girder is valued for both its structural and architectural advantages.
 Its open look is aesthetically pleasing, and it can be used to create large, unobstructed
 spaces where diagonal bracing would be an eyesore or functional obstacle. Modern computer
 analysis has also made designing these complex frames much more practical.
This makes it a popular choice for pedestrian footbridges, where it offers a sleek,
 modern aesthetic with functional open spans.
"""
#%%-------------------------------------------------  
# WIKIPEDIA: Vierendeel bridge
'https://en.wikipedia.org/wiki/Vierendeel_bridge'
'https://de.wikipedia.org/wiki/Gustav-Heinemann-Br%C3%BCcke_(Berlin)'
# PAPER: The Gustav-Heinemann-Bridge across the river Spree in the new centre of Berlin. Modern solutions for the dynamic problems of a slender superstructure.
'https://onlinelibrary.wiley.com/doi/10.1002/stab.200610064'    
# IMAGE OF BRIDGE AND MORE DETAILS:
'https://www.steelconstruction.info/File:Isfield_03.jpg'
# YOUTUBE: Trusses in Vierendeel Design
'https://www.youtube.com/watch?v=Q5X2i2qf5yM'
# YOUTUBE: Approximate analysis of frames (1 of 7). Vierendeel girder - single point load
'https://www.youtube.com/watch?v=XYqzSWX_NwQ'    
# YOUTUBE: 4 Vierendeel Truss
'https://www.youtube.com/watch?v=U-LnMXqgF5E&pp=ygUUVmllcmVuZGVlbCBzdHJ1Y3R1cmU%3D'    
# PAPER: Numerical simulation analysis of Steel-concrete Vierendeel Sandwich Plate Composite Bridge under a close-in blast loading
'https://www.sciencedirect.com/science/article/pii/S2352012423016508'
# PAPER: THE VIERENDEEL TRUSS: PAST AND PRESENT OF AN INNOVATIVE TYPOLOGY
'https://www.redalyc.org/journal/1936/193660402011/html/'
#%%----------------------------------------------------
import numpy as np
import openseespy.opensees as ops
import matplotlib.pyplot as plt 
import PLOT_2D as S01
import ANALYSIS_FUNCTION as S02
import PERIOD_FUN as S03
import DAMPING_RATIO_FUN as S04
import EIGENVALUE_ANALYSIS_FUN as S05
import RAYLEIGH_DAMPING_FUN as S06
import BILINEAR_CURVE as S07
import SECTION_CURVATURE_FUN as S09
import SECTION_STRESS_STRAIN_FUN as S10
import I_STEEL_PLATE_SECTION_DIFF_WIDTH_HEIGHT_FUN_EXTRA as S11
#%%----------------------------------------------------
def BRIDGE_VIERENDEEL_GIRDERS(LENGTH, HEIGHT, MAT_TYPE, ELE_TYPE, ANAL_TYPE, TOTAL_MASS):
    # Initialize model
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    
    MAX_ITERATIONS = 1000    # Maximum number of iterations
    MAX_TOLERANCE = 1.0e-8   # Specified tolerance for convergence
    
    #%% GEOMETRY DEFINITION
    NUM_POINTS = 30                  # Structural Node counts
    ε = 0.0033 * LENGTH              # [mm] Bridge deck camber (Imperfection amplitude)
    # Camber is the slight upward curvature intentionally built into a bridge deck so that the surface becomes level after the applied loads.
    dy = LENGTH / (NUM_POINTS-1)     # Element length
    for i in range(NUM_POINTS):
        x = i * dy
        y1 = ε * np.sin(np.pi * x / LENGTH)           # Bridge deck camber in y-direction
        y2 = HEIGHT + ε * np.sin(np.pi * x / LENGTH)  # Bridge deck camber in y-direction
        ops.node(i+1, x, y1)  
        ops.node(i+NUM_POINTS+1, x, y2)   
        #print(i+NUM_POINTS+1)             
    #%% SUPPORTS
    ops.fix(1, 1, 1, 0)       # left bottom - pinned
    ops.fix(30, 1, 1, 0)      # right bottom - pinned
    #ops.fix(31, 1, 1, 0)      # right top - pinned
    #ops.fix(60, 1, 1, 0)      # right top - pinned

    #%% COMPOSITE BEAM SECTION
    fc, Kc = 25, 1.25
    STEEL_DENSITY = 7850/1e9         # [kg/m^3] -> [kg/mm^3] Steel Material Density
    CONCRETE_DENSITY = 2500/1e9      # [kg/m^3] -> [kg/mm^3] Concrete Material Density

    STEEL_TYPE = 'INELASTIC'
    SEC_TAG = 1
    Depth, Ele_Mass = S11.I_STEEL_PLATE_SECTION_DIFF_WIDTH_HEIGHT_FUN_EXTRA(SEC_TAG, STEEL_TYPE, STEEL_DENSITY, plot=True)    
    #%% CREATE ELEMENTS
    ele_tag = 1
    transfTag = 1
    #ops.geomTransf('Linear', transfTag)
    #ops.geomTransf('PDelta', transfTag)
    ops.geomTransf('Corotational', transfTag)
    numIntgrPts = 5
    if ELE_TYPE == 'elasticBeamColumn':
        for J in range(NUM_POINTS-1):
            ops.element('elasticBeamColumn', J+1, J+1, J+2, SEC_TAG, transfTag,'-mass', Ele_Mass)
            ops.element('elasticBeamColumn', NUM_POINTS+J+1, NUM_POINTS+J+1, NUM_POINTS+J+2, SEC_TAG, transfTag,'-mass', Ele_Mass)
            ops.element('elasticBeamColumn', int(2*NUM_POINTS+J+1), J+1, NUM_POINTS+J+1, SEC_TAG, transfTag,'-mass', Ele_Mass)            
        ops.element('elasticBeamColumn', int(3*NUM_POINTS+1), int(NUM_POINTS), int(2*NUM_POINTS), SEC_TAG, transfTag,'-mass', Ele_Mass)
    if ELE_TYPE == 'nonlinearBeamColumn':
        for J in range(NUM_POINTS-1):
            ops.element('nonlinearBeamColumn', J+1, J+1, J+2, numIntgrPts, SEC_TAG, transfTag,'-mass', Ele_Mass) 
            ops.element('nonlinearBeamColumn', NUM_POINTS+J+1, NUM_POINTS+J+1, NUM_POINTS+J+2, numIntgrPts, SEC_TAG, transfTag,'-mass', Ele_Mass)
            ops.element('nonlinearBeamColumn', int(2*NUM_POINTS+J+1), J+1, NUM_POINTS+J+1, numIntgrPts, SEC_TAG, transfTag,'-mass', Ele_Mass) 
        ops.element('nonlinearBeamColumn', int(3*NUM_POINTS+1), int(NUM_POINTS), int(2*NUM_POINTS), numIntgrPts, SEC_TAG, transfTag,'-mass', Ele_Mass)       
    ele_tag = 30+30+1
    
    #%%--------------------------------------
    # Recorder for element forces
    """
    ops.recorder(
        'Element',
        '-file', f'ElementForces_{ANAL_TYPE}.txt',
        '-time',
        '-ele', *range(1, ele_tag),
        'localForce'
    )
    """
    """
    # Local Force Time History for Each Element 
    for i in range(1, ele_tag):
        ops.recorder('Element', '-file', f"LOCALFORCE_{ANAL_TYPE}_{i}.txt",'-time', '-ele', i,'localForce')
    # Disp. & Velocity & Acceleration, Time History for Each Node
    if ANAL_TYPE in ['DYNAMIC_EXTERNAL_TIME-DEPENDENT_LOADING', 'FREE-VIBRATION', 'SEISMIC']:
        for i in range(1, ele_tag+1):
            ops.recorder('Node', '-file', f"DISP_{ANAL_TYPE}_{i}.txt",'-time', '-node', i,'disp', '-dof', 1, 2, 3)        
            ops.recorder('Node', '-file', f"VELO_{ANAL_TYPE}_{i}.txt",'-time', '-node', i,'vel', '-dof', 1, 2, 3)
            ops.recorder('Node', '-file', f"ACCE_{ANAL_TYPE}_{i}.txt",'-time', '-node', i,'accel', '-dof', 1, 2, 3)
    """
    #%%--------------------------------------
    tag = ele_tag
    GMfact = 9810.0           # [mm/s²] standard acceleration of gravity or standard acceleration
    TOTAL_WEIGHT = TOTAL_MASS * GMfact
    ENW = TOTAL_WEIGHT / tag  # EACH NODE WEIGHT
    ENM = TOTAL_MASS / tag    # EACH NODE MASS
    #%% LOAD & MASS & ANALYSIS
    center_node = 15
    TS_TAG = 1
    PATT_TAG = 1
    ops.timeSeries('Linear', TS_TAG)
    ops.pattern('Plain', PATT_TAG, TS_TAG)
    for ii in range(1, tag):
        ops.load(ii, 0.0, -ENW, 0.0)   # [N] downward
        ops.mass(ii, ENM, ENM, 0.0)
        
    TS_TAG = 2
    PATT_TAG = 2
    ops.timeSeries('Linear', TS_TAG)
    ops.pattern('Plain', PATT_TAG, TS_TAG)    
    #ops.timeSeries('Linear', 2)
    #ops.pattern('Plain', 2, 1)    
    ops.load(center_node, 0.0, -1.0, 0.0)   # [N] downward
    
    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.system('BandGeneral')
    
    if ELE_TYPE == 'ELASTIC':
        ops.algorithm('Linear')
    if ELE_TYPE == 'INELASTIC':
        ops.algorithm('Newton')
        
    time = []
    dispX, dispY = [], []
    veloX, veloY = [], []
    accX, accY = [], []
    reaction, disp_mid = [], []
    stiffness = []
    OMEGA, PERIOD = [], []
    PERIOD_MIN, PERIOD_MAX = [], []
    STRESSb = []; STRAINb = []; 
    STRESSt = []; STRAINt = [];
    CURVATURE = [];
    
    # Initialize lists for each element's axial force, shear and moment
    num_elements = NUM_POINTS-1   # change as needed
    num_nodes = NUM_POINTS
    
    ele_axialforce = {tag: [] for tag in range(1, num_elements + 1)}
    ele_shearforce = {tag: [] for tag in range(1, num_elements + 1)}
    ele_momentforce = {tag: [] for tag in range(1, num_elements + 1)}
    # Initialize lists for each node's displacement
    node_displacementsX = {tag: [] for tag in range(1, num_nodes + 1)}
    node_displacementsY = {tag: [] for tag in range(1, num_nodes + 1)}
    if ANAL_TYPE == 'PERIOD': 
      #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(3, 0.5*DR, DR, 0, 1)
      PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(3, PLOT=True)
      # Compute modal properties
      ops.modalProperties("-print", "-file", f"SALAR_ModalReport_{ANAL_TYPE}.txt", "-unorm")
      return PERIODmin, PERIODmax
  
    if ANAL_TYPE == 'STATIC': 
        ops.integrator('LoadControl', 1.0)
        ops.analysis('Static')
        OK = ops.analyze(1)
        S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS)
        ops.reactions()
        reaction.append(ops.nodeReaction(1, 2)+ops.nodeReaction(30, 2)) # SHEAR BASE REACTION
        disp_mid.append(ops.nodeDisp(center_node, 2))                   # DISPLACEMENT NODE 15 
        dispX.append(ops.nodeDisp(center_node, 1))                      # DISPLACEMENT NODE 15 IN X DIR 
        dispY.append(ops.nodeDisp(center_node, 2))                      # DISPLACEMENT NODE 15 IN Y DIR 
        # Store axial forces, strain and stress
        for ele_id in ele_axialforce.keys(): 
            ele_axialforce[ele_id].append(ops.eleResponse(ele_id, 'force')[0])        # [N] ELEMENT AXIAL FORCE
            ele_shearforce[ele_id].append(ops.eleResponse(ele_id, 'force')[1])        # [N] ELEMENT SHEAR FORCE
            ele_momentforce[ele_id].append(ops.eleResponse(ele_id, 'force')[2])       # [N] ELEMENT MOMENT FORCE                
        # Store displacements
        for node_id in node_displacementsX.keys():    
            node_displacementsX[node_id].append(ops.nodeDisp(node_id, 1))
            node_displacementsY[node_id].append(ops.nodeDisp(node_id, 2))
        else:
            print('\n\nSTATIC ANALYSIS DONE.\n\n')  
            
        DATA = (reaction, disp_mid,
                ele_axialforce, ele_shearforce, ele_momentforce,
                node_displacementsX, node_displacementsY,
                dispX, dispY)
        
        return  DATA
    
    if ANAL_TYPE == 'PUSHOVER': # STATIC TIME-HISTORY ANALYSIS
        IDctrlDOF = 2   # 1: Horizental Dispalcement - 2: Vertical Dispalcement
        DINCR = -0.5    # [mm] Incremental Vertical Displacement
        DMAX = -1200.0   # [mm] Max. Vertical Displacement
        ops.integrator('DisplacementControl', center_node, IDctrlDOF, DINCR)
        ops.analysis('Static')
        Nsteps =  int(np.abs(DMAX/ DINCR)) 
        STEP = 0.0
        for step in range(Nsteps):
            OK = ops.analyze(1)
            S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS)
            ops.reactions()
            reaction.append(ops.nodeReaction(1, 2)+ops.nodeReaction(30, 2)) # SHEAR BASE REACTION
            disp_mid.append(ops.nodeDisp(center_node, 2))                   # DISPLACEMENT NODE 15 
            dispX.append(ops.nodeDisp(center_node, 1))                      # DISPLACEMENT NODE 15 IN X DIR 
            dispY.append(ops.nodeDisp(center_node, 2))                      # DISPLACEMENT NODE 15 IN Y DIR 
            # Store elements force, strain and stress
            for ele_id in ele_axialforce.keys(): 
                ele_axialforce[ele_id].append(ops.eleResponse(ele_id, 'force')[0])        # [N] ELEMENT AXIAL FORCE
                ele_shearforce[ele_id].append(ops.eleResponse(ele_id, 'force')[1])        # [N] ELEMENT SHEAR FORCE
                ele_momentforce[ele_id].append(ops.eleResponse(ele_id, 'force')[2])       # [N] ELEMENT MOMENT FORCE
            # Store displacements
            for node_id in node_displacementsX.keys():    
                node_displacementsX[node_id].append(ops.nodeDisp(node_id, 1))
                node_displacementsY[node_id].append(ops.nodeDisp(node_id, 2))      
            # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
            #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(3, 0.5*DR, DR, 0, 1)
            PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(3, PLOT=True)
            PERIOD_MIN.append(PERIODmin)
            PERIOD_MAX.append(PERIODmax) 
            # SECTION STRESS-STRAIN FOR bOTTOM FIBER
            stressB, strainB = S10.SECTION_STRESS_STRAIN(center_node-1, SEC_TAG, -Depth*0.5, 0.0)
            # Store in separate lists
            STRESSb.append(stressB)
            STRAINb.append(strainB)
            # SECTION STRESS-STRAIN FOR TOP FIBER
            stressT, strainT = S10.SECTION_STRESS_STRAIN(center_node-1, SEC_TAG, +Depth*0.5, 0.0)
            # Store in separate lists
            STRESSt.append(stressT)
            STRAINt.append(strainT)
            # SECTION CURVATURE
            CURVATURE.append(S09.SECTION_CURVATURE(center_node-1, SEC_TAG))
            STEP += 1
            #print(STEP, dispY[-1], reaction[-1])
            print(f"Step: {STEP}, Displacement: {dispY[-1]:.4f} mm, Reaction: {reaction[-1]:.2f} N")     
        else:
            print('\n\nPUSHOVER ANALYSIS DONE.\n\n')
            
        DATA = (reaction, disp_mid,
                ele_axialforce, ele_shearforce, ele_momentforce,
                node_displacementsX, node_displacementsY,
                dispX, dispY,
                np.array(PERIOD_MIN), np.array(PERIOD_MAX),
                STRESSb, STRAINb, STRESSt, STRAINt, CURVATURE)
    
        return  DATA
    
    if ANAL_TYPE == 'CYCLIC_DISPLACEMENT': # STATIC TIME-HISTORY ANALYSIS
        IDctrlDOF = 2   # 1: Horizental Dispalcement - 2: Vertical Dispalcement
        DMAX = -1200.0   # [mm] Max. Vertical Displacement
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
        
        # Generate 10000-point displacement protocol
        disp_protocol = np.interp(
            np.linspace(0, len(key_disp) - 1, n_points),
            np.arange(len(key_disp)),
            key_disp
        )
        
        ops.analysis('Static')
        STEP = 0.0
        for target_disp in disp_protocol:
            current_disp = ops.nodeDisp(center_node, 2) # DISPALCEMENT APPLIED IN MIDDLE NOD IN Y DIR.
            dU = target_disp - current_disp
            ops.integrator('DisplacementControl', center_node, IDctrlDOF, dU)
            OK = ops.analyze(1)
            S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS)
            ops.reactions()
            reaction.append(ops.nodeReaction(1, 2)+ops.nodeReaction(30, 2)) # SHEAR BASE REACTION
            disp_mid.append(ops.nodeDisp(center_node, 2))                   # DISPLACEMENT NODE 15 
            dispX.append(ops.nodeDisp(center_node, 1))                      # DISPLACEMENT NODE 15 IN X DIR 
            dispY.append(ops.nodeDisp(center_node, 2))                      # DISPLACEMENT NODE 15 IN Y DIR 
            # Store elements force, strain and stress
            for ele_id in ele_axialforce.keys(): 
                ele_axialforce[ele_id].append(ops.eleResponse(ele_id, 'force')[0])        # [N] ELEMENT AXIAL FORCE
                ele_shearforce[ele_id].append(ops.eleResponse(ele_id, 'force')[1])        # [N] ELEMENT SHEAR FORCE
                ele_momentforce[ele_id].append(ops.eleResponse(ele_id, 'force')[2])       # [N] ELEMENT MOMENT FORCE
            # Store displacements
            for node_id in node_displacementsX.keys():    
                node_displacementsX[node_id].append(ops.nodeDisp(node_id, 1))
                node_displacementsY[node_id].append(ops.nodeDisp(node_id, 2))      
            # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
            #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(3, 0.5*DR, DR, 0, 1)
            PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(3, PLOT=True)
            PERIOD_MIN.append(PERIODmin)
            PERIOD_MAX.append(PERIODmax) 
            # SECTION STRESS-STRAIN FOR bOTTOM FIBER
            stressB, strainB = S10.SECTION_STRESS_STRAIN(center_node-1, SEC_TAG, -Depth*0.5, 0.0)
            # Store in separate lists
            STRESSb.append(stressB)
            STRAINb.append(strainB)
            # SECTION STRESS-STRAIN FOR TOP FIBER
            stressT, strainT = S10.SECTION_STRESS_STRAIN(center_node-1, SEC_TAG, +Depth*0.5, 0.0)
            # Store in separate lists
            STRESSt.append(stressT)
            STRAINt.append(strainT)
            # SECTION CURVATURE
            CURVATURE.append(S09.SECTION_CURVATURE(center_node-1, SEC_TAG))
            STEP += 1
            #print(STEP, dispY[-1], reaction[-1])
            print(f"Step: {STEP}, Displacement: {dispY[-1]:.4f} mm, Reaction: {reaction[-1]:.2f} N")     
        else:
            print('\n\nCYCLIC DISPLAEMENT ANALYSIS DONE.\n\n')
            
        DATA = (reaction, disp_mid,
                ele_axialforce, ele_shearforce, ele_momentforce,
                node_displacementsX, node_displacementsY,
                dispX, dispY,
                np.array(PERIOD_MIN), np.array(PERIOD_MAX),
                STRESSb, STRAINb, STRESSt, STRAINt, CURVATURE)
    
        return  DATA 

    if ANAL_TYPE == 'STATIC_EXTERNAL_TIME-DEPENDENT_LOADING': # STATIC TIME-HISTORY ANALYSIS
        #%% DEFINE EXTERNAL TIME-DEPENDENT LOADING PROPERTIES
        # IN HERE ANALYSIS TIME AND DURATION ARE LOAD STEPS
        duration = 20.0             # [s] Analysis duration
        dt = 0.01                   # [s] Time step
        DT = dt                     # [s] Time step
        DT_time = 5.0               # [s] Total external Load Analysis Durations [*******]
        force_amplitude = 19600.0   # [N] Amplitude Force
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
        ops.load(center_node, 0.0, -1.0, 0.0)
        
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
            reaction.append(ops.nodeReaction(1, 2)+ops.nodeReaction(30, 2)) # SHEAR BASE REACTION
            disp_mid.append(ops.nodeDisp(center_node, 2))                   # DISPLACEMENT NODE 15 
            dispX.append(ops.nodeDisp(center_node, 1))                      # DISPLACEMENT NODE 15 IN X DIR 
            dispY.append(ops.nodeDisp(center_node, 2))                      # DISPLACEMENT NODE 15 IN Y DIR 
            # Store elements force, strain and stress
            for ele_id in ele_axialforce.keys(): 
                ele_axialforce[ele_id].append(ops.eleResponse(ele_id, 'force')[0])        # [N] ELEMENT AXIAL FORCE
                ele_shearforce[ele_id].append(ops.eleResponse(ele_id, 'force')[1])        # [N] ELEMENT SHEAR FORCE
                ele_momentforce[ele_id].append(ops.eleResponse(ele_id, 'force')[2])       # [N] ELEMENT MOMENT FORCE
            # Store displacements
            for node_id in node_displacementsX.keys():    
                node_displacementsX[node_id].append(ops.nodeDisp(node_id, 1))
                node_displacementsY[node_id].append(ops.nodeDisp(node_id, 2))      
            # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
            #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(3, 0.5*DR, DR, 0, 1)
            PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(3, PLOT=True)
            PERIOD_MIN.append(PERIODmin)
            PERIOD_MAX.append(PERIODmax) 
            # SECTION STRESS-STRAIN FOR BOTTOM FIBER
            stressB, strainB = S10.SECTION_STRESS_STRAIN(center_node-1, SEC_TAG, -Depth*0.5, 0.0)
            # Store in separate lists
            STRESSb.append(stressB)
            STRAINb.append(strainB)
            # SECTION STRESS-STRAIN FOR TOP FIBER
            stressT, strainT = S10.SECTION_STRESS_STRAIN(center_node-1, SEC_TAG, +Depth*0.5, 0.0)
            # Store in separate lists
            STRESSt.append(stressT)
            STRAINt.append(strainT)
            # SECTION CURVATURE
            CURVATURE.append(S09.SECTION_CURVATURE(center_node-1, SEC_TAG))
            STEP += 1
            #print(STEP, dispY[-1], reaction[-1])
            print(f"Step: {STEP}, Displacement: {dispY[-1]:.4f} mm, Reaction: {reaction[-1]:.2f} N")     
        else:
            print('\n\nSTATIC EXTERNAL TIME-DEPENDENT LOADING ANALYSIS DONE.\n\n')
            
        DATA = (reaction, disp_mid,
                ele_axialforce, ele_shearforce, ele_momentforce,
                node_displacementsX, node_displacementsY,
                dispX, dispY,
                np.array(PERIOD_MIN), np.array(PERIOD_MAX),
                STRESSb, STRAINb, STRESSt, STRAINt, CURVATURE)
    
        return  DATA
    
    if ANAL_TYPE == 'DYNAMIC_EXTERNAL_TIME-DEPENDENT_LOADING': # DYNAMIC TIME-HISTORY ANALYSIS
        #%% DEFINE EXTERNAL TIME-DEPENDENT LOADING PROPERTIES
        duration = 20.0             # [s] Analysis duration
        dt = 0.001                  # [s] Time step
        DT = dt                     # [s] Time step
        DT_time = 5.0               # [s] Total external Load Analysis Durations [*******]
        force_amplitude = 19600.0   # [N] Amplitude Force
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
        ops.load(center_node, 0.0, -1.0, 0.0)
        
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
            reaction.append(ops.nodeReaction(1, 2)+ops.nodeReaction(30, 2)) # SHEAR BASE REACTION
            disp_mid.append(ops.nodeDisp(center_node, 2))                   # DISPLACEMENT NODE 15 
            dispX.append(ops.nodeDisp(center_node, 1))                      # DISPLACEMENT NODE 15 IN X DIR 
            dispY.append(ops.nodeDisp(center_node, 2))                      # DISPLACEMENT NODE 15 IN Y DIR 
            veloX.append(ops.nodeVel(center_node, 1))                       # VELOCITY NODE 15
            veloY.append(ops.nodeVel(center_node, 2))                       # VELOCITY NODE 15
            accX.append(ops.nodeAccel(center_node, 1))                      # ACCELERATION NODE 15
            accY.append(ops.nodeAccel(center_node, 2))                      # ACCELERATION NODE 15
            stiffness.append(np.abs(reaction[-1] / disp_mid[-1]))
            OMEGA.append(np.sqrt(stiffness[-1]/TOTAL_MASS))
            PERIOD.append((np.pi * 2) / OMEGA[-1])
            # Store elements force, strain and stress
            for ele_id in ele_axialforce.keys(): 
                ele_axialforce[ele_id].append(ops.eleResponse(ele_id, 'force')[0])        # [N] ELEMENT AXIAL FORCE
                ele_shearforce[ele_id].append(ops.eleResponse(ele_id, 'force')[1])        # [N] ELEMENT SHEAR FORCE
                ele_momentforce[ele_id].append(ops.eleResponse(ele_id, 'force')[2])       # [N] ELEMENT MOMENT FORCE
            # Store displacements
            for node_id in node_displacementsX.keys():    
                node_displacementsX[node_id].append(ops.nodeDisp(node_id, 1))
                node_displacementsY[node_id].append(ops.nodeDisp(node_id, 2))   
            # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
            #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(3, 0.5*DR, DR, 0, 1)
            PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(3, PLOT=True)
            PERIOD_MIN.append(PERIODmin)
            PERIOD_MAX.append(PERIODmax)
            # SECTION STRESS-STRAIN FOR bOTTOM FIBER
            stressB, strainB = S10.SECTION_STRESS_STRAIN(center_node-1, SEC_TAG, -Depth*0.5, 0.0)
            # Store in separate lists
            STRESSb.append(stressB)
            STRAINb.append(strainB)
            # SECTION STRESS-STRAIN FOR TOP FIBER
            stressT, strainT = S10.SECTION_STRESS_STRAIN(center_node-1, SEC_TAG, +Depth*0.5, 0.0)
            # Store in separate lists
            STRESSt.append(stressT)
            STRAINt.append(strainT)
            # SECTION CURVATURE
            CURVATURE.append(S09.SECTION_CURVATURE(center_node-1, SEC_TAG))
            #print(time[-1], dispY[-1], veloY[-1])
            print(f"Time: {time[-1]:.4f}, Displacement: {dispY[-1]:.4f} mm, Reaction: {reaction[-1]:.2f} N")      
        else:
            print('\n\nDYNAMIC EXTERNAL TIME-DEPENDENT LOADING ANALYSIS DONE.\n\n')  
        # Calculating Damping Ratio and Period Using Logarithmic Decrement Analysis 
        damping_ratio = S04.DAMPING_RATIO(dispY)  
        
        # Compute modal properties
        ops.modalProperties("-print", "-file", f"SALAR_ModalReport_{ANAL_TYPE}.txt", "-unorm") 
        
        DATA = (time, reaction, disp_mid,
                ele_axialforce, ele_shearforce, ele_momentforce,
                node_displacementsX, node_displacementsY,
                dispX, dispY,
                veloX, veloY,
                accX, accY,
                stiffness, PERIOD, damping_ratio,
                np.array(PERIOD_MIN), np.array(PERIOD_MAX),
                STRESSb, STRAINb, STRESSt, STRAINt, CURVATURE)
        
        return  DATA
        
    if ANAL_TYPE == 'FREE-VIBRATION': # DYNAMIC TIME-HISTORY ANALYSIS
        #%% DEFINE PARAMETERS FOR FREE-VIBRATION ANALYSIS
        u0 = -0.00015                      # [mm] Initial displacement
        v0 = 5.15                          # [mm/s] Initial velocity
        a0 = 0.065                         # [mm/s^2] Initial acceleration
        IU = True                          # Free Vibration with Initial Displacement
        IV = True                          # Free Vibration with Initial Velocity
        IA = True                          # Free Vibration with Initial Acceleration
        duration = 20.0                    # [s] Analysis duration
        dt = 0.001                         # [s] Time step
        DR = 0.03                          # Damping Ratio
        
        if IU == True:
            # Define initial displacment
            ops.setNodeDisp(center_node, 2, u0, '-commit') # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/setNodeDisp.html
        if IV == True:
            # Define initial velocity
            ops.setNodeVel(center_node, 2, v0, '-commit')  # INFO LINK: https://openseespydoc.readthedocs.io/en/stable/src/setNodeVel.html
        if IA == True:
            # Define initial  acceleration
            ops.setNodeAccel(center_node, 2, a0, '-commit') # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/setNodeAccel.html
            
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
            reaction.append(ops.nodeReaction(1, 2)+ops.nodeReaction(30, 2)) # SHEAR BASE REACTION
            disp_mid.append(ops.nodeDisp(center_node, 2))                   # DISPLACEMENT NODE 15 
            dispX.append(ops.nodeDisp(center_node, 1))                      # DISPLACEMENT NODE 15 IN X DIR 
            dispY.append(ops.nodeDisp(center_node, 2))                      # DISPLACEMENT NODE 15 IN Y DIR 
            veloX.append(ops.nodeVel(center_node, 1))                       # VELOCITY NODE 15
            veloY.append(ops.nodeVel(center_node, 2))                       # VELOCITY NODE 15
            accX.append(ops.nodeAccel(center_node, 1))                      # ACCELERATION NODE 15
            accY.append(ops.nodeAccel(center_node, 2))                      # ACCELERATION NODE 15
            stiffness.append(np.abs(reaction[-1] / disp_mid[-1]))
            OMEGA.append(np.sqrt(stiffness[-1]/TOTAL_MASS))
            PERIOD.append((np.pi * 2) / OMEGA[-1])
            # Store elements force, strain and stress
            for ele_id in ele_axialforce.keys(): 
                ele_axialforce[ele_id].append(ops.eleResponse(ele_id, 'force')[0])        # [N] ELEMENT AXIAL FORCE
                ele_shearforce[ele_id].append(ops.eleResponse(ele_id, 'force')[1])        # [N] ELEMENT SHEAR FORCE
                ele_momentforce[ele_id].append(ops.eleResponse(ele_id, 'force')[2])       # [N] ELEMENT MOMENT FORCE
            # Store displacements
            for node_id in node_displacementsX.keys():    
                node_displacementsX[node_id].append(ops.nodeDisp(node_id, 1))
                node_displacementsY[node_id].append(ops.nodeDisp(node_id, 2))   
            # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
            #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(3, 0.5*DR, DR, 0, 1)
            PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(3, PLOT=True)
            PERIOD_MIN.append(PERIODmin)
            PERIOD_MAX.append(PERIODmax)
            # SECTION STRESS-STRAIN FOR BOTTOM FIBER
            stressB, strainB = S10.SECTION_STRESS_STRAIN(center_node-1, SEC_TAG, -Depth*0.5, 0.0)
            # Store in separate lists
            STRESSb.append(stressB)
            STRAINb.append(strainB)
            # SECTION STRESS-STRAIN FOR TOP FIBER
            stressT, strainT = S10.SECTION_STRESS_STRAIN(center_node-1, SEC_TAG, +Depth*0.5, 0.0)
            # Store in separate lists
            STRESSt.append(stressT)
            STRAINt.append(strainT)
            # SECTION CURVATURE
            CURVATURE.append(S09.SECTION_CURVATURE(center_node-1, SEC_TAG))
            #print(time[-1], dispY[-1], veloY[-1])
            print(f"Time: {time[-1]:.4f}, Displacement: {dispY[-1]:.4f} mm, Reaction: {reaction[-1]:.2f} N")      
        else:
            print('\n\nFREE-VIBRATION ANALYSIS DONE.\n\n')    
        # Calculating Damping Ratio and Period Using Logarithmic Decrement Analysis 
        damping_ratio = S04.DAMPING_RATIO(dispY)  
        
        # Compute modal properties
        ops.modalProperties("-print", "-file", f"SALAR_ModalReport_{ANAL_TYPE}.txt", "-unorm") 
        
        DATA = (time, reaction, disp_mid,
                ele_axialforce, ele_shearforce, ele_momentforce,
                node_displacementsX, node_displacementsY,
                dispX, dispY,
                veloX, veloY,
                accX, accY,
                stiffness, PERIOD, damping_ratio,
                np.array(PERIOD_MIN), np.array(PERIOD_MAX),
                STRESSb, STRAINb, STRESSt, STRAINt, CURVATURE)
        
        return DATA
    
    if ANAL_TYPE == 'SEISMIC': # DYNAMIC TIME-HISTORY ANALYSIS
        #%% DEFINE PARAMETERS FOR SEISMIC ANALYSIS
        duration = 5.0                     # [s] Analysis duration
        dt = 0.01                          # [s] Time step
        GMfact = 9810                      # [mm/s²]standard acceleration of gravity or standard acceleration
        SSF_X = 0.5                        # Seismic Acceleration Scale Factor in X Direction
        SSF_Y = 0.5                        # Seismic Acceleration Scale Factor in Y Direction
        iv0_X = 0.005                      # [mm/s] Initial velocity applied to the node  in X Direction
        iv0_Y = 0.005                      # [mm/s] Initial velocity applied to the node  in Y Direction
        st_iv0 = 0.0                       # [s] Initial velocity applied starting time
        SEI = 'Y'                          # Seismic Direction
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
            reaction.append(ops.nodeReaction(1, 2)+ops.nodeReaction(30, 2)) # SHEAR BASE REACTION
            disp_mid.append(ops.nodeDisp(center_node, 2))                   # DISPLACEMENT NODE 15 
            dispX.append(ops.nodeDisp(center_node, 1))                      # DISPLACEMENT NODE 15 IN X DIR 
            dispY.append(ops.nodeDisp(center_node, 2))                      # DISPLACEMENT NODE 15 IN Y DIR 
            veloX.append(ops.nodeVel(center_node, 1))                       # VELOCITY NODE 15
            veloY.append(ops.nodeVel(center_node, 2))                       # VELOCITY NODE 15
            accX.append(ops.nodeAccel(center_node, 1))                      # ACCELERATION NODE 15
            accY.append(ops.nodeAccel(center_node, 2))                      # ACCELERATION NODE 15
            stiffness.append(np.abs(reaction[-1] / disp_mid[-1]))
            OMEGA.append(np.sqrt(stiffness[-1]/TOTAL_MASS))
            PERIOD.append((np.pi * 2) / OMEGA[-1])
            # Store elements force, strain and stress
            for ele_id in ele_axialforce.keys(): 
                ele_axialforce[ele_id].append(ops.eleResponse(ele_id, 'force')[0])        # [N] ELEMENT AXIAL FORCE
                ele_shearforce[ele_id].append(ops.eleResponse(ele_id, 'force')[1])        # [N] ELEMENT SHEAR FORCE
                ele_momentforce[ele_id].append(ops.eleResponse(ele_id, 'force')[2])       # [N] ELEMENT MOMENT FORCE
            # Store displacements
            for node_id in node_displacementsX.keys():    
                node_displacementsX[node_id].append(ops.nodeDisp(node_id, 1))
                node_displacementsY[node_id].append(ops.nodeDisp(node_id, 2))
            # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
            #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(3, 0.5*DR, DR, 0, 1)
            PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(3, PLOT=True)
            PERIOD_MIN.append(PERIODmin)
            PERIOD_MAX.append(PERIODmax)
            # SECTION STRESS-STRAIN FOR BOTTOM FIBER
            stressB, strainB = S10.SECTION_STRESS_STRAIN(center_node-1, SEC_TAG, -Depth*0.5, 0.0)
            # Store in separate lists
            STRESSb.append(stressB)
            STRAINb.append(strainB)
            # SECTION STRESS-STRAIN FOR TOP FIBER
            stressT, strainT = S10.SECTION_STRESS_STRAIN(center_node-1, SEC_TAG, +Depth*0.5, 0.0)
            # Store in separate lists
            STRESSt.append(stressT)
            STRAINt.append(strainT)
            # SECTION CURVATURE
            CURVATURE.append(S09.SECTION_CURVATURE(center_node-1, SEC_TAG))
            #print(time[-1], dispY[-1], veloY[-1])
            print(f"Time: {time[-1]:.4f}, Displacement: {dispY[-1]:.4f} mm, Reaction: {reaction[-1]:.2f} N")

        else:
            print('\n\nSEISMIC ANALYSIS DONE.\n\n')     
        # Calculating Damping Ratio and Period Using Logarithmic Decrement Analysis 
        damping_ratio = S04.DAMPING_RATIO(dispY)  
        
        # Compute modal properties
        ops.modalProperties("-print", "-file", f"SALAR_ModalReport_{ANAL_TYPE}.txt", "-unorm") 
        
        DATA = (time, reaction, disp_mid,
                ele_axialforce, ele_shearforce, ele_momentforce,
                node_displacementsX, node_displacementsY,
                dispX, dispY,
                veloX, veloY,
                accX, accY,
                stiffness, PERIOD, damping_ratio,
                np.array(PERIOD_MIN), np.array(PERIOD_MAX),
                STRESSb, STRAINb, STRESSt, STRAINt, CURVATURE)
        
        return DATA    

#%%-------------------------------------------------------
def PLOT_TIME_HISTORY(time,
                      reaction, disp_mid,
                      dispX, dispY,
                      veloX, veloY,
                      accX, accY):

    import matplotlib.pyplot as plt

    data_dict = {
        "Reaction [N]": reaction,
        "Mid Displacement [mm]": disp_mid,
        "Disp X [mm]": dispX,
        "Disp Y [mm]": dispY,
        "Velocity X [mm/s]": veloX,
        "Velocity Y [mm/s]": veloY,
        "Acceleration X [mm/s²]": accX,
        "Acceleration Y [mm/s²]": accY,
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
    plt.ylabel('Displacement [mm]')
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
LENGTH = 66000.0      # [mm] Bridge Horizental Length
HEIGHT = 1600.0       # [mm] Bridge Height
TOTAL_MASS = 0.0      # [kg] Total Mass of Structure
#%%----------------------------------------------------
# PERIOD ANALYSIS
#ELE_TYPE = 'elasticBeamColumn'      # MATERIAL LINEARITY
ELE_TYPE = 'nonlinearBeamColumn'     # MATERIAL AND GEOMETRIC NONLINEARITIES
MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'
ANAL_TYPE = 'PERIOD'

DATA = BRIDGE_VIERENDEEL_GIRDERS(LENGTH, HEIGHT, MAT_TYPE, ELE_TYPE, ANAL_TYPE, TOTAL_MASS)
(PERIOD_MIN_X, PERIOD_MAX_X) = DATA
print('Structure First Period:  ', PERIOD_MIN_X)
print('Structure Second Period: ', PERIOD_MAX_X) 

S01.PLOT_2D_FRAME(deformed_scale=1.0)
#%%----------------------------------------------------
# STATIC ANALYSIS
#ELE_TYPE = 'elasticBeamColumn'      # MATERIAL LINEARITY
ELE_TYPE = 'nonlinearBeamColumn'     # MATERIAL AND GEOMETRIC NONLINEARITIES
MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'
ANAL_TYPE = 'STATIC'

DATA = BRIDGE_VIERENDEEL_GIRDERS(LENGTH, HEIGHT, MAT_TYPE, ELE_TYPE, ANAL_TYPE, TOTAL_MASS)
(reaction, disp_mid, 
 ele_axialforce, ele_shearforce, ele_momentforce,
 node_displacementsX, node_displacementsY,
 dispX, dispY) = DATA

S01.PLOT_2D_FRAME(deformed_scale=1.0)
#%%----------------------------------------------------
# PUSHOVER ANALYSIS (STATIC TIME-HISTORY ANALYSIS)
#ELE_TYPE = 'elasticBeamColumn'      # MATERIAL LINEARITY
ELE_TYPE = 'nonlinearBeamColumn'     # MATERIAL AND GEOMETRIC NONLINEARITIES
MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'
ANAL_TYPE = 'PUSHOVER'

DATA = BRIDGE_VIERENDEEL_GIRDERS(LENGTH, HEIGHT, MAT_TYPE, ELE_TYPE, ANAL_TYPE, TOTAL_MASS)
(reaction_PUSH, disp_mid_PUSH, 
 ele_axialforce_PUSH, ele_shearforce_PUSH, ele_momentforce_PUSH,
 node_displacementsX_PUSH, node_displacementsY_PUSH,
  dispX_PUSH, dispY_PUSH,
 PERIOD_MIN_PUSH, PERIOD_MAX_PUSH,
 STRESSb_PUSH, STRAINb_PUSH, STRESSt_PUSH, STRAINt_PUSH, CURVATURE_PUSH) = DATA


XDATA = disp_mid_PUSH
YDATA = reaction_PUSH
XLABEL = 'Displacement in Middle Span [mm]'
YLABEL = 'Base Reaction [N]'
TITLE = 'Base Reaction and Displacement of Structure During Pushover Analysis'
COLOR = 'black'
SEMILOGY = False
PLOT(XDATA, YDATA, TITLE, XLABEL, YLABEL, COLOR, SEMILOGY)

DATA = S07.BILNEAR_CURVE(np.abs(disp_mid_PUSH), np.abs(reaction_PUSH), SLOPE_NODE=10)
(X_PUSH, Y_PUSH, Elastic_ST, Plastic_ST, Tangent_ST, Ductility_Rito, Over_Strength_Factor) = DATA

# PLOT STRUCTURAL PERIOD DURING THE ANALYSIS
plt.figure(0, figsize=(12, 8))
plt.plot(disp_mid_PUSH, PERIOD_MIN_PUSH, linewidth=3)
plt.plot(disp_mid_PUSH, PERIOD_MAX_PUSH, linewidth=3)
plt.title('Period of Structure During Pushover Analysis')
plt.ylabel('Structural Period [s]')
plt.xlabel('Displacement [mm]')
#plt.semilogy()
plt.grid()
plt.legend([f'PERIOD - MIN VALUES: Min: {np.min(PERIOD_MIN_PUSH):.3f} (s) - Mean: {np.mean(PERIOD_MIN_PUSH):.3f} (s) - Max: {np.max(PERIOD_MIN_PUSH):.3f} (s)', 
            f'PERIOD - MAX VALUES:  Min: {np.min(PERIOD_MAX_PUSH):.3f} (s) - Mean: {np.mean(PERIOD_MAX_PUSH):.3f} (s) - Max: {np.max(PERIOD_MAX_PUSH):.3f} (s)',
            ])
plt.show()

# Plot Stress-strain of top and bottom fibers
plt.figure(44, figsize=(12, 8))
plt.plot(STRAINb_PUSH, STRESSb_PUSH, linewidth=4, label=f"BOT FIBER - Stress: {np.abs(STRESSb_PUSH[-1]) : .3f}")
plt.plot(STRAINt_PUSH, STRESSt_PUSH, linewidth=4, label=f"TOP FIBER - Stress: {np.abs(STRESSt_PUSH[-1]) : .3f}")
plt.xlabel('Strain (mm/mm)')
plt.ylabel('Stress (N/mm^2)')
plt.title('Behavior of beam - Stress-strain of Top and Bottom Fibers Curve')
plt.grid(True)
plt.legend()
plt.show()

S01.PLOT_2D_FRAME(deformed_scale=1.0)
#%%----------------------------------------------------
# CYCLIC DISPLACEMENT ANALYSIS (STATIC TIME-HISTORY ANALYSIS)
#ELE_TYPE = 'elasticBeamColumn'      # MATERIAL LINEARITY
ELE_TYPE = 'nonlinearBeamColumn'     # MATERIAL AND GEOMETRIC NONLINEARITIES
MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'
ANAL_TYPE = 'CYCLIC_DISPLACEMENT'

DATA = BRIDGE_VIERENDEEL_GIRDERS(LENGTH, HEIGHT, MAT_TYPE, ELE_TYPE, ANAL_TYPE, TOTAL_MASS)
(reaction_CP, disp_mid_CP, 
 ele_axialforce_CP, ele_shearforce_CP, ele_momentforce_CP,
 node_displacementsX_CP, node_displacementsY_CP,
  dispX_CP, dispY_CP,
 PERIOD_MIN_CP, PERIOD_MAX_CP,
 STRESSb_CP, STRAINb_CP, STRESSt_CP, STRAINt_CP, CURVATURE_CP) = DATA


XDATA = disp_mid_CP
YDATA = reaction_CP
XLABEL = 'Displacement in Middle Span [mm]'
YLABEL = 'Base Reaction [N]'
TITLE = 'Base Reaction and Displacement of Structure During Cyclic-Displacement Analysis'
COLOR = 'black'
SEMILOGY = False
PLOT(XDATA, YDATA, TITLE, XLABEL, YLABEL, COLOR, SEMILOGY)

# PLOT STRUCTURAL PERIOD DURING THE ANALYSIS
plt.figure(0, figsize=(12, 8))
plt.plot(disp_mid_CP, PERIOD_MIN_CP, linewidth=3)
plt.plot(disp_mid_CP, PERIOD_MAX_CP, linewidth=3)
plt.title('Period of Structure During Cyclic Displacement Analysis')
plt.ylabel('Structural Period [s]')
plt.xlabel('Displacement [mm]')
#plt.semilogy()
plt.grid()
plt.legend([f'PERIOD - MIN VALUES: Min: {np.min(PERIOD_MIN_CP):.3f} (s) - Mean: {np.mean(PERIOD_MIN_CP):.3f} (s) - Max: {np.max(PERIOD_MIN_CP):.3f} (s)', 
            f'PERIOD - MAX VALUES:  Min: {np.min(PERIOD_MAX_CP):.3f} (s) - Mean: {np.mean(PERIOD_MAX_CP):.3f} (s) - Max: {np.max(PERIOD_MAX_CP):.3f} (s)',
            ])
plt.show()

# Plot Stress-strain of top and bottom fibers
plt.figure(44, figsize=(12, 8))
plt.plot(STRAINb_CP, STRESSb_CP, linewidth=4, label=f"BOT FIBER - Stress: {np.abs(STRESSb_CP[-1]) : .3f}")
plt.plot(STRAINt_CP, STRESSt_CP, linewidth=4, label=f"TOP FIBER - Stress: {np.abs(STRESSt_CP[-1]) : .3f}")
plt.xlabel('Strain (mm/mm)')
plt.ylabel('Stress (N/mm^2)')
plt.title('Behavior of beam - Stress-strain of Top and Bottom Fibers Curve')
plt.grid(True)
plt.legend()
plt.show()

S01.PLOT_2D_FRAME(deformed_scale=1.0)
#%%----------------------------------------------------
# EXTERNAL TIME-DEPENDENT LOADING ANALYSIS (STATIC TIME-HISTORY ANALYSIS)
#ELE_TYPE = 'elasticBeamColumn'      # MATERIAL LINEARITY
ELE_TYPE = 'nonlinearBeamColumn'     # MATERIAL AND GEOMETRIC NONLINEARITIES
MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'
ANAL_TYPE = 'STATIC_EXTERNAL_TIME-DEPENDENT_LOADING'

DATA = BRIDGE_VIERENDEEL_GIRDERS(LENGTH, HEIGHT, MAT_TYPE, ELE_TYPE, ANAL_TYPE, TOTAL_MASS)

(reaction_ETDLS, disp_mid_ETDLS, 
 ele_axialforce_ETDLS, ele_shearforce_ETDLS, ele_momentforce_ETDLS,
 node_displacementsX_ETDLS, node_displacementsY_ETDLS,
  dispX_ETDLS, dispY_ETDLS,
 PERIOD_MIN_ETDLS, PERIOD_MAX_ETDLS,
 STRESSb_ETDLS, STRAINb_ETDLS, STRESSt_ETDLS, STRAINt_ETDLS, CURVATURE_ETDLS) = DATA


XDATA = disp_mid_ETDLS
YDATA = reaction_ETDLS
XLABEL = 'Displacement in Middle Span [mm]'
YLABEL = 'Base Reaction [N]'
TITLE = 'Base Reaction and Displacement of Structure During Static External Time-dependent Loading Analysis'
COLOR = 'black'
SEMILOGY = False
PLOT(XDATA, YDATA, TITLE, XLABEL, YLABEL, COLOR, SEMILOGY)

# PLOT STRUCTURAL PERIOD DURING THE ANALYSIS
plt.figure(0, figsize=(12, 8))
plt.plot(disp_mid_ETDLS, PERIOD_MIN_ETDLS, linewidth=3)
plt.plot(disp_mid_ETDLS, PERIOD_MAX_ETDLS, linewidth=3)
plt.title('Period of Structure During Static External Time-dependent Loading Analysis')
plt.ylabel('Structural Period [s]')
plt.xlabel('Displacement [mm]')
#plt.semilogy()
plt.grid()
plt.legend([f'PERIOD - MIN VALUES: Min: {np.min(PERIOD_MIN_ETDLS):.3f} (s) - Mean: {np.mean(PERIOD_MIN_ETDLS):.3f} (s) - Max: {np.max(PERIOD_MIN_ETDLS):.3f} (s)', 
            f'PERIOD - MAX VALUES:  Min: {np.min(PERIOD_MAX_ETDLS):.3f} (s) - Mean: {np.mean(PERIOD_MAX_ETDLS):.3f} (s) - Max: {np.max(PERIOD_MAX_ETDLS):.3f} (s)',
            ])
plt.show()

# Plot Stress-strain of top and bottom fibers
plt.figure(44, figsize=(12, 8))
plt.plot(STRAINb_ETDLS, STRESSb_ETDLS, linewidth=4, label=f"BOT FIBER - Stress: {np.abs(STRESSb_ETDLS[-1]) : .3f}")
plt.plot(STRAINt_ETDLS, STRESSt_ETDLS, linewidth=4, label=f"TOP FIBER - Stress: {np.abs(STRESSt_ETDLS[-1]) : .3f}")
plt.xlabel('Strain (mm/mm)')
plt.ylabel('Stress (N/mm^2)')
plt.title('Behavior of beam - Stress-strain of Top and Bottom Fibers Curve')
plt.grid(True)
plt.legend()
plt.show()

S01.PLOT_2D_FRAME(deformed_scale=1.0)
#%%----------------------------------------------------
# EXTERNAL TIME-DEPENDENT LOADING ANALYSIS (DYNAMIC TIME-HISTORY ANALYSIS)
#ELE_TYPE = 'elasticBeamColumn'      # MATERIAL LINEARITY
ELE_TYPE = 'nonlinearBeamColumn'     # MATERIAL AND GEOMETRIC NONLINEARITIES
MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'
ANAL_TYPE = 'DYNAMIC_EXTERNAL_TIME-DEPENDENT_LOADING'

DATA = BRIDGE_VIERENDEEL_GIRDERS(LENGTH, HEIGHT, MAT_TYPE, ELE_TYPE, ANAL_TYPE, TOTAL_MASS)

(time_ETDLD, reaction_ETDLD, disp_mid_ETDLD,
ele_axialforce_ETDLD, ele_shearforce_ETDLD, ele_momentforce_ETDLD,
node_displacementsX_ETDLD, node_displacementsY_ETDLD,
dispX_ETDLD, dispY_ETDLD,
veloX_ETDLD, veloY_ETDLD,
accX_ETDLD, accY_ETDLD,
stiffness_ETDLD, PERIOD_ETDLD, damping_ratio_ETDLD,
PERIOD_MIN_ETDLD, PERIOD_MAX_ETDLD,
STRESSb_ETDLD, STRAINb_ETDLD, STRESSt_ETDLD, STRAINt_ETDLD, CURVATURE_ETDLD) = DATA


XDATA = disp_mid_ETDLD
YDATA = reaction_ETDLD
XLABEL = 'Displacement in Middle Span [mm]'
YLABEL = 'Base Reaction [N]'
TITLE = 'Base Reaction and Displacement of Structure During Dynamic External Time-dependent Loading Analysis'
COLOR = 'black'
SEMILOGY = False
PLOT(XDATA, YDATA, TITLE, XLABEL, YLABEL, COLOR, SEMILOGY)

PLOT_TIME_HISTORY(time_ETDLD, reaction_ETDLD, disp_mid_ETDLD,
                          dispX_ETDLD, dispY_ETDLD,
                          veloX_ETDLD, veloY_ETDLD,
                          accX_ETDLD, accY_ETDLD)

# PLOT ELEMENTS AXIAL FORCE
YLABEL = 'Element Axial Force [N]'  
TITLE = "elements  force in X Dir. vs Time for Beam Element During Dynamic External Time-dependent Loading Analysis"
PLOT_FORCES(time_ETDLD, ele_axialforce_ETDLD, YLABEL, TITLE) # Beam Element - ELEMENTS AXIAL FORCE
# PLOT ELEMENTS SHEAR FORCE
YLABEL = 'Element Shear Force [N]'  
TITLE = "elements  force in Y Dir. vs Time for Beam Element During Dynamic External Time-dependent Loading Analysis"
PLOT_FORCES(time_ETDLD, ele_shearforce_ETDLD, YLABEL, TITLE) # Beam Element - ELEMENTS SHEAR FORCE
# PLOT ELEMENTS MOMENT FORCE
YLABEL = 'Element Moment Force [N.mm]'  
TITLE = "elements  force in Z Dir. vs Time for Beam Element During Dynamic External Time-dependent Loading Analysis"
PLOT_FORCES(time_ETDLD, ele_momentforce_ETDLD, YLABEL, TITLE) # Beam Element - ELEMENTS MOMENT FORCE

# PLOT NODAL DISPALEMENTS
PLOT_DISPLAEMENTS(time_ETDLD, node_displacementsX_ETDLD, TITLE = "Node Displacements in X Dir. vs Time for Beam Element During Dynamic External Time-dependent Loading Analysis") # Beam Element IN X DIR.
PLOT_DISPLAEMENTS(time_ETDLD, node_displacementsY_ETDLD, TITLE = "Node Displacements in Y Dir. vs Time for Beam Element During Dynamic External Time-dependent Loading Analysis") # Beam Element IN Y DIR.


# PLOT STRUCTURAL PERIOD DURING THE ANALYSIS
plt.figure(0, figsize=(12, 8))
plt.plot(disp_mid_ETDLD, PERIOD_MIN_ETDLD, linewidth=3)
plt.plot(disp_mid_ETDLD, PERIOD_MAX_ETDLD, linewidth=3)
plt.title('Period of Structure During Dynamic External Time-dependent Loading Analysis')
plt.ylabel('Structural Period [s]')
plt.xlabel('Displacement [mm]')
#plt.semilogy()
plt.grid()
plt.legend([f'PERIOD - MIN VALUES: Min: {np.min(PERIOD_MIN_ETDLD):.3f} (s) - Mean: {np.mean(PERIOD_MIN_ETDLD):.3f} (s) - Max: {np.max(PERIOD_MIN_ETDLD):.3f} (s)', 
            f'PERIOD - MAX VALUES:  Min: {np.min(PERIOD_MAX_ETDLD):.3f} (s) - Mean: {np.mean(PERIOD_MAX_ETDLD):.3f} (s) - Max: {np.max(PERIOD_MAX_ETDLD):.3f} (s)',
            ])
plt.show()

# Plot Stress-strain of top and bottom fibers
plt.figure(44, figsize=(12, 8))
plt.plot(STRAINb_ETDLD, STRESSb_ETDLD, linewidth=4, label=f"BOT FIBER - Stress: {np.abs(STRESSb_ETDLD[-1]) : .3f}")
plt.plot(STRAINt_ETDLD, STRESSt_ETDLD, linewidth=4, label=f"TOP FIBER - Stress: {np.abs(STRESSt_ETDLD[-1]) : .3f}")
plt.xlabel('Strain (mm/mm)')
plt.ylabel('Stress (N/mm^2)')
plt.title('Behavior of beam - Stress-strain of Top and Bottom Fibers Curve')
plt.grid(True)
plt.legend()
plt.show()

S01.PLOT_2D_FRAME(deformed_scale=1.0)
#%%----------------------------------------------------
# FREE-VIBRATION ANALYSIS (DYNAMIC TIME-HISTORY ANALYSIS)
#ELE_TYPE = 'elasticBeamColumn'      # MATERIAL LINEARITY
ELE_TYPE = 'nonlinearBeamColumn'     # MATERIAL AND GEOMETRIC NONLINEARITIES
MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'
ANAL_TYPE = 'FREE-VIBRATION'
DATA = BRIDGE_VIERENDEEL_GIRDERS(LENGTH, HEIGHT, MAT_TYPE, ELE_TYPE, ANAL_TYPE, TOTAL_MASS)

(time_FV, reaction_FV, disp_mid_FV,
ele_axialforce_FV, ele_shearforce_FV, ele_momentforce_FV,
node_displacementsX_FV, node_displacementsY_FV,
dispX_FV, dispY_FV,
veloX_FV, veloY_FV,
accX_FV, accY_FV,
stiffness_FV, PERIOD_FV, damping_ratio_FV,
PERIOD_MIN_FV, PERIOD_MAX_FV,
STRESSb_FV, STRAINb_FV, STRESSt_FV, STRAINt_FV, CURVATURE_FV) = DATA

XDATA = disp_mid_FV
YDATA = reaction_FV
XLABEL = 'Displacement in Middle Span [mm]'
YLABEL = 'Base Reaction [N]'
TITLE = 'Base Reaction and Dispalcement of Structure During Free-vibration Analysis'
COLOR = 'black'
SEMILOGY = False
PLOT(XDATA, YDATA, TITLE, XLABEL, YLABEL, COLOR, SEMILOGY)

PLOT_TIME_HISTORY(time_FV, reaction_FV, disp_mid_FV,
                          dispX_FV, dispY_FV,
                          veloX_FV, veloY_FV,
                          accX_FV, accY_FV)

# PLOT ELEMENTS AXIAL FORCE
YLABEL = 'Element Axial Force [N]'
TITLE = "elements  force in X Dir. vs Time for Beam Element During Free-vibration Analysis"
PLOT_FORCES(time_FV, ele_axialforce_FV, YLABEL, TITLE) # Beam Element - ELEMENTS AXIAL FORCE
# PLOT ELEMENTS SHEAR FORCE
YLABEL = 'Element Shear Force [N]'  
TITLE = "elements  force in Y Dir. vs Time for Beam Element During Free-vibration Analysis"
PLOT_FORCES(time_FV, ele_shearforce_FV, YLABEL, TITLE) # Beam Element - ELEMENTS SHEAR FORCE
# PLOT ELEMENTS MOMENT FORCE
YLABEL = 'Element Moment Force [N.mm]'  
TITLE = "elements  force in Z Dir. vs Time for Beam Element During Free-vibration Analysis"
PLOT_FORCES(time_FV, ele_momentforce_FV, YLABEL, TITLE) # Beam Element - ELEMENTS MOMENT FORCE

# PLOT NODAL DISPALEMENTS
PLOT_DISPLAEMENTS(time_FV, node_displacementsX_FV, TITLE = "Node Displacements in X Dir. vs Time for Beam Element During Free-vibration Analysis") # Beam Element IN X DIR.
PLOT_DISPLAEMENTS(time_FV, node_displacementsY_FV, TITLE = "Node Displacements in Y Dir. vs Time for Beam Element During Free-vibration Analysis") # Beam Element IN Y DIR.

# PLOT STRUCTURAL PERIOD DURING THE ANALYSIS
plt.figure(0, figsize=(12, 8))
plt.plot(time_FV, PERIOD_MIN_FV, linewidth=3)
plt.plot(time_FV, PERIOD_MAX_FV, linewidth=3)
plt.title('Period of Structure During Free-vibration Analysis')
plt.ylabel('Structural Period [s]')
plt.xlabel('Time [s]')
#plt.semilogy()
plt.grid()
plt.legend([f'PERIOD - MIN VALUES: Min: {np.min(PERIOD_MIN_FV):.3f} (s) - Mean: {np.mean(PERIOD_MIN_FV):.3f} (s) - Max: {np.max(PERIOD_MIN_FV):.3f} (s)', 
            f'PERIOD - MAX VALUES:  Min: {np.min(PERIOD_MAX_FV):.3f} (s) - Mean: {np.mean(PERIOD_MAX_FV):.3f} (s) - Max: {np.max(PERIOD_MAX_FV):.3f} (s)',
            ])
plt.show()

# Plot Stress-strain of top and bottom fibers
plt.figure(44, figsize=(12, 8))
plt.plot(STRAINb_FV, STRESSb_FV, linewidth=4, label=f"BOT FIBER - Stress: {np.abs(STRESSb_FV[-1]) : .3f}")
plt.plot(STRAINt_FV, STRESSt_FV, linewidth=4, label=f"TOP FIBER - Stress: {np.abs(STRESSt_FV[-1]) : .3f}")
plt.xlabel('Strain (mm/mm)')
plt.ylabel('Stress (N/mm^2)')
plt.title('Behavior of beam - Stress-strain of Top and Bottom Fibers Curve')
plt.grid(True)
plt.legend()
plt.show()

S01.PLOT_2D_FRAME(deformed_scale=1.0)
#%%----------------------------------------------------
# SEISMIC ANALYSIS (DYNAMIC TIME-HISTORY ANALYSIS)
#ELE_TYPE = 'elasticBeamColumn'      # MATERIAL LINEARITY
ELE_TYPE = 'nonlinearBeamColumn'     # MATERIAL AND GEOMETRIC NONLINEARITIES
MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'
ANAL_TYPE = 'SEISMIC'

DATA = BRIDGE_VIERENDEEL_GIRDERS(LENGTH, HEIGHT, MAT_TYPE, ELE_TYPE, ANAL_TYPE, TOTAL_MASS)

(time_SEI, reaction_SEI, disp_mid_SEI,
 ele_axialforce_SEI, ele_shearforce_SEI, ele_momentforce_SEI,
 node_displacementsX_SEI, node_displacementsY_SEI,
 dispX_SEI, dispY_SEI,
 veloX_SEI, veloY_SEI,
 accX_SEI, accY_SEI,
 stiffness_SEI, PERIOD_SEI, damping_ratio_SEI,
 PERIOD_MIN_SEI, PERIOD_MAX_SEI,
 STRESSb_SEI, STRAINb_SEI, STRESSt_SEI, STRAINt_SEI, CURVATURE_SEI) = DATA


XDATA = disp_mid_SEI
YDATA = reaction_SEI
XLABEL = 'Displacement in Middle Span [mm]'
YLABEL = 'Base Reaction [N]'
TITLE = 'Base Reaction and Dispalcement of Structure During Seismic Analysis'
COLOR = 'black'
SEMILOGY = False
PLOT(XDATA, YDATA, TITLE, XLABEL, YLABEL, COLOR, SEMILOGY)

PLOT_TIME_HISTORY(time_SEI, reaction_SEI, disp_mid_SEI,
                          dispX_SEI, dispY_SEI,
                          veloX_SEI, veloY_SEI,
                          accX_SEI, accY_SEI)

# PLOT ELEMENTS AXIAL FORCE
YLABEL = 'Element Axial Force [N]'
TITLE = "elements  force in X Dir. vs Time for Beam Element During Seismic Analysis"
PLOT_FORCES(time_SEI, ele_axialforce_SEI, YLABEL, TITLE) # Beam Element - ELEMENTS AXIAL FORCE
# PLOT ELEMENTS SHEAR FORCE
YLABEL = 'Element Shear Force [N]'  
TITLE = "elements  force in Y Dir. vs Time for Beam Element During Seismic Analysis"
PLOT_FORCES(time_SEI, ele_shearforce_SEI, YLABEL, TITLE) # Beam Element - ELEMENTS SHEAR FORCE
# PLOT ELEMENTS MOMENT FORCE
YLABEL = 'Element Moment Force [N.mm]'  
TITLE = "elements  force in Z Dir. vs Time for Beam Element During Seismic Analysis"
PLOT_FORCES(time_SEI, ele_momentforce_SEI, YLABEL, TITLE) # Beam Element - ELEMENTS MOMENT FORCE

# PLOT NODAL DISPALEMENTS
PLOT_DISPLAEMENTS(time_SEI, node_displacementsX_SEI, TITLE = "Node Displacements in X Dir. vs Time for Beam Element During Seismic Analysis") # Beam Element IN X DIR.
PLOT_DISPLAEMENTS(time_SEI, node_displacementsY_SEI, TITLE = "Node Displacements in Y Dir. vs Time for Beam Element During Seismic Analysis") # Beam Element IN Y DIR.

# PLOT STRUCTURAL PERIOD DURING THE ANALYSIS
plt.figure(0, figsize=(12, 8))
plt.plot(time_SEI, PERIOD_MIN_SEI, linewidth=3)
plt.plot(time_SEI, PERIOD_MAX_SEI, linewidth=3)
plt.title('Period of Structure During Seismic Analysis')
plt.ylabel('Structural Period [s]')
plt.xlabel('Time [s]')
#plt.semilogy()
plt.grid()
plt.legend([f'PERIOD - MIN VALUES: Min: {np.min(PERIOD_MIN_SEI):.3f} (s) - Mean: {np.mean(PERIOD_MIN_SEI):.3f} (s) - Max: {np.max(PERIOD_MIN_SEI):.3f} (s)', 
            f'PERIOD - MAX VALUES:  Min: {np.min(PERIOD_MAX_SEI):.3f} (s) - Mean: {np.mean(PERIOD_MAX_SEI):.3f} (s) - Max: {np.max(PERIOD_MAX_SEI):.3f} (s)',
            ])
plt.show()

# Plot Stress-strain of top and bottom fibers
plt.figure(44, figsize=(12, 8))
plt.plot(STRAINb_SEI, STRESSb_SEI, linewidth=4, label=f"BOT FIBER - Stress: {np.abs(STRESSb_SEI[-1]) : .3f}")
plt.plot(STRAINt_SEI, STRESSt_SEI, linewidth=4, label=f"TOP FIBER - Stress: {np.abs(STRESSt_SEI[-1]) : .3f}")
plt.xlabel('Strain (mm/mm)')
plt.ylabel('Stress (N/mm^2)')
plt.title('Behavior of beam - Stress-strain of Top and Bottom Fibers Curve')
plt.grid(True)
plt.legend()
plt.show()

S01.PLOT_2D_FRAME(deformed_scale=1.0)
#%%----------------------------------------------------
# --------------------------------------
#  Plot BaseShear-Displacement Analysis 
# --------------------------------------
XX = np.abs(disp_mid_PUSH); YY = np.abs(reaction_PUSH); # ABSOLUTE VALUE
SLOPE_NODE = 10

DATA = S07.BILNEAR_CURVE(XX, YY, SLOPE_NODE)
X, Y, Elastic_ST, Plastic_ST, Tangent_ST, Ductility_Rito, Over_Strength_Factor = DATA

XLABEL = 'Displacement in Y [mm]'
YLABEL = 'Base-Shear Reaction [N]'
LEGEND01 = 'Curve'
LEGEND02 = 'Bilinear Fitted'
LEGEND03 = 'Undefined'
TITLE = f'Last Data of Base Axial-Displacement Analysis - Ductility Ratio: {X[2]/X[1]:.4f} - Over Strength Factor: {Y[2]/Y[1]:.4f}'
COLOR = 'black'
PLOT_2D(np.abs(disp_mid_PUSH), np.abs(reaction_PUSH), X, Y, X, Y, XLABEL, YLABEL, TITLE, LEGEND01, LEGEND02, LEGEND03, COLOR='black', Z=2) 
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
Dd = np.max(np.abs(disp_mid_SEI))
DIy = (Dd - X[1]) /(X[2] - X[1])
print(f'Structural Ductility Damage Index in Y Direction: {100*DIy:.4f} (%)')
#%%----------------------------------------------------
# EVALUATION OF DISSIPATED ENERGY CAPACITY INDEX
def DISSIPATED_ENERGY_FUN_WITH_PLOT(displacement, base_shear, title="Hysteresis Curve"):
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
    ax.set_title(title)
    ax.set_xlabel("Displacement (mm)")
    ax.set_ylabel("Base Shear (N)")
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.legend()

    return area, fig

Ed_SEI, fig_SEI = DISSIPATED_ENERGY_FUN_WITH_PLOT(
    disp_mid_SEI, reaction_SEI, 
    title="Earthquake Response – Dissipated Energy (Convex Hull)"
)
fig_SEI.show()

print(f"Dissipated Energy from Earthquake= {Ed_SEI:.2f} N·mm")

Ed_CP, fig_CP = DISSIPATED_ENERGY_FUN_WITH_PLOT(
    disp_mid_CP, reaction_CP,
    title="Cyclic Loading – Dissipated Energy (Convex Hull)"
)
fig_CP.show()

print(f"Dissipated Energy from Cyclic Displacement= {Ed_CP:.2f} N·mm")


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
#%%----------------------------------------------------
# %% FRAGILITY ANALYSIS
#exec(open("FRAGILITY_FUN.py").read())
import FRAGILITY_FUN as FF

damage_states = {
    'Minor Damage Level': (20, 20),# Median DI=20%, β=40%
    'Moderate Damage Level': (40, 20),
    'Severe Damage Level': (60, 20),
    'Failure Level': (100, 20)
}

# Intensity Measure (IM) values from 0.0 to 100.0
Ddi = np.abs(disp_mid_SEI)
DIyi = 100*(Ddi - X[1]) /(X[2] - X[1])
for KK in range(len(DIyi)):
    if DIyi[KK] <= 0.0: # DAMAGE INDEX -> 0.0 <= DI <=100 
        DIyi[KK] = 0.0
    if DIyi[KK] >= 100: 
        DIyi[KK] = 100.0    
        
im_values = DIyi # Structural Ductility Damage Index
TITLE = 'Structural Ductility Damage Index (%)  [IM]'
FF.FRAGILITY_ANALYSIS(damage_states, im_values, TITLE, SCATTER='True', SEMI_LOG='False')
#%%----------------------------------------------------
# EXCEL OUTPUT
import pandas as pd

# Create DataFrame function
def create_df(reaction, disp_mid, dispX, dispY, PERIOD_MIN, PERIOD_MAX):
    df = pd.DataFrame({
        "reaction": reaction,
        "disp_mid": disp_mid,
        "dispX": dispX,
        "dispY": dispY,
        "PERIOD_MIN": PERIOD_MIN,
        "PERIOD_MAX": PERIOD_MAX,        
    })
    return df


# Save to Excel
with pd.ExcelWriter("P(t)_PEDESTRIAN_BRIDGE_VIERENDEEL_GIRDERS_8_ANA_OUTPUT.xlsx", engine='openpyxl') as writer:
    
    # PUSHOVER
    df1 = create_df(reaction_PUSH, disp_mid_PUSH, dispX_PUSH, dispY_PUSH, PERIOD_MIN_PUSH, PERIOD_MAX_PUSH)
    df1.to_excel(writer, sheet_name="PUSHOVER", index=False)
                 
    # CYCLIC DISPLACEMENT
    df1 = create_df(reaction_CP, disp_mid_CP, dispX_CP, dispY_CP, PERIOD_MIN_CP, PERIOD_MAX_CP)
    df1.to_excel(writer, sheet_name="CYCLIC_DISPLACEMENT", index=False)
    
    # STATIC EXTERNAL TIME-DEPENDENT LOADING
    df2 = create_df(reaction_ETDLS, disp_mid_ETDLS, dispX_ETDLS, dispY_ETDLS, PERIOD_MIN_ETDLS, PERIOD_MAX_ETDLS)
    df2.to_excel(writer, sheet_name="STATIC_EXTERNAL_TIME-DEPENDENT_LOADING", index=False)

    # DYNAMIC EXTERNAL TIME-DEPENDENT LOADING
    df3 = create_df(reaction_ETDLD, disp_mid_ETDLD, dispX_ETDLD, dispY_ETDLD, PERIOD_MIN_ETDLD, PERIOD_MAX_ETDLD)
    df3.to_excel(writer, sheet_name="DYNAMIC_EXTERNAL_TIME-DEPENDENT_LOADING", index=False)
    
    # FREE-VIBRATION
    df3 = create_df(reaction_FV, disp_mid_FV, dispX_FV, dispY_FV, PERIOD_MIN_FV, PERIOD_MAX_FV)
    df3.to_excel(writer, sheet_name="FREE-VIBRATION", index=False)

    # SEISMIC
    df4 = create_df(reaction_SEI, disp_mid_SEI, dispX_SEI, dispY_SEI, PERIOD_MIN_SEI, PERIOD_MAX_SEI)
    df4.to_excel(writer, sheet_name="SEISMIC", index=False)
#%%----------------------------------------------------