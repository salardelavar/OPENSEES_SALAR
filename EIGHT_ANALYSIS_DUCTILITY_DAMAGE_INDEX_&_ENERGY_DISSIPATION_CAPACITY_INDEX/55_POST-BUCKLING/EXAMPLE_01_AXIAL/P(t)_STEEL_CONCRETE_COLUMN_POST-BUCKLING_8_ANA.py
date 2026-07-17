###########################################################################################################
#                   >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                      #
#COMPREHENSIVE NONLINEAR SEISMIC ASSESSMENT OF MULTI-MODE POST-BUCKLING PHENOMENA A COLUMN WITH SEMI-RIGID#
#  CONNECTIONS: AN OPENSEES FRAMEWORK FOR STATIC PUSHOVER, CYCLIC DEGRADATION, STATIC TIME-HISTORY AND    #
#                                 DYNAMIC TIME-HISTORY ANALYSIS                                           #
#---------------------------------------------------------------------------------------------------------#
#                                                   P(t) = P0 sin(wt)                                     #
#                                           P(t) = P0 exp(-0.05wt) sin(wt)                                #
#---------------------------------------------------------------------------------------------------------#
#         ASSESSMENT OF DUCTILITY DAMAGE INDICES FOR STRUCTURAL ELEMENTS AND SYSTEMS AND EVALUATION       #
#                                        OF ENERGY DISSIPATION CAPACITY INDEX                             #
#---------------------------------------------------------------------------------------------------------#
#                    THIS PYTHON SCRIPT WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                     #
#                                   EMAIL: salar.d.ghashghaei@gmail.com                                   #
###########################################################################################################
"""
Nonlinear Seismic Performance Assessment of a column:
An OpenSeesPy Framework for Material and Geometric Nonlinearity Under Static, Cyclic, and Earthquake Loading

This OpenSeesPy script performs rigorous nonlinear static and dynamic analysis of a column
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
"""
# BOOK: Dynamics of Structure and Foundation; A Unified Approach-1. Fundamentals-Indrajit Chowdhury-CRC Press-2009
# BOOK: Dynamics of structure and foundation; A Uniﬁed Approach-2. Applications-Indrajit Chowdhury-CRC Press-2009
# BOOK: Plastic Design and Second-Order Analysis of Steel Frames-W. F. Chen, I. Sohal-Springer New York-1995
# BOOK: Structural Stability-W. F. Chen'
# BOOK: Differential Equations for Engineers-Wei-Chau Xie-CAMBRIDGE-2010
# BOOK: Structural dynamics THEORY AND COMPUTATION MARIO BAZ-5 EDITION   
# BOOK: Dynamics of Structures in SI Units -  Anil Kumar Chopra
'https://share.google/TDN5O4eWmw5pH8zUH' 
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
import COMPOSITE_REC_CROSS_SHAPED_STEEL_SECTION_FUN_EXTRA as S11
import COMPOSITE_SECTION_DIFF_WIDTH_HEIGHT_FUN_EXTRA_TWO as S12
import CONCRETE_REC_SECTION_FUN_EXTRA as S13
import CROSS_SHAPED_STEEL_SECTION_DIFF_WIDTH_HEIGHT_FUN_EXTRA as S14
import STEEL_CORNER_BOX_SECTION_DIFF_WIDTH_HEIGHT_FUN_EXTRA as S15
import PLATE_DOUBLE_I_STEEL_SECTION_DIFF_WIDTH_HEIGHT_FUN_EXTRA as S16
import SECTION_CRACK_DEPTH_FUN as S18
import OPENSEEES_HYSTERETICSM_FORCE_DISP_FUN as S033
import OPENSEEES_HYSTERETIC_FORCE_DISP_FUN as S044
import EQULIVALENT_VISCOUS_DAMPING_RATIO_FUN as S055
#%%----------------------------------------------------
def POST_BUCKLING_COLUMN(LENGTH, MAT_TYPE, ELE_TYPE, ANAL_TYPE, TOTAL_MASS, MODE_NUM, SEMI_RIGID):
    # Initialize model
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    
    MAX_ITERATIONS = 5000   # Maximum number of iterations
    MAX_TOLERANCE = 1.0e-9  # Specified tolerance for convergence
    
    N = 100             # Number of elements
    dy = LENGTH / N     # Element length
    ε = 0.001 * LENGTH  # [mm] Imperfection amplitude
    
    #%% GEOMETRY DEFINITION
    # Define nodes with initial imperfection ((1/MODE_NUM) sine wave)
    for i in range(N+1):
        y = i * dy
        if MODE_NUM == 0:
          x = 0.0
        else:  
            x = ε * np.sin(MODE_NUM * np.pi * y / LENGTH)  # Imperfection in x-direction
        ops.node(i+1, x, y)    

    #------------------------------------------------------------------------------------------------
    # Define Semi-rigid connection with rotational spring
    if SEMI_RIGID == True:
        ops.node(N+2, 0, 0)      # Rotational spring in the bottom of column   
        ops.fix(N+2, 1, 1, 1)
        ops.node(N+3, 0, LENGTH) # Rotational spring in the top of column   
        ops.fix(N+3, 1, 1, 1)
    #------------------------------------------------------------------------------------------------
    # Define boundary conditions for simply supported column
    ops.fix(1, 1, 1, 0)    # Bottom node: fix u_x, u_y, free theta
    ops.fix(N+1, 1, 0, 0)  # Top node: fix u_x, free u_y, free theta
    #------------------------------------------------------------------------------------------------
    # Define Nonlinear Connection in Bottom and Top Rotational Springs
    if SEMI_RIGID == True:
        matTag_S = 500
        MY = 12430.0           # [N.mm] Yield Moment of springs
        TY = 0.001474          # [Rad] Yield Rotation of springs
        MU = 1.12 * MY         # [N.mm] Utimate Moment of springs
        TSU = 0.0541           # [Rad] Utimate Rotation of springs
        #ROTATIONAl_STIFFNESS = 1e14 # [N/mm] Rotational stiffness of springs
        pinchX = 0.4           # Pinching factor in X direction
        pinchY = 0.2           # Pinching factor in Y direction
        damage1 = 0.0          # Damage due to ductility
        damage2 = 0.0          # Damage due to energy
        beta = 0.1             # Stiffness degradation parameter
        # MOMENT-ROTATION RELATIONSHIP OF ROTATIONAL SPRING AND PLOT 
        DP = [0, 0, 0, 0]
        FP = [0, 0, 0, 0]
        DN = [0, 0, 0, 0]
        FN = [0, 0, 0, 0]
        DP[0], FP[0] = TY, MY
        DP[1], FP[1] = TSU, MU 
        DP[2], FP[2] = 1.13*TSU, 0.43*MU
        DP[3], FP[3] = 1.35*TSU, 0.14*MU
        DN[0], FN[0] = -TY, -MY 
        DN[1], FN[1] = -TSU, -MU
        DN[2], FN[2] = -1.07*TSU, -0.43*MU
        DN[3], FN[3] = -1.24*TSU, -0.12*MU   
        S033.OPENSEEES_HYSTERETICSM_FORCE_DISP_FUN(matTag_S, DP, FP, DN, FN, 
                                            pinchX, pinchY,
                                            damage1, damage2, beta,
                                            PLOT = True, X_LABEL='Rotation (rad)', Y_LABEL='Moment [N.mm]', TITLE='MOMENT-ROTATION CURVE')
    
        #ops.uniaxialMaterial('Elastic', matTag_S, ROTATIONAl_STIFFNESS)
        ops.element('zeroLength', N+200, N+2, 1, '-mat', matTag_S, '-dir', 3)   # ROTATIONAL SPRING AT BOTTOM OF COLUMN
        ops.element('zeroLength', N+300, N+1, N+3, '-mat', matTag_S, '-dir', 3) # ROTATIONAL SPRING AT TOP OF COLUMN  
    #------------------------------------------------------------------------------------------------
    
    #%% CONCRETE SECTION
    fc = 30.0        # [N/mm²] Concrete Compressive Strength
    Kfc = 1.24       # ratio of confined to unconfined concrete strength
    SEC_TAG = 1
    nFib = 40
    STEEL_DENSITY = 7850/1e9         # [kg/m^3] -> [kg/mm^3] Steel Material Density
    CONCRETE_DENSITY = 2500/1e9      # [kg/m^3] -> [kg/mm^3] Concrete Material Density
    STEEL_TYPE = 'INELASTIC'
    """
    Bsec, Hsec, cover = 400.0, 400.0, 50.0 # [mm]
    Depth, Ele_Mass = S13.CONCRETE_REC_SECTION_FUN(SEC_TAG, STEEL_TYPE, fc, Kfc,
                               Bsec, Hsec, cover,
                               nFib, CONCRETE_DENSITY,
                               plot=True)
    """
    """
    #%% COMPOSITE SECTION
    HSec, BSec, cover = 400.0, 400.0, 50.0
    Depth, Ele_Mass = S11.COMPOSITE_REC_CROSS_SHAPED_STEEL_SECTION_FUN_EXTRA(SEC_TAG, HSec, BSec, cover,
                              STEEL_TYPE, fc, Kfc,
                              STEEL_DENSITY,
                              CONCRETE_DENSITY, PLOT=True)
   
    """
    """
    Depth, Ele_Mass = S12.COMPOSITE_SECTION_DIFF_WIDTH_HEIGHT_FUN_EXTRA_TWO(SEC_TAG, STEEL_TYPE, fc, Kfc,
                                                                            CONCRETE_DENSITY, plot=True)
    """
    #%% STEEL SECTION
    #Depth, Ele_Mass = S14.STEEL_SECTION_DIFF_WIDTH_HEIGHT_FUN_EXTRA(SEC_TAG, STEEL_TYPE, STEEL_DENSITY, plot=True)
    #Depth, Ele_Mass = S15.STEEL_CORNER_BOX_SECTION_DIFF_WIDTH_HEIGHT_FUN_EXTRA(SEC_TAG, STEEL_TYPE, STEEL_DENSITY, plot=True)
    Depth, Ele_Mass = S16.STEEL_SECTION_DIFF_WIDTH_HEIGHT_FUN_EXTRA(SEC_TAG, STEEL_TYPE, STEEL_DENSITY, plot=True)
    #%% CREATE ELEMENTS
    ele_tag = 1
    transfTag = 1
    #ops.geomTransf('Linear', transfTag)
    #ops.geomTransf('PDelta', transfTag)
    ops.geomTransf('Corotational', transfTag)
    numIntgrPts = 5
    if ELE_TYPE == 'elasticBeamColumn':
        for i in range(1, N+1):
            ops.element('elasticBeamColumn', i, i, i+1, SEC_TAG, transfTag,'-mass', Ele_Mass)
    if ELE_TYPE == 'nonlinearBeamColumn': 
        for i in range(1, N+1):
            ops.element('nonlinearBeamColumn', i, i, i+1, numIntgrPts, SEC_TAG, transfTag,'-mass', Ele_Mass) 
    ele_tag = N
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
    center_node = N+1
    TS_TAG = 1
    PATT_TAG = 1
    ops.timeSeries('Linear', TS_TAG)
    ops.pattern('Plain', PATT_TAG, TS_TAG)
    for ii in range(1, tag):
        ops.load(ii, 0.0, -ENW, 0.0)   # [N] downward
        ops.mass(ii, ENM, ENM, 0.0)
        
    #ops.mass(7, ENM, ENM, 0.0)
    TS_TAG = 2
    PATT_TAG = 2
    ops.timeSeries('Linear', TS_TAG)
    ops.pattern('Plain', PATT_TAG, TS_TAG)       
    #for ii in range(1, tag-1):
    #    ops.eleLoad('-ele', ii,'-type', '-beamUniform', -0.0003, 0.0) # [N/mm] Uniformly-distributed load  
    ops.load(center_node, 0.0, -1.0, 0.0)
    
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
    reactionX, reactionY, reactionZ, disp_mid = [], [], [], []
    stiffness = []
    OMEGA, PERIOD = [], []
    PERIOD_MIN, PERIOD_MAX = [], []
    STRESSb = []; STRAINb = []; 
    STRESSt = []; STRAINt = [];
    CURVATURE = [];
    NA, CD = [], [];
    # Initialize lists for each element's axial force, shear force and moment force
    ele_axialforce = {tag: [] for tag in range(1, N - 1)}
    ele_shearforce = {tag: [] for tag in range(1, N - 1)}
    ele_momentforce = {tag: [] for tag in range(1, N - 1)}
    # Initialize lists for each node's displacement
    node_displacementsX = {tag: [] for tag in range(1, N + 1)}
    node_displacementsY = {tag: [] for tag in range(1, N + 1)}
    
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
        reactionX.append(ops.nodeReaction(1, 1))                        # SHEAR BASE REACTION
        reactionY.append(ops.nodeReaction(1, 2))                        # AXIAL BASE REACTION
        reactionZ.append(ops.nodeReaction(1, 3))                        # MOMENT BASE REACTION
        disp_mid.append(ops.nodeDisp(center_node, 2))                   # DISPLACEMENT NODE 07 
        dispX.append(ops.nodeDisp(center_node, 1))                      # DISPLACEMENT NODE 07 IN X DIR 
        dispY.append(ops.nodeDisp(center_node, 2))                      # DISPLACEMENT NODE 07 IN Y DIR 
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
            
        DATA = (reactionX, reactionY, reactionZ, disp_mid,
                ele_axialforce, ele_shearforce, ele_momentforce,
                node_displacementsX, node_displacementsY,
                dispX, dispY)
        
        return  DATA
    
    if ANAL_TYPE == 'PUSHOVER': # STATIC TIME-HISTORY ANALYSIS
        IDctrlDOF = 2   # 1: Horizental Dispalcement - 2: Vertical Dispalcement
        DINCR = -0.01   # [mm] Incremental Displacement
        DMAX = -25.0    # [mm] Max. Displacement
        ops.integrator('DisplacementControl', center_node, IDctrlDOF, DINCR)
        ops.analysis('Static')
        Nsteps =  int(np.abs(DMAX / DINCR)) 
        STEP = 0.0
        for step in range(Nsteps):
            OK = ops.analyze(1)
            S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS)
            ops.reactions()
            reactionX.append(ops.nodeReaction(1, 1))                        # SHEAR BASE REACTION
            reactionY.append(ops.nodeReaction(1, 2))                        # AXIAL BASE REACTION
            reactionZ.append(ops.nodeReaction(1, 3))                        # MOMENT BASE REACTION
            disp_mid.append(ops.nodeDisp(center_node, 2))                   # DISPLACEMENT NODE 07 
            dispX.append(ops.nodeDisp(center_node, 1))                      # DISPLACEMENT NODE 07 IN X DIR 
            dispY.append(ops.nodeDisp(center_node, 2))                      # DISPLACEMENT NODE 07 IN Y DIR 
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
            if ELE_TYPE == 'nonlinearBeamColumn': 
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
                # SECTION NEUTRAL AXIS & CRACK DEPTH
                na, cd = S18.SECTION_CRACK_DEPTH_FUN(center_node-1, SEC_TAG,
                                        -0.15,          # TOP FIBER MATERIAL ULTIMATE STRAIN FOR CRACK DEPTH EVALUTION
                                        -Depth*0.5, 0.0, # BOTTOM FIBER FROM NEUTRAL AXIS
                                        +Depth*0.5, 0.0, # TOP FIBER FROM NEUTRAL AXIS 
                                        )
                NA.append(na); CD.append(cd);
            else:
                STRESSb.append(0.0)
                STRAINb.append(0.0)
                STRESSt.append(0.0)
                STRAINt.append(0.0)
                CURVATURE.append(0.0)
                NA.append(0.0); CD.append(0.0);
            STEP += 1
            #print(STEP, dispY[-1], reaction[-1])
            print(f"Step: {STEP}, Displacement: {dispY[-1]:.4f} mm, Reaction: {reactionY[-1]:.2f} N")     
        else:
            print('\n\nPUSHOVER ANALYSIS DONE.\n\n')
            
        DATA = (reactionX, reactionY, reactionZ, disp_mid,
                ele_axialforce, ele_shearforce, ele_momentforce,
                node_displacementsX, node_displacementsY,
                dispX, dispY,
                np.array(PERIOD_MIN), np.array(PERIOD_MAX),
                STRESSb, STRAINb, STRESSt, STRAINt, CURVATURE, NA, CD)
    
        return  DATA
    
    if ANAL_TYPE == 'CYCLIC_DISPLACEMENT': # STATIC TIME-HISTORY ANALYSIS
        IDctrlDOF = 2   # 1: Horizental Dispalcement - 2: Vertical Dispalcement
        DMAX = -25.0    # [mm] Max. Displacement
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
            current_disp = ops.nodeDisp(center_node, 2) # DISPALCEMENT APPLIED IN NODE IN X DIR.
            dU = target_disp - current_disp
            ops.integrator('DisplacementControl', center_node, IDctrlDOF, dU)
            OK = ops.analyze(1)
            S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS)
            ops.reactions()
            reactionX.append(ops.nodeReaction(1, 1))                        # SHEAR BASE REACTION
            reactionY.append(ops.nodeReaction(1, 2))                        # AXIAL BASE REACTION
            reactionZ.append(ops.nodeReaction(1, 3))                        # MOMENT BASE REACTIO
            disp_mid.append(ops.nodeDisp(center_node, 2))                   # DISPLACEMENT NODE 07
            dispX.append(ops.nodeDisp(center_node, 1))                      # DISPLACEMENT NODE 07 IN X DIR 
            dispY.append(ops.nodeDisp(center_node, 2))                      # DISPLACEMENT NODE 07 IN Y DIR 
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
            if ELE_TYPE == 'nonlinearBeamColumn': 
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
                # SECTION NEUTRAL AXIS & CRACK DEPTH
                na, cd = S18.SECTION_CRACK_DEPTH_FUN(center_node-1, SEC_TAG,
                                        -0.15,          # TOP FIBER MATERIAL ULTIMATE STRAIN FOR CRACK DEPTH EVALUTION
                                        -Depth*0.5, 0.0, # BOTTOM FIBER FROM NEUTRAL AXIS
                                        +Depth*0.5, 0.0, # TOP FIBER FROM NEUTRAL AXIS 
                                        )
                NA.append(na); CD.append(cd);
            else:
                STRESSb.append(0.0)
                STRAINb.append(0.0)
                STRESSt.append(0.0)
                STRAINt.append(0.0)
                CURVATURE.append(0.0)
                NA.append(0.0); CD.append(0.0);
            STEP += 1
            #print(STEP, dispY[-1], reaction[-1])
            print(f"Step: {STEP}, Displacement: {dispY[-1]:.4f} mm, Reaction: {reactionY[-1]:.2f} N")     
        else:
            print('\n\nCYCLIC DISPLAEMENT ANALYSIS DONE.\n\n')
            
        DATA = (reactionX, reactionY, reactionZ, disp_mid,
                ele_axialforce, ele_shearforce, ele_momentforce,
                node_displacementsX, node_displacementsY,
                dispX, dispY,
                np.array(PERIOD_MIN), np.array(PERIOD_MAX),
                STRESSb, STRAINb, STRESSt, STRAINt, CURVATURE, NA, CD)
    
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
        #for ii in range(1, tag-1):# W(x,t) = W0 exp(-0.05wt) sin(wt) 
        #    ops.eleLoad('-ele', ii,'-type', '-beamUniform', -0.0003, 0.0) # [N/mm] Uniformly-distributed load
        
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
            reactionX.append(ops.nodeReaction(1, 1))                        # SHEAR BASE REACTION
            reactionY.append(ops.nodeReaction(1, 2))                        # AXIAL BASE REACTION
            reactionZ.append(ops.nodeReaction(1, 3))                        # MOMENT BASE REACTION
            disp_mid.append(ops.nodeDisp(center_node, 2))                   # DISPLACEMENT NODE 07 
            dispX.append(ops.nodeDisp(center_node, 1))                      # DISPLACEMENT NODE 07 IN X DIR 
            dispY.append(ops.nodeDisp(center_node, 2))                      # DISPLACEMENT NODE 07 IN Y DIR 
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
            if ELE_TYPE == 'nonlinearBeamColumn': 
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
                # SECTION NEUTRAL AXIS & CRACK DEPTH
                na, cd = S18.SECTION_CRACK_DEPTH_FUN(center_node-1, SEC_TAG,
                                        -0.15,          # TOP FIBER MATERIAL ULTIMATE STRAIN FOR CRACK DEPTH EVALUTION
                                        -Depth*0.5, 0.0, # BOTTOM FIBER FROM NEUTRAL AXIS
                                        +Depth*0.5, 0.0, # TOP FIBER FROM NEUTRAL AXIS 
                                        )
                NA.append(na); CD.append(cd);
            else:
                STRESSb.append(0.0)
                STRAINb.append(0.0)
                STRESSt.append(0.0)
                STRAINt.append(0.0)
                CURVATURE.append(0.0)
                NA.append(0.0); CD.append(0.0);
            STEP += 1
            #print(STEP, dispY[-1], reaction[-1])
            print(f"Step: {STEP}, Displacement: {dispY[-1]:.4f} mm, Reaction: {reactionY[-1]:.2f} N")     
        else:
            print('\n\nSTATIC EXTERNAL TIME-DEPENDENT LOADING ANALYSIS DONE.\n\n')
            
        DATA = (reactionX, reactionY, reactionZ, disp_mid,
                ele_axialforce, ele_shearforce, ele_momentforce,
                node_displacementsX, node_displacementsY,
                dispX, dispY,
                np.array(PERIOD_MIN), np.array(PERIOD_MAX),
                STRESSb, STRAINb, STRESSt, STRAINt, CURVATURE, NA, CD)
    
        return  DATA
    
    if ANAL_TYPE == 'DYNAMIC_EXTERNAL_TIME-DEPENDENT_LOADING': # DYNAMIC TIME-HISTORY ANALYSIS
        #%% DEFINE EXTERNAL TIME-DEPENDENT LOADING PROPERTIES
        duration = 20.0             # [s] Analysis duration
        dt = 0.01                   # [s] Time step
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
        #for ii in range(1, tag-1):# W(x,t) = W0 exp(-0.05wt) sin(wt) 
        #    ops.eleLoad('-ele', ii,'-type', '-beamUniform', -0.0003, 0.0) # [N/mm] Uniformly-distributed load
        
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
            reactionX.append(ops.nodeReaction(1, 1))              # SHEAR BASE REACTION
            reactionY.append(ops.nodeReaction(1, 2))              # AXIAL BASE REACTION
            reactionZ.append(ops.nodeReaction(1, 3))              # MOMENT BASE REACTION
            disp_mid.append(ops.nodeDisp(center_node, 2))         # DISPLACEMENT NODE 07 
            dispX.append(ops.nodeDisp(center_node, 1))            # DISPLACEMENT NODE 07 IN X DIR 
            dispY.append(ops.nodeDisp(center_node, 2))            # DISPLACEMENT NODE 07 IN Y DIR 
            veloX.append(ops.nodeVel(center_node, 1))             # VELOCITY NODE 07
            veloY.append(ops.nodeVel(center_node, 2))             # VELOCITY NODE 07
            accX.append(ops.nodeAccel(center_node, 1))            # ACCELERATION NODE 07
            accY.append(ops.nodeAccel(center_node, 2))            # ACCELERATION NODE 07
            stiffness.append(np.abs(reactionX[-1]) / np.abs(disp_mid[-1]))
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
            if ELE_TYPE == 'nonlinearBeamColumn': 
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
                # SECTION NEUTRAL AXIS & CRACK DEPTH
                na, cd = S18.SECTION_CRACK_DEPTH_FUN(center_node-1, SEC_TAG,
                                        -0.15,          # TOP FIBER MATERIAL ULTIMATE STRAIN FOR CRACK DEPTH EVALUTION
                                        -Depth*0.5, 0.0, # BOTTOM FIBER FROM NEUTRAL AXIS
                                        +Depth*0.5, 0.0, # TOP FIBER FROM NEUTRAL AXIS 
                                        )
                NA.append(na); CD.append(cd);
            else:
                STRESSb.append(0.0)
                STRAINb.append(0.0)
                STRESSt.append(0.0)
                STRAINt.append(0.0)
                CURVATURE.append(0.0)
                NA.append(0.0); CD.append(0.0);
            #print(time[-1], dispY[-1], veloY[-1])
            print(f"Time: {time[-1]:.4f}, Displacement: {dispY[-1]:.4f} mm, Reaction: {reactionY[-1]:.2f} N")      
        else:
            print('\n\nDYNAMIC EXTERNAL TIME-DEPENDENT LOADING ANALYSIS DONE.\n\n')  
        # Calculating Damping Ratio and Period Using Logarithmic Decrement Analysis 
        damping_ratio = S04.DAMPING_RATIO(dispY)  
        
        # Compute modal properties
        ops.modalProperties("-print", "-file", f"SALAR_ModalReport_{ANAL_TYPE}.txt", "-unorm") 
        
        DATA = (time, reactionX, reactionY, reactionZ, disp_mid,
                ele_axialforce, ele_shearforce, ele_momentforce,
                node_displacementsX, node_displacementsY,
                dispX, dispY,
                veloX, veloY,
                accX, accY,
                stiffness, PERIOD, damping_ratio,
                np.array(PERIOD_MIN), np.array(PERIOD_MAX),
                STRESSb, STRAINb, STRESSt, STRAINt, CURVATURE, NA, CD)
        
        return  DATA
        
    if ANAL_TYPE == 'FREE-VIBRATION': # DYNAMIC TIME-HISTORY ANALYSIS
        #%% DEFINE PARAMETERS FOR FREE-VIBRATION ANALYSIS
        u0 = -0.0035                       # [mm] Initial displacement
        v0 = 0.015                         # [mm/s] Initial velocity
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
            reactionX.append(ops.nodeReaction(1, 1))              # SHEAR BASE REACTION
            reactionY.append(ops.nodeReaction(1, 2))              # AXIAL BASE REACTION
            reactionZ.append(ops.nodeReaction(1, 3))              # MOMENT BASE REACTION
            disp_mid.append(ops.nodeDisp(center_node, 2))         # DISPLACEMENT NODE 07 
            dispX.append(ops.nodeDisp(center_node, 1))            # DISPLACEMENT NODE 07 IN X DIR 
            dispY.append(ops.nodeDisp(center_node, 2))            # DISPLACEMENT NODE 07 IN Y DIR 
            veloX.append(ops.nodeVel(center_node, 1))             # VELOCITY NODE 07
            veloY.append(ops.nodeVel(center_node, 2))             # VELOCITY NODE 07
            accX.append(ops.nodeAccel(center_node, 1))            # ACCELERATION NODE 07
            accY.append(ops.nodeAccel(center_node, 2))            # ACCELERATION NODE 07
            stiffness.append(np.abs(reactionX[-1]) / np.abs(disp_mid[-1]))
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
            if ELE_TYPE == 'nonlinearBeamColumn': 
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
                # SECTION NEUTRAL AXIS & CRACK DEPTH
                na, cd = S18.SECTION_CRACK_DEPTH_FUN(center_node-1, SEC_TAG,
                                        -0.15,          # TOP FIBER MATERIAL ULTIMATE STRAIN FOR CRACK DEPTH EVALUTION
                                        -Depth*0.5, 0.0, # BOTTOM FIBER FROM NEUTRAL AXIS
                                        +Depth*0.5, 0.0, # TOP FIBER FROM NEUTRAL AXIS 
                                        )
                NA.append(na); CD.append(cd);
            else:
                STRESSb.append(0.0)
                STRAINb.append(0.0)
                STRESSt.append(0.0)
                STRAINt.append(0.0)
                CURVATURE.append(0.0)
                NA.append(0.0); CD.append(0.0);
            #print(time[-1], dispY[-1], veloY[-1])
            print(f"Time: {time[-1]:.4f}, Displacement: {dispY[-1]:.4f} mm, Reaction: {reactionY[-1]:.2f} N")      
        else:
            print('\n\nFREE-VIBRATION ANALYSIS DONE.\n\n')    
        # Calculating Damping Ratio and Period Using Logarithmic Decrement Analysis 
        damping_ratio = S04.DAMPING_RATIO(dispY)  
        
        # Compute modal properties
        ops.modalProperties("-print", "-file", f"SALAR_ModalReport_{ANAL_TYPE}.txt", "-unorm") 
        
        DATA = (time, reactionX, reactionY, reactionZ, disp_mid,
                ele_axialforce, ele_shearforce, ele_momentforce,
                node_displacementsX, node_displacementsY,
                dispX, dispY,
                veloX, veloY,
                accX, accY,
                stiffness, PERIOD, damping_ratio,
                np.array(PERIOD_MIN), np.array(PERIOD_MAX),
                STRESSb, STRAINb, STRESSt, STRAINt, CURVATURE, NA, CD)
        
        return DATA
    if ANAL_TYPE == 'SEISMIC': # DYNAMIC TIME-HISTORY ANALYSIS
        #%% DEFINE PARAMETERS FOR SEISMIC ANALYSIS
        duration = 20.0                    # [s] Analysis duration
        dt = 0.01                          # [s] Time step
        GMfact = 9810                      # [mm/s²]standard acceleration of gravity or standard acceleration
        SSF_X = 2.50                       # Seismic Acceleration Scale Factor in X Direction
        SSF_Y = 2.50                       # Seismic Acceleration Scale Factor in Y Direction
        iv0_X = 0.0005                     # [mm/s] Initial velocity applied to the node  in X Direction
        iv0_Y = 0.0005                     # [mm/s] Initial velocity applied to the node  in Y Direction
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
            reactionX.append(ops.nodeReaction(1, 1))              # SHEAR BASE REACTION
            reactionY.append(ops.nodeReaction(1, 2))              # AXIAL BASE REACTION
            reactionZ.append(ops.nodeReaction(1, 3))              # MOMENT BASE REACTION
            disp_mid.append(ops.nodeDisp(center_node, 2))         # DISPLACEMENT NODE 07 
            dispX.append(ops.nodeDisp(center_node, 1))            # DISPLACEMENT NODE 07 IN X DIR 
            dispY.append(ops.nodeDisp(center_node, 2))            # DISPLACEMENT NODE 07 IN Y DIR 
            veloX.append(ops.nodeVel(center_node, 1))             # VELOCITY NODE 07
            veloY.append(ops.nodeVel(center_node, 2))             # VELOCITY NODE 07
            accX.append(ops.nodeAccel(center_node, 1))            # ACCELERATION NODE 07
            accY.append(ops.nodeAccel(center_node, 2))            # ACCELERATION NODE 07
            stiffness.append(np.abs(reactionX[-1]) / np.abs(disp_mid[-1]))
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
            if ELE_TYPE == 'nonlinearBeamColumn': 
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
                # SECTION NEUTRAL AXIS & CRACK DEPTH
                na, cd = S18.SECTION_CRACK_DEPTH_FUN(center_node-1, SEC_TAG,
                                        -0.15,          # TOP FIBER MATERIAL ULTIMATE STRAIN FOR CRACK DEPTH EVALUTION
                                        -Depth*0.5, 0.0, # BOTTOM FIBER FROM NEUTRAL AXIS
                                        +Depth*0.5, 0.0, # TOP FIBER FROM NEUTRAL AXIS 
                                        )
                NA.append(na); CD.append(cd);
            else:
                STRESSb.append(0.0)
                STRAINb.append(0.0)
                STRESSt.append(0.0)
                STRAINt.append(0.0)
                CURVATURE.append(0.0)
                NA.append(0.0); CD.append(0.0);
            #print(time[-1], dispY[-1], veloY[-1])
            print(f"Time: {time[-1]:.4f}, Displacement: {dispY[-1]:.4f} mm, Reaction: {reactionY[-1]:.2f} N")

        else:
            print('\n\nSEISMIC ANALYSIS DONE.\n\n')     
        # Calculating Damping Ratio and Period Using Logarithmic Decrement Analysis 
        damping_ratio = S04.DAMPING_RATIO(dispY)  
        
        # Compute modal properties
        ops.modalProperties("-print", "-file", f"SALAR_ModalReport_{ANAL_TYPE}.txt", "-unorm") 
        
        DATA = (time, reactionX, reactionY, reactionZ, disp_mid,
                ele_axialforce, ele_shearforce, ele_momentforce,
                node_displacementsX, node_displacementsY,
                dispX, dispY,
                veloX, veloY,
                accX, accY,
                stiffness, PERIOD, damping_ratio,
                np.array(PERIOD_MIN), np.array(PERIOD_MAX),
                STRESSb, STRAINb, STRESSt, STRAINt, CURVATURE, NA, CD)
        
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
LENGTH = 3000.0       # [mm] Column Length
TOTAL_MASS = 0.0      # [kg] Total Mass of Structure
MODE_NUM = 1
SEMI_RIGID = True
#%%----------------------------------------------------
# PERIOD ANALYSIS
#ELE_TYPE = 'elasticBeamColumn'      # MATERIAL LINEARITY
ELE_TYPE = 'nonlinearBeamColumn'     # MATERIAL AND GEOMETRIC NONLINEARITIES
MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'
ANAL_TYPE = 'PERIOD'

DATA = POST_BUCKLING_COLUMN(LENGTH, MAT_TYPE, ELE_TYPE, ANAL_TYPE, TOTAL_MASS, MODE_NUM, SEMI_RIGID)
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

DATA = POST_BUCKLING_COLUMN(LENGTH, MAT_TYPE, ELE_TYPE, ANAL_TYPE, TOTAL_MASS, MODE_NUM, SEMI_RIGID)
(reactionX, reactionY, reactionZ, disp_mid, 
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

DATA = POST_BUCKLING_COLUMN(LENGTH, MAT_TYPE, ELE_TYPE, ANAL_TYPE, TOTAL_MASS, MODE_NUM, SEMI_RIGID)
(reactionX_PUSH, reactionY_PUSH, reactionZ_PUSH, disp_mid_PUSH, 
 ele_axialforce_PUSH, ele_shearforce_PUSH, ele_momentforce_PUSH,
 node_displacementsX_PUSH, node_displacementsY_PUSH,
  dispX_PUSH, dispY_PUSH,
 PERIOD_MIN_PUSH, PERIOD_MAX_PUSH,
 STRESSb_PUSH, STRAINb_PUSH, STRESSt_PUSH, STRAINt_PUSH, CURVATURE_PUSH, NA_PUSH, CD_PUSH) = DATA


XDATA = disp_mid_PUSH
YDATA = reactionY_PUSH
XLABEL = 'Displacement [mm]'
YLABEL = 'Base Reaction [N]'
TITLE = 'Base Reaction and Displacement of Structure During Pushover Analysis'
COLOR = 'black'
SEMILOGY = False
PLOT(XDATA, YDATA, TITLE, XLABEL, YLABEL, COLOR, SEMILOGY)

DATA = S07.BILNEAR_CURVE(np.abs(disp_mid_PUSH), np.abs(reactionX_PUSH), SLOPE_NODE=10)
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
plt.title('Behavior of element - Stress-strain of Top and Bottom Fibers Curve')
plt.grid(True)
plt.legend()
plt.show()

# Plot Neutral Axis During Analysis
plt.figure(55, figsize=(12, 8))
plt.plot(disp_mid_PUSH, NA_PUSH, linewidth=4, color='purple', label=f"TOP FIBER CRACK DEPTH: {np.abs(NA_PUSH[-1]) : .3f}")
plt.xlabel('Displacment (mm)')
plt.ylabel('Neutral Axis (mm)')
plt.title('Neutral Axis Curve During Analysis')
plt.grid(True)
plt.legend()
plt.show()

# Plot Crack Depth During Analysis
plt.figure(66, figsize=(12, 8))
plt.plot(disp_mid_PUSH, CD_PUSH, linewidth=4, color='black', label=f"NEUTRAL AXIS: {np.abs(CD_PUSH[-1]) : .3f}")
plt.xlabel('Displacment (mm)')
plt.ylabel('Crack Depth (mm)')
plt.title('Crack Depth Curve During Analysis')
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

DATA = POST_BUCKLING_COLUMN(LENGTH, MAT_TYPE, ELE_TYPE, ANAL_TYPE, TOTAL_MASS, MODE_NUM, SEMI_RIGID)
(reactionX_CP, reactionY_CP, reactionZ_CP, disp_mid_CP, 
 ele_axialforce_CP, ele_shearforce_CP, ele_momentforce_CP,
 node_displacementsX_CP, node_displacementsY_CP,
  dispX_CP, dispY_CP,
 PERIOD_MIN_CP, PERIOD_MAX_CP,
 STRESSb_CP, STRAINb_CP, STRESSt_CP, STRAINt_CP, CURVATURE_CP, NA_CP, CD_CP) = DATA


XDATA = disp_mid_CP
YDATA = reactionY_CP
XLABEL = 'Displacement [mm]'
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
plt.title('Behavior of element - Stress-strain of Top and Bottom Fibers Curve')
plt.grid(True)
plt.legend()
plt.show()

# Plot Neutral Axis During Analysis
plt.figure(55, figsize=(12, 8))
plt.plot(disp_mid_CP, NA_CP, linewidth=4, color='purple', label=f"NEUTRAL AXIS: {np.abs(NA_CP[-1]) : .3f}")
plt.xlabel('Displacment (mm)')
plt.ylabel('Neutral Axis (mm)')
plt.title('Neutral Axis Curve During Analysis')
plt.grid(True)
plt.legend()
plt.show()

# Plot Crack Depth During Analysis
plt.figure(66, figsize=(12, 8))
plt.plot(disp_mid_CP, CD_CP, linewidth=4, color='black', label=f"TOP FIBER CRACK DEPTH: {np.abs(CD_CP[-1]) : .3f}")
plt.xlabel('Displacment (mm)')
plt.ylabel('Crack Depth (mm)')
plt.title('Crack Depth Curve During Analysis')
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

DATA = POST_BUCKLING_COLUMN(LENGTH, MAT_TYPE, ELE_TYPE, ANAL_TYPE, TOTAL_MASS, MODE_NUM, SEMI_RIGID)

(reactionX_ETDLS, reactionY_ETDLS, reactionZ_ETDLS, disp_mid_ETDLS, 
 ele_axialforce_ETDLS, ele_shearforce_ETDLS, ele_momentforce_ETDLS,
 node_displacementsX_ETDLS, node_displacementsY_ETDLS,
  dispX_ETDLS, dispY_ETDLS,
 PERIOD_MIN_ETDLS, PERIOD_MAX_ETDLS,
 STRESSb_ETDLS, STRAINb_ETDLS, STRESSt_ETDLS, STRAINt_ETDLS, CURVATURE_ETDLS, NA_ETDLS, CD_ETDLS) = DATA


XDATA = disp_mid_ETDLS
YDATA = reactionY_ETDLS
XLABEL = 'Displacement [mm]'
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
plt.title('Behavior of element - Stress-strain of Top and Bottom Fibers Curve')
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

DATA = POST_BUCKLING_COLUMN(LENGTH, MAT_TYPE, ELE_TYPE, ANAL_TYPE, TOTAL_MASS, MODE_NUM, SEMI_RIGID)

(time_ETDLD, reactionX_ETDLD, reactionY_ETDLD, reactionZ_ETDLD, disp_mid_ETDLD,
ele_axialforce_ETDLD, ele_shearforce_ETDLD, ele_momentforce_ETDLD,
node_displacementsX_ETDLD, node_displacementsY_ETDLD,
dispX_ETDLD, dispY_ETDLD,
veloX_ETDLD, veloY_ETDLD,
accX_ETDLD, accY_ETDLD,
stiffness_ETDLD, PERIOD_ETDLD, damping_ratio_ETDLD,
PERIOD_MIN_ETDLD, PERIOD_MAX_ETDLD,
STRESSb_ETDLD, STRAINb_ETDLD, STRESSt_ETDLD, STRAINt_ETDLD, CURVATURE_ETDLD, CD_ETDLD, CD_ETDLD) = DATA


XDATA = disp_mid_ETDLD
YDATA = reactionY_ETDLD
XLABEL = 'Displacement [mm]'
YLABEL = 'Base Reaction [N]'
TITLE = 'Base Reaction and Displacement of Structure During Dynamic External Time-dependent Loading Analysis'
COLOR = 'black'
SEMILOGY = False
PLOT(XDATA, YDATA, TITLE, XLABEL, YLABEL, COLOR, SEMILOGY)

PLOT_TIME_HISTORY(time_ETDLD, reactionX_ETDLD, disp_mid_ETDLD,
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
YLABEL = 'Element Moment Force [N]'  
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
plt.title('Behavior of element - Stress-strain of Top and Bottom Fibers Curve')
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
DATA = POST_BUCKLING_COLUMN(LENGTH, MAT_TYPE, ELE_TYPE, ANAL_TYPE, TOTAL_MASS, MODE_NUM, SEMI_RIGID)

(time_FV, reactionX_FV, reactionY_FV, reactionZ_FV, disp_mid_FV,
ele_axialforce_FV, ele_shearforce_FV, ele_momentforce_FV,
node_displacementsX_FV, node_displacementsY_FV,
dispX_FV, dispY_FV,
veloX_FV, veloY_FV,
accX_FV, accY_FV,
stiffness_FV, PERIOD_FV, damping_ratio_FV,
PERIOD_MIN_FV, PERIOD_MAX_FV,
STRESSb_FV, STRAINb_FV, STRESSt_FV, STRAINt_FV, CURVATURE_FV, NA_FV, CD_FV) = DATA

XDATA = disp_mid_FV
YDATA = reactionY_FV
XLABEL = 'Displacement [mm]'
YLABEL = 'Base Reaction [N]'
TITLE = 'Base Reaction and Dispalcement of Structure During Free-vibration Analysis'
COLOR = 'black'
SEMILOGY = False
PLOT(XDATA, YDATA, TITLE, XLABEL, YLABEL, COLOR, SEMILOGY)

PLOT_TIME_HISTORY(time_FV, reactionX_FV, disp_mid_FV,
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
YLABEL = 'Element Moment Force [N]'  
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
plt.title('Behavior of element - Stress-strain of Top and Bottom Fibers Curve')
plt.grid(True)
plt.legend()
plt.show()

# Plot Crack Depth During Analysis
plt.figure(55, figsize=(12, 8))
plt.plot(disp_mid_FV, CD_FV, linewidth=4, color='purple', label=f"TOP FIBER CRACK DEPTH: {np.abs(CD_FV[-1]) : .3f}")
plt.xlabel('Displacment (mm)')
plt.ylabel('Crack Depth (mm)')
plt.title('Crack Depth Curve During Analysis')
plt.grid(True)
plt.legend()
plt.show()

# Plot Neutral Axis During Analysis
plt.figure(66, figsize=(12, 8))
plt.plot(disp_mid_FV, NA_FV, linewidth=4, color='black', label=f"NEUTRAL AXIS: {np.abs(NA_FV[-1]) : .3f}")
plt.xlabel('Displacment (mm)')
plt.ylabel('Neutral Axis (mm)')
plt.title('Neutral Axis Curve During Analysis')
plt.grid(True)
plt.legend()
plt.show()

XLABEL = "Displacement Y [mm]"
YLABEL = "Base Axial [N]"
TITLE = "Equivalent viscous damping ratio - Free-vibration Hysteresis"
method = 1 # 
zeta = S055.EQULIVALENT_VISCOUS_DAMPING_RATIO_FUN(dispY_FV, reactionY_FV, method, XLABEL, YLABEL, TITLE)
print(f"Equivalent viscous damping ratio = {zeta:.4f}")

S01.PLOT_2D_FRAME(deformed_scale=1.0)
#%%----------------------------------------------------
# SEISMIC ANALYSIS (DYNAMIC TIME-HISTORY ANALYSIS)
#ELE_TYPE = 'elasticBeamColumn'      # MATERIAL LINEARITY
ELE_TYPE = 'nonlinearBeamColumn'     # MATERIAL AND GEOMETRIC NONLINEARITIES
MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'
ANAL_TYPE = 'SEISMIC'

DATA = POST_BUCKLING_COLUMN(LENGTH, MAT_TYPE, ELE_TYPE, ANAL_TYPE, TOTAL_MASS, MODE_NUM, SEMI_RIGID)

(time_SEI, reactionX_SEI, reactionY_SEI, reactionZ_SEI, disp_mid_SEI,
 ele_axialforce_SEI, ele_shearforce_SEI, ele_momentforce_SEI,
 node_displacementsX_SEI, node_displacementsY_SEI,
 dispX_SEI, dispY_SEI,
 veloX_SEI, veloY_SEI,
 accX_SEI, accY_SEI,
 stiffness_SEI, PERIOD_SEI, damping_ratio_SEI,
 PERIOD_MIN_SEI, PERIOD_MAX_SEI,
 STRESSb_SEI, STRAINb_SEI, STRESSt_SEI, STRAINt_SEI, CURVATURE_SEI, NA_SEI, CD_SEI) = DATA


XDATA = disp_mid_SEI
YDATA = reactionY_SEI
XLABEL = 'Displacement [mm]'
YLABEL = 'Base Reaction [N]'
TITLE = 'Base Reaction and Dispalcement of Structure During Seismic Analysis'
COLOR = 'black'
SEMILOGY = False
PLOT(XDATA, YDATA, TITLE, XLABEL, YLABEL, COLOR, SEMILOGY)

PLOT_TIME_HISTORY(time_SEI, reactionX_SEI, disp_mid_SEI,
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
YLABEL = 'Element Moment Force [N]'  
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
plt.title('Behavior of element - Stress-strain of Top and Bottom Fibers Curve')
plt.grid(True)
plt.legend()
plt.show()

# Plot Crack Depth During Analysis
plt.figure(55, figsize=(12, 8))
plt.plot(disp_mid_SEI, CD_SEI, linewidth=4, color='purple', label=f"TOP FIBER CRACK DEPTH: {np.abs(CD_SEI[-1]) : .3f}")
plt.xlabel('Displacment (mm)')
plt.ylabel('Crack Depth (mm)')
plt.title('Crack Depth Curve During Analysis')
plt.grid(True)
plt.legend()
plt.show()

# Plot Neutral Axis During Analysis
plt.figure(66, figsize=(12, 8))
plt.plot(disp_mid_SEI, NA_SEI, linewidth=4, color='black', label=f"NEUTRAL AXIS: {np.abs(NA_SEI[-1]) : .3f}")
plt.xlabel('Displacment (mm)')
plt.ylabel('Neutral Axis (mm)')
plt.title('Neutral Axis Curve During Analysis')
plt.grid(True)
plt.legend()
plt.show()

S01.PLOT_2D_FRAME(deformed_scale=1.0)
#%%----------------------------------------------------
# --------------------------------------
#  Plot BaseShear-Displacement Analysis 
# --------------------------------------
XX = np.abs(disp_mid_PUSH); YY = np.abs(reactionY_PUSH); # ABSOLUTE VALUE
SLOPE_NODE = 10

DATA = S07.BILNEAR_CURVE(XX, YY, SLOPE_NODE)
X, Y, Elastic_ST, Plastic_ST, Tangent_ST, Ductility_Rito, Over_Strength_Factor = DATA

XLABEL = 'Displacement in X [mm]'
YLABEL = 'Base-Shear Reaction [N]'
LEGEND01 = 'Curve'
LEGEND02 = 'Bilinear Fitted'
LEGEND03 = 'Undefined'
TITLE = f'Last Data of Base Shear-Displacement Analysis - Ductility Ratio: {X[2]/X[1]:.4f} - Over Strength Factor: {Y[2]/Y[1]:.4f}'
COLOR = 'black'
PLOT_2D(np.abs(disp_mid_PUSH), np.abs(reactionY_PUSH), X, Y, X, Y, XLABEL, YLABEL, TITLE, LEGEND01, LEGEND02, LEGEND03, COLOR='black', Z=2) 
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
print(f'Structural Ductility Damage Index in X Direction: {100*DIy:.4f} (%)')
#%%----------------------------------------------------
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
        ax.set_xlabel("Displacement (mm)")
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
        ax.set_xlabel("Displacement (mm)")
        ax.set_ylabel("Base Shear (N)")
        ax.grid(True, linestyle='--', alpha=0.5)
        ax.legend(loc='lower right')
        fig.tight_layout()
    
        return area, fig

Ed_SEI, fig_SEI = DISSIPATED_ENERGY_FUN_WITH_PLOT(
    dispY_SEI, reactionY_SEI, method = 2, 
    title="Earthquake Response – Dissipated Energy"
)
fig_SEI.show()

print(f"Dissipated Energy from Earthquake= {Ed_SEI:.2f} N·m")

Ed_CP, fig_CP = DISSIPATED_ENERGY_FUN_WITH_PLOT(
    dispY_CP, reactionY_CP, method = 2,
    title="Cyclic Loading – Dissipated Energy"
)
fig_CP.show()

print(f"Dissipated Energy from Cyclic Displacement= {Ed_CP:.2f} N·m")


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
    'Minor Damage Level': (20, 40),# Median DI=20%, β=40%
    'Moderate Damage Level': (40, 40),
    'Severe Damage Level': (60, 50),
    'Failure Level': (100, 50)
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

# Add interpretation
results = FF.INTERPRET_FRAGILITY(damage_states, im_values)
print(f"\nProbability of Failure at max displacement: {results['probabilities']['Failure Level']*100:.1f}%")
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
with pd.ExcelWriter("W(x,t)_STEEL_CONCRETE_COLUMN_POST-BUCKLING_8_ANA_OUTPUT.xlsx", engine='openpyxl') as writer:
    
    # PUSHOVER
    df1 = create_df(reactionX_PUSH, disp_mid_PUSH, dispX_PUSH, dispY_PUSH, PERIOD_MIN_PUSH, PERIOD_MAX_PUSH)
    df1.to_excel(writer, sheet_name="PUSHOVER", index=False)
                 
    # CYCLIC DISPLACEMENT
    df1 = create_df(reactionX_CP, disp_mid_CP, dispX_CP, dispY_CP, PERIOD_MIN_CP, PERIOD_MAX_CP)
    df1.to_excel(writer, sheet_name="CYCLIC_DISPLACEMENT", index=False)
    
    # STATIC EXTERNAL TIME-DEPENDENT LOADING
    df2 = create_df(reactionX_ETDLS, disp_mid_ETDLS, dispX_ETDLS, dispY_ETDLS, PERIOD_MIN_ETDLS, PERIOD_MAX_ETDLS)
    df2.to_excel(writer, sheet_name="STATIC_EXTERNAL_TIME-DEPENDENT_LOADING", index=False)

    # DYNAMIC EXTERNAL TIME-DEPENDENT LOADING
    df3 = create_df(reactionX_ETDLD, disp_mid_ETDLD, dispX_ETDLD, dispY_ETDLD, PERIOD_MIN_ETDLD, PERIOD_MAX_ETDLD)
    df3.to_excel(writer, sheet_name="DYNAMIC_EXTERNAL_TIME-DEPENDENT_LOADING", index=False)
    
    # FREE-VIBRATION
    df3 = create_df(reactionX_FV, disp_mid_FV, dispX_FV, dispY_FV, PERIOD_MIN_FV, PERIOD_MAX_FV)
    df3.to_excel(writer, sheet_name="FREE-VIBRATION", index=False)

    # SEISMIC
    df4 = create_df(reactionX_SEI, disp_mid_SEI, dispX_SEI, dispY_SEI, PERIOD_MIN_SEI, PERIOD_MAX_SEI)
    df4.to_excel(writer, sheet_name="SEISMIC", index=False)
#%%----------------------------------------------------