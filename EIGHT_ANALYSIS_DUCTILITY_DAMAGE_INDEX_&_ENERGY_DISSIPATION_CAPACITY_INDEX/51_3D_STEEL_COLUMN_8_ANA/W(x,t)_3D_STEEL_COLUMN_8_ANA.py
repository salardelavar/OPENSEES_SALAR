###########################################################################################################
#                   >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                      #
#           COMPREHENSIVE NONLINEAR SEISMIC ASSESSMENT OF A 3D STEEL COLUMN : AN OPENSEES                 #
#FRAMEWORK FOR STATIC PUSHOVER, CYCLIC DEGRADATION, STATIC TIME-HISTORY AND DYNAMIC TIME-HISTORY ANALYSIS #
#---------------------------------------------------------------------------------------------------------#
#                                                 W(x,t) = W0 sin(wt)                                     #
#                                          W(x,t) = W0 exp(-0.05wt) sin(wt)                               #
#---------------------------------------------------------------------------------------------------------#
#         ASSESSMENT OF DUCTILITY DAMAGE INDICES FOR STRUCTURAL ELEMENTS AND SYSTEMS AND EVALUATION       #
#                                        OF ENERGY DISSIPATION CAPACITY INDEX                             #
#---------------------------------------------------------------------------------------------------------#
#                       THIS PYTHON SCRIPT WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                  #
#                                   EMAIL: salar.d.ghashghaei@gmail.com                                   #
###########################################################################################################
"""
Nonlinear Seismic Performance Assessment of a steel column:
An OpenSeesPy Framework for Material and Geometric Nonlinearity Under Static, Cyclic, and Earthquake Loading

This OpenSeesPy script performs rigorous nonlinear static and dynamic analysis of a steel column
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
#%%----------------------------------------------------
import numpy as np
import openseespy.opensees as ops
import matplotlib.pyplot as plt 
import PLOT_3D as S01
import ANALYSIS_FUNCTION as S02
import PERIOD_FUN as S03
import DAMPING_RATIO_FUN as S04
import EIGENVALUE_ANALYSIS_FUN as S05
import RAYLEIGH_DAMPING_FUN as S06
import BILINEAR_CURVE as S07
import SECTION_CURVATURE_FUN as S09
import SECTION_STRESS_STRAIN_FUN as S10
import I_STEEL_PLATE_SECTION_DIFF_WIDTH_HEIGHT_FUN_EXTRA_3D as S11
import PLATE_DOUBLE_I_STEEL_SECTION_DIFF_WIDTH_HEIGHT_FUN_EXTRA_3D as S12
import PLATE_I_STEEL_SECTION_DIFF_WIDTH_HEIGHT_FUN_EXTRA_3D as S13
import STEEL_BOX_SECTION_DIFF_WIDTH_HEIGHT_FUN_EXTRA_3D as S14
import STEEL_CORNER_BOX_SECTION_DIFF_WIDTH_HEIGHT_FUN_EXTRA_3D as S15
import STEEL_UNP_BOX_SECTION_DIFF_WIDTH_HEIGHT_FUN_EXTRA_3D as S16
import SECTION_CRACK_DEPTH_FUN as S18
import FRAGILITY_ANALYSIS_TWO_PARAMETERS_FUN as S19
import DISSIPATED_ENERGY_WITH_PLOT_FUN as S25
#%%----------------------------------------------------
def STEEL_COLUMN_3D(LENGTH, MAT_TYPE, ELE_TYPE, ANAL_TYPE, TOTAL_MASS):
    # Initialize model
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 6)
    
    MAX_ITERATIONS = 5000   # Maximum number of iterations
    MAX_TOLERANCE = 1.0e-6  # Specified tolerance for convergence
    
    #%% GEOMETRY DEFINITION
    ops.node(1, 0.0, 0.0, (0/6)*LENGTH)
    ops.node(2, 0.0, 0.0, (1/6)*LENGTH)
    ops.node(3, 0.0, 0.0, (2/6)*LENGTH)
    ops.node(4, 0.0, 0.0, (3/6)*LENGTH)
    ops.node(5, 0.0, 0.0, (4/6)*LENGTH)
    ops.node(6, 0.0, 0.0, (5/6)*LENGTH)
    ops.node(7, 0.0, 0.0, (6/6)*LENGTH)
    #%% SUPPORTS
    ops.fix(1, 1, 1, 1, 1, 1, 1)      # Bot  - fixed
    ops.fix(7, 0, 0, 0, 1, 1, 0)      # Top  - rotation fixed
    
    #%% STEEL SECTION
    STEEL_TYPE = 'INELASTIC'
    STEEL_DENSITY = 7850/1e9      # [kg/m^3] -> [kg/mm^3] Steel Material Density
    SEC_TAG = 1
    #Depth, Ele_Mass = S11.I_STEEL_PLATE_SECTION_DIFF_WIDTH_HEIGHT_FUN_EXTRA_3D(SEC_TAG, STEEL_TYPE, STEEL_DENSITY, plot=True)
    #Depth, Ele_Mass = S12.PLATE_DOUBLE_I_STEEL_SECTION_DIFF_WIDTH_HEIGHT_FUN_EXTRA_3D(SEC_TAG, STEEL_TYPE, STEEL_DENSITY, plot=True)
    #Depth, Ele_Mass = S13.PLATE_I_STEEL_SECTION_DIFF_WIDTH_HEIGHT_FUN_EXTRA_3D(SEC_TAG, STEEL_TYPE, STEEL_DENSITY, plot=True)
    #Depth, Ele_Mass = S14.STEEL_BOX_SECTION_DIFF_WIDTH_HEIGHT_FUN_EXTRA_3D(SEC_TAG, STEEL_TYPE, STEEL_DENSITY, plot=True)
    Depth, Ele_Mass = S15.STEEL_CORNER_BOX_SECTION_DIFF_WIDTH_HEIGHT_FUN_EXTRA_3D(SEC_TAG, STEEL_TYPE, STEEL_DENSITY, plot=True)
    #Depth, Ele_Mass = S16.STEEL_UNP_BOX_SECTION_DIFF_WIDTH_HEIGHT_FUN_EXTRA_3D(SEC_TAG, STEEL_TYPE, STEEL_DENSITY, plot=True)
    
    #%% CREATE ELEMENTS
    ele_tag = 1
    transfTag = 1
    ops.geomTransf('Linear', transfTag, 1, 1, 0)
    #ops.geomTransf('PDelta', transfTag, 1, 1, 0)
    #ops.geomTransf('Corotational', transfTag, 1, 1, 0)
    """
    X, Y, and Z components of vecxz, the vector used to define the local x-z plane of the local-coordinate
    system. The local y-axis is defined by taking the cross product of the vecxz vector and the x-axis.
    These components are specified in the global-coordinate system X,Y,Z and define a vector that is in
    a plane parallel to the x-z plane of the local-coordinate system. These items need to be specified for
    the three-dimensional problem.
    """
    numIntgrPts = 5
    if ELE_TYPE == 'elasticBeamColumn':
        ops.element('elasticBeamColumn', 1, 1, 2, SEC_TAG, transfTag,'-mass', Ele_Mass)
        ops.element('elasticBeamColumn', 2, 2, 3, SEC_TAG, transfTag,'-mass', Ele_Mass)
        ops.element('elasticBeamColumn', 3, 3, 4, SEC_TAG, transfTag,'-mass', Ele_Mass)
        ops.element('elasticBeamColumn', 4, 4, 5, SEC_TAG, transfTag,'-mass', Ele_Mass)
        ops.element('elasticBeamColumn', 5, 5, 6, SEC_TAG, transfTag,'-mass', Ele_Mass)
        ops.element('elasticBeamColumn', 6, 6, 7, SEC_TAG, transfTag,'-mass', Ele_Mass)
    if ELE_TYPE == 'nonlinearBeamColumn': 
        ops.element('nonlinearBeamColumn', 1, 1, 2, numIntgrPts, SEC_TAG, transfTag,'-mass', Ele_Mass) 
        ops.element('nonlinearBeamColumn', 2, 2, 3, numIntgrPts, SEC_TAG, transfTag,'-mass', Ele_Mass) 
        ops.element('nonlinearBeamColumn', 3, 3, 4, numIntgrPts, SEC_TAG, transfTag,'-mass', Ele_Mass)
        ops.element('nonlinearBeamColumn', 4, 4, 5, numIntgrPts, SEC_TAG, transfTag,'-mass', Ele_Mass) 
        ops.element('nonlinearBeamColumn', 5, 5, 6, numIntgrPts, SEC_TAG, transfTag,'-mass', Ele_Mass) 
        ops.element('nonlinearBeamColumn', 6, 6, 7, numIntgrPts, SEC_TAG, transfTag,'-mass', Ele_Mass) 
    if ELE_TYPE == 'dispBeamColumnInt': 
        # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/dispBeamColumnInt.html
        cRot = 0.40 # Fraction of the height is distance from bottom to the center of rotation (0 to 1)
        ops.element('dispBeamColumnInt', 1, 1, 2, numIntgrPts, SEC_TAG, transfTag, cRot, '-mass', Ele_Mass) 
        ops.element('dispBeamColumnInt', 2, 2, 3, numIntgrPts, SEC_TAG, transfTag, cRot, '-mass', Ele_Mass) 
        ops.element('dispBeamColumnInt', 3, 3, 4, numIntgrPts, SEC_TAG, transfTag, cRot, '-mass', Ele_Mass)
        ops.element('dispBeamColumnInt', 4, 4, 5, numIntgrPts, SEC_TAG, transfTag, cRot, '-mass', Ele_Mass) 
        ops.element('dispBeamColumnInt', 5, 5, 6, numIntgrPts, SEC_TAG, transfTag, cRot, '-mass', Ele_Mass) 
        ops.element('dispBeamColumnInt', 6, 6, 7, numIntgrPts, SEC_TAG, transfTag, cRot, '-mass', Ele_Mass)
    ele_tag = 6+1
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
    center_node = 7
    TS_TAG = 1
    PATT_TAG = 1
    ops.timeSeries('Linear', TS_TAG)
    ops.pattern('Plain', PATT_TAG, TS_TAG)
    for ii in range(1, tag):
        ops.load(ii, 0.0, 0.0, -ENW, 0.0, 0.0, 0.0)   # [N] downward
        ops.mass(ii, ENM, ENM, ENM, 0.0, 0.0, 0.0)
        
    #ops.mass(7, ENM, ENM, 0.0)
    TS_TAG = 2
    PATT_TAG = 2
    ops.timeSeries('Linear', TS_TAG)
    ops.pattern('Plain', PATT_TAG, TS_TAG)       
    #for ii in range(1, tag-1):
    #    ops.eleLoad('-ele', ii,'-type', '-beamUniform', -0.0003, 0.0) # [N/mm] Uniformly-distributed load  
    ops.load(center_node, 1.0, 1.0, 0.0,
             0.0, 0.0, 0.0)
    
    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.system('BandGeneral')
    
    if ELE_TYPE == 'ELASTIC':
        ops.algorithm('Linear')
    if ELE_TYPE == 'INELASTIC':
        ops.algorithm('Newton')
        
    time = []
    dispX, dispY, dispZ = [], [], []
    veloX, veloY, veloZ = [], [], []
    accX, accY, accZ = [], [], []
    reactionX, reactionY, reactionZ = [], [], []
    reactionMX, reactionMY, reactionMZ = [], [], []
    disp_mid = []
    stiffness = []
    OMEGA, PERIOD = [], []
    PERIOD_MIN, PERIOD_MAX = [], []
    STRESSb = []; STRAINb = []; 
    STRESSt = []; STRAINt = [];
    CURVATURE = [];
    NA, CD = [], [];
    # Initialize lists for each element's axial force, shear force and SHEAR FORCE
    ele_axialforce = {
        1: [],  # AXIALFORCE-01
        2: [],  # AXIALFORCE-02
        3: [],  # AXIALFORCE-03
        4: [],  # AXIALFORCE-04
        5: [],  # AXIALFORCE-05
        6: [],  # AXIALFORCE-06
        }  
    ele_shearforceX = {
        1: [],  # SHEARFORCE-01
        2: [],  # SHEARFORCE-02
        3: [],  # SHEARFORCE-03
        4: [],  # SHEARFORCE-04
        5: [],  # SHEARFORCE-05
        6: [],  # SHEARFORCE-06
        }  
    ele_shearforceY = {
        1: [],  # SHEARFORCE-01
        2: [],  # SHEARFORCE-02
        3: [],  # SHEARFORCE-03
        4: [],  # SHEARFORCE-04
        5: [],  # SHEARFORCE-05
        6: [],  # SHEARFORCE-06
        }     
    # Initialize lists for each node's displacement
    node_displacementsX = {
        1: [],  # DISP01
        2: [],  # DISP02
        3: [],  # DISP03
        4: [],  # DISP04
        5: [],  # DISP05
        6: [],  # DISP06
        7: [],  # DISP07
        }
    node_displacementsY = {
        1: [],  # DISP01
        2: [],  # DISP02
        3: [],  # DISP03
        4: [],  # DISP04
        5: [],  # DISP05
        6: [],  # DISP06
        7: [],  # DISP07
        }
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
        reactionX.append(ops.nodeReaction(1, 1))                        # SHEAR BASE REACTION IN X
        reactionY.append(ops.nodeReaction(1, 2))                        # SHEAR BASE REACTION IN Y
        reactionZ.append(ops.nodeReaction(1, 3))                        # AXIAL BASE REACTION IN Z
        reactionMX.append(ops.nodeReaction(1, 4))                       # MOMENT REACTION IN X
        reactionMY.append(ops.nodeReaction(1, 5))                       # MOMENT BASE REACTION IN Y
        reactionMZ.append(ops.nodeReaction(1, 6))                       # TORSION BASE REACTION IN Z
        disp_mid.append(ops.nodeDisp(center_node, 1))                   # DISPLACEMENT NODE 07 
        dispX.append(ops.nodeDisp(center_node, 1))                      # DISPLACEMENT NODE 07 IN X DIR 
        dispY.append(ops.nodeDisp(center_node, 2))                      # DISPLACEMENT NODE 07 IN Y DIR 
        dispZ.append(ops.nodeDisp(center_node, 3))                      # DISPLACEMENT NODE 07 IN Z DIR
        # Store axial forces, strain and stress
        for ele_id in ele_axialforce.keys(): 
            ele_axialforce[ele_id].append(ops.eleResponse(ele_id, 'force')[0])        # [N] ELEMENT AXIAL FORCE
            ele_shearforceX[ele_id].append(ops.eleResponse(ele_id, 'force')[1])        # [N] ELEMENT SHEAR FORCE
            ele_shearforceY[ele_id].append(ops.eleResponse(ele_id, 'force')[2])       # [N] ELEMENT SHEAR FORCE                
        # Store displacements
        for node_id in node_displacementsX.keys():    
            node_displacementsX[node_id].append(ops.nodeDisp(node_id, 1))
            node_displacementsY[node_id].append(ops.nodeDisp(node_id, 2))
        else:
            print('\n\nSTATIC ANALYSIS DONE.\n\n')  
            
        DATA = (reactionX, reactionY, reactionZ,
                reactionMX, reactionMY, reactionMZ,
                disp_mid,
                ele_axialforce, ele_shearforceX, ele_shearforceY,
                node_displacementsX, node_displacementsY,
                dispX, dispY)
        
        return  DATA
    
    if ANAL_TYPE == 'PUSHOVER': # STATIC TIME-HISTORY ANALYSIS
        IDctrlDOF = 1   # 1: Horizental Dispalcement - 2: Vertical Dispalcement
        DINCR = -0.01   # [mm] Incremental Displacement
        DMAX = -300.0   # [mm] Max. Displacement
        ops.integrator('DisplacementControl', center_node, IDctrlDOF, DINCR)
        ops.integrator('DisplacementControl', center_node, 2, DINCR)
        ops.analysis('Static')
        Nsteps =  int(np.abs(DMAX/ DINCR)) 
        STEP = 0.0
        for step in range(Nsteps):
            OK = ops.analyze(1)
            S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS)
            ops.reactions()
            reactionX.append(ops.nodeReaction(1, 1))                        # SHEAR BASE REACTION IN X
            reactionY.append(ops.nodeReaction(1, 2))                        # SHEAR BASE REACTION IN Y
            reactionZ.append(ops.nodeReaction(1, 3))                        # AXIAL BASE REACTION IN Z
            reactionMX.append(ops.nodeReaction(1, 4))                       # MOMENT REACTION IN X
            reactionMY.append(ops.nodeReaction(1, 5))                       # MOMENT BASE REACTION IN Y
            reactionMZ.append(ops.nodeReaction(1, 6))                       # TORSION BASE REACTION IN Z
            disp_mid.append(ops.nodeDisp(center_node, 1))                   # DISPLACEMENT NODE 07 
            dispX.append(ops.nodeDisp(center_node, 1))                      # DISPLACEMENT NODE 07 IN X DIR 
            dispY.append(ops.nodeDisp(center_node, 2))                      # DISPLACEMENT NODE 07 IN Y DIR 
            dispZ.append(ops.nodeDisp(center_node, 3))                      # DISPLACEMENT NODE 07 IN Z DIR
            # Store elements force, strain and stress
            for ele_id in ele_axialforce.keys(): 
                ele_axialforce[ele_id].append(ops.eleResponse(ele_id, 'force')[0])        # [N] ELEMENT AXIAL FORCE
                ele_shearforceX[ele_id].append(ops.eleResponse(ele_id, 'force')[1])        # [N] ELEMENT SHEAR FORCE
                ele_shearforceY[ele_id].append(ops.eleResponse(ele_id, 'force')[2])       # [N] ELEMENT SHEAR FORCE
            # Store displacements
            for node_id in node_displacementsX.keys():    
                node_displacementsX[node_id].append(ops.nodeDisp(node_id, 1))
                node_displacementsY[node_id].append(ops.nodeDisp(node_id, 2))      
            # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
            #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(3, 0.5*DR, DR, 0, 1)
            PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(3, PLOT=True)
            PERIOD_MIN.append(PERIODmin)
            PERIOD_MAX.append(PERIODmax) 
            if ELE_TYPE == 'nonlinearBeamColumn' or ELE_TYPE == 'dispBeamColumnInt':
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
                                        -240.0/200000.0,          # TOP FIBER MATERIAL ULTIMATE STRAIN FOR CRACK DEPTH EVALUTION
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
            print(f"Step: {STEP}, Displacement: {dispX[-1]:.4f} mm, Reaction: {reactionX[-1]:.2f} N")     
        else:
            print('\n\nPUSHOVER ANALYSIS DONE.\n\n')
            
        DATA = (reactionX, reactionY, reactionZ,
                reactionMX, reactionMY, reactionMZ,
                disp_mid,
                ele_axialforce, ele_shearforceX, ele_shearforceY,
                node_displacementsX, node_displacementsY,
                dispX, dispY,
                np.array(PERIOD_MIN), np.array(PERIOD_MAX),
                STRESSb, STRAINb, STRESSt, STRAINt, CURVATURE, NA, CD)
    
        return  DATA
    
    if ANAL_TYPE == 'CYCLIC_DISPLACEMENT': # STATIC TIME-HISTORY ANALYSIS
        IDctrlDOF = 1   # 1: Horizental Dispalcement - 2: Vertical Dispalcement
        DMAX = -300.0   # [mm] Max. Displacement
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
            current_disp = ops.nodeDisp(center_node, 1) # DISPALCEMENT APPLIED IN NODE 
            dU = target_disp - current_disp
            ops.integrator('DisplacementControl', center_node, IDctrlDOF, dU)
            OK = ops.analyze(1)
            S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS)
            ops.reactions()
            reactionX.append(ops.nodeReaction(1, 1))                        # SHEAR BASE REACTION IN X
            reactionY.append(ops.nodeReaction(1, 2))                        # SHEAR BASE REACTION IN Y
            reactionZ.append(ops.nodeReaction(1, 3))                        # AXIAL BASE REACTION IN Z
            reactionMX.append(ops.nodeReaction(1, 4))                       # MOMENT REACTION IN X
            reactionMY.append(ops.nodeReaction(1, 5))                       # MOMENT BASE REACTION IN Y
            reactionMZ.append(ops.nodeReaction(1, 6))                       # TORSION BASE REACTION IN Z                  # DISPLACEMENT NODE 07 
            disp_mid.append(ops.nodeDisp(center_node, 1))                   # DISPLACEMENT NODE 07
            dispX.append(ops.nodeDisp(center_node, 1))                      # DISPLACEMENT NODE 07 IN X DIR 
            dispY.append(ops.nodeDisp(center_node, 2))                      # DISPLACEMENT NODE 07 IN Y DIR 
            dispZ.append(ops.nodeDisp(center_node, 3))                      # DISPLACEMENT NODE 07 IN Z DIR
            # Store elements force, strain and stress
            for ele_id in ele_axialforce.keys(): 
                ele_axialforce[ele_id].append(ops.eleResponse(ele_id, 'force')[0])        # [N] ELEMENT AXIAL FORCE
                ele_shearforceX[ele_id].append(ops.eleResponse(ele_id, 'force')[1])        # [N] ELEMENT SHEAR FORCE
                ele_shearforceY[ele_id].append(ops.eleResponse(ele_id, 'force')[2])       # [N] ELEMENT SHEAR FORCE
            # Store displacements
            for node_id in node_displacementsX.keys():    
                node_displacementsX[node_id].append(ops.nodeDisp(node_id, 1))
                node_displacementsY[node_id].append(ops.nodeDisp(node_id, 2))      
            # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
            #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(3, 0.5*DR, DR, 0, 1)
            PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(3, PLOT=True)
            PERIOD_MIN.append(PERIODmin)
            PERIOD_MAX.append(PERIODmax) 
            if ELE_TYPE == 'nonlinearBeamColumn' or ELE_TYPE == 'dispBeamColumnInt':
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
                                        -240.0/200000.0,          # TOP FIBER MATERIAL ULTIMATE STRAIN FOR CRACK DEPTH EVALUTION
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
            print(f"Step: {STEP}, Displacement: {dispX[-1]:.4f} mm, Reaction: {reactionX[-1]:.2f} N")     
        else:
            print('\n\nCYCLIC DISPLAEMENT ANALYSIS DONE.\n\n')
            
        DATA = (reactionX, reactionY, reactionZ,
                reactionMX, reactionMY, reactionMZ,
                disp_mid,
                ele_axialforce, ele_shearforceX, ele_shearforceY,
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
        #ops.load(center_node, 1.0, 0.0, 0.0)
        for ii in range(1, tag-1):# W(x,t) = W0 exp(-0.05wt) sin(wt) 
            ops.eleLoad('-ele', ii,'-type', '-beamUniform', -0.0003, 0.0) # [N/mm] Uniformly-distributed load
        
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
            reactionX.append(ops.nodeReaction(1, 1))                        # SHEAR BASE REACTION IN X
            reactionY.append(ops.nodeReaction(1, 2))                        # SHEAR BASE REACTION IN Y
            reactionZ.append(ops.nodeReaction(1, 3))                        # AXIAL BASE REACTION IN Z
            reactionMX.append(ops.nodeReaction(1, 4))                       # MOMENT REACTION IN X
            reactionMY.append(ops.nodeReaction(1, 5))                       # MOMENT BASE REACTION IN Y
            reactionMZ.append(ops.nodeReaction(1, 6))                       # TORSION BASE REACTION IN Z
            disp_mid.append(ops.nodeDisp(center_node, 1))                   # DISPLACEMENT NODE 07 
            dispX.append(ops.nodeDisp(center_node, 1))                      # DISPLACEMENT NODE 07 IN X DIR 
            dispY.append(ops.nodeDisp(center_node, 2))                      # DISPLACEMENT NODE 07 IN Y DIR 
            dispZ.append(ops.nodeDisp(center_node, 3))                      # DISPLACEMENT NODE 07 IN Z DIR
            # Store elements force, strain and stress
            for ele_id in ele_axialforce.keys(): 
                ele_axialforce[ele_id].append(ops.eleResponse(ele_id, 'force')[0])        # [N] ELEMENT AXIAL FORCE
                ele_shearforceX[ele_id].append(ops.eleResponse(ele_id, 'force')[1])        # [N] ELEMENT SHEAR FORCE
                ele_shearforceY[ele_id].append(ops.eleResponse(ele_id, 'force')[2])       # [N] ELEMENT SHEAR FORCE
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
                                        -240.0/200000.0,          # TOP FIBER MATERIAL ULTIMATE STRAIN FOR CRACK DEPTH EVALUTION
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
            print(f"Step: {STEP}, Displacement: {dispX[-1]:.4f} mm, Reaction: {reactionX[-1]:.2f} N")     
        else:
            print('\n\nSTATIC EXTERNAL TIME-DEPENDENT LOADING ANALYSIS DONE.\n\n')
            
        DATA = (reactionX, reactionY, reactionZ,
                reactionMX, reactionMY, reactionMZ,
                disp_mid,
                ele_axialforce, ele_shearforceX, ele_shearforceY,
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
        #ops.load(center_node, 1.0, 0.0, 0.0)
        for ii in range(1, tag-1):# W(x,t) = W0 exp(-0.05wt) sin(wt) 
            ops.eleLoad('-ele', ii,'-type', '-beamUniform', -0.0003, 0.0) # [N/mm] Uniformly-distributed load
        
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
            reactionX.append(ops.nodeReaction(1, 1))              # SHEAR BASE REACTION IN X
            reactionY.append(ops.nodeReaction(1, 2))              # SHEAR BASE REACTION IN Y
            reactionZ.append(ops.nodeReaction(1, 3))              # AXIAL BASE REACTION IN Z
            reactionMX.append(ops.nodeReaction(1, 4))             # MOMENT REACTION IN X
            reactionMY.append(ops.nodeReaction(1, 5))             # MOMENT BASE REACTION IN Y
            reactionMZ.append(ops.nodeReaction(1, 6))             # TORSION BASE REACTION IN Z
            disp_mid.append(ops.nodeDisp(center_node, 1))         # DISPLACEMENT NODE 07 
            dispX.append(ops.nodeDisp(center_node, 1))            # DISPLACEMENT NODE 07 IN X DIR 
            dispY.append(ops.nodeDisp(center_node, 2))            # DISPLACEMENT NODE 07 IN Y DIR 
            dispZ.append(ops.nodeDisp(center_node, 3))            # DISPLACEMENT NODE 07 IN Z DIR
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
                ele_shearforceX[ele_id].append(ops.eleResponse(ele_id, 'force')[1])        # [N] ELEMENT SHEAR FORCE
                ele_shearforceY[ele_id].append(ops.eleResponse(ele_id, 'force')[2])       # [N] ELEMENT SHEAR FORCE
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
                                        -240.0/200000.0,          # TOP FIBER MATERIAL ULTIMATE STRAIN FOR CRACK DEPTH EVALUTION
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
            print(f"Time: {time[-1]:.4f}, Displacement: {dispX[-1]:.4f} mm, Reaction: {reactionX[-1]:.2f} N")      
        else:
            print('\n\nDYNAMIC EXTERNAL TIME-DEPENDENT LOADING ANALYSIS DONE.\n\n')  
        # Calculating Damping Ratio and Period Using Logarithmic Decrement Analysis 
        damping_ratio = S04.DAMPING_RATIO(dispX)  
        
        # Compute modal properties
        ops.modalProperties("-print", "-file", f"SALAR_ModalReport_{ANAL_TYPE}.txt", "-unorm") 
        
        DATA = (time,
                reactionX, reactionY, reactionZ,
                reactionMX, reactionMY, reactionMZ,
                disp_mid,
                ele_axialforce, ele_shearforceX, ele_shearforceY,
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
        u0 = -1.35                         # [mm] Initial displacement
        v0 = 0.15                          # [mm/s] Initial velocity
        a0 = 0.065                         # [mm/s^2] Initial acceleration
        IU = True                          # Free Vibration with Initial Displacement
        IV = True                          # Free Vibration with Initial Velocity
        IA = True                          # Free Vibration with Initial Acceleration
        duration = 20.0                    # [s] Analysis duration
        dt = 0.001                         # [s] Time step
        DR = 0.03                          # Damping Ratio
        
        if IU == True:
            # Define initial displacment
            ops.setNodeDisp(center_node, 1, u0, '-commit') # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/setNodeDisp.html
            ops.setNodeDisp(center_node, 2, u0, '-commit')
        if IV == True:
            # Define initial velocity
            ops.setNodeVel(center_node, 1, v0, '-commit')  # INFO LINK: https://openseespydoc.readthedocs.io/en/stable/src/setNodeVel.html
            ops.setNodeVel(center_node, 2, v0, '-commit')
        if IA == True:
            # Define initial  acceleration
            ops.setNodeAccel(center_node, 1, a0, '-commit') # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/setNodeAccel.html
            ops.setNodeAccel(center_node, 2, a0, '-commit')
            
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
            reactionX.append(ops.nodeReaction(1, 1))              # SHEAR BASE REACTION IN X
            reactionY.append(ops.nodeReaction(1, 2))              # SHEAR BASE REACTION IN Y
            reactionZ.append(ops.nodeReaction(1, 3))              # AXIAL BASE REACTION IN Z
            reactionMX.append(ops.nodeReaction(1, 4))             # MOMENT REACTION IN X
            reactionMY.append(ops.nodeReaction(1, 5))             # MOMENT BASE REACTION IN Y
            reactionMZ.append(ops.nodeReaction(1, 6))             # TORSION BASE REACTION IN Z
            disp_mid.append(ops.nodeDisp(center_node, 1))         # DISPLACEMENT NODE 07 
            dispX.append(ops.nodeDisp(center_node, 1))            # DISPLACEMENT NODE 07 IN X DIR 
            dispY.append(ops.nodeDisp(center_node, 2))            # DISPLACEMENT NODE 07 IN Y DIR 
            dispZ.append(ops.nodeDisp(center_node, 3))            # DISPLACEMENT NODE 07 IN Z DIR
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
                ele_shearforceX[ele_id].append(ops.eleResponse(ele_id, 'force')[1])        # [N] ELEMENT SHEAR FORCE
                ele_shearforceY[ele_id].append(ops.eleResponse(ele_id, 'force')[2])       # [N] ELEMENT SHEAR FORCE
            # Store displacements
            for node_id in node_displacementsX.keys():    
                node_displacementsX[node_id].append(ops.nodeDisp(node_id, 1))
                node_displacementsY[node_id].append(ops.nodeDisp(node_id, 2))   
            # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
            #PERIODmin, PERIODmax = S06.RAYLEIGH_DAMPING(3, 0.5*DR, DR, 0, 1)
            PERIODmin, PERIODmax = S05.EIGENVALUE_ANALYSIS(3, PLOT=True)
            PERIOD_MIN.append(PERIODmin)
            PERIOD_MAX.append(PERIODmax)
            if ELE_TYPE == 'nonlinearBeamColumn' or ELE_TYPE == 'dispBeamColumnInt':
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
                                        -240.0/200000.0,          # TOP FIBER MATERIAL ULTIMATE STRAIN FOR CRACK DEPTH EVALUTION
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
            print(f"Time: {time[-1]:.4f}, Displacement: {dispX[-1]:.4f} mm, Reaction: {reactionX[-1]:.2f} N")      
        else:
            print('\n\nFREE-VIBRATION ANALYSIS DONE.\n\n')    
        # Calculating Damping Ratio and Period Using Logarithmic Decrement Analysis 
        damping_ratioX = S04.DAMPING_RATIO(dispX) 
        print('Damping Ratio in X Dir.: {damping_ratioX:.3f} [%]')
        damping_ratioY = S04.DAMPING_RATIO(dispY)
        print('Damping Ratio in Y Dir.: {damping_ratioY:.3f} [%]')
        
        # Compute modal properties
        ops.modalProperties("-print", "-file", f"SALAR_ModalReport_{ANAL_TYPE}.txt", "-unorm") 
        
        DATA = (time,
                reactionX, reactionY, reactionZ,
                reactionMX, reactionMY, reactionMZ,
                disp_mid,
                ele_axialforce, ele_shearforceX, ele_shearforceY,
                node_displacementsX, node_displacementsY,
                dispX, dispY,
                veloX, veloY,
                accX, accY,
                stiffness, PERIOD, damping_ratioX,
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
        SEI = 'XY'                         # Seismic Direction
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
            reactionX.append(ops.nodeReaction(1, 1))              # SHEAR BASE REACTION IN X
            reactionY.append(ops.nodeReaction(1, 2))              # SHEAR BASE REACTION IN Y
            reactionZ.append(ops.nodeReaction(1, 3))              # AXIAL BASE REACTION IN Z
            reactionMX.append(ops.nodeReaction(1, 4))             # MOMENT REACTION IN X
            reactionMY.append(ops.nodeReaction(1, 5))             # MOMENT BASE REACTION IN Y
            reactionMZ.append(ops.nodeReaction(1, 6))             # TORSION BASE REACTION IN Z
            disp_mid.append(ops.nodeDisp(center_node, 1))         # DISPLACEMENT NODE 07 
            dispX.append(ops.nodeDisp(center_node, 1))            # DISPLACEMENT NODE 07 IN X DIR 
            dispY.append(ops.nodeDisp(center_node, 2))            # DISPLACEMENT NODE 07 IN Y DIR 
            dispZ.append(ops.nodeDisp(center_node, 3))            # DISPLACEMENT NODE 07 IN Z DIR
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
                ele_shearforceX[ele_id].append(ops.eleResponse(ele_id, 'force')[1])        # [N] ELEMENT SHEAR FORCE
                ele_shearforceY[ele_id].append(ops.eleResponse(ele_id, 'force')[2])       # [N] ELEMENT SHEAR FORCE
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
                                        -240.0/200000.0,          # TOP FIBER MATERIAL ULTIMATE STRAIN FOR CRACK DEPTH EVALUTION
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
            print(f"Time: {time[-1]:.4f}, Displacement: {dispX[-1]:.4f} mm, Reaction: {reactionX[-1]:.2f} N")

        else:
            print('\n\nSEISMIC ANALYSIS DONE.\n\n')     
        # Calculating Damping Ratio and Period Using Logarithmic Decrement Analysis 
        damping_ratio = S04.DAMPING_RATIO(dispX)  
        
        # Compute modal properties
        ops.modalProperties("-print", "-file", f"SALAR_ModalReport_{ANAL_TYPE}.txt", "-unorm") 
        
        DATA = (time,
                reactionX, reactionY, reactionZ,
                reactionMX, reactionMY, reactionMZ,
                disp_mid,
                ele_axialforce, ele_shearforceX, ele_shearforceY,
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
    plt.show()
    
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
#%%----------------------------------------------------
# PERIOD ANALYSIS
#ELE_TYPE = 'elasticBeamColumn'      # MATERIAL LINEARITY
ELE_TYPE = 'nonlinearBeamColumn'     # MATERIAL AND GEOMETRIC NONLINEARITIES
MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'
ANAL_TYPE = 'PERIOD'

DATA = STEEL_COLUMN_3D(LENGTH, MAT_TYPE, ELE_TYPE, ANAL_TYPE, TOTAL_MASS)
(PERIOD_MIN_X, PERIOD_MAX_X) = DATA
print('Structure First Period:  ', PERIOD_MIN_X)
print('Structure Second Period: ', PERIOD_MAX_X) 

S01.PLOT_3D_FRAME(deformed_scale=1.0)
#%%----------------------------------------------------
# STATIC ANALYSIS
#ELE_TYPE = 'elasticBeamColumn'      # MATERIAL LINEARITY
ELE_TYPE = 'nonlinearBeamColumn'     # MATERIAL AND GEOMETRIC NONLINEARITIES
MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'
ANAL_TYPE = 'STATIC'

DATA = STEEL_COLUMN_3D(LENGTH, MAT_TYPE, ELE_TYPE, ANAL_TYPE, TOTAL_MASS)
(reactionX, reactionY, reactionZ,
 reactionMX, reactionMY, reactionMZ,
 disp_mid, 
 ele_axialforce, ele_shearforceX, ele_shearforceY,
 node_displacementsX, node_displacementsY,
 dispX, dispY) = DATA

S01.PLOT_3D_FRAME(deformed_scale=1.0)
#%%----------------------------------------------------
# PUSHOVER ANALYSIS (STATIC TIME-HISTORY ANALYSIS)
#ELE_TYPE = 'elasticBeamColumn'      # MATERIAL LINEARITY
ELE_TYPE = 'nonlinearBeamColumn'     # MATERIAL AND GEOMETRIC NONLINEARITIES
MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'
ANAL_TYPE = 'PUSHOVER'

DATA = STEEL_COLUMN_3D(LENGTH, MAT_TYPE, ELE_TYPE, ANAL_TYPE, TOTAL_MASS)
(reactionX_PUSH, reactionY_PUSH, reactionZ_PUSH,
 reactionMX_PUSH, reactionMY_PUSH, reactionMZ_PUSH,
 disp_mid_PUSH, 
 ele_axialforce_PUSH, ele_shearforce_PUSH, ele_momentforce_PUSH,
 node_displacementsX_PUSH, node_displacementsY_PUSH,
  dispX_PUSH, dispY_PUSH,
 PERIOD_MIN_PUSH, PERIOD_MAX_PUSH,
 STRESSb_PUSH, STRAINb_PUSH, STRESSt_PUSH, STRAINt_PUSH, CURVATURE_PUSH, NA_PUSH, CD_PUSH) = DATA


XDATA = dispX_PUSH
YDATA = reactionX_PUSH
XLABEL = 'Displacement in X Dir. [mm]'
YLABEL = 'Base Reaction in X Dir. [N]'
TITLE = 'Base Reaction and Displacement in X Dir. of Structure During Pushover Analysis'
COLOR = 'black'
SEMILOGY = False
PLOT(XDATA, YDATA, TITLE, XLABEL, YLABEL, COLOR, SEMILOGY)

XDATA = dispY_PUSH
YDATA = reactionY_PUSH
XLABEL = 'Displacement in Y Dir. [mm]'
YLABEL = 'Base Reaction in Y Dir. [N]'
TITLE = 'Base Reaction and Displacement in Y Dir. of Structure During Pushover Analysis'
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
plt.plot(disp_mid_PUSH, NA_PUSH, linewidth=4, color='purple', label=f"NEUTRAL AXIS: {np.abs(NA_PUSH[-1]) : .3f}")
plt.xlabel('Displacment (mm)')
plt.ylabel('Neutral Axis (mm)')
plt.title('Neutral Axis Curve During Analysis')
plt.grid(True)
plt.legend()
plt.show()

# Plot Crack Depth During Analysis
plt.figure(66, figsize=(12, 8))
plt.plot(disp_mid_PUSH, CD_PUSH, linewidth=4, color='black', label=f"TOP FIBER CRACK DEPTH: {np.abs(CD_PUSH[-1]) : .3f}")
plt.xlabel('Displacment (mm)')
plt.ylabel('Crack Depth (mm)')
plt.title('Crack Depth Curve During Analysis')
plt.grid(True)
plt.legend()
plt.show()

S01.PLOT_3D_FRAME(deformed_scale=1.0)
#%%----------------------------------------------------
# CYCLIC DISPLACEMENT ANALYSIS (STATIC TIME-HISTORY ANALYSIS)
#ELE_TYPE = 'elasticBeamColumn'      # MATERIAL LINEARITY
ELE_TYPE = 'nonlinearBeamColumn'     # MATERIAL AND GEOMETRIC NONLINEARITIES
MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'
ANAL_TYPE = 'CYCLIC_DISPLACEMENT'

DATA = STEEL_COLUMN_3D(LENGTH, MAT_TYPE, ELE_TYPE, ANAL_TYPE, TOTAL_MASS)
(reactionX_CP, reactionY_CP, reactionZ_CP,
 reactionMX_CP, reactionMY_CP, reactionMZ_CP,
 disp_mid_CP, 
 ele_axialforce_CP, ele_shearforce_CP, ele_momentforce_CP,
 node_displacementsX_CP, node_displacementsY_CP,
  dispX_CP, dispY_CP,
 PERIOD_MIN_CP, PERIOD_MAX_CP,
 STRESSb_CP, STRAINb_CP, STRESSt_CP, STRAINt_CP, CURVATURE_CP, NA_CP, CD_CP) = DATA


XDATA = dispX_CP
YDATA = reactionX_CP
XLABEL = 'Displacement in X Dir. [mm]'
YLABEL = 'Base Reaction in X Dir. [N]'
TITLE = 'Base Reaction and Displacement in X Dir. of Structure During Cyclic-Displacement Analysis'
COLOR = 'black'
SEMILOGY = False
PLOT(XDATA, YDATA, TITLE, XLABEL, YLABEL, COLOR, SEMILOGY)

XDATA = dispY_CP
YDATA = reactionY_CP
XLABEL = 'Displacement in Y Dir. [mm]'
YLABEL = 'Base Reaction in Y Dir. [N]'
TITLE = 'Base Reaction and Displacement in Y Dir. of Structure During Cyclic-Displacement Analysis'
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
plt.plot(disp_mid_CP, CD_CP, linewidth=4, color='black',label=f"TOP FIBER CRACK DEPTH: {np.abs(CD_CP[-1]) : .3f}")
plt.xlabel('Displacment (mm)')
plt.ylabel('Crack Depth (mm)')
plt.title('Crack Depth Curve During Analysis')
plt.grid(True)
plt.legend()
plt.show()

S01.PLOT_3D_FRAME(deformed_scale=1.0)
#%%----------------------------------------------------
# EXTERNAL TIME-DEPENDENT LOADING ANALYSIS (STATIC TIME-HISTORY ANALYSIS)
#ELE_TYPE = 'elasticBeamColumn'      # MATERIAL LINEARITY
ELE_TYPE = 'nonlinearBeamColumn'     # MATERIAL AND GEOMETRIC NONLINEARITIES
MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'
ANAL_TYPE = 'STATIC_EXTERNAL_TIME-DEPENDENT_LOADING'

DATA = STEEL_COLUMN_3D(LENGTH, MAT_TYPE, ELE_TYPE, ANAL_TYPE, TOTAL_MASS)

(reactionX_ETDLS, reactionY_ETDLS, reactionZ_ETDLS,
 reactionMX_ETDLS, reactionMY_ETDLS, reactionMZ_ETDLS,
 disp_mid_ETDLS, 
 ele_axialforce_ETDLS, ele_shearforce_ETDLS, ele_momentforce_ETDLS,
 node_displacementsX_ETDLS, node_displacementsY_ETDLS,
  dispX_ETDLS, dispY_ETDLS,
 PERIOD_MIN_ETDLS, PERIOD_MAX_ETDLS,
 STRESSb_ETDLS, STRAINb_ETDLS, STRESSt_ETDLS, STRAINt_ETDLS, CURVATURE_ETDLS, NA_ETDLS, CD_ETDLS) = DATA


XDATA = disp_mid_ETDLS
YDATA = reactionX_ETDLS
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

S01.PLOT_3D_FRAME(deformed_scale=1.0)
#%%----------------------------------------------------
# EXTERNAL TIME-DEPENDENT LOADING ANALYSIS (DYNAMIC TIME-HISTORY ANALYSIS)
#ELE_TYPE = 'elasticBeamColumn'      # MATERIAL LINEARITY
ELE_TYPE = 'nonlinearBeamColumn'     # MATERIAL AND GEOMETRIC NONLINEARITIES
MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'
ANAL_TYPE = 'DYNAMIC_EXTERNAL_TIME-DEPENDENT_LOADING'

DATA = STEEL_COLUMN_3D(LENGTH, MAT_TYPE, ELE_TYPE, ANAL_TYPE, TOTAL_MASS)

(time_ETDLD, reactionX_ETDLD, reactionY_ETDLD, reactionZ_ETDLD,
 reactionMX_ETDLD, reactionMY_ETDLS, reactionMZ_ETDLD,
 disp_mid_ETDLD,
ele_axialforce_ETDLD, ele_shearforce_ETDLD, ele_momentforce_ETDLD,
node_displacementsX_ETDLD, node_displacementsY_ETDLD,
dispX_ETDLD, dispY_ETDLD,
veloX_ETDLD, veloY_ETDLD,
accX_ETDLD, accY_ETDLD,
stiffness_ETDLD, PERIOD_ETDLD, damping_ratio_ETDLD,
PERIOD_MIN_ETDLD, PERIOD_MAX_ETDLD,
STRESSb_ETDLD, STRAINb_ETDLD, STRESSt_ETDLD, STRAINt_ETDLD, CURVATURE_ETDLD, NA_ETDLD, CD_ETDLD) = DATA


XDATA = disp_mid_ETDLD
YDATA = reactionX_ETDLD
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
TITLE = "elements  force  vs Time for Beam Element During Dynamic External Time-dependent Loading Analysis"
PLOT_FORCES(time_ETDLD, ele_axialforce_ETDLD, YLABEL, TITLE) # Beam Element - ELEMENTS AXIAL FORCE
# PLOT ELEMENTS SHEAR FORCE
YLABEL = 'Element Shear Force [N]'  
TITLE = "elements  force  vs Time for Beam Element During Dynamic External Time-dependent Loading Analysis"
PLOT_FORCES(time_ETDLD, ele_shearforce_ETDLD, YLABEL, TITLE) # Beam Element - ELEMENTS SHEAR FORCE
# PLOT ELEMENTS SHEAR FORCE
YLABEL = 'Element Shear Force [N]'  
TITLE = "elements  force  vs Time for Beam Element During Dynamic External Time-dependent Loading Analysis"
PLOT_FORCES(time_ETDLD, ele_momentforce_ETDLD, YLABEL, TITLE) # Beam Element - ELEMENTS SHEAR FORCE

# PLOT NODAL DISPALEMENTS
PLOT_DISPLAEMENTS(time_ETDLD, node_displacementsX_ETDLD, TITLE = "Node Displacements in X Dir. vs Time for Beam Element During Dynamic External Time-dependent Loading Analysis") # Beam Element 
PLOT_DISPLAEMENTS(time_ETDLD, node_displacementsY_ETDLD, TITLE = "Node Displacements in Y Dir. vs Time for Beam Element During Dynamic External Time-dependent Loading Analysis") # Beam Element 


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

S01.PLOT_3D_FRAME(deformed_scale=1.0)
#%%----------------------------------------------------
# FREE-VIBRATION ANALYSIS (DYNAMIC TIME-HISTORY ANALYSIS)
#ELE_TYPE = 'elasticBeamColumn'      # MATERIAL LINEARITY
ELE_TYPE = 'nonlinearBeamColumn'     # MATERIAL AND GEOMETRIC NONLINEARITIES
MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'
ANAL_TYPE = 'FREE-VIBRATION'
DATA = STEEL_COLUMN_3D(LENGTH, MAT_TYPE, ELE_TYPE, ANAL_TYPE, TOTAL_MASS)

(time_FV, reactionX_FV, reactionY_FV, reactionZ_FV,
 reactionMX_FV, reactionMY_FV, reactionMZ_FV,
 disp_mid_FV,
ele_axialforce_FV, ele_shearforce_FV, ele_momentforce_FV,
node_displacementsX_FV, node_displacementsY_FV,
dispX_FV, dispY_FV,
veloX_FV, veloY_FV,
accX_FV, accY_FV,
stiffness_FV, PERIOD_FV, damping_ratio_FV,
PERIOD_MIN_FV, PERIOD_MAX_FV,
STRESSb_FV, STRAINb_FV, STRESSt_FV, STRAINt_FV, CURVATURE_FV, NA_FV, CD_FV) = DATA

XDATA = dispX_FV
YDATA = reactionX_FV
XLABEL = 'Displacement in X Dir. [mm]'
YLABEL = 'Base Reaction in X Dir. [N]'
TITLE = 'Base Reaction and Displacement in X Dir. of Structure During Free-vibration Analysis'
COLOR = 'black'
SEMILOGY = False
PLOT(XDATA, YDATA, TITLE, XLABEL, YLABEL, COLOR, SEMILOGY)

XDATA = dispY_FV
YDATA = reactionY_FV
XLABEL = 'Displacement in Y Dir. [mm]'
YLABEL = 'Base Reaction in Y Dir. [N]'
TITLE = 'Base Reaction and Displacement in Y Dir. of Structure During Free-vibration Analysis'
COLOR = 'black'
SEMILOGY = False
PLOT(XDATA, YDATA, TITLE, XLABEL, YLABEL, COLOR, SEMILOGY)

PLOT_TIME_HISTORY(time_FV, reactionX_FV, disp_mid_FV,
                          dispX_FV, dispY_FV,
                          veloX_FV, veloY_FV,
                          accX_FV, accY_FV)

# PLOT ELEMENTS AXIAL FORCE
YLABEL = 'Element Axial Force [N]'
TITLE = "elements  force  vs Time for Beam Element During Free-vibration Analysis"
PLOT_FORCES(time_FV, ele_axialforce_FV, YLABEL, TITLE) # Beam Element - ELEMENTS AXIAL FORCE
# PLOT ELEMENTS SHEAR FORCE
YLABEL = 'Element Shear Force [N]'  
TITLE = "elements  force  vs Time for Beam Element During Free-vibration Analysis"
PLOT_FORCES(time_FV, ele_shearforce_FV, YLABEL, TITLE) # Beam Element - ELEMENTS SHEAR FORCE
# PLOT ELEMENTS SHEAR FORCE
YLABEL = 'Element Shear Force [N]'  
TITLE = "elements  force  vs Time for Beam Element During Free-vibration Analysis"
PLOT_FORCES(time_FV, ele_momentforce_FV, YLABEL, TITLE) # Beam Element - ELEMENTS SHEAR FORCE

# PLOT NODAL DISPALEMENTS
PLOT_DISPLAEMENTS(time_FV, node_displacementsX_FV, TITLE = "Node Displacements in X Dir. vs Time for Beam Element During Free-vibration Analysis") # Beam Element 
PLOT_DISPLAEMENTS(time_FV, node_displacementsY_FV, TITLE = "Node Displacements in Y Dir. vs Time for Beam Element During Free-vibration Analysis") # Beam Element 

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

S01.PLOT_3D_FRAME(deformed_scale=1.0)
#%%----------------------------------------------------
# SEISMIC ANALYSIS (DYNAMIC TIME-HISTORY ANALYSIS)
#ELE_TYPE = 'elasticBeamColumn'      # MATERIAL LINEARITY
ELE_TYPE = 'nonlinearBeamColumn'     # MATERIAL AND GEOMETRIC NONLINEARITIES
MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'
ANAL_TYPE = 'SEISMIC'

DATA = STEEL_COLUMN_3D(LENGTH, MAT_TYPE, ELE_TYPE, ANAL_TYPE, TOTAL_MASS)

(time_SEI, reactionX_SEI, reactionY_SEI, reactionZ_SEI,
 reactionMX_SEI, reactionMY_SEI, reactionMZ_SEI,
 disp_mid_SEI,
 ele_axialforce_SEI, ele_shearforce_SEI, ele_momentforce_SEI,
 node_displacementsX_SEI, node_displacementsY_SEI,
 dispX_SEI, dispY_SEI,
 veloX_SEI, veloY_SEI,
 accX_SEI, accY_SEI,
 stiffness_SEI, PERIOD_SEI, damping_ratio_SEI,
 PERIOD_MIN_SEI, PERIOD_MAX_SEI,
 STRESSb_SEI, STRAINb_SEI, STRESSt_SEI, STRAINt_SEI, CURVATURE_SEI, NA_SEI, CD_SEI) = DATA

XDATA = dispX_SEI
YDATA = reactionX_SEI
XLABEL = 'Displacement in X Dir. [mm]'
YLABEL = 'Base Reaction in X Dir. [N]'
TITLE = 'Base Reaction and Displacement in X Dir. of Structure During Seismic Analysis'
COLOR = 'black'
SEMILOGY = False
PLOT(XDATA, YDATA, TITLE, XLABEL, YLABEL, COLOR, SEMILOGY)

XDATA = dispY_SEI
YDATA = reactionY_SEI
XLABEL = 'Displacement in Y Dir. [mm]'
YLABEL = 'Base Reaction in Y Dir. [N]'
TITLE = 'Base Reaction and Displacement in Y Dir. of Structure During Seismic Analysis'
COLOR = 'black'
SEMILOGY = False
PLOT(XDATA, YDATA, TITLE, XLABEL, YLABEL, COLOR, SEMILOGY)

PLOT_TIME_HISTORY(time_SEI, reactionX_SEI, disp_mid_SEI,
                          dispX_SEI, dispY_SEI,
                          veloX_SEI, veloY_SEI,
                          accX_SEI, accY_SEI)

# PLOT ELEMENTS AXIAL FORCE
YLABEL = 'Element Axial Force [N]'
TITLE = "elements  force  vs Time for Beam Element During Seismic Analysis"
PLOT_FORCES(time_SEI, ele_axialforce_SEI, YLABEL, TITLE) # Beam Element - ELEMENTS AXIAL FORCE
# PLOT ELEMENTS SHEAR FORCE
YLABEL = 'Element Shear Force [N]'  
TITLE = "elements  force  vs Time for Beam Element During Seismic Analysis"
PLOT_FORCES(time_SEI, ele_shearforce_SEI, YLABEL, TITLE) # Beam Element - ELEMENTS SHEAR FORCE
# PLOT ELEMENTS SHEAR FORCE
YLABEL = 'Element Shear Force [N]'  
TITLE = "elements  force  vs Time for Beam Element During Seismic Analysis"
PLOT_FORCES(time_SEI, ele_momentforce_SEI, YLABEL, TITLE) # Beam Element - ELEMENTS SHEAR FORCE

# PLOT NODAL DISPALEMENTS
PLOT_DISPLAEMENTS(time_SEI, node_displacementsX_SEI, TITLE = "Node Displacements in X Dir. vs Time for Beam Element During Seismic Analysis") # Beam Element 
PLOT_DISPLAEMENTS(time_SEI, node_displacementsY_SEI, TITLE = "Node Displacements in Y Dir. vs Time for Beam Element During Seismic Analysis") # Beam Element 

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

S01.PLOT_3D_FRAME(deformed_scale=1.0)
#%%----------------------------------------------------
# --------------------------------------
#  Plot BaseShear-Displacement Analysis 
# --------------------------------------
XX = np.abs(dispX_PUSH); YY = np.abs(reactionX_PUSH); # ABSOLUTE VALUE
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
PLOT_2D(np.abs(disp_mid_PUSH), np.abs(reactionX_PUSH), X, Y, X, Y, XLABEL, YLABEL, TITLE, LEGEND01, LEGEND02, LEGEND03, COLOR='black', Z=2) 
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
S25.ENERGY_DISSIPATION_CAPACITY_INDEX(dispX_SEI, reactionX_SEI, dispX_CP, reactionX_CP)
S25.ENERGY_DISSIPATION_CAPACITY_INDEX(dispY_SEI, reactionY_SEI, dispY_CP, reactionY_CP)
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
# IN X DIR.
Ddi = np.abs(dispX_SEI)
DIyi = 100*(Ddi - X[1]) /(X[2] - X[1])
for KK in range(len(DIyi)):
    if DIyi[KK] <= 0.0: # DAMAGE INDEX -> 0.0 <= DI <=100 
        DIyi[KK] = 0.0
    if DIyi[KK] >= 100: 
        DIyi[KK] = 100.0    
        
im_valuesX = DIyi # Structural Ductility Damage Index
TITLE = 'Structural Ductility Damage Index in X Dir. (%)  [IM]'
FF.FRAGILITY_ANALYSIS(damage_states, im_valuesX, TITLE, SCATTER='True', SEMI_LOG='False')
# IN Y DIR.
Ddi = np.abs(dispY_SEI)
DIyi = 100*(Ddi - X[1]) /(X[2] - X[1])
for KK in range(len(DIyi)):
    if DIyi[KK] <= 0.0: # DAMAGE INDEX -> 0.0 <= DI <=100 
        DIyi[KK] = 0.0
    if DIyi[KK] >= 100: 
        DIyi[KK] = 100.0    
        
im_valuesY = DIyi # Structural Ductility Damage Index
TITLE = 'Structural Ductility Damage Index in Y Dir. (%)  [IM]'
FF.FRAGILITY_ANALYSIS(damage_states, im_valuesY, TITLE, SCATTER='True', SEMI_LOG='False')

#%%----------------------------------------------------
# BIVARIATE FRAGILITY MODELING
damage = {
    'Minor Damage Level': {'median1': 20.0, 'beta1': 40.0, 'median2': 20.0, 'beta2': 40.0, 'rho': 0.3},
    'Moderate Damage Level': {'median1': 50.0, 'beta1': 40.0, 'median2': 50.0, 'beta2': 40.0, 'rho': 0.3},
    'Severe Damage Level': {'median1': 70.0, 'beta1': 40.0, 'median2': 70.0, 'beta2': 40.0, 'rho': 0.3},
    'Failure Level': {'median1': 100.0, 'beta1': 40.0, 'median2': 100.0, 'beta2': 40.0, 'rho': 0.3},
}

im1_vals = im_valuesX 
im2_vals = im_valuesY

S19.FRAGILITY_ANALYSIS_TWO_PARAMETERS_FUN(damage, im1_vals, im2_vals, 'Structural Damage Index vs Duration', CONTOUR='True', SURFACE='True')
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
with pd.ExcelWriter("W(x,t)_3D_STEEL_COLUMN_8_ANA_OUTPUT.xlsx", engine='openpyxl') as writer:
    
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