###########################################################################################################
#                   >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                      #
#     COMPREHENSIVE NONLINEAR SEISMIC ASSESSMENT OF A CONCRETE T SECTION SLAB SUPERSTRUCTURE SIMPILY      #
#       SUPPORTED SHORT SPAN BRIDGE : AN OPENSEES FRAMEWORK FOR MOMENT-ROTATION ANALYSIS                  #
#---------------------------------------------------------------------------------------------------------#
#         ASSESSMENT OF DUCTILITY DAMAGE INDICES FOR STRUCTURAL ELEMENTS AND SYSTEMS AND EVALUATION       #
#                                        OF ENERGY DISSIPATION CAPACITY INDEX                             #
#---------------------------------------------------------------------------------------------------------#
#                          THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                     #
#                                   EMAIL: salar.d.ghashghaei@gmail.com                                   #
###########################################################################################################
"""
Nonlinear Seismic Performance Assessment of a simpily supported concrete slab short span bridge:
An OpenSeesPy Framework for Material and Geometric Nonlinearity Under Static, Cyclic, and Earthquake Loading

This OpenSeesPy script performs rigorous nonlinear static analysis of a simpily supported concrete slab short span bridge
 for performance-based earthquake engineering. The 2D model incorporates both material nonlinearity
 (elastic-perfectly plastic or hysteretic steel with strain hardening)
 and geometric nonlinearity to capture P-delta effects
 and large displacements—essential for collapse assessment.
 
Eight analysis protocols are implemented:    
(1) [PERIOD] : Structural Period
(2) [STATIC] : Gravity load analysis establishing dead load state
(3) [PUSHOVER] : Displacement-controlled pushover generating full capacity curves
 and plastic mechanism identification
 
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
#%%-------------------------------------------------
# BRIDGE IMAGES AND DATA
'https://www.shortspansteelbridges.org/wp-content/uploads/2020/07/Audrain-County-Concrete.png'    
# CONCRETE SOLID SLAB BRIDGE ANALYSIS WITH MIDAS SOFTWARE
'https://resource.midasuser.com/en/blog/bridge/case-study-analysis-rc-solid-slab-bridge#:~:text=Solid%20slab%20bridges%20are%20essentially,of%20slab%20are%20shown%20below.'
# BOOK: Bridge Deck Analysis
'https://www.taylorfrancis.com/books/mono/10.1201/b17475/bridge-deck-analysis-eugene-obrien-alan-connor-damien-keogh'
# BOOK: Finite Element Analysis and Design of Steel and Steel–Concrete Composite Bridges
'https://www.sciencedirect.com/book/9780124172470/finite-element-analysis-and-design-of-steel-and-steel-concrete-composite-bridges'
# BOOK: Steel-concrete composite bridge design guide
'https://www.nzta.govt.nz/assets/resources/research/reports/525/docs/525.pdf' 
# PAPER: Seismic damage prediction by deterministic methods: Concepts and procedures
'https://onlinelibrary.wiley.com/doi/10.1002/eqe.4290160507'
# PAPER: RC Slab Bridge Design Overview
'https://www.scribd.com/document/345712603/Design-of-Slab-Bridge'
# YOUTUBE: Mar 10, 2022 Bridges 07 Seismic Design of Highway Bridges
'https://www.youtube.com/watch?v=9wDZtnyHd-4'    
# YOUTUBE: Mar 2, 2022 Bridges 03 Bridge Deck Design AASHTO LRFD 2017
'https://www.youtube.com/watch?v=eEbWxb7o_lk'

# BOOK: Highway Bridge Superstructure Engineering; LRFD Approaches to Design and Analysis - Narendra Taly - CRC Press 
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
import CONCRETE_SLAB_T_SECTION_FUN_EXTRA as S14
import CONCRETE_SLAB_U_SECTION_FUN_EXTRA as S16
import CONCRETE_SLAB_E_SECTION_FUN_EXTRA as S17
#%%----------------------------------------------------
def SHORT_SPAN_T_SLAB_CONCRETE_BRIDGE(LENGTH, MAT_TYPE, ELE_TYPE, ANAL_TYPE, TOTAL_MASS):
    # Initialize model
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    
    MAX_ITERATIONS = 5000   # Maximum number of iterations
    MAX_TOLERANCE = 1.0e-10 # Specified tolerance for convergence
    
    #%% GEOMETRY DEFINITION
    NUM_POINTS = 30                  # Structural Node counts
    ε = 0.001 * LENGTH               # [mm] Bridge deck camber (Imperfection amplitude)
    # Camber is the slight upward curvature intentionally built into a bridge deck so that the surface becomes level after the applied loads.
    dy = LENGTH / (NUM_POINTS-1)     # Element length
    for i in range(NUM_POINTS):
        x = i * dy
        y = 0.0
       # y = ε * np.sin(np.pi * x / LENGTH)  # Bridge deck camber in y-direction
        ops.node(i+1, x, y)        
    #%% SUPPORTS
    ops.fix(1, 1, 1, 0)       # left  - pinned
    ops.fix(30, 1, 1, 0)      # right - pinned
    

    STEEL_TYPE = 'INELASTIC'
    
    # SLAB CONCRETE T SECTION
    Bsec = 4000.0    # [mm] width
    Hsec = 500.0     # [mm] height
    cover = 50.0     # [mm]      
    fc = 30.0        # [N/mm²] Concrete Compressive Strength
    Kfc = 1.36       # ratio of confined to unconfined concrete strength
    SEC_TAG = 1
    nFib = 100
    CONCRETE_DENSITY = 2500/1e9      # [kg/m^3] -> [kg/mm^3] Concrete Material Density
    Depth, Ele_Mass = S14.CONCRETE_SLAB_T_SECTION_FUN(SEC_TAG, STEEL_TYPE, fc, Kfc,
                               Bsec, Hsec, cover,
                               nFib, CONCRETE_DENSITY,
                               plot=True)
    
    """
    # SLAB CONCRETE E SECTION
    Bsec = 4000.0    # [mm] width
    Hsec = 500.0     # [mm] height
    cover = 50.0     # [mm]      
    fc = 30.0        # [N/mm²] Concrete Compressive Strength
    Kfc = 1.36       # ratio of confined to unconfined concrete strength
    SEC_TAG = 1
    nFib = 100
    CONCRETE_DENSITY = 2500/1e9      # [kg/m^3] -> [kg/mm^3] Concrete Material Density
    Depth, Ele_Mass = S17.CONCRETE_SLAB_E_SECTION_FUN(SEC_TAG, STEEL_TYPE, fc, Kfc,
                               Bsec, Hsec, cover,
                               nFib, CONCRETE_DENSITY,
                               plot=True)
    """
    """
    # SLAB CONCRETE U SECTION
    Bsec = 4000.0    # [mm] width
    Hsec = 500.0     # [mm] height
    cover = 50.0     # [mm]      
    fc = 30.0        # [N/mm²] Concrete Compressive Strength
    Kfc = 1.36       # ratio of confined to unconfined concrete strength
    SEC_TAG = 1
    nFib = 100
    CONCRETE_DENSITY = 2500/1e9      # [kg/m^3] -> [kg/mm^3] Concrete Material Density
    Depth, Ele_Mass = S16.CONCRETE_SLAB_U_SECTION_FUN(SEC_TAG, STEEL_TYPE, fc, Kfc,
                               Bsec, Hsec, cover,
                               nFib, CONCRETE_DENSITY,
                               plot=True)
    """
    #%% CREATE ELEMENTS
    ele_tag = 1
    transfTag = 1
    ops.geomTransf('Linear', transfTag)
    #ops.geomTransf('PDelta', transfTag)
    #ops.geomTransf('Corotational', transfTag)
    numIntgrPts = 5
    if ELE_TYPE == 'elasticBeamColumn':
        for J in range(NUM_POINTS-1):
            ops.element('elasticBeamColumn', J+1, J+1, J+2, numIntgrPts, SEC_TAG, transfTag,'-mass', Ele_Mass)

    if ELE_TYPE == 'nonlinearBeamColumn':
        for J in range(NUM_POINTS-1):
            ops.element('nonlinearBeamColumn', J+1, J+1, J+2, numIntgrPts, SEC_TAG, transfTag,'-mass', Ele_Mass) 
    ele_tag = 29+1
    
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
    ops.load(1, 0.0, 0.0, -1.0)             # [N.mm]    
    ops.load(center_node, 0.0, 0.0, -1.0)   # [N.mm]
    
    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.system('BandGeneral')
    
    if ELE_TYPE == 'ELASTIC':
        ops.algorithm('Linear')
    if ELE_TYPE == 'INELASTIC':
        ops.algorithm('Newton')
        
    time = []
    dispZ, dispY = [], []
    veloX, veloY = [], []
    accX, accY = [], []
    MOMENT, disp_mid = [], []
    stiffness = []
    OMEGA, PERIOD = [], []
    PERIOD_MIN, PERIOD_MAX = [], []
    STRESSb = []; STRAINb = []; 
    STRESSt = []; STRAINt = [];
    CURVATURE = [];
    # Initialize lists for each element's axial force, shear force and moment force
    ele_axialforce = {
        1: [],  # AXIALFORCE-01
        2: [],  # AXIALFORCE-02
        3: [],  # AXIALFORCE-03
        4: [],  # AXIALFORCE-04
        5: [],  # AXIALFORCE-05
        6: [],  # AXIALFORCE-06
        7: [],  # AXIALFORCE-07
        8: [],  # AXIALFORCE-08
        9: [],  # AXIALFORCE-09
        10: [],  # AXIALFORCE-10
        11: [],  # AXIALFORCE-11
        12: [],  # AXIALFORCE-12
        13: [],  # AXIALFORCE-13
        14: [],  # AXIALFORCE-14
        15: [],  # AXIALFORCE-15
        16: [],  # AXIALFORCE-16
        17: [],  # AXIALFORCE-17
        18: [],  # AXIALFORCE-18
        19: [],  # AXIALFORCE-19
        20: [],  # AXIALFORCE-20
        21: [],  # AXIALFORCE-21
        22: [],  # AXIALFORCE-22
        23: [],  # AXIALFORCE-23
        24: [],  # AXIALFORCE-24
        25: [],  # AXIALFORCE-25
        26: [],  # AXIALFORCE-26
        27: [],  # AXIALFORCE-27
        28: [],  # AXIALFORCE-28
        29: [],  # AXIALFORCE-29
        }  
    ele_shearforce = {
        1: [],  # SHEARFORCE-01
        2: [],  # SHEARFORCE-02
        3: [],  # SHEARFORCE-03
        4: [],  # SHEARFORCE-04
        5: [],  # SHEARFORCE-05
        6: [],  # SHEARFORCE-06
        7: [],  # SHEARFORCE-07
        8: [],  # SHEARFORCE-08
        9: [],  # SHEARFORCE-09
        10: [],  # SHEARFORCE-10
        11: [],  # SHEARFORCE-11
        12: [],  # SHEARFORCE-12
        13: [],  # SHEARFORCE-13
        14: [],  # SHEARFORCE-14
        15: [],  # SHEARFORCE-15
        16: [],  # SHEARFORCE-16
        17: [],  # SHEARFORCE-17
        18: [],  # SHEARFORCE-18
        19: [],  # SHEARFORCE-19
        20: [],  # SHEARFORCE-20
        21: [],  # SHEARFORCE-21
        22: [],  # SHEARFORCE-22
        23: [],  # SHEARFORCE-23
        24: [],  # SHEARFORCE-24
        25: [],  # SHEARFORCE-25
        26: [],  # SHEARFORCE-26
        27: [],  # SHEARFORCE-27
        28: [],  # SHEARFORCE-28
        29: [],  # SHEARFORCE-29
        }  
    ele_momentforce = {
        1: [],  # MOMENTFORCE-01
        2: [],  # MOMENTFORCE-02
        3: [],  # MOMENTFORCE-03
        4: [],  # MOMENTFORCE-04
        5: [],  # MOMENTFORCE-05
        6: [],  # MOMENTFORCE-06
        7: [],  # MOMENTFORCE-07
        8: [],  # MOMENTFORCE-08
        9: [],  # MOMENTFORCE-09
        10: [],  # MOMENTFORCE-10
        11: [],  # MOMENTFORCE-11
        12: [],  # MOMENTFORCE-12
        13: [],  # MOMENTFORCE-13
        14: [],  # MOMENTFORCE-14
        15: [],  # MOMENTFORCE-15
        16: [],  # MOMENTFORCE-16
        17: [],  # MOMENTFORCE-17
        18: [],  # MOMENTFORCE-18
        19: [],  # MOMENTFORCE-19
        20: [],  # MOMENTFORCE-20
        21: [],  # MOMENTFORCE-21
        22: [],  # MOMENTFORCE-22
        23: [],  # MOMENTFORCE-23
        24: [],  # MOMENTFORCE-24
        25: [],  # MOMENTFORCE-25
        26: [],  # MOMENTFORCE-26
        27: [],  # MOMENTFORCE-27
        28: [],  # MOMENTFORCE-28
        29: [],  # MOMENTFORCE-29
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
        8: [],  # DISP08
        9: [],  # DISP09
        10: [],  # DISP10
        11: [],  # DISP11
        12: [],  # DISP12
        13: [],  # DISP13
        14: [],  # DISP14
        15: [],  # DISP15
        16: [],  # DISP16
        17: [],  # DISP17
        18: [],  # DISP18
        19: [],  # DISP19
        20: [],  # DISP20
        21: [],  # DISP21
        22: [],  # DISP22
        23: [],  # DISP23
        24: [],  # DISP24
        25: [],  # DISP25
        26: [],  # DISP26
        27: [],  # DISP27
        28: [],  # DISP28
        29: [],  # DISP29
        30: [],  # DISP30
        }
    node_displacementsY = {
        1: [],  # DISP01
        2: [],  # DISP02
        3: [],  # DISP03
        4: [],  # DISP04
        5: [],  # DISP05
        6: [],  # DISP06
        7: [],  # DISP07
        8: [],  # DISP08
        9: [],  # DISP09
        10: [],  # DISP10
        11: [],  # DISP11
        12: [],  # DISP12
        13: [],  # DISP13
        14: [],  # DISP14
        15: [],  # DISP15
        16: [],  # DISP16
        17: [],  # DISP17
        18: [],  # DISP18
        19: [],  # DISP19
        20: [],  # DISP20
        21: [],  # DISP21
        22: [],  # DISP22
        23: [],  # DISP23
        24: [],  # DISP24
        25: [],  # DISP25
        26: [],  # DISP26
        27: [],  # DISP27
        28: [],  # DISP28
        29: [],  # DISP29
        30: [],  # DISP30
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
        MOMENT.append(np.abs(ops.eleResponse(1, 'force')[2]))           # MOMENT FORCE
        disp_mid.append(ops.nodeDisp(center_node, 2))                   # DISPLACEMENT NODE 15 
        dispZ.append(ops.nodeDisp(center_node, 3))                      # ROTATION NODE 15 IN Z DIR 
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
            
        DATA = (MOMENT, disp_mid,
                ele_axialforce, ele_shearforce, ele_momentforce,
                node_displacementsX, node_displacementsY,
                dispZ, dispY)
        
        return  DATA
    
    if ANAL_TYPE == 'PUSHOVER': # STATIC TIME-HISTORY ANALYSIS
        IDctrlDOF = 3     # 1: Horizental Dispalcement - 2: Vertical Dispalcement - 3: Rotation
        DINCR = -0.000001 # [rad] Incremental Rotation
        DMAX = -0.00040   # [rad] Max. Rotation
        ops.integrator('DisplacementControl', center_node, IDctrlDOF, DINCR)
        ops.analysis('Static')
        Nsteps =  int(np.abs(DMAX/ DINCR)) 
        STEP = 0.0
        for step in range(Nsteps):
            OK = ops.analyze(1)
            S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS)
            MOMENT.append(np.abs(ops.eleResponse(1, 'force')[2]))           # MOMENT FORCE
            disp_mid.append(ops.nodeDisp(center_node, 2))                   # DISPLACEMENT NODE 15 
            dispZ.append(ops.nodeDisp(center_node, 3))                      # ROTATION NODE 15 IN Z DIR 
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
            print(f"Step: {STEP}, Rotation: {dispZ[-1]:.8f} rad, Moment: {MOMENT[-1]:.3f} N.mm")     
        else:
            print('\n\nPUSHOVER ANALYSIS DONE.\n\n')
            
        DATA = (MOMENT, disp_mid,
                ele_axialforce, ele_shearforce, ele_momentforce,
                node_displacementsX, node_displacementsY,
                dispZ, dispY,
                np.array(PERIOD_MIN), np.array(PERIOD_MAX),
                STRESSb, STRAINb, STRESSt, STRAINt, CURVATURE)
    
        return  DATA
#%%-------------------------------------------------------
def PLOT_TIME_HISTORY(time,
                      MOMENT, disp_mid,
                      dispZ, dispY,
                      veloX, veloY,
                      accX, accY):

    import matplotlib.pyplot as plt

    data_dict = {
        "Reaction [N]": MOMENT,
        "Mid Displacement [mm]": disp_mid,
        "Disp X [mm]": dispZ,
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
LENGTH = 7000.0       # [mm] Bridge Horizental Length
TOTAL_MASS = 0.0      # [kg] Total Mass of Structure
#%%----------------------------------------------------
# PERIOD ANALYSIS
#ELE_TYPE = 'elasticBeamColumn'      # MATERIAL LINEARITY
ELE_TYPE = 'nonlinearBeamColumn'     # MATERIAL AND GEOMETRIC NONLINEARITIES
MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'
ANAL_TYPE = 'PERIOD'

DATA = SHORT_SPAN_T_SLAB_CONCRETE_BRIDGE(LENGTH, MAT_TYPE, ELE_TYPE, ANAL_TYPE, TOTAL_MASS)
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

DATA = SHORT_SPAN_T_SLAB_CONCRETE_BRIDGE(LENGTH, MAT_TYPE, ELE_TYPE, ANAL_TYPE, TOTAL_MASS)
(MOMENT, disp_mid, 
 ele_axialforce, ele_shearforce, ele_momentforce,
 node_displacementsX, node_displacementsY,
 dispZ, dispY) = DATA

S01.PLOT_2D_FRAME(deformed_scale=1.0)
#%%----------------------------------------------------
# PUSHOVER ANALYSIS (STATIC TIME-HISTORY ANALYSIS)
#ELE_TYPE = 'elasticBeamColumn'      # MATERIAL LINEARITY
ELE_TYPE = 'nonlinearBeamColumn'     # MATERIAL AND GEOMETRIC NONLINEARITIES
MAT_TYPE = 'INELASTIC'   # 'ELASTIC' OR 'INELASTIC'
ANAL_TYPE = 'PUSHOVER'

DATA = SHORT_SPAN_T_SLAB_CONCRETE_BRIDGE(LENGTH, MAT_TYPE, ELE_TYPE, ANAL_TYPE, TOTAL_MASS)
(reaction_PUSH, disp_mid_PUSH, 
 ele_axialforce_PUSH, ele_shearforce_PUSH, ele_momentforce_PUSH,
 node_displacementsX_PUSH, node_displacementsY_PUSH,
  dispZ_PUSH, dispY_PUSH,
 PERIOD_MIN_PUSH, PERIOD_MAX_PUSH,
 STRESSb_PUSH, STRAINb_PUSH, STRESSt_PUSH, STRAINt_PUSH, CURVATURE_PUSH) = DATA


XDATA = dispZ_PUSH
YDATA = reaction_PUSH
XLABEL = 'Rotation [rad]'
YLABEL = 'Moment [N.mm]'
TITLE = 'Moment and Rotation of Structure During Pushover Analysis'
COLOR = 'black'
SEMILOGY = False
PLOT(XDATA, YDATA, TITLE, XLABEL, YLABEL, COLOR, SEMILOGY)

# PLOT STRUCTURAL PERIOD DURING THE ANALYSIS
plt.figure(0, figsize=(12, 8))
plt.plot(dispZ_PUSH, PERIOD_MIN_PUSH, linewidth=3)
plt.plot(dispZ_PUSH, PERIOD_MAX_PUSH, linewidth=3)
plt.title('Period of Structure During Pushover Analysis')
plt.ylabel('Structural Period [s]')
plt.xlabel('Rotation [rad]')
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
# --------------------------------------
#  Plot BaseShear-Displacement Analysis 
# --------------------------------------
XX = np.abs(disp_mid_PUSH); YY = np.abs(reaction_PUSH); # ABSOLUTE VALUE
SLOPE_NODE = 10

DATA = S07.BILNEAR_CURVE(XX, YY, SLOPE_NODE)
X, Y, Elastic_ST, Plastic_ST, Tangent_ST, Ductility_Rito, Over_Strength_Factor = DATA

XLABEL = 'Displacement [mm]'
YLABEL = 'Base-reaction [N]'
LEGEND01 = 'Curve'
LEGEND02 = 'Bilinear Fitted'
LEGEND03 = 'Undefined'
TITLE = f'Last Data of BaseReaction-Displacement Analysis - Ductility Ratio: {X[2]/X[1]:.4f} - Over Strength Factor: {Y[2]/Y[1]:.4f}'
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
#%%----------------------------------------------------
# EXCEL OUTPUT
import pandas as pd

# Create DataFrame function
def create_df(MOMENT, disp_mid, dispZ, dispY, PERIOD_MIN, PERIOD_MAX):
    df = pd.DataFrame({
        "reaction": MOMENT,
        "disp_mid": disp_mid,
        "dispZ": dispZ,
        "dispY": dispY,
        "PERIOD_MIN": PERIOD_MIN,
        "PERIOD_MAX": PERIOD_MAX,        
    })
    return df


# Save to Excel
with pd.ExcelWriter("P(t)_BRIDGE_CONCRETE_SLAB_T_SUPERSTRUCTURE_SIMPLY_SUPPORTED_OUTPUT.xlsx", engine='openpyxl') as writer:
    
    # PUSHOVER
    df1 = create_df(reaction_PUSH, disp_mid_PUSH, dispZ_PUSH, dispY_PUSH, PERIOD_MIN_PUSH, PERIOD_MAX_PUSH)
    df1.to_excel(writer, sheet_name="PUSHOVER", index=False)
#%%----------------------------------------------------