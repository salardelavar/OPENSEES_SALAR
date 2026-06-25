######################################################################################################################
#                        >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                            #
#                              SOIL-FOUNDATION-STRUCTURE INTERACTION USING OPENSEES                                  #
#--------------------------------------------------------------------------------------------------------------------#
#                                                      EXAMPLE 02:                                                   #
#                     SOIL SIMULATED WITH SIMPLE SPRINGS AND THE  SOIL SPRINGS VALUES ARE NOT EXACT.                 #
#--------------------------------------------------------------------------------------------------------------------#
#                              THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                            #
#                                       EMAIL: salar.d.ghashghaei@gmail.com                                          #
######################################################################################################################
"""
[1] Nonlinear Frame Modeling: 2D RC frame with distributed plasticity (fiber sections) using `nonlinearBeamColumn` elements.  
[2] Material Laws:  
   - *Concrete*: `Concrete01` with confined (core) and unconfined (cover) properties.  
   - *Steel*: `Hysteretic` model with pinching, hardening, and cyclic degradation.  
[3] Seismic Loads:  
   - Pushover: Displacement-controlled lateral loading to failure.  
   - Dynamic: Uniform excitation with user-defined ground motions (X/Y components).  
[4] Damping: Rayleigh damping (a0, a1) calibrated via eigenvalue analysis (modes 1–2).  
[5] Performance Metrics:  
   - Ductility Ratio (μ): Derived from bilinearized pushover curves.  
   - Overstrength (Ω₀): Yield vs. ultimate capacity.  
   - Damage Index (DI): Normalized displacement demand/capacity.  
[6] Advanced Solver: HHT-α integrator (unconditionally stable) with Newton-Raphson iterations.  
[7] Outputs:  
   - Hysteretic responses (P-M, V-Δ, M-θ).  
   - Time-history plots (displacement, base shear).  
   - Stiffness degradation tracking.  
[8] Validation: Logarithmic decrement method for damping ratio verification.  
[9] Ductility Damage Index (DDI) Implementation:
    DDI quantifies structural damage via normalized displacement demand.
    # After bilinear fit (X[1] = Δ_y, X[2] = Δ_u):
    Dd = max(abs(DISP_Xd))  # Max dynamic displacement demand
    DI = (Dd - Dy) / (Du - Dy)  # Ductility Damage Index (X-dir)
    DI ≈ 0: Elastic response (no damage).
    DI ≥ 1: Collapse (demand ≥ ultimate capacity).
Key Innovations:  
- Combines static (pushover) and dynamic analyses for comprehensive fragility assessment.  
- Automated bilinear curve fitting for rapid performance evaluation.  
- Exportable results for FEMA P-58 compliance checks.

Utilizing OpenSees for Soil-Structure Interaction via the Spring Method:
Modeling soil-structure interaction (SSI) with lateral and rotational springs in OpenSees involves defining
the properties and behavior of the springs to simulate the interaction between the soil and structure.
Considering the large variation of stiffness across various deposits or within the same soil,
its required to conduct site or lab-tests to estimate a reliable value of stiffness (direct or coupled stiffness)
under proper drainage condition. An SSI study with an incorrect stiffness would spoil the purpose of the study.
Further, the designer need to have a good estimate of strain levels to assume a suitable value of damping. 
For soils exhibiting elastic or elastoplastic behavior, the stiffness and damping ratio can vary significantly. 
    
Here are some general estimates:
Stiffness (k in N/m)
Clay soils: Stiffness can range from 5 MPa to 20 MPa (MegaPascals, equivalent to N/m²).
Sandy soils: Stiffness can range from 10 MPa to 50 MPa.
Gravelly soils: Stiffness can range from 20 MPa to 100 MPa.
Stony soils: Stiffness can range from 200 MPa to 1500 MPa.
Damping Ratio (ζ)
Clay soils: Damping ratio typically ranges from 5% to 15%.
Sandy soils: Damping ratio typically ranges from 2% to 10%.
Gravelly soils: Damping ratio typically ranges from 1% to 5%.
Stony soils: Damping ratio typically ranges from 0.5% to 3%.
These values are approximate and can vary based on factors such as moisture content, compaction, and the presence of organic material.
 For accurate assessments, it's best to conduct site-specific geotechnical investigations.

Contact problems in finite element analysis (FEA) are quite fascinating and complex! 
They involve the interaction between two or more bodies that come into contact with each other. 

Here's a brief overview:
What are Contact Problems?
Contact problems occur when two or more bodies come into contact and interact with each other. This interaction can involve forces,
 displacements, and stresses at the contact interface. Examples include metal forming, vehicle crashes, and gear systems.

Lateral Springs in Contact Problems
Lateral springs are used in FEA to model the contact behavior between bodies. They provide a way to simulate the stiffness and damping properties
 of the contact interface. 
 
Here's how they work:
Representation: Lateral springs are represented as elements that connect the nodes of the contacting bodies.
Stiffness: The stiffness of the lateral springs determines how much resistance they provide against lateral movements.
Damping: Damping properties help in absorbing energy and reducing oscillations at the contact interface.
Challenges in Modeling Contact Problems
Contact problems are inherently nonlinear and can be challenging to solve. Some of the key challenges include:

Nonlinearity: The contact force-displacement relationship is highly nonlinear.
Convergence: Ensuring the numerical model converges can be difficult due to the complex nature of contact interactions.
Computational Cost: Tracking contact and separation can be computationally expensive.
Applications
Contact problems are crucial in various engineering fields, such as:

Automotive Engineering: For crash simulations and component interactions.
Mechanical Engineering: For gear and bearing analysis.
Manufacturing: For metal forming and stamping processes.
"""
#%% ---------------------------------------------------------------
# BOOK: 
'https://www.fema.gov/sites/default/files/documents/fema-p-2091-soil-structure-interaction.pdf'    
# BOOK: Soil-Foundation-Structure Interaction
'https://www.routledge.com/Soil-Foundation-Structure-Interaction/Orense-Chouw-Pender/p/book/9781138114128'
# REPORT: Soil-Structure Interaction for Building Structures (NEHRP Consultants)
'https://www.nehrp.gov/pdf/nistgcr12-917-21.pdf'    
# YOUTUBE: Soil Structure Interaction
'https://www.youtube.com/watch?v=klgFAtHJyFE'
# YOUTUBE: 3rd Kenji Ishihara Colloquium Series on Earthquake Engineering: Part 3 - Soil-Structure Interaction
'https://www.youtube.com/watch?v=G0z4Eq0r99Q'    
# YOUTUBE: FEMA P-2091, Webinar on A Practical Guide to Soil-Structure Interaction
'https://www.youtube.com/watch?v=Z_6rRtVgcME'
#%% ---------------------------------------------------------------

import openseespy.opensees as ops
import matplotlib.pyplot as plt
import numpy as np
import time as TI
import ANALYSIS_FUNCTION as S02
import CONCRETE_SECTION_FUN_THREE as S03
import PLOT_2D as S04
import CONCRETE_FOUNDATION_SECTION_FUN as S05
import EIGENVALUE_ANALYSIS_FUN as S07
import RAYLEIGH_DAMPING_FUN as S08
import DAMPING_RATIO_FUN as S09

# Define materials for nonlinear columns and beam
# Define parameters (units: mm, N)
# ------------------------------------------
# CONCRETE      tag   f'c        ec0   f'cu        ecu
# Cover concrete (unconfined)
fc = 25.0                # [N/mm²] Concrete Compressive Strength

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

# Foundation Section
Bf = 1200                # [mm] Depth of the Section 
Hf = 700                 # [mm] Height of the Section  
coverF = 50              # [mm] Concrete Section Cover
DIAf = 10                # [mm] # Rebar Size Diameter
AsF = np.pi*(DIAf**2)/4  # [mm²] Area of Rebar

#%% DEFINE THE ELEMENTS LENGTH
LENGTH_COL = 3000          # [mm] Column Length 
LENGTH_BM = 7000           # [mm] Beam Length 
FOUNDATION_LENGTH = 10000  # [mm] Founation Length 
FOUNDATION_NODES = 200     # Number of Nodes

#%% DEFINE PARAMEETRS FOR NONLINEAR DYNAMIC ANALYSIS
GMfact = 9810    # [mm/s²]standard acceleration of gravity or standard acceleration
SSF_X = 0.1      # Seismic Acceleration Scale Factor in X Direction
SSF_Y = 0.1      # Seismic Acceleration Scale Factor in Y Direction
iv0_X = 0.00005  # [mm/s] Initial velocity applied to the node  in X Direction
iv0_Y = 0.00005  # [mm/s] Initial velocity applied to the node  in Y Direction
st_iv0 = 0.0     # [s] Initial velocity applied starting time
SEI = 'X'        # Seismic Direction
DR = 0.05        # Intial Guess for Damping ratio
duration = 15.0  # [s] Total simulation duration
dt = 0.01        # [s] Time step
MASS = 12000     # [kg] Mass on the each column
DL_SOIL = -5.0   # [N/mm] Distribued Load from Compacted Soil On The Foundation

#%% DEFINE PARAMEETRS FOR NONLINEAR STATIC ANALYSIS 
DMAX = 175.0     # [mm] Maximum Displacement
DINCR = 0.1      # [mm] Incremental Displacement

#%% DEFINE ANALYSIS PROPERTIES
MAX_ITERATIONS = 5000    # Convergence iteration for test
MAX_TOLERANCE = 1.0e-6   # Convergence tolerance for test
#STEEL_TYPE: 1 -> WITHOUT HARDENING AND ULTIMATE STRAIN
#STEEL_TYPE: 2 -> WITH HARDENING AND ULTIMATE STRAIN
#%%------------------------------------------------------------------------------
def PD_ANALYSIS(ELE_TYPE, STEEL_TYPE, CONCRETE_TYPE, ANA_TYPE, SOIL_TYPE):
    # Reset model
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)

    # CORNER NODES
    ops.node(3, 0.0, LENGTH_COL)
    ops.node(4, LENGTH_BM, LENGTH_COL)
    
    # LEFT COLUMN MIDDLE NODES
    ops.node(5, 0.0, (1/3)*LENGTH_COL)
    ops.node(6, 0.0, (2/3)*LENGTH_COL)
    
    # RIGHT COLUMN MIDDLE NODES
    ops.node(7, LENGTH_BM, (1/3)*LENGTH_COL)
    ops.node(8, LENGTH_BM, (2/3)*LENGTH_COL)
    
    # BEAM MIDDLE NODES
    ops.node(9, (1/7)*LENGTH_BM, LENGTH_COL)
    ops.node(10, (6/7)*LENGTH_BM, LENGTH_COL)
    

    # Calculate foundation parameters
    P = (FOUNDATION_LENGTH - LENGTH_BM) * 0.5    # Overhang length on each side
    FL = (FOUNDATION_LENGTH / FOUNDATION_NODES)  # Spacing between foundation nodes

    # Create foundation nodes (at bottom)
    for II in range(1, FOUNDATION_NODES+1):
        F = -P + FL * II
        ops.node(1000+II, F, 0.0)  # Foundation nodes start at ID 5
    #%%--------------------------
    if SOIL_TYPE == 'FIXED':
        # Fix foundation nodes (prevent movement)
        for II in range(1, FOUNDATION_NODES+1):
            ops.fix(1000+II, 1, 1, 1)  # Fix all foundation nodes (x, y, rotation)
    #%%--------------------------
    if SOIL_TYPE == 'SPRING':
        # Spring constants for soil-structure interaction
        k_horizontal = 210.0e10  # [N/mm] Horizontal spring constant
        k_vertical = 150.0e10    # [N/mm] Vertical spring constant
        k_rotational = 550.0e6   # [N.mm/Rad] Rotational spring constant
        # Elastic-No Tension Material:
        ops.uniaxialMaterial('ENT', 4000, k_horizontal)     # Horizontal spring
        ops.uniaxialMaterial('ENT', 5000, k_vertical)       # Vertical spring
        ops.uniaxialMaterial('ENT', 6000, k_rotational)     # Rotational spring
        # Sprigs foundation nodes (prevent movement)
        for II in range(1, FOUNDATION_NODES+1):
            F = -P + FL * II
            ops.node(4000+II, F, 0.0)  # Foundation nodes start at ID 5
            ops.fix(4000+II, 1, 1, 1)  # Fix all foundation nodes (x, y, rotation)
            ops.element('zeroLength', 40000+II, 4000+II, 1000+II, '-mat', 4000, '-dir', 1) # Horizontal spring for soil stiffness
            ops.element('zeroLength', 50000+II, 4000+II, 1000+II, '-mat', 5000, '-dir', 2) # Vertical spring for soil stiffness
            ops.element('zeroLength', 60000+II, 4000+II, 1000+II, '-mat', 6000, '-dir', 3) # Rotational spring for soil stiffness
    #%%--------------------------
    # Create straight columns between beam and foundation
    # First find the foundation nodes that align with beam nodes
    left_col_node = int((P / FOUNDATION_LENGTH) * FOUNDATION_NODES) + 1000
    right_col_node = int(((P + LENGTH_BM) / FOUNDATION_LENGTH) * FOUNDATION_NODES) + 1000
    print(left_col_node, right_col_node)
    ops.fix(left_col_node, 1, 1, 0)   # Fix Left Column
    ops.fix(right_col_node, 1, 1, 0)  # Fix Left Column


    #%%--------------------------
    secTagC01, secTagC02 = 10, 20
    secTagB01, secTagB02 = 30, 40
    secTagF = 50000
    # COLUMN SECTION FOR ELEMENT 01, 03, 04, 06
    Kc = 1.25        # Confinement Enhancement Ratio
    #CONCRETE_TYPE = 2
    DepthC, Ele_MassC = S03.CONFINED_CONCRETE_SECTION(secTagC01, Hc, Bc, coverC, DIAc, STEEL_TYPE, CONCRETE_TYPE, fc, Kc, COL=True, DENSITY=2500/1e9)
    # COLUMN SECTION FOR ELEMENT 02, 05
    DepthC, Ele_MassC = S03.CONFINED_CONCRETE_SECTION(secTagC02, Hc, Bc, coverC, DIAc, STEEL_TYPE, CONCRETE_TYPE, fc, Kc-0.15, COL=True, DENSITY=2500/1e9)
    # BEAM SECTION FOR ELEMENT 07, 09
    DepthB, Ele_MassB = S03.CONFINED_CONCRETE_SECTION(secTagB01, Hb, Bb, coverB, DIAb, STEEL_TYPE, CONCRETE_TYPE, fc, Kc-0.1, COL=False, DENSITY=2500/1e9)
    # BEAM SECTION FOR ELEMENT 08
    DepthB, Ele_MassB = S03.CONFINED_CONCRETE_SECTION(secTagB02, Hb, Bb, coverB, DIAb, STEEL_TYPE, CONCRETE_TYPE, fc, Kc-0.2, COL=False, DENSITY=2500/1e9)
    # FOUNDATION SECTION
    DepthF, Ele_MassF = S05.CONFINED_CONCRETE_CONCRETE_FOUNDATION_SECTION_FUN(secTagF, Hf, Bf, coverF, DIAf, STEEL_TYPE, CONCRETE_TYPE, fc, Kc+0.2, DENSITY=2500/1e9)
    #%%--------------------------
    # Define geometric transformation (corotational for large displacements)
    transfTag = 1
    #ops.geomTransf('Linear', transfTag)
    #ops.geomTransf('PDelta', transfTag)
    ops.geomTransf('Corotational', transfTag)
    numIntgrPts = 5
    CONCRETE_DENSITY = 2500/1e9      # [kg/m^3] -> [kg/mm^3] Concrete Material Density
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/nonlinearBeamColumn.html
    
    if ELE_TYPE == 'INELASTIC':
        numIntgrPts = 10
        # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/nonlinearBeamColumn.html
        # LEFT COLUMN
        ops.element('nonlinearBeamColumn', 1, left_col_node, 5, numIntgrPts, secTagC01, transfTag, '-mass', Ele_MassC) # COLUMN - ELEMENT 01
        ops.element('nonlinearBeamColumn', 2, 5, 6, numIntgrPts, secTagC02, transfTag, '-mass', Ele_MassC)             # COLUMN - ELEMENT 02
        ops.element('nonlinearBeamColumn', 3, 6, 3, numIntgrPts, secTagC01, transfTag, '-mass', Ele_MassC)             # COLUMN - ELEMENT 03
        # RIGHT COLUMN
        ops.element('nonlinearBeamColumn', 4, right_col_node, 7, numIntgrPts, secTagC01, transfTag, '-mass', Ele_MassC) # COLUMN - ELEMENT 04
        ops.element('nonlinearBeamColumn', 5, 7, 8, numIntgrPts, secTagC02, transfTag, '-mass', Ele_MassC)              # COLUMN - ELEMENT 05
        ops.element('nonlinearBeamColumn', 6, 8, 4, numIntgrPts, secTagC01, transfTag, '-mass', Ele_MassC)              # COLUMN - ELEMENT 06
        # BEAM
        ops.element('nonlinearBeamColumn', 7, 3, 9, numIntgrPts, secTagB01, transfTag, '-mass', Ele_MassB)               # BEAM - ELEMENT 07
        ops.element('nonlinearBeamColumn', 8, 9, 10, numIntgrPts, secTagB02, transfTag, '-mass', Ele_MassB)              # BEAM - ELEMENT 08
        ops.element('nonlinearBeamColumn', 9, 10, 4, numIntgrPts, secTagB01, transfTag, '-mass', Ele_MassB)              # BEAM - ELEMENT 09
        
    if ELE_TYPE == 'ELASTIC':
        """
        Ec = 4700.0 * np.sqrt(fc) # [N/mm^2] Concrete Elastic Modulus
        AREAc = Bc*Hc
        IZc = Bc*(Hc**3) / 12
        AREAb = Bb*Hb
        IZb = Bb*(Hb**3) / 12
        # element elasticBeamColumn $eleTag $iNode $jNode $A $E $Iz $transfTag <-mass $massDens> <-cMass>
        # LEFT COLUMN
        ops.element('elasticBeamColumn', 1, left_col_node, 5, AREAc, Ec, IZc, transfTag, '-mass', Bc*Hc*CONCRETE_DENSITY) # COLUMN - ELEMENT 01
        ops.element('elasticBeamColumn', 2, 5, 6, AREAc, Ec, IZc, transfTag, '-mass', Bc*Hc*CONCRETE_DENSITY) # COLUMN - ELEMENT 02
        ops.element('elasticBeamColumn', 3, 6, 3, AREAc, Ec, IZc, transfTag, '-mass', Bc*Hc*CONCRETE_DENSITY) # COLUMN - ELEMENT 03
        # RIGHT COLUMN
        ops.element('elasticBeamColumn', 4, right_col_node, 7, AREAc, Ec, IZc, transfTag, '-mass', Bc*Hc*CONCRETE_DENSITY) # COLUMN - ELEMENT 04
        ops.element('elasticBeamColumn', 5, 7, 8, AREAc, Ec, IZc, transfTag, '-mass', Bc*Hc*CONCRETE_DENSITY) # COLUMN - ELEMENT 05
        ops.element('elasticBeamColumn', 6, 8, 4, AREAc, Ec, IZc, transfTag, '-mass', Bc*Hc*CONCRETE_DENSITY) # COLUMN - ELEMENT 06
        # BEAM
        ops.element('elasticBeamColumn', 7, 3, 9, AREAb, Ec, IZb, transfTag, '-mass', Bb*Hb*CONCRETE_DENSITY)  # BEAM - ELEMENT 07
        ops.element('elasticBeamColumn', 8, 9, 10, AREAb, Ec, IZb, transfTag, '-mass', Bb*Hb*CONCRETE_DENSITY) # BEAM - ELEMENT 08
        ops.element('elasticBeamColumn', 9, 10, 4, AREAb, Ec, IZb, transfTag, '-mass', Bb*Hb*CONCRETE_DENSITY) # BEAM - ELEMENT 09
        """
        # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/elasticBeamColumn.html
        # LEFT COLUMN
        ops.element('elasticBeamColumn', 1, left_col_node, 5, secTagC01, transfTag, '-mass', Ele_MassC) # COLUMN - ELEMENT 01
        ops.element('elasticBeamColumn', 2, 5, 6, secTagC02, transfTag, '-mass', Ele_MassC)             # COLUMN - ELEMENT 02
        ops.element('elasticBeamColumn', 3, 6, 3, secTagC01, transfTag, '-mass', Ele_MassC)             # COLUMN - ELEMENT 03
        # RIGHT COLUMN
        ops.element('elasticBeamColumn', 4, right_col_node, 7, secTagC01, transfTag, '-mass', Ele_MassC) # COLUMN - ELEMENT 04
        ops.element('elasticBeamColumn', 5, 7, 8, secTagC02, transfTag, '-mass', Ele_MassC)              # COLUMN - ELEMENT 05
        ops.element('elasticBeamColumn', 6, 8, 4, secTagC01, transfTag, '-mass', Ele_MassC)              # COLUMN - ELEMENT 06
        # BEAM
        ops.element('elasticBeamColumn', 7, 3, 9, secTagB01, transfTag, '-mass', Ele_MassB)              # BEAM - ELEMENT 07
        ops.element('elasticBeamColumn', 8, 9, 10, secTagB02, transfTag, '-mass', Ele_MassB)             # BEAM - ELEMENT 08
        ops.element('elasticBeamColumn', 9, 10, 4, secTagB01, transfTag, '-mass', Ele_MassB)             # BEAM - ELEMENT 09
    #%%--------------------------
    # FOUNDATION
    transfTagF = 2
    #ops.geomTransf('Linear', transfTagF)
    #ops.geomTransf('PDelta', transfTagF)
    ops.geomTransf('Corotational', transfTagF)
    
    ops.timeSeries('Linear', 0)
    ops.pattern('Plain', 0, 0)
    for II in range(1, 101):
        if ELE_TYPE == 'INELASTIC':
            ops.element('nonlinearBeamColumn', 999+II, 1000+II, 1000+II+1, numIntgrPts, secTagF, transfTagF, '-mass', Ele_MassF) # FOUNDATION
        if ELE_TYPE == 'ELASTIC':
            ops.element('elasticBeamColumn', 999+II, 1000+II, 1000+II+1, secTagF, transfTagF, '-mass', Ele_MassF) # FOUNDATION
        # Apply distributed load to all beams in the structure
        # eleLoad -ele $eleTag1 <$eleTag2 ....> -type -beamUniform $Wz <$Wx>
        ops.eleLoad('-ele', 999+II, '-type', '-beamUniform', DL_SOIL) 
    #%%--------------------------
    if ANA_TYPE == 'PUSHOVER':
        WEIGHT = MASS * GMfact
        # Data storage
        reactionX, reactionY, reactionZ = [], [], []
        DISP_X, DISP_Y, ROT = [], [], []
        KA, KS, KI, STEP = [], [], [], []
        PERIOD_MIN, PERIOD_MAX = [], []
        
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
            # NOTICE: YOU USED FONDATION ELEMENT, SO WE MUST OUTPUT COLUMN BOTTOM NODE ELEMNT FORCE
            ops.reactions()
            reactionX.append(ops.nodeReaction(left_col_node, 1) + ops.nodeReaction(right_col_node, 1))       # SHEAR BASE REACTION
            reactionY.append(ops.nodeReaction(left_col_node, 2) + ops.nodeReaction(right_col_node, 2))       # AXIAL BASE REACTION
            reactionZ.append(ops.nodeReaction(left_col_node, 3) + ops.nodeReaction(right_col_node, 3))       # MOMENT BASE REACTION
            disp_X = ops.nodeDisp(3, 1) # LATERAL DISPLACEMENT IN X FOR NODE 3
            disp_Y = ops.nodeDisp(3, 2) # LATERAL DISPLACEMENT IN Y FOR NODE 3
            rot = ops.nodeDisp(3, 3)    # ROTATION IN Z FOR NODE 3
            DISP_X.append(disp_X)
            DISP_Y.append(disp_Y)
            ROT.append(rot)
            KS.append(np.abs(reactionX)/np.abs(disp_X)) # LATERAL STIFFNESS IN X
            KA.append(np.abs(reactionY)/np.abs(disp_Y)) # LATERAL STIFFNESS IN Y
            KI.append(np.abs(reactionZ)/np.abs(rot))    # ROTATIONAL STIFFNESS IN Z
            STEP.append(step)
            # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
            #PERIODmin, PERIODmax = S08.RAYLEIGH_DAMPING(3, 0.5*DR, DR, 0, 1)
            PERIODmin, PERIODmax = S07.EIGENVALUE_ANALYSIS(3, PLOT=True)
            PERIOD_MIN.append(PERIODmin)
            PERIOD_MAX.append(PERIODmax)
            print(f"Step: {STEP[-1]}, Displacement: {DISP_X[-1]:.4f} mm, Reaction: {reactionX[-1]:.2f} N")     
            
        # Compute modal properties
        ops.modalProperties("-print", "-file", f"SALAR_ModalReport_{ANA_TYPE}.txt", "-unorm") 
        
        DATA  = (reactionX, reactionY, reactionZ,
                 DISP_X, DISP_Y, ROT,
                 KA, KS, KI,
                 STEP, PERIOD_MIN, PERIOD_MAX)
        return DATA     
    
    if ANA_TYPE == 'DYNAMIC':
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
        ops.integrator('Newmark', 0.5, 0.25) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/newmark.html
        #alpha = 1;gamma=1.5-alpha; beta=((2-alpha)**2)/4;
        #ops.integrator('HHT', alpha, gamma, beta) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/hht.html
        ops.algorithm('Newton') # INFO LINK: https://openseespydoc.readthedocs.io/en/stable/src/algorithm.html
        ops.analysis('Transient')
        
        # Calculate Rayleigh damping factors
        #Lambda01 = ops.eigen('-fullGenLapack', 2)  # eigenvalue mode 2
        Lambda01 = ops.eigen('-genBandArpack', 2) # eigenvalue mode 2
        Omega01 = np.power(max(Lambda01), 0.5)
        Omega02 = np.power(min(Lambda01), 0.5)
        a0 = (2 * Omega01 * Omega02 * DR) / (Omega01 + Omega02) # c = a0 * m : Mass-proportional damping
        a1 = (DR * 2) / (Omega01 + Omega02)   # c = a1 * k : Stiffness-proportional damping
        # Apply Rayleigh damping
        ops.rayleigh(a0, a1, 0, 0)   # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/reyleigh.html
        #ops.rayleigh(0, 0, 2 * DR * Omega01, 0) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/reyleigh.html
        PERIOD_01 = (np.pi*2) / Omega01 # Structure First Period
        PERIOD_02 = (np.pi*2) / Omega02 # Structure Second Period
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
        print('Seismic Defined Done.')
            
        # Data storage
        reactionX, reactionY, reactionZ = [], [], []
        DISP_X, DISP_Y, ROT = [], [], []
        KA, KS, KI, STEP = [], [], [], []
        time = []
        displacement = []
        velocity_X, velocity_Y = [], []
        acceleration_X, acceleration_Y = [], []
        PERIOD_MIN, PERIOD_MAX = [], []
        
        stable = 0
        current_time = 0.0
        
        while stable == 0 and current_time < duration:
            stable = ops.analyze(1, dt)
            #S02.ANALYSIS(stable, 1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS
            current_time = ops.getTime()
            time.append(current_time)
            # NOTICE: YOU USED FONDATION ELEMENT, SO WE MUST OUTPUT COLUMN BOTTOM NODE ELEMNT FORCE
            ops.reactions()
            reactionX.append(ops.nodeReaction(left_col_node, 1) + ops.nodeReaction(right_col_node, 1))       # SHEAR BASE REACTION
            reactionY.append(ops.nodeReaction(left_col_node, 2) + ops.nodeReaction(right_col_node, 2))       # AXIAL BASE REACTION
            reactionZ.append(ops.nodeReaction(left_col_node, 3) + ops.nodeReaction(right_col_node, 3))       # MOMENT BASE REACTION
            disp_X = ops.nodeDisp(3, 1) # LATERAL DISPLACEMENT IN X FOR NODE 3
            disp_Y = ops.nodeDisp(3, 2) # LATERAL DISPLACEMENT IN Y FOR NODE 3
            rot = ops.nodeDisp(3, 3)    # ROTATION IN Z FOR NODE 3
            velocity_X.append(ops.nodeVel(3, 1))       # LATERAL VELOCITY IN X FOR NODE 3
            acceleration_X.append(ops.nodeAccel(3, 1)) # LATERAL ACCELERATION IN X FOR NODE 3
            velocity_Y.append(ops.nodeVel(3, 2))       # LATERAL VELOCITY IN Y FOR NODE 3
            acceleration_Y.append(ops.nodeAccel(3, 2)) # LATERAL ACCELERATION IN Y FOR NODE 3
            DISP_X.append(disp_X)
            DISP_Y.append(disp_Y)
            ROT.append(rot)
            KS.append(np.abs(reactionX)/np.abs(disp_X)) # LATERAL STIFFNESS IN X
            KA.append(np.abs(reactionY)/np.abs(disp_Y)) # LATERAL STIFFNESS IN Y
            KI.append(np.abs(reactionZ)/np.abs(rot))    # ROTATIONAL STIFFNESS IN Z
            print(f"Time: {time[-1]:.4f}, Displacement: {DISP_X[-1]:.4f} mm, Reaction: {reactionX[-1]:.2f} N") 
        
        damping_ratio = S09.DAMPING_RATIO(DISP_X)  
        
        # Compute modal properties
        ops.modalProperties("-print", "-file", f"SALAR_ModalReport_{ANA_TYPE}.txt", "-unorm") 
        
        DATA = (reactionX, reactionY, reactionZ,
                DISP_X, DISP_Y, ROT,
                KA, KS, KI,
                time, velocity_X, velocity_Y,
                acceleration_X, acceleration_Y,
                PERIOD_01, PERIOD_02, damping_ratio,
                PERIOD_MIN, PERIOD_MAX)
        
        return DATA
#%%------------------------------------------------------------------------------
# Analysis Durations for Nonlinear Static Analysis:
starttime = TI.process_time()

# RUN NONLINEAR STATIC ANALYSIS
DATA = PD_ANALYSIS(ELE_TYPE='INELASTIC', STEEL_TYPE=2, CONCRETE_TYPE=2, ANA_TYPE='PUSHOVER', SOIL_TYPE = 'SPRING')
(FORCE_Sp, FORCE_Ap, MOMENTp,
 DISP_Xp, DISP_Yp, ROTp,
 KAp, KSp, KIp, STEPp,
 PERIOD_MINp, PERIOD_MAXp) = DATA

totaltime = TI.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n') 

# %% Plot 2D Frame Shapes for Nonlinear Static Analysis
S04.PLOT_2D_FRAME(deformed_scale=1)  # Adjust scale factor as needed  
#%%------------------------------------------------------------------------------
# Analysis Durations for Nonlinear Dynamic Analysis:
starttime = TI.process_time()

# RUN NONLINEAR DYNAMIC ANALYSIS
DATA = PD_ANALYSIS(ELE_TYPE='INELASTIC', STEEL_TYPE=2, CONCRETE_TYPE=2, ANA_TYPE='DYNAMIC', SOIL_TYPE = 'SPRING')
(FORCE_Sd, FORCE_Ad, MOMENTd,
 DISP_Xd, DISP_Yd, ROTd,
 KAd, KSd, KId,
 timed, velocity_Xd, velocity_Yd,
 acceleration_Xd, acceleration_Yd,
 PERIOD_01d, PERIOD_02d, damping_ratio,
 PERIOD_MINd, PERIOD_MAXd) = DATA

print(f"\n Period 01: {PERIOD_01d:.4e}  - Period 02: {PERIOD_02d:.4e}")

totaltime = TI.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n') 

# %% Plot 2D Frame Shapes for Nonlinear Dynamic Analysis
S04.PLOT_2D_FRAME(deformed_scale=1000)  # Adjust scale factor as needed  
#%%------------------------------------------------------------------------------
plt.figure(1, figsize=(12, 8))
plt.plot(MOMENTp, FORCE_Ap, color='black')
plt.plot(MOMENTd, FORCE_Ad, color='purple')
#plt.scatter(MOMENTp, FORCE_Ap, color='black', linewidth=2)
#plt.scatter(MOMENTd, FORCE_Ad, color='purple', linewidth=2)
plt.title('P-M Interaction')
plt.ylabel('Axial Force [N]')
plt.xlabel('Bending Moment [N.mm]')
plt.legend(['PUSHOVER', 'DYNAMIC'])
plt.grid()
plt.show()

plt.figure(2, figsize=(12, 8))
plt.plot(DISP_Xp, FORCE_Sp, color='green', linewidth=2)
plt.plot(DISP_Xd, FORCE_Sd, color='lime', linewidth=2)
#plt.scatter(DISP_Xp, FORCE_Sp, color='green', linewidth=2)
#plt.scatter(DISP_Xd, FORCE_Sd, color='lime', linewidth=2)
plt.title('SHEAR FORCE-DISPLACEMENT DIAGRAM')
plt.ylabel('Shear Force [N]')
plt.xlabel('Displacement  in X [mm]')
plt.legend(['PUSHOVER', 'DYNAMIC'])
plt.grid()
plt.show()

plt.figure(3, figsize=(12, 8))
plt.plot(DISP_Yp, FORCE_Ap, color='purple', linewidth=2)
plt.plot(DISP_Yd, FORCE_Ad, color='green', linewidth=2)
#plt.scatter(DISP_Yp, FORCE_Ap, color='purple', linewidth=2)
#plt.scatter(DISP_Yd, FORCE_Ad, color='green', linewidth=2)
plt.title('AXIAL FORCE-DISPLACEMENT DIAGRAM')
plt.ylabel('Axial Force [N]')
plt.xlabel('Displacement in Y [mm]')
plt.legend(['PUSHOVER', 'DYNAMIC'])
plt.grid()
plt.show()

plt.figure(4, figsize=(12, 8))
plt.plot(ROTp, MOMENTp, color='red', linewidth=2)
plt.plot(ROTd, MOMENTd, color='pink', linewidth=2)
#plt.scatter(ROTp, MOMENTp, color='red', linewidth=2)
#plt.scatter(ROTd, MOMENTd, color='pink', linewidth=2)
plt.title('MOMENT-ROTATION DIAGRAM')
plt.ylabel('Moment [kN.mm]')
plt.xlabel('Rotation [rad]')
plt.legend(['PUSHOVER', 'DYNAMIC'])

plt.grid()
plt.show()
"""
plt.figure(5, figsize=(12, 8))
#plt.plot(KIp, KSp, color='black', linewidth=2)
#plt.plot(KId, KSd, color='grey', linewidth=2)
plt.scatter(KIp, KSp, color='black', linewidth=2)
plt.scatter(KId, KSd, color='grey', linewidth=2)
plt.title('ROTATIONAL STIFFNESS-LATERAL STIFFNESS DIAGRAM')
plt.ylabel('Rotational Stiffness [N.mm/Rad]')
plt.xlabel('Lateral Stiffness in X Dir. [N/mm]')
plt.semilogx()
plt.semilogy()
plt.legend(['PUSHOVER', 'DYNAMIC'])
plt.grid()
plt.show()

plt.figure(6, figsize=(12, 8))
#plt.plot(KIp, KAp, color='black', linewidth=2)
#plt.plot(KId, KAd, color='grey', linewidth=2)
plt.scatter(KIp, KAp, color='black', linewidth=2)
plt.scatter(KId, KAd, color='grey', linewidth=2)
plt.title('ROTATIONAL STIFFNESS-LATERAL STIFFNESS DIAGRAM')
plt.ylabel('Rotational Stiffness [N.mm/Rad]')
plt.xlabel('Lateral Stiffness in Y Dir. [N/mm]')
plt.semilogx()
plt.semilogy()
plt.legend(['PUSHOVER', 'DYNAMIC'])
plt.grid()
plt.show()
"""
plt.figure(7, figsize=(12, 8))
plt.plot(timed, FORCE_Ad, color='brown', linewidth=2)
#plt.scatter(timed, FORCE_Ad, color='brown', linewidth=2)
plt.title('Axial Force During the Analysis')
plt.ylabel('Axial Force [N]')
plt.xlabel('Times')
plt.grid()
plt.show()

plt.figure(8, figsize=(12, 8))
plt.plot(timed, FORCE_Sd, color='purple', linewidth=2)
#plt.scatter(timed, FORCE_Sd, color='purple', linewidth=2)
plt.title('Shear Force During the Analysis')
plt.ylabel('Shear Force [N]')
plt.xlabel('Times')
plt.grid()
plt.show()

plt.figure(9, figsize=(12, 8))
plt.plot(timed, MOMENTd, color='green', linewidth=2)
#plt.scatter(timed, MOMENTd, color='green', linewidth=2)
plt.title('Moment During the Analysis')
plt.ylabel('Moment [N.mm]')
plt.xlabel('Times')
plt.grid()
plt.show()

plt.figure(10, figsize=(12, 8))
plt.plot(timed, DISP_Xd, color='brown', linewidth=2)
#plt.scatter(timed, DISP_Xd, color='brown', linewidth=2)
plt.title('Displacement During the Analysis')
plt.ylabel('Displacement - X [mm]')
plt.xlabel('Times')
plt.grid()
plt.show()

plt.figure(11, figsize=(12, 8))
plt.plot(timed, DISP_Yd, color='blue', linewidth=2)
#plt.scatter(timed, DISP_Xd, color='brown', linewidth=2)
plt.title('Displacement During the Analysis')
plt.ylabel('Displacement - Y [mm]')
plt.xlabel('Times')
plt.grid()
plt.show()

plt.figure(12, figsize=(12, 8))
plt.plot(timed, ROTd, color='black', linewidth=2)
#plt.scatter(timed, ROTd, color='black', linewidth=2)
plt.title('Rotation During the Analysis')
plt.ylabel('Rotation [rad]')
plt.xlabel('Times')
plt.grid()
plt.show()
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

DISP_ZX = MAX_ABS(DISP_Xd)  
DISP_ZY = MAX_ABS(DISP_Yd) 
VELO_Z = MAX_ABS(velocity_Xd) 
ACCE_Z = MAX_ABS(acceleration_Xd) 
BASE_Z = MAX_ABS(FORCE_Sd) 

plt.figure(1, figsize=(8, 6))
plt.plot(timed, DISP_Xd, color='blue', linewidth=2)
plt.plot(timed, DISP_ZX, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Displacement in X [mm]')
plt.title(f'Time vs Displacement - MAX. ABS: {DISP_ZX[-1]} | ξ (Calculated): {damping_ratio:.5e} %')
plt.grid()
plt.show()

plt.figure(2, figsize=(8, 6))
plt.plot(timed, DISP_Yd, color='blue', linewidth=2)
plt.plot(timed, DISP_ZY, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Displacement in Y [mm]')
plt.title(f'Time vs Displacement - MAX. ABS: {DISP_ZY[-1]}')
plt.grid()
plt.show()

plt.figure(3, figsize=(8, 6))
plt.plot(timed, velocity_Xd, color='blue', linewidth=2)
plt.plot(timed, VELO_Z, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Velocity in X [mm/s]')
plt.title(f'Time vs Velocity - MAX. ABS: {VELO_Z[-1]}')
plt.grid()
plt.show()

plt.figure(4, figsize=(8, 6))
plt.plot(timed, acceleration_Xd, color='blue', linewidth=2)
plt.plot(timed, ACCE_Z, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Acceleration in X [mm/s^2]')
plt.title(f'Time vs Acceleration - MAX. ABS: {ACCE_Z[-1]}')
plt.grid()
plt.show()

plt.figure(5, figsize=(8, 6))
plt.plot(timed, FORCE_Sd, color='blue', linewidth=2)
plt.plot(timed, BASE_Z, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Base-reaction [N]')
plt.title(f'Time vs Base-reaction - MAX. ABS: {BASE_Z[-1]}')
plt.grid()
plt.show()

plt.figure(6, figsize=(8, 6))
plt.plot(DISP_Xd, FORCE_Sd, color='black', linewidth=2)
plt.xlabel('Displacement in X [mm]')
plt.ylabel('Shear Base-reaction [N]')
plt.title(f'Displacement vs Shear Base-reaction')
plt.grid()
plt.show()

plt.figure(7, figsize=(8, 6))
plt.plot(DISP_Yd, FORCE_Ad, color='brown', linewidth=2)
plt.xlabel('Displacement in Y [mm]')
plt.ylabel('Axial Base-reaction [N]')
plt.title(f'Displacement vs Axial Base-reaction')
plt.grid()
plt.show()

#%%------------------------------------------------------------------------------
import BILINEAR_CURVE as BC

# --------------------------------------
#  Plot BaseShear-Displacement Analysis 
# --------------------------------------
XX = np.abs(DISP_Xp); YY = np.abs(FORCE_Sp); # ABSOLUTE VALUE
SLOPE_NODE = 10

DATA = BC.BILNEAR_CURVE(XX, YY, SLOPE_NODE)
X, Y, Elastic_ST, Plastic_ST, Tangent_ST, Ductility_Rito, Over_Strength_Factor = DATA

XLABEL = 'Displacement in X [mm]'
YLABEL = 'Base-Shear Reaction [N]'
LEGEND01 = 'Curve'
LEGEND02 = 'Bilinear Fitted'
LEGEND03 = 'Undefined'
TITLE = f'Last Data of BaseShear-Displacement Analysis - Ductility Ratio: {X[2]/X[1]:.4f} - Over Strength Factor: {Y[2]/Y[1]:.4f}'
COLOR = 'black'
BC.PLOT_2D(np.abs(DISP_Xp), np.abs(FORCE_Sp), X, Y, X, Y, XLABEL, YLABEL, TITLE, LEGEND01, LEGEND02, LEGEND03, COLOR='black', Z=2) 
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
Dd = np.max(np.abs(DISP_Xd))
DIx = (100*Dd - X[1]) /(X[2] - X[1])
print(f'Structural Ductility Damage Index in X Direction: {DIx:.4f} [%]')

# ---------------------------------------
#  Plot BaseAxial-Displacement Analysis
# ---------------------------------------
XX = np.abs(DISP_Yp); YY = np.abs(FORCE_Ap); # ABSOLUTE VALUE
SLOPE_NODE = 10

DATA = BC.BILNEAR_CURVE(XX, YY, SLOPE_NODE)
X, Y, Elastic_ST, Plastic_ST, Tangent_ST, Ductility_Rito, Over_Strength_Factor = DATA

XLABEL = 'Displacement in Y [mm]'
YLABEL = 'Base-Axial Reaction [N]'
LEGEND01 = 'Curve'
LEGEND02 = 'Bilinear Fitted'
LEGEND03 = 'Undefined'
TITLE = f'Last Data of BaseAxial-Displacement Analysis - Ductility Ratio: {X[2]/X[1]:.4f} - Over Strength Factor: {Y[2]/Y[1]:.4f}'
COLOR = 'black'
BC.PLOT_2D(np.abs(DISP_Yp), np.abs(FORCE_Ap), X, Y, X, Y, XLABEL, YLABEL, TITLE, LEGEND01, LEGEND02, LEGEND03, COLOR='black', Z=2) 
#print(f'\t\t Ductility Ratio: {YY[2]/YY[1]:.4f}')

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
Dd = np.max(np.abs(DISP_Yd))
DIy = (100*Dd - X[1]) /(X[2] - X[1])
print(f'Structural Ductility Damage Index in Y Direction: {DIy:.4f} [%]')
#%%------------------------------------------------------------------------------
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
Ddi = np.abs(DISP_Xd)
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
# EXPORT DATA TO EXCEL
import pandas as pd
DATA_TOTAL = { 
    'DISP_X_WO': DISP_Xd,
    'DISP_Y_WO': DISP_Yd,
    'ROTATION_WO': ROTd,
    'VELO_WO': velocity_Xd,
    'ACCEL_WO': acceleration_Xd,
    'AXIAL_FORCE_WO': FORCE_Ad,
    'SHEAR_FORCE_WO': FORCE_Sd,
    'MOMENT_WO': MOMENTd,
    'AXIAL_RIGIDITY_WO': np.abs(FORCE_Ad),
    'ROTATIONAL_ST_WO': KId,
    'LATERAL_ST_Y_WO': KAd,
    'LATERAL_ST_X_WO': KSd,
}
# Convert to DataFrame
results_df = pd.DataFrame(DATA_TOTAL)
# Export the DataFrame to an Excel file
results_df.to_excel('CONCRETE_FRAME_SOIL-FOUNDATION-STRUCTURE_RESULTS.xlsx', index=False) 
#%%------------------------------------------------------------------------------
# Print out the state of nodes 3 and 4
ops.printModel("node",3, 4)
# Print out the state of element 1 , 2 and 3
ops.printModel("ele", 1, 2 , 3)
# Print the Model
#printModel()
ops.printModel("-JSON", "-file", "CONCRETE_FRAME_SOIL-FOUNDATION-STRUCTURE.json")
#%%------------------------------------------------------------------------------   

    
