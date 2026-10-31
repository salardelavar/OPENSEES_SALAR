###########################################################################################################
#                    >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                     #
#  NONLINEAR PUSHOVER ANALYSIS OF A ONE-BAY, ONE-STORY TIMBER FRAME WITH SEMI-RIGID AND RIGID CONNECTIONS #
#                                              IN OPENSEES                                                #
#---------------------------------------------------------------------------------------------------------#
# Key Features:                                                                                           #
# - Nonlinear finite-element modeling of a one-bay, one-story timber frame using OpenSeesPy.              #
# - Simulation of hysteretic beam-to-column and column-to-foundation connections using Pinching4.         #
# - Consideration of semi-rigid and rigid beam-to-column connection configurations.                       #
# - Modeling of elastic timber beams and columns with nonlinear connection behavior.                      #
# - Displacement-controlled nonlinear static pushover analysis.                                           #
# - Evaluation of global lateral response and base shear capacity.                                        #
# - Recording of roof displacement and base reactions during analysis.                                    #
# - Visualization of the deformed frame configuration and pushover capacity curve.                        #
# - Export of numerical results and graphical outputs for further analysis.                               #
#---------------------------------------------------------------------------------------------------------#
#        THIS PROGRAM WRITTEN BY DR. KASHIF SALMAN AND MODIFIED BY SALAR DELAVAR GHASHGHAEI (QASHQAI)     #
#                                       EMAIL: kashifsalman96@gmail.com                                   #   
#                                            salar.d.ghashghaei@gmail.com                                 #
########################################################################################################### 

"""
This Python script develops a nonlinear finite-element model of a one-bay, one-story (1B1S) timber frame with
 hysteretic beam-to-column and column-to-foundation connections using OpenSeesPy. The frame has a single 4,100
 mm bay and a story height of 2,740 mm, with elastic beam and column members connected through zero-length elements
 representing nonlinear rotational behavior. The connection response is defined using the `Pinching4` uniaxial material
 model, which captures pinching effects, stiffness and strength degradation, and energy-based deterioration under cyclic
 loading. Alternatively, the beam-to-column connections can be assigned rigid elastic behavior by changing `BCJType` to 
 `"Rigid"`. The model uses elastic material properties for the timber structural members and highly stiff elastic materials
 to represent shear constraints at the connections. After defining the nodal coordinates, elements, boundary conditions, and
 connection properties, the script establishes nodal recorders to capture roof displacement and base reactions.
 A displacement-controlled nonlinear static pushover analysis is then performed by incrementally increasing the horizontal
 displacement of the roof control node up to 260 mm, with a displacement increment of 1 mm per analysis step.
 The resulting roof displacement and accumulated base shear are extracted from the recorder files and used to generate the
 global pushover capacity curve. Finally, the `PLOT_2D_FRAME` function is called to visualize the undeformed and deformed
 frame configurations. Overall, this model provides a computational framework for investigating the nonlinear lateral response,
 hysteretic connection behavior, and global seismic capacity of a one-bay, one-story timber frame with semi-rigid or rigid connection assumptions.
"""
# PAPER: Seismic performance of glulam timber post and beam structures with and without light frame timber shear wall infill 
'https://www.sciencedirect.com/science/article/abs/pii/S2352710222009767'    
# BOOK: Advanced Analysis And Design Of Steel Frames-Guo-Qiang Li-WILEY-2007
# BOOK: Analysis and Design of Steel and Composite Structures-Qing Quan Liang-CRC Press-2015
# BOOK: Handbook of strcutural steel Connection Design and Detail_Akbar R Tamboli
# BOOK: Non-Linear Static and Cyclic Analysis of Steel Frames with Semi-Rigid Connections-S.L. Chan and P.P.T. Chui-Elsevier-2000
# BOOK: Plastic Design and Second-Order Analysis of Steel Frames-W. F. Chen, I. Sohal-Springer New York-1995
# BOOK: Seismic Design of Steel Structures-Victor Gioncu & Federico Mazzolani-CRC Press-2014    
# BOOK: Dynamics of Structure and Foundation; A Unified Approach-1. Fundamentals-Indrajit Chowdhury-CRC Press-2009
# BOOK: Dynamics of structure and foundation; A Uniﬁed Approach-2. Applications-Indrajit Chowdhury-CRC Press-2009
# BOOK: Plastic Design and Second-Order Analysis of Steel Frames-W. F. Chen, I. Sohal-Springer New York-1995
# BOOK: Structural Stability-W. F. Chen'
# BOOK: Differential Equations for Engineers-Wei-Chau Xie-CAMBRIDGE-2010
# BOOK: Structural dynamics THEORY AND COMPUTATION MARIO BAZ-5 EDITION   
# BOOK: Dynamics of Structures in SI Units -  Anil Kumar Chopra
'https://share.google/TDN5O4eWmw5pH8zUH'
# EXTRACT DATA:
'https://web.eecs.utk.edu/~dcostine/personal/PowerDeviceLib/DigiTest/index.html'
#%% ------------------------------------------
import openseespy.opensees as ops
import os
import numpy as np
import matplotlib.pyplot as plt
import PLOT_2D as S01
import ANALYSIS_FUNCTION as S02
#%% ------------------------------------------
# General design information

BCJType = "Hysteresis"  # "Hysteresis" or "Rigid"
Nstory = 1
BayTotal = 1
Story1Height = 2740.0  # [mm] STORY 1 HEIGHT
StoryHeight = 2740.0   # [mm] STORY HEIGHT
FrameLength = [4100.0] # [mm] BAY LENGTH


# Define Analysis Properties
MAX_ITERATIONS = 1000000   # Convergence iteration for test
MAX_TOLERANCE = 1.0e-6     # Convergence tolerance for test



#%% ------------------------------------------
# Build model
ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 3)

# Frame location (in bays) and starting coordinate point
FrameStarts = 0.0
XCoords = FrameStarts
BayXCoords = [FrameStarts]
for y in range(BayTotal):
    XCoords += FrameLength[y]
    BayXCoords.append(XCoords)

#%% ------------------------------------------
# Basic material models

# Stiff material
StiffMatTagBCJ = 150000001
ops.uniaxialMaterial('Elastic', StiffMatTagBCJ, 10000000000000.0)  # stiffness

# Beam-to-column connection materials (Pinching4)
BC_tag = 115
for story in range(1, Nstory + 1):
    for grid in range(1, len(BayXCoords) + 1):
        for bay in range(1, len(BayXCoords)):
            # Determine if this combination is used; only certain combos exist in original Tcl
            # We'll replicate the exact conditions from the Tcl if/elif chain
            if (story == 1 and grid == 1 and bay == 1) or \
               (story == 1 and grid == 2 and bay == 1) or \
               (story == 1 and grid == 2 and bay == 2) or \
               (story == 1 and grid == 3 and bay == 2) or \
               (story == 2 and grid == 1 and bay == 1) or \
               (story == 2 and grid == 2 and bay == 1) or \
               (story == 2 and grid == 2 and bay == 2) or \
               (story == 2 and grid == 3 and bay == 2):
                matTag = 1000000 * BC_tag + 10000 * story + 100 * grid + bay
                # Common parameters
                xMy_1P = 785200.0
                xMy_2P = 7.31241e6
                xMy_3P = 2.3432e7
                xMy_4P = 2.4691e7
                xth_1P = 0.000241
                xth_2P = 0.0521
                xth_3P = 0.1822
                xth_4P = 0.2901

                xMy_1N = -xMy_1P
                xMy_2N = -xMy_2P
                xMy_3N = -xMy_3P
                xMy_4N = -xMy_4P
                xth_1N = -xth_1P
                xth_2N = -xth_2P
                xth_3N = -xth_3P
                xth_4N = -xth_4P

                rDispP = 0.7
                rForceP = 0.155
                uForceP = 0.01012
                rDispN = -rDispP
                rForceN = -rForceP
                uForceN = -uForceP

                gK1 = 1.0
                gK2 = 0.5
                gK3 = 0.5
                gK4 = 0.5
                gKLim = 0.01
                gD1 = 0.5
                gD2 = 0.5
                gD3 = 1.0
                # gD4 depends on specific case
                if story == 2 and grid == 2 and bay == 2:
                    gD4 = 0.82
                elif story == 2 and grid == 3 and bay == 2:
                    gD4 = 0.82
                else:
                    gD4 = 0.8
                gDLim = 0.2
                gF1 = 1.0
                gF2 = 0.0
                gF3 = 1.0
                gF4 = 1.0
                gFLim = 0.01
                gE = 10.0
                damage = "energy"

                ops.uniaxialMaterial('Pinching4', matTag,
                                     xMy_1P, xth_1P, xMy_2P, xth_2P, xMy_3P, xth_3P, xMy_4P, xth_4P,
                                     xMy_1N, xth_1N, xMy_2N, xth_2N, xMy_3N, xth_3N, xMy_4N, xth_4N,
                                     rDispP, rForceP, uForceP,
                                     rDispN, rForceN, uForceN,
                                     gK1, gK2, gK3, gK4, gKLim,
                                     gD1, gD2, gD3, gD4, gDLim,
                                     gF1, gF2, gF3, gF4, gFLim,
                                     gE, damage)
                # INFO LINK: https://opensees.berkeley.edu/OpenSees/manuals/usermanual/1152.htm

# Column-to-base connection materials (Pinching4)
CF_tagCode = 116
for story in range(0, Nstory):
    for grid in range(1, len(BayXCoords) + 1):
        # Original Tcl only defined for certain (story, grid) combos
        if (story == 0 and grid in [1, 2, 3]) or (story == 1 and grid in [1, 2, 3]):
            matTag = 1000000 * CF_tagCode + 10000 * story + 100 * grid + 1
            # Parameters
            xMy_1P = 786500.0
            xMy_2P = 6.12074e7
            xMy_3P = 6.3224e7
            xMy_4P = 5.4712e7
            xth_1P = 0.000137
            xth_2P = 0.0416
            xth_3P = 0.1709
            xth_4P = 0.2595

            xMy_1N = -xMy_1P
            xMy_2N = -xMy_2P
            xMy_3N = -xMy_3P
            xMy_4N = -xMy_4P
            xth_1N = -xth_1P
            xth_2N = -xth_2P
            xth_3N = -xth_3P
            xth_4N = -xth_4P

            rDispP = 0.7
            rForceP = 0.155
            uForceP = 0.01012
            rDispN = -rDispP
            rForceN = -rForceP
            uForceN = -uForceP

            gK1 = 1.0
            gK2 = 0.5
            gK3 = 0.5
            gK4 = 0.5
            gKLim = 0.01
            gD1 = 0.5
            gD2 = 0.5
            gD3 = 1.0
            gD4 = 0.8
            gDLim = 0.28
            gF1 = 1.0
            gF2 = 0.0
            gF3 = 1.0
            gF4 = 1.0
            gFLim = 0.01
            gE = 10.0
            damage = "energy"

            ops.uniaxialMaterial('Pinching4', matTag,
                                 xMy_1P, xth_1P, xMy_2P, xth_2P, xMy_3P, xth_3P, xMy_4P, xth_4P,
                                 xMy_1N, xth_1N, xMy_2N, xth_2N, xMy_3N, xth_3N, xMy_4N, xth_4N,
                                 rDispP, rForceP, uForceP,
                                 rDispN, rForceN, uForceN,
                                 gK1, gK2, gK3, gK4, gKLim,
                                 gD1, gD2, gD3, gD4, gDLim,
                                 gF1, gF2, gF3, gF4, gFLim,
                                 gE, damage)

# Shear connection material (stiff)
Sheari_tag_par = 117000011
ops.uniaxialMaterial('Elastic', Sheari_tag_par, 10000000000000.0)

#%% ------------------------------------------
# Model Columns

ops.geomTransf('Linear', 1)

CNODE_ID = []
NODE_X = []
NODE_Y = []
CSupportNode = []
GridID = []

for story in range(0, Nstory + 1):
    for grid in range(1, len(BayXCoords) + 1):
        if story == 0:
            GridID.append(grid)
        nodeTag = 101000000 + story * 10000 + grid * 100 + 1
        X = BayXCoords[grid - 1]
        if story == 0:
            Y = 0.0
        else:
            Y = Story1Height + (story - 1) * StoryHeight
        ops.node(nodeTag, X, Y)
        CNODE_ID.append(nodeTag)
        NODE_X.append(X)
        NODE_Y.append(Y)

        if story > 0:
            Ec = 9542.0      # [N/mm^2]
            Ac = 64400.0     # [mm^2]
            Izc = 2.839e8    # [mm^4]
            eleTag = 111000000 + story * 10000 + grid * 100 + 1
            iNode = 101000000 + (story - 1) * 10000 + grid * 100 + 1
            jNode = 101000000 + story * 10000 + grid * 100 + 1
            ops.element('elasticBeamColumn', eleTag, iNode, jNode, Ac, Ec, Izc, 1)
        else:
            CSupportNode.append(nodeTag)

StoryNodes = []
for story in range(1, Nstory + 1):
    StoryNodes.append(101000000 + story * 10000 + 1 * 100 + 1)

# Define foundation connections (zeroLength)
NODE_ID = []
NODE_X2 = []
NODE_Y2 = []
C2SupportNode = []
CFJoint_ID = []
CFJoint_I = []
CFJoint_J = []

for i in range(len(CSupportNode)):
    grid_id = GridID[i]
    eleTag = 114000000 + 0 * 10000 + grid_id * 100 + 1
    jNode = 101000000 + 0 * 10000 + grid_id * 100 + 1
    iNode = 101000000 + 0 * 10000 + grid_id * 100 + 0
    coords = ops.nodeCoord(jNode)
    ops.node(iNode, coords[0], coords[1])
    NODE_ID.append(iNode)
    NODE_X2.append(coords[0])
    NODE_Y2.append(coords[1])
    ops.fix(iNode, 1, 1, 1)
    C2SupportNode.append(iNode)
    # material tags for dir 1,2: Sheari_tag_par; dir 3: Pinching4 for CF
    matCF = 1000000 * CF_tagCode + 10000 * 0 + 100 * (i + 1) + 1
    ops.element('zeroLength', eleTag, iNode, jNode,
                '-mat', Sheari_tag_par, Sheari_tag_par, matCF,
                '-dir', 1, 2, 3)
    CFJoint_ID.append(eleTag)
    CFJoint_I.append(iNode)
    CFJoint_J.append(jNode)

iSupportNode = C2SupportNode

#%% ------------------------------------------
# Model Beams

BLNODE_ID = []
BRNODE_ID = []
NODE_X3 = []
NODE_Y3 = []

for story in range(1, Nstory + 1):
    for grid in range(1, len(BayXCoords) + 1):
        X = BayXCoords[grid - 1]
        Y = Story1Height + (story - 1) * StoryHeight
        if grid == 1:
            nodeTag = 101000000 + story * 10000 + grid * 100 + 11
            ops.node(nodeTag, X, Y)
            BLNODE_ID.append(nodeTag)
            NODE_X3.append(X)
            NODE_Y3.append(Y)
        elif grid == len(BayXCoords):
            nodeTag = 101000000 + story * 10000 + (grid - 1) * 100 + 12
            ops.node(nodeTag, X, Y)
            BRNODE_ID.append(nodeTag)
            NODE_X3.append(X)
            NODE_Y3.append(Y)
        else:
            nodeTag = 101000000 + story * 10000 + grid * 100 + 11
            ops.node(nodeTag, X, Y)
            BLNODE_ID.append(nodeTag)
            NODE_X3.append(X)
            NODE_Y3.append(Y)
            nodeTag2 = 101000000 + story * 10000 + (grid - 1) * 100 + 12
            ops.node(nodeTag2, X, Y)
            BRNODE_ID.append(nodeTag2)
            NODE_X3.append(X)
            NODE_Y3.append(Y)

        if grid > 1:
            Eb = 9542.0
            Ab = 50400.0
            Izb = 1.361e8
            eleTag = 112000000 + story * 10000 + (grid - 1)
            iNode = 101000000 + story * 10000 + (grid - 1) * 100 + 11
            jNode = 101000000 + story * 10000 + (grid - 1) * 100 + 12
            ops.element('elasticBeamColumn', eleTag, iNode, jNode, Ab, Eb, Izb, 1)

# Beam-to-column connections (zeroLength)
BCJ_ID = []
BCJ_I = []
BCJ_J = []

if BCJType == "Hysteresis":
    for story in range(1, Nstory + 1):
        r = len(BayXCoords)
        for grid in range(1, r + 1):
            for bay in range(1, r):
                if grid == 1 and bay == 1:
                    eleTag = 113000000 + story * 10000 + bay * 100 + 2
                    iNode = 101000000 + story * 10000 + grid * 100 + 1
                    jNode = 101000000 + story * 10000 + bay * 100 + 11
                    matTag = 1000000 * BC_tag + 10000 * story + 100 * grid + bay
                    ops.element('zeroLength', eleTag, iNode, jNode,
                                '-mat', StiffMatTagBCJ, StiffMatTagBCJ, matTag,
                                '-dir', 1, 2, 3)
                    BCJ_ID.append(eleTag)
                    BCJ_I.append(iNode)
                    BCJ_J.append(jNode)
                elif grid == r and bay == r - 1:
                    eleTag = 113000000 + story * 10000 + bay * 100 + 3
                    iNode = 101000000 + story * 10000 + bay * 100 + 12
                    jNode = 101000000 + story * 10000 + grid * 100 + 1
                    matTag = 1000000 * BC_tag + 10000 * story + 100 * grid + bay
                    ops.element('zeroLength', eleTag, iNode, jNode,
                                '-mat', StiffMatTagBCJ, StiffMatTagBCJ, matTag,
                                '-dir', 1, 2, 3)
                    BCJ_ID.append(eleTag)
                    BCJ_I.append(iNode)
                    BCJ_J.append(jNode)
                elif r > 2 and grid == 2 and bay == 1:
                    eleTag = 113000000 + story * 10000 + bay * 100 + 3
                    iNode = 101000000 + story * 10000 + bay * 100 + 12
                    jNode = 101000000 + story * 10000 + grid * 100 + 1
                    matTag = 1000000 * BC_tag + 10000 * story + 100 * grid + bay
                    ops.element('zeroLength', eleTag, iNode, jNode,
                                '-mat', StiffMatTagBCJ, StiffMatTagBCJ, matTag,
                                '-dir', 1, 2, 3)
                    BCJ_ID.append(eleTag)
                    BCJ_I.append(iNode)
                    BCJ_J.append(jNode)
                elif r > 2 and grid == 2 and bay == 2:
                    eleTag = 113000000 + story * 10000 + bay * 100 + 2
                    iNode = 101000000 + story * 10000 + grid * 100 + 1
                    jNode = 101000000 + story * 10000 + bay * 100 + 11
                    matTag = 1000000 * BC_tag + 10000 * story + 100 * grid + bay
                    ops.element('zeroLength', eleTag, iNode, jNode,
                                '-mat', StiffMatTagBCJ, StiffMatTagBCJ, matTag,
                                '-dir', 1, 2, 3)
                    BCJ_ID.append(eleTag)
                    BCJ_I.append(iNode)
                    BCJ_J.append(jNode)
                elif r > 3 and grid == 3 and bay == 2:
                    eleTag = 113000000 + story * 10000 + bay * 100 + 3
                    iNode = 101000000 + story * 10000 + bay * 100 + 12
                    jNode = 101000000 + story * 10000 + grid * 100 + 1
                    matTag = 1000000 * BC_tag + 10000 * story + 100 * grid + bay
                    ops.element('zeroLength', eleTag, iNode, jNode,
                                '-mat', StiffMatTagBCJ, StiffMatTagBCJ, matTag,
                                '-dir', 1, 2, 3)
                    BCJ_ID.append(eleTag)
                    BCJ_I.append(iNode)
                    BCJ_J.append(jNode)
else:  # Rigid
    for story in range(1, Nstory + 1):
        r = len(BayXCoords)
        for grid in range(1, r + 1):
            for bay in range(1, r):
                if grid == 1 and bay == 1:
                    eleTag = 113000000 + story * 10000 + bay * 100 + 2
                    iNode = 101000000 + story * 10000 + grid * 100 + 1
                    jNode = 101000000 + story * 10000 + bay * 100 + 11
                    ops.element('zeroLength', eleTag, iNode, jNode,
                                '-mat', StiffMatTagBCJ, StiffMatTagBCJ, StiffMatTagBCJ,
                                '-dir', 1, 2, 3)
                    BCJ_ID.append(eleTag)
                    BCJ_I.append(iNode)
                    BCJ_J.append(jNode)
                elif grid == r and bay == r - 1:
                    eleTag = 113000000 + story * 10000 + bay * 100 + 3
                    iNode = 101000000 + story * 10000 + bay * 100 + 12
                    jNode = 101000000 + story * 10000 + grid * 100 + 1
                    ops.element('zeroLength', eleTag, iNode, jNode,
                                '-mat', StiffMatTagBCJ, StiffMatTagBCJ, StiffMatTagBCJ,
                                '-dir', 1, 2, 3)
                    BCJ_ID.append(eleTag)
                    BCJ_I.append(iNode)
                    BCJ_J.append(jNode)
                elif r > 2 and grid == 2 and bay == 1:
                    eleTag = 113000000 + story * 10000 + bay * 100 + 3
                    iNode = 101000000 + story * 10000 + bay * 100 + 12
                    jNode = 101000000 + story * 10000 + grid * 100 + 1
                    ops.element('zeroLength', eleTag, iNode, jNode,
                                '-mat', StiffMatTagBCJ, StiffMatTagBCJ, StiffMatTagBCJ,
                                '-dir', 1, 2, 3)
                    BCJ_ID.append(eleTag)
                    BCJ_I.append(iNode)
                    BCJ_J.append(jNode)
                elif r > 2 and grid == 2 and bay == 2:
                    eleTag = 113000000 + story * 10000 + bay * 100 + 2
                    iNode = 101000000 + story * 10000 + grid * 100 + 1
                    jNode = 101000000 + story * 10000 + bay * 100 + 11
                    ops.element('zeroLength', eleTag, iNode, jNode,
                                '-mat', StiffMatTagBCJ, StiffMatTagBCJ, StiffMatTagBCJ,
                                '-dir', 1, 2, 3)
                    BCJ_ID.append(eleTag)
                    BCJ_I.append(iNode)
                    BCJ_J.append(jNode)
                elif r > 3 and grid == 3 and bay == 2:
                    eleTag = 113000000 + story * 10000 + bay * 100 + 3
                    iNode = 101000000 + story * 10000 + bay * 100 + 12
                    jNode = 101000000 + story * 10000 + grid * 100 + 1
                    ops.element('zeroLength', eleTag, iNode, jNode,
                                '-mat', StiffMatTagBCJ, StiffMatTagBCJ, StiffMatTagBCJ,
                                '-dir', 1, 2, 3)
                    BCJ_ID.append(eleTag)
                    BCJ_I.append(iNode)
                    BCJ_J.append(jNode)

# Control node
IDctrlNode = 101000000 + Nstory * 10000 + 1 * 100 + 1

#%% ------------------------------------------
# Display (optional) - not necessary for analysis but we can keep a simple recorder view
# In openseespy, display is not directly available; we skip or use ops.plotModel()

#%% ------------------------------------------
# Recorders

baseDir = os.getcwd()
dataDir = "SALAR_RESULTS"
if not os.path.exists(dataDir):
    os.mkdir(dataDir)
resultsDir = os.path.join(baseDir, dataDir)

# Global recorders
ops.recorder('Node', '-file', os.path.join(resultsDir, 'DRoof.txt'),
             '-node', IDctrlNode, '-dof', 1, 'disp')
# Record reactions at the two base nodes (foundation nodes)
# The base nodes are 101000100 and 101000200 (grid 1 and 2, story 0)
base_nodes = [101000100, 101000200] if len(BayXCoords) == 2 else [101000100, 101000200, 101000300]  # adjust for more bays
ops.recorder('Node', '-file', os.path.join(resultsDir, 'RBase.txt'),
             '-node', *base_nodes, '-dof', 1, 'reaction')

#%% ------------------------------------------
# Analysis

# Initial gravity analysis? Tcl does analyze 10 steps with LoadControl 0.1 (no loads defined)
# This is typical for gravity, but no loads are applied; we can do a dummy analysis.
ops.constraints('Transformation')
ops.numberer('RCM')
ops.system('FullGeneral')
ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS, 1)
ops.algorithm('ModifiedNewton')
ops.integrator('LoadControl', 0.1)
ops.analysis('Static')
ops.analyze(10)

ops.loadConst('-time', 0.0)

# Pushover analysis
TS_TAG = 2
PATT_TAG = 2
ops.timeSeries('Linear', TS_TAG)
ops.pattern('Plain', PATT_TAG, TS_TAG)    
ops.load(IDctrlNode, 1.0, 0.0, 0.0)

ops.constraints('Transformation')
ops.numberer('Plain')
ops.system('FullGeneral')
ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS, 1)
ops.algorithm('ModifiedNewton')

initialDisp = 0.0
targetDisp = Nstory * 260.0  # [mm]
dispIncr = 1.0
numSteps = int((targetDisp - initialDisp) / dispIncr)

for step in range(1, numSteps + 1):
    ops.integrator('DisplacementControl', IDctrlNode, 1, dispIncr)
    ops.analysis('Static')
    OK = ops.analyze(1)
    S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS

print("Pushover analysis completed.")

#%% ------------------------------------------
# Plot results

# Read roof displacement
disp_file = os.path.join(resultsDir, 'DRoof.txt')
if os.path.exists(disp_file):
    disp_data = np.loadtxt(disp_file)
    if disp_data.ndim == 1:
        disp_data = disp_data.reshape(-1, 1)
    roof_disp = disp_data[:, 0]  # first column (only dof 1)
else:
    roof_disp = None
    print("DRoof.out not found.")

# Read base reactions
reac_file = os.path.join(resultsDir, 'RBase.txt')
if os.path.exists(reac_file):
    reac_data = np.loadtxt(reac_file)
    if reac_data.ndim == 1:
        reac_data = reac_data.reshape(1, -1)
    # Sum reactions from all base nodes (each node contributes one column)
    base_shear = np.sum(reac_data, axis=1)  # sum across columns (nodes)
else:
    base_shear = None
    print("RBase.out not found.")

# Plot if data available
if roof_disp is not None and base_shear is not None:
    # Ensure same length (might have one extra step due to recorder timing)
    min_len = min(len(roof_disp), len(base_shear))
    roof_disp = roof_disp[:min_len]
    base_shear = base_shear[:min_len]

    plt.figure(figsize=(8, 6))
    plt.plot(roof_disp, base_shear, 'b-', linewidth=2)
    plt.xlabel('Roof Displacement (mm)')
    plt.ylabel('Base Shear (N)')
    plt.title('Pushover Curve')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(resultsDir, 'pushover_curve.png'), dpi=300)
    plt.show()
else:
    print("Could not plot due to missing data.")
    
    
S01.PLOT_2D_FRAME(deformed_scale=1.0) 


import pandas as pd
import matplotlib.pyplot as plt

# LOAD EXCEL DATA
FILE_PATH = "RESULTS.xlsx"

df = pd.read_excel(FILE_PATH, sheet_name="SALAR")

# EXTRACT DATA
POINT = df["Point"]
DISPLACEMENT = df["Displacement (mm)"]         # mm
LOAD = df["Load (kN)"] * 1000                  # kN → N

# DISPLAY DATA
print(df.head())
print(f"Number of points: {len(df)}")

# PLOT HYSTERESIS RESULTS
plt.figure(figsize=(12, 8))

plt.plot(+1 *  np.array(roof_disp),
         -1 *  np.array(base_shear),
         'r-', linewidth=5, label='Numerical Data - Hystersis Curve')

plt.plot(
    DISPLACEMENT,
    LOAD,
    color="blue",
    linewidth=1.5,
    label="Experimental Data - Outer Boundary Hystersis Curve"
)

plt.axhline(0, color="gray", linewidth=0.8)
plt.axvline(0, color="gray", linewidth=0.8)

plt.xlabel("Displacement (m)")
plt.ylabel("Load (N)")
plt.title("NUMERICAL AND EXPERIMENAL DATA")

plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()   