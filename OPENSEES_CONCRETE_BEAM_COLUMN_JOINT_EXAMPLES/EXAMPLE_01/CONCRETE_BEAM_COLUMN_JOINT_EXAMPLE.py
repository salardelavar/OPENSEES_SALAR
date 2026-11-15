#######################################################################################################
#                    >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                 #
# NONLINEAR CYCLIC ANALYSIS OF A REINFORCED CONCRETE BEAM-COLUMN JOINT SUB-ASSEMBLAGE USING OPENSEESPY#
#-----------------------------------------------------------------------------------------------------#
#   ORIGINAL TCL EXAMPLE:                                                                             #
#   https://opensees.berkeley.edu/wiki/index.php/BeamColumnJointExample                               #
#   created : NM (nmitra@u.washington.edu) dated: Dec. 2002                                           #
#   Rewritted by Salar Delavar Ghashghaei (salar.d.ghashghaei@gmail.com) dated: Dec. 2026             #     
#   UNITS: N, mm                                                                                      #                                        
#######################################################################################################
"""
 This script is an OpenSeesPy version of the OpenSees wiki example BeamColumnJointExample
 (https://opensees.berkeley.edu/wiki/index.php/BeamColumnJointExample). The original is written 
 in Tcl as `PR1.tcl`, `procMKPC.tcl`, `procUniaxialPinching.tcl` and `procRC.tcl`.
 The example models a reinforced-concrete beam-column cruciform. The wiki describes it as carrying
 constant gravity load at nodes 4 and 7 and cyclic lateral displacement control at node 10. The beams
 and columns are fiber-based nonlinear elements (elements 1–6), and the joint (element 7) is
 a `beamColumnJoint` element with 13 components, 9 of which use the Pinching4 material. In the script,
 `MKPC_FUN` computes the confined-concrete properties with the Modified Kent–Park model. It converts to
 psi and inches internally and returns MPa and mm/mm. `FUN_uniaxial_pinching` wraps the `Pinching4` command,
 which is used for the joint shear panel. `MODEL_FUN` builds the whole model in N, mm and MPa. This includes
 the nodes, the concrete and steel materials, the column and beam fiber sections, the force-based elements,
 the 8 bar-slip materials, the shear-panel material and the supports.

The second part of `MODEL_FUN` runs the analysis. It applies the gravity pattern, but the load increment
 `GRAVITY_DLAMBDA` is 0.0, so the 55 kN loads are never actually applied, exactly as in the Tcl. It then
 adds a unit lateral load at node 10 and drives the column top through the reversed-cyclic displacement
 history in `peakpts`, using displacement control. If a step fails, it tries KrylovNewton, NewtonLineSearch
 and Newton in turn, and uses ArcLength only as a last resort. All fallbacks use the same tolerance of
 1e-6 with 500 iterations. `new_record`, `record_step` and `save_record` store the load, displacement,
 shear-panel and bar-slip histories and write them to `.txt` files. `PLOTS_FUN` saves and shows the
 load–displacement, shear-panel and bar-slip plots in N and mm, and `main()` reads the command-line
 options and runs everything.
"""

import argparse
import os
import time as TI
import numpy as np
import matplotlib.pyplot as plt
import openseespy.opensees as ops
import opsvis as opsv
import PLOT_2D as S01
import ANALYSIS_FUNCTION as S02
# ---------------------------------------------
# PEER REPORT: A Beam-Column Joint Model for Simulating the Earthquake Response of Reinforced Concrete Frames
'https://peer.berkeley.edu/sites/default/files/0310_l._lowes_n._mitra_a._altoontash.pdf'  
# PhD Thesis: An analytical study of reinforced concrete beam-column joint behavior under seismic loading
'https://faculty.washington.edu/lowes/lowes/studentTheses/2007_NilanjanMitra_PHDthesis.pdf'  

# ---------------------------------------------
#%% Modified Kent-Park confined concrete envelope
# --------------------------
def MKPC_FUN(CUnconfFc, CUnconfEc, Y, Z, Cov, TSspace, TSlength, TSFy, TSarea,
              Strfactor, Lenfactor):
    """
    MKPC_FUN: Computes Modified Kent–Park confined and unconfined concrete envelope
    parameters (peak and ultimate stress and strain) from the section geometry and
    hoop reinforcement, converting to psi and inches internally and returning MPa and mm/mm.
    """
    
    CUnconfEcu = -0.004

    SecWid = Lenfactor * Z
    SecDep = Lenfactor * Y
    cover = Lenfactor * Cov

    UFc = -Strfactor * CUnconfFc          # psi, positive
    Ue0 = -CUnconfEc
    Uecu = -CUnconfEcu

    hoopSpc = Lenfactor * TSspace
    hoopLngth = Lenfactor * TSlength
    hoopFy = Strfactor * TSFy
    hoopArea = TSarea * Lenfactor * Lenfactor

    # ratio of hoop steel volume to confined core volume (to outside of hoops)
    rhoS = (hoopLngth * hoopArea) / ((SecWid - 2 * cover) * (SecDep - 2 * cover) * hoopSpc)

    b = SecWid - 2 * cover
    temp = b / hoopSpc

    e50u = (3 + 0.002 * UFc) / (UFc - 1000)
    e50h = 3 * rhoS * temp ** 0.5 / 4

    Zm = 0.5 * (UFc - 1000) / (3 + 0.002 * UFc)
    Zc = 0.5 / (e50u + e50h - Ue0)

    K = 1 + rhoS * hoopFy / UFc

    UFcu = -UFc * (1 - Zm * (Uecu - Ue0)) / Strfactor      # [MPa] Unconfined crushing stress
    Ce0 = -K * Ue0                                         # [mm/mm] Confined peak strain
    CFc = -K * UFc / Strfactor                             # [MPa] Confined peak stress
    CFcu = 0.2 * CFc                                       # [MPa] Confined residual stress
    Cecu = -(0.8 / Zc - Ce0)                               # [mm/mm] Confined ultimate strain
    
    # Returns [Fc, Ec, Fcu, Ecu (unconfined), Fc, Ec, Fcu, Ecu (confined)] in MPa (N/mm2).
    DATA = [CUnconfFc, CUnconfEc, UFcu, CUnconfEcu, CFc, Ce0, CFcu, Cecu]
    return DATA

#%% Pinching4 wrapper
# --------------------------
def FUN_uniaxial_pinching(tag, pEnvStress, nEnvStress, pEnvStrain, nEnvStrain,
                           rDisp, rForce, uForce, gammaK, gammaD, gammaF, gammaE,
                           damage):
    """
    FUN_uniaxial_pinching: Wraps the OpenSees Pinching4 command, packing the envelope,
    pinching, and damage arguments into one call to define the joint shear-panel material.
    """
    args = []
    for i in range(4):
        args += [pEnvStress[i], pEnvStrain[i]]
    for i in range(4):
        args += [nEnvStress[i], nEnvStrain[i]]
    args += [rDisp[0], rForce[0], uForce[0], rDisp[1], rForce[1], uForce[1]]
    args += list(gammaK) + list(gammaD) + list(gammaF) + [gammaE, damage]
    ops.uniaxialMaterial("Pinching4", tag, *args)
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/Pinching4.html



#%% Model
# --------------------------
def MODEL_FUN(incre, node, dof, peakpts, rec, RECORD_SPRINGS = True , PLOT_SECTIONS = True):
    """
    MODEL_FUN: Builds the whole beam-column joint model (nodes, materials, fiber sections,
    elements, bar-slip springs, shear panel, supports), then runs the gravity step and the
    displacement-controlled reversed-cyclic analysis.
    """
    MAX_ITERATIONS = 500    # Maximum number of iterations
    MAX_TOLERANCE = 1.0e-6  # Specified tolerance for convergence
    
    # Reset model
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)

    Strfactor = 145.0
    Lenfactor = 1 / 25.4

    # Geometry for Columns 
    colY, colZ, colCov = 406.0, 305.0, 43.0 # [mm] Column Depth, Width, Cover
    # Geometry for Beams  (mm)
    bmY, bmZ = 457.0, 229.0                 # [mm] Bolumn Depth, Width
    bmCov1, bmCov2 = 42.0, 33.0
    bmCov = bmCov1

    JointWidth, JointHeight, JointDepth = colY, bmY, colZ # Joint size: width = column depth, height = beam depth, depth (out of plane) = column width.
    BeamLengthIn, BeamLengthOut, ColumnLengthClear = 645.0, 1271.0, 1008.0 # Beam length near the joint (inner), beam length beyond it (outer), and clear column length, all in mm.
    JointVolume = JointWidth * JointHeight * JointDepth   # Joint volume (mm³), used later to convert shear-panel stresses into moments.

    # Column section materials
    CUnconfFc, CUnconfEc = -45.9, -0.002   # Unconfined concrete peak stress (MPa, negative = compression) and its strain.
    CTSspace, CTSlength, CTSFy, CTSarea = 60.0, 1853.53, 282.0, 28.3 # Column hoop data: spacing, length of one hoop, yield stress, and bar area.
    CFy, CEs, CsHratio, CAs = 498.0, 196600.0, 0.004216, 201.06 # Column longitudinal steel: yield stress, elastic modulus, strain-hardening ratio, and area per bar.
    cp = MKPC_FUN(CUnconfFc, CUnconfEc, colY, colZ, colCov,
                   CTSspace, CTSlength, CTSFy, CTSarea, Strfactor, Lenfactor) # Runs the Modified Kent–Park model to get the column’s confined and unconfined concrete parameters.
    CUnconfFcu, CUnconfEcu = cp[2], cp[3] # Unconfined crushing stress and ultimate strain from the result list.
    CConfFc, CConfEc, CConfFcu, CConfEcu = cp[4], cp[5], cp[6], cp[7] # Confined peak stress and strain, plus residual stress and ultimate strain.

    # Beam section materials
    BUnconfFc, BUnconfEc = -45.9, -0.002  # Beam unconfined concrete peak stress and strain (same values as the column).
    BTSspace, BTSlength, BTSFy, BTSarea = 80.0, 1036.0, 282.0, 28.3 # Beam stirrup data: spacing, length of one stirrup, yield stress, and bar area.
    BFy, BEs, BAs, BsHratio = 294.0, 210400.0, 201.06, 0.002322     # Beam stirrup data: spacing, length of one stirrup, yield stress, and bar area.
    bp = MKPC_FUN(BUnconfFc, BUnconfEc, bmY, bmZ, bmCov,
                   BTSspace, BTSlength, BTSFy, BTSarea, Strfactor, Lenfactor) # Runs the Modified Kent–Park model for the beam section.
    BUnconfFcu, BUnconfEcu = bp[2], bp[3] # Beam unconfined crushing stress and ultimate strain.
    BConfFc, BConfEc, BConfFcu, BConfEcu = bp[4], bp[5], bp[6], bp[7] # Beam confined peak stress and strain, plus residual stress and ultimate strain.

    # Bar-slip data
    bs_fc, bs_fs, bs_es, bs_fsu, bs_dbar = -BUnconfFc, BFy, BEs, 434.0, 16.0 # Beam bar-slip inputs: concrete strength (made positive), bar yield stress, elastic modulus, ultimate stress (434 MPa), and bar diameter (16 mm).
    bs_esh = BsHratio * BEs            # Beam bar strain-hardening modulus (ratio × elastic modulus).
    bs_wid, bs_dep = colZ, bmY         # Width and depth of the member the beam bars are anchored in: the joint width (column width) and the beam depth.
    bsT_nbars, bsB_nbars = 5, 2        # Number of top and bottom beam bars.
    bs_ljoint = colY                   # Anchorage length of the beam bars through the joint, which is the column depth.

    cs_fc, cs_fs, cs_es, cs_fsu, cs_dbar = -CUnconfFc, CFy, CEs, 660.0, 16.0 # Column bar-slip inputs: concrete strength (positive), bar yield stress, elastic modulus, ultimate stress (660 MPa), and bar diameter (16 mm).
    cs_esh = CsHratio * CEs            # Column bar strain-hardening modulus.
    cs_wid, cs_dep = colZ, colY        # Column width and depth, used as the section dimensions for the column bar-slip springs.
    cs_nbars = 3                       # Number of column bars on each side that slip.
    cs_ljoint = bmY                    # Anchorage length of the column bars through the joint, which is the beam depth.

    # Nodes
    yb = ColumnLengthClear + JointHeight / 2 # Height of the beam centerline above the base: the clear column length plus half the joint height (beam depth), used as the y-coordinate of nodes 3 to 8.
    ops.node(1, 0.0, 0.0)
    ops.node(2, 0.0, ColumnLengthClear)
    ops.node(3, -BeamLengthOut - BeamLengthIn - JointWidth / 2, yb)
    ops.node(4, -BeamLengthIn - JointWidth / 2, yb)
    ops.node(5, -JointWidth / 2, yb)
    ops.node(6, JointWidth / 2, yb)
    ops.node(7, BeamLengthIn + JointWidth / 2, yb)
    ops.node(8, BeamLengthOut + BeamLengthIn + JointWidth / 2, yb)
    ops.node(9, 0.0, ColumnLengthClear + JointHeight)
    ops.node(10, 0.0, 2 * ColumnLengthClear + JointHeight)

    # Materials
    ops.uniaxialMaterial("Elastic", 1, 1.0e10)

    # Beam
    ops.uniaxialMaterial("Concrete01", 10, BUnconfFc, BUnconfEc, BUnconfFcu, BUnconfEcu)
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/Concrete01.html
    ops.uniaxialMaterial("Concrete01", 20, BConfFc, BConfEc, BConfFcu, BConfEcu)
    ops.uniaxialMaterial("Steel02", 30, BFy, BEs, BsHratio, 18.5, 0.925, 0.15, 0.0, 0.4, 0.0, 0.5)
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/steel02.html
    # Column
    ops.uniaxialMaterial("Concrete01", 40, CUnconfFc, CUnconfEc, CUnconfFcu, CUnconfEcu)
    ops.uniaxialMaterial("Concrete01", 50, CConfFc, CConfEc, CConfFcu, CConfEcu)
    ops.uniaxialMaterial("Steel02", 60, CFy, CEs, CsHratio, 18.5, 0.925, 0.15, 0.0, 0.4, 0.0, 0.5)
    COL_TAG, BE_TAG = 1, 2
    """
    # Fiber section: column (tag COL_TAG)
    z, y = colZ / 2.0, colY / 2.0
    ops.section("Fiber", COL_TAG)
    ops.patch("rect", 50, 8, 1, colCov - y, colCov - z, y - colCov, z - colCov)   # core
    ops.patch("rect", 40, 2, 1, -y, colCov - z, colCov - y, z - colCov)           # cover
    ops.patch("rect", 40, 2, 1, y - colCov, colCov - z, y, z - colCov)
    ops.patch("rect", 40, 8, 1, -y, -z, y, colCov - z)
    ops.patch("rect", 40, 8, 1, -y, z - colCov, y, z)
    ops.layer("straight", 60, 3, CAs, y - colCov, colCov - z, y - colCov, z - colCov)
    ops.layer("straight", 60, 2, CAs, 0.0, colCov - z, 0.0, z - colCov)
    ops.layer("straight", 60, 3, CAs, colCov - y, colCov - z, colCov - y, z - colCov)

    # Fiber section: beam (tag BE_TAG)
    z, y = bmZ / 2.0, bmY / 2.0
    ops.section("Fiber", 2)
    ops.patch("rect", 20, 8, 1, bmCov1 - y, bmCov1 - z, y - bmCov1, z - bmCov1)
    ops.patch("rect", 10, 2, 1, -y, bmCov1 - z, bmCov1 - y, z - bmCov1)
    ops.patch("rect", 10, 2, 1, y - bmCov1, bmCov1 - z, y, z - bmCov1)
    ops.patch("rect", 10, 8, 1, -y, -z, y, bmCov1 - z)
    ops.patch("rect", 10, 8, 1, -y, z - bmCov1, y, z)
    ops.layer("straight", 30, 3, BAs, y - bmCov1, bmCov1 - z, y - bmCov1, z - bmCov1)
    ops.layer("straight", 30, 2, BAs, y - bmCov1 - bmCov2, bmCov1 - z, y - bmCov1 - bmCov2, z - bmCov1)
    ops.layer("straight", 30, 2, BAs, bmCov1 - y, bmCov1 - z, bmCov1 - y, z - bmCov1)
    """
    # Fiber sections defined as opsvis lists
    # ---------------------------------------
    # Column section (tag COL_TAG)
    z, y = colZ / 2.0, colY / 2.0
    col_sec = [
        ['section', 'Fiber', COL_TAG, '-Gj', 1e6],
        # core (confined)
        ['patch', 'rect', 50, 8, 1, colCov - y, colCov - z, y - colCov, z - colCov],
        # cover (unconfined)
        ['patch', 'rect', 40, 2, 1, -y, colCov - z, colCov - y, z - colCov],
        ['patch', 'rect', 40, 2, 1, y - colCov, colCov - z, y, z - colCov],
        ['patch', 'rect', 40, 8, 1, -y, -z, y, colCov - z],
        ['patch', 'rect', 40, 8, 1, -y, z - colCov, y, z],
        # reinforcement
        ['layer', 'straight', 60, 3, CAs, y - colCov, colCov - z, y - colCov, z - colCov],
        ['layer', 'straight', 60, 2, CAs, 0.0, colCov - z, 0.0, z - colCov],
        ['layer', 'straight', 60, 3, CAs, colCov - y, colCov - z, colCov - y, z - colCov],
    ]
    
    # Beam section (tag BE_TAG)
    z, y = bmZ / 2.0, bmY / 2.0
    beam_sec = [
        ['section', 'Fiber', BE_TAG, '-Gj', 1e6],
        # core (confined)
        ['patch', 'rect', 20, 8, 1, bmCov1 - y, bmCov1 - z, y - bmCov1, z - bmCov1],
        # cover (unconfined)
        ['patch', 'rect', 10, 2, 1, -y, bmCov1 - z, bmCov1 - y, z - bmCov1],
        ['patch', 'rect', 10, 2, 1, y - bmCov1, bmCov1 - z, y, z - bmCov1],
        ['patch', 'rect', 10, 8, 1, -y, -z, y, bmCov1 - z],
        ['patch', 'rect', 10, 8, 1, -y, z - bmCov1, y, z],
        # reinforcement
        ['layer', 'straight', 30, 3, BAs, y - bmCov1, bmCov1 - z, y - bmCov1, z - bmCov1],
        ['layer', 'straight', 30, 2, BAs, y - bmCov1 - bmCov2, bmCov1 - z, y - bmCov1 - bmCov2, z - bmCov1],
        ['layer', 'straight', 30, 2, BAs, bmCov1 - y, bmCov1 - z, bmCov1 - y, z - bmCov1],
    ]
    
    # Create the sections in OpenSees from the lists
    opsv.fib_sec_list_to_cmds(col_sec)
    opsv.fib_sec_list_to_cmds(beam_sec)
    
    #%% Optional plotting
    # ----------------------------------
    if PLOT_SECTIONS:
        # plot_fiber_section picks colors by (matTag - 1), so remap the
        # model tags (10, 20, 30...) to 1 = cover, 2 = core, 3 = steel
        def _plot_copy(sec, cover, core, steel):
            tagmap = {cover: 1, core: 2, steel: 3}
            out = []
            for item in sec:
                item = list(item)
                if item[0] in ('patch', 'layer'):
                    item[2] = tagmap[item[2]]
                out.append(item)
            return out
    
        mcolors = ['lightgrey', 'tan', 'k']   # cover, core, steel
    
        opsv.plot_fiber_section(_plot_copy(col_sec, 40, 50, 60), matcolor=mcolors)
        plt.title('Column section (mm)')
        plt.axis('equal')
    
        opsv.plot_fiber_section(_plot_copy(beam_sec, 10, 20, 30), matcolor=mcolors)
        plt.title('Beam section (mm)')
        plt.axis('equal')
    
        plt.show()
    
    # Geometric transformations and beam-column elements
    ops.geomTransf("Linear", 1)
    ops.geomTransf("Linear", 2)
    #ops.geomTransf('PDelta', 1)
    #ops.geomTransf('Corotational', 1)
    # nonlinearBeamColumn (tag, iNode, jNode, nIP, secTag, transfTag) -> Lobatto force-based
    ops.beamIntegration("Lobatto", 1, 1, 5)
    ops.beamIntegration("Lobatto", 2, 2, 3)
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/Lobatto.html
    ops.element("forceBeamColumn", 1, 1, 2, 2, 1)     # column lower
    ops.element("forceBeamColumn", 2, 9, 10, 2, 1)    # column upper
    ops.element("forceBeamColumn", 3, 3, 4, 1, 2)     # beam left outer
    ops.element("forceBeamColumn", 4, 4, 5, 1, 2)     # beam left inner
    ops.element("forceBeamColumn", 5, 6, 7, 1, 2)     # beam right inner
    ops.element("forceBeamColumn", 6, 7, 8, 1, 2)     # beam right outer
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/ForceBeamColumn.html

    # Bar-slip materials (12 springs of the joint element)
    def beam_slip(tag, nbars, kind):
        # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/BarSlip.html
        ops.uniaxialMaterial("BarSlip", tag, bs_fc, bs_fs, bs_es, bs_fsu, bs_esh, bs_dbar,
                             bs_ljoint, nbars, bs_wid, bs_dep, "strong", kind)
    def col_slip(tag):
        ops.uniaxialMaterial("BarSlip", tag, cs_fc, cs_fs, cs_es, cs_fsu, cs_esh, cs_dbar,
                             cs_ljoint, cs_nbars, cs_wid, cs_dep, "strong", "column")
        
    beam_slip(21, bsB_nbars, "beamBot")
    beam_slip(22, bsB_nbars, "beamBot")
    beam_slip(31, bsT_nbars, "beamTop")
    beam_slip(32, bsT_nbars, "beamTop")
    for t in (41, 42, 43, 44):
        col_slip(t)

    # Shear panel material (Pinching4); stresses (MPa) x joint volume -> N.mm
    p1, p2, p3 = 2.1932, 4.0872, 4.4862
    p4 = p3 * 1e-3
    pEnvStrsp = [p * JointVolume for p in (p1, p2, p3, p4)]
    nEnvStrsp = [-p for p in pEnvStrsp]
    pEnvStnsp = [0.0002, 0.004465, 0.0131, 0.0269]
    nEnvStnsp = [-e for e in pEnvStnsp]
    FUN_uniaxial_pinching(
        5, pEnvStrsp, nEnvStrsp, pEnvStnsp, nEnvStnsp,
        rDisp=[0.25, 0.25], rForce=[0.15, 0.15], uForce=[0.0, 0.0],
        gammaK=[1.13364492409642, 0.0, 0.10111033064469, 0.0, 0.91652498468618],
        gammaD=[0.12, 0.0, 0.23, 0.0, 0.95],
        gammaF=[1.11, 0.0, 0.319, 0.0, 0.125],
        gammaE=10.0, damage="energy")

    # Beam-column joint element (anticlockwise nodes 2,6,9,5; 12 springs + shear panel)
    ops.element("beamColumnJoint", 7, 2, 6, 9, 5,
                41, 42, 1, 21, 31, 1, 43, 44, 1, 22, 32, 1, 5)
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/beamColumnJoint.html

    # Boundary conditions
    ops.fix(1, 1, 1, 0)
    ops.fix(2, 0, 0, 0)
    ops.fix(3, 0, 1, 0)
    ops.fix(4, 0, 0, 0)
    ops.fix(5, 0, 0, 0)
    ops.fix(6, 0, 0, 0)
    ops.fix(7, 0, 0, 0)
    ops.fix(8, 0, 1, 0)
    ops.fix(9, 0, 0, 0)
    ops.fix(10, 0, 0, 0)
    
    #%% Gravity stage and lateral load pattern
    # --------------------------
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    ops.load(4, 0.0, -55000.0, 0.0)
    ops.load(7, 0.0, -55000.0, 0.0)

    ops.system("ProfileSPD")
    ops.constraints("Plain")
    GRAVITY_DLAMBDA = 0.0 
    ops.integrator("LoadControl", GRAVITY_DLAMBDA, 1, 0.0, 0.0)
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/loadControl.html
    ops.test("NormDispIncr", MAX_TOLERANCE, MAX_ITERATIONS)
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/normDispIncr.html
    ops.algorithm("Newton")
    ops.numberer("RCM")
    ops.analysis("Static")
    ops.analyze(1)
    
    ops.loadConst("-time", 0.0)
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/loadConst.html
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(10, 1.0, 0.0, 0.0)
    
    numTimes = 150
    dU = peakpts[0] / incre
    dU0 = dU / 10000.0


    fallbacks = [("KrylovNewton", (), MAX_TOLERANCE),
                 ("NewtonLineSearch", (0.8,), MAX_TOLERANCE),
                 ("Newton", (), MAX_TOLERANCE)]

    def step(sign):
        # One displacement-control step with fallbacks, ArcLength last.
        ops.integrator("DisplacementControl", node, dof, sign * dU, 1, sign * dU, sign * dU)
        ok = ops.analyze(1)
        if ok != 0:
            for name, extra, tol in fallbacks:
                ops.algorithm(name, *extra)
                ops.test("NormDispIncr", MAX_TOLERANCE, MAX_ITERATIONS, 0)
                ok = ops.analyze(1)
                if ok == 0:
                    break
            ops.algorithm("Newton")
            ops.test("NormDispIncr", MAX_TOLERANCE, MAX_ITERATIONS, 0)
        if ok != 0:
            print("Displacement control failed ..... trying Arc-Length control")
            print("Current Displacement is", ops.nodeDisp(node, dof))
            ops.test("NormDispIncr", 1e-6, numTimes, 0)
            ops.integrator("ArcLength", dU0, 1.0)
            ops.analyze(1)
            ops.test("NormDispIncr", MAX_TOLERANCE, MAX_ITERATIONS, 0)
        record_step(rec, node, dof, RECORD_SPRINGS)
        return ok

    ops.integrator("DisplacementControl", node, dof, dU, 1, dU, dU)
    ops.analysis("Static")

    # first peak: +x, -2x, +x
    for n, sign in ((incre, +1), (2 * incre, -1), (incre, +1)):
        ops.integrator("DisplacementControl", node, dof, sign * dU, 1, sign * dU, sign * dU)
        for _ in range(n):
            ops.analyze(1)
            record_step(rec, node, dof, RECORD_SPRINGS)

    # remaining peaks
    for y in peakpts[1:]:
        dS = int(y / dU)
        ops.test("NormDispIncr", MAX_TOLERANCE, MAX_ITERATIONS, 0)
        ops.algorithm("Newton")
        ok = 0
        for n, sign in ((dS, +1), (2 * dS, -1), (dS, +1)):
            for _ in range(n):
                ok = step(sign)
        if ok == 0:
            print(f"analysis successful at {y} mm displacement")
        else:
            print(f"analysis could not proceed fine beyond {y} mm displacement")     
      
    
#%% Result recorder
# --------------------------
SLIP_RESPONSES = ["node2BarSlipB", "node2BarSlipT", "node4BarSlipB", "node4BarSlipT",
                  "node1BarSlipL", "node1BarSlipR", "node3BarSlipL", "node3BarSlipR"]
def new_record():
    """
    new_record: Returns an empty results dictionary for the load-displacement, shear-panel,
    and bar-slip histories.
    """
    DATA = {"top": [],                                 # [load_factor, disp]  (N, mm)
            "panel": [],                               # [shear panel moment, rotation]  (N.mm, rad)
            "slip": {k: [] for k in SLIP_RESPONSES}}   # [stress MPa, strain]
    return DATA


def record_step(rec, node, dof, record_springs=True):
    """
    record_step: Appends the current load factor and top displacement, plus the shear-panel
    and bar-slip responses if enabled, to the results dictionary.
    """
    rec["top"].append([ops.getLoadFactor(1), ops.nodeDisp(node, dof)])
    if record_springs:
        rec["panel"].append(ops.eleResponse(7, "shearpanel", "stressStrain"))
        # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/eleResponse.html
        for k in SLIP_RESPONSES:
            rec["slip"][k].append(ops.eleResponse(7, k, "stressStrain"))


def save_record(rec, prefix="PR1", folder="."):
    """
    save_record: Writes the recorded histories to .txt files with np.savetxt.
    """
    np.savetxt(os.path.join(folder, prefix + ".txt"), np.array(rec["top"]))
    np.savetxt(os.path.join(folder, prefix + "_Sp.txt"), np.array(rec["panel"]))
    for k, v in rec["slip"].items():
        np.savetxt(os.path.join(folder, f"{prefix}_{k}.txt"), np.array(v))


# Plots (N, mm, MPa)
# --------------------------
def PLOTS_FUN(rec, outdir="."):
    """
    PLOTS_FUN: Plots and saves the load-displacement, shear-panel, and bar-slip hysteresis
    curves, then displays them.
    """
    top = np.array(rec["top"])
    d, F = top[:, 1], top[:, 0]

    fig, ax = plt.subplots(figsize=(7, 5.5))
    ax.plot(d, F, lw=0.8, color="tab:blue")
    ax.axhline(0, color="k", lw=0.5)
    ax.axvline(0, color="k", lw=0.5)
    ax.set_xlabel("Column-top displacement [mm]")
    ax.set_ylabel("Applied lateral load [N]")
    ax.set_title("Lateral load vs top displacement (N, mm)")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "Load_displacement_N.png"), dpi=200)

    sp = np.array(rec["panel"])
    if sp.size:
        fig, ax = plt.subplots(figsize=(7, 5.5))
        ax.plot(sp[:, 1], sp[:, 0], lw=0.8, color="tab:red")
        ax.axhline(0, color="k", lw=0.5)
        ax.axvline(0, color="k", lw=0.5)
        ax.set_xlabel("Shear-panel deformation [rad]")
        ax.set_ylabel("Shear-panel moment [N.mm]")
        ax.set_title("Joint shear panel (Pinching4)")
        ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, "Shear_panel.png"), dpi=200)

        fig, axs = plt.subplots(2, 4, figsize=(15, 7))
        for ax, k in zip(axs.ravel(), SLIP_RESPONSES):
            v = np.array(rec["slip"][k])
            ax.plot(v[:, 1], v[:, 0], lw=0.7)
            ax.set_title(k, fontsize=9)
            ax.grid(alpha=0.3)
            ax.set_xlabel("Slip strain")
            ax.set_ylabel("Bar stress [MPa]")
        fig.suptitle("Bar-slip springs (stress in N/mm2)")
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, "Bar_slip.png"), dpi=200)

    plt.show()      # opens all figures in interactive windows


#%% Main Function
# --------------------------
def main():
    """
    main: Parses the command-line options, runs MODEL_FUN, saves the results,
    plots them, and prints the peak loads.
    """
    ap = argparse.ArgumentParser()
    ap.add_argument("--increment", type=int, default=10,
                    help="steps in the first 0.1 mm peak; dU = 0.1/increment mm")
    ap.add_argument("--peaks", type=float, nargs="+",
                    default=[0.1, 10, 10, 30, 30, 45, 45, 60, 60, 75, 75, 90, 90, 105, 105])
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--no-springs", action="store_true",
                    help="do not record shear-panel / bar-slip histories (faster)")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    rec = new_record()
    t0 = TI.time()
    
    RECORD_SPRINGS = True
    PLOT_SECTIONS = True
    MODEL_FUN(a.increment, 10, 1, a.peaks, rec, RECORD_SPRINGS, PLOT_SECTIONS)
    print(f"analysis time: {TI.time() - t0:.1f} s, steps recorded: {len(rec['top'])}")

    save_record(rec, "SALAR", a.outdir)
    PLOTS_FUN(rec, a.outdir)
    top = np.array(rec["top"])
    print(f"peak +load {top[:,0].max():.0f} N at {top[top[:,0].argmax(),1]:.2f} mm ; "
          f"peak -load {top[:,0].min():.0f} N at {top[top[:,0].argmin(),1]:.2f} mm")
    
    #S01.PLOT_2D_FRAME(deformed_scale=1.0)


if __name__ == "__main__":
    main()
