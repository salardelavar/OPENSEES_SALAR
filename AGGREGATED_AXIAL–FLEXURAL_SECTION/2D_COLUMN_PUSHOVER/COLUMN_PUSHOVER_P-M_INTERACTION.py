######################################################################################################################
#                          >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                          #
#                    AXIAL–FLEXURAL SECTION MODELING FOR NONLINEAR COLUMN ANALYSIS USING OPENSEES                    #
#--------------------------------------------------------------------------------------------------------------------#
#                              THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                            #
#                                       EMAIL: salar.d.ghashghaei@gmail.com                                          #
######################################################################################################################
"""
In this example, both the axial behavior (typically elastic)
 and the flexural behavior (moment curvature) are defined
 indepenently and are then "aggregated" into a section.
 This is a characteristic of the uniaxial section: there
 is no coupling of behaviors.
 
1. A nonlinear 2D cantilever column model was developed using an aggregated axial–flexural section to capture uncoupled P–Δ and moment–curvature behavior.
2. Flexural response was modeled via a calibrated **Hysteretic** material to represent cracking, yielding, pinching, and post-yield degradation.
3. Axial behavior was kept elastic to permit pure column-stability interaction during pushover.
4. A nonlinearBeamColumn element with 5 integration points was adopted to ensure accurate curvature distribution.
5. Gravity load was applied incrementally prior to lateral displacement-controlled pushover.
6. The pushover incorporated displacement-control at the top node to capture softening and post-buckling response.
7. Reaction forces were extracted to build axial–displacement, shear–drift, and moment–rotation curves.
8. Results clearly show stiffness degradation, strength loss, and geometric-nonlinearity dominance at large drifts.
9. The axial load–displacement plot indicates the onset of instability and progressive post-buckling softening.
10. Deformed-shape visualization confirms localized curvature concentration and instability forming near the base hinge.
 
OPENSEES EXAMPLE: 
https://opensees.berkeley.edu/OpenSees/manuals/ExamplesManual/HTML/3600.htm

"""
import openseespy.opensees as ops
import time as TI
import os
import numpy as np
import matplotlib.pyplot as plt
import ANALYSIS_FUNCTION as S02
import SECTION_ANALYSIS_FUN as S03
#------------------------------------------------------------------------------------------------
ops.wipe()
if not os.path.exists("Data"):
    os.mkdir("Data")
#------------------------------------------------------------------------------------------------
ops.model('basic', '-ndm', 2, '-ndf', 3)
#------------------------------------------------------------------------------------------------
# Define parameters (units: mm, N)
fy = -240.0           # [N/mm^2] Steel Yield Strength
Es = 21e4             # [N/mm^2] Steel Modulus of Elasticity
LCol = 5000.0         # [mm] Column length
Weight = 200000.0     # [kg] Ttop Column Weight
Dmax = 0.35 * LCol    # [mm] Ultimate Displacement
Dincr = 0.001 * LCol  # [mm] Displacement Increment 
#------------------------------------------------------------------------------------------------
# section geometry
x_c, y_c, A_total, Ix_total, Iy_total = S03.SECTION_ANALYSIS()
#------------------------------------------------------------------------------------------------
# calculated parameters
PCol = Weight
g = 9810.0
Mass = PCol / g
#------------------------------------------------------------------------------------------------
# calculated geometry
ACol = A_total
IzCol = Ix_total
#------------------------------------------------------------------------------------------------
# Define Analysis Properties
MAX_ITERATIONS = 1000      # Convergence iteration for test
MAX_TOLERANCE = 1.0e-6     # Convergence tolerance for test
#------------------------------------------------------------------------------------------------
# nodal coordinates:
N = 1             # Number of elements
node_tag = N + 1  # Top node for displacement control
dy = LCol / N     # Element length    
for i in range(N+1):
    y = i * dy
    x = 0.0
    ops.node(i+1, x, y)
#------------------------------------------------------------------------------------------------
# Define boundary conditions for simply supported column
ops.fix(1, 1, 1, 1)    # Bottom node: fix u_x, u_y, free theta
#ops.fix(N+1, 1, 0, 0)  # Top node: fix u_x, free u_y, free theta
#------------------------------------------------------------------------------------------------
# nodal masses:
ops.mass(N+1, Mass, Mass, 0.0)
#------------------------------------------------------------------------------------------------
ColMatTagFlex = 2
ColMatTagAxial = 3
ColSecTag = 1

EICol = Es * IzCol
EACol = Es * ACol
#------------------------------------------------------------------------------------------------
# Define Nonlinear Sectional Moment-Curvature
matTag_S = 500
FiY = 0.000001474      # [1/mm] Yield Curvature of element
MY = EICol * FiY       # [N.mm] Yield Moment of element
MU = 1.36 * MY         # [N.mm] Utimate Moment of element
FiSU = 0.000541        # [1/mm] Utimate Curvature of element
EIColCrack = MY / FiY  # [N.mm^2] Elasic Flextural Rigidity
#ROTATIONAl_STIFFNESS = 1e14 # [N/mm] Rotational stiffness of springs
pinchX = 0.4           # Pinching factor in X direction
pinchY = 0.2           # Pinching factor in Y direction
damage1 = 0.0          # Damage due to ductility
damage2 = 0.0          # Damage due to energy
beta = 0.1             # Stiffness degradation 

# MOMENT-CURVATURE RELATION FOR STEEL SECTION
ops.uniaxialMaterial('Hysteretic', ColMatTagFlex, MY, FiY, MU, FiSU, 0.23*MU, 1.13*FiSU, -MY, -FiY, -MU, -FiSU, -0.23*MU, -1.07*FiSU, pinchX, pinchY, damage1, damage2, beta)
# INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Hysteretic_Material
#------------------------------------------------------------------------------------------------
ops.uniaxialMaterial("Elastic", ColMatTagAxial, EACol)
#------------------------------------------------------------------------------------------------
ops.section("Aggregator", ColSecTag, ColMatTagAxial, "P", ColMatTagFlex, "Mz")
#------------------------------------------------------------------------------------------------
# Define Geometric Transformation
ColTransfTag = 1
ops.geomTransf("Linear", ColTransfTag)
#ops.geomTransf('Corotational', ColTransfTag)

# Define Element
numIntgrPts = 5
ops.element("nonlinearBeamColumn", 1, 1, 2, numIntgrPts, ColSecTag, ColTransfTag)
#------------------------------------------------------------------------------------------------
ops.recorder("Node", "-file", "Data/DFree.out", "-time", "-node", 2, "-dof", 1, 2, 3, "disp")
ops.recorder("Node", "-file", "Data/DBase.out", "-time", "-node", 1, "-dof", 1, 2, 3, "disp")
ops.recorder("Node", "-file", "Data/RBase.out", "-time", "-node", 1, "-dof", 1, 2, 3, "reaction")
ops.recorder("Element", "-file", "Data/FCol.out", "-time", "-ele", 1, "globalForce")
ops.recorder("Element", "-file", f"Data/ForceColSec1.out", "-time", "-ele", 1, "section", 1, "force")
ops.recorder("Element", "-file", f"Data/DefoColSec1.out", "-time", "-ele", 1, "section", 1, "deformation")
ops.recorder("Element", "-file", f"Data/ForceColSec{numIntgrPts}.out", "-time", "-ele", 1, "section", numIntgrPts, "force")
ops.recorder("Element", "-file", f"Data/DefoColSec{numIntgrPts}.out", "-time", "-ele", 1, "section", numIntgrPts, "deformation")

#------------------------------------------------------------------------------------------------ops.timeSeries('Linear', 1)
# Define Gravity
ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)
ops.load(N+1, 0.0, -PCol, 0.0)

ops.constraints("Plain")
ops.numberer("Plain")
ops.system("BandGeneral")
ops.test("NormDispIncr", MAX_TOLERANCE, MAX_ITERATIONS)
ops.algorithm("Newton")
NstepGravity = 10
ops.integrator("LoadControl", 1.0/NstepGravity)
ops.analysis("Static")
ops.analyze(NstepGravity)

ops.loadConst("-time", 0.0)
print("Model Built")
#------------------------------------------------------------------------------------------------
# Run Static Pushover
IDctrlNode = node_tag
IDctrlDOF = 1

ops.timeSeries('Linear', 100)
ops.pattern('Plain', 100, 1)
ops.load(N+1, 1.0, -1.0, 0.0)

ops.constraints("Plain")
ops.numberer("Plain")
ops.system("BandGeneral")

testType = "EnergyIncr"
ops.test(testType, MAX_TOLERANCE, MAX_ITERATIONS, 2)

ops.algorithm("Newton")
ops.integrator("DisplacementControl", IDctrlNode, IDctrlDOF, Dincr)
ops.analysis("Static")
#------------------------------------------------------------------------------------------------
# Run analysis and record results
mid_node = (N//2) + 1  # Mid-height node for lateral displacement
axial_disp = []        # Axial displacement at top node
rotation = []          # Rotation of column
lateral_disp = []      # Lateral displacement at mid-height
axial_load = []        # Axial compressive load (reaction force)
shear_load = []        # Shear load (reaction force)
moment_load = []       # Moment load (reaction force)

Nsteps = int(Dmax / Dincr) # Number of analysis steps

# Analysis Durations:
starttime = TI.process_time()

for step in range(Nsteps):
    OK = ops.analyze(1)  # Perform one analysis step
    S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS)  # CHECK THE ANALYSIS
    axial_disp.append(ops.nodeDisp(node_tag, 2))         # Record axial displacement
    rotation.append(ops.nodeDisp(node_tag, 3))           # Record rotation
    lateral_disp.append(ops.nodeDisp(node_tag, 1))       # Record lateral displacement
    axial_load.append(-ops.eleResponse(1, 'force')[1])   # Record reaction force at bottom node (Axial compressive load)
    shear_load.append(-ops.eleResponse(1, 'force')[0])   # Record reaction force at bottom node (shear load)
    moment_load.append(-ops.eleResponse(1, 'force')[2])  # Record reaction force at bottom node (moment load)
    print(f'STEP {step+1} DONE')
else:
    print('Analysis completed successfully')
    
totaltime = TI.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n')

#------------------------------------------------------------------------------------------------    
# Plot axial load vs. axial displacement
plt.figure(1, figsize=(8, 6))
plt.plot(axial_disp, np.abs(axial_load), 'black', linewidth=4)
# ADD CRITICAL BUCKLING LOAD AS A HORIZONTAL LINE
plt.xlabel('Axial displacement at top-height (mm)')
plt.ylabel('Axial compressive force (N)')
plt.title('Axial force-displacement behavior of column')
plt.grid(True)
plt.show()    
#------------------------------------------------------------------------------------------------
# Plot axial load vs. lateral displacement
plt.figure(2, figsize=(8, 6))
plt.plot(lateral_disp, np.abs(shear_load), 'black', linewidth=4)
plt.xlabel('Lateral displacement at top-height (mm)')
plt.ylabel('Shear force (N)')
plt.title('Shear force-displacement behavior of column')
plt.grid(True)
plt.show()
#------------------------------------------------------------------------------------------------
# Plot moment vs, rotation
plt.figure(3, figsize=(8, 6))
plt.plot(np.abs(rotation), np.abs(moment_load), 'black', linewidth=4)
plt.xlabel('Rotation at top-height (rad)')
plt.ylabel('Moment (N.mm)')
plt.title('Post-buckling behavior of column')
plt.grid(True)
plt.show()
#------------------------------------------------------------------------------------------------
# Define a function to plot the frame shapes
def PLOT_FRAME(deformed_scale=1.0):
    fig, ax = plt.subplots(1, figsize=(20, 16))

    # Extract node coordinates
    nodes = ops.getNodeTags()
    node_coords = {node: ops.nodeCoord(node) for node in nodes}

    # Plot undeformed shape
    for ele in ops.getEleTags():
        node1, node2 = ops.eleNodes(ele)
        x1, y1 = node_coords[node1]
        x2, y2 = node_coords[node2]
        ax.plot([x1, x2], [y1, y2], 'k-', label='Undeformed' if ele == 1 else "")  # Black line for undeformed

    # Plot deformed shape
    for ele in ops.getEleTags():
        node1, node2 = ops.eleNodes(ele)
        x1, y1 = node_coords[node1]
        x2, y2 = node_coords[node2]

        ux1, uy1, _ = ops.nodeDisp(node1)  # Displacement at node1
        ux2, uy2, _ = ops.nodeDisp(node2)  # Displacement at node2

        ax.plot([x1 + deformed_scale * ux1, x2 + deformed_scale * ux2],
                [y1 + deformed_scale * uy1, y2 + deformed_scale * uy2],
                'r--', label='Deformed' if ele == 1 else "")  # Red dashed line for deformed
                
    # Annotate nodes with their tags
    for node, (x, y) in node_coords.items():
        ux, uy, _ = ops.nodeDisp(node)  # Displacement at node
        ax.text(x, y, f"{node}", color='blue', fontsize=12, ha='center', label='Node Tags' if node == 1 else "")  # Undeformed
        ax.text(x + deformed_scale * ux, y + deformed_scale * uy, f"{node}", color='purple', fontsize=12, ha='center')  # Deformed            

    #ax.set_aspect('equal', 'box')
    ax.set_xlabel('X [mm]')
    ax.set_ylabel('Y [mm]')
    ax.set_title(f'Undeformed and Deformed Shapes - SCALE: {deformed_scale:.3f}')
    ax.legend()
    ax.grid()
    plt.show()
    
# Plot frame shapes
PLOT_FRAME(deformed_scale=1)  # Adjust scale factor as needed
#------------------------------------------------------------------------------------------------
