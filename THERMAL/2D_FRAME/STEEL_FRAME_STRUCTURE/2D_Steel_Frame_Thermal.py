####################################################################################
#                              IN THE NAME OF ALLAH                                #
#     THERMAL LOADING STRUCTURAL ANALYSIS OF A 2D STEEL FRAME USING OPENSEES       #
#----------------------------------------------------------------------------------#
# In all the beams of the floors, the extended load due to the dead load has been  #
# applied, and only in the beams of the first floor, the load due to the heat of   #
# the fire has been applied.                                                       #
#----------------------------------------------------------------------------------#
#       THIS PROGRAM IS WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)              #
#                    EMAIL: SALAR.D.GHASHGHAEI@GMAIL.COM                           #
####################################################################################

"""
Models and analyzes a 2D steel frame subjected to thermal and distributed loads using OpenSees. 
Key points:  
[1] Model Definition: The 2D frame has specified node coordinates for stories and bays, with fixed supports at the base.
 Material properties for steel (with thermal effects) and section geometries (I-section, circular, and rectangular tubes) are defined using fiber elements.  
[2] Element and Load Setup: Beam-column elements with corotational geometric transformation and Lobatto beam integration are created.
 Distributed loads are applied to beams, and a thermal gradient is applied to the first story beams.  
[3] Analysis Setup: The analysis uses static load control with thermal increments, and the Newton-Raphson algorithm ensures convergence.
 Convergence tolerances and maximum iterations are defined.  
[4] Output and Post-processing: Displacements, reactions, and deformations are recorded during the analysis. Data is extracted from output
 files for plotting base reactions (axial, shear, moment) and node displacements against temperature or applied load.  
[5] Visualization: The frame's undeformed and deformed shapes are plotted, and results like temperature-displacement relationships
 and base reactions are visualized using Matplotlib.
"""

import time as TI
import numpy as np
import matplotlib.pyplot as plt
import openseespy.opensees as ops
from Analysis_Function import ANALYSIS
from STEEL_FIBERTHERMAL_SECTION import I_SECTION, DOUBLE_I_SECTION, C_TUBE_SECTION, R_TUBE_SECTION, R_RECTANGULAR_STEEL_SECTION

#--------------------------------------------------------------------
# Define parameters (units: mm, kN)
# Define structure parameters
num_stories = 3        # Number of Stories
num_bays = 4           # Number of Bays
story_height = 3000.0  # [mm] Height of each story
bay_width = 7000.0     # [mm] Width of each bay
#--------------------------------------------------------------------
# DefineSteel Section Material Properties (Steel01Thermal)
fy = 0.24              # [kN/mm²] Yield strength of steel section
fu = 1.5 * fy          # [kN/mm²] Ultimate strength of steel section
Es = 200               # [kN/mm²] Modulus of elasticity of steel section
ey = fy / Es           # [mm/mm] Yield steel strain
esu = 0.35             # [mm/mm] Ultimate steel strain
Esh = (fu - fy) / (esu - ey)  # [kN/mm^2] Strain hardening modulus
b = Esh / Es                  # Strain hardening ratio
#--------------------------------------------------------------------
# Define Thermal and Distributed load
Max_Thermal = 800.0        # [°C] Temperature
distributed_load = -0.003  # [kN/mm] Uniform Distributed Loads
#--------------------------------------------------------------------
# Define Analysis Properties
MAX_ITERATIONS = 10000    # Convergence iteration for test
MAX_TOLERANCE = 1.0e-10   # Convergence tolerance for test

Nstep = 400           # Number of incremental steps
Incr_Temp = 1/Nstep   # Incremental temperature step
#--------------------------------------------------------------------
# Clear pervious model
ops.wipe()
# Define model
ops.model('basic', '-ndm', 2, '-ndf', 3)
# Node coordinates
node_id = 1
for i in range(num_stories + 1):  # Including ground level
    for j in range(num_bays + 1):  # Including leftmost column
        ops.node(node_id, j * bay_width, i * story_height)
        if i == 0:  # Fix supports at ground level
            ops.fix(node_id, 1, 1, 1)
        else:  # Free to move for higher levels
            ops.fix(node_id, 0, 0, 0)
        node_id += 1
        
print(' Structure Coordinates Done.') 
#--------------------------------------------------------------------
# Materials (MATERIAL NONLINEARITY)
matTag = 1
ops.uniaxialMaterial('Steel01Thermal', matTag, fy, Es, b) # Steel with bilinear kinematic hardening Material with thermaal effect
#--------------------------------------------------------------------
# BEAMS SECTION
secTag01 = 1
Depth01 = I_SECTION(secTag01, matTag, PLOT=True) 
# COLUMNS SECTION
secTag02 = 2
Depth02 = DOUBLE_I_SECTION(secTag02, matTag, PLOT=True)
#Depth = C_TUBE_SECTION(secTag, matTag, PLOT=True)
#Depth = R_TUBE_SECTION(secTag, matTag, PLOT=True)
#Depth = R_RECTANGULAR_STEEL_SECTION(secTag, matTag, PLOT=True)

print(' Thermal Section Done.')
#--------------------------------------------------------------------
# Define geometric transformation
transfTag = 1
ops.geomTransf('Corotational', transfTag)
#--------------------------------------------------------------------
# Define beam integration (Lobatto integration)
numIntegrationPoints = 5
biTag01 = 1
ops.beamIntegration('Lobatto', biTag01, secTag01, numIntegrationPoints)

# Define column integration (Lobatto integration)
numIntegrationPoints = 5
biTag02 = 2
ops.beamIntegration('Lobatto', biTag02, secTag02, numIntegrationPoints)
#--------------------------------------------------------------------
# Define beam-column elements
# Thermo-mechaical beam-column Elements
# INFO LINK: https://openseesforfire.github.io/Subpages/Elecmds.html
element_id = 1
beam_id = [] # Initialize the beam ID list
for i in range(num_stories):
    for j in range(num_bays + 1):  # Vertical elements (columns)
        node1 = i * (num_bays + 1) + j + 1
        node2 = (i + 1) * (num_bays + 1) + j + 1
        # element dispBeamColumnThermal $eleID $NodeTag0 $NodeTag1 $NumInts $secTag $TransfTag
        ops.element('dispBeamColumnThermal', element_id, node1, node2, transfTag, biTag02)
        # element forceBeamColumnThermal $eleID $NodeTag0 $NodeTag1 $NumInts $secTag $TransfTag
        #ops.element('forceBeamColumnThermal', element_id, node1, node2, secTag, transfTag)
        element_id += 1

    for j in range(num_bays):  # Horizontal elements (beams)
        node1 = (i + 1) * (num_bays + 1) + j + 1
        node2 = (i + 1) * (num_bays + 1) + j + 2
        ops.element('dispBeamColumnThermal', element_id, node1, node2, transfTag, biTag01)
        if i == 0:  # Only store beams from the first story
            beam_id.append(element_id)
        element_id += 1
        
print(' Thermal Element Done.')        
#--------------------------------------------------------------------
# Define time series
Tag = 1
ops.timeSeries('Linear', Tag)
# Define load pattern
patternTag = 1
ops.pattern('Plain', patternTag, Tag)
#--------------------------------------------------------------------
# Apply distributed load to all beams in the structure
for ele_id in beam_id:
    # eleLoad -ele $eleTag1 <$eleTag2 ....> -type -beamUniform $Wz <$Wx>
    ops.eleLoad('-ele', ele_id, '-type', '-beamUniform', distributed_load)
#--------------------------------------------------------------------
"""
# Define uniform distributed load analysis
ops.system('BandGeneral')
ops.constraints('Plain')
ops.numberer('Plain')
ops.test('NormDispIncr', TOLERANCE, MAX_ITERATIONS) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/normDispIncr.html
ops.algorithm('Newton')
ops.integrator('LoadControl', 0.1)
ops.analysis('Static')
ops.analyze(10)

print(' Uniform Distributed Load Analysis Properties Done.')   
ops.wipeAnalysis()   
"""  
#--------------------------------------------------------------------    
# A linear thermal gradient is applied to the elements in the first story.
# The temperature varies from -DD to DD across the section height.
# INFO LINK: https://openseesforfire.github.io/Subpages/ThermalActionCmds.html

# Apply thermal load to all elements
#for ele_id in range(1, element_id):
#    ops.eleLoad('-ele', ele_id, '-type', '-beamThermal', Max_Thermal, -DD, Max_Thermal, DD)

# Apply thermal load only to the identified beam elements
DD = 0.5 * Depth01
for ele_id in beam_id:
    if ele_id == 6 or ele_id == 7 or ele_id == 8 or ele_id == 9: # Apply Thermal Loads only to the Beams of Story 1
        # eleLoad -ele $eleTag -type -beamThermal $T1 $y1 $T2 $Y2
        ops.eleLoad('-ele', ele_id, '-type', '-beamThermal', Max_Thermal, -DD, Max_Thermal, DD)  

print(' Thermal Force Done.')      
#--------------------------------------------------------------------
# Define analysis
ops.system('BandGeneral')
ops.constraints('Plain')
ops.numberer('Plain')
ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS, 2) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/normDispIncr.html
ops.algorithm('Newton')
ops.integrator('LoadControl', Incr_Temp)
ops.analysis('Static')

print(' Thermal Analysis Properties Done.')
#--------------------------------------------------------------------
# Output Data
ops.recorder('Node', '-file', "DTH_PUSH.txt",'-time', '-node', 6, '-dof', 1,2,3, 'disp')        # Displacement Time History Node 6
ops.recorder('Node', '-file', "BTH_PUSH_01.txt",'-time', '-node', 1, '-dof', 1,2,3, 'reaction') # Base Reaction Time History Node 1
ops.recorder('Node', '-file', "BTH_PUSH_02.txt",'-time', '-node', 2, '-dof', 1,2,3, 'reaction') # Base Reaction Time History Node 2
ops.recorder('Node', '-file', "BTH_PUSH_03.txt",'-time', '-node', 3, '-dof', 1,2,3, 'reaction') # Base Reaction Time History Node 3
ops.recorder('Node', '-file', "BTH_PUSH_04.txt",'-time', '-node', 4, '-dof', 1,2,3, 'reaction') # Base Reaction Time History Node 4
ops.recorder('Node', '-file', "BTH_PUSH_05.txt",'-time', '-node', 5, '-dof', 1,2,3, 'reaction') # Base Reaction Time History Node 5
#ops.recorder('Element', '-file', 'STRESS_STRAIN_BOT.txt', '-time', '-ele', 6, 'section', secTag, 'fiber', -DD, 0, 'stressStrainTangent')
#ops.recorder('Element', '-file', 'STRESS_STRAIN_TOP.txt', '-time', '-ele', 6, 'section', secTag, 'fiber', +DD, 0, 'stressStrainTangent')
ops.recorder('Element', '-file', 'STRESS_STRAIN_BOT.txt', '-time', '-ele', 6, 'section', secTag01, 'fiber', -DD, 0, 'stressStrain')
ops.recorder('Element', '-file', 'STRESS_STRAIN_TOP.txt', '-time', '-ele', 6, 'section', secTag01, 'fiber', +DD, 0, 'stressStrain')
#--------------------------------------------------------------------
# Perform analysis
temp = []
dispX = []
dispY = []

mid_node = (num_stories * (num_bays + 1) + num_bays // 2 + 1)  # Middle node at top

# Analysis Durations:
starttime = TI.process_time()

for i in range(Nstep):
    print('STEP: ', i + 1)
    OK = ops.analyze(1)
    ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS

    temp.append(ops.getLoadFactor(patternTag) * Max_Thermal)
    dispX.append(ops.nodeDisp(mid_node, 1)) # X Displacement
    dispY.append(ops.nodeDisp(mid_node, 2)) # Y Displacement

totaltime = TI.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n')

#--------------------------------------------------------------------

# Define a function to plot the frame shapes
def plot_frame(deformed_scale=1.0):
    fig, ax = plt.subplots(figsize=(20, 16))

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
    ax.set_title('Undeformed and Deformed Shapes')
    ax.legend()
    ax.grid()
    plt.show()
    
# Plot frame shapes
plot_frame(deformed_scale=500)  # Adjust scale factor as needed

#--------------------------------------------------------------------

def OUTPUT_SECOND_COLUMN(X, COLUMN):
    import numpy as np
    # Time History
    filename = f"{X}.txt"
    data_collected = np.loadtxt(filename)
    X = data_collected[:, COLUMN]   
    return X 
    
disp_X = OUTPUT_SECOND_COLUMN('DTH_PUSH', 1) # Reading Disp from Text file - X Direaction
disp_Y = OUTPUT_SECOND_COLUMN('DTH_PUSH', 2) # Reading Disp from Text file - Y Direaction
disp_Z = OUTPUT_SECOND_COLUMN('DTH_PUSH', 3) # Reading Disp from Text file - Z Direaction (Rotation)
# AXIAL
base01_Y = OUTPUT_SECOND_COLUMN('BTH_PUSH_01', 1) # Reading base reaction from Text file - NODE 1
base02_Y = OUTPUT_SECOND_COLUMN('BTH_PUSH_02', 1) # Reading base reaction from Text file - NODE 2
base03_Y = OUTPUT_SECOND_COLUMN('BTH_PUSH_03', 1) # Reading base reaction from Text file - NODE 3
base04_Y = OUTPUT_SECOND_COLUMN('BTH_PUSH_04', 1) # Reading base reaction from Text file - NODE 4
base05_Y = OUTPUT_SECOND_COLUMN('BTH_PUSH_05', 1) # Reading base reaction from Text file - NODE 5
BASES_AXIAL = base01_Y + base02_Y + base03_Y + base04_Y + base05_Y
# SHEAR
base01_X = OUTPUT_SECOND_COLUMN('BTH_PUSH_01', 2) # Reading base reaction from Text file - NODE 1
base02_X = OUTPUT_SECOND_COLUMN('BTH_PUSH_02', 2) # Reading base reaction from Text file - NODE 2
base03_X = OUTPUT_SECOND_COLUMN('BTH_PUSH_03', 2) # Reading base reaction from Text file - NODE 3
base04_X = OUTPUT_SECOND_COLUMN('BTH_PUSH_04', 2) # Reading base reaction from Text file - NODE 4
base05_X = OUTPUT_SECOND_COLUMN('BTH_PUSH_05', 2) # Reading base reaction from Text file - NODE 5
BASES_SHEAR = base01_X + base02_X + base03_X + base04_X + base05_X
# MOMENT
base01_Z = OUTPUT_SECOND_COLUMN('BTH_PUSH_01', 3) # Reading base reaction from Text file - NODE 1
base02_Z = OUTPUT_SECOND_COLUMN('BTH_PUSH_02', 3) # Reading base reaction from Text file - NODE 2
base03_Z = OUTPUT_SECOND_COLUMN('BTH_PUSH_03', 3) # Reading base reaction from Text file - NODE 3
base04_Z = OUTPUT_SECOND_COLUMN('BTH_PUSH_04', 3) # Reading base reaction from Text file - NODE 4
base05_Z = OUTPUT_SECOND_COLUMN('BTH_PUSH_05', 3) # Reading base reaction from Text file - NODE 5
BASES_MOMENT = base01_Z + base02_Z + base03_Z + base04_Z + base05_Z
# STRESS AND STRAIN OF ELEMENT 6
strain_B = OUTPUT_SECOND_COLUMN('STRESS_STRAIN_BOT', 2) # Reading bottom strain from Text file
stress_B = OUTPUT_SECOND_COLUMN('STRESS_STRAIN_BOT', 1) # Reading bottom stress from Text file
strain_T = OUTPUT_SECOND_COLUMN('STRESS_STRAIN_TOP', 2) # Reading top strain from Text file
stress_T = OUTPUT_SECOND_COLUMN('STRESS_STRAIN_TOP', 1) # Reading top stress from Text file

ops.wipe()
#--------------------------------------------------------------------

# Plot results
plt.figure(1, figsize=(8, 6))
plt.plot(disp_Y, BASES_AXIAL, color='blue', linewidth=2)
plt.xlabel('Node 6 Displacement Y (mm)')
plt.ylabel('Base Axial Reaction (kN)')
plt.title('Base Axial Reaction vs Displacement Y')
plt.grid()
plt.show()

plt.figure(2, figsize=(8, 6))
plt.plot(disp_X, BASES_SHEAR, color='red', linewidth=2)
plt.xlabel('Node 6 Displacement X (mm)')
plt.ylabel('Base Shear Reaction (kN)')
plt.title('Base Shear Reaction vs Displacement X')
plt.grid()
plt.show()

plt.figure(3, figsize=(8, 6))
plt.plot(disp_Z, BASES_MOMENT, color='green', linewidth=2)
plt.xlabel('Node 6 Rotation (rad)')
plt.ylabel('Base Moment Reaction (kN.mm)')
plt.title('Base Moment Reaction vs Rotation')
plt.grid()
plt.show()

plt.figure(4, figsize=(8, 6))
plt.plot(temp, dispX, color='purple', linewidth=2)
plt.xlabel('Temperature (°C)')
plt.ylabel(f'Middle Node {mid_node} Displacement X (mm)')
plt.title('Temperature vs Displacement X')
plt.grid()
plt.show()

plt.figure(5, figsize=(8, 6))
plt.plot(temp, dispY, color='black', linewidth=2)
plt.xlabel('Temperature (°C)')
plt.ylabel(f'Middle Node {mid_node} Displacement Y (mm)')
plt.title('Temperature vs Displacement Y')
plt.grid()
plt.show()

plt.figure(6, figsize=(8, 6))
plt.plot(strain_B, stress_B, color='blue', label='Bottom Fiber', linewidth=2)
plt.plot(strain_T, stress_T, color='red', label='Top Fiber', linewidth=2)
plt.xlabel('Strain (mm/mm)')
plt.ylabel('Stress (kN/mm^2)')
plt.title(f'Stress-Strain Relation of Element {6} Top & Bottom Fibers')
plt.grid()
plt.legend()
plt.show()

#--------------------------------------------------------------------