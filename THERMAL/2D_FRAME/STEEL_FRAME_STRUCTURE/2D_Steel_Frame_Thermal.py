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

import time as TI
import numpy as np
import matplotlib.pyplot as plt
import openseespy.opensees as ops
from Analysis_Function import ANALYSIS

"""
Models and analyzes a 2D steel frame subjected to thermal and distributed loads using OpenSees. 
Key points:  
1. Model Definition: The 2D frame has specified node coordinates for stories and bays, with fixed supports at the base.
 Material properties for steel (with thermal effects) and section geometries (I-section, circular, and rectangular tubes) are defined using fiber elements.  
2. Element and Load Setup: Beam-column elements with corotational geometric transformation and Lobatto beam integration are created.
 Distributed loads are applied to beams, and a thermal gradient is applied to the first story beams.  
3. Analysis Setup: The analysis uses static load control with thermal increments, and the Newton-Raphson algorithm ensures convergence.
 Convergence tolerances and maximum iterations are defined.  
4. Output and Post-processing: Displacements, reactions, and deformations are recorded during the analysis. Data is extracted from output
 files for plotting base reactions (axial, shear, moment) and node displacements against temperature or applied load.  
5. Visualization: The frame's undeformed and deformed shapes are plotted, and results like temperature-displacement relationships
 and base reactions are visualized using Matplotlib.
"""

# Analysis Durations:
starttime = TI.process_time()
#--------------------------------------------------------------------
# Define parameters (units: mm, kN)
num_stories = 3        # Number of Stories
num_bays = 4           # Number of Bays
story_height = 3000.0  # [mm] Height of each story
bay_width = 7000.0     # [mm] Width of each bay
#--------------------------------------------------------------------
# DefineSteel Section Material Properties (Steel01Thermal)
fy = 0.24              # [kN/mm^2] Yield strength of steel section
fu = 1.5 * fy          # [kN/mm^2] Ultimate strength of steel section
Es = 200               # [kN/mm^2] Modulus of elasticity of steel section
ey = fy / Es           # [mm/mm] Yield steel strain
esu = 0.35             # [mm/mm] Ultimate steel strain
Esh = (fu - fy) / (esu - ey)  # [kN/mm^2] Strain hardening modulus
b = Esh / Es                  # Strain hardening ratio
#--------------------------------------------------------------------
# Define Thermal and Distributed load
Max_Thermal = 300.0       # [°C] Temperature
distributed_load = -0.3   # [kN/mm] Distributed load
#--------------------------------------------------------------------
# Define Analysis Properties
MAX_ITERATIONS = 10000    # Convergence iteration for test
TOLERANCE = 1.0e-10       # Convergence tolerance for test

Incr_Temp = 0.0005                 # Incremental temperature step
#Nstep = int(Thermal/Incr_Temp)    # Number of incremental steps
Nstep = 600                        # Number of incremental steps
#--------------------------------------------------------------------
# Define model
ops.model('basic', '-ndm', 2, '-ndf', 3)
ops.wipe()
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
#----------------------------
# I SECTION
#----------------------------
def I_SECTION(secTag):
    # Define geometric properties of the steel tube
    #secTag = 1
    # Define geometric properties of the steel I section
    ops.section('FiberThermal', secTag)
    # Define section (FiberThermal)
    bf = 300               # [mm] Flange width
    tf = 20                # [mm] Flange thickness
    tw = 10                # [mm] Web thickness
    hw = 400               # [mm] Web height
    NUM = 100              # Number of fibers for web
    d = 2 * tf + hw
    # Top flange fibers
    ops.fiber(bf / 2, hw / 2 + tf / 2, bf * tf / 4, matTag)
    ops.fiber(-bf / 2, hw / 2 + tf / 2, bf * tf / 4, matTag)
    # Bottom flange fibers
    ops.fiber(bf / 2, -hw / 2 - tf / 2, bf * tf / 4, matTag)
    ops.fiber(-bf / 2, -hw / 2 - tf / 2, bf * tf / 4, matTag)
    # Web fibers
    for i in range(NUM):
        yLoc = hw / 2 - i * (hw / NUM)
        ops.fiber(0.0, yLoc, tw * hw / NUM, matTag) 
    
    #-------------------
    # PLOT THE SECTION
    #-------------------
    import matplotlib.pyplot as plt
    import numpy as np
    fig, ax = plt.subplots(figsize=(8, 8))

    # Plot the top flange fibers
    flange_area = bf * tf / 4
    ax.add_patch(plt.Rectangle((-bf / 2, hw / 2), bf, tf, color='lightgrey', alpha=0.7))

    # Plot the bottom flange fibers
    ax.add_patch(plt.Rectangle((-bf / 2, -hw / 2 - tf), bf, tf, color='lightgrey', alpha=0.7))

    # Plot the web fibers
    fiber_height = hw / NUM
    for i in range(NUM):
        yLoc = hw / 2 - i * fiber_height
        ax.add_patch(plt.Rectangle((-tw / 2, yLoc - fiber_height / 2), tw, fiber_height, color='lightgrey', alpha=0.7))

    # Add section outline
    #ax.add_patch(plt.Rectangle((-bf / 2, -hw / 2 - tf), bf, hw + 2 * tf, edgecolor='black', fill=False, linewidth=1.5))

    # Set plot limits
    ax.set_xlim([-bf / 2 - 20, bf / 2 + 20])
    ax.set_ylim([-hw / 2 - tf - 20, hw / 2 + tf + 20])
    ax.set_aspect('equal', adjustable='box')

    # Labels and title
    plt.xlabel('Width (bf)', fontsize=12)
    plt.ylabel('Height (hw + 2 * tf)', fontsize=12)
    plt.title('I-Section with Fibers', fontsize=14)

    # Show plot
    plt.grid(True)
    plt.show()
    
    return d # Return Section Height
    
#----------------------------    
# CIRCULAR TUBE SECTION  
#---------------------------- 
def C_TUBE_SECTION(secTag):
    # Define geometric properties of the steel tube
    #secTag = 1
    D = 500          # [mm] Outer diameter of the tube
    t = 20           # [mm] Wall thickness of the tube
    r = D / 2        # [mm] Outer radius of the tube (m)
    r_inner = r - t  # [mm] Inner radius of the tube

    # Define material tag
    matTag = 1

    # Number of fibers along the radial and circumferential directions
    NUM_R = 10  # Number of fibers radially
    NUM_C = 20  # Number of fibers circumferentially

    # Define the steel tube section using fibers
    ops.section('FiberThermal', secTag)

    # Loop over radial and circumferential directions to define fibers
    import numpy as np
    for i in range(NUM_R):
        # Compute the radial location of the fiber
        r_loc = r_inner + (i + 0.5) * (t / NUM_R)
        for j in range(NUM_C):
            # Compute the angular location of the fiber
            theta = j * (2 * np.pi / NUM_C)
            x_loc = r_loc * np.cos(theta)
            y_loc = r_loc * np.sin(theta)
            # Compute the area of the fiber
            fiber_area = (2 * np.pi * r_loc / NUM_C) * (t / NUM_R)
            # Add the fiber to the section
            ops.fiber(x_loc, y_loc, fiber_area, matTag)
    
    #-------------------
    # PLOT THE SECTION
    #-------------------
    import matplotlib.pyplot as plt
    # Calculate radii
    r = D / 2  # Outer radius
    r_inner = r - t  # Inner radius

    # Initialize the plot
    fig, ax = plt.subplots(figsize=(8, 8))

    # Plot fibers
    for i in range(NUM_R):
        # Radial location of fibers
        r_loc = r_inner + (i + 0.5) * (t / NUM_R)
        for j in range(NUM_C):
            # Angular location of fibers
            theta = j * (2 * np.pi / NUM_C)
            x_loc = r_loc * np.cos(theta)
            y_loc = r_loc * np.sin(theta)
            ax.add_patch(plt.Circle((x_loc, y_loc), 2, color='lightgrey', alpha=0.7))  # Fibers as small circles

    # Add the outer and inner circles for reference
    outer_circle = plt.Circle((0, 0), r, color='black', fill=False, linewidth=1.5)
    inner_circle = plt.Circle((0, 0), r_inner, color='black', fill=False, linewidth=1.5)
    ax.add_patch(outer_circle)
    ax.add_patch(inner_circle)

    # Set aspect ratio and limits
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlim([-r - 20, r + 20])
    ax.set_ylim([-r - 20, r + 20])

    # Labels and title
    plt.xlabel('X (mm)', fontsize=12)
    plt.ylabel('Y (mm)', fontsize=12)
    plt.title('Circular Tube Section with Fibers', fontsize=14)

    # Show the plot
    plt.grid(True)
    plt.show()
    
    return D # Return Section Height

#----------------------------
# RECTANGULAR TUBE SECTION 
#----------------------------
def R_TUBE_SECTION(secTag):
    # Define section tag for rectangular tube steel section
    #secTag = 1
    # Define geometric properties of the rectangular tube
    B = 300  # [mm] Outer width of the tube
    H = 500  # [mm] Outer height of the tube
    t = 20   # [mm] Wall thickness of the tube
    # Define material tag
    matTag = 1
    # Number of fibers along each wall direction
    NUM_B = 10  # Number of fibers along the width
    NUM_H = 20  # Number of fibers along the height
    # Define the rectangular tube section using fibers
    ops.section('FiberThermal', secTag)
    # Outer top wall fibers
    for i in range(NUM_B):
        x_loc = -B / 2 + (i + 0.5) * (B / NUM_B)
        ops.fiber(x_loc, H / 2 - t / 2, (B / NUM_B) * t, matTag)
    # Outer bottom wall fibers
    for i in range(NUM_B):
        x_loc = -B / 2 + (i + 0.5) * (B / NUM_B)
        ops.fiber(x_loc, -H / 2 + t / 2, (B / NUM_B) * t, matTag)
    # Outer left wall fibers
    for i in range(NUM_H):
        y_loc = -H / 2 + (i + 0.5) * (H / NUM_H)
        ops.fiber(-B / 2 + t / 2, y_loc, t * (H / NUM_H), matTag)
    # Outer right wall fibers
    for i in range(NUM_H):
        y_loc = -H / 2 + (i + 0.5) * (H / NUM_H)
        ops.fiber(B / 2 - t / 2, y_loc, t * (H / NUM_H), matTag)
        
    
    #-------------------
    # PLOT THE SECTION
    #-------------------
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 8))

    # Plot outer top wall fibers
    fiber_width = B / NUM_B
    for i in range(NUM_B):
        x_loc = -B / 2 + i * fiber_width
        ax.add_patch(plt.Rectangle((x_loc, H / 2 - t), fiber_width, t, color='lightgrey', alpha=0.7))

    # Plot outer bottom wall fibers
    for i in range(NUM_B):
        x_loc = -B / 2 + i * fiber_width
        ax.add_patch(plt.Rectangle((x_loc, -H / 2), fiber_width, t, color='lightgrey', alpha=0.7))

    # Plot outer left wall fibers
    fiber_height = H / NUM_H
    for i in range(NUM_H):
        y_loc = -H / 2 + i * fiber_height
        ax.add_patch(plt.Rectangle((-B / 2, y_loc), t, fiber_height, color='lightgrey', alpha=0.7))

    # Plot outer right wall fibers
    for i in range(NUM_H):
        y_loc = -H / 2 + i * fiber_height
        ax.add_patch(plt.Rectangle((B / 2 - t, y_loc), t, fiber_height, color='lightgrey', alpha=0.7))

    # Add section outline
    outer_outline = plt.Rectangle((-B / 2, -H / 2), B, H, edgecolor='black', fill=False, linewidth=1.5)
    inner_outline = plt.Rectangle((-B / 2 + t, -H / 2 + t), B - 2 * t, H - 2 * t, edgecolor='black', fill=False, linewidth=1.5)
    ax.add_patch(outer_outline)
    ax.add_patch(inner_outline)

    # Set plot limits
    ax.set_xlim([-B / 2 - 20, B / 2 + 20])
    ax.set_ylim([-H / 2 - 20, H / 2 + 20])
    ax.set_aspect('equal', adjustable='box')

    # Labels and title
    plt.xlabel('Width (B)', fontsize=12)
    plt.ylabel('Height (H)', fontsize=12)
    plt.title('Rectangular Tube Section with Fibers', fontsize=14)

    # Show plot
    plt.grid(True)
    plt.show()
    
    return H # Return Section Height
#---------------------------
# RECTANGULAR SECTION 
#---------------------------
def R_RECTANGULAR_STEEL_SECTION(secTag):
    # Define a rectangular steel section using fibers.
    B = 400       # [mm] Width of the rectangular section
    H = 500       # [mm] Height of the rectangular section
    NUM_B = 10   # Number of fibers along the width of the section
    NUM_H = 10   # Number of fibers along the height of the section
    secTag = 1    # Section tag identifier
    matTag = 1    # Material tag for the concrete
    ops.section('FiberThermal', secTag)
    # Create fibers for the entire section
    for i in range(NUM_B):
        for j in range(NUM_H):
            # Calculate the x and y location of the fiber centroid
            x_loc = -B / 2 + (i + 0.5) * (B / NUM_B)
            y_loc = -H / 2 + (j + 0.5) * (H / NUM_H)
            # Define the fiber area
            fiber_area = (B / NUM_B) * (H / NUM_H)
            # Add the fiber to the section
            ops.fiber(x_loc, y_loc, fiber_area, matTag)
            
    #-------------------
    # PLOT THE SECTION
    #------------------- 
    fiber_width = B / NUM_B
    fiber_height = H / NUM_H
    
    # Prepare the plot
    fig, ax = plt.subplots()
    ax.set_aspect('equal', 'box')
    ax.set_xlim(-B / 2, B / 2)
    ax.set_ylim(-H / 2, H / 2)
    
    # Plot the fibers
    for i in range(NUM_B):
        for j in range(NUM_H):
            # Calculate the x and y location of the fiber centroid
            x_loc = -B / 2 + i * fiber_width
            y_loc = -H / 2 + j * fiber_height
            # Draw the rectangle representing the fiber
            rect = plt.Rectangle((x_loc, y_loc), fiber_width, fiber_height, color='lightgrey', edgecolor='black', alpha=0.7)
            ax.add_patch(rect)
    
    # Add labels and grid
    ax.set_xlabel("Width [mm]")
    ax.set_ylabel("Height [mm]")
    ax.grid(True)
    ax.set_title("Rectangular Steel Section with Fibers")
    
    # Show the plot
    plt.show()           

    return H # Return Section Height
    
#--------------------------------------------------------------------
secTag = 1
Depth = I_SECTION(secTag) 
#Depth = C_TUBE_SECTION(secTag)
#Depth = R_TUBE_SECTION(secTag)
#Depth = R_RECTANGULAR_STEEL_SECTION(secTag)

print(' Thermal Section Done.')
#--------------------------------------------------------------------
# Define geometric transformation
transfTag = 1
ops.geomTransf('Corotational', transfTag)
#--------------------------------------------------------------------
# Define beam integration
NP = 3
biTag = 1
ops.beamIntegration('Lobatto', biTag, secTag, NP)
#--------------------------------------------------------------------
# Define beam-column elements
# Thermo-mechaical beam-column Elements
# LINK: https://openseesforfire.github.io/Subpages/Elecmds.html
element_id = 1
beam_id = [] # Initialize the beam ID list
for i in range(num_stories):
    for j in range(num_bays + 1):  # Vertical elements (columns)
        node1 = i * (num_bays + 1) + j + 1
        node2 = (i + 1) * (num_bays + 1) + j + 1
        # element dispBeamColumnThermal $eleID $NodeTag0 $NodeTag1 $NumInts $secTag $TransfTag
        ops.element('dispBeamColumnThermal', element_id, node1, node2, transfTag, biTag)
        element_id += 1

    for j in range(num_bays):  # Horizontal elements (beams)
        node1 = (i + 1) * (num_bays + 1) + j + 1
        node2 = (i + 1) * (num_bays + 1) + j + 2
        ops.element('dispBeamColumnThermal', element_id, node1, node2, transfTag, biTag)
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
# A linear thermal gradient is applied to the elements in the first story.
# The temperature varies from -DD to DD across the section height.
#LINK: https://openseesforfire.github.io/Subpages/ThermalActionCmds.html

# Apply thermal load to all elements
#for ele_id in range(1, element_id):
#    ops.eleLoad('-ele', ele_id, '-type', '-beamThermal', Max_Thermal, -DD, Max_Thermal, DD)

# Apply thermal load only to the identified beam elements
DD = 0.5 * Depth
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
ops.test('NormDispIncr', TOLERANCE, MAX_ITERATIONS, 1)
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
#--------------------------------------------------------------------
# Perform analysis
temp = []
dispX = []
dispY = []

mid_node = (num_stories * (num_bays + 1) + num_bays // 2 + 1)  # Middle node at top

for i in range(Nstep):
    print('STEP: ', i + 1)
    OK = ops.analyze(1)
    ANALYSIS(OK, 1, TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS

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
plot_frame(deformed_scale=60)  # Adjust scale factor as needed

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
base01_X = OUTPUT_SECOND_COLUMN('BTH_PUSH_01', 1) # Reading base reaction from Text file - NODE 1
base02_X = OUTPUT_SECOND_COLUMN('BTH_PUSH_02', 1) # Reading base reaction from Text file - NODE 2
base03_X = OUTPUT_SECOND_COLUMN('BTH_PUSH_03', 1) # Reading base reaction from Text file - NODE 3
base04_X = OUTPUT_SECOND_COLUMN('BTH_PUSH_04', 1) # Reading base reaction from Text file - NODE 4
base05_X = OUTPUT_SECOND_COLUMN('BTH_PUSH_05', 1) # Reading base reaction from Text file - NODE 5
BASES_AXIAL = base01_X + base02_X + base03_X + base04_X + base05_X
# SHEAR
base01_Y = OUTPUT_SECOND_COLUMN('BTH_PUSH_01', 2) # Reading base reaction from Text file - NODE 1
base02_Y = OUTPUT_SECOND_COLUMN('BTH_PUSH_02', 2) # Reading base reaction from Text file - NODE 2
base03_Y = OUTPUT_SECOND_COLUMN('BTH_PUSH_03', 2) # Reading base reaction from Text file - NODE 3
base04_Y = OUTPUT_SECOND_COLUMN('BTH_PUSH_04', 2) # Reading base reaction from Text file - NODE 4
base05_Y = OUTPUT_SECOND_COLUMN('BTH_PUSH_05', 2) # Reading base reaction from Text file - NODE 5
BASES_SHEAR = base01_Y + base02_Y + base03_Y + base04_Y + base05_Y
# MOMENT
base01_Z = OUTPUT_SECOND_COLUMN('BTH_PUSH_01', 3) # Reading base reaction from Text file - NODE 1
base02_Z = OUTPUT_SECOND_COLUMN('BTH_PUSH_02', 3) # Reading base reaction from Text file - NODE 2
base03_Z = OUTPUT_SECOND_COLUMN('BTH_PUSH_03', 3) # Reading base reaction from Text file - NODE 3
base04_Z = OUTPUT_SECOND_COLUMN('BTH_PUSH_04', 3) # Reading base reaction from Text file - NODE 4
base05_Z = OUTPUT_SECOND_COLUMN('BTH_PUSH_05', 3) # Reading base reaction from Text file - NODE 5
BASES_MOMENT = base01_Z + base02_Z + base03_Z + base04_Z + base05_Z

ops.wipe()
#--------------------------------------------------------------------

# Plot results
plt.figure(figsize=(8, 6))
plt.plot(disp_Y, BASES_AXIAL, color='blue', linewidth=2)
plt.xlabel('Node 6 Displacement Y (mm)')
plt.ylabel('Base Axial Reaction (kN)')
plt.title('Base Axial Reaction vs Displacement Y')
plt.grid()
plt.show()

plt.figure(figsize=(8, 6))
plt.plot(disp_X, BASES_SHEAR, color='red', linewidth=2)
plt.xlabel('Node 6 Displacement X (mm)')
plt.ylabel('Base Shear Reaction (kN)')
plt.title('Base Shear Reaction vs Displacement X')
plt.grid()
plt.show()

plt.figure(figsize=(8, 6))
plt.plot(disp_Z, BASES_MOMENT, color='green', linewidth=2)
plt.xlabel('Node 6 Rotation (rad)')
plt.ylabel('Base Moment Reaction (kN.mm)')
plt.title('Base Moment Reaction vs Rotation')
plt.grid()
plt.show()

plt.figure(figsize=(8, 6))
plt.plot(temp, dispX, color='purple', linewidth=2)
plt.xlabel('Temperature (°C)')
plt.ylabel(f'Middle Node {mid_node} Displacement X (mm)')
plt.title('Temperature vs Displacement X')
plt.grid()
plt.show()

plt.figure(figsize=(8, 6))
plt.plot(temp, dispY, color='orange', linewidth=2)
plt.xlabel('Temperature (°C)')
plt.ylabel(f'Middle Node {mid_node} Displacement Y (mm)')
plt.title('Temperature vs Displacement Y')
plt.grid()
plt.show()

#--------------------------------------------------------------------