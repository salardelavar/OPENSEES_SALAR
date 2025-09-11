#          ############################################################################################
#          #              >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<            #
#          #                                       CABLE SUSPENSION BRIDGE 04                         #
#          #------------------------------------------------------------------------------------------#
#          #                                        President Ibáñez Bridge                           #
#          #------------------------------------------------------------------------------------------#
#          #                  THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)              #
#          #                            EMAIL: salar.d.ghashghaei@gmail.com                           #
#          ############################################################################################
"""
Suspension Bridge:
A suspension bridge is a bridge type that uses strong steel cables and tall towers to support the bridge deck. It is designed to carry vehicles, people, and other loads across large natural obstacles such as wide rivers, deep valleys, or rugged terrain.
Unlike simple beam bridges, which rest on multiple supports, suspension bridges can cover very long distances with only a few main supports.

It has a total length of about 210 meters, two steel towers about 25 meters tall,
 and sets of eight main steel cables on each side that support the concrete deck
 through vertical hangers. This structure shows how suspension bridges can cross 
 large rivers even in remote regions.
-------------------------------------------------
Key Features:
- Structural Concept
  Suspension bridges have steel cables stretched between tall towers.
  The cables carry the deck’s weight and transfer the load to the towers, which send the forces down into the ground through strong anchorages.
  The deck itself is usually stiffened with longitudinal and transverse beams to resist bending.

- Design Purpose
  They are designed to carry vehicles, bicycles, and pedestrians across wide gaps such as rivers, bays, or mountain gorges.
  Suspension bridges are often used where long spans are required and building many supports is not possible.

- Load Characteristics
  These bridges carry heavy loads, including car and truck traffic.
  The loads are transferred from the deck to vertical hangers, then to the main cables, and finally to the towers and anchorages.
-------------------------------------------------
Advantages:
- Longer Spans — They can cross very long distances without intermediate supports.
- Overcoming Natural Obstacles — Ideal for spanning rivers, fjords, or deep valleys.
- Aesthetic Appeal — Their tall towers and sweeping cables make them visually impressive landmarks.
-------------------------------------------------
Construction and Use:
- Building suspension bridges requires precise engineering and heavy equipment to install towers, anchorages, and cables.
- They are used on major transport routes where long-span crossings are needed, especially in areas with difficult geography like wide rivers or canyons.
-------------------------------------------------
Differences from Other Bridges
- Load Capacity
  Suspension bridges are designed for heavy loads such as vehicles and trucks, unlike small beam or arch bridges.
- Scale and Dimensions
  They are much longer and taller, with main spans reaching hundreds or even thousands of meters.
-------------------------------------------------
Summary:
For a simplified modeling approach, it is assumed that the bridge deck, which is originally truss-shaped,
 is modeled as a simple element having the same cross-sectional area as the deck.
A suspension bridge uses the cable-and-tower system to carry heavy loads over long distances.
The  is an example, showing how this design allows safe vehicle crossings over wide rivers while reducing the number of supports needed.
"""
#%%------------------------------------------------------------------
# Load the image
def PLOT_IMAGE(image):
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg
    image = mpimg.imread(image_path)

    # Display the image
    plt.figure(figsize=(15, 8))
    plt.imshow(image)
    plt.axis('off')  # Hide axes
    plt.show()
    
image_path = 'OPENSEES_CABLE_SUSPENSION_BRIDGE_04.png'    
PLOT_IMAGE(image_path)
#%%------------------------------------------------------------------
# WIKIPEDIA: Simple suspension bridge
'https://es.wikipedia.org/wiki/Puente_Presidente_Ib%C3%A1%C3%B1ez'
# IMAGE:
'https://www.monumentos.gob.cl/monumentos/monumentos-historicos/puente-presidente-ibanez'
'https://www.minube.com/rincon/puente-presidente-ibanez-a461551'
# PAPER: Bridge Map and Data
'https://www.monumentos.gob.cl/sites/default/files/planos/2025-06/00818_MH_11201_POd_D138_2002.pdf'
#%%------------------------------------------------------------------
#import the os module
import os
import time
import numpy as np
import openseespy.opensees as op
import matplotlib.pyplot as plt
import ANALYSIS_FUNCTION as S02
import DAMPING_RATIO_FUN as S05
import RAYLEIGH_DAMPING_FUN as S06
import BILINEAR_CURVE as BC
import EIGENVALUE_ANALYSIS_FUN as S07
#%%------------------------------------------------------------------
# Create a directory at specified path with name 'directory_path'
import os
directory_path = 'C:\\OPENSEESPY_SALAR'

# Check if the directory already exists
if not os.path.exists(directory_path):
    os.mkdir(directory_path)
    print(f"Directory '{directory_path}' created successfully.")
else:
    print(f"Directory '{directory_path}' already exists. Skipping creation.")
#-------------------------------------------------------------
# Create folder name
FOLDER_NAME = 'OPENSEES_CABLE_SUSPENSION_BRIDGE'
dir = f"C:\\OPENSEESPY_SALAR\\{FOLDER_NAME}\\"
if not os.path.exists(dir):
    os.makedirs(dir)    
#-------------------------------------------------------------
# OUTPUT DATA ADDRESS:
SALAR_DIR = f'C://OPENSEESPY_SALAR//{FOLDER_NAME}//';
#-------------------------------------------------------------
## DELETE ALL FILES IN DIRECTORY 
def DELETE_FOLDER_CONTANTS(folder_path):
    import os
    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
        except Exception as e:
            print(f"Failed to delete {file_path}. Reason: {e}")
    print("Deletion done")
   
FOLDER_PATH = f'C:\\OPENSEESPY_SALAR\\{FOLDER_NAME}'  # Specify the folder path
#DELETE_FOLDER_CONTANTS(FOLDER_PATH)   

#%%------------------------------------------------------------------
def CURRENT_TIME():
    import time
    t = time.localtime()
    current_time = time.strftime("%H:%M:%S", t)
    print(f"Current time (HH:MM:SS): {current_time}\n\n")
    
# ---------------------
def plot_shapes(initial_coords, displacements):
    import numpy as np
    import matplotlib.pyplot as plt
    """
    Plot the initial and deformed shapes of the arch.

    Parameters:
    - initial_coords: Initial coordinates of the nodes.
    - displacements: List of displacements for each node.
    """
    # Plot initial and deformed shape
    plt.figure(figsize=(12, 8))
    plt.plot(initial_coords[:num_nodes, 0], initial_coords[:num_nodes, 1], 'bo-', label='Initial Cable Shape')
    plt.plot(deformed_coords[:num_nodes, 0], deformed_coords[:num_nodes, 1], 'ro-', label='Deformed Cable Shape')
    plt.plot(initial_coords[num_nodes:, 0], initial_coords[num_nodes:, 1], 'go-', label='Initial Deck Shape')
    plt.plot(deformed_coords[num_nodes:, 0], deformed_coords[num_nodes:, 1], 'yo-', label='Deformed Deck Shape')
    
    # Plot connecting elements between the cable and the beam
    for i in range(num_nodes):
        plt.plot([initial_coords[i, 0], initial_coords[num_nodes + i, 0]],
                 [initial_coords[i, 1], initial_coords[num_nodes + i, 1]], 'c--', label='Initial Connecting Elements' if i == 0 else "")
        plt.plot([deformed_coords[i, 0], deformed_coords[num_nodes + i, 0]],
                 [deformed_coords[i, 1], deformed_coords[num_nodes + i, 1]], 'm--', label='Deformed Connecting Elements' if i == 0 else "")
    
    # Mark simply supported nodes with red triangles
    support_nodes = [0, num_nodes - 1, num_nodes, 2*num_nodes - 1]
    plt.plot(initial_coords[support_nodes, 0], initial_coords[support_nodes, 1], 'r^', markersize=20, label='Simply Supported Nodes')
    
    plt.xlabel('X [mm]')
    plt.ylabel('Y [mm]')
    plt.legend()
    plt.title('Pushover Analysis: Cable and Deck Deformed Shapes')
    plt.show()
    
# ---------------------

def plot_reactions(disp_x, reaction_x, disp_y, reaction_y):
    import matplotlib.pyplot as plt
    """
    Plot the base reaction versus displacement.

    Parameters:
    - disp_x: List of x displacements of the middle node.
    - reaction_x: List of base reactions in x.
    - disp_y: List of y displacements of the middle node.
    - reaction_y: List of base reactions in y.
    """
    plt.figure(figsize=(10, 6))
    plt.plot(disp_x, reaction_x, 'b-')
    plt.xlabel('Displacement X [mm]')
    plt.ylabel('Base Reaction X [N]')
    plt.title('Base Reaction X vs. Displacement X')
    plt.show()

    plt.figure(figsize=(10, 6))
    plt.plot(disp_y, reaction_y, 'r-')
    plt.xlabel('Displacement Y [mm]')
    plt.ylabel('Base Reaction Y [N]')
    plt.title('Base Reaction Y vs. Displacement Y')
    plt.show()
    
# ---------------------
#%%------------------------------------------------------------------
#          ############################################################################################
#          #                                          IN THE NAME OF ALLAH                            #
#          #                         PUSHOVER ANALYSIS OF CABLE SUSPENSION BRIDGE 04                  #
#          #------------------------------------------------------------------------------------------#
#          #                  THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)              #
#          #                            EMAIL: salar.d.ghashghaei@gmail.com                           #
#          ############################################################################################
#%%------------------------------------------------------------------
# --------------------
#   PUSHOVER ANALYSIS
# --------------------
def PUSHOVER_ANALYSIS(LINEAR, L, H1, arc_depth, E_cable, Cable_Dia_01, Cable_Dia_02, Cable_Dia_03, num_nodes, MAX_DISP, disp_incr, MAX_ITERATIONS, MAX_TOLERANCE):
    import openseespy.opensees as ops
    import numpy as np
    import matplotlib.pyplot as plt
    
    A_cable_01 = (np.pi * Cable_Dia_01 **2) / 4  # [mm^2] Horizontal Longitudinal Cable Area Bottom
    A_cable_02 = (np.pi * Cable_Dia_02 **2) / 4  # [mm^2] Horizontal Longitudinal Cable Area Top
    A_cable_03 = (np.pi * Cable_Dia_03 **2) / 4  # [mm^2] Vertical Longitudinal Cable Area Top
    
    # Define model builder
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 2)
    dx = L / (num_nodes - 1)
    
    # Create nodes
    for i in range(num_nodes):
        x = i * dx
        y = H1 - arc_depth * np.sin(np.pi * x / L)
        ops.node(i + 1, x, y)   
        ops.node(num_nodes + i + 1, i * dx, 0.0)

    # Define boundary conditions (fixed at both ends)
    ops.fix(1, 1, 1)              # TOP CABLE
    ops.fix(num_nodes, 1, 1)      # TOP CABLE
    ops.fix(num_nodes + 1, 1, 1)  # BOTTOM CABLE - DECK
    ops.fix(2 * num_nodes, 1, 1)  # BOTTOM CABLE - DECK
    
    # Define material properties
    if LINEAR == True:
        ops.uniaxialMaterial('Elastic', 1, E_cable)
        #                                  TENSION    COMPRESSION
        #ops.uniaxialMaterial('Elastic', 1, E_cable, 0, 0.5 * E_cable)
    if LINEAR == False:
        #E_cable = 210e5   # [N/mm²] Modulus of Elasticity
        #Fy_cable = 3550   # [N/mm^2] Yield strength of the cable
        #b0 = 0.01
        #ops.uniaxialMaterial('Steel01', 1, Fy_cable, E_cable, b0)
        fy = 3550                                 # [N/mm²] Steel Rebar Yield Strength   
        Es = 210e5                                # [N/mm²] Modulus of Elasticity
        ey = fy/Es                                # [mm/mm] Steel Rebar Yield Strain
        fu = 1.1818*fy                            # [N/mm²] Steel Rebar Ultimate Strength
        esu = ey*75.2                             # [mm/mm] Steel Rebar Ultimate Strain
        Esh = (fu - fy)/(esu - ey)
        Bs = Esh / Es
        pinchX = 0.8   # Pinching factor in X direction
        pinchY = 0.5   # Pinching factor in Y direction
        damage1 = 0.0  # Damage due to ductility
        damage2 = 0.0  # Damage due to energy
        beta = 0.1     # Stiffness degradation parameter
        ops.uniaxialMaterial('Hysteretic', 1, fy, ey, fu, esu, 0.2*fu, 1.1*esu, -fy, -ey, -fu, -esu, -0.2*fu, -1.1*esu, pinchX, pinchY, damage1, damage2, beta)
        # INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Hysteretic_Material
        
    # Define truss elements
    for i in range(num_nodes - 1):
        ops.element('corotTruss', i + 1, i + 1, i + 2, A_cable_02, 1)  # Cable element
        ops.element('corotTruss', num_nodes + i + 1, num_nodes + i + 1, num_nodes + i + 2, A_cable_01, 1)  # Deck element

    # Connect each node of the cable with the corresponding node of the deck using truss elements
    for i in range(num_nodes):
        ops.element('corotTruss', 2 * num_nodes + i + 1, i + 1, num_nodes + i + 1, A_cable_03, 1)
    
    #mid_node = num_nodes // 2 + 1
    mid_node = int(num_nodes + 0.5 * num_nodes)
    # Define load pattern with displacement at the middle span
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(mid_node, 0.0, -1)  # Apply a horizontal load at node
    n_steps = int(np.abs(MAX_DISP / disp_incr))  # Analysis Steps
    
    # Define analysis parameters
    ops.system('BandGeneral')
    ops.numberer('Plain')
    ops.constraints('Plain')
    ops.integrator('DisplacementControl', mid_node, 2, disp_incr)
    ops.test('EnergyIncr', MAX_TOLERANCE, MAX_ITERATIONS)
    ops.algorithm('ModifiedNewton')
    ops.analysis('Static')

    # OUTPUT DATA
    ops.recorder('Node', '-file', f"{SALAR_DIR}DTH_PUSH.txt",'-time', '-node', mid_node, '-dof', 1,2, 'disp')# Displacement Time History
    ops.recorder('Node', '-file', f"{SALAR_DIR}BTH_PUSH_01.txt",'-time', '-node', 1, '-dof', 1,2, 'reaction')# Base Shear Time History
    ops.recorder('Node', '-file', f"{SALAR_DIR}BTH_PUSH_02.txt",'-time', '-node', num_nodes, '-dof', 1,2, 'reaction')# Base Shear Time History
    ops.recorder('Node', '-file', f"{SALAR_DIR}BTH_PUSH_03.txt",'-time', '-node', num_nodes + 1, '-dof', 1,2, 'reaction')# Base Shear Time History
    ops.recorder('Node', '-file', f"{SALAR_DIR}BTH_PUSH_04.txt",'-time', '-node', 2*num_nodes, '-dof', 1,2, 'reaction')# Base Shear Time History

    DISPLACEMENTS = []
    DISP_X , DISP_Y = [], []
    BASE_REACTION_X, BASE_REACTION_Y = [], []
    
    # Perform the analysis with increments
    for i in range(n_steps):
        OK = ops.analyze(1)
        S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS)
        DISPLACEMENTS.append([ops.nodeDisp(j + 1) for j in range(num_nodes * 2)])
        DISP_X.append(ops.nodeDisp(mid_node, 1))
        DISP_Y.append(ops.nodeDisp(mid_node, 2))
        x1 = ops.nodeResponse(1, 1, 6)
        x2 = ops.nodeResponse(num_nodes, 1, 6)
        x3 = ops.nodeResponse(num_nodes + 1, 1, 6)
        x4 = ops.nodeResponse(2*num_nodes, 1, 6)
        
        y1 = ops.nodeResponse(1, 2, 6)
        y2 = ops.nodeResponse(num_nodes, 2, 6)
        y3 = ops.nodeResponse(num_nodes + 1, 2, 6)
        y4 = ops.nodeResponse(2*num_nodes, 2, 6)
        BASE_REACTION_X.append(x1+x2+x3+x4) # CABLE REACION-X
        BASE_REACTION_Y.append(y1+y2+y3+y4) # CABLE REACION-Y
        #print(f'STEP: {i+1}')
    
    # Get initial and final node coordinates for plotting
    initial_coords = np.array([ops.nodeCoord(i + 1) for i in range(num_nodes * 2)])
    deformed_coords = initial_coords + np.array(DISPLACEMENTS[-1])
    
    # Output all unconstrained node displacements
    for i in range(num_nodes * 2):
        disp = ops.nodeDisp(i + 1)
        print(f"Node {i + 1}: X Displacement = {disp[0]}, Y Displacement = {disp[1]}")
    
    DATA = (initial_coords, deformed_coords, DISPLACEMENTS,
            DISP_X, DISP_Y, BASE_REACTION_X, BASE_REACTION_Y)
    
    return DATA
#%%------------------------------------------------------------------
# --------------------
#   PUSHOVER ANALYSIS
# --------------------
# Parameters for the analysis
L = 210000.0          # [mm] Bridge span lentgh
H1 = 25000.0          # [mm] Height of Top Cable
arc_depth = 20000.0   # [mm]
E_cable = 210e5       # [N/mm^2] Modulus of elasticity Cable
Cable_Dia_01 = 145600 # [mm] Horizontal Longitudinal Cable Diameter Bottom
Cable_Dia_02 = 300    # [mm] Horizontal Longitudinal Cable Diameter Top
Cable_Dia_03 = 50     # [mm] Vertical Longitudinal Cable Diameter Top 
num_nodes = 22        # Cable Arc Number of nodes
MAX_DISP = -2500.0    # [mm] Maximum Displacment
disp_incr = -0.01     # [mm] Incremental Displacement

MAX_ITERATIONS = 50000  # Maximum number of iterations
MAX_TOLERANCE = 1.0e-6  # Tolerance for convergence

LINEAR = True # False: Cable Nonlinear Materials Properties

starttime = time.process_time()
# Run the analysis
results = PUSHOVER_ANALYSIS(LINEAR, L, H1, arc_depth, E_cable, Cable_Dia_01, Cable_Dia_02, Cable_Dia_03, num_nodes, MAX_DISP, disp_incr, MAX_ITERATIONS, MAX_TOLERANCE)
initial_coords, deformed_coords, DISPLACEMENTS, DISP_X, DISP_Y, BASE_REACTION_X, BASE_REACTION_Y = results
totaltime = time.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n')
#%%------------------------------------------------------------------
# Plot results
plot_shapes(initial_coords, DISPLACEMENTS)
plot_reactions(DISP_X, BASE_REACTION_X, DISP_Y, BASE_REACTION_Y)
#%%------------------------------------------------------------------
#          ############################################################################################
#          #                                          IN THE NAME OF ALLAH                            #
#          #                     FREE VIBRATION ANALYSIS OF CABLE SUSPENSION BRIDGE 04                #
#          #------------------------------------------------------------------------------------------#
#          #                  THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)              #
#          #                            EMAIL: salar.d.ghashghaei@gmail.com                           #
#          ############################################################################################
#%%------------------------------------------------------------------
# --------------------------
#   FREE VIBRATION ANALYSIS
# --------------------------
def FREE_VIBRATION_ANALYSIS(damping, damping_ratio, LINEAR, duration, dt, TOTAL_MASS, u0, L, H1, arc_depth, E_cable, Cable_Dia_01, Cable_Dia_02, Cable_Dia_03, num_nodes, MAX_ITERATIONS, MAX_TOLERANCE):
    import openseespy.opensees as ops
    import numpy as np
    import matplotlib.pyplot as plt
    
    A_cable_01 = (np.pi * Cable_Dia_01 **2) / 4  # [mm^2] Horizontal Longitudinal Cable Area Bottom
    A_cable_02 = (np.pi * Cable_Dia_02 **2) / 4  # [mm^2] Horizontal Longitudinal Cable Area Top
    A_cable_03 = (np.pi * Cable_Dia_03 **2) / 4  # [mm^2] Vertical Longitudinal Cable Area Top
    
    # Define model builder
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 2)
    KE = 1;
    dx = L / (num_nodes - 1)
    MASS = TOTAL_MASS / num_nodes
    # Create nodes
    for i in range(num_nodes):
        X = i * dx
        Y1 = H1 - arc_depth * np.sin(np.pi * X / L)
        ops.node(i + 1, X, Y1)
        ops.node(num_nodes + i + 1, i * dx, 0.0)
        # Define mass
        ops.mass(num_nodes + i + 1, MASS, MASS, 0)

    # Define boundary conditions (fixed at both ends)
    ops.fix(1, 1, 1)              # TOP CABLE
    ops.fix(num_nodes, 1, 1)      # TOP CABLE
    ops.fix(num_nodes + 1, 1, 1)  # BOTTOM CABLE
    ops.fix(2 * num_nodes, 1, 1)  # BOTTOM CABLE
    
    # Define material properties
    if LINEAR == True:
        ops.uniaxialMaterial('Elastic', 1, E_cable)
        #                                  TENSION    COMPRESSION
        #ops.uniaxialMaterial('Elastic', 1, E_cable, 0, 0.5 * E_cable)
    if LINEAR == False:
        #E_cable = 210e5   # [N/mm²] Modulus of Elasticity
        #Fy_cable = 3550   # [N/mm^2] Yield strength of the cable
        #b0 = 0.01
        #ops.uniaxialMaterial('Steel01', 1, Fy_cable, E_cable, b0)
        fy = 3550                                 # [N/mm²] Steel Rebar Yield Strength   
        Es = 210e5                                # [N/mm²] Modulus of Elasticity
        ey = fy/Es                                # [mm/mm] Steel Rebar Yield Strain
        fu = 1.1818*fy                            # [N/mm²] Steel Rebar Ultimate Strength
        esu = ey*75.2                             # [mm/mm] Steel Rebar Ultimate Strain
        Esh = (fu - fy)/(esu - ey)
        Bs = Esh / Es
        pinchX = 0.8   # Pinching factor in X direction
        pinchY = 0.5   # Pinching factor in Y direction
        damage1 = 0.0  # Damage due to ductility
        damage2 = 0.0  # Damage due to energy
        beta = 0.1     # Stiffness degradation parameter
        ops.uniaxialMaterial('Hysteretic', 1, fy, ey, fu, esu, 0.2*fu, 1.1*esu, -fy, -ey, -fu, -esu, -0.2*fu, -1.1*esu, pinchX, pinchY, damage1, damage2, beta)
        # INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Hysteretic_Material
 
    # Define truss elements
    for i in range(num_nodes - 1):
        ops.element('corotTruss', i + 1, i + 1, i + 2, A_cable_02, 1)  # Top Cable element
        ops.element('corotTruss', num_nodes + i + 1, num_nodes + i + 1, num_nodes + i + 2, A_cable_01, 1)  # Bottom Cable Element

    # Connect each node of the cable with the corresponding node of the deck using truss elements
    for i in range(num_nodes):
        ops.element('corotTruss', 2 * num_nodes + i + 1, i + 1, num_nodes + i + 1, A_cable_03, 1) # Vertical Cable Element
    
    #mid_node = num_nodes // 2 + 1
    mid_node = int(num_nodes + 0.5 * num_nodes)
    # Define load pattern with displacement at the middle span
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(mid_node, 0.0, -1)  # Apply a horizontal load at node
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('BandGeneral')
    ops.algorithm('Newton')
    ops.test('EnergyIncr', MAX_TOLERANCE, MAX_ITERATIONS)
    ops.integrator('DisplacementControl', mid_node, 2, u0)
    ops.analysis('Static')
    ops.analyze(1)

    ops.setTime(0.0)

    # Wipe analysis and reset time
    ops.wipeAnalysis()
    ops.remove('loadPattern', 1)
    ops.system('UmfPack')

    # Dynamic analysis
    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.system('UmfPack')
    ops.test('EnergyIncr', MAX_TOLERANCE, MAX_ITERATIONS)
    #ops.integrator('CentralDifference')
    ops.integrator('Newmark', 0.5, 0.25) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/newmark.html
    #alpha = 1;gamma=1.5-alpha; beta=((2-alpha)**2)/4;
    #ops.integrator('HHT', alpha, gamma, beta) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/hht.html
    ops.algorithm('Newton')

    if damping == True:
        # Calculate Rayleigh damping factors
        PERIOD_01, PERIOD_02 = S06.RAYLEIGH_DAMPING(2, damping_ratio, 0.6*damping_ratio, 0, 1)

    ops.analysis('Transient')

    # OUTPUT DATA
    ops.recorder('Node', '-file', f"{SALAR_DIR}DTH_DYN.txt",'-time', '-node', mid_node, '-dof', 1,2, 'disp')# Displacement Time History
    ops.recorder('Node', '-file', f"{SALAR_DIR}VTH_DYN.txt",'-time', '-node', mid_node, '-dof', 1,2, 'vel')  # Velocity Time History
    ops.recorder('Node', '-file', f"{SALAR_DIR}ATH_DYN.txt",'-time', '-node', mid_node, '-dof', 1,2, 'accel') # Acceleration Time History
    ops.recorder('Node', '-file', f"{SALAR_DIR}BTH_DYN_01.txt",'-time', '-node', 1, '-dof', 1,2, 'reaction')# Base Shear Time History
    ops.recorder('Node', '-file', f"{SALAR_DIR}BTH_DYN_02.txt",'-time', '-node', num_nodes, '-dof', 1,2, 'reaction')# Base Shear Time History
    ops.recorder('Node', '-file', f"{SALAR_DIR}BTH_DYN_03.txt",'-time', '-node', num_nodes + 1, '-dof', 1,2, 'reaction')# Base Shear Time History
    ops.recorder('Node', '-file', f"{SALAR_DIR}BTH_DYN_04.txt",'-time', '-node', 2*num_nodes, '-dof', 1,2, 'reaction')# Base Shear

    DISPLACEMENTS = []
    DISP_X , DISP_Y, VELOCITY, ACCELERATION = [], [], [], []
    BASE_REACTION_X, BASE_REACTION_Y = [], []

    stable = 0
    current_time = 0.0
     # Perform the analysis with increments
    while stable == 0 and current_time < duration:
        OK = ops.analyze(1, dt)
        S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS)
        current_time = ops.getTime()
        DISPLACEMENTS.append([ops.nodeDisp(j + 1) for j in range(num_nodes * 2)])
        DISP_X.append(ops.nodeDisp(mid_node, 1))
        DISP_Y.append(ops.nodeDisp(mid_node, 2))
        VELOCITY.append(ops.nodeVel(mid_node, 2))
        ACCELERATION.append(ops.nodeAccel(mid_node, 2))
        x1 = ops.nodeResponse(1, 1, 6)
        x2 = ops.nodeResponse(num_nodes, 1, 6)
        x3 = ops.nodeResponse(num_nodes + 1, 1, 6)
        x4 = ops.nodeResponse(2*num_nodes, 1, 6)
        y1 = ops.nodeResponse(1, 2, 6)
        y2 = ops.nodeResponse(num_nodes, 2, 6)
        y3 = ops.nodeResponse(num_nodes + 1, 2, 6)
        y4 = ops.nodeResponse(2*num_nodes, 2, 6)
        BASE_REACTION_X.append(x1+x2+x3+x4) # CABLE REACION-X
        BASE_REACTION_Y.append(y1+y2+y3+y4) # CABLE REACION-Y
        #print(f'time: {current_time}')
    
    # Get initial and final node coordinates for plotting
    initial_coords = np.array([ops.nodeCoord(i + 1) for i in range(num_nodes * 2)])
    deformed_coords = initial_coords + np.array(DISPLACEMENTS[-1])
    
    # Output all unconstrained node displacements
    for i in range(num_nodes * 2):
        disp = ops.nodeDisp(i + 1)
        print(f"Node {i + 1}: X Displacement = {disp[0]}, Y Displacement = {disp[1]}")
    
    DATA = (initial_coords, deformed_coords, DISPLACEMENTS,
            VELOCITY, ACCELERATION, DISP_X, DISP_Y,
            BASE_REACTION_X, BASE_REACTION_Y)
    
    return DATA
#%%------------------------------------------------------------------
# Define the plotting function
def plot_time_history(TIME, DISP_X_undamped, DISP_X_damped, DISP_Y_undamped, DISP_Y_damped,
                      VELOCITY_undamped, VELOCITY_damped, ACCELERATION_undamped, ACCELERATION_damped,
                      BASE_REACTION_X_undamped, BASE_REACTION_X_damped, BASE_REACTION_Y_undamped, BASE_REACTION_Y_damped):
    import matplotlib.pyplot as plt
    import numpy as np
    
    plt.figure(figsize=(14, 26))

    # Displacement X
    plt.subplot(6, 1, 1)
    P1 = np.max(np.abs(np.array(DISP_X_undamped)))
    P2 = np.max(np.abs(np.array(DISP_X_damped)))
    plt.plot(TIME, DISP_X_undamped, color='black', label=f'Undamped: {P1:.4e}')
    plt.plot(TIME, DISP_X_damped, color='red', label=f'Damped: {P2:.4e}')
    plt.xlabel('Time [s]')
    plt.ylabel('Displacement X [mm]')
    plt.legend()

    # Displacement Y
    plt.subplot(6, 1, 2)
    P1 = np.max(np.abs(np.array(DISP_Y_undamped)))
    P2 = np.max(np.abs(np.array(DISP_Y_damped)))
    plt.plot(TIME, DISP_Y_undamped, color='black', label=f'Undamped: {P1:.4e}')
    plt.plot(TIME, DISP_Y_damped, color='red', label=f'Damped: {P2:.4e}')
    plt.xlabel('Time [s]')
    plt.ylabel('Displacement Y [mm]')
    plt.legend()

    # Velocity
    plt.subplot(6, 1, 3)
    P1 = np.max(np.abs(np.array(VELOCITY_damped)))
    P2 = np.max(np.abs(np.array(VELOCITY_undamped)))
    plt.plot(TIME, VELOCITY_undamped, color='black', label=f'Undamped: {P1:.4e}')
    plt.plot(TIME, VELOCITY_damped, color='red', label=f'Damped: {P2:.4e}')
    plt.xlabel('Time [s]')
    plt.ylabel('Velocity [mm/s]')
    plt.legend()

    # Acceleration
    plt.subplot(6, 1, 4)
    P1 = np.max(np.abs(np.array(ACCELERATION_undamped)))
    P2 = np.max(np.abs(np.array(ACCELERATION_damped)))
    plt.plot(TIME, ACCELERATION_undamped, color='black', label=f'Undamped: {P1:.4e}')
    plt.plot(TIME, ACCELERATION_damped, color='red', label=f'Damped: {P2:.4e}')
    plt.xlabel('Time [s]')
    plt.ylabel('Acceleration [mm/s^2]')
    plt.legend()

    # Base Reaction X
    plt.subplot(6, 1, 5)
    P1 = np.max(np.abs(np.array(BASE_REACTION_X_undamped)))
    P2 = np.max(np.abs(np.array(BASE_REACTION_X_damped)))
    plt.plot(TIME, BASE_REACTION_X_undamped, color='black', label=f'Undamped: {P1:.4e}')
    plt.plot(TIME, BASE_REACTION_X_damped, color='red', label=f'Damped: {P2:.4e}')
    plt.xlabel('Time [s]')
    plt.ylabel('Base Reaction X [N]')
    plt.legend()

    # Base Reaction Y
    plt.subplot(6, 1, 6)
    P1 = np.max(np.abs(np.array(BASE_REACTION_Y_undamped)))
    P2 = np.max(np.abs(np.array(BASE_REACTION_Y_damped)))
    plt.plot(TIME, BASE_REACTION_Y_undamped, color='black', label=f'Undamped: {P1:.4e}')
    plt.plot(TIME, BASE_REACTION_Y_damped, color='red', label=f'Damped: {P2:.4e}')
    plt.xlabel('Time [s]')
    plt.ylabel('Base Reaction Y [N]')
    plt.legend()
    
    plt.suptitle('Time History of Free Vibration Analysis: Damped vs Undamped')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()
#%%------------------------------------------------------------------
# --------------------------
#   FREE VIBRATION ANALYSIS
# --------------------------
""" OMITTED, BECAUSE THEY ARE AS SAME AS TOP ANALYSIS PARAMETERS
# Parameters for the analysis
L = 210000.0        # [mm] Bridge span lentgh
H1 = 25000.0        # [mm] Height of Top Cable
arc_depth = 20000.0 # [mm]
E_cable = 210e5     # [N/mm^2] Modulus of elasticity Cable
Cable_Dia_01 = 145600 # [mm] Horizontal Longitudinal Cable Diameter Bottom
Cable_Dia_02 = 300  # [mm] Horizontal Longitudinal Cable Diameter Top
Cable_Dia_03 = 50   # [mm] Vertical Longitudinal Cable Diameter Top 
num_nodes = 22      # Cable Arc Number of nodes
"""
TOTAL_MASS = 50_000_000.0 # [kg] Total Mass of Structure
u0 = 20.0              # [mm] Initial displacement
damping_ratio = 0.01   # Damping ratio
duration = 50.0        # [s] Duration of the analysis in seconds
dt = 0.01              # Time step in seconds

MAX_ITERATIONS = 10000  # Maximum number of iterations
MAX_TOLERANCE = 1.0e-6  # Tolerance for convergence

import time
starttime = time.process_time()

# Run the undamped analysis
damping = False
LINEAR = True # False: Cable Nonlinear Materials Properties
results_undamped = FREE_VIBRATION_ANALYSIS(damping, damping_ratio, LINEAR, duration, dt, TOTAL_MASS, u0, L, H1, arc_depth, E_cable, Cable_Dia_01, Cable_Dia_02, Cable_Dia_03, num_nodes, MAX_ITERATIONS, MAX_TOLERANCE)
initial_coords, deformed_coords, DISPLACEMENTS_undamped, VELOCITY_undamped, ACCELERATION_undamped, DISP_X_undamped, DISP_Y_undamped, BASE_REACTION_X_undamped, BASE_REACTION_Y_undamped = results_undamped

# Run the damped analysis
damping = True
LINEAR = True # False: Cable Nonlinear Materials Properties
results_damped = FREE_VIBRATION_ANALYSIS(damping, damping_ratio, LINEAR, duration, dt, TOTAL_MASS, u0, L, H1, arc_depth, E_cable, Cable_Dia_01, Cable_Dia_02, Cable_Dia_03, num_nodes, MAX_ITERATIONS, MAX_TOLERANCE)
initial_coords, deformed_coords, DISPLACEMENTS_damped, VELOCITY_damped, ACCELERATION_damped, DISP_X_damped, DISP_Y_damped, BASE_REACTION_X_damped, BASE_REACTION_Y_damped = results_damped

totaltime = time.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n')
#%%------------------------------------------------------------------
# Plotting the time history
TIME = np.arange(dt, duration+2*dt, dt)

plot_time_history(TIME, DISP_X_undamped, DISP_X_damped, DISP_Y_undamped, DISP_Y_damped,
                  VELOCITY_undamped, VELOCITY_damped, ACCELERATION_undamped, ACCELERATION_damped,
                  BASE_REACTION_X_undamped, BASE_REACTION_X_damped, BASE_REACTION_Y_undamped, BASE_REACTION_Y_damped)
#%%------------------------------------------------------------------
# Plot damped results
plot_shapes(initial_coords, DISPLACEMENTS_damped)
plot_reactions(DISP_X_damped, BASE_REACTION_X_damped, DISP_Y_damped, BASE_REACTION_Y_damped)
#%%------------------------------------------------------------------
# Plot undamped results
plot_shapes(initial_coords, DISPLACEMENTS_undamped)
plot_reactions(DISP_X_undamped, BASE_REACTION_X_undamped, DISP_Y_undamped, BASE_REACTION_Y_undamped)
#%%------------------------------------------------------------------
def HISROGRAM_BOXPLOT(X, HISTO_COLOR, LABEL):
    import numpy as np
    import matplotlib.pyplot as plt
    X = np.array(X)
    print("-------------------------")
    from scipy.stats import skew, kurtosis
    MINIMUM = np.min(X)
    MAXIMUM = np.max(X)
    #MODE = max(set(X), key=list(X).count)
    MEDIAN = np.quantile(X, .50)#q2
    MEAN = np.mean(X)
    STD = np.std(X)
    q1 = np.quantile(X, .25)
    q3 = np.quantile(X, .75)
    SKEW = skew(X)
    KURT = kurtosis(X)
    #SKEW = (MEAN - MODE) / STD
    #KURT = (np.mean((X - MEAN)**4) / STD**4)
    # Estimate confidence intervals of the output variable
    lower_bound = np.quantile(X, .05)
    upper_bound = np.quantile(X, .95)
    print("Box-Chart Datas: ")
    print(f'Minimum: {MINIMUM:.4e}')
    print(f'First quantile: {q1:.4e}')
    #print(f'Mode: {MODE:.4e}')
    print(f'Median: {MEDIAN:.4e}')
    print(f'Mean: {MEAN:.4e}')
    print(f'Std: {STD:.4e}')
    print(f'Third quantile: {q3:.4e}')
    print(f'Maximum: {MAXIMUM :.4e}')
    print(f'Skewness: {skew(X) :.4e}')
    print(f'kurtosis: {kurtosis(X) :.4e}')
    print(f"90% Confidence Interval: ({lower_bound:.4e}, {upper_bound:.4e})")
    print("-------------------------")

    plt.figure(figsize=(10,6))
    # Plot histogram of data
    count, bins, ignored = plt.hist(X, bins=100, color=HISTO_COLOR, density=True, align='mid')#, edgecolor="black"
    
    # Plot lognormal PDF
    x = np.linspace(min(bins), max(bins), 10000)
    pdf = (np.exp(-(x - MEAN)**2 / (2 * STD**2)) / (STD * np.sqrt(2 * np.pi)))
    plt.plot(x, pdf, linewidth=2, color='r', label="Normal PDF")
    
    # Plot vertical lines for risk measures
    plt.axvline(q1, color="black", linestyle="--", label=f"Quantile 0.25: {q1:.4e}")
    plt.axvline(MEDIAN, color="green", linestyle="--", label=f"Median: {MEDIAN:.4e}")
    plt.axvline(q3, color="black", linestyle="--", label=f"Quantile 0.75: {q3:.4e}")
    #plt.axvline(MODE, color="purple", linestyle="--", label=f"Mode: {MODE:.4e}")
    plt.axvline(MEAN, color="red", linestyle="--", label=f"Mean: {MEAN:.4e}")
    plt.axvline(MEAN-STD, color="blue", linestyle="--", label=f"Mean-Std: {MEAN-STD:.4e}")
    plt.axvline(MEAN+STD, color="blue", linestyle="--", label=f"Mean+Std: {MEAN+STD:.4e}")
    plt.xlabel(LABEL)
    plt.ylabel("Frequency")
    prob = np.sum(X > 0) / len(X)
    plt.title(f"Histogram - Probability of Positive {LABEL} is {100*prob:.2f} %")
    plt.legend()
    #plt.grid()
    plt.show()

    #Plot boxplot with outliers
    plt.figure(figsize=(10,6))
    plt.boxplot(X, vert=0)
    # Write the quartile data on the chart
    plt.text(q1, 1.05, f" Q1: {q1:.4e}")
    plt.text(MEDIAN, 1.1, f" Q2: {MEDIAN:.4e}")
    plt.text(q3, 1.05, f" Q3: {q3:.4e}")
    #plt.text(MODE, 1.15, f" Mode: {MODE:.4e}")
    
    #plt.text(MEAN, 0.9, f" Mean: {MEAN:.4e}")
    #plt.text(MEAN-STD, 0.9, f" Mean-Std: {MEAN-STD:.4e}")
    #plt.text(MEAN+STD, 0.9, f" Mean+Std: {MEAN+STD:.4e}")
    plt.scatter(MEAN, 1, color="red", marker="+", s=200, label=f"Mean: {MEAN:.4e}")
    plt.scatter(MEAN-STD, 1, color="green", marker="X", s=200, label=f"Mean-Std: {MEAN-STD:.4e}")
    plt.scatter(MEAN+STD, 1, color="blue", marker="*", s=200, label=f"Mean+Std:  {MEAN+STD:.4e}")
    plt.xlabel(LABEL)
    plt.ylabel("Data")
    plt.title(f"Boxplot of {LABEL}")
    plt.legend()
    plt.grid()
    plt.show()
# ----------------------------------------------- 
def HISTOGRAM_BOXPLOT_PLOTLY( DATA, XLABEL='X', TITLE='A', COLOR='cyan'):
    # Plotting histogram and boxplot
    import plotly.express as px
    fig = px.histogram(x=DATA, marginal="box", color_discrete_sequence=[COLOR])
    fig.update_layout(title=TITLE, xaxis_title=XLABEL, yaxis_title="Frequency")
    fig.show()
#%%------------------------------------------------------------------
HISROGRAM_BOXPLOT(DISP_X_undamped, HISTO_COLOR='black', LABEL='Displacement X (Undamped)')
HISROGRAM_BOXPLOT(DISP_X_damped, HISTO_COLOR='grey', LABEL='Displacement X (Damped)')
#HISTOGRAM_BOXPLOT_PLOTLY(DISP_X_undamped, XLABEL='Displacement X (Undamped)', TITLE='Displacement X (Undamped)', COLOR='black')
#HISTOGRAM_BOXPLOT_PLOTLY(DISP_X_damped, XLABEL='Displacement X (Damped)', TITLE='Displacement X (Damped)', COLOR='grey')
#%%------------------------------------------------------------------
HISROGRAM_BOXPLOT(DISP_Y_undamped, HISTO_COLOR='orange', LABEL='Displacement Y (Undamped)')
HISROGRAM_BOXPLOT(DISP_Y_damped, HISTO_COLOR='brown', LABEL='Displacement Y (Damped)')
#HISTOGRAM_BOXPLOT_PLOTLY(DISP_Y_undamped, XLABEL='Displacement Y (Undamped)', TITLE='Displacement Y (Undamped)', COLOR='orange')
#HISTOGRAM_BOXPLOT_PLOTLY(DISP_Y_damped, XLABEL='Displacement Y (Damped)', TITLE='Displacement Y (Damped)', COLOR='brown')
#%%------------------------------------------------------------------
HISROGRAM_BOXPLOT(VELOCITY_undamped, HISTO_COLOR='red', LABEL='Velocity (Undamped)')
HISROGRAM_BOXPLOT(VELOCITY_damped, HISTO_COLOR='pink', LABEL='Velocity (Damped)')
#HISTOGRAM_BOXPLOT_PLOTLY(VELOCITY_undamped, XLABEL='Velocity (Undamped)', TITLE='Velocity (Undamped)', COLOR='red')
#HISTOGRAM_BOXPLOT_PLOTLY(VELOCITY_damped, XLABEL='Velocity (Damped)', TITLE='Velocity (Damped)', COLOR='pink')
#%%------------------------------------------------------------------
HISROGRAM_BOXPLOT(ACCELERATION_undamped, HISTO_COLOR='gold', LABEL='Acceleration (Undamped)')
HISROGRAM_BOXPLOT(ACCELERATION_damped, HISTO_COLOR='brown', LABEL='Acceleration (Damped)')
#HISTOGRAM_BOXPLOT_PLOTLY(ACCELERATION_undamped, XLABEL='Acceleration (Undamped)', TITLE='Acceleration (Undamped)', COLOR='gold')
#HISTOGRAM_BOXPLOT_PLOTLY(ACCELERATION_damped, XLABEL='Acceleration (Damped)', TITLE='Acceleration (Damped)', COLOR='brown')
#%%------------------------------------------------------------------
HISROGRAM_BOXPLOT(BASE_REACTION_X_undamped, HISTO_COLOR='blue', LABEL='Base Reaction X (Undamped)')
HISROGRAM_BOXPLOT(BASE_REACTION_X_damped, HISTO_COLOR='cyan', LABEL='Base Reaction X (Damped)')
#HISTOGRAM_BOXPLOT_PLOTLY(BASE_REACTION_X_undamped, XLABEL='Base Reaction X (Undamped)', TITLE='Base Reaction X (Undamped)', COLOR='blue')
#HISTOGRAM_BOXPLOT_PLOTLY(BASE_REACTION_X_damped, XLABEL='Base Reaction X (Damped)', TITLE='Base Reaction X (Damped)', COLOR='cyan')
#%%------------------------------------------------------------------
HISROGRAM_BOXPLOT(BASE_REACTION_Y_undamped, HISTO_COLOR='green', LABEL='Base Reaction Y (Undamped)')
HISROGRAM_BOXPLOT(BASE_REACTION_Y_damped, HISTO_COLOR='lime', LABEL='Base Reaction Y (Damped)')
#HISTOGRAM_BOXPLOT_PLOTLY(BASE_REACTION_Y_undamped, XLABEL='Base Reaction Y (Undamped)', TITLE='Base Reaction Y (Undamped)', COLOR='green')
#HISTOGRAM_BOXPLOT_PLOTLY(BASE_REACTION_Y_damped, XLABEL='Base Reaction Y (Damped)', TITLE='Base Reaction Y (Damped)', COLOR='lime')
#%%------------------------------------------------------------------
#          ############################################################################################
#          #                                          IN THE NAME OF ALLAH                            #
#          #                          DYNAMIC ANALYSIS OF CABLE SUSPENSION BRIDGE 04                  #
#          #------------------------------------------------------------------------------------------#
#          #                  THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)              #
#          #                            EMAIL: salar.d.ghashghaei@gmail.com                           #
#          ############################################################################################
#%%------------------------------------------------------------------
# -------------------
#   DYNAMIC ANALYSIS
# -------------------
def DYNAMIC_ANALYSIS(damping, damping_ratio, LINEAR, periodTF, duration, dt, TOTAL_MASS, L, H1, arc_depth, E_cable, Cable_Dia_01, Cable_Dia_02, Cable_Dia_03, num_nodes, MAX_ITERATIONS, MAX_TOLERANCE):
    import openseespy.opensees as ops
    import numpy as np
    import matplotlib.pyplot as plt
    
    #GMfact = 1;
    GMfact = 9810    # [mm/s²]standard acceleration of gravity or standard acceleration
    SSF_X = 1.0      # Seismic Acceleration Scale Factor in X Direction
    SSF_Y = 1.0      # Seismic Acceleration Scale Factor in Y Direction
    iv0_X = 0.0005   # [mm/s] Initial velocity applied to the node  in X Direction
    iv0_Y = 0.0005   # [mm/s] Initial velocity applied to the node  in Y Direction
    st_iv0 = 0.0     # [s] Initial velocity applied starting time
    SEI = 'Y'        # Seismic Direction
    
    A_cable_01 = (np.pi * Cable_Dia_01 **2) / 4  # [mm^2] Horizontal Longitudinal Cable Area Bottom
    A_cable_02 = (np.pi * Cable_Dia_02 **2) / 4  # [mm^2] Horizontal Longitudinal Cable Area Top
    A_cable_03 = (np.pi * Cable_Dia_03 **2) / 4  # [mm^2] Vertical Longitudinal Cable Area Top
    
    # Define model builder
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 2)
    dx = L / (num_nodes - 1)
    MASS = TOTAL_MASS / num_nodes
    # Create nodes
    for i in range(num_nodes):
        X = i * dx
        Y1 = H1 - arc_depth * np.sin(np.pi * X / L)
        ops.node(i + 1, X, Y1)
        ops.node(num_nodes + i + 1, i * dx, 0.0)
        # Define mass
        ops.mass(num_nodes + i + 1, MASS, MASS, 0)

    # Define boundary conditions (fixed at both ends)
    ops.fix(1, 1, 1)              # TOP CABLE
    ops.fix(num_nodes, 1, 1)      # TOP CABLE
    ops.fix(num_nodes + 1, 1, 1)  # BOTTOM CABLE
    ops.fix(2 * num_nodes, 1, 1)  # BOTTOM CABLE
    
    # Define material properties
    if LINEAR == True:
        ops.uniaxialMaterial('Elastic', 1, E_cable)
        #                                  TENSION    COMPRESSION
        #ops.uniaxialMaterial('Elastic', 1, E_cable, 0, 0.5 * E_cable)
    if LINEAR == False:
        #E_cable = 210e5   # [N/mm²] Modulus of Elasticity
        #Fy_cable = 3550   # [N/mm^2] Yield strength of the cable
        #b0 = 0.01
        #ops.uniaxialMaterial('Steel01', 1, Fy_cable, E_cable, b0)
        fy = 3550                                 # [N/mm²] Steel Rebar Yield Strength   
        Es = 210e5                                # [N/mm²] Modulus of Elasticity
        ey = fy/Es                                # [mm/mm] Steel Rebar Yield Strain
        fu = 1.1818*fy                            # [N/mm²] Steel Rebar Ultimate Strength
        esu = ey*75.2                             # [mm/mm] Steel Rebar Ultimate Strain
        Esh = (fu - fy)/(esu - ey)
        Bs = Esh / Es
        pinchX = 0.8   # Pinching factor in X direction
        pinchY = 0.5   # Pinching factor in Y direction
        damage1 = 0.0  # Damage due to ductility
        damage2 = 0.0  # Damage due to energy
        beta = 0.1     # Stiffness degradation parameter
        ops.uniaxialMaterial('Hysteretic', 1, fy, ey, fu, esu, 0.2*fu, 1.1*esu, -fy, -ey, -fu, -esu, -0.2*fu, -1.1*esu, pinchX, pinchY, damage1, damage2, beta)
        # INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Hysteretic_Material
 
    # Define truss elements
    for i in range(num_nodes - 1):
        ops.element('corotTruss', i + 1, i + 1, i + 2, A_cable_02, 1)  # Top Cable element
        ops.element('corotTruss', num_nodes + i + 1, num_nodes + i + 1, num_nodes + i + 2, A_cable_01, 1)  # Bottom Cable Element

    # Connect each node of the cable with the corresponding node of the deck using truss elements
    for i in range(num_nodes):
        ops.element('corotTruss', 2 * num_nodes + i + 1, i + 1, num_nodes + i + 1, A_cable_03, 1) # Vertical Cable Element
    
    #mid_node = num_nodes // 2 + 1
    mid_node = int(num_nodes + 0.5 * num_nodes)

    # Dynamic analysis setup
    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.system('UmfPack')
    ops.test('EnergyIncr', MAX_TOLERANCE, MAX_ITERATIONS)
    #ops.integrator('CentralDifference')
    #ops.integrator('Newmark', 0.5, 0.25) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/newmark.html
    alpha = 1;gamma=1.5-alpha; beta=((2-alpha)**2)/4;
    ops.integrator('HHT', alpha, gamma, beta) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/hht.html
    ops.algorithm('ModifiedNewton') # INFO LINK: https://openseespydoc.readthedocs.io/en/stable/src/algorithm.html

    # Define analysis type
    ops.analysis('Transient')
    
    
        
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
    
    # Perform eigenvalue analysis to determine modal periods
    if periodTF == True:
        # Calculate Rayleigh damping factors
        PERIOD_01, PERIOD_02 = S06.RAYLEIGH_DAMPING(2, damping_ratio, 0.6*damping_ratio, 0, 1)
        modal_period = PERIOD_01        # [Second]  
        frequency = 1 / modal_period    # [Hertz]

    if periodTF == False:
        modal_period = 0.0            # [Second]  
        frequency = 0.0               # [Hertz]


    # OUTPUT DATA
    ops.recorder('Node', '-file', f"{SALAR_DIR}DTH_DYN.txt",'-time', '-node', mid_node, '-dof', 1,2, 'disp')# Displacement Time History
    ops.recorder('Node', '-file', f"{SALAR_DIR}VTH_DYN.txt",'-time', '-node', mid_node, '-dof', 1,2, 'vel')  # Velocity Time History
    ops.recorder('Node', '-file', f"{SALAR_DIR}ATH_DYN.txt",'-time', '-node', mid_node, '-dof', 1,2, 'accel') # Acceleration Time History
    ops.recorder('Node', '-file', f"{SALAR_DIR}BTH_DYN_01.txt",'-time', '-node', 1, '-dof', 1,2, 'reaction')# Base Shear Time History
    ops.recorder('Node', '-file', f"{SALAR_DIR}BTH_DYN_02.txt",'-time', '-node', num_nodes, '-dof', 1,2, 'reaction')# Base Shear Time History
    ops.recorder('Node', '-file', f"{SALAR_DIR}BTH_DYN_03.txt",'-time', '-node', num_nodes + 1, '-dof', 1,2, 'reaction')# Base Shear Time History
    ops.recorder('Node', '-file', f"{SALAR_DIR}BTH_DYN_04.txt",'-time', '-node', 2*num_nodes, '-dof', 1,2, 'reaction')# Base Shear

    DISPLACEMENTS = []
    DISP_X , DISP_Y, VELOCITY, ACCELERATION = [], [], [], []
    BASE_REACTION_X, BASE_REACTION_Y = [], []
    PERIOD_MIN, PERIOD_MAX = [], []

    stable = 0
    current_time = 0.0
     # Perform the analysis with increments
    while stable == 0 and current_time < duration:
        OK = ops.analyze(1, dt)
        S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS)
        current_time = ops.getTime()
        DISPLACEMENTS.append([ops.nodeDisp(j + 1) for j in range(num_nodes * 2)])
        DISP_X.append(ops.nodeDisp(mid_node, 1))
        DISP_Y.append(ops.nodeDisp(mid_node, 2))
        VELOCITY.append(ops.nodeVel(mid_node, 2))
        ACCELERATION.append(ops.nodeAccel(mid_node, 2))
        x1 = ops.nodeResponse(1, 1, 6)
        x2 = ops.nodeResponse(num_nodes, 1, 6)
        x3 = ops.nodeResponse(num_nodes + 1, 1, 6)
        x4 = ops.nodeResponse(2*num_nodes, 1, 6)
        y1 = ops.nodeResponse(1, 2, 6)
        y2 = ops.nodeResponse(num_nodes, 2, 6)
        y3 = ops.nodeResponse(num_nodes + 1, 2, 6)
        y4 = ops.nodeResponse(2*num_nodes, 2, 6)
        BASE_REACTION_X.append(x1+x2+x3+x4) # CABLE REACION-X
        BASE_REACTION_Y.append(y1+y2+y3+y4) # CABLE REACION-Y
        # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
        Period_MIN, Period_MAX = S07.EIGENVALUE_ANALYSIS(4, PLOT=True)
        PERIOD_MIN.append(Period_MIN) 
        PERIOD_MAX.append(Period_MAX)
        #print(f'time: {current_time}')
    
    # Get initial and final node coordinates for plotting
    initial_coords = np.array([ops.nodeCoord(i + 1) for i in range(num_nodes * 2)])
    deformed_coords = initial_coords + np.array(DISPLACEMENTS[-1])
    
    # Output all unconstrained node displacements
    for i in range(num_nodes * 2):
        disp = ops.nodeDisp(i + 1)
        print(f"Node {i + 1}: X Displacement = {disp[0]}, Y Displacement = {disp[1]}")
        
    #ops.wipe()
    damping_ratioX = S05.DAMPING_RATIO(DISP_X)
    print(f' Damping Structure Dynamic Analysis Done.') 
    print(f'Structural Damping Ratio in X Dir. : {damping_ratioX:.3f}')
    damping_ratioY = S05.DAMPING_RATIO(DISP_Y) 
    print(f'Structural Damping Ratio in Y Dir. : {damping_ratioY:.3f}')
    """
    if damping == True:
        # Calculate Structure Damping Ratio Based on Lateral Displacement
        damping_ratioX = S05.DAMPING_RATIO(DISP_X)
        print(f' Damping Structure Dynamic Analysis Done.') 
        print(f'Structural Damping Ratio in X Dir. : {damping_ratioX:.3f}')
        damping_ratioY = S05.DAMPING_RATIO(DISP_Y) 
        print(f'Structural Damping Ratio in Y Dir. : {damping_ratioY:.3f}')
    """    
    #print("Period: ", modal_period)
    #print("Frequency: ", frequency)

    DATA = (initial_coords, deformed_coords, DISPLACEMENTS, VELOCITY,
            ACCELERATION, DISP_X, DISP_Y, BASE_REACTION_X, BASE_REACTION_Y, damping_ratioX, damping_ratioY)
    
    return DATA
#%%------------------------------------------------------------------
def plot_time_history(DISP_X_undamped, DISP_X_damped, DISP_Y_undamped, DISP_Y_damped, 
                      VELOCITY_undamped, VELOCITY_damped, ACCELERATION_undamped, ACCELERATION_damped, 
                      BASE_REACTION_X_undamped, BASE_REACTION_X_damped, BASE_REACTION_Y_undamped, BASE_REACTION_Y_damped):
    import matplotlib.pyplot as plt
    import numpy as np
    # Assuming you have 'time' array representing the time steps
    time = [i * dt for i in range(len(DISPLACEMENTS_undamped))]
    plt.figure(figsize=(14, 20))

    # Plot Displacements in X direction
    plt.subplot(6, 1, 1)
    P1 = DISP_X_undamped
    P2 = DISP_X_damped
    P11 = np.max(np.abs(P1))
    P22 = np.max(np.abs(P2))
    P3 = P11 / P22
    plt.plot(time, P1, color='black', label=f'Undamped X: {P11:.5e}')
    plt.plot(time, P2, color='red', label=f'Damped X: {P22:.5e}')
    plt.title(f'Displacement Time History (X direction) - Amplification Factor: {P3:.3f}')
    plt.xlabel('Time [s]')
    plt.ylabel('Displacement [m] (X direction)')
    plt.legend()

    # Plot Displacements in Y direction
    plt.subplot(6, 1, 2)
    P1 = DISP_Y_undamped
    P2 = DISP_Y_damped
    P11 = np.max(np.abs(P1))
    P22 = np.max(np.abs(P2))
    P3 = P11 / P22
    plt.plot(time, P1, color='black', label=f'Undamped Y: {P11:.5e}')
    plt.plot(time, P2, color='red', label=f'Damped Y: {P22:.5e}')
    plt.title(f'Displacement Time History (Y direction) - Amplification Factor: {P3:.3f}')
    plt.xlabel('Time [s]')
    plt.ylabel('Displacement [m] (Y direction)')
    plt.legend()

    # Plot Velocities
    plt.subplot(6, 1, 3)
    P1 = VELOCITY_undamped
    P2 = VELOCITY_damped
    P11 = np.max(np.abs(P1))
    P22 = np.max(np.abs(P2))
    P3 = P11 / P22
    plt.plot(time, P1, color='black', label=f'Undamped: {P11:.5e}')
    plt.plot(time, P2, color='red', label=f'Damped: {P22:.5e}')
    plt.title(f'Velocity Time History - Amplification Factor: {P3:.3f}')
    plt.xlabel('Time [s]')
    plt.ylabel('Velocity [m/s] (Y direction)')
    plt.legend()

    # Plot Accelerations
    plt.subplot(6, 1, 4)
    P1 = ACCELERATION_undamped
    P2 = ACCELERATION_damped
    P11 = np.max(np.abs(P1))
    P22 = np.max(np.abs(P2))
    P3 = P11 / P22
    plt.plot(time, P1, color='black', label=f'Undamped: {P11:.5e}')
    plt.plot(time, P2, color='red', label=f'Damped: {P22:.5e}')
    plt.title(f'Acceleration Time History - Amplification Factor: {P3:.3f}')
    plt.xlabel('Time [s]')
    plt.ylabel('Acceleration [m/s²] (Y direction)')
    plt.legend()

    # Plot Base Reactions in X direction
    plt.subplot(6, 1, 5)
    P1 = BASE_REACTION_X_undamped
    P2 = BASE_REACTION_X_damped
    P11 = np.max(np.abs(P1))
    P22 = np.max(np.abs(P2))
    P3 = P11 / P22
    plt.plot(time, P1, color='black', label=f'Undamped: {P11:.5e}')
    plt.plot(time, P2, color='red', label=f'Damped: {P22:.5e}')
    plt.title(f'Base Reaction Time History (X direction) - Amplification Factor: {P3:.3f}')
    plt.xlabel('Time [s]')
    plt.ylabel('Base Reaction [N] (X direction)')
    plt.legend()

    # Plot Base Reactions in Y direction
    plt.subplot(6, 1, 6)
    P1 = BASE_REACTION_Y_undamped
    P2 = BASE_REACTION_Y_damped
    P11 = np.max(np.abs(P1))
    P22 = np.max(np.abs(P2))
    P3 = P11 / P22
    plt.plot(time, P1, color='black', label=f'Undamped: {P11:.5e}')
    plt.plot(time, P2, color='red', label=f'Damped: {P22:.5e}')
    plt.title(f'Base Reaction Time History (Y direction) - Amplification Factor: {P3:.3f}')
    plt.xlabel('Time [s]')
    plt.ylabel('Base Reaction [N] (Y direction)')
    plt.legend()

    plt.tight_layout()
    plt.show()
#%%------------------------------------------------------------------
# -------------------
#   DYNAMIC ANALYSIS
# -------------------
""" OMITTED, BECAUSE THEY ARE AS SAME AS TOP ANALYSIS PARAMETERS
# Parameters for the analysis
L = 210000.0        # [mm] Bridge span lentgh
H1 = 25000.0        # [mm] Height of Top Cable
arc_depth = 20000.0 # [mm]
E_cable = 210e5     # [N/mm^2] Modulus of elasticity Cable
Cable_Dia_01 = 145600 # [mm] Horizontal Longitudinal Cable Diameter Bottom
Cable_Dia_02 = 300  # [mm] Horizontal Longitudinal Cable Diameter Top
Cable_Dia_03 = 50   # [mm] Vertical Longitudinal Cable Diameter Top 
num_nodes = 22      # Cable Arc Number of nodes
"""
TOTAL_MASS = 50_000_000.0 # [kg] Total Mass of Structure
damping_ratio = 0.01   # Damping ratio
duration = 50.0        # [s] Duration of the analysis in seconds
dt = 0.01              # Time step in seconds

MAX_ITERATIONS = 10000  # Maximum number of iterations
MAX_TOLERANCE = 1.0e-6  # Tolerance for convergence

import time
starttime = time.process_time()

# Run the undamped analysis
damping = False
LINEAR = True # False: Cable Nonlinear Materials Properties
periodTF = False
results_undamped = DYNAMIC_ANALYSIS(damping, damping_ratio, LINEAR, periodTF, duration, dt, TOTAL_MASS, L, H1, arc_depth, E_cable, Cable_Dia_01, Cable_Dia_02, Cable_Dia_03, num_nodes, MAX_ITERATIONS, MAX_TOLERANCE)
initial_coords, deformed_coords, DISPLACEMENTS_undamped, VELOCITY_undamped, ACCELERATION_undamped, DISP_X_undamped, DISP_Y_undamped, BASE_REACTION_X_undamped, BASE_REACTION_Y_undamped, damping_ratioX_undamped, damping_ratioY_undamped = results_undamped

# Run the damped analysis
damping = True
LINEAR = True # False: Cable Nonlinear Materials Properties
periodTF = False
results_damped = DYNAMIC_ANALYSIS(damping, damping_ratio, LINEAR, periodTF, duration, dt, TOTAL_MASS, L, H1, arc_depth, E_cable, Cable_Dia_01, Cable_Dia_02, Cable_Dia_03, num_nodes, MAX_ITERATIONS, MAX_TOLERANCE)
initial_coords, deformed_coords, DISPLACEMENTS_damped, VELOCITY_damped, ACCELERATION_damped, DISP_X_damped, DISP_Y_damped, BASE_REACTION_X_damped, BASE_REACTION_Y_damped, damping_ratioX_damped, damping_ratioY_damped = results_damped

totaltime = time.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n')
#%%------------------------------------------------------------------
plot_time_history(DISP_X_undamped, DISP_X_damped, DISP_Y_undamped, DISP_Y_damped, VELOCITY_undamped, VELOCITY_damped, ACCELERATION_undamped, ACCELERATION_damped, BASE_REACTION_X_undamped, BASE_REACTION_X_damped, BASE_REACTION_Y_undamped, BASE_REACTION_Y_damped)
#%%------------------------------------------------------------------
def PLOT_2D(X1, Y1, X2, Y2, XLABEL, YLABEL, TITLE):
    plt.figure(figsize=(10, 6))
    plt.plot(X1, Y1, label='Undamped', color='black')
    plt.plot(X2, Y2, label='Damped', color='red')
    plt.xlabel(XLABEL)
    plt.ylabel(YLABEL)
    plt.title(TITLE)
    plt.grid(True)
    #plt.semilogy()
    plt.legend()
    plt.show()
#%%------------------------------------------------------------------
### BASE REACTION & DISPALCEMENT IN X:
X1 =  DISP_X_undamped
X2 =  DISP_X_damped
Y1 = BASE_REACTION_X_undamped
Y2 = BASE_REACTION_X_damped
XLABEL = 'Displacement X'
YLABEL = 'Base-Reaction X'
TITLE = 'Base-Reaction and Displacement X'

PLOT_2D(X1, Y1, X2, Y2, XLABEL, YLABEL, TITLE)
#%%------------------------------------------------------------------
### BASE REACTION & DISPALCEMENT IN Y:
X1 =  DISP_Y_undamped
X2 =  DISP_Y_damped
Y1 = BASE_REACTION_Y_undamped
Y2 = BASE_REACTION_Y_damped
XLABEL = 'Displacement Y'
YLABEL = 'Base-Reaction Y'
TITLE = 'Base-Reaction and Displacement Y'

PLOT_2D(X1, Y1, X2, Y2, XLABEL, YLABEL, TITLE)
#%%------------------------------------------------------------------
# Plot damped results
plot_shapes(initial_coords, DISPLACEMENTS_damped)
#%%------------------------------------------------------------------
"""
Plotting a histogram for the results of dynamic analysis of a structure, such as displacements, accelerations, and support reactions, helps analyze the statistical distribution and trends in the data. The following insights can be drawn from these histograms:

Range Identification The histogram shows the range within which displacements, accelerations, or support reactions vary. This is crucial for assessing the stability of the structure and determining the allowable design limits.
Central Tendency From the histogram, measures like mean, median, or mode can be extracted. This helps identify the most likely values and provides an overview of the dominant behavior of the structure.
Data Dispersion The histogram illustrates the spread of the data (e.g., through its width and shape). A wide histogram may indicate nonlinear behavior or unusual effects in the structure.
Anomaly Detection Unexpected peaks or asymmetric distributions may indicate abnormal structural behavior (e.g., weaknesses or inelastic behavior). It can also reveal specific phenomena like resonances or amplification effects.
Probability of Specific Values The histogram's shape can indicate the probability distribution of the data (e.g., normal, uniform, or Poisson). This is useful for simulations and predicting future structural behavior.
Behavior Under Specific Conditions Analyzing histograms of acceleration or displacement can reveal dominant modes of the structure and their impact on overall response. These insights are helpful in designing vibration control systems (e.g., dampers).
Impact of External Forces or Loading An unusual histogram shape might suggest the influence of extraordinary dynamic loads (e.g., severe earthquakes or impact forces). Conclusion: Histograms are a simple yet effective tool for analyzing dynamic data. They help engineers and designers:
Understand the overall behavior of a structure. Assess the likelihood of extreme or unexpected responses. Define appropriate design or retrofitting criteria. For deeper analysis, histograms can be combined with other statistical tools, such as standard deviation or spectral analysis.
"""
#%%------------------------------------------------------------------
HISROGRAM_BOXPLOT(DISP_X_undamped, HISTO_COLOR='black', LABEL='Displacement X (Undamped)')
HISROGRAM_BOXPLOT(DISP_X_damped, HISTO_COLOR='grey', LABEL='Displacement X (Damped)')
#HISTOGRAM_BOXPLOT_PLOTLY(DISP_X_undamped, XLABEL='Displacement X (Undamped)', TITLE='Displacement X (Undamped)', COLOR='black')
#HISTOGRAM_BOXPLOT_PLOTLY(DISP_X_damped, XLABEL='Displacement X (Damped)', TITLE='Displacement X (Damped)', COLOR='grey')
#%%------------------------------------------------------------------
HISROGRAM_BOXPLOT(DISP_Y_undamped, HISTO_COLOR='orange', LABEL='Displacement Y (Undamped)')
HISROGRAM_BOXPLOT(DISP_Y_damped, HISTO_COLOR='brown', LABEL='Displacement Y (Damped)')
#HISTOGRAM_BOXPLOT_PLOTLY(DISP_Y_undamped, XLABEL='Displacement Y (Undamped)', TITLE='Displacement Y (Undamped)', COLOR='orange')
#HISTOGRAM_BOXPLOT_PLOTLY(DISP_Y_damped, XLABEL='Displacement Y (Damped)', TITLE='Displacement Y (Damped)', COLOR='brown')
#%%------------------------------------------------------------------
HISROGRAM_BOXPLOT(VELOCITY_undamped, HISTO_COLOR='red', LABEL='Velocity (Undamped)')
HISROGRAM_BOXPLOT(VELOCITY_damped, HISTO_COLOR='pink', LABEL='Velocity (Damped)')
#HISTOGRAM_BOXPLOT_PLOTLY(VELOCITY_undamped, XLABEL='Velocity (Undamped)', TITLE='Velocity (Undamped)', COLOR='red')
#HISTOGRAM_BOXPLOT_PLOTLY(VELOCITY_damped, XLABEL='Velocity (Damped)', TITLE='Velocity (Damped)', COLOR='pink')
#%%------------------------------------------------------------------
HISROGRAM_BOXPLOT(ACCELERATION_undamped, HISTO_COLOR='gold', LABEL='Acceleration (Undamped)')
HISROGRAM_BOXPLOT(ACCELERATION_damped, HISTO_COLOR='brown', LABEL='Acceleration (Damped)')
#HISTOGRAM_BOXPLOT_PLOTLY(ACCELERATION_undamped, XLABEL='Acceleration (Undamped)', TITLE='Acceleration (Undamped)', COLOR='gold')
#HISTOGRAM_BOXPLOT_PLOTLY(ACCELERATION_damped, XLABEL='Acceleration (Damped)', TITLE='Acceleration (Damped)', COLOR='brown')
#%%------------------------------------------------------------------
HISROGRAM_BOXPLOT(BASE_REACTION_X_undamped, HISTO_COLOR='blue', LABEL='Base Reaction X (Undamped)')
HISROGRAM_BOXPLOT(BASE_REACTION_X_damped, HISTO_COLOR='cyan', LABEL='Base Reaction X (Damped)')
#HISTOGRAM_BOXPLOT_PLOTLY(BASE_REACTION_X_undamped, XLABEL='Base Reaction X (Undamped)', TITLE='Base Reaction X (Undamped)', COLOR='blue')
#HISTOGRAM_BOXPLOT_PLOTLY(BASE_REACTION_X_damped, XLABEL='Base Reaction X (Damped)', TITLE='Base Reaction X (Damped)', COLOR='cyan')
#%%------------------------------------------------------------------
HISROGRAM_BOXPLOT(BASE_REACTION_Y_undamped, HISTO_COLOR='green', LABEL='Base Reaction Y (Undamped)')
HISROGRAM_BOXPLOT(BASE_REACTION_Y_damped, HISTO_COLOR='lime', LABEL='Base Reaction Y (Damped)')
#HISTOGRAM_BOXPLOT_PLOTLY(BASE_REACTION_Y_undamped, XLABEL='Base Reaction Y (Undamped)', TITLE='Base Reaction Y (Undamped)', COLOR='green')
#HISTOGRAM_BOXPLOT_PLOTLY(BASE_REACTION_Y_damped, XLABEL='Base Reaction Y (Damped)', TITLE='Base Reaction Y (Damped)', COLOR='lime')
#%%------------------------------------------------------------------