#          #####################################################################################
#          #                                 IN THE NAME OF ALLAH                              #
#          # PROGRESSIVE COLLAPSE ANALYSIS OF CONCRETE 2-STORY FRAME WITH DISPLACEMENT CONTROL #
#          #-----------------------------------------------------------------------------------#
#          #              THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)           #
#          #                       EMAIL: salar.d.ghashghaei@gmail.com                         #
#          #####################################################################################
"""
Progressive collapse of reinforced concrete frames occurs when a local failure—due to accidental
 actions such as impact, explosion or fire—triggers a chain reaction of element removals, leading
 to partial or total structural loss. Advanced assessment hinges on capturing nonlinear material
 behavior, geometric effects, and load‐redistribution mechanisms that dictate whether alternative
 load paths can sustain the imposed demands.

[1] Modeling Philosophy:
- Fiber based sections discretize concrete and steel across the cross-section, enabling accurate
 stress–strain representation under combined axial, bending and shear demands. Cover, core concrete,
 and rebar layouts are modeled with uniaxial constitutive laws that include confinement, cracking,
 strain hardening and ultimate strain limits.
- Nonlinear beam–column elements employ Gauss integration points along member length, paired with
 corotational kinematics to account for large displacements and P-Δ effects in a fully consistent 2D formulation.

[2] Analysis Strategy:
- Alternate load‐path method: deliberately remove one or more columns (or beams) after applying
 gravity loads, then trace the static response under incremental displacement control at a critical
 location. The structural response captures bending yielding, shear failure, catenary action and
 eventual loss of load‐bearing capacity.
- Pushover framework: displacement control at a predefined “attack” node (e.g., mid-height of a key column)
 simulates the increasing drift demands after element removal. Reaction forces at the base yield
 a capacity curve relating force vs. displacement, from which reserve strength and ductility can be assessed.

[3] Key Response Mechanisms:
- Flexural yielding and plastic hinge formation in adjacent beams and columns allow moment redistribution.
 Hinge rotation capacity depends on reinforcement ratio, concrete confinement and strain‐hardening characteristics of steel.
- P-Δ instability magnifies demands when large drifts develop; corotational transforms ensure equilibrium
 accounts for geometric nonlinearity.
- Catenary action engages once flexural capacity is exhausted and members deform significantly, mobilizing
 tensile forces in reinforcement. Accurate modeling of ultimate tendon strain (eult) is critical to predict post-peak response.
- Shear failure remains brittle; its prevention through detailing (stirrups, confinement) is vital to allow
 ductile mechanisms to develop.

[4] Collapse Criteria and Robustness:
- Vertical and lateral drift limits define collapse thresholds. Exceeding the ultimate drift at a control
 node triggers element deletion, simulating fracture or buckling.
- Progressive removal tests on different locations probe system robustness, verifying that the structure
 retains sufficient redundancy and alternative load paths.

[5] Practical Implications:
- Design against progressive collapse requires enforcing continuity (tie forces), detailing for ductility
 (strong-column–weak-beam hierarchy), and redundancy (multiple load paths).
- Nonlinear analyses—both static pushover and dynamic removal simulations—inform code provisions
 (e.g., UFC 4-023-03, GSA Guidelines) by quantifying reserve strength margins and post-failure behavior.

In workflow, the combination of fiber section definitions, corotational beam–column elements,
 displacement-controlled pushover and element removal routines captures the critical phases of
 progressive collapse: initial yielding, redistribution, catenary action and final instability.
 Interpreting capacity curves and post-peak degradation provides insights on member detailing and
 overall frame robustness under accidental collapse scenarios.

 I Really appreciate  Dr.Barham Ali for his best constructive suggestions and sending some real damaged concrete frame structure pictures
"""
#%%-----------------------------------------------------------------------------------
# wikipedia: Progressive collapse
'https://en.wikipedia.org/wiki/Progressive_collapse'
# BOOK: Progressive Collapse Resilience of Concrete Structures: Mechanisms, Simulations and Experiments
'https://link.springer.com/book/10.1007/978-981-99-0772-4'
# REPORT: Computational Modeling of Progressive Collapse in Reinforced Concrete Frame Structures
'https://peer.berkeley.edu/sites/default/files/webpeer710_mohamed_m._talaat_khalid_m._mosalam.pdf'
# PAPER: The Performance of Resistance Progressive Collapse Analysis for High-Rise Frame-Shear Structure Based on OpenSees
'https://onlinelibrary.wiley.com/doi/10.1155/2017/3518232'
# PAPER: Benchmark Numerical Model for Progressive Collapse Analysis of RC Beam-Column Sub-Assemblages
'https://www.mdpi.com/2075-5309/12/2/122?type=check_update&version=1'
# PAPER: A computationally efficient numerical model for progressive collapse analysis of reinforced concrete structures
'https://journals.sagepub.com/doi/10.1177/2041419619854768?icid=int.sj-full-text.similar-articles.2'
# PAPER: The Performance of Resistance Progressive Collapse Analysis for High-Rise Frame-Shear Structure Based on OpenSees
'https://onlinelibrary.wiley.com/doi/10.1155/2017/3518232'
# PAPER: Refined dynamic progressive collapse analysis of RC structures
'https://www.researchgate.net/publication/321948788_Refined_dynamic_progressive_collapse_analysis_of_RC_structures'
#%%-----------------------------------------------------------------------------------
#import the os module
import os
import time
import openseespy.opensees as ops
import opsvis as opsv
import numpy as np
import matplotlib.pyplot as plt
import CONCRETE_SECTION_FUN as S02
import ANALYSIS_FUNCTION as S03
import PLOT_2D as S04
#%%-----------------------------------------------------------------------------------
# Load the image
def PLOT_IMAGE(image):
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg
    image = mpimg.imread(image_path)

    # Display the image
    plt.figure(figsize=(12, 8))
    plt.imshow(image)
    plt.axis('off')  # Hide axes
    plt.show()
    
image_path = 'OPENSEES_PROGRESSIVE_COLLAPSE.png'    
PLOT_IMAGE(image_path)

#%%-----------------------------------------------------------------------------------
image_path = 'OPENSEES_PROGRESSIVE_COLLAPSE_CONCRETE_FRAME.png'    
PLOT_IMAGE(image_path)
#%%-----------------------------------------------------------------------------------
def CURRENT_TIME():
    import time
    t = time.localtime()
    current_time = time.strftime("%H:%M:%S", t)
    print(f"Current time (HH:MM:SS): {current_time}\n\n")
#%%-----------------------------------------------------------------------------------
import matplotlib.pyplot as plt

# Define node coordinates (in mm)
H1 = 3000 # [mm] 1st floor length
H2 = 3000 # [mm] 2nd floor length
L1 = 5000 # [mm] 1st bay length
L2 = 7000 # [mm] 2nd bay length

# Node coordinates (in mm)
node_coords = [
    (0, 0), (0, H1), (0, H1+H2),  # Left column (base, 1st floor, 2nd floor)
    (L1, 0), (L1, H1), (L1, H1+H2),  # Middle column (base, 1st floor, 2nd floor)
    (L1+L2, 0), (L1+L2, H1), (L1+L2, H1+H2)  # Right column (base, 1st floor, 2nd floor)
]

# Define nonlinear beam-column elements for columns and beams
elementsZ = [
    (0, 1), (1, 2),  # Left column
    (3, 4), (4, 5),  # Middle column
    (6, 7), (7, 8),  # Right column
    (1, 4), (4, 7),  # Bottom beam
    (2, 5), (5, 8)   # Top beam
]

# Define fixed nodes (base nodes)
fixed_nodes = [0, 3, 6]

# Extract node coordinates
x_coords, y_coords = zip(*node_coords)

# Plot the 2-story, 2-bay structure
plt.figure(figsize=(8, 6))
for element in elementsZ:
    x_vals = [x_coords[element[0]], x_coords[element[1]]]
    y_vals = [y_coords[element[0]], y_coords[element[1]]]
    plt.plot(x_vals, y_vals, 'b-o', markersize=8)

# Annotate node numbers
for i, (x, y) in enumerate(node_coords):
    plt.text(x, y, f'{i+1}', fontsize=12, ha='right')

# Plot red rectangles to show fixed constraints
for node in fixed_nodes:
    plt.plot(x_coords[node], y_coords[node], 's', color='red', markersize=15)

plt.xlabel('X (mm)')
plt.ylabel('Y (mm)')
plt.title('2-Story, 2-Bay Structure with Fixed Constraints')
plt.grid(True)
plt.axis('equal')
plt.show()


#%%-----------------------------------------------------------------------------------
nFibCoverZ, nFibCoverY = 1 , 1
nFibCoreZ, nFibCoreY = 1, 1
FS01, FS02 = S02.CONCRETE_SECTION(600, 600, 600, 400, 50, 25,nFibCoverZ, nFibCoverY, nFibCoreZ, nFibCoreY, PLOT=1)
#%%-----------------------------------------------------------------------------------


# Initialize the model
ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 3)

# Define nonlinear beam-column elements for columns and beams
elements = [
    (1, 2), (2, 3), # Left column
    (4, 5), (5, 6), # Middle column
    (7, 8), (8, 9), # Right column
    (2, 5), (5, 8), # Bottom beam
    (3, 6), (6, 9)  # Top beam
]

# Define nodes
for i, coord in enumerate(node_coords):
    ops.node(i + 1, *coord)

# Fix base nodes
ops.fix(1, 1, 1, 1)
ops.fix(4, 1, 1, 1)
ops.fix(7, 1, 1, 1)

#ops.uniaxialMaterial('ReinforcingSteel', 1, fy_steel, fu_steel, E_steel, Esh, esh, eult)
# LINK: https://opensees.berkeley.edu/wiki/index.php?title=Reinforcing_Steel_Material

# Define fiber section for I-section
Bcol, Hcol, Bbeam, Hbeam = 600, 600, 600, 400; # [mm] Column & Beam Section Diamenstion Properties
COVER = 50        # [mm] Concrete Cover
REBAR_DIA = 25    # [mm] Steel Rebar Diameter

MAX_ITERAIONS = 5000
TOLERANCE = 1.0e-12

# Concrete Sections for Beams and Columns
nFibCoverZ, nFibCoverY, nFibCoreZ, nFibCoreY = 3, 120, 3, 120
SECTION01, SECTION02 = S02.CONCRETE_SECTION(Bcol, Hcol, Bbeam, Hbeam, COVER, REBAR_DIA,
                                             nFibCoverZ, nFibCoverY, nFibCoreZ, nFibCoreY, PLOT=0)
opsv.fib_sec_list_to_cmds(SECTION01) # COLUMNS SECTION
opsv.fib_sec_list_to_cmds(SECTION02) # BEAMS SECTION

UDL = -0.001   # [N/mm] Uniform Distributed Loads
PY = -12000    # [N] Verictal Constant Load in Node [5]

Max_Disp = -2500        # [mm] Maximum Vertical Displacement
Collapse_Disp = -2500   # [mm] Absolute Value Collapse Vertical Displacement
disp_incr = -0.5        # [mm] Displacement Increment

# Define geometric transformation
# Linear:
# Small displacement assumptions in local to basic transformation
# Linear transformation of forces and displacements
# ops.geomTransf('Linear', 1)

# P-Delta:
# Small displacement assumption transformation of displacements
# Account for transverse displacement of axial load in equilibrium relationship
# ops.geomTransf('PDelta', 1)

# Corotational:
# Fully nonlinear transformation of displacements and force
# Exact in 2D but some approximations in 3D
ops.geomTransf('Corotational', 1) 

for i, (iNode, jNode) in enumerate(elements):
    #                                 $eleTag $iNode $jNode $numIntgrPts $secTag $transfTag
    if i <= 5: # COLUMNS
        ops.element('nonlinearBeamColumn', i + 1, iNode, jNode, 5, 1, 1)
    else:      # BEAMS
        ops.element('nonlinearBeamColumn', i + 1, iNode, jNode, 5, 2, 1)
    

                
# Apply load pattern for pushover analysis on the middle column (node 5)
ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)
#ops.load(2, 0.0, PY, 0.0)  # Vertical load
ops.load(5, 0.0, PY, 0.0)  # Vertical load
#ops.load(8, 0.0, PY, 0.0)  # Vertical load
# Uniform Distributed Load
for i in range(6, 10):
        # mag of uniformily distributed ref load acting in local y direction of element
        ops.eleLoad('-ele', i + 1,'-type', '-beamUniform', UDL, 0.0)
    
#ops.recorder('Collapse', '-ele', 4, '-node', 5, '-file', 'Collapse.txt')
# Define analysis parameters
ops.system('BandGeneral')
ops.numberer('RCM')
ops.constraints('Plain')
ops.test('EnergyIncr', TOLERANCE, MAX_ITERAIONS)
ops.algorithm('ModifiedNewton')
ops.integrator('DisplacementControl', 5, 2, disp_incr)
ops.analysis('Static')

print('Model Done.')

# Perform pushover analysis
n_steps = int(np.abs(Max_Disp / disp_incr)) # Analysis Steps
displacementsX = []
displacementsY = []
rotations = []
forcesH = []
forcesV = []
forcesM = []
S_col01, S_col02, S_col03 = [], [], []
A_col01, A_col02, A_col03 = [], [], []
M_col01, M_col02, M_col03 = [], [], []
KA, KS, KI = [], [], []

# PLOT CURRENT TIME
CURRENT_TIME()
delete_element = False

for step in range(n_steps):
    #print(step + 1)
    OK = ops.analyze(1)
    S03.ANALYSIS(OK, 1,TOLERANCE, MAX_ITERAIONS)
    dispX = ops.nodeDisp(5, 1) # HORIZENTAL DISPLACEMENT FOR NODE 5
    dispY = ops.nodeDisp(5, 2) # VERTICAL DISPLACEMENT FOR
    rotat = ops.nodeDisp(5, 3) # ROTATION IN NODE 5
    ops.reactions()
    if abs(dispY) < abs(Collapse_Disp): # DISPLACEMENT CONTROL FOR COLUMN 3 UNTIL COLUMN IS NOT DELETED
        S = ops.nodeReaction(1, 1) + ops.nodeReaction(4, 1) + ops.nodeReaction(7, 1) # SHEAR BASE REACTION
        A = ops.nodeReaction(1, 2) + ops.nodeReaction(4, 2) + ops.nodeReaction(7, 2) # AXIAL BASE REACTION
        M = ops.nodeReaction(1, 3) + ops.nodeReaction(4, 3) + ops.nodeReaction(7, 3) # MOMENT BASE REACTION
    if abs(dispY) == abs(Collapse_Disp) and not delete_element: # DISPLACEMENT CONTROL FOR COLUMN 3 AND DELETE ELEMENT
        print(f"Displacement exceeds {Collapse_Disp} mm. Removing element 5 at step {step + 1}.") 
        ops.remove('element', 3) # REMOVE ELEMENT
        ops.remove('sp', 4) # REMOVE FIX SUPPORT
        delete_element = True 
        S = ops.nodeReaction(1, 1) + ops.nodeReaction(7, 1) # SHEAR BASE REACTION
        A = ops.nodeReaction(1, 2) + ops.nodeReaction(7, 2) # AXIAL BASE REACTION
        M = ops.nodeReaction(1, 3) + ops.nodeReaction(7, 3) # MOMENT BASE REACTION
    if abs(dispY) > abs(Collapse_Disp): # DISPLACEMENT CONTROL FOR COLUMN 3 UNTIL COLUMN IS NOT DELETED
        S = ops.nodeReaction(1, 1) + ops.nodeReaction(7, 1) # SHEAR BASE REACTION
        A = ops.nodeReaction(1, 2) + ops.nodeReaction(7, 2) # AXIAL BASE REACTION
        M = ops.nodeReaction(1, 3) + ops.nodeReaction(7, 3) # MOMENT BASE REACTION
        
    #force =  ops.eleResponse(1, 'force')[1] + ops.eleResponse(3, 'force')[1] + ops.eleResponse(5, 'force')[1]    
    #print(force)
    displacementsX.append(np.abs(dispX)) # DISPLACEMENT IN X FOR  NODE [5]
    displacementsY.append(np.abs(dispY)) # DISPLACEMENT IN Y FOR  NODE [5]
    rotations.append(rotat)              # ROTAION NODE[5]
    forcesH.append(S)                    # BASE-SHEAR REACTION 
    forcesV.append(A)                    # BASE-AXIAL REACTION 
    forcesM.append(M)                    # BASE-MOMENT REACTION 
    KS.append(np.abs(S)/np.abs(dispX))   # LATERAL STIFFNESS IN X
    KA.append(np.abs(A)/np.abs(dispY))   # LATERAL STIFFNESS IN Y
    KI.append(np.abs(M)/np.abs(rotat))   # ROTATIONAL STIFFNESS IN Z
    # COLUMN AXIAL FORCE
    A_col01.append(ops.eleResponse(1, 'force')[1])
    if abs(dispY) < abs(Collapse_Disp):
        A_col02.append(ops.eleResponse(3, 'force')[1])
    else:
        A_col02.append(0)  
    A_col03.append(ops.eleResponse(5, 'force')[1])
    # COLUMN SHEAR FORCE
    S_col01.append(ops.eleResponse(1, 'force')[0])
    if abs(dispY) < abs(Collapse_Disp):
        S_col02.append(ops.eleResponse(3, 'force')[0])
    else:
        S_col02.append(0)     
    S_col03.append(ops.eleResponse(5, 'force')[0])
    # COLUMN MOMENT FORCE
    M_col01.append(ops.eleResponse(1, 'force')[2])
    if abs(dispY) < abs(Collapse_Disp):
        M_col02.append(ops.eleResponse(3, 'force')[2])
    else:
        M_col02.append(0)    
    M_col03.append(ops.eleResponse(5, 'force')[2])

    print(step + 1, dispY, A)


#ops.wipe() -> IF IT WILL BE ACTIVE, SO WE CAN'T SEE FRAME SHAPE

print('Analysis is Done.')
# PLOT CURRENT TIME
CURRENT_TIME()
#%%-----------------------------------------------------------------------------------
# Plot the Pushover urve
def PLOT_2D(X, Y, XLABEL, YLABEL, COLOR):
    plt.figure(figsize=(10, 6))
    plt.plot(X, Y,  color=COLOR, linewidth=4)
    plt.xlabel(XLABEL)
    plt.ylabel(YLABEL)
    plt.title('Pushover Analysis of Structure')
    plt.grid(True)
    #plt.semilogy()
    plt.show()
#%%-----------------------------------------------------------------------------------    
X = displacementsX 
Y = forcesH
XLABEL = 'Horizontal Displacement (mm) NODE[5]'
YLABEL = 'Shear Base Reaction (N)'
COLOR = 'blue'
PLOT_2D(X, Y, XLABEL, YLABEL, COLOR)  
#%%-----------------------------------------------------------------------------------
X = displacementsY
Y = forcesV
XLABEL = 'Vertical Displacement (mm) NODE[5]'
YLABEL = 'Axial Base Reaction (N)'
COLOR = 'green'
PLOT_2D(X, Y, XLABEL, YLABEL, COLOR)  
#%%-----------------------------------------------------------------------------------
X = rotations 
Y = forcesM
XLABEL = 'Rotation (rad) NODE[5]'
YLABEL = 'Moment Base Reaction (N.mm)'
COLOR = 'purple'
PLOT_2D(X, Y, XLABEL, YLABEL, COLOR)
#%%-----------------------------------------------------------------------------------
X = displacementsX 
Y = displacementsY
XLABEL = 'Horizontal Displacement (mm) in X Dir. NODE[5]'
YLABEL = 'Vertical Displacement (mm) in Y Dir. NODE[5]'
COLOR = 'brown'
PLOT_2D(X, Y, XLABEL, YLABEL, COLOR)
#%%-----------------------------------------------------------------------------------
plt.figure(figsize=(10, 6))
plt.plot(M_col01, A_col01, color='black', linewidth=4)
plt.plot(M_col02, A_col02, color='red', linewidth=4)
plt.plot(M_col03, A_col03, color='green', linewidth=4)
plt.xlabel('MOMENT FORCE')
plt.ylabel('AXIAL FORCE')
plt.title('Column Force - Pushover Analysis of Structure')
plt.grid(True)
plt.legend(['COL 1', 'COL 2', 'COL 3'])
#plt.semilogy()
plt.show()
#%%-----------------------------------------------------------------------------------
plt.figure(figsize=(10, 6))
plt.plot(S_col01, A_col01, color='black', linewidth=4)
plt.plot(S_col02, A_col02, color='red', linewidth=4)
plt.plot(S_col03, A_col03, color='green', linewidth=4)
plt.xlabel('SHEAR FORCE')
plt.ylabel('AXIAL FORCE')
plt.title('Column Force - Pushover Analysis of Structure')
plt.grid(True)
plt.legend(['COL 1', 'COL 2', 'COL 3'])
#plt.semilogy()
plt.show()
#%%-----------------------------------------------------------------------------------
plt.figure(figsize=(10, 6))
plt.plot(S_col01, M_col01, color='black', linewidth=4)
plt.plot(S_col02, M_col02, color='red', linewidth=4)
plt.plot(S_col03, M_col03, color='green', linewidth=4)
plt.xlabel('SHEAR FORCE')
plt.ylabel('MOMENT FORCE')
plt.title('Column Force - Pushover Analysis of Structure')
plt.grid(True)
plt.legend(['COL 1', 'COL 2', 'COL 3'])
#plt.semilogy()
plt.show()
#%%-----------------------------------------------------------------------------------
plt.figure(figsize=(12, 8))
#plt.plot(KI, KS, color='black', linewidth=2)
#plt.plot(KI02, KS02, color='grey', linestyle='--', linewidth=2)
plt.scatter(KI, KS, color='black', linewidth=2)
plt.title('ROTATIONAL STIFFNESS-LATERAL STIFFNESS IN X DIR. FOR STORY-1 DIAGRAM')
plt.ylabel('Rotational Stiffness [N.mm/Rad]')
plt.xlabel('Lateral Stiffness in X [N/mm]')
plt.semilogx()
plt.semilogy()
plt.grid()
plt.show()
#%%-----------------------------------------------------------------------------------
plt.figure(figsize=(12, 8))
#plt.plot(KI, KA, color='black', linewidth=2)
#plt.plot(KI02, KA02, color='grey', linestyle='--', linewidth=2)
plt.scatter(KI, KA, color='black', linewidth=2)
plt.title('ROTATIONAL STIFFNESS-LATERAL STIFFNESS IN Y DIR. FOR STORY-1 DIAGRAM')
plt.ylabel('Rotational Stiffness [N.mm/Rad]')
plt.xlabel('Lateral Stiffness in Y [N/mm]')
plt.semilogx()
plt.semilogy()
plt.grid()
plt.show()
#%%-----------------------------------------------------------------------------------
# %% Plot 2D Frame Shapes
S04.PLOT_2D_FRAME(deformed_scale=1)  # Adjust scale factor as needed
#%%-----------------------------------------------------------------------------------   
# EXPORT DATA TO EXCEL
import pandas as pd
DATA_TOTAL = {
    'DISP-X': displacementsX, 
    'DISP-Y': displacementsY, 
    'ROTATION': rotations,
    'BASE-SHEAR': forcesH,
    'BASE-AXIAL': forcesV,
    'BASE-MOMENT': forcesM,
    'COL-01-AXIAL': A_col01,
    'COL-01-SHEAR': S_col01,
    'COL-01-MOMENT': M_col01,
    'COL-02-AXIAL': A_col02,
    'COL-02-SHEAR': S_col02,
    'COL-02-MOMENT': M_col02,
    'COL-03-AXIAL': A_col03,
    'COL-03-SHEAR': S_col03,
    'COL-03-MOMENT': M_col03,
}
# Convert to DataFrame
results_df = pd.DataFrame(DATA_TOTAL)
# Export the DataFrame to an Excel file
results_df.to_excel('PROGRESSIVE_COLLAPSE_CONCRETE_2STORY_RESULTS.xlsx', index=False) 
#%%-----------------------------------------------------------------------------------


