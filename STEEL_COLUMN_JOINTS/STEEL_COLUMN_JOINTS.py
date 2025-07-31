#          #####################################################################################
#          #                                  IN THE NAME OF ALLAH                             #
#          #      PUSHOVER ANALYSIS OF STEEL COLUMNS WITH AXIAL AND ROTATIONAL SPRINGS         #
#          #                         FOR MODELING BEAM COLUMN JOINTS                           #
#          #-----------------------------------------------------------------------------------#
#          #              THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)           #
#          #                       EMAIL: salar.d.ghashghaei@gmail.com                         #
#          #####################################################################################
#%%-------------------------------------------------------------
# PAPER: Study on Structural Behaviour of Fully Bolted Beam Column Joints in Modular Steel Buildings
'https://link.springer.com/chapter/10.1007/978-3-031-12011-4_63'
# MASTER THESIS: MODELING OF INTERIOR BEAM-COLUMN JOINTS FOR NONLINEAR ANALYSIS OF REINFORCED CONCRETE FRAMES
'https://tspace.library.utoronto.ca/bitstream/1807/75836/3/Pan_Zhangcheng_201611_MAS_thesis.pdf'
# YOUTUBE: Beam column joint modelling in ABAQUS
'https://www.youtube.com/watch?v=j33VlhZ9oZo'
# YOUTUBE: ABAQUS Tutorial, Reinforced Concrete Beam-Column Joint Modeling, Analysis and behavior
'https://www.youtube.com/watch?v=AfeYRbJQfzw'
#%%-------------------------------------------------------------
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
    
image_path = 'STEEL_COLUMN_JOINTS.jpg'       
PLOT_IMAGE(image_path)
#%%-------------------------------------------------------------

"""
When OK equals -1, it generally indicates that the command or operation was not executed
because it was already in progress or had already been completed. This can happen if you
try to run a command that is already running or has been completed in a previous step.

When OK equals -2, it typically indicates that the command or operation was not executed
because it was not recognized or not implemented. This could mean that the command
is either misspelled, not available in the current version of OpenSees, or not applicable to the current context.

When OK equals -3, it typically means that the command or operation failed.
This could be due to various reasons, such as incorrect input parameters,
syntax errors, or issues with the model setup.
"""
def ANALYSIS(OK, INCREMENT, TOLERANCE, MAX_ITERAIONS):
    import openseespy.opensees as op
    test = {1:'NormDispIncr', 2: 'RelativeEnergyIncr', 4: 'RelativeNormUnbalance',5: 'RelativeNormDispIncr', 6: 'NormUnbalance'}
    algorithm = {1:'KrylovNewton', 2: 'SecantNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'BFGS', 7: 'Broyden', 8: 'NewtonLineSearch'}

    for i in test:
        for j in algorithm:
            if OK != 0:
                if j < 4:
                    op.algorithm(algorithm[j], '-initial')

                else:
                    op.algorithm(algorithm[j])

                op.test(test[i], TOLERANCE, MAX_ITERAIONS) 
                OK = op.analyze(INCREMENT)                            
                print(test[i], algorithm[j], OK)             
                if OK == 0:
                    break
            else:
                continue
#%%-------------------------------------------------------------
import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt

COL = 3000        # [mm] Column Length
PX = -1.0         # [N] HORIZENTAL FORCE - DOF[10]
PY = -5.0         # [N] VERTICAL FORCE - DOF[11]
MZ = 0.0          # [N.mm] MOMENT FORCE - DOF[12]
max_disp = 550.5  # [mm] MAXIMUM DISPLCEMENT - DOF[12]
disp_incr = 0.01  # [mm] EACH DISPLCEMENT INCREMENT - DOF[12]

# Define Analysis Properties
MAX_ITERATIONS = 20000     # Convergence iteration for test
MAX_TOLERANCE = 1.0e-10    # Convergence tolerance for test

# Create nodes
num_nodes = 4
node_coords = [(0, 0), (0, COL/3), (0, 2 * (COL/3)), (0, COL)]  # [mm] Column coordinates

# Define the model
ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 3)

# Structure Coordinate 
for i, coord in enumerate(node_coords):
    ops.node(i + 1, *coord)

# Fix the base node
ops.fix(1, 0, 1, 0) # Constraint DOF[1,2,3]

# Springs Coordinate 
ops.node(100, 0.0, 0.0)
ops.node(400, 0.0, COL)

# Fix the base Springs node
ops.fix(100, 1, 1, 1) # Springs Constraint DOF[1,2,3]
ops.fix(400, 1, 1, 1) # Springs Constraint DOF[10,11,12]

# Define material and section properties
Es = 210e3        # [N/mm^2] Young's modulus 
fy = 355          # [N/mm^2] Yield strength
ey = fy/Es        # [mm/mm] Steel Yield Strain
fu = 1.1818*fy    # [N/mm²] Steel Ultimate Strength
esu = ey*75.2     # [mm/mm] Steel Ultimate Strain
Esh = (fu - fy)/(esu - ey)
b = 0.01     # Strain hardening ratio
R0 = 17.0				# control the transition from elastic to plastic branches
cR1 = 0.925				# control the transition from elastic to plastic branches
cR2 = 0.15				# control the transition from elastic to plastic branches

# Define fiber section for nonlinear behavior using I-section
#ops.uniaxialMaterial('Steel01', 1, fy, Es, b)
ops.uniaxialMaterial('Steel02', 1, fy, Es, b, R0,  cR1, cR2)
ops.section('Fiber', 1)
web_thickness = 10  # [mm] Web thickness of the I-section
flange_width = 100  # [mm] Flange width of the I-section
flange_thickness = 20  # [mm] Flange thickness of the I-section
d = 300  # [mm] Total depth of the I-section
plate_thickness = 5  # [mm] Thickness of additional plates

# Define the fibers for the I-section
#           $matTag $numSubdivY $numSubdivZ $yI $zI $yJ $zJ
ops.patch('rect', 1, 20, 3, -flange_width / 2, d / 2, flange_width / 2, d / 2 - flange_thickness); # TOP FLANGE
ops.patch('rect', 1, 20, 3, -web_thickness / 2, d / 2 - flange_thickness, web_thickness / 2, -d / 2 + flange_thickness); # WEB
ops.patch('rect', 1, 20, 3, -flange_width / 2, -d / 2 + flange_thickness, flange_width / 2, -d / 2); # BOTTOM FLANGE
# Additional plates on top and bottom
ops.patch('rect', 1, 20, 3, -flange_width / 2, d / 2 + plate_thickness, flange_width / 2, d / 2)# TOP PLATE ON FLANGE
ops.patch('rect', 1, 20, 3, -flange_width / 2, -d / 2, flange_width / 2, -d / 2 - plate_thickness)# BOTTOM PLATE ON FLANGE
# Additional plates on left and right web
ops.patch('rect', 1, 20, 3, -flange_width / 2 - plate_thickness, d / 2 - flange_thickness, -flange_width / 2, -d / 2 + flange_thickness) # LEFT PLATE ON WEB
ops.patch('rect', 1, 20, 3, flange_width / 2, d / 2 - flange_thickness, flange_width / 2 + plate_thickness, -d / 2 + flange_thickness) # RIGHT PLATE ON WEB

# Define geometric transformation
ops.geomTransf('Linear', 1)
# ALL COLUMN HAS ONE SECTION AND 12 DEGREES OF FREEDOM
# Define elements as nonlinear beam-column elements
for i in range(num_nodes - 1):
#                                    $eleTag $iNode $jNode $numIntgrPts $secTag $transfTag    
    ops.element('nonlinearBeamColumn', i + 1, i + 1, i + 2, 5, 1, 1)

# Define nonlinear springs
# ALPHA: AXIAL SPRING 
# BETA: ROTATIONAL SPRING 
#ops.uniaxialMaterial('Steel01', 2, fy, Es, b)
#ops.uniaxialMaterial('Steel02', 2, fy, Es, b, R0,  cR1, cR2)
pinchX = 0.8   # Pinching factor in X direction
pinchY = 0.5   # Pinching factor in Y direction
damage1 = 0.0  # Damage due to ductility
damage2 = 0.0  # Damage due to energy
beta = 0.1     # Stiffness degradation parameter
ops.uniaxialMaterial('Hysteretic', 2, fy, ey, fu, esu, 0.2*fu, 1.1*esu, -fy, -ey, -fu, -0.5*esu, -0.2*fu, -1.1*esu, pinchX, pinchY, damage1, damage2, beta)
# INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Hysteretic_Material
#                   $eleTag $iNode $jNode -mat $matTag -dir $dir
ops.element('zeroLength', 4, 100, 1, '-mat', 2, '-dir', 1)  # DOF [1]
ops.element('zeroLength', 5, 100, 1, '-mat', 2, '-dir', 2)  # DOF [2]
ops.element('zeroLength', 6, 100, 1, '-mat', 2, '-dir', 3)  # DOF [3]
ops.element('zeroLength', 7, 400, 4, '-mat', 2, '-dir', 10) # DOF [10]
ops.element('zeroLength', 8, 400, 4, '-mat', 2, '-dir', 11) # DOF [11]
ops.element('zeroLength', 9, 400, 4, '-mat', 2, '-dir', 12) # DOF [12]

# Define load pattern for pushover analysis
ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)
ops.load(num_nodes, PX, PY, MZ)  # [N] Load at the top node (node 4)

# Define analysis parameters
ops.system('BandGeneral')
ops.numberer('RCM')
ops.constraints('Plain')
ops.test('NormUnbalance', MAX_TOLERANCE, MAX_ITERATIONS)
ops.algorithm('Newton')
ops.integrator('DisplacementControl', num_nodes, 1, disp_incr)
ops.analysis('Static')

# Perform pushover analysis
n_steps = int(max_disp / disp_incr)
displacements_x = []
displacements_y = []
drift = []
rotations = []
base_axials = []
base_shears = []
base_moments = []

for step in range(n_steps):
    OK = ops.analyze(1)
    ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS)
    disp_x = ops.nodeDisp(num_nodes, 1)
    drift_x = ((disp_x - ops.nodeDisp(1, 1)) / COL) * 100
    disp_y = ops.nodeDisp(num_nodes, 2)
    rotat = ops.nodeDisp(num_nodes, 3)
    axial_V_force = ops.eleResponse(1, 'force')[0] # SHEAR - Axial force in the zeroLength spring DOF (1)
    axial_A_force = ops.eleResponse(1, 'force')[1] # AXIAL - Axial force in the zeroLength spring DOF (2)
    axial_M_force = ops.eleResponse(1, 'force')[2] # MOMENT - Axial force in the zeroLength spring DOF (3)
    displacements_x.append(np.abs(disp_x))
    displacements_y.append(np.abs(disp_y))
    drift.append(np.abs(drift_x))
    rotations.append(np.abs(rotat))
    base_shears.append(np.abs(axial_V_force))
    base_axials.append(np.abs(axial_A_force))
    base_moments.append(np.abs(axial_M_force))

print('Analysis Done.')
#ops.wipe()    


#%%-------------------------------------------------------------
# Plot results
plt.figure(1, figsize=(10, 6))
plt.plot(displacements_y, base_axials, color='blue')
plt.xlabel('Displacement (mm) DOF[11]')
plt.ylabel('Base Axial Force (N) DOF[2]')
plt.title('Bottom Base-Axial Force & Top Displacement Absolute Data Diagram')
plt.grid(True)
plt.show()
#%%-------------------------------------------------------------
plt.figure(2 ,figsize=(10, 6))
plt.plot(displacements_x, base_shears, color='blue')
plt.xlabel('Displacement (mm) DOF[10]')
plt.ylabel('Base Shear (N) DOF[1]')
plt.title('Bottom Base-Shear & Top Displacement Absolute Data Diagram')
plt.grid(True)
plt.show()
#%%-------------------------------------------------------------
plt.figure(3 ,figsize=(10, 6))
plt.plot(rotations, base_moments, color='blue')
plt.xlabel('Rotation (Rad) DOF[12]')
plt.ylabel('Base Moment (N.mm) DOF[3]')
plt.title('Bottom Base-Moment & Bottom Rotation Absolute Data Diagram')
plt.grid(True)
plt.show()
#%%-------------------------------------------------------------
plt.figure(4 ,figsize=(10, 6))
plt.plot(drift, base_shears, color='blue')
plt.xlabel('Drift (%)')
plt.ylabel('Base Shear (N) DOF[1]')
plt.title('Bottom Base-Shear & Drift Absolute Data Diagram')
plt.grid(True)
plt.show()
#%%-------------------------------------------------------------



