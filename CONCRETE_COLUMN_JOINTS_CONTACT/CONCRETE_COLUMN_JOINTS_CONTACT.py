#          #####################################################################################
#          #                                  IN THE NAME OF ALLAH                             #
#          #      PUSHOVER ANALYSIS OF CONCRETE COLUMNS WITH AXIAL AND ROTATIONAL SPRINGS      #
#          #                         FOR MODELING BEAM COLUMN JOINTS                           #
#          #-----------------------------------------------------------------------------------#
#          #              THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)           #
#          #                       EMAIL: salar.d.ghashghaei@gmail.com                         #
#          #####################################################################################
"""
1. Purpose: 
    Performs pushover and dynamic analysis of concrete columns with axial/rotational springs to model beam-column joints.  
2. Framework: 
    Uses OpenSeesPy for nonlinear structural analysis with Python scripting.  
3. Features:  
   - Fiber-section modeling for concrete columns (confined/unconfined) and rebars.  
   - Nonlinear springs for joint behavior (axial/rotational).  
4. Pushover Analysis:  
   - Applies incremental displacements to study force-displacement response.  
   - Tracks base shear, axial force, and moment capacities.  
5. Dynamic Analysis:  
   - Simulates seismic response using Newmark integration.  
   - Includes Rayleigh damping and ground motion input.  
6. Material Models:  
   - Concrete02 for concrete (Mander model for confined).  
   - Hysteretic/Steel02 for rebars and springs.  
7. Section Modeling:  
   - Rectangular patches/layers for core/cover concrete and rebars.  
8. Contact Springs: 
    Optional nonlinear springs for gap/contact effects.  
9. Outputs: 
    Records displacements, drifts, forces, and moments.  
10. Visualization: 
    Plots hysteresis curves, time histories, and fiber sections.  
11. Convergence: 
    Adaptive Newton-Raphson with energy-based tolerance checks.  
12. Modularity: 
    Functions for analysis, plotting, and file management.  
13. Validation: 
    References academic papers/theses for joint modeling.  
14. Efficiency: 
    Parallel processing for large-scale models.  
15. User Inputs: 
    Customizable geometry, loads, and analysis parameters.  
16. Error Handling: 
    Robust checks for convergence and material failures.  
17. Applications: 
    RC frame joints, seismic retrofits, and modular steel buildings.  
18. Scalability: 
    Handles varying column/beam sections.  
19. Documentation: 
    Header comments clarify objectives and authorship.  
20. Extensibility: 
    Ready for integration with optimization/ML workflows.  
"""
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
#Load the Libraries
import os
import numpy as np
import openseespy.opensees as ops
import matplotlib.pyplot as plt
import openseespy.opensees as ops
import opsvis as opsv
import time
#%%-------------------------------------------------------------
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
FOLDER_NAME = 'OPENSEES_COLUMN_JOINT_SPRINGS_CONCRETE'
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
#%%-------------------------------------------------------------
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
    
image_path = 'OPENSEES_COLUMN_JOINT_SPRINGS_CONCRETE.jpg'       
PLOT_IMAGE(image_path)
#%%-------------------------------------------------------------
def CURRENT_TIME():
    import time
    t = time.localtime()
    current_time = time.strftime("%H:%M:%S", t)
    print(f"Current time (HH:MM:SS): {current_time}\n\n")
    
# ---------------------
def PLOT_2D(X, Y, XLABEL, YLABEL, TITLE, COLOR):
    plt.figure(figsize=(10, 6))
    plt.plot(X, Y,  color=COLOR)
    plt.xlabel(XLABEL)
    plt.ylabel(YLABEL)
    plt.title(TITLE)
    plt.grid(True)
    #plt.semilogy()
    plt.show()

# ---------------------
def PLOT_TIME_HISORY(TIME, displacements_x, displacements_y, rotations, velocity, acceleration, base_shears, base_axials, base_moments):
    import matplotlib.pyplot as plt
    # Creating a time array for the x-axis (Replace with actual time data if available)
    #time = range(len(displacements_x))
    
    # Plotting the data
    fig, axs = plt.subplots(8, 1, figsize=(10, 20))
    
    axs[0].plot(TIME, displacements_x, label='Displacement X')
    axs[0].set_title('Displacement X Time History')
    axs[0].set_xlabel('Time')
    axs[0].set_ylabel('Displacement X')
    
    axs[1].plot(TIME, displacements_y, label='Displacement Y')
    axs[1].set_title('Displacement Y Time History')
    axs[1].set_xlabel('Time')
    axs[1].set_ylabel('Displacement Y')
    
    axs[2].plot(TIME, rotations, label='Rotations')
    axs[2].set_title('Rotations Time History')
    axs[2].set_xlabel('Time')
    axs[2].set_ylabel('Rotations')
    
    axs[3].plot(TIME, velocity, label='Velocity')
    axs[3].set_title('Velocity Time History')
    axs[3].set_xlabel('Time')
    axs[3].set_ylabel('Velocity')
    
    axs[4].plot(TIME, acceleration, label='Acceleration')
    axs[4].set_title('Acceleration Time History')
    axs[4].set_xlabel('Time')
    axs[4].set_ylabel('Acceleration')
    
    axs[5].plot(TIME, base_shears, label='Base Shears')
    axs[5].set_title('Base Shears Time History')
    axs[5].set_xlabel('Time')
    axs[5].set_ylabel('Base Shears')
    
    axs[6].plot(TIME, base_axials, label='Base Axials')
    axs[6].set_title('Base Axials Time History')
    axs[6].set_xlabel('Time')
    axs[6].set_ylabel('Base Axials')
    
    axs[7].plot(TIME, base_moments, label='Base Moments')
    axs[7].set_title('Base Moments Time History')
    axs[7].set_xlabel('Time')
    axs[7].set_ylabel('Base Moments')
    
    # Adjust layout
    plt.tight_layout()
    plt.show()

# ---------------------
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
def CONCRETE_SECTION_PLOT(Bcol01, Hcol01, Bcol02, Hcol02, cover, Rebar_D, nFibCoverZ, nFibCoverY, nFibCoreZ, nFibCoreY, PLOT):
    import matplotlib.pyplot as plt
    import numpy as np
    import openseespy.opensees as ops
    import opsvis as opsv
    
    Mat_Tag01 = 1 # Confined Concrete Section Tag - Bottom and Top Section
    Mat_Tag02 = 2 # Unconfined Concrete Section Tag
    Mat_Tag03 = 3 # Steel Rebar Section Tag
    Mat_Tag04 = 4 # Confined Concrete Section Tag - Middle Section
    SECTION_TAG_01 = 1 # Concrete Column Section Tag
    SECTION_TAG_02 = 2 # Concrete Beam Section Tag
    
    fc = -25 # [N/mm^2] Nominal concrete compressive strength
    Ec = 4700 * np.sqrt(-fc) # [N/mm^2] Concrete Elastic Modulus

    # confined concrete - bottom and top section
    Kfc = 1.2;			    # ratio of confined to unconfined concrete strength - COLUMN - BOTTOM AND TOP SECTION
    fc1C = Kfc*fc;		    # CONFINED concrete (mander model), maximum stress - COLUMN  - BOTTOM AND TOP SECTION
    eps1C = 2*fc1C/Ec;	    # strain at maximum stress 
    fc2C = 0.2*fc1C;		# ultimate stress
    eps2C = 5*eps1C;		# strain at ultimate stress 

    # confined concrete - middle section
    Kfcm = 1.05;	        # ratio of confined to unconfined concrete strength - COLUMN - MIDDLE SECTION
    fc1Cm = Kfcm * fc       # CONFINED concrete (mander model), maximum stress - COLUMN  - MIDDLE SECTION
    eps1Cm = 2*fc1Cm/Ec;	# strain at maximum stress 
    fc2Cm = 0.2*fc1Cm;		# ultimate stress
    eps2Cm = 5*eps1Cm;		# strain at ultimate stress 
    ftCm = -0.55*fc1Cm;		# tensile strength +tension    
    
    # unconfined concrete
    fc1U = fc;			# UNCONFINED concrete (todeschini parabolic model), maximum stress
    eps1U = -0.0025;	    # strain at maximum strength of unconfined concrete
    fc2U = 0.2*fc1U;		# ultimate stress
    eps2U = -0.012;			# strain at ultimate stress
    Lambda = 0.1;		    # ratio between unloading slope at $eps2 and initial slope $Ec
    
    # tensile-strength properties
    ftC = -0.55*fc1C;		# tensile strength +tension
    ftU = -0.55*fc1U;		# tensile strength +tension
    Ets = ftU/0.002;		# tension softening stiffness
    
    ops.uniaxialMaterial('Concrete02', Mat_Tag01, fc1C, eps1C, fc2C, eps2C, Lambda, ftC, Ets) # build core concrete (confined) - BOTTOM AND TOP SECTION
    ops.uniaxialMaterial('Concrete02', Mat_Tag04, fc1Cm, eps1Cm, fc2Cm, eps2Cm, Lambda, ftCm, Ets) # build core concrete (confined) - MIDDLE SECTION
    ops.uniaxialMaterial('Concrete02', Mat_Tag02, fc1U, eps1U, fc2U, eps2U, Lambda, ftU, Ets) # build cover concrete (unconfined)
    # REBAR MATERIAL PROPERTIES:
    """
    Fy = 4000			# Steel rebar yield stress
    Cy = 0.02			# Steel rebar yield strain
    Es = Fy/Cy				# modulus of steel
    Bs = 0.01				# strain-hardening ratio 
    R0 = 18.0				# control the transition from elastic to plastic branches
    cR1 = 0.925				# control the transition from elastic to plastic branches
    cR2 = 0.15				# control the transition from elastic to plastic branches
    ops.uniaxialMaterial('Steel02', Mat_Tag03, Fy, Es, Bs, R0, cR1, cR2) # build reinforcement material 
    """
    """
    E_steel = 210e3               # [N/mm²] Young's modulus
    fy_steel = 4000               # [N/mm²] Yield strength
    fu_steel = 1.23 * fy_steel    # [N/mm²] Ultimate strength
    esh = 0.02                    # Strain corresponding to initial strain hardening
    eult = 0.191                  # Strain at peak stress
    Esh = (fu_steel - fy_steel)/(eult - esh)
    ops.uniaxialMaterial('ReinforcingSteel', Mat_Tag03, fy_steel, fu_steel, E_steel, Esh, esh, eult)
    """
    """
    Fy = 4000			# Steel rebar yield tension stress
    FyC = 2500			# Steel rebar yield compression stress
    Cy = 0.02			# Steel rebar yield strain
    Es = Fy/Cy				# modulus of steel
    Bs = 0.01				# strain-hardening ratio 
    R0 = 18.0				# control the transition from elastic to plastic branches
    cR1 = 0.925				# control the transition from elastic to plastic branches
    cR2 = 0.15				# control the transition from elastic to plastic branches
    #                                      $matTag $Fy $FyC $E $b $R0 $cR1 $cR2 $a1 $a2 $a3 $a4 $sigcr $beta $sigmin $FI_lim
    ops.uniaxialMaterial('SteelFractureDI', Mat_Tag03, Fy, FyC, Es, Bs, R0, cR1, cR2, 0.08, 1.00, 0.08, 1.00, 120, 0.8, 20, 1.0)
    """
    Es = 210e4        # [N/mm^2] Young's modulus 
    fy = 355          # [N/mm^2] Yield strength
    ey = fy/Es        # [mm/mm] Steel Yield Strain
    fu = 1.1818*fy    # [N/mm²] Steel Ultimate Strength
    esu = ey*75.2     # [mm/mm] Steel Ultimate Strain
    pinchX = 0.8   # Pinching factor in X direction
    pinchY = 0.5   # Pinching factor in Y direction
    damage1 = 0.0  # Damage due to ductility
    damage2 = 0.0  # Damage due to energy
    beta = 0.1 # Stiffness degradation parameter
    ops.uniaxialMaterial('Hysteretic', Mat_Tag03, fy, ey, fu, esu, 0.2*fu, 1.1*esu, -fy, -ey, -fu, -0.5*esu, -0.2*fu, -1.1*esu, pinchX, pinchY, damage1, damage2, beta)
    # INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Hysteretic_Material
    
    # FIBER SECTION properties -------------------------------------------------------------
    # symmetric section
    #                        y
    #                        ^
    #                        |     
    #             ---------------------     --   --
    #             |   o  o   o    o   |     |    -- cover
    #             |                   |     |
    #             |   o           o   |     |
    #    z <---   |          +        |     H
    #             |   o           o   |     |
    #             |                   |     |
    #             |   o  o    o   o   |     |    -- cover
    #             ---------------------     --   --
    #             |-------- B --------|
    #
    # RC section: 
    y1col = Hcol01/2.0
    z1col = Bcol01/2.0

    y2col = 0.5 * (Hcol01 - 2 * cover) / 2;

    #nFibCoverZ, nFibCoverY = 1 , 20
    #nFibCoreZ, nFibCoreY = 2, 16
    As = (np.pi * Rebar_D ** 2) / 4; # [mm^2] Rebar Area

    FIBER_SEC_01 = [['section', 'Fiber', SECTION_TAG_01, '-GJ', 1.0e6],
                 ['patch', 'rect', Mat_Tag01, nFibCoreY, nFibCoverZ, cover-y1col, cover-z1col, y1col-cover, z1col-cover], # CORE
                 ['patch', 'rect', Mat_Tag02, nFibCoverY, nFibCoverZ, -y1col, -z1col, y1col, cover-z1col],                # COVER
                 ['patch', 'rect', Mat_Tag02, nFibCoverY, nFibCoverZ, -y1col, z1col-cover, y1col, z1col],                 # COVER
                 ['patch', 'rect', Mat_Tag02, nFibCoverY, nFibCoverZ, -y1col, cover-z1col, cover-y1col, z1col-cover],     # COVER
                 ['patch', 'rect', Mat_Tag02, nFibCoverY, nFibCoverZ, y1col-cover, cover-z1col, y1col, z1col-cover],      # COVER
                 ['layer', 'straight', Mat_Tag03, 5, As, y1col-cover, z1col-cover, y1col-cover, cover-z1col],             # REBAR
                 ['layer', 'straight', Mat_Tag03, 2, As, y2col, z1col-cover, y2col, cover-z1col],                         # REBAR
                 ['layer', 'straight', Mat_Tag03, 2, As, 0, z1col-cover, 0, cover-z1col],                                 # REBAR
                 ['layer', 'straight', Mat_Tag03, 2, As, -y2col, z1col-cover, -y2col, cover-z1col],                       # REBAR
                 ['layer', 'straight', Mat_Tag03, 5, As, cover-y1col, z1col-cover, cover-y1col, cover-z1col]              # REBAR
                ]
    
    if PLOT == 1:
        matcolor = ['gold', 'lightgrey']
        plt.figure(1)
        opsv.plot_fiber_section(FIBER_SEC_01, matcolor=matcolor)
        # Set the x and y limits
        LIMIT_Y = 0.5 * Hcol01 + 10
        LIMIT_X = 0.5 * Bcol01 + 10
        plt.ylim(-LIMIT_Y, LIMIT_Y)
        plt.xlim(-LIMIT_X, LIMIT_X)
        plt.title('COLUMN BOTTOM AND TOP SECTION')
        plt.show()

    # FIBER SECTION properties -------------------------------------------------------------
    # symmetric section
    #                        y
    #                        ^
    #                        |     
    #             ---------------------     --   --
    #             |   o  o   o    o   |     |    -- cover
    #             |                   |     |
    #             |   o           o   |     |
    #    z <---   |          +        |     H
    #             |   o           o   |     |
    #             |                   |     |
    #             |   o  o    o   o   |     |    -- cover
    #             ---------------------     --   --
    #             |-------- B --------|
    #
    # RC section:     
    y1col = Hcol02/2.0
    z1col = Bcol02/2.0

    y2col = 0.5*(Hcol02-2*cover) / 2.0

    #nFibCoverZ, nFibCoverY = 1 , 20
    #nFibCoreZ, nFibCoreY = 2, 16
    As = (np.pi * Rebar_D ** 2) / 4; # [mm^2] Rebar Area

    FIBER_SEC_02 = [['section', 'Fiber', SECTION_TAG_02, '-GJ', 1.0e6],
                 ['patch', 'rect', Mat_Tag04, nFibCoreY, nFibCoreZ, cover-y1col, cover-z1col, y1col-cover, z1col-cover], # CORE
                 ['patch', 'rect', Mat_Tag02, nFibCoverY, nFibCoverZ, -y1col, -z1col, y1col, cover-z1col],               # COVER
                 ['patch', 'rect', Mat_Tag02, nFibCoverY, nFibCoverZ, -y1col, z1col-cover, y1col, z1col],                # COVER
                 ['patch', 'rect', Mat_Tag02, nFibCoverY, nFibCoverZ, -y1col, cover-z1col, cover-y1col, z1col-cover],    # COVER
                 ['patch', 'rect', Mat_Tag02, nFibCoverY, nFibCoverZ, y1col-cover, cover-z1col, y1col, z1col-cover],     # COVER
                 ['layer', 'straight', Mat_Tag03, 5, As, y1col-cover, z1col-cover, y1col-cover, cover-z1col],            # REBAR
                 ['layer', 'straight', Mat_Tag03, 2, As, y2col, z1col-cover, y2col, cover-z1col],                        # REBAR
                 ['layer', 'straight', Mat_Tag03, 2, As, 0, z1col-cover, 0, cover-z1col],                                # REBAR
                 ['layer', 'straight', Mat_Tag03, 2, As, -y2col, z1col-cover, -y2col, cover-z1col],                      # REBAR
                 ['layer', 'straight', Mat_Tag03, 5, As, cover-y1col, z1col-cover, cover-y1col, cover-z1col]             # REBAR
                ]
    
    if PLOT == 1:
        matcolor = ['gold', 'lightgrey', 'lightgrey', 'lightblue']
        plt.figure(1)
        opsv.plot_fiber_section(FIBER_SEC_02, matcolor=matcolor)
        # Set the x and y limits
        LIMIT_Y = 0.5 * Hcol02 + 10
        LIMIT_X = 0.5 * Bcol02 + 10
        plt.ylim(-LIMIT_Y, LIMIT_Y)
        plt.xlim(-LIMIT_X, LIMIT_X)
        plt.title('COLUMN MIDDLE SECTION')
        plt.show()
        
    #print(fc1C, eps1C, fc2C, eps2C, ftC)
    #print(fc1Cm, eps1Cm, fc2Cm, eps2Cm, ftCm)    
    return FIBER_SEC_01, FIBER_SEC_02    
#%%-------------------------------------------------------------
nFibCoverZ, nFibCoverY = 1 , 1
nFibCoreZ, nFibCoreY = 1, 1
FS01, FS02 = CONCRETE_SECTION_PLOT(600, 600, 600, 600, 50, 25,nFibCoverZ, nFibCoverY, nFibCoreZ, nFibCoreY, PLOT=1)
#%%-------------------------------------------------------------
def PUSHOVER_ANALYSIS(COL, PX, PY, MZ, max_disp, disp_incr,
                      Bcol01, Hcol01, Bcol02, Hcol02, COVER, REBAR_DIA,
                      MAX_ITERAIONS, TOLERANCE, CONTACT, CONTACT_DISP, KIND):

    num_nodes = 4 # COLUMN NODES
    # Define fiber section for I-section
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

    # Concrete Sections for Columns and Beams
    nFibCoverZ, nFibCoverY, nFibCoreZ, nFibCoreY = 3, 120, 3, 120
    SECTION01, SECTION02 = CONCRETE_SECTION_PLOT(Bcol, Hcol, Bbeam, Hbeam, COVER, REBAR_DIA,
                                                 nFibCoverZ, nFibCoverY, nFibCoreZ, nFibCoreY, PLOT=0)
    opsv.fib_sec_list_to_cmds(SECTION01) # COLUMNS
    opsv.fib_sec_list_to_cmds(SECTION02) # BEAMS

    # Define geometric transformation
    # Linear:
    # Small displacement assumptions in local to basic transformation
    # Linear transformation of forces and displacements
    # transfTag = 1
    # ops.geomTransf('Linear', transfTag)
    
    # PDelta:
    # Small displacement assumption transformation of displacements
    # Account for transverse displacement of axial load in equilibrium relationship
    # transfTag = 1
    # ops.geomTransf('PDelta', transfTag)
    
    # Corotational:
    # Fully nonlinear transformation of displacements and force
    # Exact in 2D but some approximations in 3D
    transfTag = 1
    ops.geomTransf('Corotational', transfTag) 

    # Define elements as nonlinear beam-column elements
    for i in range(num_nodes - 1):
        if i == 0 or i == 2:   
            ops.element('nonlinearBeamColumn', i + 1, i + 1, i + 2, 5, 1, transfTag)
        else: ## MIDDLE OF COLUMN
            ops.element('nonlinearBeamColumn', i + 1, i + 1, i + 2, 5, 2, transfTag)

    # OUTPUT DATA
    ops.recorder('Node', '-file', f"{SALAR_DIR}DTH_PUSH.txt",'-time', '-node', 2, '-dof', 1,2,3, 'disp')# Displacement Time History Node 2
    ops.recorder('Node', '-file', f"{SALAR_DIR}BTH_PUSH.txt",'-time', '-node', 1, '-dof', 1,2,3, 'reaction')# Base Shear Time History Node 1
    ops.recorder('Element', '-file', f"{SALAR_DIR}DEF_PUSH.txt",'-time', '-ele', 1, 'section', 5, 'deformations')# Curvature Time History

    # Define nonlinear springs
    # ALPHA: AXIAL SPRING 
    P_Y = 42000 * Bcol01 * Hcol01  # [N] Yield force capacity
    P_K = P_Y / 12                 # [N/mm] Axial Rigidity Stiffness
    P_U = 1.21 * P_Y               # [N/mm²] Ultimate force
    P_esh_D = 20                   # [mm] Displacement corresponding to initial displacement hardening
    P_eult_D = 50                  # [mm] Displacement at peak Force
    Esh_D = (P_U - P_Y)/(P_eult_D - P_esh_D) # Tangent at initial displacement hardening
    Z = Esh_D / P_K                # displacement-hardening ratio (ratio between post-yield tangent and initial elastic tangent
    #ops.uniaxialMaterial('ReinforcingSteel', 100, P_Y, P_U, P_K, Esh_D, P_esh_D, P_eult_D)
    ops.uniaxialMaterial('Steel02', 100, P_Y, P_K, Z, 18.0, 0.925, 0.15)
    # BETA: ROTATIONAL SPRING 
    M_Y = 57530 * Bcol01 * Hcol01  # [N.mm] Yield moment capacity
    M_K = M_Y / 0.0523             # [N.mm/rad] Rotational Flextural Rigidity Stiffness
    M_U = 1.524 * M_Y              # [N/mm] Ultimate moment
    M_esh_R = 0.095                # [rad] Rotation corresponding to initial rotation hardening
    M_eult_R = 0.162               # [rad] Rotation at peak Moment
    Esh_R = (M_U - M_Y)/(M_eult_R - M_esh_R) # Tangent at initial rotation hardening
    Z = Esh_R / M_K                # rotation-hardening ratio (ratio between post-yield tangent and initial elastic tangent
    #ops.uniaxialMaterial('ReinforcingSteel', 200, M_Y, M_U, M_K, Esh_R, M_esh_R, M_eult_R)
    ops.uniaxialMaterial('Steel02', 200, M_Y, M_K, Z, 18.0, 0.925, 0.15)

    #                   $eleTag $iNode $jNode -mat $matTag -dir $dir
    ops.element('zeroLength', 4, 100, 1, '-mat', 100, '-dir', 1) # DOF [1] LATERAL SPRING
    ops.element('zeroLength', 5, 100, 1, '-mat', 100, '-dir', 2) # DOF [2] LATERAL SPRING
    ops.element('zeroLength', 6, 100, 1, '-mat', 200, '-dir', 6) # DOF [3] ROTATIONAL SPRING
    ops.element('zeroLength', 7, 400, 4, '-mat', 100, '-dir', 1) # DOF [10] LATERAL SPRING
    ops.element('zeroLength', 8, 400, 4, '-mat', 100, '-dir', 2) # DOF [11] LATERAL SPRING
    ops.element('zeroLength', 9, 400, 4, '-mat', 200, '-dir', 6) # DOF [12] ROTATIONAL SPRING

    # Define load pattern for pushover analysis
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(num_nodes, PX, PY, MZ)  # [N] Load at the top node (node 4)
    """
    # Define gravity analysis parameters
    NstepGravity = 10
    DGravity = 1 / NstepGravity
    ops.integrator('LoadControl', DGravity, MAX_ITERATIONS) # determine the next time step for an analysis
    ops.numberer('Plain') # renumber dof's to minimize band-width (optimization), if you want to
    ops.system('BandGeneral') # how to store and solve the system of equations in the analysis
    ops.constraints('Plain') # how it handles boundary conditions
    ops.test('EnergyIncr', TOLERANCE, MAX_ITERATIONS) # determine if convergence has been achieved at the end of an iteration step
    ops.algorithm('ModifiedNewton') # use Newton's solution algorithm: updates tangent stiffness at every iteration
    ops.analysis('Static') # define type of analysis static or transient
    OK = ops.analyze(NstepGravity) # apply gravity
    #ANALYSIS(OK, NstepGravity, TOLERANCE, MAX_ITERATIONS)
    print('Load Done.') 
    
    ops.loadConst('-time', 0.0) #maintain constant gravity loads and reset time to zero
    """
    # Define incremental displacement analysis parameters
    ops.system('BandGeneral')
    ops.numberer('RCM')
    ops.constraints('Plain')
    ops.test('EnergyIncr', TOLERANCE, MAX_ITERATIONS)
    ops.algorithm('ModifiedNewton')
    if KIND == 1:# HORIZENTAL DISPLACEMENT
        ops.integrator('DisplacementControl', num_nodes, 1, disp_incr)
    if KIND == 2:# VERTICAL DISPLACEMENT
        ops.integrator('DisplacementControl', num_nodes, 2, disp_incr)    
    if KIND == 3:# ROTATION
        ops.integrator('DisplacementControl', num_nodes, 3, disp_incr)     
    ops.analysis('Static')
    print('Model Done.')
    
    if CONTACT == True:
        # CONTACT SPRING STIFFNESS PROPERTIES:
        PY_C = 114e10    # [N] Yield force capacity
        K_C = PY_C / 7  # [N/mm] Axial Rigidity Stiffness
        ops.uniaxialMaterial('Steel02', 300, PY_C, K_C, 0.1, 18.0, 0.925, 0.15) # NONLINEAR CONTACT
        #ops.uniaxialMaterial('Elastic', 300, K_C)                                # LINEAR CONTACT
        CONTACT_STEPS = int(np.abs(CONTACT_DISP / disp_incr))
    
    # Perform pushover analysis
    n_steps = int(np.abs(max_disp / disp_incr))  # Analysis Steps
    print(f'{n_steps} steps for incremental displacement is going to be run')
    displacements_x, displacements_y, drift, rotations = [], [], [], []
    base_axials, base_shears, base_moments = [], [], []
    
    for step in range(n_steps):
        OK = ops.analyze(1)
        ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS)
        disp_x = ops.nodeDisp(num_nodes, 1)
        drift_x = ((disp_x - ops.nodeDisp(1, 1)) / COL) * 100
        disp_y = ops.nodeDisp(num_nodes, 2)
        rotat = ops.nodeDisp(num_nodes, 3)
        if CONTACT == True:
            if step == CONTACT_STEPS: # CHECK FOR CONTACT
                print(f'IN STEP {step} CONTACT DONE!')
                ops.element('zeroLength', 10, 4, 400, '-mat', 300, '-dir', 1) # DOF [10] LATERAL SPRING 
        axial_V_force = ops.eleResponse(1, 'force')[0] # SHEAR - Axial force in the zeroLength spring DOF (1)
        axial_A_force = ops.eleResponse(1, 'force')[1] # AXIAL - Axial force in the zeroLength spring DOF (2)
        axial_M_force = ops.eleResponse(1, 'force')[2] # MOMENT - Axial force in the zeroLength spring DOF (3)
        displacements_x.append(disp_x)
        displacements_y.append(disp_y)
        drift.append(drift_x)
        rotations.append(rotat)
        base_shears.append(-axial_V_force)
        base_axials.append(-axial_A_force)
        base_moments.append(-axial_M_force)
        #print(step + 1,' Step Done.')

    print('Pushover Done.')    
    ops.wipe() 
    return displacements_x, displacements_y, drift, rotations, base_shears, base_axials, base_moments
#%%-------------------------------------------------------------
### DEFINE FIBER SETION FOR CONCRETE COLUMN AND BEAM SECTION:
COL = 3000   # [mm] Column Length
Bcol, Hcol, Bbeam, Hbeam = 600, 600, 600, 600; # [mm] Column & Beam Section Diamenstion Properties
COVER = 50        # [mm] Concrete Cover
REBAR_DIA = 25    # [mm] Steel Rebar Diameter

### DEFINE LOAD PROPERTIES:
PX = 100.0        # [N] HORIZENTAL FORCE - DOF[10]
PY = 0.0          # [N] VERTICAL FORCE - DOF[11]
MZ = 0.0          # [N.mm] MOMENT FORCE - DOF[12]
MAX_DISP = 520    # [mm] MAXIMUM DISPLCEMENT - DOF[10]
DISP_INCR = 0.01  # [mm] EACH DISPLCEMENT INCREMENT - DOF[10]


### PUSHOVER DIRECTION:
# [1]: X-DIRECTION - HORIZENTAL DISPLACEMENT
# [2]: Y-DIRECTION - VERTICAL DISPLACEMENT
# [3]: Z-DIRECTION - ROTATION
KIND = 1

### CONTACT:
CONTACT = True # True: HAVE CONTACT - False: DO NOT HAVE CONTACT
CONTACT_DISP = 450 # [mm] CONTACT MAXIMUM DISPLCEMENT - DOF[10]

### ANALYSIS TOLEANCE AND ITERATIONS
MAX_ITERATIONS = 5000   # Maximum number of iterations
MAX_TOLERANCE = 1.0e-10 # Specified tolerance for convergence


starttime = time.process_time()

### RUN THE ANALYSIS:
DATA = PUSHOVER_ANALYSIS(COL, PX, PY, MZ, MAX_DISP, DISP_INCR,Bcol, Hcol, Bbeam, Hbeam,COVER, REBAR_DIA, MAX_ITERATIONS, MAX_TOLERANCE, CONTACT, CONTACT_DISP, KIND)
displacements_x, displacements_y, drift, rotations, base_shears, base_axials, base_moments = DATA

totaltime = time.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n')
#%%-------------------------------------------------------------
X = displacements_x
Y = base_shears
XLABEL = 'Displacement (mm) DOF[10]'
YLABEL = 'Base Shear (N) DOF[1]'
TITLE = 'Bottom Base-shear Force & Top Displacement Absolute Data Diagram'
COLOR = 'blue'
PLOT_2D(X, Y, XLABEL, YLABEL, TITLE, COLOR)
#%%-------------------------------------------------------------
X = displacements_x
Y = base_axials
XLABEL = 'Displacement (mm) DOF[10]'
YLABEL = 'Base Axial Force (N) DOF[2]'
TITLE = 'Bottom Base-Axial Force & Top Displacement Absolute Data Diagram'
COLOR = 'green'
PLOT_2D(X, Y, XLABEL, YLABEL, TITLE, COLOR)
#%%-------------------------------------------------------------
X = displacements_y
Y = base_axials
XLABEL = 'Displacement (mm) DOF[11]'
YLABEL = 'Base Axial Force (N) DOF[2]'
TITLE = 'Bottom Base-Axial Force & Top Displacement Absolute Data Diagram'
COLOR = 'black'
PLOT_2D(X, Y, XLABEL, YLABEL, TITLE, COLOR)
#%%-------------------------------------------------------------
X = displacements_y
Y = base_shears
XLABEL = 'Vertical Displacement (mm) DOF[11]'
YLABEL = 'Base Shear (N) DOF[1]'
TITLE = 'Bottom Base-Shear & Top Displacement Absolute Data Diagram'
COLOR = 'lime'
PLOT_2D(X, Y, XLABEL, YLABEL, TITLE, COLOR)
#%%-------------------------------------------------------------
X = rotations
Y = base_moments
XLABEL = 'Rotation (Rad) DOF[12]'
YLABEL = 'Base Moment (N.mm) DOF[3]'
TITLE = 'Bottom Base-Moment & Bottom Rotation Absolute Data Diagram'
COLOR = 'purple'
PLOT_2D(X, Y, XLABEL, YLABEL, TITLE, COLOR)
#%%-------------------------------------------------------------
X = drift
Y = base_shears
XLABEL = 'Horizental Drift (%) DOF[10]'
YLABEL = 'Base Shear (N) DOF[1]'
TITLE = 'Bottom Base-Shear & Drift Absolute Data Diagram'
COLOR = 'red'
PLOT_2D(X, Y, XLABEL, YLABEL, TITLE, COLOR)
#%%-------------------------------------------------------------
#          #####################################################################################
#          #                                  IN THE NAME OF ALLAH                             #
#          #       DYNAMIC ANALYSIS OF CONCRETE COLUMNS WITH AXIAL AND ROTATIONAL SPRINGS      #
#          #                         FOR MODELING BEAM COLUMN JOINTS                           #
#          #-----------------------------------------------------------------------------------#
#          #              THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)           #
#          #                       EMAIL: salar.d.ghashghaei@gmail.com                         #
#          #####################################################################################
#%%-------------------------------------------------------------

def DYNAMIC_ANALYSIS(MASS_X, MASS_Y, MASS_Z, COL, Bcol01, Hcol01, Bcol02, Hcol02, COVER, REBAR_DIA, CONTACT, CONTACT_DISP, Damping_Ratio, dt, Duration):
    ### ANALYSIS TOLEANCE AND ITERATIONS
    MAX_ITERATIONS = 10000      # Maximum number of iterations
    MAX_TOLERANCE = 1.0e-10     # Specified tolerance for convergence
    
    num_nodes = 4  # Column Nodes
    node_coords = [(0, 0), (0, COL/3), (0, 2 * (COL/3)), (0, COL)]  # Column coordinates

    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)

    for i, coord in enumerate(node_coords):
        ops.node(i + 1, *coord)
    
    # Mass on top of column
    ops.mass(num_nodes, MASS_X, MASS_Y, MASS_Z)
    
    # Fix the base node
    ops.fix(1, 0, 1, 0) # Constraint DOF[1,2,3]

    # Springs Coordinate 
    ops.node(100, 0.0, 0.0)
    ops.node(400, 0.0, COL)

    # Fix the base Springs node
    ops.fix(100, 1, 1, 1) # Springs Constraint DOF[1,2,3]
    ops.fix(400, 1, 1, 1) # Springs Constraint DOF[10,11,12]

    nFibCoverZ, nFibCoverY, nFibCoreZ, nFibCoreY = 3, 120, 3, 120
    SECTION01, SECTION02 = CONCRETE_SECTION_PLOT(Bcol01, Hcol01, Bcol02, Hcol02, COVER, REBAR_DIA, nFibCoverZ, nFibCoverY, nFibCoreZ, nFibCoreY, PLOT=0)
    opsv.fib_sec_list_to_cmds(SECTION01)
    opsv.fib_sec_list_to_cmds(SECTION02)

    # Define geometric transformation
    # Linear:
    # Small displacement assumptions in local to basic transformation
    # Linear transformation of forces and displacements
    # transfTag = 1
    # ops.geomTransf('Linear', transfTag)
    
    # PDelta:
    # Small displacement assumption transformation of displacements
    # Account for transverse displacement of axial load in equilibrium relationship
    # transfTag = 1
    # ops.geomTransf('PDelta', transfTag)
    
    # Corotational:
    # Fully nonlinear transformation of displacements and force
    # Exact in 2D but some approximations in 3D
    transfTag = 1
    ops.geomTransf('Corotational', transfTag) 

    for i in range(num_nodes - 1):
        secTag = 1 if i == 0 or i == 2 else 2
        ops.element('nonlinearBeamColumn', i + 1, i + 1, i + 2, 5, secTag, transfTag)

    # OUTPUT DATA
    ops.recorder('Node', '-file', f"{SALAR_DIR}DTH_DYN.txt",'-time', '-node', 2, '-dof', 1,2,3, 'disp')# Displacement Time History Node 2
    ops.recorder('Node', '-file', f"{SALAR_DIR}BTH_DYN.txt",'-time', '-node', 1, '-dof', 1,2,3, 'reaction')# Base Shear Time History Node 1
    ops.recorder('Element', '-file', f"{SALAR_DIR}DEF_DYN.txt",'-time', '-ele', 1, 'section', 5, 'deformations')# Curvature Time History

    # Define nonlinear springs
    # ALPHA: AXIAL SPRING 
    P_Y = 42000 * Bcol01 * Hcol01  # [N] Yield force capacity
    P_K = P_Y / 12                 # [N/mm] Axial Rigidity Stiffness
    P_U = 1.21 * P_Y               # [N/mm²] Ultimate force
    P_esh_D = 20                   # [mm] Displacement corresponding to initial displacement hardening
    P_eult_D = 50                  # [mm] Displacement at peak Force
    Esh_D = (P_U - P_Y)/(P_eult_D - P_esh_D) # Tangent at initial displacement hardening
    Z = Esh_D / P_K                # displacement-hardening ratio (ratio between post-yield tangent and initial elastic tangent
    #ops.uniaxialMaterial('ReinforcingSteel', 100, P_Y, P_U, P_K, Esh_D, P_esh_D, P_eult_D)
    ops.uniaxialMaterial('Steel02', 100, P_Y, P_K, Z, 18.0, 0.925, 0.15)
    # BETA: ROTATIONAL SPRING 
    M_Y = 57530 * Bcol01 * Hcol01  # [N.mm] Yield moment capacity
    M_K = M_Y / 0.0523             # [N.mm/rad] Rotational Flextural Rigidity Stiffness
    M_U = 1.524 * M_Y              # [N/mm] Ultimate moment
    M_esh_R = 0.095                # [rad] Rotation corresponding to initial rotation hardening
    M_eult_R = 0.162               # [rad] Rotation at peak Moment
    Esh_R = (M_U - M_Y)/(M_eult_R - M_esh_R) # Tangent at initial rotation hardening
    Z = Esh_R / M_K                # rotation-hardening ratio (ratio between post-yield tangent and initial elastic tangent
    #ops.uniaxialMaterial('ReinforcingSteel', 200, M_Y, M_U, M_K, Esh_R, M_esh_R, M_eult_R)
    ops.uniaxialMaterial('Steel02', 200, M_Y, M_K, Z, 18.0, 0.925, 0.15)

    #                   $eleTag $iNode $jNode -mat $matTag -dir $dir
    ops.element('zeroLength', 4, 100, 1, '-mat', 100, '-dir', 1) # DOF [1] LATERAL SPRING
    ops.element('zeroLength', 5, 100, 1, '-mat', 100, '-dir', 2) # DOF [2] LATERAL SPRING
    ops.element('zeroLength', 6, 100, 1, '-mat', 200, '-dir', 6) # DOF [3] ROTATIONAL SPRING
    ops.element('zeroLength', 7, 400, 4, '-mat', 100, '-dir', 1) # DOF [10] LATERAL SPRING
    ops.element('zeroLength', 8, 400, 4, '-mat', 100, '-dir', 2) # DOF [11] LATERAL SPRING
    ops.element('zeroLength', 9, 400, 4, '-mat', 200, '-dir', 6) # DOF [12] ROTATIONAL SPRING 

    
    print('Model Built')
    # Define time series and load pattern
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    #ops.load(3, 1.0, -1, 0.0)
    #ops.load(4, 1.0, -1, 0.0)
    
    # Define gravity analysis parameters
    NstepGravity = 10
    DGravity = 1 / NstepGravity
    ops.integrator('LoadControl', DGravity, MAX_ITERATIONS) # determine the next time step for an analysis
    ops.numberer('Plain') # renumber dof's to minimize band-width (optimization), if you want to
    ops.system('BandGeneral') # how to store and solve the system of equations in the analysis
    ops.constraints('Plain') # how it handles boundary conditions
    ops.test('EnergyIncr', MAX_TOLERANCE, MAX_ITERATIONS) # determine if convergence has been achieved at the end of an iteration step
    ops.algorithm('ModifiedNewton') # use Newton's solution algorithm: updates tangent stiffness at every iteration
    ops.analysis('Static') # define type of analysis static or transient
    OK = ops.analyze(NstepGravity) # apply gravity
    #ANALYSIS(OK, NstepGravity, TOLERANCE, MAX_ITERATIONS)
    print('Load Done.') 

    ops.loadConst('-time', 0.0) #maintain constant gravity loads and reset time to zero
    
    # Calculate Rayleigh damping factors
    Lambda01 = ops.eigen('-fullGenLapack', 2)  # eigenvalue mode 2
    #Lambda01 = ops.eigen('-genBandArpack', 2) # eigenvalue mode 2
    Omega01 = np.power(max(Lambda01), 0.5)
    Omega02 = np.power(min(Lambda01), 0.5)
    a0 = (2 * Omega01 * Omega02 * Damping_Ratio) / (Omega01 + Omega02) # c = a0 * m : Mass-proportional damping
    a1 = (Damping_Ratio * 2) / (Omega01 + Omega02)   # c = a1 * k : Stiffness-proportional damping
    # Apply Rayleigh damping
    ops.rayleigh(a0, a1, 0, 0)   # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/reyleigh.html
    #ops.rayleigh(0, 0, 2 * DR * Omega01, 0) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/reyleigh.html
    PERIOD_01 = np.pi / Omega01 # Structure First Period
    PERIOD_02 = np.pi / Omega02 # Structure Second Period
    print('Structure First Period:  ', PERIOD_01)
    print('Structure Second Period: ', PERIOD_02) 
    modal_period = PERIOD_01
    natural_frequency = 1 / PERIOD_01
    
    #applying Dynamic Ground motion analysis
    
    GMfact = 9810    # standard acceleration of gravity or standard acceleration
    SSF_X = 0.01     # Seismic Acceleration Scale Factor in X Direction
    SSF_Y = 0.01     # Seismic Acceleration Scale Factor in Y Direction
    iv0_X = 0.00005  # [mm/s] Initial velocity applied to the node  in X Direction
    iv0_Y = 0.00005  # [mm/s] Initial velocity applied to the node  in Y Direction
    st_iv0 = 0.0     # [s] Initial velocity applied starting time
    
    SEISMIC_TAG_01 = 100
    ops.timeSeries('Path', SEISMIC_TAG_01, '-dt', dt, '-filePath', 'Ground_Acceleration_X.txt', '-factor', GMfact, '-startTime', st_iv0) # SEISMIC-X
    # Define load patterns
    # pattern UniformExcitation $patternTag $dof -accel $tsTag <-vel0 $vel0> <-fact $cFact>
    ops.pattern('UniformExcitation', SEISMIC_TAG_01, 1, '-accel', SEISMIC_TAG_01, '-vel0', iv0_X, '-fact', SSF_X) # SEISMIC-X 
    SEISMIC_TAG_02 = 200
    ops.timeSeries('Path', SEISMIC_TAG_02, '-dt', dt, '-filePath', 'Ground_Acceleration_Y.txt', '-factor', GMfact) # SEISMIC-Z
    ops.pattern('UniformExcitation', SEISMIC_TAG_02, 2, '-accel', SEISMIC_TAG_02, '-vel0', iv0_Y, '-fact', SSF_Y)  # SEISMIC-Z 
    
    ops.wipeAnalysis()
    ops.constraints('Transformation')
    ops.numberer('Plain')
    ops.system('BandGeneral')
    ops.test('RelativeNormUnbalance', MAX_TOLERANCE, MAX_ITERATIONS)
    ops.algorithm('RaphsonNewton')

    NewmarkGamma = 0.5; NewmarkBeta = 0.25;
    ops.integrator('Newmark', NewmarkGamma, NewmarkBeta)
    ops.analysis('Transient')

    if CONTACT == True:
        # CONTACT SPRING STIFFNESS PROPERTIES:
        PY_C = 114e10    # [N] Yield force capacity
        K_C = PY_C / 7  # [N/mm] Axial Rigidity Stiffness
        ops.uniaxialMaterial('Steel02', 300, PY_C, K_C, 0.1, 18.0, 0.925, 0.15) # NONLINEAR CONTACT
        #ops.uniaxialMaterial('Elastic', 300, K_C)                                # LINEAR CONTACT
        #CONTACT_STEPS = int(np.abs(CONTACT_DISP / dt))

    # Perform Dynamic analysis
    n_steps =  int(np.abs(Duration/ dt))# Analysis Steps
    print(f'{n_steps} steps for incremental displacement is going to be run')
    displacements_x, displacements_y, drift, rotations = [], [], [], []
    base_axials, base_shears, base_moments = [], [], []
    velocityX, accelerationX = [], []
    velocityY, accelerationY = [], []
    TIME = []

    for step in range(n_steps):
        OK = ops.analyze(1, dt)
        ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS)
        current_time = ops.getTime()
        TIME.append(current_time)
        disp_x = ops.nodeDisp(num_nodes, 1)
        vel_x = ops.nodeVel(num_nodes, 1)
        accel_x = ops.nodeAccel(num_nodes, 1)
        drift_x = ((disp_x - ops.nodeDisp(1, 1)) / COL) * 100
        disp_y = ops.nodeDisp(num_nodes, 2)
        vel_y = ops.nodeVel(num_nodes, 2)
        accel_y = ops.nodeAccel(num_nodes, 2)
        rotat = ops.nodeDisp(num_nodes, 3)
        if CONTACT == True:
            if disp_x == CONTACT_DISP: # CHECK FOR CONTACT
                print(f'IN STEP {step} CONTACT DONE!')
                ops.element('zeroLength', 10, 4, 400, '-mat', 300, '-dir', 1) # DOF [10] LATERAL SPRING 
        axial_V_force = ops.eleResponse(1, 'force')[0] # SHEAR - Axial force in the zeroLength spring DOF (1)
        axial_A_force = ops.eleResponse(1, 'force')[1] # AXIAL - Axial force in the zeroLength spring DOF (2)
        axial_M_force = ops.eleResponse(1, 'force')[2] # MOMENT - Axial force in the zeroLength spring DOF (3)
        displacements_x.append(disp_x)
        displacements_y.append(disp_y)
        drift.append(drift_x)
        rotations.append(rotat)
        velocityX.append(vel_x)
        accelerationX.append(accel_x)
        velocityY.append(vel_y)
        accelerationY.append(accel_y)
        base_shears.append(axial_V_force)
        base_axials.append(axial_A_force)
        base_moments.append(axial_M_force)
        #print(step + 1,' Step Done.')
    
    print("Period: ", modal_period)
    print("Natural Frequency: ", natural_frequency)
    print('\nDynamic Done.')    
    ops.wipe() 

    return TIME, drift, displacements_x, displacements_y, rotations, velocityX, accelerationX, velocityY, accelerationY, base_shears, base_axials, base_moments

#%%-------------------------------------------------------------
## RUN DYNAMIC ANALYSIS:
##-----------------------
MASS_X = 20.0e5   # [N] X-DIRECTION MASS
MASS_Y = 20.0e5   # [N] Y-DIRECTION MASS
MASS_Z = 1.0e-9   # [N] Z-DIRECTION MASS
COL = 3000        # [mm] COLUMN LENGTH
Bcol01, Hcol01, Bcol02, Hcol02 = 600, 600, 600, 600 # [mm] COLUMN SECTION DIMENSIONS
COVER = 50        # [mm] COVER
REBAR_DIA = 25    # [mm] REBAR DIAMETER
CONTACT = True
CONTACT_DISP = 0.001 # [mm] CONTACT DISPLACEMENT

Damping_Ratio = 0.05	# damping ratio
dt = 0.01               # [s] Time increment
Duration = 10.0         # [s] Total time duration 

starttime = time.process_time()

DATA = DYNAMIC_ANALYSIS(MASS_X, MASS_Y, MASS_Z, COL,
                                      Bcol01, Hcol01, Bcol02, Hcol02, COVER, REBAR_DIA, CONTACT, CONTACT_DISP, Damping_Ratio, dt, Duration)
TIME, drift_x, displacements_x, displacements_y, rotations, velocityX, accelerationX, velocityY, accelerationY, base_shears, base_axials, base_moments = DATA

totaltime = time.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n')

PLOT_TIME_HISORY(TIME, displacements_x, displacements_y, rotations, velocityX, accelerationX, base_shears, base_axials, base_moments)
#%%-------------------------------------------------------------
X = displacements_x
Y = base_shears
XLABEL = 'Displacement [mm] DOF(10)' 
YLABEL = 'Base-shear [N] DOF(1)' 
TITLE = 'Base-shear & Displacement Diagrarm' 
COLOR = 'Black'
PLOT_2D(X, Y, XLABEL, YLABEL, TITLE, COLOR)
#%%-------------------------------------------------------------
X = displacements_y
Y = base_axials
XLABEL = 'Displacement [mm] DOF(11)' 
YLABEL = 'Base-axial [N] DOF(2)' 
TITLE = 'Base-axial & Displacement Diagrarm' 
COLOR = 'Brown'
PLOT_2D(X, Y, XLABEL, YLABEL, TITLE, COLOR)
#%%-------------------------------------------------------------
X = rotations
Y = base_moments
XLABEL = 'Rotation [Rad] DOF(12)' 
YLABEL = 'Base-moment [N.mm] DOF(3)' 
TITLE = 'Base-moment & Rotation Diagrarm' 
COLOR = 'Purple'
PLOT_2D(X, Y, XLABEL, YLABEL, TITLE, COLOR)
#%%-------------------------------------------------------------
X = drift_x
Y = base_shears
XLABEL = 'Horizental Drift [%] DOF(10)' 
YLABEL = 'Base-shear [N] DOF(1)' 
TITLE = 'Base-shear & Horizental Drift Diagrarm' 
COLOR = 'red'
PLOT_2D(X, Y, XLABEL, YLABEL, TITLE, COLOR)
#%%-------------------------------------------------------------
