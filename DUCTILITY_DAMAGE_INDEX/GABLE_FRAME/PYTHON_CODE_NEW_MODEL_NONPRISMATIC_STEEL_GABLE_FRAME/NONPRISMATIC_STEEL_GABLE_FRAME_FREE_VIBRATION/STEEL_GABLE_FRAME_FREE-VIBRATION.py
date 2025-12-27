#  #########################################################################
#  #                      >> IN THE NAME OF ALLAH <<                       #
#  #    COMPARATIVE FREE-VIBRATION ANALYSIS OF A MDOF STRUCTURE:           #
#  #              ELASTIC VS INELASTIC RESPONSE USING OPENSEES             #
#  #-----------------------------------------------------------------------#
#  #     NONLINEAR DYNAMIC ANALYSIS OF A NONPRISMATIC STEEL GABLE FRAME    #
#  #              I SECTION COLUMN WITH FINITE PRISMATIC COLUMN            #
#  #-----------------------------------------------------------------------#
#  #    MODELING OF NONPRISMATIC ELEMENT WITH MULTI PRISMATIC ELEMENTS     #
#  #-----------------------------------------------------------------------#
#  #      THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)       #
#  #                   EMAIL: salar.d.ghashghaei@gmail.com                 #
#  #########################################################################

#%%-------------------------------------------------------------
#import the os module
import os
import time as TI
import numpy as np
import openseespy.opensees as op
import ANALYSIS_FUNCTION as S01
import PLOT_2D as S02
import DAMPING_RATIO_FUN as S03
import EIGENVALUE_ANALYSIS_FUN as S04
#%%-------------------------------------------------------------
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

# Load the image
image_path = 'STEEL_GABLE_FRAME_FREE-VIBRATION.PNG'
image = mpimg.imread(image_path)

# Display the image
plt.figure(figsize=(30, 16))
plt.imshow(image)
plt.axis('off')  # Hide axes
plt.show()
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
FOLDER_NAME = 'OPENSEES_STEEL_GABLE_FRAME'
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
def MAXABS_FUN(DATA_FILE, COLUMN):
    import numpy as np
    # Read and process displacement data
    NameFiles = DATA_FILE
    filename = f"{NameFiles}.txt"
    D = np.loadtxt(filename)
    #print(D)
    MAXABS = np.max(np.abs([D[:, COLUMN]]))
    #print("MAX. ABS. :", MAXABS)
    return MAXABS
# -----------------------------------------------
def PLOT_2D(X, Y, XLABEL, YLABEL, TITLE, COLOR):
    plt.figure(figsize=(10, 6))
    plt.plot(X, Y,  color=COLOR)
    plt.xlabel(XLABEL)
    plt.ylabel(YLABEL)
    plt.title(TITLE)
    plt.grid(True)
    #plt.semilogy()
    plt.show()
# -----------------------------------------------
def OUTPUT_SECOND_COLUMN(FOLDER_NAME, X, COLUMN):
    import numpy as np
    # Time History
    filename = f"C:\OPENSEESPY_SALAR\{FOLDER_NAME}\\{X}.txt"
    data_collected = np.loadtxt(filename)
    X = data_collected[:, COLUMN]
    return X 
# -----------------------------------------------
def PLOT_TIME_HISORY(TIME, displacements_x, displacements_y, rotations, velocity_X, acceleration_X, velocity_Y, acceleration_Y, base_shears, base_axials, base_moments):
    import matplotlib.pyplot as plt
    # Creating a time array for the x-axis (Replace with actual time data if available)
    #time = range(len(displacements_x))
    
    # Plotting the data
    fig, axs = plt.subplots(10, 1, figsize=(12, 24))
    
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
    
    axs[3].plot(TIME, velocity_X, label='Velocity')
    axs[3].set_title('Velocity Time History')
    axs[3].set_xlabel('Time')
    axs[3].set_ylabel('Velocity in X Dir.')
    
    axs[4].plot(TIME, acceleration_X, label='Acceleration')
    axs[4].set_title('Acceleration Time History')
    axs[4].set_xlabel('Time')
    axs[4].set_ylabel('Acceleration in X Dir.')
    
    axs[5].plot(TIME, velocity_Y, label='Velocity')
    axs[5].set_title('Velocity Time History')
    axs[5].set_xlabel('Time')
    axs[5].set_ylabel('Velocity in Y Dir.')
    
    axs[6].plot(TIME, acceleration_Y, label='Acceleration')
    axs[6].set_title('Acceleration Time History')
    axs[6].set_xlabel('Time')
    axs[6].set_ylabel('Acceleration in Y Dir.')
    
    axs[7].plot(TIME, base_shears, label='Base Shears')
    axs[7].set_title('Base Shears Time History')
    axs[7].set_xlabel('Time')
    axs[7].set_ylabel('Base Shears')
    
    axs[8].plot(TIME, base_axials, label='Base Axials')
    axs[8].set_title('Base Axials Time History')
    axs[8].set_xlabel('Time')
    axs[8].set_ylabel('Base Axials')
    
    axs[9].plot(TIME, base_moments, label='Base Moments')
    axs[9].set_title('Base Moments Time History')
    axs[9].set_xlabel('Time')
    axs[9].set_ylabel('Base Moments')
    
    # Adjust layout
    plt.tight_layout()
    plt.show()

# ---------------------
#%%-------------------------------------------------------------
def FREE_VIBRAION_ANALYSIS(nfCoreY, nfCoreZ, MAX_ITERATIONS, MAX_TOLERANCE,
                           u0, v0, a0, IA, IU, IV,
                           hw_COL, tw_COL, tf_COL, bf_COL, hw_BEAM, tw_BEAM,
                           tf_BEAM, bf_BEAM, H1, H2, L1, N, ND, Weight, UL,
                           Damping_Ratio, dt, Duration, ELE_KIND, DIST_LOAD):
    ### ANALYSIS TOLEANCE AND ITERATIONS
    MAX_ITERATIONS = 5000   # Maximum number of iterations
    MAX_TOLERANCE = 1.0e-6  # Specified tolerance for convergence
    
    NN = int(N / 4) # 4 ELEMENTS
    
    # Define the model
    op.wipe()
    op.model('basic', '-ndm', 2, '-ndf', 3) 
    PCol = Weight    # nodal dead-load weight per column
    GMfact = 9810    # [mm/s^2] standard acceleration of gravity or standard acceleration
    Mass =  PCol/GMfact
    h1 = H1 / (NN - 1)
    h2 = H2 / (NN - 1)
    #bl = (L1 **2 + H2 **2)**0.5  # Beam Length
    l1 = L1 / (NN) 
          
    
    # Define geometric transformation
    # Linear:
    # Small displacement assumptions in local to basic transformation
    # Linear transformation of forces and displacements
    # transfTag = 1
    # op.geomTransf('Linear', transfTag)
    
    # PDelta:
    # Small displacement assumption transformation of displacements
    # Account for transverse displacement of axial load in equilibrium relationship
    # transfTag = 1
    # op.geomTransf('PDelta', transfTag)
    
    # Corotational:
    # Fully nonlinear transformation of displacements and force
    # Exact in 2D but some approximations in 3D
    transfTag = 1
    op.geomTransf('Corotational', transfTag)
    
    # MATERIAL parameters -------------------------------------------------------------------
    Mat_Tag = 1; 				# material ID tag -- steel 
    """
    Fy = 240			# STEEL yield stress
    Cy = 0.0012			# STEEL yield stress
    Es = Fy/Cy				# modulus of steel
    Bs = 0.01				# strain-hardening ratio 
    R0 = 18.0				# control the transition from elastic to plastic branches
    cR1 = 0.925				# control the transition from elastic to plastic branches
    cR2 = 0.15				# control the transition from elastic to plastic branches
    op.uniaxialMaterial('Steel02', Mat_Tag, Fy, Es, Bs, R0,cR1,cR2) # build reinforcement material
    """
    Es = 2.10e5       # [N/mm^2] Young's modulus 
    fy = 240          # [N/mm^2] Yield strength
    ey = fy/Es        # [mm/mm] Steel Yield Strain
    fu = 1.1818*fy    # [N/mm²] Steel Ultimate Strength
    esu = 0.2         # [mm/mm] Steel Ultimate Strain
    pinchX = 0.8      # Pinching factor in X direction
    pinchY = 0.5      # Pinching factor in Y direction
    damage1 = 0.0     # Damage due to ductility
    damage2 = 0.0     # Damage due to energy
    beta = 0.1        # Stiffness degradation parameter
    if ELE_KIND == 'INELASTIC':
        op.uniaxialMaterial('HystereticSM', Mat_Tag,
                            fy, ey, fu, esu, 0.2*fu, 1.1*esu, 0.1*fu, 1.2*esu,
                            -fy, -ey, -fu, -0.5*esu, -0.2*fu, -1.1*esu, -0.1*fu, -1.2*esu,
                            pinchX, pinchY, damage1, damage2, beta)
        # INFO LINK: https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/HystereticSM.html
        #Esh = (fu-fy)/(esu-ey); b = Esh/Es; R0  = 10; CR1  = 0.925; CR2 = 0.15; 
        #op.uniaxialMaterial('Steel02', Mat_Tag, fy, Es, b, R0 ,CR1 ,CR2)
        # INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Steel02_Material_--_Giuffr%C3%A9-Menegotto-Pinto_Model_with_Isotropic_Strain_Hardening
    if ELE_KIND == 'ELASTIC':
        op.uniaxialMaterial('Elastic', Mat_Tag, Es)
        # INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Elastic_Uniaxial_Material
    # ---------------------------------------------------------------------------------------
    hwbot_COL = 0.5 * hw_COL # Bottom column web height
    
    # Nodal coordinates:
    op.node(1, 0.0, 0.0) # node X, Y
    op.node(N, L1 * 2, 0.0) # node X, Y
    
    # Constraints -- Boundary Conditions
    op.fix(1, 1, 1, 1) # node DX DY RZ
    op.fix(N, 1, 1, 1) # node DX DY RZ
    
    # ----------------------------------------------------------
    # COLUMN 01:
    hwbot_COL = 0.5 * hw_COL # Bottom column web height - SECTION 01
    for i in range(2, NN + 1, 1):
        pp = i-1
        x = 0; y = h1 * pp;
        op.node(i, x, y)
        ColSecTag = pp			# assign a tag number to the column section
        ColTransfTag = pp
        eleTag = pp

        # FIBER SECTION properties -------------------------------------------------------------
        # symmetric section
        #                        y
        #                        ^
        #                        |     
        #              _____________________    --   --
        #             |_________   ________|    |    -- tf
        #                      |  |             |
        #                      |  |             |
        #    z <---       hw   |tw|             H
        #                      |  |             |
        #              ________|  |_________    |
        #             |____________________|    |    -- tf
        #                                      --    --
        #             |-------- bf --------|
        #
        # STEEL I SECTION: 
        HW_COL =  hwbot_COL + ((hw_COL - hwbot_COL) / H1) * h1 * (i-1) # Varying web height of nonprismatic section - SECTION 02
        #print(i, x, y, i-1, i, HW_COL)
        coverY = (HW_COL + tf_COL) / 2.0
        coverZ = tw_COL / 2.0
        coreY = coverY - tf_COL
        coreZ02 = bf_COL / 2.0

        #nfCoreY = 20;			# number of fibers for steel in y-direction
        #nfCoreZ = 5;			# number of fibers for steel in z-direction


        op.section('Fiber', ColSecTag)
        # Define the core patch
        op.patch('quad', Mat_Tag, nfCoreZ, nfCoreY, -coverY, coreZ02, -coverY, -coreZ02, -coreY,-coreZ02, coreY, coreZ02) # BOTTOM FLANGE
        op.patch('quad', Mat_Tag, nfCoreZ, nfCoreY, -coreY,coverZ, -coreY,-coverZ, coreY,-coverZ, coreY, coverZ) # MIDDLE WEB
        op.patch('quad', Mat_Tag, nfCoreZ, nfCoreY, coreY, coreZ02, -coreY, coreZ02, coverY,-coreZ02, coreY, coreZ02) # TOP FLANGE

        numIntgrPts = 5
        op.element('nonlinearBeamColumn', eleTag, i-1, i, numIntgrPts, ColSecTag, transfTag)
    # ----------------------------------------------------------
    # BEAM 01:
    hwbot_BEAM = 0.5 * hw_BEAM # j section beam web height  - SECTION 04
    for i in range(NN + 1, 2 * NN + 1, 1):
        zz = i - NN
        x = l1 * zz; y = H1 + h2 * zz;
        op.node(i, x, y)
        ColSecTag = i-1			# assign a tag number to the column section
        ColTransfTag = i-1
        eleTag = i-1

        # FIBER SECTION properties -------------------------------------------------------------
        # symmetric section
        #                        y
        #                        ^
        #                        |     
        #              _____________________    --   --
        #             |_________   ________|    |    -- tf
        #                      |  |             |
        #                      |  |             |
        #    z <---       hw   |tw|             H
        #                      |  |             |
        #              ________|  |_________    |
        #             |____________________|    |    -- tf
        #                                      --    --
        #             |-------- bf --------|
        #
        # STEEL I SECTION: 
        HW_BEAM =  hw_BEAM + ((hwbot_BEAM - hw_BEAM) / L1) * l1 * zz # Varying web height of nonprismatic section - SECTION 03
        #print(i, x, y, i-1, i, HW_BEAM)
        coverY = (HW_BEAM + tf_BEAM) / 2.0
        coverZ = tw_BEAM / 2.0
        coreY = coverY - tf_BEAM
        coreZ02 = bf_BEAM / 2.0

        #nfCoreY = 20;			# number of fibers for steel in y-direction
        #nfCoreZ = 5;			# number of fibers for steel in z-direction


        op.section('Fiber', ColSecTag)
        # Define the core patch
        op.patch('quad', Mat_Tag, nfCoreZ, nfCoreY, -coverY, coreZ02, -coverY, -coreZ02, -coreY,-coreZ02, coreY, coreZ02) # BOTTOM FLANGE
        op.patch('quad', Mat_Tag, nfCoreZ, nfCoreY, -coreY,coverZ, -coreY,-coverZ, coreY,-coverZ, coreY, coverZ) # MIDDLE WEB
        op.patch('quad', Mat_Tag, nfCoreZ, nfCoreY, coreY, coreZ02, -coreY, coreZ02, coverY,-coreZ02, coreY, coreZ02) # TOP FLANGE

        numIntgrPts = 5
        op.element('nonlinearBeamColumn', eleTag, i-1, i, numIntgrPts, ColSecTag, transfTag) 
    # ----------------------------------------------------------
    # COLUMN 02:
    hwbot_COL = 0.5 * hw_COL # Bottom column web height
    for i in range(N-1, 3 * NN, -1):
        pp = i - 1
        zz = N - i
        #print(pp)
        x = L1 * 2; y = h1 * zz;
        #print(x, y)
        op.node(i, x, y)
        ColSecTag = i			# assign a tag number to the column section
        ColTransfTag = i
        eleTag = i

        # FIBER SECTION properties -------------------------------------------------------------
        # symmetric section
        #                        y
        #                        ^
        #                        |     
        #              _____________________    --   --
        #             |_________   ________|    |    -- tf
        #                      |  |             |
        #                      |  |             |
        #    z <---       hw   |tw|             H
        #                      |  |             |
        #              ________|  |_________    |
        #             |____________________|    |    -- tf
        #                                      --    --
        #             |-------- bf --------|
        #
        # STEEL I SECTION: 
        HW_COL =  hwbot_COL + ((hw_COL - hwbot_COL) / H1) * h1 * (N - i) # Varying web height of nonprismatic section
        #print(i, x, y, i+1, i, HW_COL)
        coverY = (HW_COL + tf_COL) / 2.0
        coverZ = tw_COL / 2.0
        coreY = coverY - tf_COL
        coreZ02 = bf_COL / 2.0

        #nfCoreY = 20			# number of fibers for steel in y-direction
        #nfCoreZ = 5;			# number of fibers for steel in z-direction


        op.section('Fiber', ColSecTag)
        # Define the core patch
        op.patch('quad', Mat_Tag, nfCoreZ, nfCoreY, -coverY, coreZ02, -coverY, -coreZ02, -coreY,-coreZ02, coreY, coreZ02) # BOTTOM FLANGE
        op.patch('quad', Mat_Tag, nfCoreZ, nfCoreY, -coreY,coverZ, -coreY,-coverZ, coreY,-coverZ, coreY, coverZ) # MIDDLE WEB
        op.patch('quad', Mat_Tag, nfCoreZ, nfCoreY, coreY, coreZ02, -coreY, coreZ02, coverY,-coreZ02, coreY, coreZ02) # TOP FLANGE

        numIntgrPts = 5
        op.element('nonlinearBeamColumn', eleTag, i+1, i, numIntgrPts, ColSecTag, transfTag)  
    # ----------------------------------------------------------
    # BEAM 02:
    hwbot_BEAM = 0.5 * hw_BEAM # j section beam web height
    for i in range(3 * NN, 2 * NN, -1):
        zz = 3 * NN + 1 - i
        x = 2 * L1 - l1 * zz; y = H1 + h2 * zz;
        op.node(i, x, y) 
        ColSecTag = i			# assign a tag number to the column section
        #ColTransfTag = i
        eleTag = i


        # FIBER SECTION properties -------------------------------------------------------------
        # symmetric section
        #                        y
        #                        ^
        #                        |     
        #              _____________________    --   --
        #             |_________   ________|    |    -- tf
        #                      |  |             |
        #                      |  |             |
        #    z <---       hw   |tw|             H
        #                      |  |             |
        #              ________|  |_________    |
        #             |____________________|    |    -- tf
        #                                      --    --
        #             |-------- bf --------|
        #
        # STEEL I SECTION: 
        HW_BEAM =  hw_BEAM + ((hwbot_BEAM - hw_BEAM) / L1) * l1 * zz # Varying web height of nonprismatic section
        #print(i, x, y, i+1, i, HW_BEAM)
        coverY = (HW_BEAM + tf_BEAM) / 2.0
        coverZ = tw_BEAM / 2.0
        coreY = coverY - tf_BEAM
        coreZ02 = bf_BEAM / 2.0

        #nfCoreY = 20;			# number of fibers for steel in y-direction
        #nfCoreZ = 5;			# number of fibers for steel in z-direction

        op.section('Fiber', ColSecTag)
        # Define the core patch
        op.patch('quad', Mat_Tag, nfCoreZ, nfCoreY, -coverY, coreZ02, -coverY, -coreZ02, -coreY,-coreZ02, coreY, coreZ02) # BOTTOM FLANGE
        op.patch('quad', Mat_Tag, nfCoreZ, nfCoreY, -coreY,coverZ, -coreY,-coverZ, coreY,-coverZ, coreY, coverZ) # MIDDLE WEB
        op.patch('quad', Mat_Tag, nfCoreZ, nfCoreY, coreY, coreZ02, -coreY, coreZ02, coverY,-coreZ02, coreY, coreZ02) # TOP FLANGE

        numIntgrPts = 5
        op.element('nonlinearBeamColumn', eleTag, i+1, i, numIntgrPts, ColSecTag, transfTag)  
        
        
    op.element('nonlinearBeamColumn', 200, 102, 100, numIntgrPts, ColSecTag, transfTag) ## THIS ELEMENT IS FOR CONNECTION 2 SEPARATE ELEMENTS TO EACH OTHER  
    #import InelasticFiberSection
    op.recorder('Node', '-file', f"{SALAR_DIR}DTH_DYN.txt",'-time', '-node', ND, '-dof', 1,2,3, 'disp')# Displacement Time History Node 150
    op.recorder('Node', '-file', f"{SALAR_DIR}BTH_DYN_01.txt",'-time', '-node', 1, '-dof', 1,2,3, 'reaction')# Base Shear Time History Node 1
    op.recorder('Node', '-file', f"{SALAR_DIR}BTH_DYN_200.txt",'-time', '-node', N, '-dof', 1,2,3, 'reaction')# Base Shear Time History Node 200
    
    #print('Model Built')
    
    # node#, Mx My Mz, Mass=Weight/g, neglect rotational inertia at nodes
    #defining mass
    op.mass(1*NN, Mass, Mass, 0.0)
    #op.mass(2*NN, Mass, Mass, 0.0)
    op.mass(3*NN, Mass, Mass, 0.0)
    
    #defining gravity loads
    op.timeSeries('Linear', 1)
    op.pattern('Plain', 1, 1)
    op.load(ND, 1.0, 0.0, 0.0)
    """
    if DIST_LOAD == True:
        for i in range(NN + 1, 2 * NN + 1, 1): # CREATE UNIFORM LOADS FOR BEAM 01
            op.eleLoad('-ele', i-1,'-type', '-beamUniform', UL, 0.0) # uniformly-distributed load
        for i in range(3 * NN, 2 * NN, -1):    # CREATE UNIFORM LOADS FOR BEAM 02
            op.eleLoad('-ele', i,'-type', '-beamUniform', UL, 0.0) # uniformly-distributed load
    
    #print('Unifrom Load Defined.')
    
    NstepGravity = 10
    DGravity = 1 / NstepGravity
    op.integrator('LoadControl', DGravity) # determine the next time step for an analysis
    op.numberer('Plain') # renumber dof's to minimize band-width (optimization), if you want to
    op.system('BandGeneral') # how to store and solve the system of equations in the analysis
    op.constraints('Plain') # how it handles boundary conditions
    op.test('NormDispIncr', MAX_TOLERANCE, 6) # determine if convergence has been achieved at the end of an iteration step
    op.algorithm('Newton') # use Newton's solution algorithm: updates tangent stiffness at every iteration
    op.analysis('Static') # define type of analysis static or transient
    op.analyze(NstepGravity) # apply gravity
    
    #print('Gravity Anlysis Done.')

    op.loadConst('-time', 0.0) #maintain constant gravity loads and reset time to zero
    
    #print('Period Anlysis Done.')
    
    #defining lateral loads
    op.timeSeries('Linear', 2)
    op.pattern('Plain', 2, 1)
    op.load(ND, 1.0, 0.0, 0.0)
    """
    if IU == True:
        # Define initial displacment
        op.setNodeDisp(ND, 1, u0, '-commit') # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/setNodeDisp.html
    if IV == True:
        # Define initial velocity
        op.setNodeVel(ND, 1, v0, '-commit')  # INFO LINK: https://openseespydoc.readthedocs.io/en/stable/src/setNodeVel.html
    if IA == True:
        # Define initial  acceleration
        op.setNodeAccel(ND, 1, a0, '-commit') # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/setNodeAccel.html
    #print('Free-vibration Defined.')
    
    #op.wipeAnalysis()
    op.constraints('Plain')
    op.numberer('Plain')
    op.system('BandGeneral')
    op.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS)
    NewmarkGamma = 0.5; NewmarkBeta = 0.25;
    op.integrator('Newmark', NewmarkGamma, NewmarkBeta) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/newmark.html
    #alpha = 1;gamma=1.5-alpha; beta=((2-alpha)**2)/4;
    #op.integrator('HHT', alpha, gamma, beta) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/hht.html
    op.algorithm('RaphsonNewton') # INFO LINK: https://openseespydoc.readthedocs.io/en/stable/src/algorithm.html
    op.analysis('Transient')

    n_steps =  int(np.abs(Duration/ dt))# Analysis Steps
    
    displacements_x, displacements_y, drift, rotations = [], [], [], []
    base_axials, base_shears, base_moments = [], [], []
    velocityX, accelerationX = [], []
    velocityY, accelerationY = [], []
    TIME = []
    stiffness = []
    PERIOD_MIN, PERIOD_MAX = [], []
    
    #OK = op.analyze(n_steps, dt)
    #S03.ANALYSIS(OK, Nsteps, MAX_TOLERANCE, MAX_ITERATIONS)
    
    for step in range(n_steps):
        OK = op.analyze(1, dt)
        S01.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS)
        current_time = op.getTime()
        TIME.append(current_time)
        disp_x = op.nodeDisp(ND, 1)
        vel_x = op.nodeVel(ND, 1)
        accel_x = op.nodeAccel(ND, 1)
        drift_x = ((disp_x ) / H1) * 100
        disp_y = op.nodeDisp(ND, 2)
        vel_y = op.nodeVel(ND, 2)
        accel_y = op.nodeAccel(ND, 2)
        rotat = op.nodeDisp(ND, 3)
        # Record results
        op.reactions()
        S = op.nodeReaction(1, 1) + op.nodeReaction(N, 1) # SHEAR BASE REACTION
        A = op.nodeReaction(1, 2) + op.nodeReaction(N, 2) # AXIAL BASE REACTION
        M = op.nodeReaction(1, 3) + op.nodeReaction(N, 3) # MOMENT BASE REACTION
        displacements_x.append(disp_x)
        displacements_y.append(disp_y)
        drift.append(drift_x)
        rotations.append(rotat)
        velocityX.append(vel_x)
        accelerationX.append(accel_x)
        velocityY.append(vel_y)
        accelerationY.append(accel_y)
        base_shears.append(S)
        base_axials.append(A)
        base_moments.append(M)
        stiffness.append(np.abs(base_shears[-1])/np.abs(displacements_x[-1])) # LATERAL STIFFNESS IN X
        # IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
        PERIODmin, PERIODmax = S04.EIGENVALUE_ANALYSIS(3, PLOT=True)
        PERIOD_MIN.append(PERIODmin)
        PERIOD_MAX.append(PERIODmax)
        print('DYN ', step + 1,' Step Done.')

    print('Free-vibration Done.\n')
    #op.wipe()
    # Calculating Damping Ratio and Period Using Logarithmic Decrement Analysis 
    damping_ratio = S03.DAMPING_RATIO(displacements_x)
    
    DATA = (TIME, drift, displacements_x, displacements_y, rotations, velocityX,
            accelerationX, velocityY, accelerationY, base_shears, base_axials,
            base_moments, PERIOD_MIN, PERIOD_MAX, stiffness, damping_ratio)
    return DATA

#%%-------------------------------------------------------------------------------
### --------------------------------------
###        FREE-VIBRATION ANALYSIS
### --------------------------------------

# Define Section Geometry
H1 = 3000.0 # [mm] Column length
H2 = 1000.0 # [mm] Beam height
L1 = 7000.0 # [mm] Beam length

hw_COL = 350.0 # [mm] Section Web Height 
tw_COL = 10.0  # [mm] Section Web Thickness
tf_COL = 10.0  # [mm] Section Flange Thickness
bf_COL = 110.0 # [mm] Section Flange Width

hw_BEAM = 150.0 # [mm] Section Web Hight 
tw_BEAM = 10.0  # [mm] Section Web Thickness
tf_BEAM = 10.0  # [mm] Section Flange Thickness
bf_BEAM = 110.0 # [mm] Section Flange Width

N = 200  # STRUCTURES NODES COUNT
ND = 150 # NODE NUMBER APPLIED FREE-VIBRATION

Weight = 100000.0 # [N] weight on he top of column
ul = -10          # [N/mm] Uniform Distributed Loads

# FREE-VIBRATION
IU = True         # Free Vibration with Initial Displacement
IV = True         # Free Vibration with Initial Velocity
IA = True         # Free Vibration with Initial Acceleration
u0 = 0.10         # [mm] Initial displacement
v0 = 0.0015       # [mm/s] Initial velocity
a0 = 0.00065      # [mm/s^2] Initial acceleration
Duration = 10.0   # [s] Analysis duration
dt = 0.01         # [s] Time step

DAMPING_RATIO = 0.05    # 5% Damping Ratio

### ANALYSIS TOLEANCE AND ITERATIONS
MAX_ITERATIONS = 5000   # Maximum number of iterations
MAX_TOLERANCE = 1.0e-6  # Specified tolerance for convergence

#%%-------------------------------------------------------------------------------
# Analysis Durations (monitor cpu time):
starttime = TI.process_time()

ELE_KIND = 'INELASTIC'    # Elements Material Behavior (INELASTIC or ELASTIC)
DIST_LOAD = True          # UNIFORMLY-DISTRIBURED LOAD
nfCoreY, nfCoreZ = 10, 50 # Number of fibers for steel in Y-direction and Z-direction

DATAd = FREE_VIBRAION_ANALYSIS(nfCoreY, nfCoreZ, MAX_ITERATIONS, MAX_TOLERANCE,
                               u0, v0, a0, IA, IU, IV,
                               hw_COL, tw_COL, tf_COL, bf_COL, hw_BEAM, tw_BEAM, tf_BEAM, bf_BEAM,
                               H1, H2, L1, N, ND, Weight, ul,
                               DAMPING_RATIO, dt, Duration, ELE_KIND, DIST_LOAD)

(TIME, driftE, displacements_xE, displacements_yE, rotationsE, velocityXE,
        accelerationXE, velocityYE, accelerationYE, base_shearsE, base_axialsE,
        base_momentsE, PERIOD_MINE, PERIOD_MAXE, stiffnessE, damping_ratioE) = DATAd


# %% Plot 2D Frame Shapes for Free-vibration Analysis
S02.PLOT_2D_FRAME(deformed_scale=10)  # Adjust scale factor as needed
#%%-----------------------------------------------------------------
ELE_KIND = 'ELASTIC'      # Elements Material Behavior (INELASTIC or ELASTIC)
DIST_LOAD = True          # UNIFORMLY-DISTRIBURED LOAD
nfCoreY, nfCoreZ = 5, 10  # Number of fibers for steel in Y-direction and Z-direction

DATAd = FREE_VIBRAION_ANALYSIS(nfCoreY, nfCoreZ, MAX_ITERATIONS, MAX_TOLERANCE,
                               u0, v0, a0, IA, IU, IV,
                               hw_COL, tw_COL, tf_COL, bf_COL, hw_BEAM, tw_BEAM, tf_BEAM, bf_BEAM,
                               H1, H2, L1, N, ND, Weight, ul,
                               DAMPING_RATIO, dt, Duration, ELE_KIND, DIST_LOAD)

(TIME, driftP, displacements_xP, displacements_yP, rotationsP, velocityXP,
        accelerationXP, velocityYP, accelerationYP, base_shearsP, base_axialsP,
        base_momentsP, PERIOD_MINP, PERIOD_MAXP, stiffnessP, damping_ratioP) = DATAd


# %% Plot 2D Frame Shapes for Free-vibration Analysis
S02.PLOT_2D_FRAME(deformed_scale=10)  # Adjust scale factor as needed

totaltime = TI.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n') 
#%%-------------------------------------------------------------------------------
# PLOT THE TIME-HISTORY
PLOT_TIME_HISORY(TIME, displacements_xE, displacements_yE, rotationsE, velocityXE, accelerationXE, velocityYE, accelerationYE, base_shearsE, base_axialsE, base_momentsE)
PLOT_TIME_HISORY(TIME, displacements_xP, displacements_yP, rotationsP, velocityXP, accelerationXP, velocityYP, accelerationYP, base_shearsP, base_axialsP, base_momentsP)
#%%-------------------------------------------------------------------------------
# PLOT STRUCTURAL PERIOD DURING THE ANALYSIS
plt.figure(0, figsize=(12, 8))
plt.plot(TIME, PERIOD_MINE)
plt.plot(TIME, PERIOD_MAXE)
plt.plot(TIME, PERIOD_MINP)
plt.plot(TIME, PERIOD_MAXP)
plt.title('Period of Structure')
plt.ylabel('Structural Period [s]')
plt.xlabel('Time [s]')
#plt.semilogy()
plt.grid()
plt.legend([f'ELASTIC PERIOD - MIN VALUES: Min: {np.min(PERIOD_MINE):.3f} (s) - Mean: {np.mean(PERIOD_MINE):.3f} (s) - Max: {np.max(PERIOD_MINE):.3f} (s)', 
            f'ELASTIC PERIOD - MAX VALUES:  Min: {np.min(PERIOD_MAXE):.3f} (s) - Mean: {np.mean(PERIOD_MAXE):.3f} (s) - Max: {np.max(PERIOD_MAXE):.3f} (s)',
            f'INELASTIC PERIOD - MIN VALUES: Min: {np.min(PERIOD_MINP):.3f} (s) - Mean: {np.mean(PERIOD_MINP):.3f} (s) - Max: {np.max(PERIOD_MINP):.3f} (s)', 
            f'INELASTIC PERIOD - MAX VALUES:  Min: {np.min(PERIOD_MAXP):.3f} (s) - Mean: {np.mean(PERIOD_MAXP):.3f} (s) - Max: {np.max(PERIOD_MAXP):.3f} (s)',
            ])
plt.show()
#%%-------------------------------------------------------------------------------
# Plot Results
plt.figure(2, figsize=(18, 14))

# Define colors
elastic_color = 'black'
inelastic_color = 'red'

# Reaction plot
plt.subplot(5, 1, 1)
plt.plot(TIME, base_shearsE, color=elastic_color, linewidth=1.5, label=f'Elastic Damping Ratio: {100*damping_ratioE:.3e} %')
plt.plot(TIME, base_shearsP, color=inelastic_color, linewidth=1.5, label=f'Inelastic Damping Ratio: {100*damping_ratioP:.3e} %')
plt.title('Reaction Forces vs Time', fontsize=12, pad=10)
plt.ylabel('Reaction (N)', fontsize=10)
plt.grid(alpha=0.3)
plt.legend(loc='upper right', framealpha=1)

# Displacement plot
plt.subplot(5, 1, 2)
plt.plot(TIME, displacements_xE, color=elastic_color, linewidth=1.5)
plt.plot(TIME, displacements_xP, color=inelastic_color, linewidth=1.5)
plt.title('Displacement vs Time', fontsize=12, pad=10)
plt.ylabel('Displacement (mm)', fontsize=10)
plt.grid(alpha=0.3)

# Velocity plot
plt.subplot(5, 1, 3)
plt.plot(TIME, velocityXE, color=elastic_color, linewidth=1.5)
plt.plot(TIME, velocityXP, color=inelastic_color, linewidth=1.5)
plt.title('Velocity vs Time', fontsize=12, pad=10)
plt.ylabel('Velocity (mm/s)', fontsize=10)
plt.grid(alpha=0.3)

# Acceleration plot
plt.subplot(5, 1, 4)
plt.plot(TIME, accelerationXE, color=elastic_color, linewidth=1.5)
plt.plot(TIME, accelerationXP, color=inelastic_color, linewidth=1.5)
plt.title('Acceleration vs Time', fontsize=12, pad=10)
plt.ylabel('Acceleration (mm/s²)', fontsize=10)
plt.grid(alpha=0.3)

# Stiffness plot
plt.subplot(5, 1, 5)
plt.plot(TIME, stiffnessE, color=elastic_color, linewidth=1.5)
plt.plot(TIME, stiffnessP, color=inelastic_color, linewidth=1.5)
plt.title('Stiffness vs Time', fontsize=12, pad=10)
plt.ylabel('Stiffness (N/mm)', fontsize=10)
plt.semilogy()
plt.grid(alpha=0.3)

plt.tight_layout()
plt.show()

plt.figure(3, figsize=(8, 6))
plt.plot(displacements_xE, base_shearsE, color='black', linewidth=2)
plt.plot(displacements_xP, base_shearsP, color='purple', linewidth=2)
plt.xlabel('Displacement [mm]')
plt.ylabel('Base-Shear Reaction [N]')
plt.title('Displacement vs Base-Shear Reaction')
plt.legend(['INELASTIC', 'ELASTIC'])
plt.grid()
plt.show()
#%%-------------------------------------------------------------------------------