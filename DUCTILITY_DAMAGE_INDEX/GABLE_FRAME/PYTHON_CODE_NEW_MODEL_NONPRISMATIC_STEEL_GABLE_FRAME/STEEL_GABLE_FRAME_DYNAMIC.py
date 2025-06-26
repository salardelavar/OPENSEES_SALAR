#  #########################################################################
#  #                           IN THE NAME OF ALLAH                        #
#  #     NONLINEAR DYNAMIC ANALYSIS OF A NONPRISMATIC STEEL GABLE FRAME    #
#  #              I SECTION COLUMN WITH FINITE PRISMATIC COLUMN            #
#  #-----------------------------------------------------------------------#
#  #    MODELING OF NONPRISMATIC ELEMENT WITH MULTI PRISMATIC ELEMENTS     #
#  #-----------------------------------------------------------------------#
#  #              THIS PROGRAM WRITTEN BY SALAR DELAVAR QASHQAI            #
#  #                   EMAIL: salar.d.ghashghaei@gmail.com                 #
#  #########################################################################

#%%-------------------------------------------------------------
#import the os module
import os
import time as TI
import numpy as np
import openseespy.opensees as op
import ANALYSIS_FUNCTION as S03
import PLOT_2D as S04
#%%-------------------------------------------------------------
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

# Load the image
image_path = 'STEEL_GABLE_FRAME.PNG'
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
def DYNAMIC_ANALYSIS(hw_COL, tw_COL,tf_COL, bf_COL, hw_BEAM, tw_BEAM, tf_BEAM, bf_BEAM, H1, H2, L1, N, ND, Weight, UL, Damping_Ratio, dt, Duration, ELE_KIND):
    ### ANALYSIS TOLEANCE AND ITERATIONS
    MAX_ITERATIONS = 1000   # Maximum number of iterations
    MAX_TOLERANCE = 1.0e-8  # Specified tolerance for convergence
    
    NN = int(N / 4) # 4 ELEMENTS
    
    # Define the model
    op.wipe()
    op.model('basic', '-ndm', 2, '-ndf', 3) 
    PCol = Weight    # nodal dead-load weight per column
    GMfact = 9810    # [mm/s^2] standard acceleration of gravity or standard acceleration
    Mass =  PCol/GMfact
    h1 = H1 / (NN - 1)
    h2 = H2 / (NN - 1)
    bl = (L1 **2 + H2 **2)**0.5  # Beam Length
    l1 = L1 / (NN) 
          
    IDctrlNode = ND ## INCREMENTAL DISPLACEMENT NODE
    IDctrlDOF = 1
    
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
    Es = 210e4        # [N/mm^2] Young's modulus 
    fy = 240          # [N/mm^2] Yield strength
    ey = fy/Es        # [mm/mm] Steel Yield Strain
    fu = 1.1818*fy    # [N/mm²] Steel Ultimate Strength
    esu = ey*75.2     # [mm/mm] Steel Ultimate Strain
    pinchX = 0.8   # Pinching factor in X direction
    pinchY = 0.5   # Pinching factor in Y direction
    damage1 = 0.0  # Damage due to ductility
    damage2 = 0.0  # Damage due to energy
    beta = 0.1     # Stiffness degradation parameter
    if ELE_KIND == 'INELASTIC':
        op.uniaxialMaterial('Hysteretic', Mat_Tag, fy, ey, fu, esu, 0.2*fu, 1.1*esu, -fy, -ey, -fu, -0.5*esu, -0.2*fu, -1.1*esu, pinchX, pinchY, damage1, damage2, beta)
        # INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Hysteretic_Material
    if ELE_KIND == 'ELASTIC':
        op.uniaxialMaterial('Elastic', Mat_Tag, Es)
        # INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Elastic_Uniaxial_Material
    # ---------------------------------------------------------------------------------------
    hwbot_COL = 0.5 * hw_COL # Bottom column web height
    # nodal coordinates:
    op.node(1, 0.0, 0.0) # node#, X, Y
    op.node(N, L1 * 2, 0.0) # node#, X, Y
    # Single point constraints -- Boundary Conditions
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

        nfCoreY = 15;			# number of fibers for steel in y-direction
        nfCoreZ = 5;			# number of fibers for steel in z-direction


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

        nfCoreY = 15;			# number of fibers for steel in y-direction
        nfCoreZ = 5;			# number of fibers for steel in z-direction


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

        nfCoreY = 15;			# number of fibers for steel in y-direction
        nfCoreZ = 5;			# number of fibers for steel in z-direction


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
        HW_BEAM =  hw_BEAM + ((hwbot_BEAM - hw_BEAM) / L1) * l1 * zz # Varying web height of nonprismatic section
        #print(i, x, y, i+1, i, HW_BEAM)
        coverY = (HW_BEAM + tf_BEAM) / 2.0
        coverZ = tw_BEAM / 2.0
        coreY = coverY - tf_BEAM
        coreZ02 = bf_BEAM / 2.0

        nfCoreY = 15;			# number of fibers for steel in y-direction
        nfCoreZ = 5;			# number of fibers for steel in z-direction

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
    for i in range(NN + 1, 2 * NN + 1, 1): # CREATE UNIFORM LOADS FOR BEAM 01
        op.eleLoad('-ele', i-1,'-type', '-beamUniform', UL, 0.0) # uniformly-distributed load
    for i in range(3 * NN, 2 * NN, -1):# CREATE UNIFORM LOADS FOR BEAM 02
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
    
    # Calculate Rayleigh damping factors
    Lambda01 = op.eigen('-fullGenLapack', 2)  # eigenvalue mode 2
    #Lambda01 = op.eigen('-genBandArpack', 2) # eigenvalue mode 2
    Omega01 = np.power(max(Lambda01), 0.5)
    Omega02 = np.power(min(Lambda01), 0.5)
    a0 = (2 * Omega01 * Omega02 * Damping_Ratio) / (Omega01 + Omega02) # c = a0 * m : Mass-proportional damping
    a1 = (Damping_Ratio * 2) / (Omega01 + Omega02)   # c = a1 * k : Stiffness-proportional damping
    # Apply Rayleigh damping
    op.rayleigh(a0, a1, 0, 0)   # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/reyleigh.html
    #op.rayleigh(0, 0, 2 * DR * Omega01, 0) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/reyleigh.html
    PERIOD_01 = np.pi / Omega01 # Structure First Period
    PERIOD_02 = np.pi / Omega02 # Structure Second Period
    print('Structure First Period:  ', PERIOD_01)
    print('Structure Second Period: ', PERIOD_02) 
    
    #print('Period Anlysis Done.')
    
    #applying Dynamic Ground motion analysis
    SSF_X = 1.00     # Seismic Acceleration Scale Factor in X Direction
    SSF_Y = 1.00     # Seismic Acceleration Scale Factor in Y Direction
    iv0_X = 0.00005  # [mm/s] Initial velocity applied to the node  in X Direction
    iv0_Y = 0.00001  # [mm/s] Initial velocity applied to the node  in Y Direction
    st_iv0 = 0.0     # [s] Initial velocity applied starting time
    
    SEISMIC_TAG_01 = 100
    op.timeSeries('Path', SEISMIC_TAG_01, '-dt', dt, '-filePath', 'Ground_Acceleration_X.txt', '-factor', GMfact, '-startTime', st_iv0) # SEISMIC-X
    # Define load patterns
    # pattern UniformExcitation $patternTag $dof -accel $tsTag <-vel0 $vel0> <-fact $cFact>
    op.pattern('UniformExcitation', SEISMIC_TAG_01, 1, '-accel', SEISMIC_TAG_01, '-vel0', iv0_X, '-fact', SSF_X) # SEISMIC-X 
    SEISMIC_TAG_02 = 200
    op.timeSeries('Path', SEISMIC_TAG_02, '-dt', dt, '-filePath', 'Ground_Acceleration_Y.txt', '-factor', GMfact) # SEISMIC-Z
    op.pattern('UniformExcitation', SEISMIC_TAG_02, 2, '-accel', SEISMIC_TAG_02, '-vel0', iv0_Y, '-fact', SSF_Y)  # SEISMIC-Z 
    
    #print('Seismic Acceleration Defined.')
    
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
    
    #OK = op.analyze(n_steps, dt)
    #S03.ANALYSIS(OK, Nsteps, MAX_TOLERANCE, MAX_ITERATIONS)
    
    for step in range(n_steps):
        OK = op.analyze(1, dt)
        S03.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS)
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
        print('DYN ', step + 1,' Step Done.')

    print('Dynamic Done.\n')
    #op.wipe()
    # Calculating Damping Ratio Using Logarithmic Decrement Analysis 
    displacement = np.array(displacements_x)
    peaks = np.array([displacement[i] for i in range(1, len(displacement)-1) if displacement[i] > displacement[i-1] and displacement[i] > displacement[i+1]])
    # Natural logarithm
    delta = np.log(peaks[:-1] / peaks[1:])
    return TIME, drift, displacements_x, displacements_y, rotations, velocityX, accelerationX, velocityY, accelerationY, base_shears, base_axials, base_moments, PERIOD_01, PERIOD_02, delta

#%%-------------------------------------------------------------
### --------------------------------------
###          DYNAMIC ANALYSIS
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
ND = 150 # NODE NUMBER APPLIED INCREMENTAL DISPLACEMENT

Weight = 100000.0 # [N] weight on he top of column
ul = -10          #[N/mm] Uniform Distributed Loads

DAMPING_RATIO = 0.02    # 2% damping ratio
dt = 0.01               # [s] Time increment
Duration = 10.0         # [s] Total time duration 


ELE_KIND = 'INELASTIC' # Elements Material Behavior (INELASTIC or ELASTIC)

# Analysis Durations (monitor cpu time):
starttime = TI.process_time()

DATAd = DYNAMIC_ANALYSIS(hw_COL, tw_COL, tf_COL, bf_COL, hw_BEAM, tw_BEAM, tf_BEAM, bf_BEAM, H1, H2, L1, N, ND, Weight, ul, DAMPING_RATIO, dt, Duration, ELE_KIND)
TIME, drift, displacements_x, displacements_y, rotations, velocityX, accelerationX, velocityY, accelerationY, base_shears, base_axials, base_moments, PERIOD_01, PERIOD_02, delta = DATAd

totaltime = TI.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n') 
#%%-------------------------------------------------------------
PLOT_TIME_HISORY(TIME, displacements_x, displacements_y, rotations, velocityX, accelerationX, velocityY, accelerationY, base_shears, base_axials, base_moments)
#%%-------------------------------------------------------------
X = displacements_x
Y = base_shears
XLABEL = 'DISPLACEMENT [mm] DOF(448)' 
YLABEL = 'BASE SHEAR DOF [N]  DOF(1)+(200)' 
TITLE = 'Base-shear & Displacement Diagrarm' 
COLOR = 'Black'
PLOT_2D(X, Y, XLABEL, YLABEL, TITLE, COLOR)
#%%-------------------------------------------------------------
X = displacements_y
Y = base_axials
XLABEL = 'Displacement [mm]' 
YLABEL = 'Base-axial [N]' 
TITLE = 'Base-axial & Displacement Diagrarm' 
COLOR = 'Brown'
PLOT_2D(X, Y, XLABEL, YLABEL, TITLE, COLOR)
#%%-------------------------------------------------------------
X = rotations
Y = base_moments
XLABEL = 'Rotation [Rad]' 
YLABEL = 'Base-moment [N.mm]' 
TITLE = 'Base-moment & Rotation Diagrarm' 
COLOR = 'Purple'
PLOT_2D(X, Y, XLABEL, YLABEL, TITLE, COLOR)
#%%-------------------------------------------------------------
X = drift
Y = base_shears
XLABEL = 'Horizental Drift [%]' 
YLABEL = 'Base-shear [N] DOF(1)' 
TITLE = 'Base-shear & Horizental Drift Diagrarm' 
COLOR = 'red'
PLOT_2D(X, Y, XLABEL, YLABEL, TITLE, COLOR)
#%%-------------------------------------------------------------
# EXACT SOLUTION OF DAMPING RATIO DURING SEISMIC ANALYSIS:
from scipy.optimize import fsolve

# Define the equation for natural logarithm of this ratio, called the logarithmic decrement, we denote by δ
def EQUATION(x, delta):
    if np.any(x == 0):  # Avoid division by zero
        return np.inf  
        
    # Calculate the value of the equation
    A = x**2 - 1 + ((2 * np.pi * x) / np.mean(delta)) ** 2
    #print(f"x: {x}, A: {A}")  # Debugging output
    # Return the difference (for root finding)
    return A
      
# Initial guess for root(s)
x0 = 1  # Intial Guess for Damping Ratio
# Solve for x
solution = fsolve(EQUATION, x0, args=(delta))
print(f"Exact Damping Ratio: {solution[0]:.8e}")
#%%-----------------------------------------------------------------
# %% Plot 2D Frame Shapes for Nonlinear Dynamic Analysis
S04.PLOT_2D_FRAME(deformed_scale=100)  # Adjust scale factor as needed
#%%-------------------------------------------------------------


