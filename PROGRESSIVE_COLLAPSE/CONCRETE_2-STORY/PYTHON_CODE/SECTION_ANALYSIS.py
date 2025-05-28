
#              #########################################################################
#              #                           IN THE NAME OF ALLAH                        #
#              #                    CONFINED CONCRETE SECTION ANALYSIS                 #
#              #      DOF[1]: CONFINED CONCRETE SECTION AXIAL FORCE-STRAIN ANALYSIS    #
#              #      DOF[3]: CONFINED CONCRETE SECTION MOMENT-CURVATURE ANALYSIS      #
#              #-----------------------------------------------------------------------#
#              #       THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)      #
#              #                   EMAIL: salar.d.ghashghaei@gmail.com                 #
#              #########################################################################
#%%-----------------------------------------------------------------------------------
#import the os module
import os
import time
import numpy as np
import openseespy.opensees as ops
import opsvis as opsv
import matplotlib.pyplot as plt
#%%-----------------------------------------------------------------------------------
# Create a directory at specified path with name 'directory_path'
import os
directory_path = 'C:\\OPENSEESPY_SALAR'

# Check if the directory already exists
if not os.path.exists(directory_path):
    os.mkdir(directory_path)
    print(f"Directory '{directory_path}' created successfully.")
else:
    print(f"Directory '{directory_path}' already exists. Skipping creation.")
#%%-----------------------------------------------------------------------------------
# OUTPUT DATA ADDRESS:
FOLDER_NAME = 'CONCRETE_SECTION_ANALYSIS'
SALAR_DIR = f'C://OPENSEESPY_SALAR//{FOLDER_NAME}//';  
#%%-----------------------------------------------------------------------------------
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
# -----------------------------------------------
def Check_Number(NUM):
    import sys as ss
    if NUM == 1 or NUM == 3:
        return True
    else:
        print('\t\t ERROR: DOF MUST BE JUST 1 OR 3')
        return ss.exit()
# -----------------------------------------------
def OUTPUT_SECOND_COLUMN(FOLDER, X, COLUMN, I, Z):
    import numpy as np
    # Time History
    if Z == 1:
        filename = f"C:\\OPENSEESPY_SALAR\\{FOLDER}\\{X}.txt"
        data_collected = np.loadtxt(filename)
        X = data_collected[:, COLUMN]
    if Z == 2:
        filename = f"C:\\OPENSEESPY_SALAR\\{FOLDER}\\{X}_{I}.txt"
        data_collected = np.loadtxt(filename)
        X = data_collected[:, COLUMN]    
    return X 
# -----------------------------------------------
def BILNEAR_CURVE(Cur, Mom, SLOPE_NODE):
    import numpy as np
    # bilinear fitting
    SIZE = len(Mom)
    hh = np.zeros(SIZE-1)
    Aa = np.zeros(SIZE-1)
    for i in range(SIZE-1):
        hh[i] = Cur[i+1] - Cur[i]
        Aa[i] = (Mom[i] + Mom[i+1]) * 0.5 * hh[i]

    Area = sum(Aa)
    k0 = Mom[SLOPE_NODE] / Cur[SLOPE_NODE]
    fiy = (Mom[i+1] * max(Cur) * 0.5 - Area) / (Mom[i+1] * 0.5 - k0 * max(Cur) * 0.5)
    My = k0 * fiy
    X = np.array([0, fiy, max(Cur)])
    Y = np.array([0, My, Mom[i+1]])
    
    # EI and Ductility_Rito
    Elastic_ST = Y[1] / X[1]
    Plastic_ST = Y[2] / X[2]
    Tangent_ST = (Y[2] - Y[1]) / (X[2] - X[1])
    Ductility_Rito = X[2] / X[1]
    Over_Strength_Factor = Y[2] / Y[1]
    """
    # MOMENT-CURVAVTURE ANALYSIS
    print('+==========================+')
    print('=   Analysis curve fitted =')
    print('  Curvature    Moment')
    print('----------------------------')
    print(np.column_stack((X.T, Y.T)))
    print('+==========================+')
    """
    print('+--------------------------------------------------------------------+')
    print(f' Elastic Flextural Rigidity :             {Elastic_ST:.2f}')
    print(f' Plastic Flextural Rigidity :             {Plastic_ST:.2f}')
    print(f' Tangent Flextural Rigidity :             {Tangent_ST:.2f}')
    print(f' Section Ductility Ratio :                {Ductility_Rito:.2f}')
    print(f' Section Over Strength Factor:            {Over_Strength_Factor:.2f}')
    print('+--------------------------------------------------------------------+')
    
    """
    # PUSHOVER ANALYSIS
    print('+==========================+')
    print('=   Analysis curve fitted =')
    print('     Disp       Base Shear ')
    print('----------------------------')
    print(np.column_stack((X.T, Y.T)))
    print('+==========================+')
    print('+----------------------------------------------------+')
    print(f' Structure Elastic Stiffness :     {Elastic_ST:.2f}')
    print(f' Structure Plastic Stiffness :     {Plastic_ST:.2f}')
    print(f' Structure Tangent Stiffness :     {Tangent_ST:.2f}')
    print(f' Structure Ductility Ratio :       {Ductility_Rito:.2f}')
    print(f' Structure Over Strength Factor:   {Over_Strength_Factor:.2f}')
    print('+----------------------------------------------------+')
    """
    return X, Y, Elastic_ST, Plastic_ST, Tangent_ST, Ductility_Rito, Over_Strength_Factor
# -----------------------------------------------
def PLOT_2D(X, Y, Xfit, Yfit, X2, Y2, XLABEL, YLABEL, TITLE, LEGEND01, LEGEND02, LEGEND03, COLOR, Z):
    import matplotlib.pyplot as plt
    #plt.figure(figsize=(12, 8))
    if Z == 1:
        # Plot 1 line
        plt.plot(X, Y,color=COLOR)
        plt.xlabel(XLABEL)
        plt.ylabel(YLABEL)
        plt.title(TITLE)
        plt.grid(True)
        plt.show()
    if Z == 2:
        # Plot 2 lines
        plt.plot(X, Y, Xfit, Yfit, 'r--', linewidth=3)
        plt.title(TITLE)
        plt.xlabel(XLABEL)
        plt.ylabel(YLABEL)
        plt.legend([LEGEND01, LEGEND02], loc='lower right')
        plt.grid(True)
        plt.show()
    if Z == 3:
        # Plot 3 lines
        plt.plot(X, Y, Xfit, Yfit, 'r--', X2, Y2, 'g-*', linewidth=3)
        plt.title(TITLE)
        plt.xlabel(XLABEL)
        plt.ylabel(YLABEL)
        plt.legend([LEGEND01, LEGEND02, LEGEND03], loc='lower right')
        plt.grid(True)
        plt.show() 
def CURRENT_TIME():
    import time
    t = time.localtime()
    current_time = time.strftime("%H:%M:%S", t)
    print(f"Current time (HH:MM:SS): {current_time}\n\n")
        
# -----------------------------------------------
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
             
# -----------------------------------------------
import openseespy.opensees as op    
import opsvis as opsv
#import SECTION_ANALYSIS_FUN as S01
import CONCRETE_SECTION_FUN as S02
### -------------------------------
###    SECTION ANALYSIS FUNCTION
### -------------------------------
def SECTION_ANALYSIS(axial, shear, moment, DR, STRAIN_ULT, DOF):
    Check_Number(DOF)
    op.wipe()
    op.model('basic','-ndm',2,'-ndf',3)
    
    # Define fiber section for I-section
    Bcol, Hcol, Bbeam, Hbeam = 600, 600, 600, 400; # [mm] Column & Beam Section Diamenstion Properties
    COVER = 50        # [mm] Concrete Cover
    REBAR_DIA = 25    # [mm] Steel Rebar Diameter
    
    NUM_INCR = 4000 # Number of analysis increments
    Yield_Strain = 0.002 # Yield Concrete Strain
    MAX_ITERATIONS = 5000# Maximum number of iterations
    TOLERANCE = 1.0e-12# Specified tolerance for convergence
    nFibCoverZ, nFibCoverY = 3, 120
    nFibCoreZ, nFibCoreY = 3, 120
    SECTION_01, SECTION_02 = S02.CONCRETE_SECTION(Hcol, Bcol, Hbeam, Bbeam, COVER, REBAR_DIA, 
                                                   nFibCoverZ, nFibCoverY, nFibCoreZ, nFibCoreY, PLOT=0)
    
    # ----------------
    opsv.fib_sec_list_to_cmds(SECTION_01) # SECTION 01 -> COLUMN SECTION
    #opsv.fib_sec_list_to_cmds(SECTION_02) # SECTION 02 -> BEAM SECTION
    SecTag = 1 # SECTION 1 SELECTED FOR ANALYSIS
    # ----------------
    
    # Yield Curvature
    if DOF == 1:
        Ky = Yield_Strain
    if DOF == 3:
        Ky = Yield_Strain / (0.5 * Hcol)    
    #print('Ky', Ky)
    
    # Ultimate Curvature
    Ku = Ky * DR
    #print('Ku', Ku)
    # Define two nodes at (0,0)
    op.node(1, 0.0, 0.0)
    op.node(2, 0.0, 0.0)

    # Fix all degrees of freedom except axial and bending
    op.fix(1, 1, 1, 1)
    op.fix(2, 0, 1, 0)  
    # Define element
    #                             tag ndI ndJ  secTag
    op.element('zeroLengthSection',  1, 1, 2, SecTag)
    # Create recorder
    op.recorder('Node', '-file', f"{SALAR_DIR}CUR.txt",'-time', '-node', 2, '-dof', DOF, 'disp')# Curvature Time History nodes 2
    op.recorder('Node', '-file', f"{SALAR_DIR}MOM.txt",'-time', '-node', 1, '-dof', DOF, 'reaction')# Base Shear Time History nodes 1

    # Define constant axial load
    op.timeSeries('Constant', 1)
    op.pattern('Plain', 1, 1)
    op.load(2, axial, shear, moment)
    
    # Define gravity analysis parameters
    NstepGravity = 10
    DGravity = 1 / NstepGravity
    op.integrator('LoadControl', DGravity) # determine the next time step for an analysis
    op.numberer('Plain') # renumber dof's to minimize band-width (optimization), if you want to
    op.system('BandGeneral') # how to store and solve the system of equations in the analysis
    op.constraints('Plain') # how it handles boundary conditions
    op.test('EnergyIncr', TOLERANCE, MAX_ITERATIONS) # determine if convergence has been achieved at the end of an iteration step
    op.algorithm('ModifiedNewton') # use Newton's solution algorithm: updates tangent stiffness at every iteration
    op.analysis('Static') # define type of analysis static or transient
    OK = op.analyze(NstepGravity) # apply gravity
    ANALYSIS(OK, NstepGravity, TOLERANCE, MAX_ITERATIONS)
    print('Load Done.')

    op.loadConst('-time', 0.0) #maintain constant gravity loads and reset time to zero
    
    # Define incremental displacement analysis parameters
    op.integrator('LoadControl', 0.01)
    op.system('SparseGeneral', '-piv')
    op.test('EnergyIncr', TOLERANCE, MAX_ITERATIONS)
    op.numberer('Plain')
    op.constraints('Plain')
    op.algorithm('ModifiedNewton')
    op.analysis('Static')

    # Do one analysis for constant axial load
    OK = op.analyze(1)
    ANALYSIS(OK, 1, TOLERANCE, MAX_ITERATIONS)

    # Define reference moment
    op.timeSeries('Linear', 2)
    op.pattern('Plain',2, 2)
    # Compute curvature increment
    dK = Ku / NUM_INCR
    MOMENT, CUR= [], []
    ops.reactions()
    if DOF == 1:# AXIAL FORCE - AXIAL STRAIN ANALYSIS
        op.load(2, 1.0, 0.0, 0.0)
        # Use displacement control at node 2 for section analysis
        op.integrator('DisplacementControl', 2, 1, dK, 1, dK, dK)
    if DOF == 3:# MOMENT-CURVATURE ANALYSIS
        op.load(2, 0.0, 0.0, 1.0)
        # Use displacement control at node 2 for section analysis
        op.integrator('DisplacementControl', 2, 3, dK, 1, dK, dK)

    Depth = Hcol
    # Do the section analysis
    strain = 0.0
    it = 0
    while strain > STRAIN_ULT:
        #OK = op.analyze(NUM_INCR)
        #ANALYSIS(OK, NUM_INCR, TOLERANCE, MAX_ITERATIONS)
        if DOF == 1:
            OK = op.analyze(1)
            ANALYSIS(OK, 1, TOLERANCE, MAX_ITERATIONS)
            strain = -op.nodeDisp(2,1)
            m = ops.nodeReaction(1, 3) # AXIAL FORCE
            cur = ops.nodeDisp(2, 3)   # STRAIN
            #print(f'ITERTION {it+1} STRAIN {strain} DONE')
        if DOF == 3:
            OK = op.analyze(1)
            ANALYSIS(OK, 1, TOLERANCE, MAX_ITERATIONS)
            strain = op.nodeDisp(2,1) - Depth/2 * op.nodeDisp(2,3)
            m = ops.nodeReaction(1, 3) # MOMENT FORCE
            cur = ops.nodeDisp(2, 3)   # CURVATURE
            #print(f'ITERTION {it+1} STRAIN {strain} DONE')
        it += 1 # update iterations 
        MOMENT.append(m)
        CUR.append(cur)
    print(f'SECTION ANALYSIS DONE.')
    op.wipe()
    return MOMENT, CUR
    


#%%-----------------------------------------------------------------------------------
import CONCRETE_SECTION_FUN as S02

# ------------------------------------
# Axial Force-Strain Section Analysis
# ------------------------------------

# Call the section analysis procedure:
axial = -50.0        # [N] PY - VERTICAL FORCE - DOF[11]
shear =  0.0         # [N] PX - HORIZENTAL FORCE - DOF[10]
moment = 0.0         # [N.mm] MZ - MOMENT FORCE - DOF[12]
DR = 4.5             # [mm/mm] set ductility ratio for axial force-strain
STRAIN_ULT = -0.004  # [mm/mm] ULTIMATE STRAIN (IN HERE NEGATIVE VALUE IS COMPRESSION)
DOF = 1              # AXIAL FORCE - AXIAL STRAIN ANALYSIS
_, _ = SECTION_ANALYSIS(axial, shear, moment, DR, STRAIN_ULT, DOF)

# ------------------------------------------
#  Plot Axial Force-Strain Section Analysis
# ------------------------------------------
STR = OUTPUT_SECOND_COLUMN(FOLDER_NAME,'CUR', 1, 0, 1)
FORCE = OUTPUT_SECOND_COLUMN(FOLDER_NAME,'MOM', 1, 0, 1)
xxc, yyc, _, _, _, _, _ = BILNEAR_CURVE(STR, -FORCE, 10)
xxc = np.abs(xxc); yyc = np.abs(yyc); # ABSOLUTE VALUE
XLABEL = 'Strain'
YLABEL = 'Axial Force'
LEGEND01 = 'Curve'
LEGEND02 = 'Bilinear Fitted'
TITLE = 'Axial Force and Strain Analysis'
COLOR = 'black'
PLOT_2D(STR, -FORCE, xxc, yyc, _, _, XLABEL, YLABEL, TITLE, LEGEND01, LEGEND02, _, COLOR='black', Z=2) 
print(f'\t\t Ductility Ratio: {yyc[2]/yyc[1]:.4f}')
#%%-----------------------------------------------------------------------------------
# ------------------------------------
#   Moment-Curvature Section Analysis
# ------------------------------------

# Call the section analysis procedure:
axial = -50.0           # [N] PY - VERTICAL FORCE - DOF[11]
shear =  0.0            # [N] PX - HORIZENTAL FORCE - DOF[10]
moment = 0.0            # [N.mm] MZ - MOMENT FORCE - DOF[12]
DR = 35.5               # [mm/mm] set ductility ratio for moment curvature
STRAIN_ULT = -0.071     # [mm/mm] ULTIMATE STRAIN
DOF = 3                 # AXIAL FORCE - AXIAL STRAIN ANALYSIS
_, _ = SECTION_ANALYSIS(axial, shear, moment, DR, STRAIN_ULT, DOF)

# -----------------------------------------
#   Plot Moment-Curvature Section Analysis
# -----------------------------------------
CUR = OUTPUT_SECOND_COLUMN(FOLDER_NAME,'CUR', 1, 0, 1)
MOM = OUTPUT_SECOND_COLUMN(FOLDER_NAME,'MOM', 1, 0, 1)
xxc, yyc, _, _, _, _, _ = BILNEAR_CURVE(CUR, -MOM, 2)
xxc = np.abs(xxc); yyc = np.abs(yyc); # ABSOLUTE VALUE
XLABEL = 'Curvature'
YLABEL = 'Moment'
LEGEND01 = 'Curve'
LEGEND02 = 'Bilinear Fitted'
TITLE = 'Moment and Curvature Analysis'
COLOR = 'black'
PLOT_2D(CUR, -MOM, xxc, yyc, _, _, XLABEL, YLABEL, TITLE, LEGEND01, LEGEND02, _, COLOR='black', Z=2) 
print(f'\t\t Ductility Ratio: {yyc[2]/yyc[1]:.4f}')
#%%-----------------------------------------------------------------------------------    
