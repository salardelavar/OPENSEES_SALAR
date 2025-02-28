###########################################################################################################
#                                                IN THE NAME OF ALLAH                                     #
#            MOMENT-CURVATURE ANALYSIS OF CONFINED CONCRETE SECTION IN UNCERTAINITY CONDITIONS            #
#---------------------------------------------------------------------------------------------------------#
#                          THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                     #
#                                   EMAIL: salar.d.ghashghaei@gmail.com                                   #
########################################################################################################### 

import random
import time as TI
import numpy as np
import matplotlib.pyplot as plt
import openseespy.opensees as ops
import Analysis_Function as S02

# --------------------------------------------------------------------------------------------------

# Function to perform Moment-Curvature analysis
def MOMENT_CURVATURE(secTag):
    # Define Analysis Properties
    MAX_ITERATIONS = 1000      # Convergence iteration for test
    MAX_TOLERANCE = 1.0e-10    # Convergence tolerance for test
    
    # Clear previous model
    ops.wipe()
    
    # Define model builder
    ops.model('basic', '-ndm', 2, '-ndf', 3)

    # Define materials for nonlinear columns
    # ------------------------------------------
    # CONCRETE                  tag   f'c        ec0   f'cu        ecu
    # Core concrete (confined)
    fcC = random.uniform(-27, -23)
    ec0C = random.uniform(-0.0045, -0.0035)
    fcUC = random.uniform(-21, -20)
    ecuC = random.uniform(-0.015, -0.014)
    ops.uniaxialMaterial('Concrete01', 1, fcC, ec0C, fcUC, ecuC)

    # Cover concrete (unconfined)
    fcU = random.uniform(-18, -15)
    ec0U = random.uniform(-0.0025, -0.0020) 
    fcUU = random.uniform(-2, -1)
    ecuU = random.uniform(-0.008, -0.007)
    ops.uniaxialMaterial('Concrete01', 2, fcU, ec0U, fcUU, ecuU)

    # STEEL
    # Reinforcing steel
    fy = random.uniform(3850, 4150)  # Yield stress
    ey = random.uniform(0.018, 0.022)  # Yield strain
    E = fy / ey  # Young's modulus
    ops.uniaxialMaterial('Steel01', 3, fy, E, 0.01)

    # Define cross-section for nonlinear columns
    # ------------------------------------------

    # Set some parameters
    colWidth = random.uniform(395, 405)
    colDepth = random.uniform(595, 605)

    cover = random.uniform(45, 55)    # concrete section cover
    RZD = random.uniform(15.9, 16.1)  # Rebar Size Diameter
    As = (np.pi * RZD**2) / 4         # Area of reinforcement bars

    # Some variables derived from the parameters
    y1 = colDepth / 2.0
    z1 = colWidth / 2.0
    NUMFIBERS = 20  # Number of layers for each fiber

    ops.section('Fiber', 1)
    # Create the concrete core fibers
    ops.patch('rect', 1, NUMFIBERS, 5, cover - y1, cover - z1, y1 - cover, z1 - cover)

    # Create the concrete cover fibers (top, bottom, left, right)
    ops.patch('rect', 2, NUMFIBERS, 5, -y1, z1 - cover, y1, z1)
    ops.patch('rect', 2, NUMFIBERS, 5, -y1, -z1, y1, cover - z1)
    ops.patch('rect', 2, NUMFIBERS, 5, -y1, cover - z1, cover - y1, z1 - cover)
    ops.patch('rect', 2, NUMFIBERS, 5, y1 - cover, cover - z1, y1, z1 - cover)

    # Create the reinforcing fibers (left, middle, right)
    ops.layer('straight', 3, 3, As, y1 - cover, z1 - cover, y1 - cover, cover - z1)
    ops.layer('straight', 3, 2, As, 0.0, z1 - cover, 0.0, cover - z1)
    ops.layer('straight', 3, 3, As, cover - y1, z1 - cover, cover - y1, cover - z1)

    # Estimate yield curvature
    # (Assuming no axial load and only top and bottom steel)
    d = colDepth - cover  # d -- from cover to rebar
    epsy = fy / E  # steel yield strain
    Ky = epsy / (0.7 * d)
    mu = 2
    maxK = mu * Ky

    # Print estimate to standard output
    #print(f"Estimated yield curvature: {Ky}")

    # Set axial load
    P = 0#random.uniform(-180, -160)

    mu = 15  # Target ductility for analysis
    numIncr = 100  # Number of analysis increments
        
    # Define two nodes at (0,0)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 0.0, 0.0)

    # Fix all degrees of freedom except axial and bending
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 0, 1, 0)

    # Define element
    ops.element('zeroLengthSection', 1, 1, 2, secTag)

    # Create recorder for curvature and moment
    ops.recorder('Node', '-file', 'CUR.txt', '-time', '-node', 2, '-dof', 3, 'disp')
    ops.recorder('Node', '-file', 'MOM.txt', '-time', '-node', 1, '-dof', 3, 'reaction')

    # Define constant axial load
    ops.timeSeries('Constant', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(2, P, 0.0, 0.0)

    # Define analysis parameters
    ops.integrator('LoadControl', 0.001)
    ops.system('SparseGeneral', '-piv')
    ops.test('NormUnbalance', MAX_TOLERANCE, MAX_ITERATIONS)
    ops.numberer('Plain')
    ops.constraints('Plain')
    ops.algorithm('Newton')
    ops.analysis('Static')

    # Do one analysis for constant axial load
    ops.analyze(1)

    # Define reference moment
    ops.timeSeries('Linear', 2)
    ops.pattern('Plain',2, 2)
    ops.load(2, 0.0, 0.0, 1.0)
   
    # Compute curvature increment
    dK = maxK / numIncr

    # Use displacement control at node 2 for section analysis
    ops.integrator('DisplacementControl', 2, 3, dK, 1, dK, dK)

    # Do the section analysis
    OK = ops.analyze(numIncr)
    S02.ANALYSIS(OK, numIncr, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS

# --------------------------------------------------------------------------------------------------

# Lists to store results for plotting
MAX_CUR = []
MAX_MOM = []

# Start time
starttime = TI.time()

# Number of Monte Carlo iterations
NUM_INCREMENT = 50

# Perform Monte Carlo simulation
for II in range(NUM_INCREMENT):
    # Call the MomentCurvature function
    MOMENT_CURVATURE(1)
    
    # Read results from the recorder file
    Curvatures = np.loadtxt('CUR.txt', usecols=1)  # Second column is curvature
    Moments = np.loadtxt('MOM.txt', usecols=1)     # Second column is moment
    MAX_CUR.append(np.max(np.abs(Curvatures)))
    MAX_MOM.append(np.max(np.abs(Moments)))
    print(f'STEP {II+1} DONE')
else:
     print('Analysis completed successfully')       


totaltime = TI.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n')

# --------------------------------------------------------------------------------------------------

def HISROGRAM_BOXPLOT_MATPLOTLIB(X, HISTO_COLOR, LABEL):
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
    print(f"Box-Chart Datas {LABEL} : ")
    print(f'Minimum: {MINIMUM:.6e}')
    print(f'First quartile: {q1:.6e}')
    #print(f'Mode: {MODE:.6e}')
    print(f'Median: {MEDIAN:.6e}')
    print(f'Mean: {MEAN:.6e}')
    print(f'Std: {STD:.6e}')
    print(f'Third quartile: {q3:.6e}')
    print(f'Maximum: {MAXIMUM :.6e}')
    print(f'Skewness: {skew(X) :.6e}')
    print(f'kurtosis: {kurtosis(X) :.6e}')
    print(f"90% Confidence Interval: ({lower_bound:.6e}, {upper_bound:.6e})")
    print("-------------------------")

    plt.figure(figsize=(10,6))
    # Plot histogram of data
    count, bins, ignored = plt.hist(X, bins=100, color=HISTO_COLOR, density=True, align='mid')#, edgecolor="black"
    
    # Plot lognormal PDF
    x = np.linspace(min(bins), max(bins), 1000)
    pdf = (np.exp(-(x - MEAN)**2 / (2 * STD**2)) / (STD * np.sqrt(2 * np.pi)))
    plt.plot(x, pdf, linewidth=2, color='r', label="Normal PDF")
    
    # Plot vertical lines for risk measures
    plt.axvline(q1, color="black", linestyle="--", label=f"Quantile 0.25: {q1:.6e}")
    plt.axvline(MEDIAN, color="green", linestyle="--", label=f"Median: {MEDIAN:.6e}")
    plt.axvline(q3, color="black", linestyle="--", label=f"Quantile 0.75: {q3:.6e}")
    #plt.axvline(MODE, color="purple", linestyle="--", label=f"Mode: {MODE:.6e}")
    plt.axvline(MEAN, color="red", linestyle="--", label=f"Mean: {MEAN:.6e}")
    plt.axvline(MEAN-STD, color="blue", linestyle="--", label=f"Mean-Std: {MEAN-STD:.6e}")
    plt.axvline(MEAN+STD, color="blue", linestyle="--", label=f"Mean+Std: {MEAN+STD:.6e}")
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
    plt.text(q1, 1.05, f" Q1: {q1:.6e}")
    plt.text(MEDIAN, 1.1, f" Q2: {MEDIAN:.6e}")
    plt.text(q3, 1.05, f" Q3: {q3:.6e}")
    #plt.text(MODE, 1.15, f" Mode: {MODE:.6e}")
    
    #plt.text(MEAN, 0.9, f" Mean: {MEAN:.6e}")
    #plt.text(MEAN-STD, 0.9, f" Mean-Std: {MEAN-STD:.6e}")
    #plt.text(MEAN+STD, 0.9, f" Mean+Std: {MEAN+STD:.6e}")
    plt.scatter(MEAN, 1, color="red", marker="+", s=200, label=f"Mean: {MEAN:.6e}")
    plt.scatter(MEAN-STD, 1, color="green", marker="X", s=200, label=f"Mean-Std: {MEAN-STD:.6e}")
    plt.scatter(MEAN+STD, 1, color="blue", marker="*", s=200, label=f"Mean+Std:  {MEAN+STD:.6e}")
    plt.xlabel(LABEL)
    plt.ylabel("Data")
    plt.title(f"Boxplot of {LABEL}")
    plt.legend()
    plt.grid()
    plt.show()

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
    print('+--------------------------------------------------------------------+')
    print(f' Elastic Flextural Rigidity :             {Elastic_EI:.2f}')
    print(f' Plastic Flextural Rigidity :             {Plastic_EI:.2f}')
    print(f' Tangent Flextural Rigidity :             {Tangent_EI:.2f}')
    print(f' Section Ductility Ratio :                {Ductility_Rito:.2f}')
    print(f' Section Over Stregth Factor:             {Over_Stregth_Factor:.2f}')
    print('+--------------------------------------------------------------------+')
    """
    """
    # PUSHOVER ANALYSIS
    print('+==========================+')
    print('=   Analysis curve fitted =')
    print('     Disp       Baser Shear')
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
        
# --------------------------------------------------------------------------------------------------

# ------------------------------------------
#  Plot Moment-Cuvature Section Analysis
# ------------------------------------------
xxc, yyc, _, _, _, _, _ = BILNEAR_CURVE(Curvatures, -Moments, 10)
xxc = np.abs(xxc); yyc = np.abs(yyc); # ABSOLUTE VALUE
XLABEL = 'Curvature'
YLABEL = 'Moment'
LEGEND01 = 'Curve'
LEGEND02 = 'Bilinear Fitted'
TITLE = 'Moment-Curvature Analysis'
COLOR = 'black'
PLOT_2D(Curvatures, -Moments, xxc, yyc, _, _, XLABEL, YLABEL, TITLE, LEGEND01, LEGEND02, _, COLOR='black', Z=2) 
print(f'\t\t Ductility Ratio: {yyc[2]/yyc[1]:.4f}')
# --------------------------------------------------------------------------------------------------
X = MAX_CUR
HISTO_COLOR = 'cyan'
LABEL = 'Maximum Curvature (K)'
HISROGRAM_BOXPLOT_MATPLOTLIB(X, HISTO_COLOR, LABEL)
# --------------------------------------------------------------------------------------------------
X = MAX_MOM
HISTO_COLOR = 'lime'
LABEL = 'Maximum Moments'
HISROGRAM_BOXPLOT_MATPLOTLIB(X, HISTO_COLOR, LABEL)
# --------------------------------------------------------------------------------------------------


