##############################################################################################################
#                                                   IN THE NAME OF ALLAH                                     #
# MOMENT-CURVATURE ANALYSIS OF STEEL BOLTED CONNECTION SECTION IN UNCERTAINITY CONDITIONS MONTE-CARLO METHOD #
#------------------------------------------------------------------------------------------------------------#
#                             THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                     #
#                                      EMAIL: salar.d.ghashghaei@gmail.com                                   #
############################################################################################################## 
# STEEL CONNCETION IMAGE IS GOTTEN FROM THIS PAPER:
# Paper Author: Jayarajan Puttatt  
# Email: puttatt1@gmail.com  
'https://www.researchgate.net/publication/370766113_Characterization_of_connections_in_the_seismic_assessment_of_steel_moment_resisting_frames'    

import random
import time as TI
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import openseespy.opensees as ops
import Analysis_Function as S02

# --------------------------------------------------------------------------------------------------

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
    print(f' Section Over Strength Factor:            {Over_Strength_Factor:.2f}')
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

    # Define materials for nonlinear bolted steel connection
    # --------------------------------------------------------
    # STEEL MATERIALS
    # Bolt reinforcing steel
    fy = random.uniform(890, 905)           # Yield stress for bolt
    ey = random.uniform(0.00445, 0.004525)  # Yield strain for bolt
    esu = random.uniform(0.1, 0.12)         # Ultimate strain for bolt
    fu = 1.1818 * fy                        # Ultimate stress for bolt
    E = fy / ey                             # Young's modulus for bolt
    Esh = (fu - fy) / (esu - ey)
    b = Esh / E
    
    # STEEL PALTE PROPERTIES
    ops.uniaxialMaterial('Steel01', 2, fy, E, b)     
    
    # Plate reinforcing steel M30 (Class 10.9)
    fyP = random.uniform(235, 245)             # Yield stress for plate
    eyP = random.uniform(0.001175, 0.001225)   # Yield strain for plate
    EP = fyP / eyP                             # Young's modulus for plate
    
    # STEEL BLOT PROPERTIES - INFO LINK: https://opensees.berkeley.edu/wiki/index.php?title=Elastic-Perfectly_Plastic_Material
    ops.uniaxialMaterial('ENT', 1, EP)               



    # Define cross-section for nonlinear bolted steel connection
    # ------------------------------------------------------------

    # Set some parameters
    plateWidth = 200 #random.uniform(395, 405)
    plateDepth = 670 #random.uniform(595, 605)

    cover = random.uniform(70, 75)        # Bolt section cover
    cover02 = random.uniform(130, 140)    # Bolt section in middle layer
    RZD = random.uniform(29.5, 30.5)      # Bolt Size Diameter
    As = (np.pi * RZD**2) / 4             # Area of reinforcement bolts

    # Some variables derived from the parameters
    y1 = plateDepth / 2.0
    z1 = plateWidth / 2.0
    NUMFIBERS = 20  # Number of layers for each fiber

    ops.section('Fiber', 1)
    # Create the plate core fibers
    ops.patch('rect', 2, NUMFIBERS, 5, cover - y1, cover - z1, y1 - cover, z1 - cover)

    # Create the plate cover fibers (top, bottom, left, right)
    ops.patch('rect', 2, NUMFIBERS, 5, -y1, z1 - cover, y1, z1)
    ops.patch('rect', 2, NUMFIBERS, 5, -y1, -z1, y1, cover - z1)
    ops.patch('rect', 2, NUMFIBERS, 5, -y1, cover - z1, cover - y1, z1 - cover)
    ops.patch('rect', 2, NUMFIBERS, 5, y1 - cover, cover - z1, y1, z1 - cover)

    # Create the reinforcing fibers (left, middle, right)
    ops.layer('straight', 1, 2, As, y1 - cover, z1 - cover, y1 - cover, cover - z1)         # LAYER 01 - TOP LAYER
    ops.layer('straight', 1, 2, As, y1 - cover02, z1 - cover02, y1 - cover02, cover02 - z1) # LAYER 02 - TOP MIDDLE
    ops.layer('straight', 1, 2, As, cover02 - y1, z1 - cover02, cover02 - y1, cover - z1)   # LAYER 03 - BOTTOM MIDDLE
    ops.layer('straight', 1, 2, As, cover - y1, z1 - cover, cover - y1, cover - z1)         # LAYER 04 - BOTTOM LAYER

    # Estimate yield curvature
    # (Assuming no axial load and only top and bottom steel)
    d = plateDepth - cover  # d -- from cover to rebar
    epsy = ey      # steel bolt yield strain
    Ky = epsy / (0.5 * d)
    DR = 10        # Target ductility for analysis
    maxK = DR * Ky
    numIncr = 100  # Number of analysis increments

    # Print estimate to standard output
    #print(f"Estimated yield curvature: {Ky}")

    # Set axial load
    P = 0.0 #random.uniform(-1800, -1600)

        
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
Y_CUR = []
Y_MOM = []
ELASTIC_ST = []
PLASTIC_ST = []
TANGENT_ST = []
DUCT_RATIO = []
OVER_STRENGTH = []

# Start time
starttime = TI.time()

# Number of Monte Carlo iterations
NUM_INCREMENT = 10000

# Perform Monte Carlo simulation
for II in range(NUM_INCREMENT):
    # Call the MomentCurvature function
    MOMENT_CURVATURE(1)
    
    # Read results from the recorder file
    Curvatures = np.loadtxt('CUR.txt', usecols=1)  # Second column is curvature
    Moments = np.loadtxt('MOM.txt', usecols=1)     # Second column is moment
    MAX_CUR.append(np.max(np.abs(Curvatures))) # Output Ultimate Curvatue
    MAX_MOM.append(np.max(np.abs(Moments)))    # Output Ultimate Moment
    xxc, yyc, Elastic_ST, Plastic_ST, Tangent_ST, Ductility_Rito, Over_Strength_Factor = BILNEAR_CURVE(Curvatures, -Moments, 10)
    Y_CUR.append(xxc[1]) # Output Yield Curvatue
    Y_MOM.append(yyc[1]) # Output Yield Moment
    ELASTIC_ST.append(Elastic_ST)
    PLASTIC_ST.append(Plastic_ST)
    TANGENT_ST.append(Tangent_ST)
    DUCT_RATIO.append(Ductility_Rito)
    OVER_STRENGTH.append(Over_Strength_Factor)
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
    print(f"Box-Chart Datas ({LABEL}) : ")
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

def PLOT_2D(X, Y, Xfit, Yfit, X2, Y2, XLABEL, YLABEL, TITLE, LEGEND01, LEGEND02, LEGEND03, COLOR, Z):
    import matplotlib.pyplot as plt
    plt.figure(figsize=(12, 8))
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
# EXPORT DATA TO EXCEL

DATA_TOTAL = {
    'Yield_Curvature': Y_CUR,
    'Yield_Moment': Y_MOM,
    'Ultimate_Curvature': MAX_CUR,
    'Ultimate_Moment': MAX_MOM,
    'Elastic_Flextural_Rigidity': ELASTIC_ST,
    'Plastic_Flextural_Rigidity': PLASTIC_ST,
    'Tangent_Flextural_Rigidity': TANGENT_ST,
    'Section_Ductility_Ratio': DUCT_RATIO,
    'Over_Stength_Fctor': OVER_STRENGTH
}
# Convert to DataFrame
results_df = pd.DataFrame(DATA_TOTAL)
# Export the DataFrame to an Excel file
results_df.to_excel('MOMENT_CURVAURE_UNCERTAINTY_STEEL_BOLTED_CONNECTION_RESULTS.xlsx', index=False)
#------------------------------------------------------------------------------------------------  
# ------------------------------------------
#  Plot Moment-Cuvature Section Analysis
# ------------------------------------------
xxc, yyc, _, _, _, _, _ = BILNEAR_CURVE(Curvatures, -Moments, 10)
xxc = np.abs(xxc); yyc = np.abs(yyc); # ABSOLUTE VALUE
XLABEL = 'Curvature'
YLABEL = 'Moment'
LEGEND01 = 'Curve'
LEGEND02 = 'Bilinear Fitted'
TITLE = f'Last Data of Moment-Curvature Analysis - Ductility Ratio: {yyc[2]/yyc[1]:.4f} - Over Strength Factor: {xxc[2]/xxc[1]:.4f}'
COLOR = 'black'
PLOT_2D(Curvatures, -Moments, xxc, yyc, _, _, XLABEL, YLABEL, TITLE, LEGEND01, LEGEND02, _, COLOR='black', Z=2) 
print(f'\t\t Ductility Ratio: {yyc[2]/yyc[1]:.4f}')
# --------------------------------------------------------------------------------------------------
X = Y_CUR
HISTO_COLOR = 'brown'
LABEL = 'Yield Curvature (K)'
HISROGRAM_BOXPLOT_MATPLOTLIB(X, HISTO_COLOR, LABEL)
# --------------------------------------------------------------------------------------------------
X = Y_MOM
HISTO_COLOR = 'gold'
LABEL = 'Yield Moments'
HISROGRAM_BOXPLOT_MATPLOTLIB(X, HISTO_COLOR, LABEL)
# --------------------------------------------------------------------------------------------------
X = MAX_CUR
HISTO_COLOR = 'cyan'
LABEL = 'Ultimate Curvature (K)'
HISROGRAM_BOXPLOT_MATPLOTLIB(X, HISTO_COLOR, LABEL)
# --------------------------------------------------------------------------------------------------
X = MAX_MOM
HISTO_COLOR = 'lime'
LABEL = 'Ultimate Moments'
HISROGRAM_BOXPLOT_MATPLOTLIB(X, HISTO_COLOR, LABEL)
# --------------------------------------------------------------------------------------------------
X = ELASTIC_ST
HISTO_COLOR = 'brown'
LABEL = 'Elastic Flextural Rigidity'
HISROGRAM_BOXPLOT_MATPLOTLIB(X, HISTO_COLOR, LABEL)
# --------------------------------------------------------------------------------------------------
X = PLASTIC_ST
HISTO_COLOR = 'pink'
LABEL = 'Plastic Flextural Rigidity '
HISROGRAM_BOXPLOT_MATPLOTLIB(X, HISTO_COLOR, LABEL)
# --------------------------------------------------------------------------------------------------
X = TANGENT_ST
HISTO_COLOR = 'yellow'
LABEL = 'Tangent Flextural Rigidity'
HISROGRAM_BOXPLOT_MATPLOTLIB(X, HISTO_COLOR, LABEL)
# --------------------------------------------------------------------------------------------------
X = DUCT_RATIO
HISTO_COLOR = 'purple'
LABEL = 'Section Ductility Ratio'
HISROGRAM_BOXPLOT_MATPLOTLIB(X, HISTO_COLOR, LABEL)
# --------------------------------------------------------------------------------------------------
X = OVER_STRENGTH
HISTO_COLOR = 'lightblue'
LABEL = 'Section Over Strength Factor'
HISROGRAM_BOXPLOT_MATPLOTLIB(X, HISTO_COLOR, LABEL)
# --------------------------------------------------------------------------------------------------
def PLOT_HEATMAP(df):
    import matplotlib.pyplot as plt
    import numpy as np

    # Calculate the correlation matrix
    corr_matrix = df.corr()

    # Create the figure and axis
    fig, ax = plt.subplots(figsize=(12, 12))

    # Create the heatmap
    cax = ax.matshow(corr_matrix, cmap='viridis')

    # Add the colorbar
    fig.colorbar(cax)

    # Set axis labels
    ax.set_xticks(np.arange(len(corr_matrix.columns)))
    ax.set_yticks(np.arange(len(corr_matrix.index)))
    ax.set_xticklabels(corr_matrix.columns, rotation=90)
    ax.set_yticklabels(corr_matrix.index)

    # Annotate the heatmap with the correlation values
    for i in range(len(corr_matrix.columns)):
        for j in range(len(corr_matrix.index)):
            ax.text(j, i, f'{corr_matrix.iloc[i, j]:.2f}', 
                    ha='center', va='center', color='white')

    # Set title and layout
    plt.title('Correlation Heatmap', pad=20)
    plt.xlabel('Variable')
    plt.ylabel('Variable')
    
    # Display the plot
    plt.tight_layout()
    plt.show()
    
# PLOT HEATMAP FOR CORRELATION 
PLOT_HEATMAP(results_df)
#------------------------------------------------------------------------------------------------


