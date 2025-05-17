######################################################################################################################
#                                                       IN THE NAME OF ALLAH                                         #
# PROBABILISTIC MOMENT-CURVATURE ANALYSIS OF STEEL-BOLTED CONNECTIONS UNDER UNCERTAINTY USING MONTE CARLO SIMULATION #
#--------------------------------------------------------------------------------------------------------------------#
# THIS CODE PERFORMS A PROBABILISTIC MOMENT-CURVATURE ANALYSIS OF STEEL-BOLTED CONNECTIONS UNDER MATERIAL            #
# UNCERTAINTY USING MONTE CARLO SIMULATIONS. IT QUANTIFIES VARIABILITY IN KEY PARAMETERS LIKE YIELD/ULTIMATE         #
# MOMENTS, CURVATURE, DUCTILITY RATIOS, AND STIFFNESS THROUGH 10,000 RANDOMIZED SIMULATIONS OF BOLT/PLATE MATERIAL   #
# PROPERTIES. THE ANALYSIS LEVERAGES OPENSEESPY FOR NONLINEAR FIBER-SECTION MODELING AND EXTRACTS BILINEAR           #
# APPROXIMATIONS OF MOMENT-CURVATURE RESPONSES. RESULTS ARE STATISTICALLY ANALYZED (HISTOGRAMS, BOXPLOTS             #
# , CORRELATION HEATMAPS) TO ASSESS CONFIDENCE INTERVALS, SKEWNESS, AND PARAMETER SENSITIVITY. THE OUTPUT            #
#  CHARACTERIZES STRUCTURAL RELIABILITY AND PROBABILISTIC BEHAVIOR OF CONNECTIONS FOR UNCERTAINTY-AWARE DESIGN.      #
#--------------------------------------------------------------------------------------------------------------------#
#                                   THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                       #
#                                             EMAIL: salar.d.ghashghaei@gmail.com                                    #
###################################################################################################################### 
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
    # STEEL BOLT MATERIALS
    # Bolt reinforcing steel
    fy = random.uniform(890, 905)           # Yield stress for bolt
    ey = random.uniform(0.00445, 0.004525)  # Yield strain for bolt
    esu = random.uniform(0.1, 0.12)         # Ultimate strain for bolt
    fu = 1.1818 * fy                        # Ultimate stress for bolt
    E = fy / ey                             # Young's modulus for bolt
    Esh = (fu - fy) / (esu - ey)
    b = Esh / E
    
    """
    MatTag = 2
    # uniaxialMaterial MultiLinear $matTag $u1 $f1 $u2 $f2 $u3 $f3 $u4 $f4
    kh = np.array([[ey, fy], [esu, fu], [1.1*esu, 0.2*fu], [1.2*esu, 0.1*fu]])# COMPRESSION & TENSION STRESS-STRAIN RELATION
    ops.uniaxialMaterial('MultiLinear', MatTag, *kh.flatten())  
    INFO LINK: https://opensees.berkeley.edu/wiki/index.php/MultiLinear_Material
    """
    """
    MatTag = 2
    KP = np.array([[fy , ey], [fu, esu], [0.2*fu, 1.1*esu], [0.1*fu , 1.2*esu]])           # TENSION STRESS-STRAIN RELATION
    KN = np.array([[-fy , -ey], [-fu, -esu], [-0.2*fu, -1.05*esu], [-0.0*fu , -1.1*esu]])  # COMPRESSION STRESS-STRAIN RELATION
    ops.uniaxialMaterial('HystereticSM', MatTag, '-posEnv', *KP.flatten(), '-negEnv', *KN.flatten(), '-pinch', 1, 1, '-damage', 0.1, 0.01, '-beta', 0),'printInput'
    #INFO LINK: https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/HystereticSM.html
    """
    MatTag = 2
    pinchX = 0.8           # Pinching factor in X direction
    pinchY = 0.5           # Pinching factor in Y direction
    damage1 = 0.0          # Damage due to ductility
    damage2 = 0.0          # Damage due to energy
    beta = 0.1             # Stiffness degradation parameter
    ops.uniaxialMaterial('Hysteretic', MatTag, fy, ey, fu, esu, 0.2*fu, 1.1*esu, -fy, -ey, -fu, -esu, -0.2*fu, -1.1*esu, pinchX, pinchY, damage1, damage2, beta)
    # INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Hysteretic_Material
    """
    MatTag = 2
    ops.uniaxialMaterial('Steel01', MatTag, fy, E, b)     
    """
    # STEEL PLATE PROPERTIES
    # Plate reinforcing steel M30 (Class 10.9)
    fyP = random.uniform(235, 245)             # Yield stress for plate
    eyP = random.uniform(0.001175, 0.001225)   # Yield strain for plate
    EP = fyP / eyP                             # Young's modulus for plate
    
    # STEEL PLATE PROPERTIES - INFO LINK: https://opensees.berkeley.edu/wiki/index.php?title=Elastic-Perfectly_Plastic_Material
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
    ops.patch('rect', 1, NUMFIBERS, 5, cover - y1, cover - z1, y1 - cover, z1 - cover)

    # Create the plate cover fibers (top, bottom, left, right)
    ops.patch('rect', 1, NUMFIBERS, 5, -y1, z1 - cover, y1, z1)
    ops.patch('rect', 1, NUMFIBERS, 5, -y1, -z1, y1, cover - z1)
    ops.patch('rect', 1, NUMFIBERS, 5, -y1, cover - z1, cover - y1, z1 - cover)
    ops.patch('rect', 1, NUMFIBERS, 5, y1 - cover, cover - z1, y1, z1 - cover)

    # Create the reinforcing fibers (left, middle, right)
    ops.layer('straight', 2, 2, As, y1 - cover, z1 - cover, y1 - cover, cover - z1)         # LAYER 01 - TOP LAYER
    ops.layer('straight', 2, 2, As, y1 - cover02, z1 - cover02, y1 - cover02, cover02 - z1) # LAYER 02 - TOP MIDDLE
    ops.layer('straight', 2, 2, As, cover02 - y1, z1 - cover02, cover02 - y1, cover - z1)   # LAYER 03 - BOTTOM MIDDLE
    ops.layer('straight', 2, 2, As, cover - y1, z1 - cover, cover - y1, cover - z1)         # LAYER 04 - BOTTOM LAYER

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
    print(f'First quantile: {q1:.6e}')
    #print(f'Mode: {MODE:.6e}')
    print(f'Median: {MEDIAN:.6e}')
    print(f'Mean: {MEAN:.6e}')
    print(f'Std: {STD:.6e}')
    print(f'Third quantile: {q3:.6e}')
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
    'Over_Strength_Factor': OVER_STRENGTH
}
# Convert to DataFrame
results_df = pd.DataFrame(DATA_TOTAL)
# Export the DataFrame to an Excel file
results_df.to_excel('MOMENT_CURVATURE_UNCERTAINTY_STEEL_BOLTED_CONNECTION_RESULTS.xlsx', index=False)
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
TITLE = f'Last Data of Moment-Curvature Analysis - Ductility Ratio: {xxc[2]/xxc[1]:.4f} - Over Strength Factor: {yyc[2]/yyc[1]:.4f}'
COLOR = 'black'
PLOT_2D(Curvatures, -Moments, xxc, yyc, _, _, XLABEL, YLABEL, TITLE, LEGEND01, LEGEND02, _, COLOR='black', Z=2) 
print(f'\t\t Ductility Ratio: {yyc[2]/yyc[1]:.4f}')
# --------------------------------------------------------------------------------------------------
X = Y_CUR
HISTO_COLOR = 'green'
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

def RANDOM_FOREST(df):
    from sklearn.model_selection import train_test_split
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.metrics import classification_report, mean_squared_error, r2_score
    import seaborn as sns
    # Target: Safety binary label (example: 1 if max displacement <= threshold)
    target_threshold = np.median(df['Yield_Moment'])  # Example threshold
    df['Safety_label'] = (df['Yield_Moment'] <= target_threshold).astype(int)

    # Step 1: Create a dataset with statistical features and response labels
    # Features: Extracted stats
    # Create a dataset with individual simulation results

    # Step 2: Split the data
    X = df.drop(columns=['Safety_label'])
    y = df['Safety_label']

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    # Step 3: Train a Random Forest Classifier
    clf = RandomForestClassifier(n_estimators=100, random_state=42)
    clf.fit(X_train, y_train)

    # Step 4: Evaluate the model
    y_pred = clf.predict(X_test)
    print(classification_report(y_test, y_pred))

    # Optional: Regression for predicting maximum displacement
    reg = RandomForestRegressor(n_estimators=100, random_state=42)
    reg.fit(X_train, y_train)
    y_pred_reg = reg.predict(X_test)

    # Evaluate regression performance
    print(f"Mean Squared Error: {mean_squared_error(y_test, y_pred_reg):.6f}")
    print(f"R2 Score: {r2_score(y_test, y_pred_reg):.6f}")

    # Step 5: Visualize feature importance
    import matplotlib.pyplot as plt

    feature_importance = clf.feature_importances_
    sns.barplot(x=feature_importance, y=X.columns)
    plt.title("Feature Importance")
    plt.show()
    
RANDOM_FOREST(results_df)    

#------------------------------------------------------------------------------------------------
####  FRAGILITY ANALYSIS
from scipy.stats import norm  
# ----------------------------
# Fragility Assessment
# ----------------------------
damage_states = {
    'Minor Damage Level': (0.5, 0.4),    # Median: 0.4, β=0.4
    'Moderate Damage Level': (2.30, 0.4),
    'Severe Damage Level': (3.60, 0.5),
    'Failure Level': (4.0, 0.5)
}

im_values = DUCT_RATIO
# --------------
# Visualization
# --------------
    
# Fragility curves
plt.figure(1, figsize=(10, 6))
# Calculate and plot fragility curves for each damage state
for damage_state, (median, beta) in damage_states.items():
    # Calculate log-normal probabilities
    ln_im = np.log(im_values)
    ln_median = np.log(median)
    probabilities = norm.cdf((ln_im - ln_median) / beta)
    plt.scatter(im_values, probabilities, marker='o', label=f'{damage_state} (η={median}, β={beta}')
    #plt.plot(im_values, probabilities, lw=2, label=f'{damage_state} (η={median}, β={beta})')
plt.xlabel('Ductility Ratio  [IM]')
plt.ylabel('Probability of Exceedance')
plt.title('Fragility Curves')
plt.legend()
plt.semilogy()
plt.ylim(0, 1.0)
plt.grid(True)
plt.tight_layout()
plt.show()    

#===========================================================

# Define damage state parameters: {Damage State: (median_IM, beta)}
damage_states = {
    'Minor Damage Level': (0.5, 0.4),    # Median: 0.4, β=0.4
    'Moderate Damage Level': (1.1, 0.4),
    'Severe Damage Level': (1.6, 0.5),
    'Failure Level': (2.0, 0.5)
}

# Generate intensity measure (IM) values
im_values = OVER_STRENGTH 
# --------------
# Visualization
# --------------
# Create plot
plt.figure(2, figsize=(10, 6))
# Calculate and plot fragility curves for each damage state
for damage_state, (median, beta) in damage_states.items():
    # Calculate log-normal probabilities
    ln_im = np.log(im_values)
    ln_median = np.log(median)
    probabilities = norm.cdf((ln_im - ln_median) / beta)
    plt.scatter(im_values, probabilities, marker='o', label=f'{damage_state} (η={median}, β={beta}')
    #plt.plot(im_values, probabilities, lw=2, label=f'{damage_state} (η={median}, β={beta})')

# Format plot
plt.xlabel('Over Stength Factor [IM]', fontsize=12)
plt.ylabel('Probability of Exceedance', fontsize=12)
plt.title('Fragility Curves', fontsize=14)
plt.legend(loc='lower right', fontsize=10)
plt.grid(True)
plt.semilogy()
plt.ylim(0, 1.0)
plt.tight_layout()
plt.show()    

#------------------------------------------------------------------------------------------------



