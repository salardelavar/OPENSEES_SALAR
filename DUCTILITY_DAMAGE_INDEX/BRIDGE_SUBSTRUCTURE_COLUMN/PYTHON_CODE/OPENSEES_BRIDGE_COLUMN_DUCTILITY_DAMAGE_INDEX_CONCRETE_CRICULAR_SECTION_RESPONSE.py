#          #########################################################################
#          #                           IN THE NAME OF ALLAH                        #
#          #              BRIDGE SUBSTRUCTURE INELASTIC RESPONSE SPECTRUM          #
#          #-----------------------------------------------------------------------#
#          #              THIS PROGRAM WRITTEN BY SALAR DELAVAR QASHQAI            #
#          #                   EMAIL: salar.d.ghashghaei@gmail.com                 #
#          #########################################################################

#%%------------------------------------------------------------------
#import the os module
import os
import math
import time
import numpy as np
import openseespy.opensees as op
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
FOLDER_NAME = 'RESPONSE'
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

#------------------------------------------------------------- 
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
    print(f'Minimum: {MINIMUM:.4f}')
    print(f'First quartile: {q1:.4f}')
    #print(f'Mode: {MODE:.4f}')
    print(f'Median: {MEDIAN:.4f}')
    print(f'Mean: {MEAN:.4f}')
    print(f'Std: {STD:.4f}')
    print(f'Third quartile: {q3:.4f}')
    print(f'Maximum: {MAXIMUM :.4f}')
    print(f'Skewness: {skew(X) :.4f}')
    print(f'kurtosis: {kurtosis(X) :.4f}')
    print(f"90% Confidence Interval: ({lower_bound:.4f}, {upper_bound:.4f})")
    print("-------------------------")

    plt.figure(figsize=(10,6))
    # Plot histogram of data
    count, bins, ignored = plt.hist(X, bins=100, color=HISTO_COLOR, density=True, align='mid')#, edgecolor="black"
    
    # Plot lognormal PDF
    x = np.linspace(min(bins), max(bins), 10000)
    pdf = (np.exp(-(x - MEAN)**2 / (2 * STD**2)) / (STD * np.sqrt(2 * np.pi)))
    plt.plot(x, pdf, linewidth=2, color='r', label="Normal PDF")
    
    # Plot vertical lines for risk measures
    plt.axvline(q1, color="black", linestyle="--", label=f"Quantile 0.25: {q1:.4f}")
    plt.axvline(MEDIAN, color="green", linestyle="--", label=f"Median: {MEDIAN:.4f}")
    plt.axvline(q3, color="black", linestyle="--", label=f"Quantile 0.75: {q3:.4f}")
    #plt.axvline(MODE, color="purple", linestyle="--", label=f"Mode: {MODE:.4f}")
    plt.axvline(MEAN, color="red", linestyle="--", label=f"Mean: {MEAN:.4f}")
    plt.axvline(MEAN-STD, color="blue", linestyle="--", label=f"Mean-Std: {MEAN-STD:.4f}")
    plt.axvline(MEAN+STD, color="blue", linestyle="--", label=f"Mean+Std: {MEAN+STD:.4f}")
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
    plt.text(q1, 1.05, f" Q1: {q1:.4f}")
    plt.text(MEDIAN, 1.1, f" Q2: {MEDIAN:.4f}")
    plt.text(q3, 1.05, f" Q3: {q3:.4f}")
    #plt.text(MODE, 1.15, f" Mode: {MODE:.4f}")
    
    #plt.text(MEAN, 0.9, f" Mean: {MEAN:.4f}")
    #plt.text(MEAN-STD, 0.9, f" Mean-Std: {MEAN-STD:.4f}")
    #plt.text(MEAN+STD, 0.9, f" Mean+Std: {MEAN+STD:.4f}")
    plt.scatter(MEAN, 1, color="red", marker="+", s=200, label=f"Mean: {MEAN:.4f}")
    plt.scatter(MEAN-STD, 1, color="green", marker="X", s=200, label=f"Mean-Std: {MEAN-STD:.4f}")
    plt.scatter(MEAN+STD, 1, color="blue", marker="*", s=200, label=f"Mean+Std:  {MEAN+STD:.4f}")
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
# -----------------------------------------------
def PLOT_TIME_HIS(x, xlabel, y1, y1label, y2, y2label, y3, y3label, y4, y4label, Z, LOG):
    ## PLOT THE DATA
    import numpy as np
    import matplotlib.pyplot as plt
    # Define colors for each dataset
    colors = ['b', 'g', 'r', 'c']

    # Create subplots based on the value of Z
    fig, axs = plt.subplots(Z, 1, figsize=(14, 14))

    # Plot each dataset with a different color
    for i, y_data in enumerate([y1, y2, y3, y4][:Z]):
        axs[i].plot(x, y_data, color=colors[i])
        axs[i].set_title(f"{[y1label, y2label, y3label, y4label][i]} - MAX ABS: {np.max(np.abs(y_data)):.6e}")
        axs[i].set_xlabel(xlabel)
        #axs[i].set_ylabel()
        axs[i].grid()
        if LOG == 1:
            axs[i].semilogy()

    # Adjust layout
    plt.tight_layout()
    plt.show()  
# -----------------------------------------------    
def MAXABS_FUN(DATA_FILE, COLUMN, I):
    import numpy as np
    # Read and process displacement data
    NameFiles = DATA_FILE
    filename = f"{NameFiles}_{I}.txt"
    D = np.loadtxt(filename)
    #print(D)
    MAXABS = np.max(np.abs([D[:, COLUMN]]))
    #print("MAX. ABS. :", MAXABS)
    return MAXABS
# -----------------------------------------------
def PLOT_2D(X, Y, Xfit, Yfit, XLABEL, YLABEL, TITLE, COLOR, Z):
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
        plt.plot(X, Y, Xfit, Yfit, 'r--', linewidth=3)
        plt.title(TITLE)
        plt.xlabel(XLABEL)
        plt.ylabel(YLABEL)
        plt.legend(['curve', 'bilinear fitted'], loc='lower right')
        plt.grid(True)
        plt.show()
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
    """
    print('+==========================+')
    print('=   Analysis curve fitted =')
    print('  Curvature    Moment')
    print('----------------------------')
    print(np.column_stack((X.T, Y.T)))
    print('+==========================+')
    """
    # EI and Ductility_Rito of Unconfined Section
    Elastic_EI = Y[1] / X[1]
    Plastic_EI = Y[2] / X[2]
    Tangent_EI = (Y[2] - Y[1]) / (X[2] - X[1])
    Ductility_Rito = X[2] / X[1]
    Over_Stregth_Factor = Y[2] / Y[1]
    """
    print('+--------------------------------------------------------------------+')
    print(f' Elastic Flextural Rigidity :             {Elastic_EI:.2f}')
    print(f' Plastic Flextural Rigidity :             {Plastic_EI:.2f}')
    print(f' Tangent Flextural Rigidity :             {Tangent_EI:.2f}')
    print(f' Section Ductility Ratio :                {Ductility_Rito:.2f}')
    print(f' Section Over Strength Factor:            {Over_Strength_Factor:.2f}')
    print('+--------------------------------------------------------------------+')
    """
    return X, Y, Elastic_EI, Plastic_EI, Tangent_EI, Ductility_Rito, Over_Stregth_Factor
# -----------------------------------------------
def PLOT3D(X, Y, Z, XLABEL, YLABEL, ZLABEL, TITLE):
    import plotly.graph_objects as go
    # Create a 3D scatter plot
    fig = go.Figure(data=[go.Scatter3d(x=X, y=Y, z=Z, mode='markers', marker=dict(size=5, color=Z))])
    fig.update_layout(scene=dict(xaxis_title=XLABEL, yaxis_title=YLABEL, zaxis_title=ZLABEL), title=TITLE)
    fig.show()
# -----------------------------------------------     
# Create a scatter plot
def PLOT_SCATTER(X, Y , XLABEL, YLABEL, TITLE, COLOR, LOG, ORDER):
    import matplotlib.pyplot as plt
    # Calculate linear regression parameters
    import numpy as np
    coefficients = np.polyfit(X, Y, ORDER)
    if ORDER == 1:
        a, b = coefficients
    if ORDER == 2:
        a, b, c = coefficients    
    if ORDER == 3:
        a, b, c, d = coefficients   
    if ORDER == 4:
        a, b, c, d, e = coefficients  
    if ORDER == 5:
        a, b, c, d, e, f = coefficients  
    if ORDER == 6:
        a, b, c, d, e, f, I = coefficients   
    if ORDER == 7:
        a, b, c, d, e, f, I, J = coefficients     
    y = [];yy = [];
    for i in range(len(X)):
        if ORDER == 1:
            y.append(a * X[i] + b)
        if ORDER == 2:
            y.append(a * X[i]**2 + b * X[i] + c)
        if ORDER == 3:
            y.append(a * X[i]**3 + b * X[i]**2 + c * X[i] + d)    
        if ORDER == 4:
            y.append(a * X[i]**4 + b * X[i]**3 + c * X[i]**2 + d * X[i] + e)  
        if ORDER == 5:
            y.append(a * X[i]**5 + b * X[i]**4 + c * X[i]**3 + d * X[i]**2 + e * X[i] + f)    
        if ORDER == 6:
            y.append(a * X[i]**6 + b * X[i]**5 + c * X[i]**4 + d * X[i]**3 + e * X[i]**2 + f * X[i] + I)     
        if ORDER == 7:
            y.append(a * X[i]**7 + b * X[i]**6 + c * X[i]**5 + d * X[i]**4 + e * X[i]**3 + f * X[i]**2 + I * X[i] + J)     
        yy.append(Y[i] - y[-1])
    y = np.array(y)    
    yy = np.array(yy) 
    # Calculate TSS
    Y_mean = np.mean(Y)
    TSS = np.sum((Y - Y_mean) ** 2)
    # Calculate RSS
    RSS = np.sum(yy ** 2)
    # Calculate R-squared
    R_squared = 1 - (RSS / TSS)
    #print(f"R-squared value: {R_squared:.4f}")
    plt.figure(figsize=(10,6))
    plt.scatter(X, Y, color=COLOR, marker='o', label='Data')
    # Add labels and title
    plt.xlabel(XLABEL)
    plt.ylabel(YLABEL)
    # Add the linear regression line
    if ORDER == 1:
        plt.plot(X, y, color='black', label=f'y = {a:.2f}x + {b:.2f} - R^2 = {R_squared:.3f}')
    if ORDER == 2:
        plt.plot(X, y, color='black', label=f'y = {a:.2f}x^2 + {b:.2f}x + {c:.2f} - R^2 = {R_squared:.3f}')
    if ORDER == 3:
        plt.plot(X, y, color='black', label=f'y = {a:.2f}x^3 + {b:.2f}x^2 + {c:.2f}x + {d:.2f} - R^2 = {R_squared:.3f}')  
    if ORDER == 4:
        plt.plot(X, y, color='black', label=f'y = {a:.2f}x^4 + {b:.2f}x^3 + {c:.2f}x^2 + {d:.2f}x + {e:.2f} - R^2 = {R_squared:.3f}') 
    if ORDER == 5:
        plt.plot(X, y, color='black', label=f'y = {a:.2f}x^5 + {b:.2f}x^4 + {c:.2f}x^3 + {d:.2f}x^2 + {e:.2f}x + {f:.2f} - R^2 = {R_squared:.3f}')  
    if ORDER == 6:
        plt.plot(X, y, color='black', label=f'y = {a:.2f}x^6 + {b:.2f}x^5 + {c:.2f}x^4 + {d:.2f}x^3 + {e:.2f}x^2 + {f:.2f}x + {I:.2f} - R^2 = {R_squared:.3f}')  
    if ORDER == 7:
        plt.plot(X, y, color='black', label=f'y = {a:.2f}x^7 + {b:.2f}x^6 + {c:.2f}x^5 + {d:.2f}x^4 + {e:.2f}x^3 + {f:.2f}x^2 + {I:.2f}x + {J:.2f} - R^2 = {R_squared:.3f}')               
    
    plt.title(TITLE)
    plt.grid(True)
    plt.legend()
    if LOG == 1:
        plt.semilogx();plt.semilogy();
    plt.show()

def plot_scatter_plotly(X, Y, XLABEL, YLABEL, TITLE, COLOR):
    import plotly.express as px
    fig = px.scatter(x=X, y=Y, color_discrete_sequence=[COLOR], labels={XLABEL: XLABEL, YLABEL: YLABEL})
    fig.update_layout(title=TITLE, xaxis_type='log', yaxis_type='log')
    fig.show() 
# ----------------------------------------------- 
def PLOT_HEATMAP(df):
    import plotly.figure_factory as ff
    # Calculate the correlation matrix
    corr_matrix = df.corr()

    # Create a correlation heatmap
    fig = ff.create_annotated_heatmap(
        z=corr_matrix.values,
        x=list(corr_matrix.columns),
        y=list(corr_matrix.index),
        annotation_text=corr_matrix.round(5).values,
        showscale=True,
        colorscale='Viridis'
    )

    # Update layout
    fig.update_layout(
        title='Correlation Heatmap',
        xaxis=dict(title='Variable'),
        yaxis=dict(title='Variable'),
        width=1200, height=1200
    )

    fig.show()
    
# -----------------------------------------------     
"""
Long short-term memory (LSTM) is a type
of recurrent neural network (RNN) aimed
at dealing with the vanishing gradient
problem present in traditional RNNs
"""
def PREDICT_LSTM(x, y, look_back, ITERATION):
    import numpy as np
    from keras.models import Sequential
    from keras.layers import LSTM, Dense
    # Prepare data for LSTM
    trainX, trainY = [], []
    for i in range(len(x) - look_back):
        trainX.append(x[i:i + look_back])
        trainY.append(y[i + look_back])

    trainX, trainY = np.array(trainX), np.array(trainY)

    # Build the LSTM model
    model = Sequential()
    model.add(LSTM(4, input_shape=(look_back, 1)))
    model.add(Dense(1))
    model.compile(loss='mean_squared_error', optimizer='adam')
    model.fit(trainX, trainY, epochs= ITERATION, batch_size=1, verbose=2)

    # Predict the next 'y' value
    next_x = np.array(x[-look_back:]).reshape(1, look_back, 1)
    predicted_y = model.predict(next_x)
    return predicted_y 

# -----------------------------------------------    
def Normal_CDF_Newton_Raphson(P_f, EPS=1e-3, tol=1e-6, max_iter=1000000):
    from scipy.stats import norm
    x = 0.0  # Initial guess (you can choose any value)
    
    for i in range(max_iter):
        xmin = x - EPS
        xmax = x + EPS
        f = norm.cdf(-x) - P_f
        fmin = norm.cdf(-xmin) - P_f
        fmax = norm.cdf(-xmax) - P_f
        df = (fmax - fmin) / (2 * EPS)
        dx = f / df
        f_prime_x = -norm.pdf(-x)
        
        if abs(dx) < tol:
            break
        
        x -= dx
    
    return x
# ----------------------------------------------- 
def MIX_HISTOGRAM(x, y, BINS, X, Y, TITLE):
    plt.figure(figsize=(8, 6))
    plt.hist(x, bins=BINS, alpha=0.5, label=X, color='blue')
    plt.hist(y, bins=BINS, alpha=0.5, label=Y, color='red')
    plt.legend(loc='upper right')
    plt.xlabel("Samples")
    plt.ylabel("Frequency")
    plt.title(TITLE)
    plt.show()
# -----------------------------------------------     
def BETA_PDF(min_x, max_x, a , b):
    return min_x + (max_x - min_x) * np.random.beta(a, b)    
# -----------------------------------------------       
#############################
# Check Different Analysis  #
#############################

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
				
  
#%%------------------------------------------------------------------
# OUTPUT DATA ADDRESS:
SALAR_DIR = f'C://OPENSEESPY_SALAR//{FOLDER_NAME}//';    
#%%------------------------------------------------------------------
### -------------------------------
###    MOMENT-CURVATURE FUNCTION
### -------------------------------

def MC_ANALYSIS(P, DR, numIncr, DSec, coverSec, numBarsSec, BD, fc, TOLERANCE, ITERATION, I):
    import openseespy.opensees as op
    def MomentCurvature(ColSecTag, axialLoad, maxK, numIncr=100):

        # Define two nodes at (0,0)
        op.node(1, 0.0, 0.0)
        op.node(2, 0.0, 0.0)

        # Fix all degrees of freedom except axial and bending
        op.fix(1, 1, 1, 1)
        op.fix(2, 0, 1, 0)

        # Define element
        #                             tag ndI ndJ  secTag
        op.element('zeroLengthSection',  1,   1,   2,  ColSecTag)
        # Create recorder
        op.recorder('Node', '-file', f"{SALAR_DIR}CUR_{I}.txt",'-time', '-node', 2, '-dof', 3, 'disp')# Curvature Time History nodes 2
        op.recorder('Node', '-file', f"{SALAR_DIR}MOM_{I}.txt",'-time', '-node', 1, '-dof', 3, 'reaction')# Base Shear Time History nodes 1
        #op.recorder('Element','-ele',1,'-file',f'{SALAR_DIR}fiberA_StressStrain.txt','section','fiber',DSec/2,0,'stressStrain')# Top fiber
        #op.recorder('Element','-ele',1,'-file',f'{SALAR_DIR}fiberB_StressStrain.txt','section','fiber',-DSec/2,0,'stressStrain')# Bottom fiber
        #op.recorder('Element','-ele',1,'-file',f'{SALAR_DIR}fiberC_StressStrain.txt','section','fiber',DSec/2 - coverSec,0,'stressStrain')# Top fiber - steel
        #op.recorder('Element','-ele',1,'-file',f'{SALAR_DIR}fiberD_StressStrain.txt','section','fiber',-DSec/2 + coverSec,0,'stressStrain')# Bottom fiber - steel
        
        # Define constant axial load
        op.timeSeries('Constant', 1)
        op.pattern('Plain', 1, 1)
        op.load(2, axialLoad, 0.0, 0.0)

        # Define analysis parameters
        op.integrator('LoadControl', 0.0)
        op.system('SparseGeneral', '-piv')
        op.test('NormUnbalance', TOLERANCE, ITERATION)
        op.numberer('Plain')
        op.constraints('Plain')
        op.algorithm('Newton')
        op.analysis('Static')

        # Do one analysis for constant axial load
        op.analyze(1)

        # Define reference moment
        op.timeSeries('Linear', 2)
        op.pattern('Plain',2, 2)
        op.load(2, 0.0, 0.0, 1.0)

        # Compute curvature increment
        dK = maxK / numIncr

        # Use displacement control at node 2 for section analysis
        op.integrator('DisplacementControl', 2,3,dK,1,dK,dK)

        # Do the section analysis
        op.analyze(numIncr)


    op.wipe()
    
    ####      Start Moment Curvature Analysis
    barAreaSec = (3.1415 * BD**2) / 4
    # Define model builder
    # --------------------
    op.model('basic','-ndm',2,'-ndf',3)
    ColSecTag = 1
    # MATERIAL parameters -------------------------------------------------------------------
    IDconcCore = 1; 				# material ID tag -- confined core concrete  - COLUMN
    IDconcCover = 2; 				# material ID tag -- unconfined cover concrete  - COLUMN
    IDreinf = 3; 				# material ID tag -- reinforcement  - COLUMN & BEAM
    # nominal concrete compressive strength
    Ec = 4700 * math.sqrt(-fc) # Concrete Elastic Modulus

    # confined concrete
    Kfc = 1.3;			# ratio of confined to unconfined concrete strength
    fc1C = Kfc*fc;		# CONFINED concrete (mander model), maximum stress
    eps1C = 2*fc1C/Ec;	# strain at maximum stress 
    fc2C = 0.2*fc1C;		# ultimate stress
    eps2C = 5*eps1C;		# strain at ultimate stress 
    # unconfined concrete
    fc1U = fc;			# UNCONFINED concrete (todeschini parabolic model), maximum stress
    eps1U = -0.0025;			# strain at maximum strength of unconfined concrete
    fc2U = 0.2*fc1U;		# ultimate stress
    eps2U = -0.012;			# strain at ultimate stress
    Lambda = 0.1;				# ratio between unloading slope at $eps2 and initial slope $Ec
    # tensile-strength properties
    ftC = -0.55*fc1C;		# tensile strength +tension
    ftU = -0.55*fc1U;		# tensile strength +tension
    Ets = ftU/0.002;		# tension softening stiffness
    # REBAR MATERIAL PROPERTIES:
    Fy = 4000			# Steel rebar yield stress
    Cy = 0.02			# Steel rebar yield strain
    Es = Fy/Cy				# modulus of steel
    Bs = 0.01				# strain-hardening ratio 
    R0 = 18.0				# control the transition from elastic to plastic branches
    cR1 = 0.925				# control the transition from elastic to plastic branches
    cR2 = 0.15				# control the transition from elastic to plastic branches

    op.uniaxialMaterial('Concrete02', IDconcCore, fc1C, eps1C, fc2C, eps2C, Lambda, ftC, Ets) # build cover concrete (confined)
    op.uniaxialMaterial('Concrete02', IDconcCover, fc1U, eps1U, fc2U, eps2U, Lambda, ftU, Ets) # build cover concrete (unconfined)
    op.uniaxialMaterial('Steel02', IDreinf, Fy, Es, Bs, R0,cR1,cR2) # build reinforcement material
    # FIBER SECTION properties -------------------------------------------------------------
    ri = 0.0;			# inner radius of the section, only for hollow sections
    ro = DSec/2;	# overall (outer) radius of the section
    rc = ro - coverSec;					# Core radius
    nfCoreR  = 80;		# number of radial divisions in the core (number of "rings")
    nfCoreT = 80;		# number of theta divisions in the core (number of "wedges")
    nfCoverR = 40;		# number of radial divisions in the cover
    nfCoverT = 80;		# number of theta divisions in the cover
    # -----------------------------------------------------------------------------------------------------
    # RETROFITTED SECTION: 
    op.section('Fiber', ColSecTag)
    # Define the core patch
    op.patch('circ', IDconcCore, nfCoreT, nfCoreR, 0, 0, ri, rc, 0, 360)
    # Define the four cover patches
    op.patch('circ', IDconcCover, nfCoverT, nfCoverR, 0, 0, rc, ro, 0, 360)
    # Define reinfocement layers 
    theta = 360.0/numBarsSec; # Determine angle increment between bars
    # Define the reinforcing layer
    op.layer('circ', IDreinf, numBarsSec, barAreaSec, 0, 0, rc, theta, 360)

    # -----------------------------------------------------------------------------------------------------
    # set yield  Curvature
    Ky = Cy / (0.5 * DSec)
    #print('Ky', Ky)

    # set ultimate Curvature
    Ku = Ky * DR
    #print('Ku', Ku)


    # Call the section analysis procedure:
    MomentCurvature(ColSecTag, P, Ku, numIncr)
    print(f'{I+1} MC Done.')
    op.wipe()
#%%------------------------------------------------------------------
### -----------------------
###    PUSHOVER FUNCTION
### -----------------------

def PUSHOVER_ANALYSIS(L, H, DSec, coverSec , numBarsSecC, BDcol, Hbeam, Bbeam, coverBeam, BDbeam, numBarsSecB, Weight, fc, ND, DMAX, I):
    import openseespy.opensees as op
    IDctrlNode = ND ## INCREMENTAL DISPLACEMENT NODE
    IDctrlDOF = 1
    # Define Analysis Properties
    MAX_ITERATIONS = 5000    # Convergence iteration for test
    MAX_TOLERANCE = 1.0e-10  # Convergence tolerance for test
    op.wipe()
    op.model('basic', '-ndm', 2, '-ndf', 3) 
    PCol = Weight
    #g =  9810 # mm/s^2
    #Mass =  PCol/g
    # nodal coordinates:
    op.node(1, 0.0, 0.0) # node#, X, Y
    op.node(2, 0.0, L)
    op.node(3, 0.5*H, 0.0) 
    op.node(4, 0.5*H, L)
    op.node(5, H, 0.0) 
    op.node(6, H, L)
    # Single point constraints -- Boundary Conditions
    op.fix(1, 1, 1, 1) # node DX DY RZ
    op.fix(3, 1, 1, 1)
    op.fix(5, 1, 1, 1)
    
    barAreaSecC = (3.1415 * BDcol**2) / 4  # area of longitudinal-reinforcement bars - COLUMN
    barAreaSecB = (3.1415 * BDbeam**2) / 4  # area of longitudinal-reinforcement bars - BEAM
    
    SecTag01 = 1			# Column Section
    SecTag02 = 2			# Beam Section
    # MATERIAL parameters -------------------------------------------------------------------
    IDconcCore = 1; 				# material ID tag -- confined core concrete - COLUMN
    IDconcCover = 2; 				# material ID tag -- unconfined cover concrete - COLUMN
    IDreinf = 3; 				# material ID tag -- reinforcement - COLUMN & BEAM
    IDconcCoreB = 4; 				# material ID tag -- confined core concrete  - BEAM
    IDconcCoverB = 5; 				# material ID tag -- unconfined cover concrete  - BEAM
    
    # nominal concrete compressive strength
    Ec = 4700 * math.sqrt(-fc) # Concrete Elastic Modulus (the term in sqr root needs to be in psi

    # confined concrete
    Kfc = 1.3;			# ratio of confined to unconfined concrete strength - COLUMN
    fc1C = Kfc*fc;		# CONFINED concrete (mander model), maximum stress - COLUMN
    KfcB = 1.1;			# ratio of confined to unconfined concrete strength - BEAM
    fc1CB = KfcB*fc;	# CONFINED concrete (mander model), maximum stress - BEAM
    
    eps1C = 2*fc1C/Ec;	# strain at maximum stress 
    fc2C = 0.2*fc1C;		# ultimate stress
    eps2C = 5*eps1C;		# strain at ultimate stress 
    # unconfined concrete
    fc1U = fc;			# UNCONFINED concrete (todeschini parabolic model), maximum stress
    eps1U = -0.0025;		# strain at maximum strength of unconfined concrete
    fc2U = 0.2*fc1U;		# ultimate stress
    eps2U = -0.012;			# strain at ultimate stress
    Lambda = 0.1;			# ratio between unloading slope at $eps2 and initial slope $Ec
    # tensile-strength properties
    ftC = -0.55*fc1C;		# tensile strength +tension
    ftU = -0.55*fc1U;		# tensile strength +tension
    Ets = ftU/0.002;		# tension softening stiffness
    # REBAR MATERIAL PROPERTIES:
    """    
    Fy = 4000			    # Steel rebar yield stress
    Cy = 0.02			    # Steel rebar yield strain
    Es = Fy/Cy				# modulus of steel
    Bs = 0.01				# strain-hardening ratio 
    R0 = 18.0				# control the transition from elastic to plastic branches
    cR1 = 0.925				# control the transition from elastic to plastic branches
    cR2 = 0.15				# control the transition from elastic to plastic branches
    op.uniaxialMaterial('Steel02', IDreinf, Fy, Es, Bs, R0,cR1,cR2) # build reinforcement material
    """
    fy = 4000         # [N/mm²] Steel Rebar Yield Strength   
    Es = 2e5          # [N/mm²] Modulus of Elasticity
    ey = fy/Es        # [mm/mm] Steel Rebar Yield Strain
    fu = 1.1818*fy    # [N/mm²] Steel Rebar Ultimate Strength
    esu = ey*75.2     # [mm/mm] Steel Rebar Ultimate Strain
    #Esh = (fu - fy)/(esu - ey)
    #Bs = Esh / Es
    pinchX = 0.8   # Pinching factor in X direction
    pinchY = 0.5   # Pinching factor in Y direction
    damage1 = 0.0  # Damage due to ductility
    damage2 = 0.0  # Damage due to energy
    beta = 0.1 # Stiffness degradation parameter
    op.uniaxialMaterial('Hysteretic', IDreinf, fy, ey, fu, esu, 0.2*fu, 1.1*esu, -fy, -ey, -fu, -esu, -0.2*fu, -1.1*esu, pinchX, pinchY, damage1, damage2, beta)
    # INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Hysteretic_Material
    # COLUMN
    op.uniaxialMaterial('Concrete02', IDconcCore, fc1C, eps1C, fc2C, eps2C, Lambda, ftC, Ets) # build cover concrete (confined)
    op.uniaxialMaterial('Concrete02', IDconcCover, fc1U, eps1U, fc2U, eps2U, Lambda, ftU, Ets) # build cover concrete (unconfined)
    # BEAM
    op.uniaxialMaterial('Concrete02', IDconcCoreB, fc1CB, eps1C, fc2C, eps2C, Lambda, ftC, Ets) # build cover concrete (confined)
    op.uniaxialMaterial('Concrete02', IDconcCoverB, fc1U, eps1U, fc2U, eps2U, Lambda, ftU, Ets) # build cover concrete (unconfined)
    # FIBER SECTION properties -------------------------------------------------------------
    # COLUMN
    ri = 0.0;			# inner radius of the section, only for hollow sections
    ro = (DSec)/2;	# overall (outer) radius of the section
    rc = ro - coverSec;					# Core radius
    nfCoreR  = 40;		# number of radial divisions in the core (number of "rings")
    nfCoreT = 40;		# number of theta divisions in the core (number of "wedges")
    nfCoverR = 40;		# number of radial divisions in the cover
    nfCoverT = 40;		# number of theta divisions in the cover
    # -----------------------------------------------------------------------------------------------------
    # CIRCULAR SECTION:
    op.section('Fiber', SecTag01)
    # Define the core patch
    op.patch('circ', IDconcCore, nfCoreT, nfCoreR, 0, 0, ri, rc, 0, 360)
    # Define the cover patches
    op.patch('circ', IDconcCover, nfCoverT, nfCoverR, 0, 0, rc, ro, 0, 360)
    # Define reinfocement layers 
    theta = 360.0/numBarsSecC;		# Determine angle increment between bars
    # Define the reinforcing layer
    op.layer('circ', IDreinf, numBarsSecC, barAreaSecC, 0, 0, rc, theta, 360)
    
    op.geomTransf('Linear', 1)
    numIntgrPts = 5
    op.element('nonlinearBeamColumn', 1, 1, 2, numIntgrPts, SecTag01, 1)
    op.element('nonlinearBeamColumn', 2, 3, 4, numIntgrPts, SecTag01, 1)
    op.element('nonlinearBeamColumn', 3, 5, 6, numIntgrPts, SecTag01, 1)
    # -----------------------------------------------------------------------------------------------------
    # BEAM
    coverY = Hbeam/2.0	# The distance from the section z-axis to the edge of the cover concrete -- outer edge of cover concrete
    coverZ = Bbeam/2.0	# The distance from the section y-axis to the edge of the cover concrete -- outer edge of cover concrete
    
    coreY = coverY - coverCol
    coreZ = coverZ - coverCol
    coreY02 = coreY - 120 # Middle Rebar Distance
    coreZ02 = coreZ      # Middle Rebar Distance
    
    nfCoreY = 60;			# number of fibers for concrete in y-direction -- core concrete
    nfCoreZ = 20;			# number of fibers for concrete in z-direction
    nfCoverY = 60;			# number of fibers for concrete in y-direction -- cover concrete
    nfCoverZ = 20;			# number of fibers for concrete in z-direction

    op.section('Fiber', SecTag02)
    # Define the core patch
    op.patch('quad', IDconcCoreB, nfCoreZ, nfCoreY, -coreY,coreZ, -coreY,-coreZ, coreY,-coreZ, coreY, coreZ) # Define the concrete patch
    # Define the four cover patches
    op.patch('quad', IDconcCoverB, nfCoverZ, nfCoverY, -coverY,coverZ, -coreY, coreZ, coreY, coreZ, coverY,coverZ) # Define the concrete patch
    op.patch('quad', IDconcCoverB, nfCoverZ, nfCoverY, -coreY, -coreZ, -coverY, -coverZ, coverY, -coverZ, coreY, -coreZ) # Define the concrete patch
    op.patch('quad', IDconcCoverB, nfCoverZ, nfCoverY, -coverY, coverZ, -coverY, -coverZ, -coreY, -coreZ, -coreY, coreZ) # Define the concrete patch
    op.patch('quad', IDconcCoverB, nfCoverZ, nfCoverY, coreY, coreZ, coreY, -coreZ, coverY,-coverZ, coverY,coverZ) # Define the concrete patch
    # Define reinfocement layers    
    op.layer('straight', IDreinf, numBarsSecB, barAreaSecB, coreY, coreZ, coreY, -coreZ)# top layer reinforcement
    op.layer('straight', IDreinf, 2, barAreaSecB, coreY02, coreZ02, coreY02, -coreZ02)# middle top layer reinforcement
    op.layer('straight', IDreinf, 2, barAreaSecB, -coreY02, coreZ02, -coreY02, -coreZ02)# middle bottom layer reinforcement
    op.layer('straight', IDreinf, numBarsSecB, barAreaSecB, -coreY, coreZ, -coreY, -coreZ)# bottom layer reinfocement
    
    op.geomTransf('Linear', 2)
    numIntgrPts = 5
    op.element('nonlinearBeamColumn', 4, 2, 4, numIntgrPts, SecTag02, 2)
    op.element('nonlinearBeamColumn', 5, 4, 6, numIntgrPts, SecTag02, 2)
    

    #import InelasticFiberSection
    op.recorder('Node', '-file', f"{SALAR_DIR}DTH_PUSH_{I}.txt",'-time', '-node', ND, '-dof', 1,2,3, 'disp')# Displacement Time History Node 2
    op.recorder('Node', '-file', f"{SALAR_DIR}BTH_PUSH_01_{I}.txt",'-time', '-node', 1, '-dof', 1,2,3, 'reaction')# Base Shear Time History Node 1
    op.recorder('Node', '-file', f"{SALAR_DIR}BTH_PUSH_03_{I}.txt",'-time', '-node', 3, '-dof', 1,2,3, 'reaction')# Base Shear Time History Node 3
    op.recorder('Node', '-file', f"{SALAR_DIR}BTH_PUSH_05_{I}.txt",'-time', '-node', 5, '-dof', 1,2,3, 'reaction')# Base Shear Time History Node 5
    
    #defining gravity loads
    op.timeSeries('Linear', 1)
    op.pattern('Plain', 1, 1)
    op.load(2, 0.0, -PCol, 0.0)
    op.load(4, 0.0, -PCol, 0.0)
    op.load(6, 0.0, -PCol, 0.0)

    Tol = 1e-8 # convergence tolerance for test
    Iter = 1000# convergence iteration for test
    NstepGravity = 10
    DGravity = 1 / NstepGravity
    op.integrator('LoadControl', DGravity) # determine the next time step for an analysis
    op.numberer('Plain') # renumber dof's to minimize band-width (optimization), if you want to
    op.system('BandGeneral') # how to store and solve the system of equations in the analysis
    op.constraints('Plain') # how it handles boundary conditions
    op.test('NormDispIncr', Tol, Iter) # determine if convergence has been achieved at the end of an iteration step
    op.algorithm('Newton') # use Newton's solution algorithm: updates tangent stiffness at every iteration
    op.analysis('Static') # define type of analysis static or transient
    op.analyze(NstepGravity) # apply gravity

    op.loadConst('-time', 0.0) #maintain constant gravity loads and reset time to zero
    #print('Model Built')
    
    Dincr = 0.001 * DMAX
    Hload = 1#Weight
    maxNumIter = 1000
    tol = 1e-8

    op.timeSeries('Linear', 2)
    op.pattern('Plain', 200, 2)
    op.load(ND, Hload, 0.0, 0.0)

    op.wipeAnalysis()
    op.constraints('Plain')
    op.numberer('Plain')
    op.system('BandGeneral')
    op.test('EnergyIncr', Tol, maxNumIter)
    op.algorithm('Newton')

    op.integrator('DisplacementControl', IDctrlNode, IDctrlDOF, Dincr)
    op.analysis('Static')


    Nsteps =  int(DMAX/ Dincr)

    OK = op.analyze(Nsteps)
    ANALYSIS(OK, Nsteps, MAX_TOLERANCE, MAX_ITERATIONS)
    print(f'{I+1} Pushover Done.')
    op.wipe()
    
#%%------------------------------------------------------------------
### ----------------------
###    DYNAMIC FUNCTION
### ----------------------

def DYNAMIC_ANALYSIS(L, H, DSec, coverSec , numBarsSecC, BDcol, Hbeam, Bbeam, coverBeam, BDbeam, numBarsSecB, Mass, Weight, fc, I):
    import openseespy.opensees as op
    # Define Analysis Properties
    MAX_ITERATIONS = 5000    # Convergence iteration for test
    MAX_TOLERANCE = 1.0e-10  # Convergence tolerance for test
    op.wipe()
    op.model('basic', '-ndm', 2, '-ndf', 3) 
    #PCol = Weight
    #g =  9810 # mm/s^2
    g =  1 # BEACUSE WE IMPORTED MASS
    #Mass =  PCol/g
    # nodal coordinates:
    op.node(1, 0.0, 0.0) # node#, X, Y
    op.node(2, 0.0, L)
    op.node(3, 0.5*H, 0.0) 
    op.node(4, 0.5*H, L)
    op.node(5, H, 0.0) 
    op.node(6, H, L)
    # Single point constraints -- Boundary Conditions
    op.fix(1, 1, 1, 1) # node DX DY RZ
    op.fix(3, 1, 1, 1)
    op.fix(5, 1, 1, 1)
    # node#, Mx My Mz, Mass=Weight/g, neglect rotational inertia at nodes
    op.mass(2, Mass, Mass, 0.0)
    op.mass(4, Mass, Mass, 0.0)
    op.mass(6, Mass, Mass, 0.0)
    
    barAreaSecC = (3.1415 * BDcol**2) / 4  # area of longitudinal-reinforcement bars - COLUMN
    barAreaSecB = (3.1415 * BDbeam**2) / 4  # area of longitudinal-reinforcement bars - BEAM
    
    SecTag01 = 1			# Column Section
    SecTag02 = 2			# Beam Section
    # MATERIAL parameters -------------------------------------------------------------------
    IDconcCore = 1; 				# material ID tag -- confined core concrete - COLUMN
    IDconcCover = 2; 				# material ID tag -- unconfined cover concrete - COLUMN
    IDreinf = 3; 				   # material ID tag -- reinforcement - COLUMN & BEAM
    IDconcCoreB = 4; 				# material ID tag -- confined core concrete  - BEAM
    IDconcCoverB = 5; 				# material ID tag -- unconfined cover concrete  - BEAM
    
    # nominal concrete compressive strength
    Ec = 4700 * math.sqrt(-fc) # Concrete Elastic Modulus

    # confined concrete
    Kfc = 1.3;			# ratio of confined to unconfined concrete strength - COLUMN
    fc1C = Kfc*fc;		# CONFINED concrete (mander model), maximum stress - COLUMN
    KfcB = 1.1;			# ratio of confined to unconfined concrete strength - BEAM
    fc1CB = KfcB*fc;	# CONFINED concrete (mander model), maximum stress - BEAM
    
    eps1C = 2*fc1C/Ec;	# strain at maximum stress 
    fc2C = 0.2*fc1C;		# ultimate stress
    eps2C = 5*eps1C;		# strain at ultimate stress 
    # unconfined concrete
    fc1U = fc;			    # UNCONFINED concrete (todeschini parabolic model), maximum stress
    eps1U = -0.0025;	    # strain at maximum strength of unconfined concrete
    fc2U = 0.2*fc1U;		# ultimate stress
    eps2U = -0.012;			# strain at ultimate stress
    Lambda = 0.1;				# ratio between unloading slope at $eps2 and initial slope $Ec
    # tensile-strength properties
    ftC = -0.55*fc1C;		# tensile strength +tension
    ftU = -0.55*fc1U;		# tensile strength +tension
    Ets = ftU/0.002;		# tension softening stiffness
    # REBAR MATERIAL PROPERTIES:
    """    
    Fy = 4000			    # Steel rebar yield stress
    Cy = 0.02			    # Steel rebar yield strain
    Es = Fy/Cy				# modulus of steel
    Bs = 0.01				# strain-hardening ratio 
    R0 = 18.0				# control the transition from elastic to plastic branches
    cR1 = 0.925				# control the transition from elastic to plastic branches
    cR2 = 0.15				# control the transition from elastic to plastic branches
    op.uniaxialMaterial('Steel02', IDreinf, Fy, Es, Bs, R0,cR1,cR2) # build reinforcement material
    """
    fy = 400          # [N/mm²] Steel Rebar Yield Strength   
    Es = 2e5          # [N/mm²] Modulus of Elasticity
    ey = fy/Es        # [mm/mm] Steel Rebar Yield Strain
    fu = 1.1818*fy    # [N/mm²] Steel Rebar Ultimate Strength
    esu = ey*75.2     # [mm/mm] Steel Rebar Ultimate Strain
    #Esh = (fu - fy)/(esu - ey)
    #Bs = Esh / Es
    pinchX = 0.8   # Pinching factor in X direction
    pinchY = 0.5   # Pinching factor in Y direction
    damage1 = 0.0  # Damage due to ductility
    damage2 = 0.0  # Damage due to energy
    beta = 0.1 # Stiffness degradation parameter
    op.uniaxialMaterial('Hysteretic', IDreinf, fy, ey, fu, esu, 0.2*fu, 1.1*esu, -fy, -ey, -fu, -esu, -0.2*fu, -1.1*esu, pinchX, pinchY, damage1, damage2, beta)
    # INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Hysteretic_Material
    # COLUMN
    op.uniaxialMaterial('Concrete02', IDconcCore, fc1C, eps1C, fc2C, eps2C, Lambda, ftC, Ets) # build cover concrete (confined)
    op.uniaxialMaterial('Concrete02', IDconcCover, fc1U, eps1U, fc2U, eps2U, Lambda, ftU, Ets) # build cover concrete (unconfined)
    # BEAM
    op.uniaxialMaterial('Concrete02', IDconcCoreB, fc1CB, eps1C, fc2C, eps2C, Lambda, ftC, Ets) # build cover concrete (confined)
    op.uniaxialMaterial('Concrete02', IDconcCoverB, fc1U, eps1U, fc2U, eps2U, Lambda, ftU, Ets) # build cover concrete (unconfined)
    # FIBER SECTION properties -------------------------------------------------------------
    # COLUMN
    ri = 0.0;			# inner radius of the section, only for hollow sections
    ro = (DSec)/2;	# overall (outer) radius of the section
    rc = ro - coverSec;					# Core radius
    nfCoreR  = 40;		# number of radial divisions in the core (number of "rings")
    nfCoreT = 40;		# number of theta divisions in the core (number of "wedges")
    nfCoverR = 40;		# number of radial divisions in the cover
    nfCoverT = 40;		# number of theta divisions in the cover
    # -----------------------------------------------------------------------------------------------------
    # CIRCULAR SECTION:
    op.section('Fiber', SecTag01)
    # Define the core patch
    op.patch('circ', IDconcCore, nfCoreT, nfCoreR, 0, 0, ri, rc, 0, 360)
    # Define the cover patches
    op.patch('circ', IDconcCover, nfCoverT, nfCoverR, 0, 0, rc, ro, 0, 360)
    # Define reinfocement layers 
    theta = 360.0/numBarsSecC;		# Determine angle increment between bars
    # Define the reinforcing layer
    op.layer('circ', IDreinf, numBarsSecC, barAreaSecC, 0, 0, rc, theta, 360)
    
    op.geomTransf('Linear', 1)
    numIntgrPts = 5
    op.element('nonlinearBeamColumn', 1, 1, 2, numIntgrPts, SecTag01, 1)
    op.element('nonlinearBeamColumn', 2, 3, 4, numIntgrPts, SecTag01, 1)
    op.element('nonlinearBeamColumn', 3, 5, 6, numIntgrPts, SecTag01, 1)
    # -----------------------------------------------------------------------------------------------------
    # BEAM
    coverY = Hbeam/2.0	# The distance from the section z-axis to the edge of the cover concrete -- outer edge of cover concrete
    coverZ = Bbeam/2.0	# The distance from the section y-axis to the edge of the cover concrete -- outer edge of cover concrete
    
    coreY = coverY - coverCol
    coreZ = coverZ - coverCol
    coreY02 = coreY - 120   # Middle Rebar Distance
    coreZ02 = coreZ         # Middle Rebar Distance
    
    nfCoreY = 60;			# number of fibers for concrete in y-direction -- core concrete
    nfCoreZ = 20;			# number of fibers for concrete in z-direction
    nfCoverY = 60;			# number of fibers for concrete in y-direction -- cover concrete
    nfCoverZ = 20;			# number of fibers for concrete in z-direction
    
    op.section('Fiber', SecTag02)
    # Define the core patch
    op.patch('quad', IDconcCoreB, nfCoreZ, nfCoreY, -coreY,coreZ, -coreY,-coreZ, coreY,-coreZ, coreY, coreZ) # Define the concrete patch
    # Define the four cover patches
    op.patch('quad', IDconcCoverB, nfCoverZ, nfCoverY, -coverY,coverZ, -coreY, coreZ, coreY, coreZ, coverY,coverZ) # Define the concrete patch
    op.patch('quad', IDconcCoverB, nfCoverZ, nfCoverY, -coreY, -coreZ, -coverY, -coverZ, coverY, -coverZ, coreY, -coreZ) # Define the concrete patch
    op.patch('quad', IDconcCoverB, nfCoverZ, nfCoverY, -coverY, coverZ, -coverY, -coverZ, -coreY, -coreZ, -coreY, coreZ) # Define the concrete patch
    op.patch('quad', IDconcCoverB, nfCoverZ, nfCoverY, coreY, coreZ, coreY, -coreZ, coverY,-coverZ, coverY,coverZ) # Define the concrete patch
    # Define reinfocement layers    
    op.layer('straight', IDreinf, numBarsSecB, barAreaSecB, coreY, coreZ, coreY, -coreZ)# top layer reinforcement
    op.layer('straight', IDreinf, 2, barAreaSecB, coreY02, coreZ02, coreY02, -coreZ02)# middle top layer reinforcement
    op.layer('straight', IDreinf, 2, barAreaSecB, -coreY02, coreZ02, -coreY02, -coreZ02)# middle bottom layer reinforcement
    op.layer('straight', IDreinf, numBarsSecB, barAreaSecB, -coreY, coreZ, -coreY, -coreZ)# bottom layer reinfocement
    
    op.geomTransf('Linear', 2)
    numIntgrPts = 5
    op.element('nonlinearBeamColumn', 4, 2, 4, numIntgrPts, SecTag02, 2)
    op.element('nonlinearBeamColumn', 5, 4, 6, numIntgrPts, SecTag02, 2)
    

    #import InelasticFiberSection
    op.recorder('EnvelopeNode','-file', f"{SALAR_DIR}MD_{I}.txt" ,'-time','-node',2,'-dof',1,'disp');# max. displacements of free nodes 2
    op.recorder('EnvelopeNode','-file',f"{SALAR_DIR}MV_{I}.txt" ,'-time','-node',2,'-dof',1,'vel');# max. vel of free nodes 2
    op.recorder('EnvelopeNode','-file', f"{SALAR_DIR}MA_{I}.txt" ,'-time','-node',2,'-dof',1,'accel');# max. accel of free nodes 2	
    op.recorder('Element', '-file', f"{SALAR_DIR}DEF_{I}.txt",'-time', '-ele', 1, 'section', 1, 'deformations')# Curvature Time History nodes 1
    op.recorder('Element','-ele',5,'-file',f'{SALAR_DIR}fiberCon_StressStrain_{I}.txt','section', 5,'fiber',Hbeam/2 - 5,0,'stressStrain')# concrete fiber
    op.recorder('Element','-ele',5,'-file',f'{SALAR_DIR}fiberReb_StressStrain_{I}.txt','section', 5,'fiber',coreY,coreZ,'stressStrain')# steel rebar fiber 
    op.recorder('Node', '-file', f"{SALAR_DIR}DTH_DYN_{I}.txt",'-time', '-node', 2, '-dof', 1,2,3, 'disp')# Displacement Time History Node 2
    op.recorder('Node', '-file', f"{SALAR_DIR}VTH_DYN_{I}.txt",'-time', '-node', 2, '-dof', 1,2,3, 'vel')# Velocity Time History Node 2
    op.recorder('Node', '-file', f"{SALAR_DIR}ATH_DYN_{I}.txt",'-time', '-node', 2, '-dof', 1,2,3, 'accel')# Acceleration Time History Node 2
    op.recorder('Node', '-file', f"{SALAR_DIR}BTH_DYN_01_{I}.txt",'-time', '-node', 1, '-dof', 1,2,3, 'reaction')# Base Shear Time History Node 1
    op.recorder('Node', '-file', f"{SALAR_DIR}BTH_DYN_03_{I}.txt",'-time', '-node', 3, '-dof', 1,2,3, 'reaction')# Base Shear Time History Node 3
    op.recorder('Node', '-file', f"{SALAR_DIR}BTH_DYN_05_{I}.txt",'-time', '-node', 5, '-dof', 1,2,3, 'reaction')# Base Shear Time History Node 5
    #defining gravity loads
    op.timeSeries('Linear', 1)
    op.pattern('Plain', 1, 1)
    op.load(2, 0.0, -Weight, 0.0)
    op.load(4, 0.0, -Weight, 0.0)
    op.load(6, 0.0, -Weight, 0.0)

    
    Tol = 1e-8 # convergence tolerance for test
    NstepGravity = 10
    DGravity = 1/NstepGravity
    op.integrator('LoadControl', DGravity) # determine the next time step for an analysis
    op.numberer('Plain') # renumber dof's to minimize band-width (optimization), if you want to
    op.system('BandGeneral') # how to store and solve the system of equations in the analysis
    op.constraints('Plain') # how it handles boundary conditions
    op.test('NormDispIncr', Tol, 6) # determine if convergence has been achieved at the end of an iteration step
    op.algorithm('Newton') # use Newton's solution algorithm: updates tangent stiffness at every iteration
    op.analysis('Static') # define type of analysis static or transient
    op.analyze(NstepGravity) # apply gravity

    op.loadConst('-time', 0.0) #maintain constant gravity loads and reset time to zero

    #%% DEFINE PARAMEETRS FOR NONLINEAR DYNAMIC ANALYSIS
    GMfact = 9810    # [mm/s²]standard acceleration of gravity or standard acceleration
    SSF_X = 1.0      # Seismic Acceleration Scale Factor in X Direction
    SSF_Y = 1.0      # Seismic Acceleration Scale Factor in Y Direction
    iv0_X = 0.0005   # [mm/s] Initial velocity applied to the node  in X Direction
    iv0_Y = 0.0005   # [mm/s] Initial velocity applied to the node  in Y Direction
    st_iv0 = 0.0     # [s] Initial velocity applied starting time
    SEI = 'X'        # Seismic Direction
    DR = 0.05        # Intial Guess for Damping ratio
    duration = 15.0  # [s] Total simulation duration
    dt = 0.01        # [s] Time step
                    
    # Dynamic analysis
    op.constraints('Plain')
    op.numberer('Plain')
    op.system('BandGeneral')
    op.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS)
    #ops.integrator('Newmark', 0.5, 0.25) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/newmark.html
    alpha = 1;gamma=1.5-alpha; beta=((2-alpha)**2)/4;
    op.integrator('HHT', alpha, gamma, beta) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/hht.html
    op.algorithm('Newton') # INFO LINK: https://openseespydoc.readthedocs.io/en/stable/src/algorithm.html
    op.analysis('Transient')
        
    # Calculate Rayleigh damping factors
    Lambda01 = op.eigen('-fullGenLapack', 2)  # eigenvalue mode 2
    #Lambda01 = ops.eigen('-genBandArpack', 2) # eigenvalue mode 2
    Omega01 = np.power(max(Lambda01), 0.5)
    Omega02 = np.power(min(Lambda01), 0.5)
    a0 = (2 * Omega01 * Omega02 * DR) / (Omega01 + Omega02) # c = a0 * m : Mass-proportional damping
    a1 = (DR * 2) / (Omega01 + Omega02)   # c = a1 * k : Stiffness-proportional damping
    # Apply Rayleigh damping
    op.rayleigh(a0, a1, 0, 0)   # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/reyleigh.html
    #op.rayleigh(0, 0, 2 * DR * Omega01, 0) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/reyleigh.html
    PERIOD_01 = (np.pi * 2) / Omega01 # Structure First Period
    PERIOD_02 = (np.pi * 2) / Omega02 # Structure Second Period
    print('Structure First Period:  ', PERIOD_01)
    print('Structure Second Period: ', PERIOD_02) 
        
    # Define time series for input motion (Acceleration time history)
    if SEI == 'X':
        SEISMIC_TAG_01 = 100
        op.timeSeries('Path', SEISMIC_TAG_01, '-dt', dt, '-filePath', 'Ground_Acceleration_X.txt', '-factor', GMfact, '-startTime', st_iv0) # SEISMIC-X
        # Define load patterns
        # pattern UniformExcitation $patternTag $dof -accel $tsTag <-vel0 $vel0> <-fact $cFact>
        op.pattern('UniformExcitation', SEISMIC_TAG_01, 1, '-accel', SEISMIC_TAG_01, '-vel0', iv0_X, '-fact', SSF_X) # SEISMIC-X
    if SEI == 'Y':
        SEISMIC_TAG_02 = 200
        op.timeSeries('Path', SEISMIC_TAG_02, '-dt', dt, '-filePath', 'Ground_Acceleration_Y.txt', '-factor', GMfact) # SEISMIC-Z
        op.pattern('UniformExcitation', SEISMIC_TAG_02, 2, '-accel', SEISMIC_TAG_02, '-vel0', iv0_Y, '-fact', SSF_Y) 
    if SEI == 'XY':
        SEISMIC_TAG_01 = 100
        op.timeSeries('Path', SEISMIC_TAG_01, '-dt', dt, '-filePath', 'Ground_Acceleration_X.txt', '-factor', GMfact, '-startTime', st_iv0) # SEISMIC-X
        # Define load patterns
        # pattern UniformExcitation $patternTag $dof -accel $tsTag <-vel0 $vel0> <-fact $cFact>
        op.pattern('UniformExcitation', SEISMIC_TAG_01, 1, '-accel', SEISMIC_TAG_01, '-vel0', iv0_X, '-fact', SSF_X) # SEISMIC-X 
        SEISMIC_TAG_02 = 200
        op.timeSeries('Path', SEISMIC_TAG_02, '-dt', dt, '-filePath', 'Ground_Acceleration_Y.txt', '-factor', GMfact) # SEISMIC-Z
        op.pattern('UniformExcitation', SEISMIC_TAG_02, 2, '-accel', SEISMIC_TAG_02, '-vel0', iv0_Y, '-fact', SSF_Y)  # SEISMIC-Z
    #print('Seismic Defined Done.')
    
    # Data storage
    FORCE_S, FORCE_A, MOMENT = [], [], []
    DISP_X, DISP_Y, ROT = [], [], []
    KA, KS, KI = [], [], []
    time = []
    displacement = []
    velocity_X, velocity_Y = [], []
    acceleration_X, acceleration_Y = [], []
            
    stable = 0
    current_time = 0.0
        
    while stable == 0 and current_time < duration:   
        stable = op.analyze(1, dt)
        ANALYSIS(stable, 1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS
        current_time = op.getTime()
        time.append(current_time)
        # Record results
        op.reactions()
        S = op.nodeReaction(1, 1) + op.nodeReaction(3, 1) + op.nodeReaction(5, 1)# SHEAR BASE REACTION
        A = op.nodeReaction(1, 2) + op.nodeReaction(3, 2) + op.nodeReaction(5, 2)# AXIAL BASE REACTION
        M = op.nodeReaction(1, 3) + op.nodeReaction(3, 3) + op.nodeReaction(5, 3)# MOMENT BASE REACTION
        #print(rot, M)
        disp_X = op.nodeDisp(2, 1) # LATERAL DISPLACEMENT IN X FOR NODE 2
        disp_Y = op.nodeDisp(2, 2) # LATERAL DISPLACEMENT IN Y FOR NODE 2
        rot = op.nodeDisp(2, 3)    # ROTATION IN Z FOR NODE 2
        velocity_X.append(op.nodeVel(2, 1))       # LATERAL VELOCITY IN X FOR NODE 2
        acceleration_X.append(op.nodeAccel(2, 1)) # LATERAL ACCELERATION IN X FOR NODE 2
        velocity_Y.append(op.nodeVel(2, 2))       # LATERAL VELOCITY IN Y FOR NODE 2
        acceleration_Y.append(op.nodeAccel(2, 2)) # LATERAL ACCELERATION IN Y FOR NODE 2
        FORCE_S.append(S)
        FORCE_A.append(A)
        MOMENT.append(M)
        DISP_X.append(disp_X)
        DISP_Y.append(disp_Y)
        ROT.append(rot)
        KS.append(np.abs(S/disp_X)) # LATERAL STIFFNESS IN X
        KA.append(np.abs(A/disp_Y)) # LATERAL STIFFNESS IN Y
        KI.append(np.abs(M/rot))    # ROTATIONAL STIFFNESS IN Z
        #print(current_time, disp_X, S)

    #OK = op.analyze(Nsteps, DtAnalysis)
    #ANALYSIS(OK, Nsteps, MAX_TOLERANCE, MAX_ITERATIONS)
    
    # Calculating Damping Ratio Using Logarithmic Decrement Analysis 
    displacement = np.array(DISP_X)
    peaks = np.array([displacement[i] for i in range(1, len(displacement)-1) if displacement[i] > displacement[i-1] and displacement[i] > displacement[i+1]])
    # Natural logarithm
    delta = np.log(peaks[:-1] / peaks[1:]) 
    
    print(f'{I+1} Dynamic Done.')
    op.wipe()

    print(PERIOD_01)
    return PERIOD_01
    
#%%------------------------------------------------------------------
### -----------------------------------------
###   Bridge Substructure Response Spectrum
### -----------------------------------------

# define section geometry
L = 2098.0 + 207  # [mm] Column length
H = 4138 - 2* 414 - 310 # [mm] Beam length

DR = 1.5 # set ductility ratio for moment curvature
numIncr = 100# Number of analysis increments
TOLERANCE = 1e-6 # Converganve tolerance for moment curvature
MAX_ITERATION = 1000000 # Maximum Tolerance for moment curvature

DSec = 310 # [mm] Column Diameter
coverCol = 25.0   # [mm] Column cover to reinforcing steel NA.
numBarsCol = 15  # number of longitudinal-reinforcement bars in column. (symmetric top & bot)
BDcol = 10  # [mm] Rebar Diamater

Hbeam = 414 # [mm] Column Depth
Bbeam = 414 # [mm] Column Width
coverBeam = 25.0   # [mm] Beam cover to reinforcing steel NA.
numBarsSecB = 10  # number of longitudinal-reinforcement bars in column. (symmetric top & bot)
BDbeam = 10  # [mm] Beam Rebar Diamater
# -------------------------

fc = -35.0 # [N/mm^2] Concrete Compressive Strength (+Tension, -Compression)
g =  9810 # mm/s^2

DMAX = 500 # [mm] Max. Pushover Incremental Displacement

Eef = 4700 * (-fc)** 0.5
Ief = (np.pi * DSec**4) / 64
Kef = 3 * (12*Eef*Ief) / L**3
#print("Effective Structural Lateral Stiffness: ", Kef)


DATA_FILE05 =f'C:\\OPENSEESPY_SALAR\\{FOLDER_NAME}\\T.txt' #STRUCTURE PERIOD

t = time.localtime()
current_time = time.strftime("%H:%M:%S", t)
print(f"Current time (HH:MM:SS): {current_time}\n\n")

NUM_ITERATION = 100
Tmax = 1 # Maximum Period 
Tdmax = Tmax / NUM_ITERATION
with open(DATA_FILE05, "w") as file:
    for i in range(NUM_ITERATION):
        T = (i+1) * Tdmax
        Massef = Kef * (T / (2*np.pi))**2
        PCol = Massef / 3 # Mass of each column
        PColw = PCol * g # Weight of each column
        print(f'{i+1} MASS: {PCol:.3f} WEIGHT: {PColw:.3f}')
        MC_ANALYSIS(-PColw, DR, numIncr, DSec, coverCol, numBarsCol, BDcol, fc, TOLERANCE, MAX_ITERATION, i)
        PUSHOVER_ANALYSIS(L, H, DSec, coverCol , numBarsCol, BDcol, Hbeam, Bbeam, coverBeam, BDbeam, numBarsSecB, PColw, fc, 2, DMAX, i)
        a = DYNAMIC_ANALYSIS(L, H, DSec, coverCol , numBarsCol, BDcol, Hbeam, Bbeam, coverBeam, BDbeam, numBarsSecB, PCol, PColw, fc, i)
        file.write(f"{i+1} {a:.4f} {Massef:.4f}\n")
        print(f'Real Period: {a:.4f} - Period: {T:.4f}')
    
t = time.localtime()
current_time = time.strftime("%H:%M:%S", t)
print(f"Current time (HH:MM:SS): {current_time}\n\n")
#%%------------------------------------------------------------------
## FILE ADDRESS
DATA_FILE01 =f'C:\\OPENSEESPY_SALAR\\{FOLDER_NAME}\\MD'  # DISPLACEMENT TIME HISTORY
DATA_FILE02 =f'C:\\OPENSEESPY_SALAR\\{FOLDER_NAME}\\MV'  # VELOCITY TIME HISTORY
DATA_FILE03 =f'C:\\OPENSEESPY_SALAR\\{FOLDER_NAME}\\MA'  # ACCELERATION TIME HISTORY
DATA_FILE04 =f'C:\\OPENSEESPY_SALAR\\{FOLDER_NAME}\\DEF'  # DEFORMATION TIME HISTORY
#%%------------------------------------------------------------------
### LOAD OUTPUT DATAS

pgd = [] # Peak Ground Displacement 
pgv = [] # Peak Ground Velocity
pga = [] # Peak Ground Acceleration
pgb = [] # Base shear during dynamic
period = []
diS = [] # Structure Ductility Damage Index
diC = [] # Section Ductility Damage Index
diAS = [] # Steel Rebar Axial Ductility Damage Index
diAC = [] # Concrete Axial Ductility Damage Index
momr = [] # Resistance Moment Capacity from Moment Curvature Analysis
moml = [] # Applied Moment from Dynamic Analysis
bDbP = [] # Ratio Base-shear Dynamic to Base-shear Pushover
dDdP = [] # Ratio Displacement Dynamic to Displacement Pushover

ys_strain = 0.02 # Yield strain - STEEL REBAR
us_strain = 0.1  # Ultimate strain - STEEL REBAR

yc_strain = 0.0002 # Yield strain - CONFINED CONCRETE
uc_strain = 0.008  # Ultimate strain -  CONFINED CONCRETE

Tmax = 1 # Maximum Period 
NUM_ITERATION = 100
Tdmax = Tmax / NUM_ITERATION

for i in range(NUM_ITERATION):
    T = (i+1) * Tdmax
    period.append(T)
    pgd.append(MAXABS_FUN(DATA_FILE01, 1, i))
    pgv.append(MAXABS_FUN(DATA_FILE02, 1, i))
    pga.append(MAXABS_FUN(DATA_FILE03, 1, i))
    base01 = OUTPUT_SECOND_COLUMN(FOLDER_NAME,'BTH_DYN_01', 1, i, 2) # Reading base shear from Text file - DYNAMIC - NODE 1
    base02 = OUTPUT_SECOND_COLUMN(FOLDER_NAME,'BTH_DYN_03', 1, i, 2) # Reading base shear from Text file - DYNAMIC - NODE 3
    base03 = OUTPUT_SECOND_COLUMN(FOLDER_NAME,'BTH_DYN_05', 1, i, 2) # Reading base shear from Text file - DYNAMIC - NODE 5
    baseD = max(abs(base01 + base02 + base03))
    pgb.append(baseD)
    # STRUCTURE DUCTILITY DAMAGE INDEX
    dispP = OUTPUT_SECOND_COLUMN(FOLDER_NAME,'DTH_PUSH', 1, i, 2) # Reading Disp from Text file - PUSHOVER
    base01 = OUTPUT_SECOND_COLUMN(FOLDER_NAME,'BTH_PUSH_01', 1, i, 2) # Reading base shear from Text file - PUSHOVER - NODE 1
    base02 = OUTPUT_SECOND_COLUMN(FOLDER_NAME,'BTH_PUSH_03', 1, i, 2) # Reading base shear from Text file - PUSHOVER - NODE 3
    base03 = OUTPUT_SECOND_COLUMN(FOLDER_NAME,'BTH_PUSH_05', 1, i, 2) # Reading base shear from Text file - PUSHOVER - NODE 5
    baseP = abs(base01 + base02 + base03)
    bDbP.append(pgb[-1] / max(baseP))
    dDdP.append(pgd[-1] / max(dispP))
    xx, yy, _, _, _, _, _ = BILNEAR_CURVE(dispP, baseP, 10)
    demand_disp = MAXABS_FUN(DATA_FILE01, 1, i)# DISPLACEMENT DYNAMIC ANALYSIS
    DIs = (demand_disp - xx[1]) / (xx[2] - xx[1])
    diS.append(DIs)
    # SECTION DUCTILITY DAMAGE INDEX
    CUR = OUTPUT_SECOND_COLUMN(FOLDER_NAME,'CUR', 1, i, 2)
    MOM = OUTPUT_SECOND_COLUMN(FOLDER_NAME,'MOM', 1, i, 2)
    xx, yy, _, _, _, _, _ = BILNEAR_CURVE(CUR, -MOM, 2)
    DEMAND_CURVATURE = MAXABS_FUN(DATA_FILE04, 1, i)# CURVATURE DYNAMIC ANALYSIS
    DIc = (DEMAND_CURVATURE - xx[1]) / (xx[2] - xx[1]) # Target Section Ductility Damage Index
    diC.append(DIc)
    # AXIAL DUCTILITY DAMAGE INDEX
    DEMAND_STRAIN = max(abs(OUTPUT_SECOND_COLUMN(FOLDER_NAME,'fiberReb_StressStrain', 1, i, 2)))
    DI_AXIAL_S = (DEMAND_STRAIN - ys_strain) / (us_strain - ys_strain)# Axial Ductility Damage Index for Steel Rebar Fiber
    #print('DI REBAR STRAIN AXIAL', DI_AXIAL) 
    diAS.append(DI_AXIAL_S)
    DEMAND_STRAIN = max(abs(OUTPUT_SECOND_COLUMN(FOLDER_NAME,'fiberCon_StressStrain', 1, i, 2)))
    DI_AXIAL_C = (DEMAND_STRAIN - yc_strain) / (uc_strain - yc_strain)# Axial Ductility Damage Index for Concrete Fiber
    #print('DI CONCRETE STRAIN AXIAL', DI_AXIAL)
    diAC.append(DI_AXIAL_C)
    # Load Moment from Dynamic Analysis & Ultimate Moment Capacity from Moment Curvature Analysis
    baseMOM01 = OUTPUT_SECOND_COLUMN(FOLDER_NAME, 'BTH_DYN_01', 3, i, 2) # Reading base moment from Text file - DYNAMIC - NODE 1
    moml.append(baseMOM01)# Applied Moment 
    momr.append(yy[2])
    
    print(f'{i+1} T: {T:.5e} pgd: {pgd[-1]:.5e} pgv: {pgv[-1]:.5e} pga: {pga[-1]:.5e} pgb: {pgb[-1]:.5e}')

#%%------------------------------------------------------------------
xlabel = 'Period'
y1label = 'Displacement Response Spectrum'
y2label = 'Velocity Response Spectrum'
y3label = 'Acceleration Response Spectrum'
y4label = 'Base-shear Response Spectrum'
PLOT_TIME_HIS(period, xlabel, pgd, y1label, pgv, y2label, pga, y3label, pgb, y4label, Z=4, LOG=0)
#%%------------------------------------------------------------------
xlabel = 'Period'
y1label = 'Structure Ductility Damage Index Response Spectrum'
y2label = 'Section Ductility Damage Index Response Spectrum'
y3label = 'Rebar Axial Ductility Damage Index Response Spectrum'
y4label = 'Concrete Axial Ductility Damage Index Response Spectrum'
PLOT_TIME_HIS(period, xlabel, diS, y1label, diC, y2label, diAS, y3label, diAC, y4label, Z = 4, LOG = 0)
#%%------------------------------------------------------------------
xlabel = 'Period'
y1label = 'Ratio Base-shear Dynamic to Base-shear Pushover'
y2label = 'Ratio Displacement Dynamic to Displacement Pushover'
PLOT_TIME_HIS(period, xlabel, bDbP, y1label, dDdP, y2label, _, _, _, _, Z = 2, LOG = 0)
#%%------------------------------------------------------------------
HISROGRAM_BOXPLOT(pgd, HISTO_COLOR='blue', LABEL='Displacement Response Spectrum')
HISROGRAM_BOXPLOT(pgv, HISTO_COLOR='purple', LABEL='Velocity Response Spectrum')
HISROGRAM_BOXPLOT(pga, HISTO_COLOR='green', LABEL='Acceleration Response Spectrum')
HISROGRAM_BOXPLOT(pgb, HISTO_COLOR='orange', LABEL='Base-shear Response Spectrum')
HISROGRAM_BOXPLOT(diS, HISTO_COLOR='lime', LABEL='Structure Ductility Damage Index Response Spectrum')
HISROGRAM_BOXPLOT(diC, HISTO_COLOR='pink', LABEL='Section Ductility Damage Index Response Spectrum')
HISROGRAM_BOXPLOT(diAS, HISTO_COLOR='cyan', LABEL='Rebar Axial Ductility Damage Index Response Spectrum')
HISROGRAM_BOXPLOT(diAC, HISTO_COLOR='yellow', LABEL='Concrete Axial Ductility Damage Index Response Spectrum')
#%%------------------------------------------------------------------
XLABEL = 'Displacement Response Spectrum'
YLABEL = 'Acceleration Response Spectrum'
TITLE = '2D CHART'
COLOR = 'purple'
PLOT_SCATTER(pgd, pga, XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 5) 
#plot_scatter_plotly(pgd, pga, XLABEL, YLABEL, TITLE, COLOR)
#%%------------------------------------------------------------------
XLABEL = 'Displacement Response Spectrum'
YLABEL = 'Velocity Response Spectrum'
ZLABEL = 'Acceleration Response Spectrum'
TITLE = '3D CHART'
PLOT3D(pgd, pgv, pga, XLABEL, YLABEL, ZLABEL, TITLE)
#%%------------------------------------------------------------------
XLABEL = 'Section Ductility Damage Index Response Spectrum'
YLABEL ='Structure Ductility Damage Index Response Spectrum'
TITLE = '2D CHART'
COLOR = 'lightblue'
PLOT_SCATTER(diC, diS, XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 5) 
#plot_scatter_plotly(diC, diS, XLABEL, YLABEL, TITLE, COLOR)
#%%------------------------------------------------------------------
XLABEL = 'Structure Ductility Damage Index Response Spectrum'
YLABEL = 'Displacement Response Spectrum'
TITLE = '2D CHART'
COLOR = 'gray'
PLOT_SCATTER(diS, pgd, XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 5) 
#plot_scatter_plotly(diS, pgd, XLABEL, YLABEL, TITLE, COLOR)
#%%------------------------------------------------------------------
XLABEL = 'Section Ductility Damage Index Response Spectrum'
YLABEL = 'Displacement Response Spectrum'
TITLE = '2D CHART'
COLOR = 'brown'
PLOT_SCATTER(diC, pgd, XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 5) 
#plot_scatter_plotly(diC, pgd, XLABEL, YLABEL, TITLE, COLOR)
#%%------------------------------------------------------------------
XLABEL = 'Section Ductility Damage Index Response Spectrum'
YLABEL ='Structure Ductility Damage Index Response Spectrum'
ZLABEL ='Displacement Response Spectrum'
TITLE = '3D CHART'
PLOT3D(diC, diS, pgd, XLABEL, YLABEL, ZLABEL, TITLE)
#%%------------------------------------------------------------------
### LAST LOAD DATA
import numpy as np
## PYSHOVER ANALYSIS DATA
# Displacement Time History
DTHP = OUTPUT_SECOND_COLUMN(FOLDER_NAME,'DTH_PUSH', 1, i, 2)
# Base Shear Time History 01
BTH01 = OUTPUT_SECOND_COLUMN(FOLDER_NAME,'BTH_PUSH_01', 1, i, 2)
# Base Shear Time History 03
BTH03 = OUTPUT_SECOND_COLUMN(FOLDER_NAME,'BTH_PUSH_03', 1, i, 2)
# Base Shear Time History 05
BTH05 = OUTPUT_SECOND_COLUMN(FOLDER_NAME,'BTH_PUSH_05', 1, i, 2)
# Total Base Shear
BTHP = abs(BTH01 + BTH03 + BTH05)

print(len(DTHP), len(BTHP))

## DYNAMAIC ANALYSIS DATA
# Displacement Time History
DTHD = OUTPUT_SECOND_COLUMN(FOLDER_NAME,'DTH_DYN', 1, i, 2)
# Velocity Time History
VTH = OUTPUT_SECOND_COLUMN(FOLDER_NAME,'VTH_DYN', 1, i, 2)
# Acceleration Time History
ATH = OUTPUT_SECOND_COLUMN(FOLDER_NAME,'ATH_DYN', 1, i, 2)
# Base Shear Time History 01
BTH01 = OUTPUT_SECOND_COLUMN(FOLDER_NAME,'BTH_DYN_01', 1, i, 2)
# Base Shear Time History 03
BTH03 = OUTPUT_SECOND_COLUMN(FOLDER_NAME,'BTH_DYN_03', 1, i, 2)
# Base Shear Time History 05
BTH05 = OUTPUT_SECOND_COLUMN(FOLDER_NAME,'BTH_DYN_05', 1, i, 2)
# Total Base Shear
BTHD = BTH01 + BTH03 + BTH05

print(len(DTHD), len(VTH),len(ATH),len(BTHD))

#%%------------------------------------------------------------------
xx, yy, _, _, _, _, _ = BILNEAR_CURVE(CUR, -MOM, 2)
demand_cur = MAXABS_FUN(DATA_FILE04, 1, i)# DEMAND DYNAMIC CURVATURE
XLABEL = 'Curvature'
YLABEL = 'Moment'
TITLE = 'Last data Moment and Curvature Analysis'
COLOR = 'black'
PLOT_2D(CUR, -MOM, xx, yy, XLABEL, YLABEL, TITLE, COLOR='black', Z=2) 
DIc = (demand_cur - xx[1]) / (xx[2] - xx[1])
print(f'\t\t\t Section Ductility Damage Index: {DIc:.3f}')
#%%------------------------------------------------------------------
xx, yy, _, _, _, _, _ = BILNEAR_CURVE(DTHP, BTHP, 10)
demand_disp = MAXABS_FUN(DATA_FILE01, 1, i)# DEMAND DYNAMIC DISPLACEMENT
XLABEL = 'Displacement'
YLABEL = 'Base Shear'
TITLE = 'Last data Base Shear and Displacement Pushover Analysis'
COLOR = 'black'
PLOT_2D(DTHP, BTHP, xx, yy, XLABEL, YLABEL, TITLE, COLOR='black', Z=2) 
DIs = (demand_disp - xx[1]) / (xx[2] - xx[1])
print(f'\t\t\t Structure Ductility Damage Index: {DIs:.3f}')
#%%------------------------------------------------------------------
xlabel = 'Time'
y1label = 'Displacement Time History'
y2label = 'Velocity Time History'
y3label = 'Acceleration Time History'
y4label = 'Base-shear  Time History'
MAX_TIME = 10 # DATA MUST BE GOTTEN FROM GROUND MOTION TIME HISTORY
DT  = MAX_TIME / len(DTHD)
TIME = np.arange(0, MAX_TIME , DT)
print(len(TIME))
PLOT_TIME_HIS(TIME, xlabel, DTHD, y1label, VTH, y2label, ATH, y3label, BTHD, y4label, Z = 4, LOG = 0)
#%%------------------------------------------------------------------
XLABEL = 'Displacement Time History'
YLABEL = 'Base Shear Time History'
TITLE = 'Last data Base Shear and Displacement Dynamic Analysis Time History'
COLOR = 'black'
PLOT_2D(DTHD, BTHD,_,_, XLABEL, YLABEL, TITLE, COLOR='black', Z=1) 
#%%------------------------------------------------------------------
import pandas as pd
# Create a DataFrame
df = pd.DataFrame({'Displacement Response Spectrum': pgd,
                   'Velocity Response Spectrum': pgv,
                   'Acceleration Response Spectrum': pga,
                   'Base-shear Response Spectrum': pgb,
                   'Structure Ductility Damage Index Response Spectrum': diS,
                   'Section Ductility Damage Index Response Spectrum': diC,
                   'Rebar Axial Ductility Damage Index Response Spectrum':diAS,
                   'Concrete Axial Ductility Damage Index Response Spectrum':diAC})
print(df)
# PLOT HEATMAP FOR CORRELATION 
PLOT_HEATMAP(df) 
#%%------------------------------------------------------------------
### Multiple Regression Model
def Multiple_Regression(df):
    import statsmodels.api as sm
    # Add a constant term for the intercept
    X = sm.add_constant(df[['Displacement Response Spectrum',
                            'Velocity Response Spectrum',
                            'Acceleration Response Spectrum',
                            'Base-shear Response Spectrum',
                            'Section Ductility Damage Index Response Spectrum',
                            'Rebar Axial Ductility Damage Index Response Spectrum',
                            'Concrete Axial Ductility Damage Index Response Spectrum']])

    # Fit the multiple regression model
    model = sm.OLS(df['Structure Ductility Damage Index Response Spectrum'], X).fit()

    # Print the summary
    print(model.summary())

Multiple_Regression(df)   
#%%------------------------------------------------------------------
x = diC # Section Ductility Damage Index Response Spectrum
y = diS # Structure Ductility Damage Index Response Spectrum
predicted_y = PREDICT_LSTM(x, y, look_back = 50, ITERATION = 100)

# Plot the results
import matplotlib.pyplot as plt
plt.figure(figsize=(10, 6))
plt.scatter(x, y, color='blue', marker='+', label='DDI')
plt.scatter(-0.05, predicted_y, color='red', marker='o', label='Predicted next DDI')
plt.title(f'MACHINE LEARNING: LONG SHORT-TREM MEMERY (LSTM) METHOD - Predicted {predicted_y}')
plt.xlabel('Section Ductility Damage Index Response Spectrum')
plt.ylabel('Structure Ductility Damage Index Response Spectrum')
plt.legend()
plt.grid()
plt.show()
#%%------------------------------------------------------------------
x = diS  # Structure Ductility Damage Index Response Spectrum
y = pgd  # Displacement Response Spectrum
predicted_y = PREDICT_LSTM(x, y, look_back = 50, ITERATION = 100)

# Plot the results
import matplotlib.pyplot as plt
plt.figure(figsize=(10, 6))
plt.scatter(x, y, color='blue', marker='+', label='DDI')
plt.scatter(-0.044, predicted_y, color='red', marker='o', label='Predicted next DDI')
plt.title(f'MACHINE LEARNING: LONG SHORT-TREM MEMERY (LSTM) METHOD - Predicted {predicted_y}')
plt.xlabel('Structure Ductility Damage Index Response Spectrum')
plt.ylabel('Displacement Response Spectrum')
plt.legend()
plt.grid()
plt.show()
#%%------------------------------------------------------------------
###       STRUCTURAL RELIABILITY ANALYSIS
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt

# Given data (mean and standard deviation)
mean_applied_moment = np.mean(np.abs(moml))  # Mean Applied Moment [DEMAND]
std_applied_moment = np.std(moml)    # Std Applied Moment  [DEMAND]
mean_resistance_moment = np.mean(momr) # Mean Resistance Moment [CAPACITY]
std_resistance_moment = np.std(momr) # Mean Resistance Moment [CAPACITY]

# Calculate reliability index (beta)
g_mean = mean_resistance_moment - mean_applied_moment
g_std = np.sqrt(std_applied_moment**2 + std_resistance_moment**2)
beta = g_mean / g_std

# Calculate failure probability
P_f = norm.cdf(-beta)

print(f"Mean Applied Moment: {mean_applied_moment:.4f}")
print(f"Std Applied Moment: {std_applied_moment:.4f}")
print(f"Mean Resistance Moment: {mean_resistance_moment:.4f}")
print(f"Std Resistance Moment: {std_resistance_moment:.4f}")
print(f"Reliability index (beta): {beta:.4f}")
print(f"Failure probability (P_f): {100 * P_f:.2f} ٪")

# Plot reliability histogram
x = np.random.normal(mean_applied_moment, std_applied_moment, 1000)
y = np.random.normal(mean_resistance_moment, std_resistance_moment, 1000)
MIX_HISTOGRAM(x, y, BINS=100, X='Applied Moment [DEMAND]', Y='Resistance Moment [CAPACITY]', TITLE=f'Applied & Resistance Moment PDF based on Failure probability: {100 * P_f:0.3f} %')



# Plot reliability diagram
beta_values = np.linspace(-3, 3, 100)
failure_probs = norm.cdf(-beta_values)

plt.figure(figsize=(8, 6))
plt.plot(beta_values, failure_probs, label="Reliability Diagram")
plt.xlabel("Reliability Index (beta)")
plt.ylabel("Failure Probability")
plt.title("Reliability Analysis")
plt.grid(True)
plt.legend()
#plt.semilogx();plt.semilogy();
plt.show()

#%%------------------------------------------------------------------
# If we Calculate Reliability Index and Mean Applied Moment Based on Failure Probability
failure_probability = 0.15  # Set your desired failure probability
root = Normal_CDF_Newton_Raphson(failure_probability)

print(f"Reliability Index (beta): {root:.6f}")

# Calculate Mean Applied Moment
mean_applied_moment = root  * g_std - mean_resistance_moment 
print(f"Mean Applied Moment {mean_applied_moment:.3f} Based on failure probability {100 * failure_probability:.3f} %\n\n")



import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt

# Given data (mean and standard deviation)
mean_applied_moment = abs(mean_applied_moment)  # Mean Applied Moment 
std_applied_moment = 0.15 * mean_applied_moment# Std Applied Moment 
mean_resistance_moment = np.mean(momr) # Mean Resistance Moment 
std_resistance_moment = np.std(momr) # Mean Resistance Moment 

# Calculate reliability index (beta)
g_mean = mean_resistance_moment - mean_applied_moment
g_std = np.sqrt(std_applied_moment**2 + std_resistance_moment**2)
beta = g_mean / g_std

# Calculate failure probability
P_f = norm.cdf(-beta)

print(f"Mean Applied Moment: {mean_applied_moment:.4f}")
print(f"Std Applied Moment: {std_applied_moment:.4f}")
print(f"Mean Resistance Moment: {mean_resistance_moment:.4f}")
print(f"Std Resistance Moment: {std_resistance_moment:.4f}")
print(f"Reliability index (beta): {beta:.4f}")
print(f"Failure probability (P_f): {100 * P_f:.2f} ٪")

# Plot reliability histogram
x = np.random.normal(mean_applied_moment, std_applied_moment, 1000)
y = np.random.normal(mean_resistance_moment, std_resistance_moment, 1000)
MIX_HISTOGRAM(x, y, BINS=100, X='Applied Moment [DEMAND]', Y='Resistance Moment [CAPACITY]', TITLE=f'Applied & Resistance Moment PDF based on Failure probability: {100 * P_f:0.3f} %')
#%%------------------------------------------------------------------
