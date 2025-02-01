import numpy as np
# -------------------------------------------------     
    
def BETA_PDF(MIN_X, MAX_X, a, b, n):
    return MIN_X + (MAX_X - MIN_X) * np.random.beta(a, b, n)
    
# -------------------------------------------------   
 
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

# -----------------------------------------------
     
"""
Long short-term memory (LSTM) is a type
of recurrent neural network (RNN) aimed
at dealing with the vanishing gradient
problem present in traditional RNNs
"""
def PREDICT_LSTM(x, y, look_back, ITERATION, XLABEL, YLABEL):
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
    
    # Plot the results
    import matplotlib.pyplot as plt
    plt.figure(figsize=(10, 6))
    plt.scatter(x, y, color='blue', marker='+', label='DDI')
    plt.scatter(next_x[-1], predicted_y, color='red', marker='o', label='Predicted next DDI')
    plt.title(f'MACHINE LEARNING: LONG SHORT-TREM MEMERY (LSTM) METHOD - Predicted {predicted_y}')
    plt.xlabel(XLABEL)
    plt.ylabel(YLABEL)
    plt.legend()
    plt.grid()
    plt.show()

# -----------------------------------------------  
  
def RELIABILITY_ANALYSIS(base_reaction, num_sim, mean_capacity, std_dev_capacity):
    """
    Perform reliability analysis for base reaction and element capacity.

    Parameters:
    - base_reaction (list or np.array): Array of maximum base reactions (demand) from simulations.
    - num_sim (int): Number of Monte Carlo simulations.
    - mean_capacity (float): Mean element capacity (resistance) in Newtons.
    - std_dev_capacity (float): Standard deviation of element capacity in Newtons.

    Returns:
    - probability_of_failure (float): Estimated probability of failure.
    - reliability_index (float): Calculated reliability index (β).
    """

    # Element capacity as a random variable
    element_capacity = np.random.normal(mean_capacity, std_dev_capacity, num_sim)

    # Convert base reaction to a numpy array (if not already)
    base_reaction = np.array(base_reaction)

    # Ensure base reaction matches simulation size
    if len(base_reaction) != num_sim:
        raise ValueError("Length of base_reaction must match the number of simulations (num_sim).")

    # Limit state function: R - D
    limit_state = element_capacity - base_reaction

    # Reliability analysis
    failures = np.sum(limit_state <= 0)
    probability_of_failure = failures / num_sim

    # Reliability index calculation
    mu_g = np.mean(limit_state)
    element_capacity_std = np.std(element_capacity)
    base_reaction_std = np.std(base_reaction)
    sigma_g = np.sqrt(element_capacity_std ** 2 +  base_reaction_std**2)
    reliability_index = mu_g / sigma_g
    
    # Visualization
    import matplotlib.pyplot as plt
    plt.figure(figsize=(10, 6))
    plt.hist(element_capacity, bins=50, color='blue', alpha=0.7, label="R: Element Capacity")
    plt.hist(base_reaction, bins=50, color='red', alpha=0.7, label="D: Base Reaction")
    plt.title("Reliability Analysis: Limit State Distribution")
    plt.xlabel("Data")
    plt.ylabel("Frequency")
    plt.legend()
    plt.grid(True)
    plt.show()

    # Print results
    print(f"Probability of Failure (P_f): {probability_of_failure:.6f}")
    print(f"Reliability Index (β): {reliability_index:.2f}")

    return probability_of_failure, reliability_index

# ----------------------------------------------- 

def MULTIPLE_REGRESSION(df):
    import statsmodels.api as sm
    # Add a constant term for the intercept
    X = sm.add_constant(df[['Max_velocity',
                            'Max_acceleration',
                            'Max_Base_Reaction',
                            'Damping_Ratio']])

    # Fit the multiple regression model
    model = sm.OLS(df['Max_displacement'], X).fit()

    # Print the summary
    print(model.summary())
    
# -----------------------------------------------     

def RANDOM_FOREST(df):
    from sklearn.model_selection import train_test_split
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.metrics import classification_report, mean_squared_error, r2_score
    import seaborn as sns
    # Target: Safety binary label (example: 1 if max displacement <= threshold)
    threshold_displacement = 0.01  # Example threshold
    df['Safety_label'] = (df['Max_displacement'] <= threshold_displacement).astype(int)

    # Step 1: Create a dataset with statistical features and response labels
    # Features: Extracted stats
    # Create a dataset with individual simulation results
    """
    data = {
        'Max_displacement': max_displacement,
        'Max_velocity': max_velocity,
        'Max_acceleration': max_acceleration,
        'Max_Base_Reaction': max_base_reaction,
        'Damping_Ratio': DR,
    }

    # Convert to DataFrame
    df = pd.DataFrame(data)
    """
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
    
# -----------------------------------------------   

def PLOT_TIME_HISTORY(time, displacement, velocity, acceleration, base_reaction):
    
    import matplotlib.pyplot as plt
    plt.figure(figsize=(18, 20))

    # Displacement
    plt.subplot(4, 1, 1)
    plt.plot(time, displacement, label='Displacement')
    plt.xlabel('Time [s]')
    plt.ylabel('Displacement [m]')
    plt.title('Displacement Time History')
    plt.grid(True)

    # Velocity
    plt.subplot(4, 1, 2)
    plt.plot(time, velocity, label='Velocity', color='orange')
    plt.xlabel('Time [s]')
    plt.ylabel('Velocity [m/s]')
    plt.title('Velocity Time History')
    plt.grid(True)

    # Acceleration
    plt.subplot(4, 1, 3)
    plt.plot(time, acceleration, label='Acceleration', color='green')
    plt.xlabel('Time [s]')
    plt.ylabel('Acceleration [m/s²]')
    plt.title('Acceleration Time History')
    plt.grid(True)

    # Base Reaction Force
    plt.subplot(4, 1, 4)
    plt.plot(time, base_reaction, label='Base Reaction', color='red')
    plt.xlabel('Time [s]')
    plt.ylabel('Base Reaction [N]')
    plt.title('Base Reaction Time History')
    plt.grid(True)

    plt.tight_layout()
    plt.show()
    
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
# ----------------------------------------------- 
# Cluster the Data
def CLUSTER_DATA(X, Y, XLABEL, YLABEL, MAX_CLUSTERS):
    import numpy as np
    from sklearn.preprocessing import StandardScaler
    from sklearn.cluster import KMeans
    from sklearn.linear_model import LinearRegression
    import matplotlib.pyplot as plt
    # Combine x and y into a single dataset and scale
    x = np.array(X);y = np.array(Y);
    data = np.column_stack((x, y))
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(data)

    best_r2 = -np.inf
    best_k = 0
    best_clusters = None
    global_mean = np.mean(y)  # For handling small clusters

    for k in range(1, MAX_CLUSTERS+1):
        kmeans = KMeans(n_clusters=k, random_state=42)
        clusters = kmeans.fit_predict(scaled_data)
        
        total_ss_res = 0.0
        for i in range(k):
            cluster_mask = (clusters == i)
            x_cluster = x[cluster_mask].reshape(-1, 1)  # Ensure 2D shape
            y_cluster = y[cluster_mask]
            
            if len(x_cluster) < 2:
                # Handle small clusters with global mean
                total_ss_res += np.sum((y_cluster - global_mean) ** 2)
                continue
                
            model = LinearRegression()
            model.fit(x_cluster, y_cluster)
            y_pred = model.predict(x_cluster)
            total_ss_res += np.sum((y_cluster - y_pred) ** 2)

        ss_tot = np.sum((y - global_mean) ** 2)
        r2 = 1 - (total_ss_res / ss_tot)
        
        if r2 > best_r2:
            best_r2 = r2
            best_k = k
            best_clusters = clusters

    print(f"Optimal clusters: {best_k}, R²: {best_r2:.4f}")

    # Plot clustered data
    plt.figure(figsize=(10, 6))
    for i in range(best_k):
        cluster_mask = (best_clusters == i)
        plt.scatter(x[cluster_mask], y[cluster_mask], alpha=0.5, label=f'Cluster {i+1}')
    plt.title(f"{XLABEL} vs {YLABEL} - {best_k} Clusters (R² = {best_r2:.4f})")
    plt.xlabel(XLABEL)
    plt.ylabel(YLABEL)
    plt.legend()
    plt.grid(True)
    plt.show()
# -----------------------------------------------     
    