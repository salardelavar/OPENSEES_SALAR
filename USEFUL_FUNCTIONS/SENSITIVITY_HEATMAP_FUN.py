def SENSITIVITY_HEATMAP_FUN(X, Y, x_labels=None, y_labels=None, method='pearson', cmap='coolwarm', annot=True, fmt='.2f'):
    """
    Plots a heatmap of sensitivity coefficients (correlation or SRC) between inputs X and outputs Y.

    Parameters:
    X : numpy array or DataFrame (n_samples, n_inputs=5)
    Y : numpy array or DataFrame (n_samples, n_outputs=5)
    x_labels : list, names of input variables
    y_labels : list, names of output variables
    method : str, 'pearson' for correlation, 'src' for standardized regression coefficients
    THIS PYTHON SCRIPT WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
    """
    import numpy as np
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    # Convert to DataFrames
    if not isinstance(X, pd.DataFrame):
        if x_labels is None:
            x_labels = [f'Input_{i+1}' for i in range(X.shape[1])]
        X = pd.DataFrame(X, columns=x_labels)
    if not isinstance(Y, pd.DataFrame):
        if y_labels is None:
            y_labels = [f'Output_{j+1}' for j in range(Y.shape[1])]
        Y = pd.DataFrame(Y, columns=y_labels)

    if method == 'pearson':
        # Simple correlation
        combined = pd.concat([X, Y], axis=1)
        corr_matrix = combined.corr()
        # Extract the X-Y block
        sensitivity_matrix = corr_matrix.loc[X.columns, Y.columns]
    elif method == 'src':
        # Standardized Regression Coefficients
        from sklearn.preprocessing import StandardScaler
        from sklearn.linear_model import LinearRegression
        scaler_X = StandardScaler()
        scaler_Y = StandardScaler()
        X_scaled = scaler_X.fit_transform(X)
        Y_scaled = scaler_Y.fit_transform(Y)
        
        sensitivity_matrix = pd.DataFrame(index=X.columns, columns=Y.columns)
        for j, y_name in enumerate(Y.columns):
            model = LinearRegression().fit(X_scaled, Y_scaled[:, j])
            # coefficients are the beta weights (SRC)
            sensitivity_matrix[y_name] = model.coef_
    else:
        raise ValueError("method must be 'pearson' or 'src'")

    # Plot
    plt.figure(figsize=(8, 6))
    sns.heatmap(sensitivity_matrix, annot=annot, fmt=fmt, cmap=cmap, 
                center=0, linewidths=0.5, cbar_kws={'label': 'Sensitivity Coefficient'})
    plt.title(f'Sensitivity Heatmap ({method.upper()})')
    plt.tight_layout()
    plt.show()
    
    return sensitivity_matrix
"""
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
# --- Example Usage ---
np.random.seed(42)
n_samples = 1000
X = np.random.uniform(-1, 1, size=(n_samples, 5))

# Define complex relationships for 5 outputs
Y = np.zeros((n_samples, 5))
Y[:, 0] = 2.0 * X[:, 0] + 0.5 * X[:, 1] + 0.1 * X[:, 2] - 0.2 * X[:, 3] + 0.3 * X[:, 4] + 0.5 * np.random.randn(n_samples)
Y[:, 1] = 0.8 * X[:, 3] - 0.9 * X[:, 4] + 0.1 * X[:, 0] + 0.2 * np.random.randn(n_samples)
Y[:, 2] = 1.5 * X[:, 1] + 1.5 * X[:, 2] + 0.3 * X[:, 0] - 0.3 * X[:, 3] + 0.1 * np.random.randn(n_samples)
Y[:, 3] = 3.0 * X[:, 4] + 0.2 * X[:, 1] - 0.1 * X[:, 2] + 0.5 * np.random.randn(n_samples)
Y[:, 4] = 0.5 * X[:, 0] + 0.5 * X[:, 1] + 0.5 * X[:, 2] + 0.5 * X[:, 3] + 0.5 * X[:, 4] + 0.5 * np.random.randn(n_samples)

x_labs = ['X1', 'X2', 'X3', 'X4', 'X5']
y_labs = ['Y1', 'Y2', 'Y3', 'Y4', 'Y5']

# Call the function
coeffs = PLOT_SENSITIVITY_HEATMAP_FUN(X, Y, x_labels=x_labs, y_labels=y_labs, method='pearson')
print("\nSensitivity Coefficient Matrix:\n", coeffs)
"""