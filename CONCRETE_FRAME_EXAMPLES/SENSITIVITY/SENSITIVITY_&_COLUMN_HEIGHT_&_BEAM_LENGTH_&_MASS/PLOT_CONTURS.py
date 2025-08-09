import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

# Define the plotting function
def PLOT_3D(TAG, X, Y, Z, XLABEL, YLABEL, ZLABEL):
    # Convert to NumPy arrays
    X = np.array(X)
    Y = np.array(Y)
    Z = np.array(Z)
    
    # Create grid for contour plot
    xi = np.linspace(min(X), max(X), 100)
    yi = np.linspace(min(Y), max(Y), 100)
    xi, yi = np.meshgrid(xi, yi)
    
    # Interpolate Z values on the grid
    zi = griddata((X, Y), Z, (xi, yi), method='cubic')
    
    # Plot 3D contour
    fig = plt.figure(TAG, figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')
    contour = ax.plot_surface(xi, yi, zi, cmap='viridis', edgecolor='none')
    
    ax.set_xlabel(XLABEL)
    ax.set_ylabel(YLABEL)
    ax.set_zlabel(ZLABEL)
    
    fig.colorbar(contour, ax=ax, shrink=0.5, aspect=5)
    plt.title(f'3D Contour Plot of {ZLABEL}')
    plt.show()

# Generate random test data
np.random.seed(0)
DIAc_MAX = np.random.choice([8, 10, 12, 14, 16, 18, 20, 22, 25, 28, 30, 32], 100)
Kfc_MAX = np.random.choice([1, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40, 1.45, 1.50, 1.55], 100)
DISP_X_MAX = np.random.uniform(low=0.5, high=15.0, size=100)  # Synthetic displacement values

# Labels
XLABEL = 'Rebar Diameter'
YLABEL = 'Strength Enhancement Factor'
ZLABEL = 'Structure Displacement in X Dir. [mm]'

# Plot
PLOT_3D(1, DIAc_MAX, Kfc_MAX, DISP_X_MAX, XLABEL, YLABEL, ZLABEL)
