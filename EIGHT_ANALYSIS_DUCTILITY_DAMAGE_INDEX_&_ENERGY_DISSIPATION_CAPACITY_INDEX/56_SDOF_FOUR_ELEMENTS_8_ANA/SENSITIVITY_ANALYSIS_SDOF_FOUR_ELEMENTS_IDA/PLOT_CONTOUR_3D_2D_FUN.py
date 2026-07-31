def PLOT_CONTOUR_3D_2D_FUN(TAG, X, Y, Z, XLABEL, YLABEL, ZLABEL):
    # THIS PYTHON SCRIPT WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.interpolate import griddata
    # Convert to NumPy arrays
    X = np.array(X)
    Y = np.array(Y)
    Z = np.array(Z)
    
    # Create grid for contour plot
    xi = np.linspace(min(X), max(X), 100)
    yi = np.linspace(min(Y), max(Y), 100)
    xi, yi = np.meshgrid(xi, yi)
    
    # Interpolate Z values on the grid
    #zi = griddata((X, Y), Z, (xi, yi), method='nearest')
    zi = griddata((X, Y), Z, (xi, yi), method='linear')
    #zi = griddata((X, Y), Z, (xi, yi), method='cubic')
    
    # Plot 3D contour
    fig = plt.figure(TAG, figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')
    contour = ax.plot_surface(xi, yi, zi, cmap='viridis', edgecolor='none')
    
    ax.set_xlabel(XLABEL)
    ax.set_ylabel(YLABEL)
    ax.set_zlabel(ZLABEL)
    
    fig.colorbar(contour, ax=ax, shrink=0.5, aspect=5)
    plt.title(f'3D CONTOUR: {XLABEL} - {YLABEL} - {ZLABEL}')
    plt.show()


    fig, ax = plt.subplots(figsize=(10, 7))
    cf = ax.contourf(xi, yi, zi, levels=20, cmap='viridis')
    cs = ax.contour(xi, yi, zi, levels=10, colors='white', linewidths=0.5)
    ax.clabel(cs, inline=True, fontsize=8)
    ax.set_xlabel(XLABEL)
    ax.set_ylabel(YLABEL)
    ax.set_title(f'2D CONTOUR: {XLABEL} – Topographic Map')
    fig.colorbar(cf, ax=ax, label=f'{ZLABEL}')
    plt.tight_layout()
    plt.show()

"""
import numpy as np
# --- Create sample data (replace with your own Kd_grid, DRx_grid, etc.) ---
nK, nD = 20, 30                     # grid size
kd_vals = np.linspace(1000, 5000, nK)
drx_vals = np.linspace(0.05, 0.2, nD)
Kd_grid, DRx_grid = np.meshgrid(kd_vals, drx_vals, indexing='ij')

# Example Z: a simple function of Kd and DRx (e.g., damping ratio = Kd/1e4 * exp(-DRx))
Z = (Kd_grid / 10000) * np.exp(-DRx_grid * 5)

# --- Call your function ---
PLOT_CONTOUR_3D_2D_FUN(
    Kd_grid, DRx_grid, Z,
    'Damper Stiffness [N/m] (Kd)',
    'Damper Damping Ratio (DRx)',
    'Damping Ratio'
)

# Damping ratio
PLOT_CONTOUR_3D_2D_FUN(
    Kd_grid, DRx_grid, WW_matrix,
    'Damper Elastic Stiffness [N/m] (Kd)',
    'Damper Damping Ratio (DRx)',
    'Structural Damping Ratio'
)

# Structural period
PLOT_CONTOUR_3D_2D_FUN(
    Kd_grid, DRx_grid, ZZ_matrix,
    'Damper Elastic Stiffness [N/m] (Kd)',
    'Damper Damping Ratio (DRx)',
    'Structural Period (s)'
)
"""


