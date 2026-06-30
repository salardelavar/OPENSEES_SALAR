def PLOT_CONTOUR_3D_2D_FUN(X, Y, Z, XLABEL, YLABEL, ZLABEL):
    # THIS PYTHON SCRIPT WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D  # for 3D plotting
    fig1 = plt.figure(figsize=(10, 6))
    ax1 = fig1.add_subplot(111, projection='3d')
    surf1 = ax1.plot_surface(X, Y, Z, cmap='viridis', edgecolor='none')
    ax1.set_xlabel(XLABEL)
    ax1.set_ylabel(YLABEL)
    ax1.set_zlabel(ZLABEL)
    ax1.set_title(f'3D CONTOUR: {XLABEL} - {YLABEL} - {ZLABEL}')
    fig1.colorbar(surf1, shrink=0.5, aspect=10)
    plt.show()

    fig, ax = plt.subplots(figsize=(7, 5))
    cf = ax.contourf(X, Y, Z, levels=20, cmap='viridis')
    cs = ax.contour(X, Y, Z, levels=10, colors='white', linewidths=0.5)
    ax.clabel(cs, inline=True, fontsize=8)
    ax.set_xlabel(XLABEL)
    ax.set_ylabel(YLABEL)
    ax.set_title(f'2D CONTOUR: {XLABEL} – Topographic Map')
    fig.colorbar(cf, ax=ax, label='Damping Ratio')
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


