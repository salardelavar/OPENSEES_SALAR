# THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI) 
def DISSIPATED_ENERGY_FUN(displacement, base_shear):
    import numpy as np
    from scipy.spatial import ConvexHull
    """
    Calculate area of the outer envelope (convex hull)
    of base shear - displacement curve.

    Parameters:
    displacement : array-like
        Displacement history
    base_shear : array-like
        Base shear history

    Returns:
    area : float
        Outer envelope area (pure geometric area)
    """
    points = np.column_stack((displacement, base_shear))
    hull = ConvexHull(points)

    return hull.volume   # In 2D, 'volume' means area

#%%----------------------------------------------
# Example hysteresis data

#disp = [0, 5, 10, 5, 0, -5, -10, -5, 0, 10, 15, 20]
#shear = [0, 50, 80, 60, 0, -60, -80, -50, 0, 100, 150, 200]
disp = [-5, 5, 5, -5]
shear = [-5, -5, 5, 5]
Ed = DISSIPATED_ENERGY_FUN(disp, shear)
print(f"Dissipated Energy = {Ed:.2f} kN·mm")

def PLOT(XDATA, YDATA, TITLE, XLABEL, YLABEL, COLOR, SEMILOGY):
    import matplotlib.pyplot as plt
    fig = plt.figure(-1, figsize=(12, 8))
    plt.plot(XDATA, YDATA, color=COLOR, linewidth=2)
    plt.title(TITLE)
    plt.ylabel(YLABEL)
    plt.xlabel(XLABEL)
    if SEMILOGY == True:
        plt.semilogy()
    plt.grid()
    
XDATA = disp
YDATA = shear
XLABEL = 'Displacement in Middle Span [mm]'
YLABEL = 'Base Reaction [N]'
TITLE = 'Base Reaction and Displacement of Structure During Pushover Analysis'
COLOR = 'black'
SEMILOGY = False
PLOT(XDATA, YDATA, TITLE, XLABEL, YLABEL, COLOR, SEMILOGY)    