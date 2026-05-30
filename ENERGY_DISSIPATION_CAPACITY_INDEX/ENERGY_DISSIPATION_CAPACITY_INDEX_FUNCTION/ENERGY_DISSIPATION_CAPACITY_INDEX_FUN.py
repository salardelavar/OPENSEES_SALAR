"""
Energy Dissipation Capacity Index (EDCI):
The Energy Dissipation Capacity Index is a quantitative measure used in structural
 engineering to evaluate how effectively a structural element (e.g., a beam, column,
 shear wall, or connection) can absorb and dissipate energy during seismic loading
 compared to its performance under controlled cyclic displacement loading.
It compares the actual energy absorbed during an earthquake with the maximum energy
 dissipation capacity that the component demonstrates in a laboratory‑style cyclic test.

Why This Index Is Important:
During an earthquake, structures undergo repeated cycles of deformation.
 A system with high energy dissipation capacity can withstand more damage
 without collapsing because it can convert seismic input energy into hysteretic energy, not elastic rebound.

The EDCI helps engineers understand:
Ductility performance
Hysteretic behavior
Damage tolerance
Collapse prevention capability
It is especially used in performance‑based seismic evaluation and retrofit design.

THIS PYTHON SCRIPT WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI) 
"""

#%%----------------------------------------------
# Example hysteresis data

#disp = [0, 5, 10, 5, 0, -5, -10, -5, 0, 10, 15, 20]
#shear = [0, 50, 80, 60, 0, -60, -80, -50, 0, 100, 150, 200]
disp_CP = [-5, 15, 5, -15]
reaction_CP = [-5, -5, 5, 5]

disp_SEI = [-5, 5, 5, -5]
reaction_SEI = [-5, -5, 5, 5]

# EVALUATION OF DISSIPATED ENERGY CAPACITY INDEX
def DISSIPATED_ENERGY_FUN_WITH_PLOT(displacement, base_shear, title="Hysteresis Curve"):
    """
    Compute dissipated energy using convex hull and plot the hysteresis curve
    with the outer hull area shaded.

    Parameters
    ----------
    displacement : array-like
    base_shear  : array-like
    title       : str

    Returns
    -------
    float
        Area of convex hull (dissipated energy)
    """
    import numpy as np
    from scipy.spatial import ConvexHull
    import matplotlib.pyplot as plt
    displacement = np.asarray(displacement)
    base_shear  = np.asarray(base_shear)

    if displacement.size != base_shear.size:
        raise ValueError("Displacement and base shear arrays must have equal lengths.")

    points = np.column_stack((displacement, base_shear))
    hull = ConvexHull(points)
    area = hull.volume   # 2D hull → area

    fig, ax = plt.subplots(figsize=(7, 6))

    # Plot full hysteresis
    ax.plot(displacement, base_shear, 'k-', linewidth=1, label="Hysteresis Curve")

    # Plot convex hull edges
    hull_pts = points[hull.vertices]
    ax.plot(hull_pts[:, 0], hull_pts[:, 1], 'r--', lw=2, label="Convex Hull")

    # Shade hull area
    ax.fill(hull_pts[:, 0], hull_pts[:, 1], color='red', alpha=0.25, label="Hull Area (Energy)")

    # Labels and style
    ax.set_title(title)
    ax.set_xlabel("Displacement (mm)")
    ax.set_ylabel("Base Shear (N)")
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.legend()

    return area, fig

Ed_SEI, fig_SEI = DISSIPATED_ENERGY_FUN_WITH_PLOT(
    disp_SEI, reaction_SEI, 
    title="Earthquake Response – Dissipated Energy (Convex Hull)"
)
fig_SEI.show()

print(f"Dissipated Energy from Earthquake= {Ed_SEI:.2f} N·mm")

Ed_CP, fig_CP = DISSIPATED_ENERGY_FUN_WITH_PLOT(
    disp_CP, reaction_CP,
    title="Cyclic Loading – Dissipated Energy (Convex Hull)"
)
fig_CP.show()

print(f"Dissipated Energy from Cyclic Displacement= {Ed_CP:.2f} N·mm")


DECI = 100 * Ed_SEI / Ed_CP

if DECI <= 100:
    print(f'\n\tDISSIPATED ENERGY CAPACITY INDEX: {DECI:.3f} [%]\n')
else:
    print('\n\tFOR EVALUATION OF DISSIPATED ENERGY CAPACITY INDEX:')
    print('\n\tCHECK THE CYCLIC DISPLACEMENT ANALYSIS AND IF IT IS POSSIBLE')
    print('\t\t\tINCREASE THE DISPLACEMENT.\n')

if DECI <= 0:
    print("\n\tZONE 0: NO DAMAGE\n")
elif DECI > 0 and DECI <= 10:
    print("\n\tZONE 1: VERY MINOR DAMAGE\n")
elif DECI > 10 and DECI <= 20:
    print("\n\tZONE 2: MINOR DAMAGE\n")
elif DECI > 20 and DECI <= 30:
    print("\n\tZONE 3: MODERATE–LOW DAMAGE\n")
elif DECI > 30 and DECI <= 40:
    print("\n\tZONE 4: MODERATE DAMAGE\n")
elif DECI > 40 and DECI <= 50:
    print("\n\tZONE 5: MODERATE–HIGH DAMAGE\n")
elif DECI > 50 and DECI <= 60:
    print("\n\tZONE 6: SEVERE–LOW DAMAGE\n")
elif DECI > 60 and DECI <= 70:
    print("\n\tZONE 7: SEVERE–MEDIUM DAMAGE\n")    
elif DECI > 70 and DECI <= 80:
    print("\n\tZONE 8: SEVERE–HIGH DAMAGE\n")
elif DECI > 80 and DECI <= 90:
    print("\n\tZONE 9: VERY SEVERE DAMAGE\n")
elif DECI > 90 and DECI <= 100:
    print("\n\tZONE 10: FAILURE DAMAGE\n")      
#%%----------------------------------------------------
