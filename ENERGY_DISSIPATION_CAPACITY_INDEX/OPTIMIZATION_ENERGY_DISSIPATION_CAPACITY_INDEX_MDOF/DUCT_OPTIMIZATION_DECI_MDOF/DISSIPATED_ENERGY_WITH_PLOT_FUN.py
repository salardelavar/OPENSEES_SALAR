# EVALUATION OF DISSIPATED ENERGY CAPACITY INDEX
# THIS PYTHON SCRIPT WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
def DISSIPATED_ENERGY_FUN_WITH_PLOT(displacement, base_shear, method, title="Hysteresis Curve"):
    if method == 1:
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
        ax.set_title(f"{title} - (Convex Hull)")
        ax.set_xlabel("Displacement (m)")
        ax.set_ylabel("Base Shear (N)")
        ax.grid(True, linestyle='--', alpha=0.5)
        ax.legend()
        
    if method == 2:  
        import numpy as np
        import matplotlib.pyplot as plt
        # Data preparation
        disp = np.asarray(displacement, dtype=float)
        shear = np.asarray(base_shear, dtype=float)

        if disp.size != shear.size:
            raise ValueError("Displacement and base shear arrays must have the same length.")
        if disp.size < 3:
            raise ValueError("At least 3 points are required to form a closed loop.")

        # Close the loop if not already closed (important for shoelace)
        if not (disp[0] == disp[-1] and shear[0] == shear[-1]):
            disp = np.append(disp, disp[0])
            shear = np.append(shear, shear[0])

        # Dissipated energy (E_d) via Shoelace formula
        x = disp
        y = shear
        area = 0.5 * np.abs(np.dot(x[:-1], y[1:]) - np.dot(y[:-1], x[1:]))
        
        # Plotting
        fig, ax = plt.subplots(figsize=(7, 6))
    
        # Hysteresis curve
        idx_max = np.argmax(np.abs(disp))
        ax.plot(disp, shear, 'k-', linewidth=1.2, label="Hysteresis Loop")
        ax.scatter(disp[idx_max], shear[idx_max], color='blue', s=80,
                   zorder=5, label=r"$(u_{\rm max}, F_{\rm max})$")
    
        # Fill the enclosed area
        ax.fill(disp, shear, color='red', alpha=0.25, label=f"E$_d$ = {area:.3f} N·m")
    
    
        ax.set_title(title)
        ax.set_xlabel("Displacement (m)")
        ax.set_ylabel("Base Shear (N)")
        ax.grid(True, linestyle='--', alpha=0.5)
        ax.legend(loc='lower right')
        fig.tight_layout()
    
        return area, fig

def ENERGY_DISSIPATION_CAPACITY_INDEX(dispX_SEI, reactionX_SEI, dispX_CP, reactionX_CP):
    Ed_SEI, fig_SEI = DISSIPATED_ENERGY_FUN_WITH_PLOT(
        dispX_SEI, reactionX_SEI, method = 2, 
        title="Earthquake Response – Dissipated Energy (Convex Hull)"
    )
    fig_SEI.show()
    
    print(f"Dissipated Energy from Earthquake= {Ed_SEI:.2f} N·m")
    
    Ed_CP, fig_CP = DISSIPATED_ENERGY_FUN_WITH_PLOT(
        dispX_CP, reactionX_CP, method = 2,
        title="Cyclic Loading – Dissipated Energy (Convex Hull)"
    )
    fig_CP.show()
    
    print(f"Dissipated Energy from Cyclic Displacement= {Ed_CP:.2f} N·m")
    
    
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
    
    return  DECI   