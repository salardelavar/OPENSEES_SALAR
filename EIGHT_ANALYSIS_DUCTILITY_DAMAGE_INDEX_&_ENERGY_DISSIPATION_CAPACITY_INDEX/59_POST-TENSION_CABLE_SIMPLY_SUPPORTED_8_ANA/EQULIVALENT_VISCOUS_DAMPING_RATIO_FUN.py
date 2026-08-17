#%% EQULIVALENT VISCOUS DAMPING RATIO
"""
[1] The function EQULIVALENT_VISCOUS_DAMPING_RATIO_FUN computes the dissipated energy (E_d) using either the convex‑hull area (method=1) or the shoelace formula (method=2), then calculates the equivalent viscous damping ratio as ζ = 100·E_d / (4π·E_s), where E_s = ½·F_max·u_max.

[2] Both methods produce the same styled plot: the hysteresis loop in black, the enclosed area filled in red, the point of maximum absolute displacement marked in blue, and a text box displaying ζ, E_d, and E_s.

[3] The bottom section generates a synthetic elliptical hysteresis loop (sinusoidal displacement and a phase‑shifted force) and calls the function with method=1 to print the damping ratio and display the plot.
THIS PYTHON SCRIPT WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
"""
def EQULIVALENT_VISCOUS_DAMPING_RATIO_FUN(displacement, base_shear, method, XLABEL, YLABEL, TITLE):
    import numpy as np
    from scipy.spatial import ConvexHull
    import matplotlib.pyplot as plt
    if method == 1:

        """
        Compute the equivalent viscous damping ratio from a hysteresis loop.
    
        The dissipated energy per cycle (E_d) is the area enclosed by the loop,
        computed using the convex hull (Shoelace formula on the hull). The maximum
        strain energy (E_s) is estimated as 0.5 * F_max * u_max, where u_max is the
        maximum absolute displacement and F_max is the absolute base shear at that
        same displacement.
    
        Parameters
        ----------
        displacement : array-like
            Displacement values (assumed to form a closed hysteresis loop).
        base_shear : array-like
            Corresponding base shear values.
        XLABEL : str
            Label for the x-axis.
        YLABEL : str
            Label for the y-axis.
        TITLE : str
            Title for the plot.
    
        Returns
        -------
        zeta_eq : float
            Equivalent viscous damping ratio (%).
        fig : matplotlib.figure.Figure
            Figure object with the plotted hysteresis loop and energy representation.
        
        """
        import numpy as np
        from scipy.spatial import ConvexHull
        import matplotlib.pyplot as plt
        
        # Convert to numpy arrays
        disp = np.asarray(displacement, dtype=float)
        shear = np.asarray(base_shear, dtype=float)
    
        if disp.size != shear.size:
            raise ValueError("Displacement and base shear arrays must have equal lengths.")
        if disp.size < 3:
            raise ValueError("At least 3 points are required to form a loop.")
    
        # ---- Dissipated energy via convex hull (method 1) ----
        points = np.column_stack((disp, shear))
        hull = ConvexHull(points)
        E_d = hull.volume  # area of the convex hull (2D)
    
        # ---- Maximum strain energy ----
        idx_max = np.argmax(np.abs(disp))
        u_max = np.abs(disp[idx_max])
        F_at_max_u = np.abs(shear[idx_max])
        E_s = 0.5 * F_at_max_u * u_max
    
        # ---- Equivalent viscous damping ratio ----
        zeta_eq = 100 * E_d / (4.0 * np.pi * E_s) if E_s != 0 else 0.0
    
        # ---- Plotting (style of method 2) ----
        # Close the loop for proper filling if it is not already closed
        if not (disp[0] == disp[-1] and shear[0] == shear[-1]):
            disp_closed = np.append(disp, disp[0])
            shear_closed = np.append(shear, shear[0])
        else:
            disp_closed = disp
            shear_closed = shear
    
        fig, ax = plt.subplots(figsize=(7, 6))
    
        # Hysteresis curve
        ax.plot(disp, shear, 'k-', linewidth=1.2, label="Hysteresis Loop")
        ax.scatter(disp[idx_max], shear[idx_max], color='blue', s=80,
                   zorder=5, label=r"$(u_{\rm max}, F_{\rm max})$")
    
        # Fill the enclosed area
        ax.fill(disp_closed, shear_closed, color='red', alpha=0.25,
                label=f"E$_d$ = {E_d:.3f} N·mm")
    
        # Reference axes
        ax.axhline(0, color='gray', linestyle='--', linewidth=0.8)
        ax.axvline(0, color='gray', linestyle='--', linewidth=0.8)
    
        # Text box with results
        textstr = (f"$\\xi_{{\\rm eq}} = {zeta_eq:.4f}$ [%]\n"
                   f"$E_d = {E_d:.3f}$ N·mm\n"
                   f"$E_s = {E_s:.3f}$ N·mm")
        props = dict(boxstyle='round', facecolor='white', alpha=0.8)
        ax.text(0.05, 0.95, textstr, transform=ax.transAxes,
                fontsize=11, verticalalignment='top', bbox=props)
    
        ax.set_title(TITLE)
        ax.set_xlabel(XLABEL)
        ax.set_ylabel(YLABEL)
        ax.grid(True, linestyle='--', alpha=0.5)
        ax.legend(loc='lower right')
        fig.tight_layout()
        plt.show()

    if method == 2:
        """
        Compute the equivalent viscous damping ratio from a hysteresis loop.

        The dissipated energy per cycle (E_d) is the area enclosed by the loop,
        computed using the Shoelace formula. The maximum strain energy (E_s) is
        estimated as 0.5 * F_max * u_max, where u_max is the maximum absolute
        displacement and F_max is the absolute base shear at that same displacement.
        """
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
        E_d = 0.5 * np.abs(np.dot(x[:-1], y[1:]) - np.dot(y[:-1], x[1:]))

        # Maximum strain energy (E_s)
        # Find the point of maximum absolute displacement
        idx_max = np.argmax(np.abs(disp))
        u_max = np.abs(disp[idx_max])
        F_at_max_u = np.abs(shear[idx_max])

        E_s = 0.5 * F_at_max_u * u_max

        # Equivalent viscous damping ratio
        zeta_eq = 100 * E_d / (4.0 * np.pi * E_s) if E_s != 0 else 0.0

        # Plotting
        fig, ax = plt.subplots(figsize=(7, 6))

        # Hysteresis curve
        ax.plot(disp, shear, 'k-', linewidth=1.2, label="Hysteresis Loop")
        ax.scatter(disp[idx_max], shear[idx_max], color='blue', s=80,
                   zorder=5, label=r"$(u_{\rm max}, F_{\rm max})$")

        # Fill the enclosed area
        ax.fill(disp, shear, color='red', alpha=0.25, label=f"E$_d$ = {E_d:.3f} N·m")

        # Mark origin
        ax.axhline(0, color='gray', linestyle='--', linewidth=0.8)
        ax.axvline(0, color='gray', linestyle='--', linewidth=0.8)

        # Annotations and text box
        textstr = (f"$\\xi_{{\\rm eq}} = {zeta_eq:.4f}$ [%]\n"
                   f"$E_d = {E_d:.3f}$ N·mm\n"
                   f"$E_s = {E_s:.3f}$ N·mm")
        props = dict(boxstyle='round', facecolor='white', alpha=0.8)
        ax.text(0.05, 0.95, textstr, transform=ax.transAxes,
                fontsize=11, verticalalignment='top', bbox=props)

        ax.set_title(TITLE)
        ax.set_xlabel(XLABEL)
        ax.set_ylabel(YLABEL)
        ax.grid(True, linestyle='--', alpha=0.5)
        ax.legend(loc='lower right')
        fig.tight_layout()
        plt.show()

    return zeta_eq



