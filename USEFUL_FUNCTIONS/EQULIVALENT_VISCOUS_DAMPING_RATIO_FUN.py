#%% EQULIVALENT VISCOUS DAMPING RATIO
def EQULIVALENT_VISCOUS_DAMPING_RATIO_FUN(displacement, base_shear, title="Hysteresis Curve"):
    """
    Compute the equivalent viscous damping ratio from a hysteresis loop.

    The dissipated energy per cycle (E_d) is the area enclosed by the loop,
    computed using the Shoelace formula. The maximum strain energy (E_s) is
    estimated as 0.5 * F_max * u_max, where u_max is the maximum absolute
    displacement and F_max is the absolute base shear at that same displacement.

    Parameters
    ----------
    displacement : array-like
        Displacement values (assumed to form a closed hysteresis loop).
    base_shear : array-like
        Corresponding base shear values.
    title : str, optional
        Title for the plot.

    Returns
    -------
    zeta_eq : float
        Equivalent viscous damping ratio.
    fig : matplotlib.figure.Figure
        Figure object with the plotted hysteresis loop and energy representation.
    
    THIS PYTHON SCRIPT WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
    """
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

    # Dissipated energy (E_d) via Shoelace formula ---
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
               f"$E_d = {E_d:.3f}$ N·m\n"
               f"$E_s = {E_s:.3f}$ N·m")
    props = dict(boxstyle='round', facecolor='white', alpha=0.8)
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes,
            fontsize=11, verticalalignment='top', bbox=props)

    ax.set_title(title)
    ax.set_xlabel("Displacement [m]")
    ax.set_ylabel("Base Shear [N]")
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.legend(loc='lower right')
    fig.tight_layout()

    return zeta_eq, fig

#%%---------------------------------
import numpy as np
theta = np.linspace(0, 2*np.pi, 500) 
u = 0.01 * np.sin(theta)                        # DISPLACEMENT
F = 1000 * np.sin(theta) + 200 * np.cos(theta)  # BASE-REACTION

zeta, fig = EQULIVALENT_VISCOUS_DAMPING_RATIO_FUN(u, F, title="Seismic Hysteresis")
fig.show()
print(f"Equivalent viscous damping ratio = {zeta:.4f}")

