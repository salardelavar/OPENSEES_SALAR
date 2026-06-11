def PLAN_IRREGULARITY_MASS_STIFFNESS_ANALYSIS_WITH_FLOOR_MASS_FUN(storey_height=3000.0):
    """
    Plan irregularity analysis for a one‑storey structure with mass and stiffness
    irregularities, using a stiffness‑transformation analogy.

    The method treats each vertical element (column or shear wall) as a "fibre"
    whose transformed property is its lateral stiffness.  The centre of rigidity
    (elastic centroid) is the stiffness‑weighted centre, and the torsional
    stiffness about that point is the polar moment of the stiffness distribution.

    Input:
        storey_height : height of columns/walls [mm] (default 3000)

    Internal data:
        - 12 columns
        - 4 shear walls
        - 5 additional point masses that have NO rigidity (floor‑only mass)

    Returns:
        x_CR, y_CR          : centre of rigidity coordinates [mm]
        x_CM, y_CM          : centre of mass coordinates [mm] (includes floor masses)
        Kx_tot, Ky_tot      : total translational stiffness in x and y [N/mm]
        K_tors              : torsional stiffness about CR [N·mm/rad]
        M_tot               : total mass [kg] (columns + walls + floor masses)

    THIS PYTHON SCRIPT WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    L = storey_height                    # column/wall height [mm]

    # ==================== INPUT DATA ====================
    # 12 columns: [x_cent, y_cent, width b, height h, E, mass]   (mm, N/mm², kg)
    columns = np.array([
        [   0,    0, 300, 300, 25000, 800],   # C1  corner
        [6000,    0, 300, 300, 25000, 1000],  # C2
        [12000,   0, 300, 300, 25000, 1100],  # C3
        [12000,4000, 300, 300, 25000, 900],   # C4
        [12000,8000, 300, 300, 25000, 145],   # C5
        [6000, 8000, 300, 300, 25000, 920],   # C6
        [   0, 8000, 300, 300, 25000, 850],   # C7
        [   0, 4000, 300, 300, 25000, 800],   # C8
        [3000, 2000, 250, 250, 30000,1100],   # C9  interior heavier + stiffer
        [9000, 2000, 250, 250, 30000,1200],   # C10
        [3000, 6000, 250, 250, 30000,1800],   # C11
        [9000, 6000, 250, 250, 30000,1200],   # C12
    ])

    # 4 shear walls: [x_cent, y_cent, length, thickness, E, mass, orientation]
    shear_walls = np.array([
        [3000, 8000, 6000, 300, 28000, 2000, 0],         # W1: along x, resists y
        [9000, 0.0, 6000, 400, 28000, 2000, 0],          # W2: along x, resists y
        [0.0, 6000, 4000.0, 200, 28000, 2000, 1],        # W3: along y, resists x
        [12000.0, 2000.0, 4000.0, 200, 28000, 2000, 1],  # W4: along y, resists x
    ])

    # ---- 5 ADDITIONAL FLOOR MASSES (NO STIFFNESS) ----
    # Format: [x_cent, y_cent, mass]   (mm, kg)
    floor_masses = np.array([
        [6000, 4000, 1500],   # F1 – centre of floor
        [3000, 4000, 1200],   # F2 – left interior
        [9000, 4000, 1200],   # F3 – right interior
        [6000, 2000,  800],   # F4 – lower interior
        [6000, 6000,  800],   # F5 – upper interior
    ])

    # ==================== STIFFNESS CALCULATION ====================
    # For a fixed‑fixed column/wall: k = 12 E I / L^3

    # Columns
    col_x = columns[:,0]; col_y = columns[:,1]
    col_b = columns[:,2]; col_h = columns[:,3]
    col_E = columns[:,4]; col_mass = columns[:,5]
    Ix_col = (col_b * col_h**3) / 12   # bending about x‑axis -> resists y
    Iy_col = (col_h * col_b**3) / 12   # bending about y‑axis -> resists x
    Kx_col = 12 * col_E * Iy_col / L**3
    Ky_col = 12 * col_E * Ix_col / L**3

    # Shear walls
    wall_x = shear_walls[:,0]; wall_y = shear_walls[:,1]
    wall_L = shear_walls[:,2]; wall_t = shear_walls[:,3]
    wall_E = shear_walls[:,4]; wall_mass = shear_walls[:,5]
    wall_orient = shear_walls[:,6].astype(int)

    Kx_wall = np.zeros(len(shear_walls))
    Ky_wall = np.zeros(len(shear_walls))
    for i in range(len(shear_walls)):
        Lw = wall_L[i]; t = wall_t[i]; Ew = wall_E[i]
        if wall_orient[i] == 0:          # length along x, resists y
            I = (t * Lw**3) / 12
            Ky_wall[i] = 12 * Ew * I / L**3
        else:                            # length along y, resists x
            I = (t * Lw**3) / 12
            Kx_wall[i] = 12 * Ew * I / L**3

    # Combine structural stiffness (floor masses do NOT contribute to stiffness)
    Kx_all = np.concatenate([Kx_col, Kx_wall])
    Ky_all = np.concatenate([Ky_col, Ky_wall])
    x_all  = np.concatenate([col_x, wall_x])
    y_all  = np.concatenate([col_y, wall_y])
    mass_all = np.concatenate([col_mass, wall_mass])

    # ==================== CENTRE OF RIGIDITY (CR) ====================
    Kx_tot = np.sum(Kx_all)
    Ky_tot = np.sum(Ky_all)
    y_CR = np.sum(Kx_all * y_all) / Kx_tot
    x_CR = np.sum(Ky_all * x_all) / Ky_tot

    # Torsional stiffness about CR
    K_tors = np.sum(Kx_all * (y_all - y_CR)**2) + np.sum(Ky_all * (x_all - x_CR)**2)

    # ==================== CENTRE OF MASS (CM) – now includes floor masses ====================
    # Extract floor mass coordinates
    fm_x = floor_masses[:,0]
    fm_y = floor_masses[:,1]
    fm_mass = floor_masses[:,2]

    # Combine all masses (structural + floor) for CM calculation
    total_mass_all = np.concatenate([mass_all, fm_mass])
    total_x_all    = np.concatenate([x_all, fm_x])
    total_y_all    = np.concatenate([y_all, fm_y])

    M_tot = np.sum(total_mass_all)
    x_CM = np.sum(total_mass_all * total_x_all) / M_tot
    y_CM = np.sum(total_mass_all * total_y_all) / M_tot

    e_x = x_CM - x_CR
    e_y = y_CM - y_CR

    # ==================== OUTPUT ====================
    print("\n" + "="*60)
    print("       STRUCTURAL IRREGULARITY ANALYSIS")
    print("       (12 columns + 4 shear walls + 5 floor masses)")
    print("="*60)
    print(f"Storey height                = {L:.0f} mm")
    print(f"Total mass (including floor) = {M_tot:.1f} kg")
    print(f"Centre of Mass (CM)          = ({x_CM:.1f}, {y_CM:.1f}) mm")
    print(f"Centre of Rigidity (CR)      = ({x_CR:.1f}, {y_CR:.1f}) mm")
    print(f"Eccentricity e_x (CM-CR)     = {e_x:.1f} mm")
    print(f"Eccentricity e_y (CM-CR)     = {e_y:.1f} mm")
    print(f"Total stiffness Kx           = {Kx_tot:.2e} N/mm")
    print(f"Total stiffness Ky           = {Ky_tot:.2e} N/mm")
    print(f"Torsional stiffness K_tors   = {K_tors:.2e} N·mm/rad")
    print(f"Torsional radius (√(Kt/Kx))  = {np.sqrt(K_tors/Kx_tot):.1f} mm  (x‑dir)")
    print(f"Torsional radius (√(Kt/Ky))  = {np.sqrt(K_tors/Ky_tot):.1f} mm  (y‑dir)")
    print("="*60)

    # ==================== PLOT WITH TAGS ====================
    fig, ax = plt.subplots(figsize=(10, 6))

    # Draw columns and tag them C1..C12
    for i in range(len(columns)):
        xc, yc, b, h = columns[i,0], columns[i,1], columns[i,2], columns[i,3]
        rect = Rectangle((xc - b/2, yc - h/2), b, h,
                         edgecolor='black', facecolor='lightblue', alpha=0.7)
        ax.add_patch(rect)
        # Tag column number
        ax.text(xc, yc, f'C{i+1}', ha='center', va='center', fontsize=8,
                fontweight='bold', color='black')

    # Draw shear walls and tag them W1..W4
    for i in range(len(shear_walls)):
        xw, yw, Lw, tw, orient = (shear_walls[i,0], shear_walls[i,1],
                                   shear_walls[i,2], shear_walls[i,3],
                                   int(shear_walls[i,6]))
        if orient == 0:          # length along x
            rect = Rectangle((xw - Lw/2, yw - tw/2), Lw, tw,
                             edgecolor='red', facecolor='salmon', alpha=0.6, linewidth=2)
        else:                    # length along y
            rect = Rectangle((xw - tw/2, yw - Lw/2), tw, Lw,
                             edgecolor='red', facecolor='salmon', alpha=0.6, linewidth=2)
        ax.add_patch(rect)
        ax.text(xw, yw, f'W{i+1}', ha='center', va='center', fontsize=9,
                fontweight='bold', color='darkred')

    # Draw the 5 floor‑mass points (no rigidity) as grey squares with labels
    for i in range(len(floor_masses)):
        xf, yf, mf = floor_masses[i]
        # Draw a small square to represent the mass (size not structural)
        side = 300  # visual size, not related to real dimensions
        rect = Rectangle((xf - side/2, yf - side/2), side, side,
                         edgecolor='grey', facecolor='lightgrey', alpha=0.5,
                         linestyle='--')
        ax.add_patch(rect)
        ax.text(xf, yf, f'F{i+1}\n({mf} kg)', ha='center', va='center', fontsize=7,
                fontweight='bold', color='dimgrey')

    # Mark centroids
    ax.plot(x_CR, y_CR, 'r+', markersize=14, markeredgewidth=2, label='Centre of Rigidity (CR)')
    ax.plot(x_CM, y_CM, 'gx', markersize=12, markeredgewidth=2, label='Centre of Mass (CM)')

    ax.set_aspect('equal')
    ax.set_xlim(-1000, 13000)
    ax.set_ylim(-1000, 9000)
    ax.set_title('Structural Plan: 12 Columns + 4 Shear Walls + 5 Floor Masses\n(irregular mass & stiffness)')
    ax.set_xlabel('x [mm]')
    ax.set_ylabel('y [mm]')
    ax.legend()
    ax.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.show()

    return x_CR, y_CR, x_CM, y_CM, Kx_tot, Ky_tot, K_tors, M_tot

PLAN_IRREGULARITY_MASS_STIFFNESS_ANALYSIS_WITH_FLOOR_MASS_FUN(storey_height=3000.0)