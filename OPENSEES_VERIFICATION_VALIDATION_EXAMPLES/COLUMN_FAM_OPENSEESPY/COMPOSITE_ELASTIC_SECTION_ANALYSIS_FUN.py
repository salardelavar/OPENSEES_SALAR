def COMPOSITE_ELASTIC_SECTION_ANALYSIS_FUN():
    """
    Composite elastic section analysis using the transformed section method.
    
    Input format for each rectangle:
        [x_centroid, y_centroid, width (b), height (h), modulus_of_elasticity (E)]
    
    Returns:
        x_el, y_el        : elastic centroid coordinates (transformed section)
        A_trans           : total transformed area
        Ix_trans, Iy_trans: transformed moments of inertia about the elastic centroid
        E_ref             : reference modulus used for transformation
        
    THIS PYTHON SCRIPT WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)    
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    # ==================== INPUT DATA ====================
    # x_centroid [mm], y_centroid [mm], width [mm], height [mm], Modulus of Elasticity [N/mm²]
    rects = np.array([
        [-55, 250, 10, 580, 200000.0],   # Left Web Plate
        [250, 520, 500, 20, 500000.0],   # Top Plate
        [250, 500, 600, 20, 200000.0],   # Top flange
        [450, 250, 20, 480, 300000.0],   # Right web
        [250, 250, 20, 480, 400000.0],   # middle web
        [50, 250, 20, 480, 200000.0],    # Left web
        [250, 0, 600, 20, 200000.0],     # Bottom flange
        [250, -20, 500, 20, 100000.0],   # Bottom Plate
        [555, 250, 10, 580, 300000.0],   # Right Web Plate
    ])

    # Extract data
    x_cent = rects[:, 0]
    y_cent = rects[:, 1]
    widths = rects[:, 2]
    heights = rects[:, 3]
    E_mod = rects[:, 4]

    # Reference modulus (largest Module of Elasticity, gives transformed areas ≤ actual areas)
    E_ref = np.max(E_mod)
    n = E_mod / E_ref          # modular ratios

    # Actual geometric areas (for reference)
    A_geom = widths * heights

    # Transformed areas (account for different Modulus of Elasticity)
    A_trans = A_geom * n

    # ============== ELASTIC CENTROID ==============
    A_trans_total = np.sum(A_trans)
    x_el = np.sum(A_trans * x_cent) / A_trans_total
    y_el = np.sum(A_trans * y_cent) / A_trans_total

    # ============== TRANSFORMED MOMENTS OF INERTIA ==============
    Ix_trans_total = 0.0
    Iy_trans_total = 0.0

    for i in range(len(rects)):
        b = widths[i]
        h = heights[i]
        A_t = A_trans[i]
        # Transformed local moments of inertia (about part's own centroid)
        Ix_loc = (b * h**3) / 12 * n[i]
        Iy_loc = (h * b**3) / 12 * n[i]

        # Parallel axis contribution to global elastic centroid
        Ix_trans_total += Ix_loc + A_t * (y_cent[i] - y_el)**2
        Iy_trans_total += Iy_loc + A_t * (x_cent[i] - x_el)**2

    # ============== OUTPUT ==============
    print("\n" + "="*50)
    print("       COMPOSITE ELASTIC SECTION ANALYSIS")
    print("="*50)
    print(f"Reference modulus (E_ref) = {E_ref:.0f} N/mm²")
    print(f"Elastic centroid: x̄ = {x_el:.2f} mm,  ȳ = {y_el:.2f} mm")
    print(f"Transformed area: A_trans = {A_trans_total:.2f} mm²")
    print(f"Transformed Ix (about elastic centroid) = {Ix_trans_total:.2e} mm⁴")
    print(f"Transformed Iy (about elastic centroid) = {Iy_trans_total:.2e} mm⁴")
    print("="*50)

    # Optional: geometric centroid and I (if all materials were the same)
    A_geom_total = np.sum(A_geom)
    x_g = np.sum(A_geom * x_cent) / A_geom_total
    y_g = np.sum(A_geom * y_cent) / A_geom_total
    print(f"\nGeometric centroid (if monolithic): x̄ = {x_g:.2f}, ȳ = {y_g:.2f}")
    print(f"Geometric area: {A_geom_total:.2f} mm²")

    # ============== PLOT ==============
    fig, ax = plt.subplots(figsize=(8, 6))

    # Normalize E for coloring (darker = stiffer)
    norm_E = (E_mod - E_mod.min()) / (E_mod.max() - E_mod.min() + 1e-12)
    cmap = plt.cm.Blues

    for i, (x, y, b, h, E_val) in enumerate(rects):
        # Rectangle defined by bottom‑left corner
        rect = Rectangle((x - b/2, y - h/2), b, h,
                         edgecolor='black',
                         facecolor=cmap(norm_E[i]),
                         alpha=0.7,
                         linewidth=1.5)
        ax.add_patch(rect)
        # Optional: annotate Modulus of Elasticity value inside large rectangles
        if b * h > 1000:   # only for sufficiently large parts
            ax.text(x, y, f'{E_val/1000:.0f}GPa', ha='center', va='center', fontsize=8)

    # Mark elastic centroid
    ax.plot(x_el, y_el, 'r+', markersize=12, markeredgewidth=2, label='Elastic centroid')
    # Mark geometric centroid (for comparison)
    ax.plot(x_g, y_g, 'gx', markersize=10, markeredgewidth=2, label='Geometric centroid')

    ax.set_aspect('equal')
    ax.set_title('Composite Section (color = relative stiffness)', fontsize=12)
    ax.set_xlabel('x (mm)')
    ax.set_ylabel('y (mm)')
    ax.legend()
    ax.grid(True, linestyle='--', alpha=0.5)

    plt.tight_layout()
    plt.show()

    return x_el, y_el, A_trans_total, Ix_trans_total, Iy_trans_total, E_ref

#COMPOSITE_ELASTIC_SECTION_ANALYSIS_FUN()