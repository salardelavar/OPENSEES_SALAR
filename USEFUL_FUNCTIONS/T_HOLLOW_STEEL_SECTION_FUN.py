def T_HOLLOW_STEEL_SECTION_FUN(secTag, matTag, B, H, tf, tw, t, nFib, DENSITY_STEEL, plot=True):
    """
    Fiber section definition for a hollow T‑section (constant wall thickness)
    and optional visualization of the section geometry.

    Parameters
    ----------
    secTag : int
        Section tag for OpenSees
    matTag : int
        Material tag for steel
    B : float
        Flange width (mm)
    H : float
        Total depth (mm) – from top of flange to bottom of web
    tf : float
        Flange thickness (mm)
    tw : float
        Web thickness (mm)
    t : float
        Wall thickness (mm) – must be less than half of flange/web dimensions
    nFib : int
        Number of fibers for discretisation (used in both directions)
    DENSITY_STEEL : float
        Steel density (mass per volume) to compute mass per length
    plot : bool, optional
        If True, draw the cross‑section (default True)

    Returns
    -------
    AREA : float
        Steel cross‑sectional area (mm²)
    ELE_MASS : float
        Mass per unit length (consistent units)
    """
    import openseespy.opensees as ops
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.patches import Polygon, Rectangle
    from matplotlib.collections import PatchCollection

    # -------------------------------------------------------------------------
    # Geometric checks
    # -------------------------------------------------------------------------
    if t >= B/2 or t >= tw/2:
        raise ValueError("Wall thickness t must be smaller than B/2 and tw/2")
    if t >= tf:
        raise ValueError("Wall thickness t must be smaller than flange thickness tf")
    if tf >= H:
        raise ValueError("Flange thickness tf must be smaller than total depth H")
    if tw >= B:
        raise ValueError("Web thickness tw must be smaller than flange width B")

    # -------------------------------------------------------------------------
    # Section area calculation (outer T minus inner T)
    # -------------------------------------------------------------------------
    # Outer T area
    A_outer_flange = B * tf
    A_outer_web   = tw * (H - tf)
    A_outer = A_outer_flange + A_outer_web

    # Inner T area (hollow region)
    A_inner_flange = (B - 2*t) * (tf - t)
    A_inner_web    = (tw - 2*t) * (H - tf - t)
    A_inner = A_inner_flange + A_inner_web

    AREA = A_outer - A_inner   # mm²

    # -------------------------------------------------------------------------
    # Coordinate definitions (centroid at (0,0))
    # -------------------------------------------------------------------------
    # Outer T
    # Flange
    y_flange_left   = -B/2
    y_flange_right  =  B/2
    z_flange_bottom =  H/2 - tf
    z_flange_top    =  H/2

    # Web
    y_web_left   = -tw/2
    y_web_right  =  tw/2
    z_web_bottom = -H/2
    z_web_top    =  H/2 - tf

    # Inner T (hollow)
    # Flange hollow
    y_flange_hollow_left   = -B/2 + t
    y_flange_hollow_right  =  B/2 - t
    z_flange_hollow_bottom =  H/2 - tf + t
    z_flange_hollow_top    =  H/2 - t

    # Web hollow
    y_web_hollow_left   = -tw/2 + t
    y_web_hollow_right  =  tw/2 - t
    z_web_hollow_bottom = -H/2 + t
    z_web_hollow_top    =  H/2 - tf - t

    # -------------------------------------------------------------------------
    # OpenSees fiber section
    # -------------------------------------------------------------------------
    ops.section('Fiber', secTag)

    # ---- Flange steel patches (8 rectangles) ----
    # 1. Left vertical strip
    ops.patch('quad', matTag, nFib, 1,
              y_flange_left,                z_flange_bottom,
              y_flange_left + t,             z_flange_bottom,
              y_flange_left + t,             z_flange_top,
              y_flange_left,                z_flange_top)

    # 2. Right vertical strip
    ops.patch('quad', matTag, nFib, 1,
              y_flange_right - t,            z_flange_bottom,
              y_flange_right,                z_flange_bottom,
              y_flange_right,                z_flange_top,
              y_flange_right - t,            z_flange_top)

    # 3. Top horizontal strip (above hollow)
    ops.patch('quad', matTag, nFib, 1,
              y_flange_hollow_left,          z_flange_hollow_top,
              y_flange_hollow_right,         z_flange_hollow_top,
              y_flange_hollow_right,         z_flange_top,
              y_flange_hollow_left,          z_flange_top)

    # 4. Bottom horizontal strip (below hollow)
    ops.patch('quad', matTag, nFib, 1,
              y_flange_hollow_left,          z_flange_bottom,
              y_flange_hollow_right,         z_flange_bottom,
              y_flange_hollow_right,         z_flange_hollow_bottom,
              y_flange_hollow_left,          z_flange_hollow_bottom)

    # ---- Web steel patches (also 4 rectangles) ----
    # 5. Left vertical strip
    ops.patch('quad', matTag, nFib, 1,
              y_web_left,                    z_web_bottom,
              y_web_left + t,                 z_web_bottom,
              y_web_left + t,                 z_web_top,
              y_web_left,                    z_web_top)

    # 6. Right vertical strip
    ops.patch('quad', matTag, nFib, 1,
              y_web_right - t,                z_web_bottom,
              y_web_right,                    z_web_bottom,
              y_web_right,                    z_web_top,
              y_web_right - t,                z_web_top)

    # 7. Bottom horizontal strip (below hollow)
    ops.patch('quad', matTag, nFib, 1,
              y_web_hollow_left,              z_web_bottom,
              y_web_hollow_right,             z_web_bottom,
              y_web_hollow_right,             z_web_hollow_bottom,
              y_web_hollow_left,              z_web_hollow_bottom)

    # 8. Top horizontal strip (above hollow, below flange)
    ops.patch('quad', matTag, nFib, 1,
              y_web_hollow_left,              z_web_hollow_top,
              y_web_hollow_right,             z_web_hollow_top,
              y_web_hollow_right,             z_web_top,
              y_web_hollow_left,              z_web_top)

    # -------------------------------------------------------------------------
    # Plot section geometry (optional)
    # -------------------------------------------------------------------------
    if plot:
        fig, ax = plt.subplots(figsize=(8, 8))

        # Collect all steel patches for visualization
        patches = []

        # Flange patches (same coordinates as above)
        patches.append(Polygon([
            [y_flange_left, z_flange_bottom],
            [y_flange_left + t, z_flange_bottom],
            [y_flange_left + t, z_flange_top],
            [y_flange_left, z_flange_top]
        ], closed=True))
        patches.append(Polygon([
            [y_flange_right - t, z_flange_bottom],
            [y_flange_right, z_flange_bottom],
            [y_flange_right, z_flange_top],
            [y_flange_right - t, z_flange_top]
        ], closed=True))
        patches.append(Polygon([
            [y_flange_hollow_left, z_flange_hollow_top],
            [y_flange_hollow_right, z_flange_hollow_top],
            [y_flange_hollow_right, z_flange_top],
            [y_flange_hollow_left, z_flange_top]
        ], closed=True))
        patches.append(Polygon([
            [y_flange_hollow_left, z_flange_bottom],
            [y_flange_hollow_right, z_flange_bottom],
            [y_flange_hollow_right, z_flange_hollow_bottom],
            [y_flange_hollow_left, z_flange_hollow_bottom]
        ], closed=True))

        # Web patches
        patches.append(Polygon([
            [y_web_left, z_web_bottom],
            [y_web_left + t, z_web_bottom],
            [y_web_left + t, z_web_top],
            [y_web_left, z_web_top]
        ], closed=True))
        patches.append(Polygon([
            [y_web_right - t, z_web_bottom],
            [y_web_right, z_web_bottom],
            [y_web_right, z_web_top],
            [y_web_right - t, z_web_top]
        ], closed=True))
        patches.append(Polygon([
            [y_web_hollow_left, z_web_bottom],
            [y_web_hollow_right, z_web_bottom],
            [y_web_hollow_right, z_web_hollow_bottom],
            [y_web_hollow_left, z_web_hollow_bottom]
        ], closed=True))
        patches.append(Polygon([
            [y_web_hollow_left, z_web_hollow_top],
            [y_web_hollow_right, z_web_hollow_top],
            [y_web_hollow_right, z_web_top],
            [y_web_hollow_left, z_web_top]
        ], closed=True))

        # Add steel patches to plot
        collection = PatchCollection(patches, facecolor='lightgray',
                                     edgecolor='darkgray', linewidth=1, alpha=0.7)
        ax.add_collection(collection)

        # Draw outer boundary (T shape)
        # Flange top edge
        ax.plot([y_flange_left, y_flange_right], [z_flange_top, z_flange_top], 'k-', linewidth=2)
        # Flange bottom edge (left part, then web top, then right part)
        ax.plot([y_flange_left, y_web_left], [z_flange_bottom, z_flange_bottom], 'k-', linewidth=2)
        ax.plot([y_web_left, y_web_right], [z_web_top, z_web_top], 'k-', linewidth=2)  # web top (hidden by flange)
        ax.plot([y_web_right, y_flange_right], [z_flange_bottom, z_flange_bottom], 'k-', linewidth=2)
        # Web bottom edge
        ax.plot([y_web_left, y_web_right], [z_web_bottom, z_web_bottom], 'k-', linewidth=2)
        # Left vertical edges
        ax.plot([y_flange_left, y_flange_left], [z_flange_bottom, z_flange_top], 'k-', linewidth=2)
        ax.plot([y_web_left, y_web_left], [z_web_bottom, z_web_top], 'k-', linewidth=2)
        # Right vertical edges
        ax.plot([y_flange_right, y_flange_right], [z_flange_bottom, z_flange_top], 'k-', linewidth=2)
        ax.plot([y_web_right, y_web_right], [z_web_bottom, z_web_top], 'k-', linewidth=2)

        # Draw inner hollow boundary (T shape, dashed)
        # Flange hollow top edge
        ax.plot([y_flange_hollow_left, y_flange_hollow_right],
                [z_flange_hollow_top, z_flange_hollow_top], 'r--', linewidth=2)
        # Flange hollow bottom edge
        ax.plot([y_flange_hollow_left, y_flange_hollow_right],
                [z_flange_hollow_bottom, z_flange_hollow_bottom], 'r--', linewidth=2)
        # Web hollow top edge
        ax.plot([y_web_hollow_left, y_web_hollow_right],
                [z_web_hollow_top, z_web_hollow_top], 'r--', linewidth=2)
        # Web hollow bottom edge
        ax.plot([y_web_hollow_left, y_web_hollow_right],
                [z_web_hollow_bottom, z_web_hollow_bottom], 'r--', linewidth=2)
        # Left vertical edges of hollow
        ax.plot([y_flange_hollow_left, y_flange_hollow_left],
                [z_flange_hollow_bottom, z_flange_hollow_top], 'r--', linewidth=2)
        ax.plot([y_web_hollow_left, y_web_hollow_left],
                [z_web_hollow_bottom, z_web_hollow_top], 'r--', linewidth=2)
        # Right vertical edges of hollow
        ax.plot([y_flange_hollow_right, y_flange_hollow_right],
                [z_flange_hollow_bottom, z_flange_hollow_top], 'r--', linewidth=2)
        ax.plot([y_web_hollow_right, y_web_hollow_right],
                [z_web_hollow_bottom, z_web_hollow_top], 'r--', linewidth=2)

        # Centroid marker
        ax.plot(0, 0, 'r+', markersize=10, markeredgewidth=2, label='Centroid')

        # Axis settings
        ax.set_aspect('equal')
        ax.set_xlabel('y (mm)')
        ax.set_ylabel('z (mm)')
        ax.set_title(f'Hollow T‑Section\nB={B}mm, H={H}mm, tf={tf}mm, tw={tw}mm, t={t}mm')
        ax.grid(True, alpha=0.3)
        margin = max(B, H) * 0.15
        ax.set_xlim([y_flange_left - margin, y_flange_right + margin])
        ax.set_ylim([z_web_bottom - margin, z_flange_top + margin])

        # Text box with area
        props_text = f'Steel area = {AREA:.0f} mm²'
        ax.text(0.02, 0.98, props_text, transform=ax.transAxes,
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

        plt.tight_layout()
        plt.show()

    # -------------------------------------------------------------------------
    # Mass per unit length
    # -------------------------------------------------------------------------
    ELE_MASS = DENSITY_STEEL * AREA   # (mass/length)

    return AREA, ELE_MASS