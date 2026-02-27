def T_STEEL_SECTION_FUN(secTag, matTag, Bf, tf, H, tw, 
                        nBf, nTf, nHw, nTw, 
                        DENSITY_STEEL, plot=True):
    """
    Fiber section definition for a solid T-shaped steel section
    and optional visualization of the section geometry.

    Parameters
    ----------
    secTag : int
        Section tag for OpenSees.
    matTag : int
        Material tag for steel.
    Bf : float
        Flange width (mm).
    tf : float
        Flange thickness (mm).
    H : float
        Total height of the section (from top of flange to bottom of web) (mm).
    tw : float
        Web thickness (mm).
    nBf : int
        Number of fiber divisions along the flange width.
    nTf : int
        Number of fiber divisions through the flange thickness.
    nHw : int
        Number of fiber divisions along the web depth.
    nTw : int
        Number of fiber divisions through the web thickness.
    DENSITY_STEEL : float
        Steel density (mass per volume) – used to compute mass per length.
    plot : bool, optional
        If True, draw a dimensioned cross‑section.

    Returns
    -------
    AREA : float
        Cross‑sectional area (mm²).
    ELE_MASS : float
        Mass per unit length (consistent with DENSITY_STEEL units).
    """
    import openseespy.opensees as ops
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.patches import Polygon, Rectangle
    from matplotlib.collections import PatchCollection

    # -------------------------
    # Geometric checks
    # -------------------------
    if tw >= Bf:
        raise ValueError("Web thickness tw must be smaller than flange width Bf")
    if tf >= H:
        raise ValueError("Flange thickness tf must be smaller than total height H")

    # -------------------------
    # Section Area Calculation
    # -------------------------
    # Flange area + web area
    AREA = Bf * tf + (H - tf) * tw   # [mm²]

    # -------------------------
    # OpenSees fiber section
    # -------------------------
    ops.section('Fiber', secTag)

    # Define coordinates (origin at centroid of the whole section)
    # Flange top fibre at y = 0, z = H/2   (if we place centroid at (0,0))
    # For simplicity we place the geometric centroid at (0,0).
    # First compute centroid location from top of flange.
    # A common choice: put the top of flange at z = H/2, so centroid is at z = 0.
    # That means the neutral axis is at mid‑height.  For a T‑beam this is not
    # the elastic centroid, but for fiber sections the position of the section
    # in local coordinates is arbitrary – the stiffness is unaffected.
    # Here we place the top of flange at z =  H/2 and the bottom of web at z = -H/2.
    # The centroid (0,0) then lies at mid‑height.  This is acceptable because
    # the section is defined in local coordinates and the user will later assign
    # the correct location via the element connectivity.

    # Flange boundaries
    y_flange_left  = -Bf / 2.0
    y_flange_right =  Bf / 2.0
    z_flange_top   =  H / 2.0
    z_flange_bot   =  H / 2.0 - tf

    # Web boundaries
    y_web_left   = -tw / 2.0
    y_web_right  =  tw / 2.0
    z_web_top    =  H / 2.0 - tf
    z_web_bot    = -H / 2.0

    # -------------------------
    # Flange (rectangular) – subdivided into nBf x nTf fibers
    # -------------------------
    # We use a loop over sub‑rectangles (patches).  OpenSees' 'quad' patch
    # requires four corner points.  Here we create a grid of small quads.
    dy_f = Bf / nBf
    dz_f = tf / nTf

    for i in range(nBf):
        y_left = y_flange_left + i * dy_f
        y_right = y_left + dy_f
        for j in range(nTf):
            z_bot = z_flange_bot + j * dz_f   # bottom of this sub‑rectangle
            z_top = z_bot + dz_f
            ops.patch('quad', matTag, 1, 1,   # 1x1 fibers inside this quad
                      y_left,  z_bot,
                      y_right, z_bot,
                      y_right, z_top,
                      y_left,  z_top)

    # -------------------------
    # Web (rectangular) – subdivided into nTw x nHw fibers
    # -------------------------
    dy_w = tw / nTw
    dz_w = (H - tf) / nHw

    for i in range(nTw):
        y_left = y_web_left + i * dy_w
        y_right = y_left + dy_w
        for j in range(nHw):
            z_bot = z_web_bot + j * dz_w
            z_top = z_bot + dz_w
            ops.patch('quad', matTag, 1, 1,
                      y_left,  z_bot,
                      y_right, z_bot,
                      y_right, z_top,
                      y_left,  z_top)

    # -------------------------
    # Plot section geometry
    # -------------------------
    if plot:
        fig, ax = plt.subplots(figsize=(8, 8))

        # Flange rectangle
        flange_rect = Rectangle(
            (y_flange_left, z_flange_bot),
            Bf, tf,
            linewidth=1, edgecolor='darkgray', facecolor='lightgray', alpha=0.7
        )
        ax.add_patch(flange_rect)

        # Web rectangle
        web_rect = Rectangle(
            (y_web_left, z_web_bot),
            tw, H - tf,
            linewidth=1, edgecolor='darkgray', facecolor='lightgray', alpha=0.7
        )
        ax.add_patch(web_rect)

        # Boundary lines
        # Outer flange
        ax.plot([y_flange_left, y_flange_right, y_flange_right, y_flange_left, y_flange_left],
                [z_flange_bot, z_flange_bot, z_flange_top, z_flange_top, z_flange_bot],
                'k-', linewidth=2, label='Flange')
        # Outer web
        ax.plot([y_web_left, y_web_right, y_web_right, y_web_left, y_web_left],
                [z_web_bot, z_web_bot, z_web_top, z_web_top, z_web_bot],
                'k-', linewidth=2, label='Web')

        # Centroid marker
        ax.plot(0, 0, 'r+', markersize=10, markeredgewidth=2, label='Centroid')

        # Center lines
        ax.axhline(y=0, color='gray', linestyle=':', linewidth=0.5, alpha=0.5)
        ax.axvline(x=0, color='gray', linestyle=':', linewidth=0.5, alpha=0.5)

        # Dimensions
        # Flange width
        ax.annotate('', xy=(y_flange_left, z_flange_top+5), xytext=(y_flange_right, z_flange_top+5),
                    arrowprops=dict(arrowstyle='<->', color='blue'))
        ax.text(0, z_flange_top+10, f'Bf = {Bf} mm', ha='center', va='bottom', color='blue')

        # Flange thickness
        ax.annotate('', xy=(y_flange_right+5, z_flange_bot), xytext=(y_flange_right+5, z_flange_top),
                    arrowprops=dict(arrowstyle='<->', color='green'))
        ax.text(y_flange_right+15, (z_flange_bot+z_flange_top)/2, f'tf = {tf} mm',
                ha='left', va='center', color='green', rotation=90)

        # Total height
        ax.annotate('', xy=(y_web_left-10, z_web_bot), xytext=(y_web_left-10, z_web_top+tf),
                    arrowprops=dict(arrowstyle='<->', color='blue'))
        ax.text(y_web_left-20, 0, f'H = {H} mm', ha='right', va='center', color='blue', rotation=90)

        # Web thickness
        ax.annotate('', xy=(y_web_left, z_web_bot-5), xytext=(y_web_right, z_web_bot-5),
                    arrowprops=dict(arrowstyle='<->', color='green'))
        ax.text(0, z_web_bot-15, f'tw = {tw} mm', ha='center', va='top', color='green')

        # Styling
        ax.set_aspect('equal')
        ax.set_xlabel('y (mm)')
        ax.set_ylabel('z (mm)')
        ax.set_title(f'T Steel Section\nBf={Bf}mm, tf={tf}mm, H={H}mm, tw={tw}mm')
        ax.grid(True, alpha=0.3)
        ax.legend(loc='upper right')

        # Margins
        margin = max(Bf, H) * 0.15
        ax.set_xlim([y_flange_left - margin, y_flange_right + margin])
        ax.set_ylim([z_web_bot - margin, z_flange_top + margin])

        # Section properties box
        props_text = f'Section Properties:\nArea = {AREA:.0f} mm²'
        ax.text(0.02, 0.98, props_text, transform=ax.transAxes,
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

        plt.tight_layout()
        plt.show()

    # Mass per unit length
    ELE_MASS = DENSITY_STEEL * AREA

    return AREA, ELE_MASS