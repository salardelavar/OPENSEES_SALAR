def CONCRETE_RECTANGULAR_HOLLOW_SECTION_FUN_3D(
    SEC_TAG, STEEL_TYPE, fc, Kfc,
    b_o, h_o, b_i, h_i,
    nFib, CONCRETE_DENSITY,
    nRebars=68, rebar_dia=25,
    plot=True
    ):
    """
    Create a hollow rectangular confined‑concrete fiber section for OpenSees,
    with automatic placement of reinforcing bars along the mid‑thickness line.

    Parameters
    ----------
    secTag          : int   – section identifier
    STEEL_TYPE      : str   – 'ELASTIC' or 'INELASTIC'
    fc              : float – unconfined concrete compressive strength (MPa, negative)
    Kfc             : float – ratio of confined to unconfined concrete strength
    b_o, h_o        : float – outer width and height (mm)
    b_i, h_i        : float – inner width and height (mm)
    nFib            : int   – number of fibers per side for concrete patches
    CONCRETE_DENSITY: float – density of concrete (kg/m³) – converted to N·s²/mm³ internally
    plot            : bool  – if True, draw a 2D sketch of the section
    nRebars         : int   – total number of reinforcing bars (default = 68)
    rebar_dia       : float – diameter of each rebar (mm), uniform for all bars (default = 25)

    Returns
    -------
    h_o             : float – section height (mm)
    ELE_MASS        : float – element mass per unit length (kg/mm)

    THIS PYTHON SCRIPT IS WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI) 
    """

    import openseespy.opensees as ops
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    import numpy as np

    # ------------------------------------------------------------------
    # 1. MATERIAL DEFINITIONS (units: mm, N)
    # ------------------------------------------------------------------

    # --- Unconfined cover concrete ---
    fcU = -fc                         # compressive strength (positive)
    Ec = 4700 * np.sqrt(-fcU)         # elastic modulus (N/mm²)
    ec0U = 2 * fcU / Ec               # strain at peak stress
    fcUU = 0.2 * fcU                  # ultimate compressive strength
    ecuU = 5 * ec0U                   # ultimate compressive strain

    # --- Confined core concrete ---
    fcC = Kfc * fcU                   # confined compressive strength
    EcC = 4700 * np.sqrt(-fcC)        # elastic modulus for confined concrete
    ec0C = 2 * fcC / EcC              # strain at peak stress
    fcUC = 0.65 * fcC                 # ultimate compressive strength
    ecuC = 15 * ec0C                  # ultimate compressive strain
    Lambda = 0.1                      # unloading slope ratio

    # Tensile properties
    ftC = 0.7 * np.sqrt(-fcC)         # tensile strength (confined)
    ftU = 0.7 * np.sqrt(-fcU)         # tensile strength (unconfined)
    EtsC = ftC / np.abs(ec0C)         # tension softening stiffness (confined)
    EtsU = ftU / np.abs(ec0U)         # tension softening stiffness (unconfined)

    # --- Reinforcing steel ---
    fy = 400.0                        # yield strength (N/mm²)
    Es = 2.0e5                        # elastic modulus (N/mm²)
    ey = fy / Es                      # yield strain
    fu = 1.1818 * fy                  # ultimate strength
    esu = 0.09                        # ultimate strain
    Esh = (fu - fy) / (esu - ey)      # post‑yield hardening modulus
    Bs = Esh / Es

    # ------------------------------------------------------------------
    # 2. CREATE UNIAXIAL MATERIALS IN OPENSEES
    # ------------------------------------------------------------------
    coreTag, coverTag, steelTag = secTag + 100, secTag + 200, secTag + 300

    if STEEL_TYPE == 'ELASTIC':
        ops.uniaxialMaterial('Steel01', steelTag, fy, Es, Bs)
    elif STEEL_TYPE == 'INELASTIC':
        pinchX = 0.8   # pinching factor in X
        pinchY = 0.5   # pinching factor in Y
        damage1 = 0.0  # ductility‑based damage
        damage2 = 0.0  # energy‑based damage
        beta = 0.1     # stiffness degradation
        ops.uniaxialMaterial(
            'Hysteretic', steelTag,
            fy, ey, fu, esu, 0.2*fu, 1.1*esu,
            -fy, -ey, -fu, -esu, -0.2*fu, -1.1*esu,
            pinchX, pinchY, damage1, damage2, beta
        )

    # Concrete02 materials (includes tension softening)
    ops.uniaxialMaterial('Concrete02', coreTag, fcC, ec0C, fcUC, ecuC,
                         Lambda, ftC, EtsC)
    ops.uniaxialMaterial('Concrete02', coverTag, fcU, ec0U, fcUU, ecuU,
                         Lambda, ftU, EtsU)

    # ------------------------------------------------------------------
    # 3. SECTION AREA
    # ------------------------------------------------------------------
    AREA = (b_o * h_o) - (b_i * h_i)
    print(f"Total Section Area = {AREA:.2f} mm²")

    # ------------------------------------------------------------------
    # 4. FIBER SECTION – PATCHES FOR COVER AND CORE
    # ------------------------------------------------------------------
    ops.section('Fiber', SEC_TAG, '-GJ', 1.0e6)

    def add_rect(matTag, y_bot, y_top, x_left, x_right):
        """Add a rectangular patch of fibers."""
        ops.patch('rect', matTag, nFib, nFib,
                  x_left, y_bot, x_right, y_top)

    T_h = h_o - h_i   # wall thickness in height direction
    T_b = b_o - b_i   # wall thickness in width direction

    # Outer cover – top and bottom strips
    add_rect(coverTag,  h_o/2 - T_h,  h_o/2, -b_o/2,  b_o/2)
    add_rect(coverTag, -h_o/2 + T_h, -h_o/2, -b_o/2,  b_o/2)
    # Left and right strips
    add_rect(coverTag, -h_o/2 + T_h,  h_o/2 - T_h, -b_o/2, -b_o/2 + T_b)
    add_rect(coverTag, -h_o/2 + T_h,  h_o/2 - T_h,  b_o/2,  b_o/2 - T_b)

    # ------------------------------------------------------------------
    # 5. AUTOMATIC REBAR PLACEMENT ALONG THE MID‑THICKNESS LINE
    # ------------------------------------------------------------------
    # Mid‑line dimensions
    b_mid = (b_o + b_i) / 2.0
    h_mid = (h_o + h_i) / 2.0
    half_b = b_mid / 2.0
    half_h = h_mid / 2.0
    perimeter = 2.0 * (b_mid + h_mid)

    rebars_coords = []
    for i in range(nRebars):
        s = (i / nRebars) * perimeter
        # Traverse the rectangle counter‑clockwise: bottom → right → top → left
        if s <= b_mid:                     # bottom edge
            x = -half_b + s
            y = -half_h
        elif s <= b_mid + h_mid:           # right edge
            x = half_b
            y = -half_h + (s - b_mid)
        elif s <= 2 * b_mid + h_mid:       # top edge
            x = half_b - (s - b_mid - h_mid)
            y = half_h
        else:                              # left edge
            x = -half_b
            y = half_h - (s - 2 * b_mid - h_mid)
        rebars_coords.append((x, y))

    # Add each rebar as a fiber
    for (x, y) in rebars_coords:
        area_bar = np.pi * (rebar_dia ** 2) / 4.0
        ops.fiber(x, y, area_bar, steelTag)

    # ------------------------------------------------------------------
    # 6. PLOT (if requested)
    # ------------------------------------------------------------------
    if plot:
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.set_xlabel('Width (mm)')
        ax.set_ylabel('Height (mm)')
        ax.set_title('Hollow rectangular section with rebars')
        ax.grid(True, linestyle='--', alpha=0.4)

        # Outer concrete (cover)
        outer = patches.Rectangle(
            (-b_o/2, -h_o/2), b_o, h_o,
            linewidth=1.5, edgecolor='black', facecolor='lightgray'
        )
        ax.add_patch(outer)

        # Inner hollow hole (white)
        inner = patches.Rectangle(
            (-b_i/2, -h_i/2), b_i, h_i,
            linewidth=1.5, edgecolor='black', facecolor='white'
        )
        ax.add_patch(inner)

        # Rebars with numbering
        for idx, (x, y) in enumerate(rebars_coords, start=1):
            circ = patches.Circle(
                (x, y), rebar_dia/2,
                edgecolor='red', facecolor='red'
            )
            ax.add_patch(circ)
            ax.text(
                x, y + rebar_dia/2 + 2, f'{idx}',
                color='purple', fontsize=6,
                ha='center', va='bottom', fontweight='bold'
            )

        # Axis limits
        lim = max(b_o, h_o) * 0.6 + 20
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_aspect('equal')
        plt.show()

    # ------------------------------------------------------------------
    # 7. MASS PER UNIT LENGTH
    # ------------------------------------------------------------------
    # Convert density from kg/m³ to kg/mm (1 m³ = 1e9 mm³)
    density_kg_per_mm3 = CONCRETE_DENSITY / 1e9   # kg/mm³
    ELE_MASS = density_kg_per_mm3 * AREA          # kg/mm

    return h_o, ELE_MASS