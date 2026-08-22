def CONCRETE_CIRCULAR_HOLLOW_PIPE_SECTION_FUN(
    secTag, STEEL_TYPE, Do, t,
    numSubdivCirc, numSubdivRad,
    fc, Kfc,
    CONCRETE_DENSITY,
    nRebars=68,
    rebar_dia=25,
    plot=True
    ):
    """
    Create a circular hollow confined‑concrete fiber section for OpenSees,
    with automatic placement of reinforcing bars along the mid‑thickness circle.

    Parameters
    ----------
    secTag          : int   – section identifier
    STEEL_TYPE      : str   – 'ELASTIC' or 'INELASTIC'
    Do              : float – outer diameter (mm)
    t               : float – wall thickness (mm)
    numSubdivCirc   : int   – number of circumferential subdivisions for concrete patch
    numSubdivRad    : int   – number of radial subdivisions for concrete patch
    fc              : float – unconfined concrete compressive strength (MPa, negative)
    Kfc             : float – ratio of confined to unconfined concrete strength
    CONCRETE_DENSITY: float – density of concrete (kg/m³) – converted to N·s²/mm³ internally
    plot            : bool  – if True, draw a 2D sketch of the section
    nRebars         : int   – total number of reinforcing bars (default = 68)
    rebar_dia       : float – diameter of each rebar (mm), uniform for all bars (default = 25)

    Returns
    -------
    Do              : float – outer diameter (mm)
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

    # Concrete01 materials (original choice; Concrete02 is commented out)
    ops.uniaxialMaterial('Concrete01', coreTag, fcC, ec0C, fcUC, ecuC)
    ops.uniaxialMaterial('Concrete01', coverTag, fcU, ec0U, fcUU, ecuU)
    # Alternative: Concrete02 (uncomment if needed)
    # ops.uniaxialMaterial('Concrete02', coreTag, fcC, ec0C, fcUC, ecuC, Lambda, ftC, EtsC)
    # ops.uniaxialMaterial('Concrete02', coverTag, fcU, ec0U, fcUU, ecuU, Lambda, ftU, EtsU)

    # ------------------------------------------------------------------
    # 3. SECTION GEOMETRY AND AREA
    # ------------------------------------------------------------------
    Di = Do - 2 * t
    r_out = Do / 2.0
    r_in = Di / 2.0

    AREA = np.pi * (r_out**2 - r_in**2)
    Iz = (np.pi / 4.0) * (r_out**4 - r_in**4)
    print(f"Section area A = {AREA:.6f} mm²")
    print(f"Moment of inertia Iz = {Iz:.8f} mm⁴")

    # ------------------------------------------------------------------
    # 4. FIBER SECTION – CONCRETE RING
    # ------------------------------------------------------------------
    ops.section('Fiber', secTag)
    # Circular patch for the entire wall (confined concrete)
    ops.patch('circ', coreTag, numSubdivRad, numSubdivCirc, 0.0, 0.0, r_in, r_out)

    # ------------------------------------------------------------------
    # 5. AUTOMATIC REBAR PLACEMENT ALONG THE MID‑THICKNESS CIRCLE
    # ------------------------------------------------------------------
    r_mid = (r_out + r_in) / 2.0   # radius at the centre of the wall

    rebars_coords = []
    for i in range(nRebars):
        theta = 2.0 * np.pi * i / nRebars
        x = r_mid * np.cos(theta)
        y = r_mid * np.sin(theta)
        rebars_coords.append((x, y))

    # Add each rebar as a fiber
    area_bar = np.pi * (rebar_dia ** 2) / 4.0
    for (x, y) in rebars_coords:
        ops.fiber(x, y, area_bar, steelTag)

    # ------------------------------------------------------------------
    # 6. PLOT (if requested)
    # ------------------------------------------------------------------
    if plot:
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.set_xlabel('Width (mm)')
        ax.set_ylabel('Height (mm)')
        ax.set_title('Circular Hollow Pipe Section with Rebars')
        ax.grid(True, linestyle='--', alpha=0.5)

        # Create circles for outer and inner boundaries
        theta_plot = np.linspace(0, 2 * np.pi, 100)
        x_out = r_out * np.cos(theta_plot)
        y_out = r_out * np.sin(theta_plot)
        x_in = r_in * np.cos(theta_plot)
        y_in = r_in * np.sin(theta_plot)

        # Plot outlines
        ax.plot(x_out, y_out, color='black', linewidth=1.5)
        ax.plot(x_in, y_in, color='black', linewidth=1.5)

        # Fill concrete (light gray) and hole (white)
        ax.fill_between(x_out, y_out, color='lightgray', alpha=0.5)
        ax.fill_between(x_in, y_in, color='white')

        # Plot rebars with numbering
        for idx, (x, y) in enumerate(rebars_coords, start=1):
            circ = patches.Circle(
                (x, y), radius=rebar_dia/2,
                edgecolor='red', facecolor='red', linewidth=1.5
            )
            ax.add_patch(circ)
            ax.text(
                x, y + rebar_dia/2 + 4, f'{idx}',
                color='purple', fontsize=6,
                ha='center', va='bottom', fontweight='bold'
            )

        # Axis limits
        max_dim = Do + 50
        ax.set_xlim(-max_dim/2, max_dim/2)
        ax.set_ylim(-max_dim/2, max_dim/2)
        ax.set_aspect('equal')
        plt.show()

    # ------------------------------------------------------------------
    # 7. MASS PER UNIT LENGTH
    # ------------------------------------------------------------------
    # Convert density from kg/m³ to kg/mm (1 m³ = 1e9 mm³)
    density_kg_per_mm3 = CONCRETE_DENSITY / 1e9   # kg/mm³
    ELE_MASS = density_kg_per_mm3 * AREA          # kg/mm

    return Do, ELE_MASS