def DUMBBELL_SHEAR_WALL_FUN(secTag, STEEL_TYPE, fc, Kfc,
                            L, t, Le, nFib,
                            CONCRETE_DENSITY, plot=True):
    """
    Create a dumbbell-shaped reinforced concrete shear wall fiber section
    (confined boundary elements + unconfined web) for OpenSees.

    Parameters
    ----------
    secTag          : int     – section identifier
    STEEL_TYPE      : str     – 'ELASTIC' or 'INELASTIC'
    fc              : float   – unconfined concrete compressive strength (MPa, *negative*)
    Kfc             : float   – confined / unconfined strength ratio (for end zones)
    L               : float   – total length of the wall (horizontal, mm)
    t               : float   – total thickness of the wall (vertical, mm)
    Le              : float   – length of each confined end zone (from both ends, mm)
    nFib            : int     – number of fibers along each direction in a rectangular patch
    CONCRETE_DENSITY: float   – concrete density (kg/m³)
    plot            : bool    – if True, draw a 2D sketch of the section

    Returns
    -------
    L, mass         : (float, float) – wall length (mm) and element mass (kg)
    
    THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
    """
    import numpy as np
    import openseespy.opensees as ops
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches

    # Define materials for nonlinear elements
    # Define parameters (units: mm, N)
    # ------------------------------------------
    # CONCRETE                  tag   f'c        ec0   f'cu        ecu
    # Cover concrete (unconfined)
    fcU = -fc                 # [N/mm²] Concrete Compressive Strength
    Ec = 4700 * np.sqrt(-fcU) # [N/mm^2] Concrete Elastic Modulus
    ec0U = 2*fcU/Ec           # [mm/mm] Concrete Compressive Strain
    fcUU = 0.2*fcU            # [N/mm²] Concrete Compressive Ultimate Strength
    ecuU = 5*ec0U             # [mm/mm] Concrete Compressive Ultimate Strain
    
    # Core concrete (confined)
    #Kfc = 1.3;			      # ratio of confined to unconfined concrete strength
    fcC = Kfc*fcU             # [N/mm²] Concrete Compressive Strength
    Ec = 4700 * np.sqrt(-fcC) # [N/mm^2] Concrete Elastic Modulus
    ec0C = 2*fcC/Ec           # [mm/mm] Concrete Compressive Strain
    fcUC = 0.65*fcC           # [N/mm²] Concrete Compressive Ultimate Strength
    ecuC = 15*ec0C            # [mm/mm] Concrete Compressive Ultimate Strain
    Lambda = 0.1;	          # ratio between unloading slope
    # tensile-strength properties
    ftC = 0.7 * np.sqrt(-fcC)  # [N/mm²] tensile strength +tension
    ftU = 0.7 * np.sqrt(-fcU)  # [N/mm²] tensile strength +tension
    EtsC = ftC/np.abs(ec0C)    # [N/mm²] tension softening stiffness
    EtsU = ftU/np.abs(ec0U)	   # [N/mm²] tension softening stiffness
    
    # STEEL
    # Reinforcing steel
    fy = 400          # [N/mm²] Steel Rebar Yield Strength   
    Es = 2e5          # [N/mm²] Modulus of Elasticity
    ey = fy/Es        # [mm/mm] Steel Rebar Yield Strain
    fu = 1.1818*fy    # [N/mm²] Steel Rebar Ultimate Strength
    esu = 0.09        # [mm/mm] Steel Rebar Ultimate Strain
    Esh = (fu - fy)/(esu - ey)
    Bs = Esh / Es
    
    coreTag, coverTag, steelTag = secTag + 100, secTag + 200, secTag + 300
    if STEEL_TYPE == 'ELASTIC':
        ops.uniaxialMaterial('Steel01', steelTag, fy, Es, Bs) 
    if STEEL_TYPE == 'INELASTIC':   
        pinchX = 0.8   # Pinching factor in X direction
        pinchY = 0.5   # Pinching factor in Y direction
        damage1 = 0.0  # Damage due to ductility
        damage2 = 0.0  # Damage due to energy
        beta = 0.1     # Stiffness degradation parameter
        ops.uniaxialMaterial('Hysteretic', steelTag, fy, ey, fu, esu, 0.2*fu, 1.1*esu, -fy, -ey, -fu, -esu, -0.2*fu, -1.1*esu, pinchX, pinchY, damage1, damage2, beta)
        # INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Hysteretic_Material
        
    #ops.uniaxialMaterial('Concrete01', coreTag, fcC, ec0C, fcUC, ecuC)  # Core concrete (confined)
    #ops.uniaxialMaterial('Concrete01', coverTag, fcU, ec0U, fcUU, ecuU) # Cover concrete (unconfined)
    ops.uniaxialMaterial('Concrete02', coreTag, fcC, ec0C, fcUC, ecuC, Lambda, ftC, EtsC) # build core concrete (confined)
    ops.uniaxialMaterial('Concrete02', coverTag, fcU, ec0U, fcUU, ecuU, Lambda, ftU, EtsU)

    # ---------- Section geometry and area ----------
    area_total = L * t
    print(f"Total cross-sectional area = {area_total:.2f} mm²")

    # Section centroid is at (0,0) – symmetric
    print(f"Neutral axis (horizontal) lies at x = 0 (centred)")

    # ---------- Build fiber section ----------
    ops.section('Fiber', secTag)

    # Helper: add a rectangular patch of concrete
    def add_rect(mat, x_left, x_right, y_bot, y_top):
        ops.patch('rect', mat, nFib, nFib,
                  x_left, y_bot, x_right, y_top)

    half_t = t / 2.0
    half_L = L / 2.0
    left_boundary_x1 = -half_L
    left_boundary_x2 = -half_L + Le
    right_boundary_x1 = half_L - Le
    right_boundary_x2 = half_L
    web_x1 = left_boundary_x2
    web_x2 = right_boundary_x1

    # Left confined end zone
    add_rect(coreTag, left_boundary_x1, left_boundary_x2, -half_t, half_t)
    # Right confined end zone
    add_rect(coreTag, right_boundary_x1, right_boundary_x2, -half_t, half_t)
    # Unconfined web
    if web_x2 > web_x1:
        add_rect(coverTag, web_x1, web_x2, -half_t, half_t)

    # ---------- Reinforcement (grid within each end zone) ----------
    rebar_dia = 25.0          # mm
    area_rebar = np.pi * rebar_dia**2 / 4.0   # mm²

    # Number of layers in vertical (y) and horizontal (x) direction per end zone
    n_rows_y = 3   # top and bottom
    n_cols_x = 10   # evenly spaced along the length of the end zone

    rebars = []   # store (x, y, area) for each rebar fiber

    # Left end zone rebars
    y_positions = np.linspace(-half_t + half_t/(n_rows_y+1),
                               half_t - half_t/(n_rows_y+1),
                               n_rows_y)
    x_positions = np.linspace(left_boundary_x1 + Le/(n_cols_x+1),
                              left_boundary_x2 - Le/(n_cols_x+1),
                              n_cols_x)
    for x in x_positions:
        for y in y_positions:
            rebars.append((x, y, area_rebar))

    # Right end zone rebars
    x_positions = np.linspace(right_boundary_x1 + Le/(n_cols_x+1),
                              right_boundary_x2 - Le/(n_cols_x+1),
                              n_cols_x)
    for x in x_positions:
        for y in y_positions:
            rebars.append((x, y, area_rebar))

    # Add all rebars to the section
    for x, y, area in rebars:
        ops.fiber(x, y, area, steelTag)

    # ---------- Optional plot ----------
    if plot:
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.set_xlabel('Width (horizontal, mm)')
        ax.set_ylabel('Thickness (vertical, mm)')
        ax.set_title('Dumbbell Shear Wall Section with Rebars')
        ax.grid(True, linestyle='--', alpha=0.5)

        # Draw concrete parts
        # Left confined zone
        rect_left = patches.Rectangle((left_boundary_x1, -half_t),
                                      Le, t,
                                      linewidth=1.5, edgecolor='black',
                                      facecolor='lightgray')
        ax.add_patch(rect_left)
        # Right confined zone
        rect_right = patches.Rectangle((right_boundary_x1, -half_t),
                                       Le, t,
                                       linewidth=1.5, edgecolor='black',
                                       facecolor='lightgray')
        ax.add_patch(rect_right)
        # Web (if visible)
        if web_x2 > web_x1:
            rect_web = patches.Rectangle((web_x1, -half_t),
                                         web_x2 - web_x1, t,
                                         linewidth=1.5, edgecolor='black',
                                         facecolor='whitesmoke')
            ax.add_patch(rect_web)

        # Draw rebars as red circles with numbering
        for idx, (x, y, _) in enumerate(rebars, start=1):
            circ = patches.Circle((x, y), radius=rebar_dia/2,
                                  edgecolor='red', facecolor='red', linewidth=1.5)
            ax.add_patch(circ)
            ax.text(x, y + rebar_dia/2 + 3, f'{idx}',
                    color='purple', fontsize=7, ha='center', va='bottom',
                    fontweight='bold')

        ax.set_xlim(-half_L - 20, half_L + 20)
        ax.set_ylim(-half_t - 30, half_t + 30)
        ax.set_aspect('equal')
        plt.show()

    # ---------- Mass calculation ----------
    density_per_mm3 = CONCRETE_DENSITY / 1e9   # kg/mm³ (since 1 kg/m³ = 1e-9 kg/mm³)
    mass = area_total * density_per_mm3        # kg per mm length? Actually mass per unit length
    # The return is consistent with original: (length, mass per mm)
    return L, mass