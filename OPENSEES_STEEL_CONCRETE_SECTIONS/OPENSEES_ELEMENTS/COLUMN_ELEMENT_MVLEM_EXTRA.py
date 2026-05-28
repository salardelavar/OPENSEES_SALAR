def MVLEM_COLUMN_EXTRA_FUN(ELE_TAG, COL_DENSITY,
                       NODEi, NODEj,
                       Bsec, Hsec,
                       cover=40.0, rebar_dia=25.0,
                       nfibers=5, plot=True):
    """
    Create a confined RC column using the MVLEM element (OpenSees).
    The column has four longitudinal rebars (one at each corner),
    modelled by assigning the steel area only to the top and bottom macro‑fibers.

    Parameters
    ----------
    ELE_TAG        : int      – element tag
    WALL_DENSITY   : float    – mass density (kg/mm³), e.g. 2.5e‑9
    NODEi, NODEj   : int      – end nodes (i = bottom, j = top)
    Bsec, Hsec     : float    – column width (out‑of‑plane) and depth (in‑plane) [mm]
    cover          : float    – concrete cover to rebar centre [mm]
    rebar_dia      : float    – diameter of each longitudinal rebar [mm]
    nfibers        : int      – number of macro‑fibers (vertical strips)
    plot           : bool     – if True, draw the column cross‑section
    
    THIS PYTHON SCRIPT WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
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
    fcU = -25.0                 # [N/mm²] Concrete Compressive Strength
    EcU = 4700 * np.sqrt(-fcU)  # [N/mm^2] Concrete Elastic Modulus
    ec0U = 2*fcU/EcU            # [mm/mm] Concrete Compressive Strain
    fcUU = 0.2*fcU              # [N/mm²] Concrete Compressive Ultimate Strength
    ecuU = 5*ec0U               # [mm/mm] Concrete Compressive Ultimate Strain
    rtU = 1.2;			        # shape parameter - tension
    rcU = 7.0;			        # shape parameter - compression
    # Core concrete (confined)
    Kfc = 1.3;			       # ratio of confined to unconfined concrete strength
    fcC = Kfc*fcU              # [N/mm²] Concrete Compressive Strength
    EcC = 4700 * np.sqrt(-fcC) # [N/mm^2] Concrete Elastic Modulus
    ec0C = 2*fcC/EcC           # [mm/mm] Concrete Compressive Strain
    fcUC = 0.65*fcC            # [N/mm²] Concrete Compressive Ultimate Strength
    ecuC = 15*ec0C             # [mm/mm] Concrete Compressive Ultimate Strain
    Lambda = 0.1;	           # ratio between unloading slope
    rtC = 1.2;			       # shape parameter - tension
    rcC = 7.42;			       # shape parameter - compression    
    # tensile-strength properties
    ftC = 0.7 * np.sqrt(-fcC)  # [N/mm²] tensile strength +tension
    ftU = 0.7 * np.sqrt(-fcU)  # [N/mm²] tensile strength +tension
    EtsC = ftC/np.abs(ec0C)    # [N/mm²] tension softening stiffness
    EtsU = ftU/np.abs(ec0U)	   # [N/mm²] tension softening stiffness
    etUC = -2*ec0C             # [mm/mm] cracking strain - tension	
    etUU = -2*ec0U             # [mm/mm] cracking strain - tension	
    
    # STEEL
    # Reinforcing steel
    fy = 400          # [N/mm²] Steel Rebar Yield Strength   
    Es = 2e5          # [N/mm²] Modulus of Elasticity
    ey = fy/Es        # [mm/mm] Steel Rebar Yield Strain
    fu = 1.1818*fy    # [N/mm²] Steel Rebar Ultimate Strength
    esu = 0.09        # [mm/mm] Steel Rebar Ultimate Strain
    Esh = (fu - fy)/(esu - ey)
    Bs = Esh / Es
    R0 = 20.0         # initial value of curvature parameter
    a1 = 0.925        # curvature degradation parameter
    a2 = 0.15         # curvature degradation parameter

    # Unique material tags
    tagCover = int(ELE_TAG * 1000 + 1)
    tagCore  = int(ELE_TAG * 1000 + 2)
    tagSteel = int(ELE_TAG * 1000 + 3)

    ops.uniaxialMaterial('ConcreteCM', tagCover,
                         fcU, ec0U, EcU, rcU, ecuU,
                         ftU, -ec0U, rtU, etUU)
    ops.uniaxialMaterial('ConcreteCM', tagCore,
                         fcC, ec0C, EcC, rcC, ecuC,
                         ftC, -ec0C, rtC, etUC)
    ops.uniaxialMaterial('SteelMPF', tagSteel,
                         fy, fy, Es, Bs, Bs, R0, a1, a2)

    #%% GEOMETRY OF MACRO‑FIBERS (vertical strips)
    # Column depth (in‑plane) = Hsec, out‑of‑plane thickness = Bsec.
    # Each fiber has constant thickness = Bsec, and a width (in‑plane) that
    # sums to Hsec. The fibers represent layers from top to bottom.
    if nfibers < 3:
        raise ValueError("At least 3 fibers needed to model cover + core")

    core_depth = Hsec - 2.0 * cover
    if core_depth <= 0:
        raise ValueError("Cover too large for given column depth")

    # Distribute fiber widths: two cover fibers (top & bottom) and (nfibers-2) core fibers
    cover_width = cover
    core_width_per_fiber = core_depth / (nfibers - 2)

    widths = []
    # Top cover fiber
    widths.append(cover_width)
    # Core fibers
    for _ in range(nfibers - 2):
        widths.append(core_width_per_fiber)
    # Bottom cover fiber
    widths.append(cover_width)

    # Thickness (out‑of‑plane) is constant = Bsec for all fibers
    thickness = Bsec
    thicks = [thickness] * nfibers

    #%% REINFORCEMENT RATIO per fiber (to represent 4 rebars)
    # Total area of 4 rebars
    rebar_area = np.pi * (rebar_dia / 2.0) ** 2
    total_steel_area = 4 * rebar_area

    # Assign steel only to the top and bottom cover fibers (each gets 2 rebars)
    steel_area_per_cover_fiber = 2 * rebar_area
    rho_cover = steel_area_per_cover_fiber / (thickness * cover_width)

    rho = [0.0] * nfibers
    rho[0] = rho_cover          # top fiber
    rho[-1] = rho_cover         # bottom fiber
    # Core fibers have rho = 0

    #%% SHEAR MATERIAL (elastic stiffness)
    matShearTag = int(ELE_TAG * 3000)
    AREA = Bsec * Hsec
    Iz = Bsec * (Hsec ** 3) / 12.0
    Avy = (5.0 / 6.0) * AREA
    E_mod = 4700 * np.sqrt(-fcC)          # MPa
    nu = 0.2
    G_mod = E_mod / (2.0 * (1.0 + nu))
    shear_stiffness = AREA * G_mod
    ops.uniaxialMaterial('Elastic', matShearTag, shear_stiffness)

    #%% ASSIGN MATERIALS TO FIBERS
    matConcrete = []
    matSteel = []
    for i in range(nfibers):
        if i == 0 or i == nfibers - 1:
            matConcrete.append(tagCover)   # cover fibers (unconfined)
        else:
            matConcrete.append(tagCore)    # core fibers (confined)
        matSteel.append(tagSteel)

    #%% CREATE MVLEM ELEMENT
    c = 0.4          # location of centre of rotation (recommended)
    ops.element('MVLEM', ELE_TAG, COL_DENSITY,
                NODEi, NODEj,
                nfibers, c,
                '-thick',  *thicks,
                '-width',  *widths,
                '-rho',    *rho,
                '-matConcrete', *matConcrete,
                '-matSteel',    *matSteel,
                '-matShear',    matShearTag)

    print(f"MVLEM element {ELE_TAG} created for column {Bsec}x{Hsec} mm, "
          f"{nfibers} fibers, 4 rebars dia={rebar_dia} mm, cover={cover} mm")

    #%% PLOT SECTION (if requested)
    if plot:
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.set_xlabel('Width (out‑of‑plane) [mm]')
        ax.set_ylabel('Depth (in‑plane) [mm]')
        ax.set_title(f'Confined Column Section (MVLEM)\n'
                     f'4 rebars Ø{rebar_dia} mm, cover = {cover} mm')
        ax.grid(True, ls='--', alpha=0.5)

        # Outer column boundary
        rect = patches.Rectangle((-Bsec/2, -Hsec/2), Bsec, Hsec,
                                 linewidth=1.5, edgecolor='black',
                                 facecolor='lightgray', label='Concrete')
        ax.add_patch(rect)

        # Indicate core region (confined)
        core_rect = patches.Rectangle((-Bsec/2, -Hsec/2 + cover),
                                      Bsec, Hsec - 2*cover,
                                      linewidth=1, edgecolor='gray',
                                      facecolor='gray', linestyle='--',
                                      label='Confined core')
        ax.add_patch(core_rect)

        # Plot rebar locations (four corners)
        rebar_radius = rebar_dia / 2.0
        corners = [(-Bsec/2 + cover, -Hsec/2 + cover),   # bottom‑left
                   ( Bsec/2 - cover, -Hsec/2 + cover),   # bottom‑right
                   (-Bsec/2 + cover,  Hsec/2 - cover),   # top‑left
                   ( Bsec/2 - cover,  Hsec/2 - cover)]   # top‑right
        for (x, y) in corners:
            circle = patches.Circle((x, y), rebar_radius,
                                    facecolor='red', edgecolor='red',
                                    linewidth=0.8)
            ax.add_patch(circle)

        ax.set_xlim(-Bsec/2 - 20, Bsec/2 + 20)
        ax.set_ylim(-Hsec/2 - 20, Hsec/2 + 20)
        ax.set_aspect('equal')
        ax.legend(loc='upper right')
        plt.show()

    # Optionally return the section total steel ratio
    total_rho = total_steel_area / AREA * 100.0
    return total_rho