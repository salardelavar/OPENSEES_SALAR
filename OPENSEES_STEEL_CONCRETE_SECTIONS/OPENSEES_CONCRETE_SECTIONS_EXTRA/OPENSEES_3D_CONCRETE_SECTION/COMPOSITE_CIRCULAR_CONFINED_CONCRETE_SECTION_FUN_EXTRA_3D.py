def COMPOSITE_CIRCULAR_CONFINED_CONCRETE_SECTION_FUN_EXTRA_3D(secTag, RI, RO, COVER,
                                                 fc, Kfc, THICKNESS,
                                                 STEEL_TYPE, STEEL_DENSITY, CONCRETE_DENSITY,
                                                 plot=True):
    """
    Create a Circular Composite with confined concrete section with an external steel pipe.
    Materials are generated and the section is built in OpenSees.
    
    Returns
    -------
    SECTION_HEIGHT : float   – total outer diameter (mm) = 2*(RO + THICKNESS)
    ELE_MASS       : float   – mass per unit length (kg/mm)  (can be scaled as needed)
    
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
    fcU = -fc                 # [N/mm²] Concrete Compressive Strength
    Ec = 4700 * np.sqrt(-fcU) # [N/mm^2] Concrete Elastic Modulus
    ec0U = 2*fcU/Ec           # [mm/mm] Concrete Compressive Strain
    fcUU = 0.2*fcU            # [N/mm²] Concrete Compressive Ultimate Strength
    ecuU = 5*ec0U             # [mm/mm] Concrete Compressive Ultimate Strain
    LambdaU = 0.1;	          # ratio between unloading slope    
    # Core concrete (confined)
    #Kfc = 1.3;			      # ratio of confined to unconfined concrete strength
    fcC = Kfc*fcU             # [N/mm²] Concrete Compressive Strength
    Ec = 4700 * np.sqrt(-fcC) # [N/mm^2] Concrete Elastic Modulus
    ec0C = 2*fcC/Ec           # [mm/mm] Concrete Compressive Strain
    fcUC = 0.65*fcC           # [N/mm²] Concrete Compressive Ultimate Strength
    ecuC = 15*ec0C            # [mm/mm] Concrete Compressive Ultimate Strain
    LambdaC = 0.1;	          # ratio between unloading slope
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
    # I Section
    fyI = 240          # [N/mm²] Steel I Section Yield Strength   
    EsI = 2e5          # [N/mm²] I Section Modulus of Elasticity
    eyI = fyI/EsI      # [mm/mm] Steel I Section Yield Strain
    fuI = 1.1818*fyI   # [N/mm²] Steel I Section Ultimate Strength
    esuI = 0.25        # [mm/mm] Steel I Section Ultimate Strain
    EshI = (fuI - fyI)/(esuI - eyI)
    BsI = EshI / EsI

    # Material tags
    coreTag   = secTag + 100
    coverTag  = secTag + 200
    steelTag  = secTag + 300   # not used for patches
    steelITag = secTag + 400

    # Define steel materials
    if STEEL_TYPE == 'ELASTIC':
        ops.uniaxialMaterial('Elastic', steelTag, Es)
        ops.uniaxialMaterial('Elastic', steelITag, EsI)
    elif STEEL_TYPE == 'INELASTIC':
        pinchX = 0.8
        pinchY = 0.5
        damage1 = 0.0
        damage2 = 0.0
        beta = 0.1
        ops.uniaxialMaterial('Hysteretic', steelTag,
                             fy, ey, fu, esu, 0.2*fu, 1.1*esu,
                             -fy, -ey, -fu, -esu, -0.2*fu, -1.1*esu,
                             pinchX, pinchY, damage1, damage2, beta)
        ops.uniaxialMaterial('Hysteretic', steelITag,
                             fyI, eyI, fuI, esuI, 0.2*fuI, 1.1*esuI,
                             -fyI, -eyI, -fuI, -esuI, -0.2*fuI, -1.1*esuI,
                             pinchX, pinchY, damage1, damage2, beta)

    # Concrete materials
    ops.uniaxialMaterial('Concrete02', coreTag, fcC, ec0C, fcUC, ecuC, LambdaC, ftC, EtsC)
    ops.uniaxialMaterial('Concrete02', coverTag, fcU, ec0U, fcUU, ecuU, LambdaU, ftU, EtsU)

    nfCoreR = 8
    nfCoreT = 8
    nfCoverR = 4
    nfCoverT = 8
    # -------------------- Fiber section -------------------------
    ops.section('Fiber', secTag, '-GJ', 10e6)
    # -------------------- Concrete Section -------------------------------
    # (layers depth[mm], layers width[mm], center y‑coord [mm], center x‑coord [mm], numSubdivY, numSubdivZ, tag, color)
    CD = RO # [mm] Section Depth and Width
    NFY, NFX = 10, 10
    mat_layers = [
        # STEEL SECTION
        (10.0, 150.0, 0.50*CD, 0.0, NFY, NFX, steelITag, 'lime'),                 # 01 -> CROSSED SECTION
        (240.0, 10.0, 0.0, 0.0, NFY, NFX, steelITag, 'lime'),                     # 02 -> CROSSED SECTION
        (10.0, 150.0, -0.50*CD, 0.0, NFY, NFX, steelITag, 'lime'),                # 03 -> CROSSED SECTION   
         
        (10.0, 105.0, 0.0, 0.0-57.5, NFY, NFX, steelITag, 'lime'),                # 04 -> CROSSED SECTION
        (10.0, 105.0, 0.0, 0.0+57.5, NFY, NFX, steelITag, 'lime'),                # 05 -> CROSSED SECTION
        
        (150.0, 10.0, 0.0, 0.0-115.0, NFY, NFX, steelITag, 'lime'),               # 06 -> CROSSED SECTION
        (150.0, 10.0, 0.0, 0.0+115.0, NFY, NFX, steelITag, 'lime'),               # 07 -> CROSSED SECTION
    ]

    rc = RO - COVER          # core radius
    # Concrete core
    ops.patch('circ', coreTag, nfCoreT, nfCoreR, 0.0, 0.0, RI, rc, 0.0, 360.0)
    # Concrete cover
    ops.patch('circ', coverTag, nfCoverT, nfCoverR, 0.0, 0.0, rc, RO, 0.0, 360.0)
    # Steel pipe (outer tube)
    ops.patch('circ', steelITag, nfCoverT, nfCoverR, 0.0, 0.0, RO, RO+THICKNESS, 0.0, 360.0)

    # ---------- Compute section height and mass ----------
    SECTION_HEIGHT = 2.0 * (RO + THICKNESS)   # outer diameter (mm)
    
    # Cross-sectional areas (mm²)
    CONCRETE_AREA = np.pi * (RO**2 - RI**2)                    # core + cover
    STEEL_PIPE_AREA = np.pi * ((RO + THICKNESS)**2 - RO**2)    # pipe only


    # Mass per unit length (kg/mm)
    ELE_MASS = CONCRETE_AREA * CONCRETE_DENSITY + STEEL_PIPE_AREA * STEEL_DENSITY

    print(f"Section height = {SECTION_HEIGHT} mm")
    # -------------------- Rebars -------------------------------
    # (diameter [mm], y‑coord [mm], x‑coord [mm])
    """
    num_bars = 16
    dtheta = 360.0 / num_bars  # 22.5 degrees
    angles_deg = np.arange(num_bars) * dtheta  # [0,22.5,45,...,337.5]
    angles_rad = np.deg2rad(angles_deg)
    x_coords = rc * np.cos(angles_rad)
    y_coords = rc * np.sin(angles_rad)
    
    diameters = [25.0, 16.0, 25.0, 16.0, 25.0, 16.0,
                 25.0, 16.0, 16.0, 16.0, 16.0, 16.0,
                 16.0, 16.0, 16.0, 16.0]
    
    rebars = [(diameters[i], x_coords[i], y_coords[i]) for i in range(num_bars)]
    for bar in rebars:
    print(f"    ({bar[0]}, {bar[1]:.4f}, {bar[2]:.4f}),")
    """
    rebars = [
        (25.0, 200.0000, 0.0000),      # rebar 1   (0°)
        (16.0, 184.7759, 76.5367),     # rebar 2   (22.5°)
        (25.0, 0.0000, 200.0000),      # rebar 3   (90°)
        (16.0, -184.7759, 76.5367),    # rebar 4   (112.5°)
        (25.0, -200.0000, 0.0000),     # rebar 5   (180°)
        (16.0, -184.7759, -76.5367),   # rebar 6   (202.5°)
        (25.0, 0.0000, -200.0000),     # rebar 7   (270°)
        (16.0, 184.7759, -76.5367),    # rebar 8   (292.5°)
        (16.0, 141.4214, -141.4214),   # rebar 9   (315°)
        (16.0, 76.5367, -184.7759),    # rebar 10  (337.5°)
        (16.0, -76.5367, -184.7759),   # rebar 11  
        (16.0, -141.4214, -141.4214),  # rebar 12
        (16.0, -76.5367, 184.7759),    # rebar 13
        (16.0, 141.4214, 141.4214),    # rebar 14
        (16.0, -141.4214, 141.4214),   # rebar 15
        (16.0, 76.5367, 184.7759),     # rebar 16
    ]
    
    for dia, y, x in rebars:
        area = np.pi * dia**2 / 4.0          # mm²
        ops.fiber(x, y, area, steelTag)      # add to OpenSees model    
    # -------------------- Plot (optional) -----------------------
    if plot:
        fig, ax = plt.subplots(figsize=(8,8))
        # Draw concrete regions (core and cover)
        # Core: inner radius RI to rc
        if RI > 0:
            core_inner = patches.Circle((0,0), RI, color='white', zorder=1)
            ax.add_patch(core_inner)
        core_disk = patches.Circle((0,0), rc, linewidth=1, edgecolor='black', facecolor='lightgray')
        ax.add_patch(core_disk)
        # Cover: rc to RO
        cover_ring = patches.Wedge((0,0), RO, 0, 360, width=RO-rc, edgecolor='black', facecolor='darkgray')
        ax.add_patch(cover_ring)
        # Steel pipe: RO to RO+THICKNESS
        pipe_ring = patches.Wedge((0,0), RO+THICKNESS, 0, 360, width=THICKNESS,
                                  linewidth=1, edgecolor='black', facecolor='cyan')
        ax.add_patch(pipe_ring)
        
        # Material geometry derived from mat_layers
        ii = 0
        for depth, width, center_y, center_x, _, _,_, COLOR in mat_layers:
            x_left = center_x - width/2
            x_right = center_x + width/2
            y_bot = center_y - depth/2
            y_top = center_y + depth/2
            rect = patches.Rectangle(
                (x_left, y_bot), width, depth,
                linewidth=1, edgecolor='black', facecolor=COLOR
            )
            ax.add_patch(rect)
            ii += 1
            ax.text(
                center_x, center_y, f'{ii}',
                color='blue', fontsize=8,
                ha='center', va='bottom', fontweight='bold'
            )
        
        # Rebars – red circles + numbers (using your existing rebars list)
        for i, (dia, y, x) in enumerate(rebars, start=1):
            circ = patches.Circle(
                (x, y), radius=dia/2,
                edgecolor='red', facecolor='red', linewidth=2
            )
            ax.add_patch(circ)
            ax.text(
                x, y + dia/2 + 5, f'{i}',
                color='purple', fontsize=8,
                ha='center', va='bottom', fontweight='bold'
            )

        # Annotations
        ax.set_xlim(-(RO+THICKNESS)*1.1, (RO+THICKNESS)*1.1)
        ax.set_ylim(-(RO+THICKNESS)*1.1, (RO+THICKNESS)*1.1)
        ax.set_aspect('equal')
        ax.set_title(f'Circular Composite Section \n Section Diameter: {SECTION_HEIGHT:.0f} mm - Section Tag: {secTag}')
        ax.set_xlabel('Width (mm)')
        ax.set_ylabel('Height (mm)')
        plt.grid(True, linestyle=':', alpha=0.5)
        plt.show()

    return SECTION_HEIGHT, ELE_MASS