def COMPOSITE_REC_CROSS_SHAPED_STEEL_SECTION_FUN_EXTRA(secTag, HSec, BSec, cover,
                              STEEL_TYPE, fc, Kc,
                              STEEL_DENSITY,
                              CONCRETE_DENSITY, PLOT):
 
    """
    Create a rectangular confined‑concrete fiber section (OpenSees) and,
    if PLOT=True, draw a simple 2‑D sketch of the section.

    Parameters
    ----------
    secTag          : int   – section identifier
    HSec, BSec      : float – height and width  (mm)
    cover           : float – concrete cover (mm)
    STEEL_TYPE      : str   – 'ELASTIC' or 'INELASTIC'
    fc              : float – unconfined concrete compressive strength (MPa, negative)
    Kfc             : float – ratio of confined to unconfined concrete strength
    CONCRETE_DENSITY: float – ρc (kg/m³) – will be converted to N·s²/mm³
    STEEL_DENSITY   : float – ρc (kg/m³) – will be converted to N·s²/mm³
    PLOT            : bool  – draw the section if True
    
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
    Kfc = 1.3;			      # ratio of confined to unconfined concrete strength
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
    
    coreTag, coverTag, steelTag, steelITag = secTag + 100, secTag + 200, secTag + 300, secTag + 400
    if STEEL_TYPE == 'ELASTIC':
        ops.uniaxialMaterial('Elastic', steelTag, Es)      # REBAR
        ops.uniaxialMaterial('Elastic', steelTag, EsI)     # I SECION
    if STEEL_TYPE == 'INELASTIC':   
        pinchX = 0.8   # Pinching factor in X direction
        pinchY = 0.5   # Pinching factor in Y direction
        damage1 = 0.0  # Damage due to ductility
        damage2 = 0.0  # Damage due to energy
        beta = 0.1     # Stiffness degradation parameter
        # REBAR
        ops.uniaxialMaterial('Hysteretic', steelTag,
                                        fy, ey,
                                        fu, esu,
                                        0.2*fu, 1.1*esu,
                                        -fy, -ey,
                                        -fu, -esu,
                                        -0.2*fu, -1.1*esu,
                                        pinchX, pinchY,
                                        damage1, damage2, beta)
        # I SECION
        ops.uniaxialMaterial('Hysteretic', steelITag,
                                        fyI, eyI,
                                        fuI, esuI,
                                        0.2*fuI, 1.1*esuI,
                                        -fyI, -eyI,
                                        -fuI, -esuI,
                                        -0.2*fuI, -1.1*esuI,
                                        pinchX, pinchY,
                                        damage1, damage2, beta)
        # INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Hysteretic_Material
        
    #ops.uniaxialMaterial('Concrete01', coreTag, fcC, ec0C, fcUC, ecuC)  # Core concrete (confined)
    #ops.uniaxialMaterial('Concrete01', coverTag, fcU, ec0U, fcUU, ecuU) # Cover concrete (unconfined)
    ops.uniaxialMaterial('Concrete02', coreTag, fcC, ec0C, fcUC, ecuC, LambdaC, ftC, EtsC) # build core concrete (confined)
    ops.uniaxialMaterial('Concrete02', coverTag, fcU, ec0U, fcUU, ecuU, LambdaU, ftU, EtsU) # build cover concrete (unconfined)
    
    # Some variables derived from the parameters
    y1 = HSec / 2.0
    z1 = BSec / 2.0
    NUMFIBERS = 20  # Number of layers for each fiber
    
    ops.section('Fiber', secTag)
    # Create the concrete core fibers
    ops.patch('rect', coreTag, NUMFIBERS, NUMFIBERS, cover - y1, cover - z1, y1 - cover, z1 - cover)
    
    # Create the concrete cover fibers (top, bottom, left, right)
    ops.patch('rect', coverTag, NUMFIBERS, 10, -y1, z1 - cover, y1, z1)
    ops.patch('rect', coverTag, NUMFIBERS, 10, -y1, -z1, y1, cover - z1)
    ops.patch('rect', coverTag, NUMFIBERS, 10, -y1, cover - z1, cover - y1, z1 - cover)
    ops.patch('rect', coverTag, NUMFIBERS, 10, y1 - cover, cover - z1, y1, z1 - cover)
    
    REBAR = 25.0 # [mm] Rebar Diameter
    PT = 10.0    # [mm] Plate Thickness
    # -------------------- Rebars -------------------------------
    # (diameter [mm], y‑coord [mm], x‑coord [mm])
    rebars = [
        (REBAR,  HSec/2 - cover, -BSec/2 + cover),    # 1
        (REBAR,  HSec/2 - cover, BSec/2 - cover),     # 2
        (REBAR,  -HSec/2 + cover, -BSec/2 + cover),   # 3
        (REBAR,  -HSec/2 + cover, +BSec/2 - cover),   # 4
        (16.0, HSec/4 - cover, -BSec/2 + cover),      # 5
        (16.0, HSec/4 - cover, +BSec/2 - cover),      # 6
        (16.0,  -HSec/4 + cover, -BSec/2 + cover),    # 7
        (16.0,  -HSec/4 + cover, +BSec/2 - cover),    # 8
        (REBAR,  0.0, -BSec/2 + cover),               # 9
        (REBAR,  0.0, +BSec/2 - cover),               # 10
        
        (REBAR,  HSec/2 - cover, 0.0),                # 11
        (REBAR,  -HSec/2 + cover, 0.0),               # 12
    ]

    for dia, y, x in rebars:
        area = np.pi * dia**2 / 4.0          # mm²
        ops.fiber(x, y, area, steelTag)      # add to OpenSees model
    
    # --------- Define Fibers for Flanges and Web ---------
    def add_rectangle_patch(y_bottom, y_top, x_left, x_right, num_y_divisions=NUMFIBERS, num_x_divisions=NUMFIBERS):
        """Helper function to add rectangular patches of fibers."""
        ops.patch('rect', steelITag, num_y_divisions, num_x_divisions, x_left, y_bottom, x_right, y_top)
    
    # -------------------- Steel Section -------------------------------
    # (layers depth[mm], layers width[mm], center y‑coord [mm], center x‑coord [mm], numSubdivY, numSubdivZ)
    steel_layers = [
        (10.0, 150.0, -125, 0.0, 10, 10),    # 1
        (240.0, 10.0, 0, 0.0, 10, 10),       # 2
        (10.0, 150.0, 125.0, 0.0, 10, 10),   # 3    
         
        (10.0, 105.0, 0.0, -57.5, 10, 10),    # 4
        (10.0, 105.0, 0.0, +57.5, 10, 10),    # 5
        
        (150.0, 10.0, 0.0, -115.0, 10, 10),   # 6
        (150.0, 10.0, 0.0, +115.0, 10, 10),   # 7
        
        (PT, BSec+2*PT, +0.5*HSec+0.5*PT, 0.0, 10, 10),  # 8 -> PLATE
        (PT, BSec+2*PT, -0.5*HSec-0.5*PT, 0.0, 10, 10),  # 9 -> PLATE
        (HSec, PT, 0.0, +0.5*BSec+0.5*PT, 10, 10),       # 10-> PLATE   
        (HSec, PT, 0.0, -0.5*BSec-0.5*PT, 10, 10),       # 11-> PLATE
        
    ]
    for depth, width, center_y, center_x, numSubdivY, numSubdivZ in steel_layers:
        x_left = center_x - width/2
        x_right = center_x + width/2
        y_bot = center_y - depth/2
        y_top = center_y + depth/2
        ops.patch('rect', steelITag, numSubdivY, numSubdivZ, x_left, y_bot, x_right, y_top)

    # -------------------- Section area -------------------------
    # Initialize
    min_bottom = float('inf')
    max_top = -float('inf')
    TOTAL_AREA = 0.0
    
    for depth, width, center_y, center_x,_ ,_ in steel_layers:
        bottom = center_y - depth/2
        top = center_y + depth/2
        min_bottom = min(min_bottom, bottom)
        max_top = max(max_top, top)
        TOTAL_AREA += depth * width
    
    SECTION_HEIGHT = max_top - min_bottom
    
    print(f"Steel Section height = {SECTION_HEIGHT} mm")
    print(f"Steel Section area   = {TOTAL_AREA} mm²")
    
    # Calculate mass
    MASS = HSec * BSec * CONCRETE_DENSITY + TOTAL_AREA * STEEL_DENSITY # kg/mm
    
    if PLOT:
        fig, ax = plt.subplots(figsize=(6, 6))
        # concrete block
        rect = plt.Rectangle((-z1, -y1), BSec, HSec,
                             facecolor='lightgray', edgecolor='k', lw=1.5, label='Unconfined')
        ax.add_patch(rect)

        # cover outline (just for visual aid)
        inner = plt.Rectangle((-z1 + cover, -y1 + cover),
                              BSec - 2 * cover, HSec - 2 * cover,
                              facecolor='gray', label='Confined')#, edgecolor='k', ls='--'
        ax.add_patch(inner)

        #fig, ax = plt.subplots(figsize=(10, 8))
        ax.set_xlabel('Width (mm)')
        ax.set_ylabel('Height (mm)')
        ax.set_title('Steel section (7 layers)')
        ax.grid(True, ls='--', alpha=0.5)
        
        # Steel geometry derived from steel_layers
        ii = 0
        for depth, width, center_y, center_x,_ ,_ in steel_layers:
            x_left = center_x - width/2
            x_right = center_x + width/2
            y_bot = center_y - depth/2
            y_top = center_y + depth/2
            rect = patches.Rectangle(
                (x_left, y_bot), width, depth,
                linewidth=1, edgecolor='black', facecolor='lightblue'
            )
            ax.add_patch(rect)
            ii += 1
            ax.text(
                center_x, center_y, f'{ii}',
                color='purple', fontsize=8,
                ha='center', va='bottom', fontweight='bold'
            )

        # Rebars – red circles + numbers
        for i, (dia, y, x) in enumerate(rebars, start=1):
            # circle
            circ = patches.Circle((x, y), radius=dia/2,
                                 edgecolor='red', facecolor='red',
                                 linewidth=2)
            ax.add_patch(circ)

            # label – placed a little above the bar
            ax.text(x, y + dia/2 + 4, f'{i}',
                    color='purple', fontsize=6,
                    ha='center', va='bottom',
                    fontweight='bold')

        ax.set_aspect('equal')
        ax.set_xlabel('Width  (mm)')
        ax.set_ylabel('Height (mm)')
        ax.set_title(f'Confined Concrete Section with Cross Shaped Section - tag {secTag}')
 
        max_dim = max(BSec, HSec) + 50
        ax.set_xlim(-max_dim/2, max_dim/2)
        ax.set_ylim(-max_dim/2, max_dim/2)
        ax.grid(True, ls=':', alpha=0.5)
        ax.legend()
        plt.show()

    return HSec+2*PT, MASS
