def STEEL_BOX_SECTION_DIFF_WIDTH_HEIGHT_FUN_EXTRA(secTag, STEEL_TYPE,STEEL_DENSITY, plot=True):
    """
    Create a steel with different widths and heights, fiber section (OpenSees) and,
    if PLOT=True, draw a simple 2‑D sketch of the section.

    Parameters
    ----------
    secTag          : int   – section identifier
    STEEL_TYPE      : str   – 'ELASTIC' or 'INELASTIC'
    STEEL_DENSITY: float – ρc (kg/m³) – will be converted to N·s²/mm³
    PLOT            : bool  – draw the section if True
    
    THIS PYTHON SCRIPT WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
    """
    import numpy as np
    import openseespy.opensees as ops
    import matplotlib.pyplot as plt
    import matplotlib.pyplot as plt, matplotlib.patches as patches
    
    # Define materials for nonlinear elements
    # Define parameters (units: mm, N)
    # ------------------------------------------    
    # BOX Section
    fyI = 240          # [N/mm²] Steel I Section Yield Strength   
    EsI = 2e5          # [N/mm²] I Section Modulus of Elasticity
    eyI = fyI/EsI      # [mm/mm] Steel I Section Yield Strain
    fuI = 1.1818*fyI   # [N/mm²] Steel I Section Ultimate Strength
    esuI = 0.25        # [mm/mm] Steel I Section Ultimate Strain
    EshI = (fuI - fyI)/(esuI - eyI)
    BsI = EshI / EsI
    
    coreTag, coverTag, steelTag, steelITag = secTag + 100, secTag + 200, secTag + 300, secTag + 400
    if STEEL_TYPE == 'ELASTIC':
        ops.uniaxialMaterial('Elastic', steelITag, EsI) 
    if STEEL_TYPE == 'INELASTIC':
        pinchX = 0.8   # Pinching factor in X direction
        pinchY = 0.5   # Pinching factor in Y direction
        damage1 = 0.0  # Damage due to ductility
        damage2 = 0.0  # Damage due to energy
        beta = 0.1     # Stiffness degradation parameter
        # BOX SECION
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

    # -------------------- Fiber section -------------------------
    ops.section('Fiber', secTag)
    # -------------------- Concrete Section -------------------------------
    # (layers depth[mm], layers width[mm], center y‑coord [mm], center x‑coord [mm], numSubdivY, numSubdivZ, tag, color)

    PT = 10.0    # [mm] Plate Thickness
    CD = 200.0   # [mm] Section Depth
    NFY, NFX = 10, 10 # Number of fibers in each Division in Y and X Dir.
    mat_layers = [
        (PT, CD+2*PT, CD+0.5*PT, 0.5*CD, NFY, NFX, steelITag, 'lightgray'),              # 01 -> BOX
        (PT, CD+2*PT, 0.0-0.5*PT, 0.5*CD, NFY, NFX, steelITag, 'lightgray'),             # 02 -> BOX
        (CD, PT, 0.5*CD, CD+0.5*PT, NFY, NFX, steelITag, 'lightgray'),                   # 03 -> BOX   
        (CD, PT, 0.5*CD, -0.5*PT, NFY, NFX, steelITag, 'lightgray'),                     # 04 -> BOX
    ]
    
    for depth, width, center_y, center_x, numSubdivY, numSubdivZ, matTAG, _ in mat_layers:
        x_left = center_x - width/2
        x_right = center_x + width/2
        y_bot = center_y - depth/2
        y_top = center_y + depth/2
        ops.patch('rect', matTAG, numSubdivY, numSubdivZ, x_left, y_bot, x_right, y_top)

    # -------------------- Section area -------------------------
    # Initialize
    min_bottom = float('inf')
    max_top = -float('inf')
    TOTAL_AREA = 0.0
    
    for depth, width, center_y, center_x,_ ,_,_, _ in mat_layers:
        bottom = center_y - depth/2
        top = center_y + depth/2
        min_bottom = min(min_bottom, bottom)
        max_top = max(max_top, top)
        TOTAL_AREA += depth * width
    
    SECTION_HEIGHT = max_top - min_bottom
    
    print(f"Section height = {SECTION_HEIGHT} mm")
    print(f"Section area   = {TOTAL_AREA} mm²")

    # -------------------- Plot (optional) -----------------------
    if plot:
        
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.set_xlabel('Width (mm)')
        ax.set_ylabel('Height (mm)')
        ax.set_title(f'Steel Box section with plates (Numbers shown) - Section Tag: {secTag}')
        ax.grid(True, ls='--', alpha=0.5)
        
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
                color='purple', fontsize=8,
                ha='center', va='bottom', fontweight='bold'
            )
                
        # Compute overall limits from the actual steel layers
        all_x = []
        all_y = []
        for depth, width, center_y, center_x,_ ,_,_,_ in mat_layers:
            all_x.append(center_x - width/2)
            all_x.append(center_x + width/2)
            all_y.append(center_y - depth/2)
            all_y.append(center_y + depth/2)
        x_min, x_max = min(all_x), max(all_x)
        y_min, y_max = min(all_y), max(all_y)
        margin = 50  # mm extra space around
        ax.set_xlim(x_min - margin, x_max + margin)
        ax.set_ylim(y_min - margin, y_max + margin)
        ax.set_aspect('equal')
        plt.show()

    # Calculate mass
    ELE_MASS = STEEL_DENSITY * TOTAL_AREA   # kg/mm
    
    return SECTION_HEIGHT, ELE_MASS