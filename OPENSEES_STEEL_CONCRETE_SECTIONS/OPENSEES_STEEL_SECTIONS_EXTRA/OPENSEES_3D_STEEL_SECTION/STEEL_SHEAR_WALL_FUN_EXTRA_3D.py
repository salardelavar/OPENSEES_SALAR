def STEEL_SHEAR_WALL_FUN_EXTRA_3D(secTag, STEEL_TYPE, STEEL_DENSITY, plot=True):
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
    # STEEL
    fy = 240          # [N/mm²] Steel Yield Strength   
    Es = 2e5          # [N/mm²] Modulus of Elasticity
    ey = fy/Es        # [mm/mm] Steel Yield Strain
    fu = 1.1818*fy    # [N/mm²] Steel Ultimate Strength
    esu = 0.15        # [mm/mm] Steel Ultimate Strain
    Esh = (fu - fy)/(esu - ey)
    Bs = Esh / Es

    steelTag = secTag + 100
    if STEEL_TYPE == 'ELASTIC':
        ops.uniaxialMaterial('Elastic', steelTag, Es) 
    if STEEL_TYPE == 'INELASTIC':   
        pinchX = 0.8   # Pinching factor in X direction
        pinchY = 0.5   # Pinching factor in Y direction
        damage1 = 0.0  # Damage due to ductility
        damage2 = 0.0  # Damage due to energy
        beta = 0.1     # Stiffness degradation parameter
        ops.uniaxialMaterial('Hysteretic', steelTag,
                             fy, ey,
                             fu, esu,
                             0.2*fu, 1.1*esu,
                             -fy, -ey,
                             -fu, -esu,
                             -0.2*fu, -1.1*esu,
                             pinchX, pinchY,
                             damage1, damage2, beta)
        # INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Hysteretic_Material
        

    # -------------------- Fiber section -------------------------
    ops.section('Fiber', secTag)
    # -------------------- Steel Section -------------------------------
    # (layers depth[mm], layers width[mm], center y‑coord [mm], center x‑coord [mm], numSubdivY, numSubdivZ)
    steel_layers = [
        (10.0, 100.0, 5.0, 100.0, 10, 10),    # 1
        (180.0, 10.0, 100.0, 100.0, 10, 10),  # 2
        (10.0, 100.0, 195.0, 100.0, 10, 10),  # 3
        
        (180.0, 10.0, 290.0, 100.0, 10, 10),  # 4 
        (10.0, 100.0, 385.0, 100.0, 10, 10),  # 5
        
        (180.0, 10.0, 480.0, 100.0, 10, 10),  # 6 
        (10.0, 100.0, 575.0, 100.0, 10, 10),  # 7
        
        (180.0, 10.0, 670.0, 100.0, 10, 10),  # 8
        (10.0, 100.0, 765.0, 100.0, 10, 10),  # 9
        
        (180.0, 10.0, 860.0, 100.0, 10, 10),  # 10
        (10.0, 100.0, 955.0, 100.0, 10, 10),  # 11
        
        (180.0, 10.0, 1050.0, 100.0, 10, 10),  # 12
        (10.0, 100.0, 1145.0, 100.0, 10, 10),  # 13
        
        (180.0, 10.0, 1240.0, 100.0, 10, 10),  # 14
        (10.0, 100.0, 1335.0, 100.0, 10, 10),  # 15
        
        (180.0, 10.0, 1430.0, 100.0, 10, 10),  # 16 
        (10.0, 100.0, 1525.0, 100.0, 10, 10),  # 17
        
        (180.0, 10.0, 1620.0, 100.0, 10, 10),  # 18 
        (10.0, 100.0, 1715.0, 100.0, 10, 10),  # 19
        
        (180.0, 10.0, 1810.0, 100.0, 10, 10),  # 20
        (10.0, 100.0, 1905.0, 100.0, 10, 10),  # 21
        
        (180.0, 10.0, 2000.0, 100.0, 10, 10),  # 22 
        (10.0, 100.0, 2095.0, 100.0, 10, 10),  # 23
        
        (180.0, 10.0, 2190.0, 100.0, 10, 10),  # 24 
        (10.0, 100.0, 2285.0, 100.0, 10, 10),  # 25
        
        (180.0, 10.0, 2380.0, 100.0, 10, 10),  # 26 
        (10.0, 100.0, 2475.0, 100.0, 10, 10),  # 27
        
        (180.0, 10.0, 2570.0, 100.0, 10, 10),  # 28 
        (10.0, 100.0, 2665.0, 100.0, 10, 10),  # 29
        
        (180.0, 10.0, 2760.0, 100.0, 10, 10),  # 30 
        (10.0, 100.0, 2855.0, 100.0, 10, 10),  # 31
        
        (180.0, 10.0, 2950.0, 100.0, 10, 10),  # 32 
        (10.0, 100.0, 3045.0, 100.0, 10, 10),  # 33
        
        (180.0, 10.0, 3140.0, 100.0, 10, 10),  # 34 
        (10.0, 100.0, 3235.0, 100.0, 10, 10),  # 35
        


    ]
    for depth, width, center_y, center_x, numSubdivY, numSubdivZ in steel_layers:
        x_left = center_x - width/2
        x_right = center_x + width/2
        y_bot = center_y - depth/2
        y_top = center_y + depth/2
        ops.patch('rect', steelTag, numSubdivY, numSubdivZ, x_left, y_bot, x_right, y_top)

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
    
    print(f"Section height = {SECTION_HEIGHT} mm")
    print(f"Section area   = {TOTAL_AREA} mm²")

    # -------------------- Plot (optional) -----------------------
    if plot:
        import matplotlib.pyplot as plt
        import matplotlib.patches as patches
        
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.set_xlabel('Width (mm)')
        ax.set_ylabel('Height (mm)')
        ax.set_title('Steel section (35 layers)')
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
                linewidth=1, edgecolor='black', facecolor='lightgray'
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
        for depth, width, center_y, center_x,_ ,_ in steel_layers:
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

    ELE_MASS = STEEL_DENSITY * TOTAL_AREA   # kg/mm

    return SECTION_HEIGHT, ELE_MASS