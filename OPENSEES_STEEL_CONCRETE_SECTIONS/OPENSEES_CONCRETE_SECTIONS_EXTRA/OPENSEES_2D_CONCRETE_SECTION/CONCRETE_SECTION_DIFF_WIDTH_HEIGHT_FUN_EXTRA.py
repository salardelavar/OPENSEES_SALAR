def CONCRETE_SECTION_DIFF_WIDTH_HEIGHT_FUN_EXTRA(secTag, STEEL_TYPE, fc, Kfc, CONCRETE_DENSITY, plot=True):
    """
    Create a concrete with different widths and heights, fiber section (OpenSees) and,
    if PLOT=True, draw a simple 2‑D sketch of the section.

    Parameters
    ----------
    secTag          : int   – section identifier
    STEEL_TYPE      : str   – 'ELASTIC' or 'INELASTIC'
    fc              : float – unconfined concrete compressive strength (MPa, negative)
    Kfc             : float – ratio of confined to unconfined concrete strength
    CONCRETE_DENSITY: float – ρc (kg/m³) – will be converted to N·s²/mm³
    PLOT            : bool  – draw the section if True
    
    THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
    """
    import numpy as np
    import openseespy.opensees as ops
    import matplotlib.pyplot as plt
    import matplotlib.pyplot as plt, matplotlib.patches as patches
    
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
    

    coreTag, coverTag, steelTag = secTag + 100, secTag + 200, secTag + 300
    if STEEL_TYPE == 'ELASTIC':
        #ops.uniaxialMaterial('Steel01', steelTag, fy, Es, Bs) 
        ops.uniaxialMaterial('Elastic', steelTag, Es) 
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
    ops.uniaxialMaterial('Concrete02', coreTag, fcC, ec0C, fcUC, ecuC, LambdaC, ftC, EtsC) # build core concrete (confined)
    ops.uniaxialMaterial('Concrete02', coverTag, fcU, ec0U, fcUU, ecuU, LambdaU, ftU, EtsU) # build cover concrete (unconfined)


    # -------------------- Fiber section -------------------------
    ops.section('Fiber', secTag)
    # -------------------- Concrete Section -------------------------------
    # (layers depth[mm], layers width[mm], center y‑coord [mm], center x‑coord [mm], numSubdivY, numSubdivZ)
    concrete_layers = [
        (50.0, 700.0, 25.0, 350.0, 10, 10),    # 1
        (50.0, 680.0, 75.0, 350.0, 10, 10),    # 2
        (50.0, 660.0, 125.0, 350.0, 10, 10),   # 3
        (50.0, 640.0, 175.0, 350.0, 10, 10),   # 4
        (50.0, 620.0, 225.0, 350.0, 10, 10),   # 5
        (50.0, 600.0, 275.0, 350.0, 10, 10),   # 6
        (50.0, 580.0, 325.0, 350.0, 10, 10),   # 7
        (50.0, 560.0, 375.0, 350.0, 10, 10),   # 8
        (50.0, 540.0, 425.0, 350.0, 10, 10),   # 9
        (50.0, 520.0, 475.0, 350.0, 10, 10),   # 10
        (50.0, 500.0, 525.0, 350.0, 10, 10),   # 11
        (50.0, 520.0, 575.0, 350.0, 10, 10),   # 12 
        (50.0, 540.0, 625.0, 350.0, 10, 10),   # 13
        (50.0, 560.0, 675.0, 350.0, 10, 10),   # 14
        (50.0, 580.0, 725.0, 350.0, 10, 10),   # 15
        (50.0, 620.0, 775.0, 350.0, 10, 10),   # 16
        (50.0, 640.0, 825.0, 350.0, 10, 10),   # 17
        (50.0, 660.0, 875.0, 350.0, 10, 10),   # 18
        (50.0, 680.0, 925.0, 350.0, 10, 10),   # 19
        (50.0, 700.0, 975.0, 350.0, 10, 10),   # 20       

    ]
    for depth, width, center_y, center_x, numSubdivY, numSubdivZ in concrete_layers:
        x_left = center_x - width/2
        x_right = center_x + width/2
        y_bot = center_y - depth/2
        y_top = center_y + depth/2
        ops.patch('rect', coverTag, numSubdivY, numSubdivZ, x_left, y_bot, x_right, y_top)

    # -------------------- Section area -------------------------
    # Initialize
    min_bottom = float('inf')
    max_top = -float('inf')
    TOTAL_AREA = 0.0
    
    for depth, width, center_y, center_x,_ ,_ in concrete_layers:
        bottom = center_y - depth/2
        top = center_y + depth/2
        min_bottom = min(min_bottom, bottom)
        max_top = max(max_top, top)
        TOTAL_AREA += depth * width
    
    SECTION_HEIGHT = max_top - min_bottom
    
    print(f"Section height = {SECTION_HEIGHT} mm")
    print(f"Section area   = {TOTAL_AREA} mm²")
    # -------------------- Rebars -------------------------------
    # (diameter [mm], y‑coord [mm], x‑coord [mm])
    rebars = [
        (25.0,  50.0, 50.0),    # 1
        (25.0,  50.0, 350.0),   # 2
        (25.0,  50.0, 650.0),   # 3
        (18.0,  100.0, 70.0),   # 4
        (18.0,  100.0, 350.0),  # 5
        (18.0,  100.0, 630.0),  # 6
        (18.0,  150.0, 90.0),   # 7
        (18.0,  150.0, 350.0),  # 8
        (18.0,  150.0, 610.0),  # 9
        (18.0,  200.0, 110.0),  # 10
        (18.0,  200.0, 350.0),  # 11
        (18.0,  200.0, 590.0),  # 12
        (25.0,  950.0, 50.0),   # 13
        (25.0,  950.0, 350.0),  # 14
        (25.0,  950.0, 650.0),  # 15
        (18.0,  900.0, 70.0),   # 16
        (18.0,  900.0, 350.0),  # 17
        (18.0,  900.0, 630.0),  # 18 
        (18.0,  850.0, 90.0),   # 19
        (18.0,  850.0, 350.0),  # 20
        (18.0,  850.0, 610.0),  # 21 
        (18.0,  800.0, 110.0),  # 22
        (18.0,  800.0, 350.0),  # 23
        (18.0,  800.0, 590.0),  # 24
        
        (16.0,  500.0, 150.0),   # 25
        (16.0,  500.0, 350.0),   # 26
        (16.0,  500.0, 550.0),   # 27
    ]

    for dia, y, x in rebars:
        area = np.pi * dia**2 / 4.0          # mm²
        ops.fiber(x, y, area, steelTag)      # add to OpenSees model

    # -------------------- Plot (optional) -----------------------
    if plot:
        import matplotlib.pyplot as plt
        import matplotlib.patches as patches
        
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.set_xlabel('Width (mm)')
        ax.set_ylabel('Height (mm)')
        ax.set_title('Concrete section (20 layers) with rebars')
        ax.grid(True, ls='--', alpha=0.5)
        
        # Concrete geometry derived from concrete_layers
        ii = 0
        for depth, width, center_y, center_x, _, _ in concrete_layers:
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
                color='green', fontsize=8,
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
        
        # Compute overall limits from the actual concrete layers
        all_x = []
        all_y = []
        for depth, width, center_y, center_x,_ ,_ in concrete_layers:
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

    ELE_MASS = CONCRETE_DENSITY * TOTAL_AREA   # kg/mm

    return SECTION_HEIGHT, ELE_MASS