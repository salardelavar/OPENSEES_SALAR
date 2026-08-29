def AFTER_FIRE_04_I_STEEL_PLATE_SECTION_DIFF_WIDTH_HEIGHT_FUN_EXTRA_3D(secTag, STEEL_TYPE, STEEL_DENSITY, plot=True):
    """
    Create a steel with different widths and heights, fiber section (OpenSees) and,
    if PLOT=True, draw a simple 2‑D sketch of the section.

    Parameters
    ----------
    secTag          : int   – section identifier
    STEEL_TYPE      : str   – 'ELASTIC' or 'INELASTIC'
    STEEL_DENSITY: float – ρc (kg/m³) – will be converted to N·s²/mm³
    PLOT            : bool  – draw the section if True
    
    Thermal Effects of Material Stres-Strain relation
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

    def STEEL_FIRE_PROPERTIES(T):
        # THERMAL EFFECTS OF STEEL MATERIAL STRESS-STRAIN
        # WRITTEN BY SALAR DElAVAR GHASHGHAEI (QASHQAI)
        # PAPER: Mechanical properties of structural steel at elevated temperatures and after cooling down
        'https://www.researchgate.net/publication/227780727_Mechanical_properties_of_structural_steel_at_elevated_temperatures_and_after_cooling_down'    
        if T <= 20:
            kE  = 1.00
            kFy = 1.00
            kFu = 1.00
            esu = 0.25
            
        elif 20 < T <= 100:
            kE  = 0.90
            kFy = 0.90
            kFu = 0.90
            esu = 0.277
            
        elif 100 < T <= 200:
            kE  = 0.80
            kFy = 0.80
            kFu = 0.80
            esu = 0.3125
            
        elif 200 < T <= 300:
            kE  = 0.70
            kFy = 0.70
            kFu = 0.70
            esu = 0.3571
            
        elif 300 < T <= 400:
            kE  = 0.60
            kFy = 0.60
            kFu = 0.60
            esu = 0.3814 
            
        elif 400 < T <= 500:
            kE  = 0.30
            kFy = 0.30
            kFu = 0.30
            esu = 0.4012              
        
        elif 500 < T <= 700:
            kE  = 0.13
            kFy = 0.23
            kFu = 0.35
            esu = 0.4124
    
        elif T > 700:
            kE  = 0.09
            kFy = 0.11
            kFu = 0.20
            esu = 0.45
    
        return kE, kFy, kFu, esu

    # Plate Section 01
    kE, kFy, kFu, esu = STEEL_FIRE_PROPERTIES(100) # THERMAL EFFECTS 100°C 
    EsI01 = kE*200000        # [N/mm²] I Section Modulus of Elasticity 
    fyI01 = kFy*240          # [N/mm²] Steel I Section Yield Strength
    fuI01 = kFu*(1.18*240)   # [N/mm²] Steel I Section Ultimate Strength
    #fyI01 = 240                
    #EsI01 = 2e5              
    eyI01 = fyI01/EsI01      # [mm/mm] Steel I Section Yield Strain
    #fuI01 = 1.1818*fyI01     # [N/mm²] Steel I Section Ultimate Strength
    esuI01 = esu            # [mm/mm] Steel I Section Ultimate Strain
    EshI01 = (fuI01 - fyI01)/(esuI01 - eyI01)
    BsI01 = EshI01 / EsI01
    
    # Plate Section 02
    kE, kFy, kFu, esu = STEEL_FIRE_PROPERTIES(300) # THERMAL EFFECTS 300°C
    EsI02 = kE*200000        # [N/mm²] I Section Modulus of Elasticity 
    fyI02 = kFy*240          # [N/mm²] Steel I Section Yield Strength
    fuI02 = kFu*(1.18*240)   # [N/mm²] Steel I Section Ultimate Strength
    #fyI02 = 140              # [N/mm²] Steel I Section Yield Strength   
    #EsI02 = 2e5              # [N/mm²] I Section Modulus of Elasticity
    eyI02 = fyI02/EsI02      # [mm/mm] Steel I Section Yield Strain
    #fuI02 = 1.1818*fyI02     # [N/mm²] Steel I Section Ultimate Strength
    esuI02 = esu            # [mm/mm] Steel I Section Ultimate Strain
    EshI02 = (fuI02 - fyI02)/(esuI02 - eyI02)
    BsI02 = EshI02 / EsI02
    
    # Plate Section 03
    kE, kFy, kFu, esu = STEEL_FIRE_PROPERTIES(500) # THERMAL EFFECTS 500°C
    EsI03 = kE*200000        # [N/mm²] I Section Modulus of Elasticity 
    fyI03 = kFy*240          # [N/mm²] Steel I Section Yield Strength
    fuI03 = kFu*(1.18*240)   # [N/mm²] Steel I Section Ultimate Strength
    #EsI03 = kE*200000        # [N/mm²] I Section Modulus of Elasticity 
    #fyI03 = kFy*240          # [N/mm²] Steel I Section Yield Strength
    #fuI3 = kFu*(1.18*240)   # [N/mm²] Steel I Section Ultimate Strength
    #fyI03 = 95              # [N/mm²] Steel I Section Yield Strength   
    #EsI03 = 2e5              # [N/mm²] I Section Modulus of Elasticity
    eyI03 = fyI03/EsI03      # [mm/mm] Steel I Section Yield Strain
    #fuI03 = 1.1818*fyI03     # [N/mm²] Steel I Section Ultimate Strength
    esuI03 = esu            # [mm/mm] Steel I Section Ultimate Strain
    EshI03 = (fuI03 - fyI03)/(esuI03 - eyI03)
    BsI03 = EshI03 / EsI03
    
    steelITag01, steelITag02, steelITag03 = secTag + 400, secTag + 500, secTag + 600
    if STEEL_TYPE == 'ELASTIC':
        ops.uniaxialMaterial('Elastic', steelITag01, EsI01) 
        ops.uniaxialMaterial('Elastic', steelITag02, EsI02) 
        ops.uniaxialMaterial('Elastic', steelITag03, EsI03) 
    if STEEL_TYPE == 'INELASTIC':   
        pinchX = 0.8   # Pinching factor in X direction
        pinchY = 0.5   # Pinching factor in Y direction
        damage1 = 0.0  # Damage due to ductility
        damage2 = 0.0  # Damage due to energy
        beta = 0.1     # Stiffness degradation parameter
        # INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Hysteretic_Material
        # PLATE SECION 01
        ops.uniaxialMaterial('Hysteretic', steelITag01,
                                        fyI01, eyI01,
                                        fuI01, esuI01,
                                        0.2*fuI01, 1.1*esuI01,
                                        -fyI01, -eyI01,
                                        -fuI01, -esuI01,
                                        -0.2*fuI01, -1.1*esuI01,
                                        pinchX, pinchY,
                                        damage1, damage2, beta)
        # PLATE SECION 02
        ops.uniaxialMaterial('Hysteretic', steelITag02,
                                        fyI02, eyI02,
                                        fuI02, esuI02,
                                        0.2*fuI02, 1.1*esuI02,
                                        -fyI02, -eyI02,
                                        -fuI02, -esuI02,
                                        -0.2*fuI02, -1.1*esuI02,
                                        pinchX, pinchY,
                                        damage1, damage2, beta)
        # PLATE SECION 03
        ops.uniaxialMaterial('Hysteretic', steelITag03,
                                        fyI03, eyI03,
                                        fuI03, esuI03,
                                        0.2*fuI03, 1.1*esuI03,
                                        -fyI03, -eyI03,
                                        -fuI03, -esuI03,
                                        -0.2*fuI03, -1.1*esuI03,
                                        pinchX, pinchY,
                                        damage1, damage2, beta)                                        
        # INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Hysteretic_Material
        
    # -------------------- Fiber section -------------------------
    ops.section('Fiber', secTag, '-GJ', 1.0e6)
    # -------------------- Steel Section -------------------------------
    # (layers depth[mm], layers width[mm], center y‑coord [mm], center x‑coord [mm], numSubdivY, numSubdivZ)
    NFY, NFX = 10, 10 # Number of fibers in each Division in Y and X Dir.
    steel_layers = [
        (10.0, 100.0, 5.0, 100.0, NFY, NFX, steelITag03, 'pink'),    # 1
        (160.0, 10.0, 90.0, 100.0, NFY, NFX, steelITag01, 'lime'),   # 2
        (10.0, 100.0, 175.0, 100.0, NFY, NFX, steelITag01, 'lime'),  # 3    
         
        (10.0, 75.0, -5.0, 100.0, NFY, NFX, steelITag03, 'pink'),    # 4
        (10.0, 75.0, 185.0, 100.0, NFY, NFX, steelITag01, 'lime'),   # 5
        
        (100.0, 10.0, 90.0, 90.0, NFY, NFX, steelITag02, 'cyan'),    # 6
        (100.0, 10.0, 90.0, 110.0, NFY, NFX, steelITag01, 'lime'),   # 7
        
    ]
    for depth, width, center_y, center_x, numSubdivY, numSubdivZ, matTAG, _ in steel_layers:
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
    
    for depth, width, center_y, center_x, _, _,_, COLOR in steel_layers:
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
        ax.set_title('After Fire Effects Steel Section (7 layers)')
        ax.grid(True, ls='--', alpha=0.5)
        
        # Steel geometry derived from steel_layers
        ii = 0
        for depth, width, center_y, center_x, _, _,_, COLOR in steel_layers:
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
        for depth, width, center_y, center_x,_, _,_, COLOR in steel_layers:
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