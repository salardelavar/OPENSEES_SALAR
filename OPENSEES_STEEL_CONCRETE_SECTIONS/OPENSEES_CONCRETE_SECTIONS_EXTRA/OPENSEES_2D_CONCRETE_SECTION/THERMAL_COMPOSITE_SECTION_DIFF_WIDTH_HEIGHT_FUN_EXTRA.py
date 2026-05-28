def THERMAL_COMPOSITE_SECTION_DIFF_WIDTH_HEIGHT_FUN_EXTRA(secTag, fc, Kfc, CONCRETE_DENSITY, plot=True):
    """
    Create a composite with different widths and heights, fiber thermal section (OpenSees) and,
    if PLOT=True, draw a simple 2‑D sketch of the section.

    Parameters
    ----------
    secTag          : int   – section identifier
    STEEL_TYPE      : str   – 'ELASTIC' or 'INELASTIC'
    fc              : float – unconfined concrete compressive strength (MPa, negative)
    Kfc             : float – ratio of confined to unconfined concrete strength
    CONCRETE_DENSITY: float – ρc (kg/m³) – will be converted to N·s²/mm³
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
    
    coreTag, coverTag, steelTag, steelITag = secTag + 100, secTag + 200, secTag + 300, secTag + 400
    ops.uniaxialMaterial('Steel01Thermal', steelTag, fy, Es, Bs)     # REBARS
    ops.uniaxialMaterial('Steel01Thermal', steelITag, fyI, EsI, BsI) # PLATES

    #ops.uniaxialMaterial('Concrete01Thermal', coverTag, fcU, ec0U, fcUU, ecuU)    
    #ops.uniaxialMaterial('Concrete01Thermal', coreTag, fcC, ec0C, fcUC, ecuC)     
    ops.uniaxialMaterial('Concrete02Thermal', coverTag, fcU, ec0U, fcUU, ecuU, LambdaU, ftU, EtsU)
    ops.uniaxialMaterial('Concrete02Thermal', coreTag, fcC, ec0C, fcUC, ecuC, LambdaC, ftC, EtsC)

    # -------------------- Fiber section -------------------------
    ops.section('FiberThermal', secTag)
    # -------------------- Concrete Section -------------------------------
    PT = 10.0    # [mm] Plate Thickness
    # (fiber area[mm²], center y‑coord [mm], center x‑coord [mm], tag, color)

    mat_layers = [
        (100.0, 50.0,   50.0, coreTag, 'lightgray'),   # 1
        (100.0, 50.0,  150.0, coreTag, 'lightgray'),   # 2
        (100.0, 50.0,  250.0, coreTag, 'lightgray'),   # 3
        (100.0, 50.0,  350.0, coreTag, 'lightgray'),   # 4
        (100.0, 50.0,  450.0, coreTag, 'lightgray'),   # 5
        (100.0, 150.0,  50.0, coreTag, 'lightgray'),   # 6
        (100.0, 150.0, 150.0, coreTag, 'lightgray'),   # 7
        (100.0, 150.0, 250.0, coreTag, 'lightgray'),   # 8
        (100.0, 150.0, 350.0, coreTag, 'lightgray'),   # 9
        (100.0, 150.0, 450.0, coreTag, 'lightgray'),   # 10
        (100.0, 250.0,  50.0, coreTag, 'lightgray'),   # 11
        (100.0, 250.0, 150.0, coreTag, 'lightgray'),   # 12
        (100.0, 250.0, 250.0, coreTag, 'lightgray'),   # 13
        (100.0, 250.0, 350.0, coreTag, 'lightgray'),   # 14
        (100.0, 250.0, 450.0, coreTag, 'lightgray'),   # 15
        (100.0, 350.0,  50.0, coreTag, 'lightgray'),   # 16
        (100.0, 350.0, 150.0, coreTag, 'lightgray'),   # 17
        (100.0, 350.0, 250.0, coreTag, 'lightgray'),   # 18
        (100.0, 350.0, 350.0, coreTag, 'lightgray'),   # 19
        (100.0, 350.0, 450.0, coreTag, 'lightgray'),   # 20
        (100.0, 450.0,  50.0, coreTag, 'lightgray'),   # 21
        (100.0, 450.0, 150.0, coreTag, 'lightgray'),   # 22
        (100.0, 450.0, 250.0, coreTag, 'lightgray'),   # 23
        (100.0, 450.0, 350.0, coreTag, 'lightgray'),   # 24
        (100.0, 450.0, 450.0, coreTag, 'lightgray'),   # 25
    
        (PT * (500.0+2*PT), 500.0+0.5*PT, 0.5*500.0, steelITag, 'cyan'),   # 26 -> PLATE
        (PT * (500.0+2*PT), 0.0-0.5*PT, 0.5*500.0, steelITag, 'cyan'),     # 27 -> PLATE
        (500.0 * PT, 0.5*500.0, 500.0+0.5*PT, steelITag, 'cyan'),        # 28 -> PLATE   
        (500.0 * PT, 0.5*500.0, -0.5*PT, steelITag, 'cyan'),             # 29 -> PLATE

        (10.0 * 150.0, 250.0-125.0, 250.0, steelITag, 'lime'),                  # 30 -> CROSSED SECTION
        (240.0 * 10.0, 250.0, 250.0, steelITag, 'lime'),                        # 31 -> CROSSED SECTION
        (10.0 * 150.0, 250+125.0, 250.0, steelITag, 'lime'),                    # 32 -> CROSSED SECTION   
         
        (10.0 * 105.0, 250.0, 250.0-57.5, steelITag, 'lime'),                   # 33 -> CROSSED SECTION
        (10.0 * 105.0, 250.0, 250.0+57.5, steelITag, 'lime'),                   # 34 -> CROSSED SECTION
        
        (150.0 * 10.0, 250.0, 250.0-115.0, steelITag, 'lime'),                  # 35 -> CROSSED SECTION
        (150.0 * 10.0, 250.0, 250.0+115.0, steelITag, 'lime'),                  # 36 -> CROSSED SECTION     
         ]
    """
    mat_layers = [
        (20.0 * 500.0, 10.0, 250.0, coreTag, 'lightgray'),    # 1
        (20.0 * 500.0, 30.0, 250.0, coreTag, 'lightgray'),    # 2
        (20.0 * 500.0, 50.0, 250.0, coreTag, 'lightgray'),    # 3
        (20.0 * 500.0, 70.0, 250.0, coreTag, 'lightgray'),    # 4
        (20.0 * 500.0, 90.0, 250.0, coreTag, 'lightgray'),    # 5
        (20.0 * 500.0, 110.0, 250.0, coreTag, 'lightgray'),   # 6
        (20.0 * 500.0, 130.0, 250.0, coreTag, 'lightgray'),   # 7
        (20.0 * 500.0, 150.0, 250.0, coreTag, 'lightgray'),   # 8
        (20.0 * 500.0, 170.0, 250.0, coreTag, 'lightgray'),   # 9
        (20.0 * 500.0, 190.0, 250.0, coreTag, 'lightgray'),   # 10
        (20.0 * 500.0, 210.0, 250.0, coreTag, 'lightgray'),   # 11
        (20.0 * 500.0, 230.0, 250.0, coreTag, 'lightgray'),   # 12 
        (20.0 * 500.0, 250.0, 250.0, coreTag, 'lightgray'),   # 13
        (20.0 * 500.0, 270.0, 250.0, coreTag, 'lightgray'),   # 14
        (20.0 * 500.0, 290.0, 250.0, coreTag, 'lightgray'),   # 15
        (20.0 * 500.0, 310.0, 250.0, coreTag, 'lightgray'),   # 16
        (20.0 * 500.0, 330.0, 250.0, coreTag, 'lightgray'),   # 17
        (20.0 * 500.0, 350.0, 250.0, coreTag, 'lightgray'),   # 18
        (20.0 * 500.0, 370.0, 250.0, coreTag, 'lightgray'),   # 19
        (20.0 * 500.0, 390.0, 250.0, coreTag, 'lightgray'),   # 20
        (20.0 * 500.0, 410.0, 250.0, coreTag, 'lightgray'),   # 21
        (20.0 * 500.0, 430.0, 250.0, coreTag, 'lightgray'),   # 22
        (20.0 * 500.0, 450.0, 250.0, coreTag, 'lightgray'),   # 23
        (20.0 * 500.0, 470.0, 250.0, coreTag, 'lightgray'),   # 24
        (20.0 * 500.0, 490.0, 250.0, coreTag, 'lightgray'),   # 25
        
        (PT * (500.0+2*PT), 500.0+0.5*PT, 0.5*500.0, steelITag, 'cyan'),   # 26 -> PLATE
        (PT * (500.0+2*PT), 0.0-0.5*PT, 0.5*500.0, steelITag, 'cyan'),     # 27 -> PLATE
        (500.0 * PT, 0.5*500.0, 500.0+0.5*PT, steelITag, 'cyan'),          # 28 -> PLATE   
        (500.0 * PT, 0.5*500.0, -0.5*PT, steelITag, 'cyan'),               # 29 -> PLATE

        (10.0 * 150.0, 250.0-125.0, 250.0, steelITag, 'lime'),                  # 30 -> CROSSED SECTION
        (240.0 * 10.0, 250.0, 250.0, steelITag, 'lime'),                        # 31 -> CROSSED SECTION
        (10.0 * 150.0, 250+125.0, 250.0, steelITag, 'lime'),                    # 32 -> CROSSED SECTION   
         
        (10.0 * 105.0, 250.0, 250.0-57.5, steelITag, 'lime'),                   # 33 -> CROSSED SECTION
        (10.0 * 105.0, 250.0, 250.0+57.5, steelITag, 'lime'),                   # 34 -> CROSSED SECTION
        
        (150.0 * 10.0, 250.0, 250.0-115.0, steelITag, 'lime'),                  # 35 -> CROSSED SECTION
        (150.0 * 10.0, 250.0, 250.0+115.0, steelITag, 'lime'),                  # 36 -> CROSSED SECTION   
        ]
    """
    for fiber_area, center_y, center_x, matTAG, _ in mat_layers:
        ops.fiber(center_x, center_y, fiber_area, matTAG)

    # -------------------- Section area -------------------------
    # Initialize
    min_bottom = float('inf')
    max_top = -float('inf')
    SECTION_AREA = 0.0
    
    for area, y, x, *_ in mat_layers:          # unpack area, center_y, center_x; ignore the rest
        # Assume square cross‑section → side = sqrt(area)
        side = np.sqrt(area)
        half_side = side / 2.0
    
        bottom = y - half_side
        top = y + half_side
    
        min_bottom = min(min_bottom, bottom)
        max_top = max(max_top, top)
    
        SECTION_AREA += area
    
    SECTION_HEIGHT = max_top - min_bottom
    
    print(f"SECTION_AREA   = {SECTION_AREA:.2f} mm²")
    print(f"SECTION_HEIGHT = {SECTION_HEIGHT:.2f} mm")
    # -------------------- Rebars -------------------------------
    # (diameter [mm], y‑coord [mm], x‑coord [mm])
    rebars = [
        (25.0,  50.0, 50.0),    # 1
        (26.0,  50.0, 250.0),   # 2
        (25.0,  50.0, 450.0),   # 3 
        (25.0,  450.0, 50.0),   # 4
        (16.0,  450.0, 250.0),  # 5
        (25.0,  450.0, 450.0),  # 6
        
        (25.0,  250.0, 50.0),   # 7
        (25.0,  250.0, 450.0),  # 8
        
    ]

    for dia, y, x in rebars:
        area = np.pi * dia**2 / 4.0          # mm²
        ops.fiber(x, y, area, steelTag)      # add to OpenSees model

    # -------------------- Plot (optional) -----------------------
    if plot:
        
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.set_xlabel('Width (mm)')
        ax.set_ylabel('Height (mm)')
        ax.set_title(f'Thermal Composite section with rebars and plates (Numbers shown) - Section Tag: {secTag}')
        ax.grid(True, ls='--', alpha=0.5)
    
        # Material geometry derived from mat_layers (square assumption)
        ii = 0
        for area, center_y, center_x, tag, color in mat_layers:
            side = np.sqrt(area)            # assume square
            width = depth = side            # both dimensions equal
    
            x_left  = center_x - width / 2
            x_right = center_x + width / 2
            y_bot   = center_y - depth / 2
            y_top   = center_y + depth / 2
    
            rect = patches.Rectangle(
                (x_left, y_bot), width, depth,
                linewidth=1, edgecolor='black', facecolor=color
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
    
        # Compute overall limits from the actual concrete layers (square assumption)
        all_x = []
        all_y = []
        for area, center_y, center_x, *_ in mat_layers:
            side = np.sqrt(area)
            half = side / 2.0
            all_x.append(center_x - half)
            all_x.append(center_x + half)
            all_y.append(center_y - half)
            all_y.append(center_y + half)
    
        x_min, x_max = min(all_x), max(all_x)
        y_min, y_max = min(all_y), max(all_y)
        margin = 50  # mm extra space around
        ax.set_xlim(x_min - margin, x_max + margin)
        ax.set_ylim(y_min - margin, y_max + margin)
        ax.set_aspect('equal')
        plt.show()

    # Calculate mass
    ELE_MASS = CONCRETE_DENSITY * SECTION_AREA   # kg/mm
    
    return SECTION_HEIGHT, ELE_MASS