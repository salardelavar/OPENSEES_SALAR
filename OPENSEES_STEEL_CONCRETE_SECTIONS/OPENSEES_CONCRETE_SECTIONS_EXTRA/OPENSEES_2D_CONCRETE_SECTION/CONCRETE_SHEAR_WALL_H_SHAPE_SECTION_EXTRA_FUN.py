def CONCRETE_SHEAR_WALL_H_SHAPE_SECTION_EXTRA_FUN(secTag, STEEL_TYPE, fc, Kfc, CONCRETE_DENSITY, plot=True):
    """
    Create a concrete H-Shape shear wall with different widths and heights, fiber section (OpenSees) and,
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
    

    coreTag, coverTag, steelTag = secTag + 100, secTag + 200, secTag + 300
    if STEEL_TYPE == 'ELASTIC':
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
    NFY, NFX = 50, 50 # Number of fibers in each Division in Y and X Dir.
    concrete_layers = [
        (4000.0, 300.0, 2000.0, 150.0, NFY, NFX, coreTag, 'lightgray'),            # 1
        (300.0, 5000.0, 2000.0, 300.0+2500.0, NFY, NFX, coreTag, 'lightgray'),    # 2
        (4000.0, 300.0, 2000.0, 5000.0+450.0, NFY, NFX, coreTag, 'lightgray'),     # 3

    ]
    for depth, width, center_y, center_x, numSubdivY, numSubdivZ, matTAG, _  in concrete_layers:
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
    
    for depth, width, center_y, center_x,_ ,_,_ ,_ in concrete_layers:
        bottom = center_y - depth/2
        top = center_y + depth/2
        min_bottom = min(min_bottom, bottom)
        max_top = max(max_top, top)
        TOTAL_AREA += depth * width
    
    SECTION_HEIGHT = max_top - min_bottom
    
    print(f"Section height = {SECTION_HEIGHT} mm")
    print(f"Section area   = {TOTAL_AREA} mm²")
    # -------------------- Rebars -------------------------------
    RD = 16.0 # [mm] Rebar Diameter
    # (diameter [mm], y‑coord [mm], x‑coord [mm])
    rebars = [
        (RD,  50.0, 50.0),      # 1 -> SHEAR WALL 01 - DISTANCE 150.0 mm
        (RD,  50.0, 250.0),     # 2
        (RD,  3950.0, 50.0),    # 3
        (RD,  3950.0, 250.0),   # 4
        (RD,  200.0, 50.0),     # 5
        (RD,  200.0, 250.0),    # 6
        (RD,  3800.0, 50.0),    # 7
        (RD,  3800.0, 250.0),   # 8
        (RD,  350.0, 50.0),     # 9
        (RD,  350.0, 250.0),    # 10
        (RD,  3650.0, 50.0),    # 11
        (RD,  3650.0, 250.0),   # 12
        (RD,  500.0, 50.0),     # 13
        (RD,  500.0, 250.0),    # 14
        (RD,  3500.0, 50.0),    # 15
        (RD,  3500.0, 250.0),   # 16
        (RD,  650.0, 50.0),     # 17
        (RD,  650.0, 250.0),    # 18
        (RD,  3350.0, 50.0),    # 19
        (RD,  3350.0, 250.0),   # 20
        (RD,  800.0, 50.0),     # 21
        (RD,  800.0, 250.0),    # 22
        (RD,  3200.0, 50.0),    # 23
        (RD,  3200.0, 250.0),   # 24
        (RD,  950.0, 50.0),     # 25
        (RD,  950.0, 250.0),    # 26
        (RD,  3050.0, 50.0),    # 27
        (RD,  3050.0, 250.0),   # 28
        
        (RD,  1150.0, 50.0),    # 29 -> SHEAR WALL 01 - DISTANCE 200.0 mm
        (RD,  1150.0, 250.0),   # 30
        (RD,  2850.0, 50.0),    # 31
        (RD,  2850.0, 250.0),   # 32
        (RD,  1150.0, 50.0),    # 33
        (RD,  1150.0, 250.0),   # 34
        (RD,  2850.0, 50.0),    # 35
        (RD,  2850.0, 250.0),   # 36
        (RD,  1350.0, 50.0),    # 37
        (RD,  1350.0, 250.0),   # 38
        (RD,  2650.0, 50.0),    # 39
        (RD,  2650.0, 250.0),   # 40
        (RD,  1550.0, 50.0),    # 41
        (RD,  1550.0, 250.0),   # 42
        (RD,  2450.0, 50.0),    # 43
        (RD,  2450.0, 250.0),   # 44
        (RD,  1750.0, 50.0),    # 45
        (RD,  1750.0, 250.0),   # 46
        (RD,  2250.0, 50.0),    # 47
        (RD,  2250.0, 250.0),   # 48
        (RD,  1950.0, 50.0),    # 49
        (RD,  1950.0, 250.0),   # 50
        (RD,  2050.0, 50.0),    # 51
        (RD,  2050.0, 250.0),   # 52
        
        (RD,  2000.0-0.5*300.0+50.0, 300.0+100.0+0*150.0),      # 53 -> SHEAR WALL 02 - DISTANCE 150.0 mm
        (RD,  2000.0+0.5*300.0-50.0, 300.0+100.0+0*150.0),      # 54
        (RD,  2000.0-0.5*300.0+50.0, 300.0+100.0+1*150.0),      # 55
        (RD,  2000.0+0.5*300.0-50.0, 300.0+100.0+1*150.0),      # 56
        (RD,  2000.0-0.5*300.0+50.0, 300.0+100.0+2*150.0),      # 57
        (RD,  2000.0+0.5*300.0-50.0, 300.0+100.0+2*150.0),      # 58
        (RD,  2000.0-0.5*300.0+50.0, 300.0+100.0+3*150.0),      # 59
        (RD,  2000.0+0.5*300.0-50.0, 300.0+100.0+3*150.0),      # 60
        (RD,  2000.0-0.5*300.0+50.0, 300.0+100.0+4*150.0),      # 61
        (RD,  2000.0+0.5*300.0-50.0, 300.0+100.0+4*150.0),      # 62
        (RD,  2000.0-0.5*300.0+50.0, 300.0+100.0+5*150.0),      # 63
        (RD,  2000.0+0.5*300.0-50.0, 300.0+100.0+5*150.0),      # 64
        (RD,  2000.0-0.5*300.0+50.0, 300.0+100.0+6*150.0),      # 65
        (RD,  2000.0+0.5*300.0-50.0, 300.0+100.0+6*150.0),      # 66
        (RD,  2000.0-0.5*300.0+50.0, 300.0+100.0+7*150.0),      # 67
        (RD,  2000.0+0.5*300.0-50.0, 300.0+100.0+7*150.0),      # 68
        (RD,  2000.0-0.5*300.0+50.0, 300.0+100.0+8*150.0),      # 69
        (RD,  2000.0+0.5*300.0-50.0, 300.0+100.0+8*150.0),      # 70
        (RD,  2000.0-0.5*300.0+50.0, 300.0+100.0+9*150.0),      # 71
        (RD,  2000.0+0.5*300.0-50.0, 300.0+100.0+9*150.0),      # 72
        (RD,  2000.0-0.5*300.0+50.0, 300.0+100.0+10*150.0),      # 73
        (RD,  2000.0+0.5*300.0-50.0, 300.0+100.0+10*150.0),      # 74
        (RD,  2000.0-0.5*300.0+50.0, 300.0+100.0+11*150.0),      # 75
        (RD,  2000.0+0.5*300.0-50.0, 300.0+100.0+11*150.0),      # 76
        (RD,  2000.0-0.5*300.0+50.0, 300.0+100.0+12*150.0),      # 77
        (RD,  2000.0+0.5*300.0-50.0, 300.0+100.0+12*150.0),      # 78
        (RD,  2000.0-0.5*300.0+50.0, 300.0+100.0+13*150.0),      # 79
        (RD,  2000.0+0.5*300.0-50.0, 300.0+100.0+13*150.0),      # 80
        (RD,  2000.0-0.5*300.0+50.0, 300.0+100.0+14*150.0),      # 81
        (RD,  2000.0+0.5*300.0-50.0, 300.0+100.0+14*150.0),      # 82
        (RD,  2000.0-0.5*300.0+50.0, 300.0+100.0+15*150.0),      # 83
        (RD,  2000.0+0.5*300.0-50.0, 300.0+100.0+15*150.0),      # 84
        (RD,  2000.0-0.5*300.0+50.0, 300.0+100.0+16*150.0),      # 85
        (RD,  2000.0+0.5*300.0-50.0, 300.0+100.0+16*150.0),      # 86
        (RD,  2000.0-0.5*300.0+50.0, 300.0+100.0+17*150.0),      # 87
        (RD,  2000.0+0.5*300.0-50.0, 300.0+100.0+17*150.0),      # 88
        (RD,  2000.0-0.5*300.0+50.0, 300.0+100.0+18*150.0),      # 89
        (RD,  2000.0+0.5*300.0-50.0, 300.0+100.0+18*150.0),      # 90
        (RD,  2000.0-0.5*300.0+50.0, 300.0+100.0+19*150.0),      # 91
        (RD,  2000.0+0.5*300.0-50.0, 300.0+100.0+19*150.0),      # 92
        (RD,  2000.0-0.5*300.0+50.0, 300.0+100.0+20*150.0),      # 93
        (RD,  2000.0+0.5*300.0-50.0, 300.0+100.0+20*150.0),      # 94
        (RD,  2000.0-0.5*300.0+50.0, 300.0+100.0+21*150.0),      # 95
        (RD,  2000.0+0.5*300.0-50.0, 300.0+100.0+21*150.0),      # 96
        (RD,  2000.0-0.5*300.0+50.0, 300.0+100.0+22*150.0),      # 97
        (RD,  2000.0+0.5*300.0-50.0, 300.0+100.0+22*150.0),      # 98
        (RD,  2000.0-0.5*300.0+50.0, 300.0+100.0+23*150.0),      # 99
        (RD,  2000.0+0.5*300.0-50.0, 300.0+100.0+23*150.0),      # 100
        (RD,  2000.0-0.5*300.0+50.0, 300.0+100.0+24*150.0),      # 101
        (RD,  2000.0+0.5*300.0-50.0, 300.0+100.0+24*150.0),      # 102
        (RD,  2000.0-0.5*300.0+50.0, 300.0+100.0+25*150.0),      # 103
        (RD,  2000.0+0.5*300.0-50.0, 300.0+100.0+25*150.0),      # 104
        (RD,  2000.0-0.5*300.0+50.0, 300.0+100.0+26*150.0),      # 105
        (RD,  2000.0+0.5*300.0-50.0, 300.0+100.0+26*150.0),      # 106
        (RD,  2000.0-0.5*300.0+50.0, 300.0+100.0+27*150.0),      # 107
        (RD,  2000.0+0.5*300.0-50.0, 300.0+100.0+27*150.0),      # 108
        (RD,  2000.0-0.5*300.0+50.0, 300.0+100.0+28*150.0),      # 109
        (RD,  2000.0+0.5*300.0-50.0, 300.0+100.0+28*150.0),      # 110
        (RD,  2000.0-0.5*300.0+50.0, 300.0+100.0+29*150.0),      # 111
        (RD,  2000.0+0.5*300.0-50.0, 300.0+100.0+29*150.0),      # 112
        (RD,  2000.0-0.5*300.0+50.0, 300.0+100.0+30*150.0),      # 113
        (RD,  2000.0+0.5*300.0-50.0, 300.0+100.0+30*150.0),      # 114
        (RD,  2000.0-0.5*300.0+50.0, 300.0+100.0+31*150.0),      # 115
        (RD,  2000.0+0.5*300.0-50.0, 300.0+100.0+31*150.0),      # 116
        (RD,  2000.0-0.5*300.0+50.0, 300.0+100.0+32*150.0),      # 117
        (RD,  2000.0+0.5*300.0-50.0, 300.0+100.0+32*150.0),      # 118
        
        (RD,  50.0, 5350.0),      # 119 -> SHEAR WALL 03 - DISTANCE 150.0 mm
        (RD,  50.0, 5550.0),      # 120
        (RD,  3950.0, 5350.0),    # 121
        (RD,  3950.0, 5550.0),    # 122
        (RD,  200.0, 5350.0),     # 123
        (RD,  200.0, 5550.0),     # 124
        (RD,  3800.0, 5350.0),    # 125
        (RD,  3800.0, 5550.0),    # 126
        (RD,  350.0, 5350.0),     # 127
        (RD,  350.0, 5550.0),     # 128
        (RD,  3650.0, 5350.0),    # 129
        (RD,  3650.0, 5550.0),    # 130
        (RD,  500.0, 5350.0),     # 131
        (RD,  500.0, 5550.0),     # 132
        (RD,  3500.0, 5350.0),    # 133
        (RD,  3500.0, 5550.0),    # 134
        (RD,  650.0, 5350.0),     # 135
        (RD,  650.0, 5550.0),     # 136
        (RD,  3350.0, 5350.0),    # 137
        (RD,  3350.0, 5550.0),    # 138
        (RD,  800.0, 5350.0),     # 139
        (RD,  800.0, 5550.0),     # 140
        (RD,  3200.0, 5350.0),    # 141
        (RD,  3200.0, 5550.0),    # 142
        (RD,  950.0, 5350.0),     # 143
        (RD,  950.0, 5550.0),     # 144
        (RD,  3050.0, 5350.0),    # 145
        (RD,  3050.0, 5550.0),    # 146
        
        (RD,  1150.0, 5350.0),    # 147 -> SHEAR WALL 03 - DISTANCE 200.0 mm
        (RD,  1150.0, 5550.0),    # 148
        (RD,  2850.0, 5350.0),    # 149
        (RD,  2850.0, 5550.0),    # 150
        (RD,  1150.0, 5350.0),    # 151
        (RD,  1150.0, 5550.0),    # 152
        (RD,  2850.0, 5350.0),    # 153
        (RD,  2850.0, 5550.0),    # 154
        (RD,  1350.0, 5350.0),    # 155
        (RD,  1350.0, 5550.0),    # 156
        (RD,  2650.0, 5350.0),    # 157
        (RD,  2650.0, 5550.0),    # 158
        (RD,  1550.0, 5350.0),    # 159
        (RD,  1550.0, 5550.0),    # 160
        (RD,  2450.0, 5350.0),    # 161
        (RD,  2450.0, 5550.0),    # 162
        (RD,  1750.0, 5350.0),    # 163
        (RD,  1750.0, 5550.0),    # 164
        (RD,  2250.0, 5350.0),    # 165
        (RD,  2250.0, 5550.0),    # 166
        (RD,  1950.0, 5350.0),    # 167
        (RD,  1950.0, 5550.0),    # 168
        (RD,  2050.0, 5350.0),    # 169
        (RD,  2050.0, 5550.0),    # 170
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
        ax.set_title('H-shape Shear-wall Concrete section (Numbers shown) with rebars')
        ax.grid(True, ls='--', alpha=0.5)
        
        # Concrete geometry derived from concrete_layers
        ii = 0
        for depth, width, center_y, center_x, _, _, _, _ in concrete_layers:
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
        for depth, width, center_y, center_x,_ ,_, _, _ in concrete_layers:
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