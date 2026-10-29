def COMPOSITE_SHEAR_WALL_SECTION_DIFF_WIDTH_HEIGHT_FUN_EXTRA(secTag, STEEL_TYPE, fc, Kfc, CONCRETE_DENSITY, plot=True):
    """
    Create a composite with different widths and heights, fiber section (OpenSees) and,
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
    if STEEL_TYPE == 'ELASTIC':
        ops.uniaxialMaterial('Elastic', steelTag, Es) 
        ops.uniaxialMaterial('Elastic', steelITag, EsI) 
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


    # -------------------- Fiber section -------------------------
    ops.section('Fiber', secTag)
    # -------------------- Concrete Section -------------------------------
    # (layers depth[mm], layers width[mm], center y‑coord [mm], center x‑coord [mm], numSubdivY, numSubdivZ, tag, color)

    PT = 10.0    # [mm] Plate Thickness
    CD = 5000.0  # [mm] Section Depth
    CW = 500.0   # [mm] Section Width
    N = 7        # Number of Division
    stepY = CD / N   # [mm] Dimension of each fiber
    stepX = CW / N   # [mm] Dimension of each fiber
    NFY, NFX = 10, 10 # Number of fibers in each Division in Y and X Dir.
    mat_layers = [
        (stepY, stepX, (i + 0.5) * stepY, (j + 0.5) * stepX, NFY, NFX, coreTag, 'lightgray')
        for i in range(N)   # N row from bottom (i=0) to top (i=N-1)
        for j in range(N)   # N columns from left (j=0) to right (j=N-1)
        ] + [
        (PT, CW+2*PT, CD+0.5*PT, 0.5*CW, NFY, NFX, steelITag, 'cyan'),                # 50 -> PLATE
        (PT, CW+2*PT, 0.0-0.5*PT, 0.5*CW, NFY, NFX, steelITag, 'cyan'),               # 51 -> PLATE
        (CD, PT, 0.5*CD, CW+0.5*PT, NFY, NFX, steelITag, 'cyan'),                     # 52 -> PLATE   
        (CD, PT, 0.5*CD, -0.5*PT, NFY, NFX, steelITag, 'cyan'),                       # 53 -> PLATE
        # BOTTOM STEEL SECTION
        (10.0, 150.0, 0.5*CW-0.25*CW, 0.5*CW, NFY, NFX, steelITag, 'lime'),                # 54 -> CROSSED SECTION
        (240.0, 10.0, 0.5*CW, 0.5*CW, NFY, NFX, steelITag, 'lime'),                        # 55 -> CROSSED SECTION
        (10.0, 150.0, 0.5*CW+0.25*CW, 0.5*CW, NFY, NFX, steelITag, 'lime'),                # 56 -> CROSSED SECTION   
         
        (10.0, 105.0, 0.5*CW, 0.5*CW-57.5, NFY, NFX, steelITag, 'lime'),                   # 57 -> CROSSED SECTION
        (10.0, 105.0, 0.5*CW, 0.5*CW+57.5, NFY, NFX, steelITag, 'lime'),                   # 58 -> CROSSED SECTION
        
        (150.0, 10.0, 0.5*CW, 0.5*CW-115.0, NFY, NFX, steelITag, 'lime'),                  # 59 -> CROSSED SECTION
        (150.0, 10.0, 0.5*CW, 0.5*CW+115.0, NFY, NFX, steelITag, 'lime'),                  # 60 -> CROSSED SECTION
        # TOP STEEL SECTION        
        (10.0, 150.0, CD-(0.5*CW-0.25*CW), 0.5*CW, NFY, NFX, steelITag, 'lime'),                # 61 -> CROSSED SECTION
        (240.0, 10.0, CD-(0.5*CW), 0.5*CW, NFY, NFX, steelITag, 'lime'),                        # 62 -> CROSSED SECTION
        (10.0, 150.0, CD-(0.5*CW+0.25*CW), 0.5*CW, NFY, NFX, steelITag, 'lime'),                # 63 -> CROSSED SECTION   
         
        (10.0, 105.0, CD-(0.5*CW), 0.5*CW-57.5, NFY, NFX, steelITag, 'lime'),                   # 64 -> CROSSED SECTION
        (10.0, 105.0, CD-(0.5*CW), 0.5*CW+57.5, NFY, NFX, steelITag, 'lime'),                   # 65 -> CROSSED SECTION
        
        (150.0, 10.0, CD-(0.5*CW), 0.5*CW-115.0, NFY, NFX, steelITag, 'lime'),                  # 66 -> CROSSED SECTION
        (150.0, 10.0, CD-(0.5*CW), 0.5*CW+115.0, NFY, NFX, steelITag, 'lime'),                  # 67 -> CROSSED SECTION
        
        (CD-2*240-2*140.0, 20.0, 0.5*CD, 0.5*CW, NFY, NFX, steelITag, 'purple'),                  # 67 -> CROSSED SECTION
    ]
    
    """
    mat_layers = [
        (100.0, 100.0, 50.0,   50.0,  NFY, NFX, coreTag, 'lightgray'),   # 1
        (100.0, 100.0, 50.0,  150.0,  NFY, NFX, coreTag, 'lightgray'),   # 2
        (100.0, 100.0, 50.0,  250.0,  NFY, NFX, coreTag, 'lightgray'),   # 3
        (100.0, 100.0, 50.0,  350.0,  NFY, NFX, coreTag, 'lightgray'),   # 4
        (100.0, 100.0, 50.0,  450.0,  NFY, NFX, coreTag, 'lightgray'),   # 5
        (100.0, 100.0, 150.0,  50.0,  NFY, NFX, coreTag, 'lightgray'),   # 6
        (100.0, 100.0, 150.0, 150.0,  NFY, NFX, coreTag, 'lightgray'),   # 7
        (100.0, 100.0, 150.0, 250.0,  NFY, NFX, coreTag, 'lightgray'),   # 8
        (100.0, 100.0, 150.0, 350.0,  NFY, NFX, coreTag, 'lightgray'),   # 9
        (100.0, 100.0, 150.0, 450.0,  NFY, NFX, coreTag, 'lightgray'),   # 10
        (100.0, 100.0, 250.0,  50.0,  NFY, NFX, coreTag, 'lightgray'),   # 11
        (100.0, 100.0, 250.0, 150.0,  NFY, NFX, coreTag, 'lightgray'),   # 12
        (100.0, 100.0, 250.0, 250.0,  NFY, NFX, coreTag, 'lightgray'),   # 13
        (100.0, 100.0, 250.0, 350.0,  NFY, NFX, coreTag, 'lightgray'),   # 14
        (100.0, 100.0, 250.0, 450.0,  NFY, NFX, coreTag, 'lightgray'),   # 15
        (100.0, 100.0, 350.0,  50.0,  NFY, NFX, coreTag, 'lightgray'),   # 16
        (100.0, 100.0, 350.0, 150.0,  NFY, NFX, coreTag, 'lightgray'),   # 17
        (100.0, 100.0, 350.0, 250.0,  NFY, NFX, coreTag, 'lightgray'),   # 18
        (100.0, 100.0, 350.0, 350.0,  NFY, NFX, coreTag, 'lightgray'),   # 19
        (100.0, 100.0, 350.0, 450.0,  NFY, NFX, coreTag, 'lightgray'),   # 20
        (100.0, 100.0, 450.0,  50.0,  NFY, NFX, coreTag, 'lightgray'),   # 21
        (100.0, 100.0, 450.0, 150.0,  NFY, NFX, coreTag, 'lightgray'),   # 22
        (100.0, 100.0, 450.0, 250.0,  NFY, NFX, coreTag, 'lightgray'),   # 23
        (100.0, 100.0, 450.0, 350.0,  NFY, NFX, coreTag, 'lightgray'),   # 24
        (100.0, 100.0, 450.0, 450.0,  NFY, NFX, coreTag, 'lightgray'),   # 25
    
        (PT, 500.0+2*PT, 500.0+0.5*PT, 0.5*500.0, NFY, NFX, steelITag, 'cyan'),   # 26 -> PLATE
        (PT, 500.0+2*PT, 0.0-0.5*PT, 0.5*500.0, NFY, NFX, steelITag, 'cyan'),     # 27 -> PLATE
        (500.0, PT, 0.5*500.0, 500.0+0.5*PT, NFY, NFX, steelITag, 'cyan'),        # 28 -> PLATE   
        (500.0, PT, 0.5*500.0, -0.5*PT, NFY, NFX, steelITag, 'cyan'),             # 29 -> PLATE

        (10.0, 150.0, 250.0-125.0, 250.0, NFY, NFX, steelITag, 'lime'),                  # 30 -> CROSSED SECTION
        (240.0, 10.0, 250.0, 250.0, NFY, NFX, steelITag, 'lime'),                        # 31 -> CROSSED SECTION
        (10.0, 150.0, 250+125.0, 250.0, NFY, NFX, steelITag, 'lime'),                    # 32 -> CROSSED SECTION   
         
        (10.0, 105.0, 250.0, 250.0-57.5, NFY, NFX, steelITag, 'lime'),                   # 33 -> CROSSED SECTION
        (10.0, 105.0, 250.0, 250.0+57.5, NFY, NFX, steelITag, 'lime'),                   # 34 -> CROSSED SECTION
        
        (150.0, 10.0, 250.0, 250.0-115.0, NFY, NFX, steelITag, 'lime'),                  # 35 -> CROSSED SECTION
        (150.0, 10.0, 250.0, 250.0+115.0, NFY, NFX, steelITag, 'lime'),                  # 36 -> CROSSED SECTION     
         ]
    """
    """
    mat_layers = [
        (20.0, 500.0, 10.0, 250.0, NFY, NFX, coreTag, 'lightgray'),    # 1
        (20.0, 500.0, 30.0, 250.0, NFY, NFX, coreTag, 'lightgray'),    # 2
        (20.0, 500.0, 50.0, 250.0, NFY, NFX, coreTag, 'lightgray'),    # 3
        (20.0, 500.0, 70.0, 250.0, NFY, NFX, coreTag, 'lightgray'),    # 4
        (20.0, 500.0, 90.0, 250.0, NFY, NFX, coreTag, 'lightgray'),    # 5
        (20.0, 500.0, 110.0, 250.0, NFY, NFX, coreTag, 'lightgray'),   # 6
        (20.0, 500.0, 130.0, 250.0, NFY, NFX, coreTag, 'lightgray'),   # 7
        (20.0, 500.0, 150.0, 250.0, NFY, NFX, coreTag, 'lightgray'),   # 8
        (20.0, 500.0, 170.0, 250.0, NFY, NFX, coreTag, 'lightgray'),   # 9
        (20.0, 500.0, 190.0, 250.0, NFY, NFX, coreTag, 'lightgray'),   # 10
        (20.0, 500.0, 210.0, 250.0, NFY, NFX, coreTag, 'lightgray'),   # 11
        (20.0, 500.0, 230.0, 250.0, NFY, NFX, coreTag, 'lightgray'),   # 12 
        (20.0, 500.0, 250.0, 250.0, NFY, NFX, coreTag, 'lightgray'),   # 13
        (20.0, 500.0, 270.0, 250.0, NFY, NFX, coreTag, 'lightgray'),   # 14
        (20.0, 500.0, 290.0, 250.0, NFY, NFX, coreTag, 'lightgray'),   # 15
        (20.0, 500.0, 310.0, 250.0, NFY, NFX, coreTag, 'lightgray'),   # 16
        (20.0, 500.0, 330.0, 250.0, NFY, NFX, coreTag, 'lightgray'),   # 17
        (20.0, 500.0, 350.0, 250.0, NFY, NFX, coreTag, 'lightgray'),   # 18
        (20.0, 500.0, 370.0, 250.0, NFY, NFX, coreTag, 'lightgray'),   # 19
        (20.0, 500.0, 390.0, 250.0, NFY, NFX, coreTag, 'lightgray'),   # 20
        (20.0, 500.0, 410.0, 250.0, NFY, NFX, coreTag, 'lightgray'),   # 21
        (20.0, 500.0, 430.0, 250.0, NFY, NFX, coreTag, 'lightgray'),   # 22
        (20.0, 500.0, 450.0, 250.0, NFY, NFX, coreTag, 'lightgray'),   # 23
        (20.0, 500.0, 470.0, 250.0, NFY, NFX, coreTag, 'lightgray'),   # 24
        (20.0, 500.0, 490.0, 250.0, NFY, NFX, coreTag, 'lightgray'),   # 25
        
        (PT, 500.0+2*PT, 500.0+0.5*PT, 0.5*500.0, NFY, NFX, steelITag, 'cyan'),   # 26 -> PLATE
        (PT, 500.0+2*PT, 0.0-0.5*PT, 0.5*500.0, NFY, NFX, steelITag, 'cyan'),     # 27 -> PLATE
        (500.0, PT, 0.5*500.0, 500.0+0.5*PT, NFY, NFX, steelITag, 'cyan'),        # 28 -> PLATE   
        (500.0, PT, 0.5*500.0, -0.5*PT, NFY, NFX, steelITag, 'cyan'),             # 29 -> PLATE

        (10.0, 150.0, 250.0-125.0, 250.0, NFY, NFX, steelITag, 'lime'),                  # 30 -> CROSSED SECTION
        (240.0, 10.0, 250.0, 250.0, NFY, NFX, steelITag, 'lime'),                        # 31 -> CROSSED SECTION
        (10.0, 150.0, 250+125.0, 250.0, NFY, NFX, steelITag, 'lime'),                    # 32 -> CROSSED SECTION   
         
        (10.0, 105.0, 250.0, 250.0-57.5, NFY, NFX, steelITag, 'lime'),                   # 33 -> CROSSED SECTION
        (10.0, 105.0, 250.0, 250.0+57.5, NFY, NFX, steelITag, 'lime'),                   # 34 -> CROSSED SECTION
        
        (150.0, 10.0, 250.0, 250.0-115.0, NFY, NFX, steelITag, 'lime'),                  # 35 -> CROSSED SECTION
        (150.0, 10.0, 250.0, 250.0+115.0, NFY, NFX, steelITag, 'lime'),                  # 36 -> CROSSED SECTION   
        ]
    """
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
    # -------------------- Rebars -------------------------------
    # (diameter [mm], y‑coord [mm], x‑coord [mm])
    rebars = [
        (25.0,  50.0, 50.0),    # 1
        (25.0,  50.0, 250.0),   # 2
        (25.0,  50.0, 450.0),   # 3 
        (25.0,  450.0, 50.0),   # 4
        (25.0,  450.0, 450.0),  # 5
        
        (25.0,  250.0, 50.0),   # 6
        (25.0,  250.0, 450.0),  # 7
        
        (25.0,  CD-(50.0), 50.0),    # 8
        (25.0,  CD-(50.0), 250.0),   # 9
        (25.0,  CD-(50.0), 450.0),   # 10 
        (25.0,  CD-(450.0), 50.0),   # 11
        (25.0,  CD-(450.0), 450.0),  # 12
        
        (25.0,  CD-(250.0), 50.0),   # 13
        (25.0,  CD-(250.0), 450.0),  # 14
        
        (12.0,  CD-(450.0)-1*200, 50.0),   # 15
        (12.0,  CD-(450.0)-1*200, 450.0),  # 16
        
        (12.0,  CD-(450.0)-2*200, 50.0),   # 17
        (12.0,  CD-(450.0)-2*200, 450.0),  # 18
        
        (12.0,  CD-(450.0)-3*200, 50.0),   # 19
        (12.0,  CD-(450.0)-3*200, 450.0),  # 20
        
        (12.0,  CD-(450.0)-4*200, 50.0),   # 21
        (12.0,  CD-(450.0)-4*200, 450.0),  # 22
        
        (12.0,  CD-(450.0)-5*200, 50.0),   # 23
        (12.0,  CD-(450.0)-5*200, 450.0),  # 24
        
        (12.0,  CD-(450.0)-6*200, 50.0),   # 25
        (12.0,  CD-(450.0)-6*200, 450.0),  # 26
        
        (12.0,  CD-(450.0)-7*200, 50.0),   # 27
        (12.0,  CD-(450.0)-7*200, 450.0),  # 28
        
        (12.0,  CD-(450.0)-8*200, 50.0),   # 29
        (12.0,  CD-(450.0)-8*200, 450.0),  # 30
        
        (12.0,  CD-(450.0)-9*200, 50.0),   # 31
        (12.0,  CD-(450.0)-9*200, 450.0),  # 32
        
        (12.0,  CD-(450.0)-10*200, 50.0),   # 33
        (12.0,  CD-(450.0)-10*200, 450.0),  # 34
        
        (12.0,  CD-(450.0)-11*200, 50.0),   # 35
        (12.0,  CD-(450.0)-11*200, 450.0),  # 36
        
        (12.0,  CD-(450.0)-12*200, 50.0),   # 37
        (12.0,  CD-(450.0)-12*200, 450.0),  # 38
        
        (12.0,  CD-(450.0)-13*200, 50.0),   # 39
        (12.0,  CD-(450.0)-13*200, 450.0),  # 40
        
        (12.0,  CD-(450.0)-14*200, 50.0),   # 41
        (12.0,  CD-(450.0)-14*200, 450.0),  # 42
        
        (12.0,  CD-(450.0)-15*200, 50.0),   # 43
        (12.0,  CD-(450.0)-15*200, 450.0),  # 44
        
        (12.0,  CD-(450.0)-16*200, 50.0),   # 45
        (12.0,  CD-(450.0)-16*200, 450.0),  # 46
        
        (12.0,  CD-(450.0)-17*200, 50.0),   # 47
        (12.0,  CD-(450.0)-17*200, 450.0),  # 48
        
        (12.0,  CD-(450.0)-18*200, 50.0),   # 49
        (12.0,  CD-(450.0)-18*200, 450.0),  # 50
        
        (12.0,  CD-(450.0)-19*200, 50.0),   # 51
        (12.0,  CD-(450.0)-19*200, 450.0),  # 52

        
    ]

    for dia, y, x in rebars:
        area = np.pi * dia**2 / 4.0          # mm²
        ops.fiber(x, y, area, steelTag)      # add to OpenSees model

    # -------------------- Plot (optional) -----------------------
    if plot:
        
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.set_xlabel('Width (mm)')
        ax.set_ylabel('Height (mm)')
        ax.set_title(f'Composite shear-wall section \n with rebars and plates (Numbers shown) - Section Tag: {secTag}')
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
        
        # Compute overall limits from the actual concrete layers
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
    ELE_MASS = CONCRETE_DENSITY * TOTAL_AREA   # kg/mm
    
    return SECTION_HEIGHT, ELE_MASS