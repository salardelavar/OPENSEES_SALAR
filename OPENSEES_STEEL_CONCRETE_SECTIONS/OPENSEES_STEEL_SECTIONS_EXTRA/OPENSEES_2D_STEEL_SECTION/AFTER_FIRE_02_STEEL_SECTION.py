def AFTER_FIRE_02_STEEL_SECTION_DIFF_WIDTH_HEIGHT_FUN_EXTRA_TWO(secTag, STEEL_TYPE, fc, Kfc, CONCRETE_DENSITY, plot=True):
    """
    Create a steel with different widths and heights, fiber section (OpenSees) and,
    if PLOT=True, draw a simple 2‑D sketch of the section.

    Parameters
    ----------
    secTag          : int   – section identifier
    STEEL_TYPE      : str   – 'ELASTIC' or 'INELASTIC'
    fc              : float – unconfined concrete compressive strength (MPa, negative)
    Kfc             : float – ratio of confined to unconfined concrete strength
    CONCRETE_DENSITY: float – ρc (kg/m³) – will be converted to N·s²/mm³
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
    esuI01 = 0.25            # [mm/mm] Steel I Section Ultimate Strain
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
    esuI02 = 0.25            # [mm/mm] Steel I Section Ultimate Strain
    EshI02 = (fuI02 - fyI02)/(esuI02 - eyI02)
    BsI02 = EshI02 / EsI02
    
    # Plate Section 03
    kE, kFy, kFu, esu = STEEL_FIRE_PROPERTIES(500) # THERMAL EFFECTS 500°C
    EsI03 = kE*200000        # [N/mm²] I Section Modulus of Elasticity 
    fyI03 = kFy*240          # [N/mm²] Steel I Section Yield Strength
    fuI03 = kFu*(1.18*240)   # [N/mm²] Steel I Section Ultimate Strength
    #EsI01 = kE*200000        # [N/mm²] I Section Modulus of Elasticity 
    #fyI01 = kFy*240          # [N/mm²] Steel I Section Yield Strength
    fuI01 = kFu*(1.18*240)   # [N/mm²] Steel I Section Ultimate Strength
    fyI03 = 95              # [N/mm²] Steel I Section Yield Strength   
    EsI03 = 2e5              # [N/mm²] I Section Modulus of Elasticity
    eyI03 = fyI03/EsI03      # [mm/mm] Steel I Section Yield Strain
    #fuI03 = 1.1818*fyI03     # [N/mm²] Steel I Section Ultimate Strength
    esuI03 = 0.25            # [mm/mm] Steel I Section Ultimate Strain
    EshI03 = (fuI03 - fyI03)/(esuI03 - eyI03)
    BsI03 = EshI03 / EsI03
    
    coreTag, coverTag, steelTag = secTag + 100, secTag + 200, secTag + 300
    steelITag01, steelITag02, steelITag03 = secTag + 400, secTag + 500, secTag + 600
    if STEEL_TYPE == 'ELASTIC':
        ops.uniaxialMaterial('Elastic', steelTag, Es) 
        ops.uniaxialMaterial('Elastic', steelITag01, EsI01) 
        ops.uniaxialMaterial('Elastic', steelITag02, EsI02) 
        ops.uniaxialMaterial('Elastic', steelITag03, EsI03) 
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
        
    #ops.uniaxialMaterial('Concrete01', coreTag, fcC, ec0C, fcUC, ecuC)  # Core concrete (confined)
    #ops.uniaxialMaterial('Concrete01', coverTag, fcU, ec0U, fcUU, ecuU) # Cover concrete (unconfined)
    ops.uniaxialMaterial('Concrete02', coreTag, fcC, ec0C, fcUC, ecuC, LambdaC, ftC, EtsC) # build core concrete (confined)
    ops.uniaxialMaterial('Concrete02', coverTag, fcU, ec0U, fcUU, ecuU, LambdaU, ftU, EtsU) # build cover concrete (unconfined)


    # -------------------- Fiber section -------------------------
    ops.section('Fiber', secTag)
    # -------------------- Concrete Section -------------------------------
    # (layers depth[mm], layers width[mm], center y‑coord [mm], center x‑coord [mm], numSubdivY, numSubdivZ, tag, color)

    NFY, NFX = 10, 10 # Number of fibers in each Division in Y and X Dir.
    mat_layers = [
        (10.0, 100.0, 5.0, 50.0, NFY, NFX, steelITag03, 'pink'),          # 1
        (160.0, 10.0, 90.0, 50.0, NFY, NFX, steelITag03, 'pink'),         # 2
        (10.0, 100.0, 175.0, 50.0, NFY, NFX, steelITag03, 'pink'),        # 3    
        
        (10.0, 100.0, 5.0, 150.0, NFY, NFX, steelITag01, 'lime'),         # 4
        (160.0, 10.0, 90.0, 150.0, NFY, NFX, steelITag01, 'lime'),        # 5
        (10.0, 100.0, 175.0, 150.0, NFY, NFX, steelITag01, 'lime'),       # 6 
        
        (6.0, 150.0, -3.0, 100.0, NFY, NFX, steelITag01, 'lime'),         # 7
        (6.0, 150.0, 183.0, 100.0, NFY, NFX, steelITag03, 'pink'),        # 8
        
        (210.0, 8.0, 90.0, -4.0, NFY, NFX, steelITag03, 'pink'),          # 9
        (210.0, 8.0, 90.0, 204.0, NFY, NFX, steelITag01, 'lime'),         # 10

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
    """
    # (diameter [mm], y‑coord [mm], x‑coord [mm])
    rebars = [
        (25.0,  50.0, 50.0),          # 1
        (18.0,  50.0, 0.5*CD),        # 2
        (25.0,  50.0, CW-50.0),       # 3 
        (25.0,  CD-50.0, 50.0),       # 4
        (18.0,  CD-50.0, 0.5*CD),     # 5
        (25.0,  CD-50.0, CW-50.0),    # 6
        
        (16.0,  0.5*CD, 50.0),        # 7
        (16.0,  0.5*CD, CW-50.0),     # 8
        
        (16.0,  0.25*CD, 50.0),       # 9
        (16.0,  0.25*CD, CW-50.0),    # 10
        
        (16.0,  0.75*CD, 50.0),       # 11
        (16.0,  0.75*CD, CW-50.0),    # 12
        
        (18.0,  50.0, 0.25*CW),       # 13
        (18.0,  50.0, 0.75*CW),       # 14
    
        (18.0,  CD-50.0, 0.25*CW),    # 15
        (18.0,  CD-50.0, 0.75*CW),    # 16
    ]

    for dia, y, x in rebars:
        area = np.pi * dia**2 / 4.0          # mm²
        ops.fiber(x, y, area, steelTag)      # add to OpenSees model
    """    
    # -------------------- Plot (optional) -----------------------
    if plot:
        
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.set_xlabel('Width (mm)')
        ax.set_ylabel('Height (mm)')
        ax.set_title(f'After Fire Effects Steel Section - Section Tag: {secTag}')
        #ax.set_title(f'After Fire Effects Steel Section \n with rebars and plates (Numbers shown) - Section Tag: {secTag}')
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
        """
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
        """
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