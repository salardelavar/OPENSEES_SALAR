def COMPOSITE_FRP_SECTION_DIFF_WIDTH_HEIGHT_FUN_EXTRA(secTag, STEEL_TYPE, fc, Kfc, CONCRETE_DENSITY, plot=True):
    """
    Create a FRP composite with different widths and heights, fiber section (OpenSees) and,
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
    # YOUTUBE: Stress-Strain Model of Concrete Confined by FRP Laminate and Spike Anchors
    "https://www.youtube.com/watch?v=rtURZ9yBL_I"
    # YOUTUBE: FRP vs. Steel Stress Strain Curve
    "https://www.youtube.com/watch?v=WBkztwXZieQ"
    # PAPER: How Concrete Filling Fundamentally Changes Stress–Strain Curve of Angle-Ply FRP Tubes in Tension
    "https://ascelibrary.org/doi/10.1061/%28ASCE%29CC.1943-5614.0001245"
    # PAPER: Structural behavior of MRPC beams exposure to riverine simulated circumstances using GFRP and CFRP bars
    "https://www.researchgate.net/publication/349611506_Structural_behavior_of_MRPC_beams_exposure_to_riverine_simulated_circumstances_using_GFRP_and_CFRP_bars?_tp=eyJjb250ZXh0Ijp7ImZpcnN0UGFnZSI6Il9kaXJlY3QiLCJwYWdlIjoiX2RpcmVjdCJ9fQ"
    # PAPER: Experimental Investigation and Modeling of the Thermal Effect on the Mechanical Properties of Polyethylene-Terephthalate FRP Laminates
    "https://www.researchgate.net/publication/343628893_Experimental_Investigation_and_Modeling_of_the_Thermal_Effect_on_the_Mechanical_Properties_of_Polyethylene-Terephthalate_FRP_Laminates?_tp=eyJjb250ZXh0Ijp7ImZpcnN0UGFnZSI6Il9kaXJlY3QiLCJwYWdlIjoiX2RpcmVjdCJ9fQ"
    # PAPER: A State-of-the-Art Review of FRP-Confined Steel-Reinforced Concrete (FCSRC) Structural Members
    "https://www.mdpi.com/2073-4360/14/4/677"    
    
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
    
    # FRP Section
    fyF = 4800.0       # [N/mm²] FRP Section Yield Strength   
    eyF = 0.02         # [mm/mm] FRP Section Yield Strain
    EsF = fyF/eyF      # [N/mm²] FRP Section Modulus of Elasticity
    
    coreTag, coverTag, steelTag, steelITag, steelFTag = secTag + 1000, secTag + 2000, secTag + 3000, secTag + 4000, secTag + 5000
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

    # FRP Section
    ops.uniaxialMaterial('Elastic', steelFTag, EsF, 0.0, 0.1*EsF)
    #ops.uniaxialMaterial('Elastic', steelFTag, EsF)
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/ElasticUni.html
    
    # -------------------- Fiber section -------------------------
    ops.section('Fiber', secTag)
    # -------------------- Concrete Section -------------------------------
    # (layers depth[mm], layers width[mm], center y‑coord [mm], center x‑coord [mm], numSubdivY, numSubdivZ, tag, color)

    PT = 2.0     # [mm] FRP Thickness
    CD = 500.0   # [mm] Section Depth
    CW = 700.0   # [mm] Section Width
    N = 7        # Number of Division
    stepY = CD / N   # [mm] Dimension of each fiber
    stepX = CW / N   # [mm] Dimension of each fiber
    NFY, NFX = 10, 10 # Number of fibers in each Division in Y and X Dir.
    mat_layers = [
        (stepY, stepX, (i + 0.5) * stepY, (j + 0.5) * stepX, NFY, NFX, coreTag, 'lightgray')
        for i in range(N)   # N row from bottom (i=0) to top (i=N-1)
        for j in range(N)   # N columns from left (j=0) to right (j=N-1)
        ] + [
        (PT, CW+2*PT, CD+0.5*PT, 0.5*CW, NFY, NFX, steelFTag, 'cyan'),                     # 50 -> FRP
        (PT, CW+2*PT, 0.0-0.5*PT, 0.5*CW, NFY, NFX, steelFTag, 'cyan'),                    # 51 -> FRP
        (CD, PT, 0.5*CD, CW+0.5*PT, NFY, NFX, steelFTag, 'cyan'),                          # 52 -> FRP  
        (CD, PT, 0.5*CD, -0.5*PT, NFY, NFX, steelFTag, 'cyan'),                            # 53 -> FRP
        # STEEL SECTION
        (10.0, 150.0, 0.5*CD-125.0, 0.5*CW, NFY, NFX, steelITag, 'lime'),                  # 54 -> CROSSED SECTION
        (240.0, 10.0, 0.5*CD, 0.5*CW, NFY, NFX, steelITag, 'lime'),                        # 55 -> CROSSED SECTION
        (10.0, 150.0, 0.5*CD+125.0, 0.5*CW, NFY, NFX, steelITag, 'lime'),                  # 56 -> CROSSED SECTION   
         
        (10.0, 105.0, 0.5*CD, 0.5*CW-57.5, NFY, NFX, steelITag, 'lime'),                   # 57 -> CROSSED SECTION
        (10.0, 105.0, 0.5*CD, 0.5*CW+57.5, NFY, NFX, steelITag, 'lime'),                   # 68 -> CROSSED SECTION
        
        (150.0, 10.0, 0.5*CD, 0.5*CW-115.0, NFY, NFX, steelITag, 'lime'),                  # 59 -> CROSSED SECTION
        (150.0, 10.0, 0.5*CD, 0.5*CW+115.0, NFY, NFX, steelITag, 'lime'),                  # 60 -> CROSSED SECTION
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

    # -------------------- Plot (optional) -----------------------
    if plot:
        
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.set_xlabel('Width (mm)')
        ax.set_ylabel('Height (mm)')
        ax.set_title(f'Composite section \n with rebars and FRP (Numbers shown) - Section Tag: {secTag}')
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