def PRESTRESS_CONCRETE_I_SECTION_DIFF_FLANGE_FUN_3D(secTag, matTag, STEEL_TYPE, fc, Kfc,
                           bfT, tfT, bfB, tfB, h, tw,
                           nFib, CONCRETE_DENSITY,
                           plot=True):
    """
    Create a Prestressed I confined‑concrete with different flanges width, fiber section (OpenSees) and,
    if PLOT=True, draw a simple 2‑D sketch of the section.

    Parameters
    ----------
    secTag          : int   – section identifier
    h, bf           : float – height and width  (mm)
    cover           : float – concrete cover (mm)
    As              : float – area of one reinforcement bar (mm²)
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
    EtsC = ftC/np.abs(ec0C)	   # [N/mm²] tension softening stiffness
    EtsU = ftU/np.abs(ec0U)	   # [N/mm²] tension softening stiffness

    # STEEL
    # Reinforcing Steel Properties
    fy = 400          # [N/mm²] Steel Rebar Yield Strength   
    Es = 2e5          # [N/mm²] Modulus of Elasticity
    ey = fy/Es        # [mm/mm] Steel Rebar Yield Strain
    fu = 1.1818*fy    # [N/mm²] Steel Rebar Ultimate Strength
    esu = 0.09        # [mm/mm] Steel Rebar Ultimate Strain
    Esh = (fu - fy)/(esu - ey)
    Bs = Esh / Es
    pinchX = 0.8   # Pinching factor in X direction
    pinchY = 0.5   # Pinching factor in Y direction
    damage1 = 0.0  # Damage due to ductility
    damage2 = 0.0  # Damage due to energy
    beta = 0.1     # Stiffness degradation parameter
    # Prestressed Rebar or Cable Steel Properties (HIGH STRENGTH STEEL)
    fyC = 1200         # [N/mm²] Steel Cable Yield Strength   
    EsC = 2e5          # [N/mm²] Modulus of Elasticity
    eyC = fyC/EsC      # [mm/mm] Steel Cable Yield Strain
    fuC = 1.1818*fyC   # [N/mm²] Steel Cable Ultimate Strength
    esuC = 0.06        # [mm/mm] Steel Cable Ultimate Strain
    EshC = (fuC - fyC)/(esuC - eyC)
    BsC = EshC / EsC
    pinchXC = 0.8   # Pinching factor in X direction
    pinchYC = 0.5   # Pinching factor in Y direction
    damage1C = 0.0  # Damage due to ductility
    damage2C = 0.0  # Damage due to energy
    betaC = 0.1     # Stiffness degradation parameter
    
    coreTag, coverTag, steelTag, cableTag, cableITag = secTag + 100, secTag + 200, secTag + 300, secTag + 400, secTag + 500
    if STEEL_TYPE == 'ELASTIC':
        #ops.uniaxialMaterial('Steel01', steelTag, fy, Es, Bs) 
        #ops.uniaxialMaterial('Steel01', cableTag, fyC, EsC, BsC) 
        ops.uniaxialMaterial('Elastic', steelTag, Es)
        ops.uniaxialMaterial('Elastic', cableTag, EsC)
    if STEEL_TYPE == 'INELASTIC':   
        # REBAR
        ops.uniaxialMaterial('Hysteretic', steelTag, fy, ey,
                             fu, esu,
                             0.2*fu, 1.1*esu,
                             -fy, -ey,
                             -fu, -esu,
                             -0.2*fu, -1.1*esu,
                             pinchX, pinchY,
                             damage1, damage2, beta)
        # PRESTRESSED CABLE
        ops.uniaxialMaterial('Hysteretic', cableTag, fyC, eyC,
                             fuC, esuC,
                             0.2*fuC, 1.1*esuC,
                             -fyC, -eyC,
                             -fuC, -esuC,
                             -0.2*fuC, -1.1*esuC,
                             pinchXC, pinchYC,
                             damage1C, damage2C, betaC)

        # INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Hysteretic_Material
        
    #ops.uniaxialMaterial('Concrete01', coreTag, fcC, ec0C, fcUC, ecuC)  # Core concrete (confined)
    #ops.uniaxialMaterial('Concrete01', coverTag, fcU, ec0U, fcUU, ecuU) # Cover concrete (unconfined)
    ops.uniaxialMaterial('Concrete02', coreTag, fcC, ec0C, fcUC, ecuC, LambdaC, ftC, EtsC) # build core concrete (confined)
    ops.uniaxialMaterial('Concrete02', coverTag, fcU, ec0U, fcUU, ecuU, LambdaU, ftU, EtsU) # build cover concrete (unconfined)

    # Define initial strain for prestressing tendon (example prestrain value)
    # LINK: https://opensees.berkeley.edu/wiki/index.php?title=Initial_Strain_Material
    InitialStrain = -0.001  # [mm/mm] Initial strain (Compression Negaive Value - Tension Postive Value)
    # uniaxialMaterial InitStrainMaterial $matTag $otherTag $initStrain
    ops.uniaxialMaterial('InitStrainMaterial', cableITag, cableTag, InitialStrain)
    # -------------------- Section area -------------------------
    hw = h-tfT-tfB
    area_top_flange = bfT * tfT
    area_bottom_flange = bfB * tfB
    area_web        = tw * hw
    AREA = area_top_flange + area_bottom_flange + area_web
    print(f"Total Section Area = {AREA:.2f} mm²")
    
    ybar = (area_top_flange * (tfB+hw+0.5*tfT) +
            area_web * (tfB+0.5*hw) +
            area_bottom_flange * tfB) / AREA
    print(f"Neutral axis in the Y‑direction = {ybar:.2f} mm")
    # -------------------- Fiber section -------------------------
    ops.section('Fiber', secTag, '-GJ', 10e6)

    def add_rect(mat, y_bot, y_top, x_left, x_right):
        ops.patch('rect', mat, nFib, nFib,
                  x_left, y_bot, x_right, y_top)

    # Concrete patches
    add_rect(coverTag,  h - tfT,  h, -bfT/2,  bfT/2)         # Top flange
    add_rect(coverTag, 0.0, tfB, -bfB/2,  bfB/2)             # Bottom flange
    add_rect(coverTag, h - tfB,  h - tfB - tfT, -tw/2, tw/2) # Web
    # -------------------- Rebars -------------------------------
    # (diameter [mm], y‑coord [mm], x‑coord [mm])
    rebars = [
        (16.0,  h - 1*tfT/3, -bfT/2 + 0.05*bfT),   # 1
        (16.0,  h - 1*tfT/3,  bfT/2 - 0.05*bfT),   # 2
        (16.0,  h - 1*tfT/3, -bfT/2 + 0.20*bfT),   # 3
        (16.0,  h - 1*tfT/3, bfT/2 - 0.20*bfT),    # 4
        (16.0, h - 1*tfT/3, -bfT/2 + 0.60*bfT),    # 5
        (16.0, h - 1*tfT/3,  bfT/2 - 0.60*bfT),    # 6
        
        (16.0,  h - 2*tfT/3, -bfT/2 + 0.05*bfT),   # 7
        (16.0,  h - 2*tfT/3,  bfT/2 - 0.05*bfT),   # 8
        (16.0,  h - 2*tfT/3, -bfT/2 + 0.20*bfT),   # 9
        (16.0,  h - 2*tfT/3, bfT/2 - 0.20*bfT),    # 10
        (16.0, h - 2*tfT/3, -bfT/2 + 0.60*bfT),    # 11
        (16.0, h - 2*tfT/3,  bfT/2 - 0.60*bfT),    # 12
        
        (12.0,  tfB+hw*0.05, -tw/4),   # 13
        (12.0,  tfB+hw*0.05,  tw/4),   # 14
        (12.0,  tfB+hw*0.15, -tw/4),   # 15
        (12.0,  tfB+hw*0.15,  tw/4),   # 16
        (12.0,  tfB+hw*0.25, -tw/4),   # 17
        (12.0,  tfB+hw*0.25,  tw/4),   # 18
        (12.0,  tfB+hw*0.35, -tw/4),   # 19
        (12.0,  tfB+hw*0.35,  tw/4),   # 20
        
        (12.0,  tfB+hw*0.45, -tw/4),   # 21
        (12.0,  tfB+hw*0.45,  tw/4),   # 22
        (12.0,  tfB+hw*0.55, -tw/4),   # 23
        (12.0,  tfB+hw*0.55,  tw/4),   # 24
        (12.0,  tfB+hw*0.65, -tw/4),   # 25
        (12.0,  tfB+hw*0.65,  tw/4),   # 26
        (12.0,  tfB+hw*0.75, -tw/4),   # 27
        (12.0,  tfB+hw*0.75,  tw/4),   # 28
        (12.0,  tfB+hw*0.85, -tw/4),   # 29
        (12.0,  tfB+hw*0.85,  tw/4),   # 30
        (12.0,  tfB+hw*0.95, -tw/4),   # 31
        (12.0,  tfB+hw*0.95,  tw/4),   # 32
        
        (16.0, 1*tfB/4, -bfB/2 + 0.10*bfB),  # 33
        (16.0, 1*tfB/4,  bfB/2 - 0.10*bfB),  # 34
        (16.0,  1*tfB/4, -bfB/2 + 0.35*bfB), # 35
        (16.0,  1*tfB/4,  bfB/2 - 0.35*bfB), # 36
        (16.0, 2*tfB/4, -bfB/2 + 0.10*bfB),  # 37
        (16.0, 2*tfB/4,  bfB/2 - 0.10*bfB),  # 38
        (16.0,  2*tfB/4, -bfB/2 + 0.35*bfB), # 39
        (16.0,  2*tfB/4,  bfB/2 - 0.35*bfB), # 40
        (16.0, 3*tfB/4, -bfB/2 + 0.10*bfB),  # 41
        (16.0, 3*tfB/4,  bfB/2 - 0.10*bfB),  # 42
        (16.0,  3*tfB/4, -bfB/2 + 0.35*bfB), # 43
        (16.0,  3*tfB/4,  bfB/2 - 0.35*bfB), # 44
    ]
    
    cables = [
        (25.0,  1*tfB/4, 0.0),   # 1
        (40.0,  2*tfB/4, 0.0),   # 2
        (25.0,  3*tfB/4, 0.0),   # 3
    ]

    for dia, y, x in rebars:
        area = np.pi * dia**2 / 4.0          # mm²
        ops.fiber(x, y, area, steelTag)      # add to OpenSees model

    # -------------------- Plot (optional) -----------------------
    if plot:
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.set_xlabel('Width (mm)')
        ax.set_ylabel('Height (mm)')
        ax.set_title('Prestressed I Section Different flange with Rebars (numbers shown)')
        ax.grid(True, ls='--', alpha=0.5)
        
        # Concrete geometry (light gray)
        geometry = [
            {"x": -bfB/2, "y": 0.0,"w": bfB, "h": tfB},                # Bottom flange
            {"x": -tw/2, "y": tfB, "w": tw, "h": hw},                  # Web
            {"x": -bfT/2, "y": h-tfT, "w": bfT, "h": tfT},             # Top flange
        ]
        for g in geometry:
            rect = patches.Rectangle((g["x"], g["y"]), g["w"], g["h"],
                                     linewidth=1.5, edgecolor='black',
                                     facecolor='lightgray')
            ax.add_patch(rect)

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
        
        # Cables – black circles + numbers
        for i, (dia, y, x) in enumerate(cables, start=1):
            # circle
            circ = patches.Circle((x, y), radius=dia/2,
                                 edgecolor='black', facecolor='black',
                                 linewidth=2)
            ax.add_patch(circ)

            # label – placed a little above the bar
            ax.text(x, y + dia/2 + 4, f'{i}',
                    color='purple', fontsize=6,
                    ha='center', va='bottom',
                    fontweight='bold')    

        max_dim = max(bfT, bfB, h) + 50
        ax.set_xlim(-max_dim/2, max_dim/2)
        ax.set_ylim(-50.0, max_dim)
        ax.set_aspect('equal')
        plt.show()

    ELE_MASS = CONCRETE_DENSITY * AREA   # [kg/mm] Mass Per Length

    return h, ELE_MASS
