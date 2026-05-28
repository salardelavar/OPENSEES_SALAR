def CONCRETE_DUMBBELL_SHEAR_WALL_DIFF_FLANGE_FUN_EXTRA_3D(secTag, STEEL_TYPE,
                                                       fc, Kfc,
                                                       bfT, tfT, bfB, tfB, h, tw,
                                                       nFib, CONCRETE_DENSITY,
                                                       plot=True):
    """
    Create a concrete Dumbbell-shaped Shear-wall I confined‑concrete with different flanges width, fiber section (OpenSees) and,
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
    
    # Core concrete (confined)
    #Kfc = 1.3;			      # ratio of confined to unconfined concrete strength
    fcC = Kfc*fcU             # [N/mm²] Concrete Compressive Strength
    Ec = 4700 * np.sqrt(-fcC) # [N/mm^2] Concrete Elastic Modulus
    ec0C = 2*fcC/Ec           # [mm/mm] Concrete Compressive Strain
    fcUC = 0.65*fcC           # [N/mm²] Concrete Compressive Ultimate Strength
    ecuC = 15*ec0C            # [mm/mm] Concrete Compressive Ultimate Strain
    Lambda = 0.1;	          # ratio between unloading slope
    # tensile-strength properties
    ftC = 0.7 * np.sqrt(-fcC)  # [N/mm²] tensile strength +tension
    ftU = 0.7 * np.sqrt(-fcU)  # [N/mm²] tensile strength +tension
    EtsC = ftC/0.002;		   # [N/mm²] tension softening stiffness
    EtsU = ftU/0.002;		   # [N/mm²] tension softening stiffness
    
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
        ops.uniaxialMaterial('Steel01', steelTag, fy, Es, Bs) 
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
    ops.uniaxialMaterial('Concrete02', coreTag, fcC, ec0C, fcUC, ecuC, Lambda, ftC, EtsC) # build core concrete (confined)
    ops.uniaxialMaterial('Concrete02', coverTag, fcU, ec0U, fcUU, ecuU, Lambda, ftU, EtsU) # build cover concrete (unconfined)

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
    cover = 50
    RD01 = 25.0
    RD02 = 16.0
    RD03 = 20.0
    rebars = [
        (RD01,  h - cover, -bfT/2 + cover),         # 1
        (RD01,  h - cover,  bfT/2 - cover),         # 2
        (RD01,  h - cover,  0.0),                   # 3
        (RD01,  h - tfT + cover, -bfT/2 + cover),   # 4
        (RD01,  h - tfT + cover, bfT/2 - cover),    # 5
        (RD01,  h - tfT + cover, 0.0),              # 6
        (RD01, h - tfT/2, -bfT/2 + cover),          # 7
        (RD01, h - tfT/2,  bfT/2 - cover),          # 8
        
        (RD02,  tfB+hw*0.1, -tw/4),                 # 9
        (RD02,  tfB+hw*0.1,  tw/4),                 # 10
        (RD02,  tfB+hw*0.2, -tw/4),                 # 11
        (RD02,  tfB+hw*0.2,  tw/4),                 # 12
        (RD02,  tfB+hw*0.3, -tw/4),                 # 13
        (RD02,  tfB+hw*0.3,  tw/4),                 # 14
        (RD02,  tfB+hw*0.4, -tw/4),                 # 15
        (RD02,  tfB+hw*0.4,  tw/4),                 # 16
        (RD02,  tfB+hw*0.5, -tw/4),                 # 17
        (RD02,  tfB+hw*0.5,  tw/4),                 # 18
        (RD02,  tfB+hw*0.6, -tw/4),                 # 19
        (RD02,  tfB+hw*0.6,  tw/4),                 # 20
        (RD02,  tfB+hw*0.7, -tw/4),                 # 21
        (RD02,  tfB+hw*0.7,  tw/4),                 # 22
        (RD02,  tfB+hw*0.8, -tw/4),                 # 23
        (RD02,  tfB+hw*0.8,  tw/4),                 # 24
        (RD02,  tfB+hw*0.9, -tw/4),                 # 23
        (RD02,  tfB+hw*0.9,  tw/4),                 # 24
        (RD02,  tfB+hw*1.0, -tw/4),                 # 25
        (RD02,  tfB+hw*1.0,  tw/4),                 # 26
        (RD02,  tfB+hw*0.0, -tw/4),                 # 27
        (RD02,  tfB+hw*0.0,  tw/4),                 # 28
        
        (RD03,  cover, -bfB/2 + cover),             # 29
        (RD03,  cover,  bfB/2 - cover),             # 30
        (RD03,  tfB - cover,  -bfB/2 + cover),      # 31
        (RD03,  tfB - cover, bfB/2 - cover),        # 32
    ]

    for dia, y, x in rebars:
        area = np.pi * dia**2 / 4.0          # mm²
        ops.fiber(x, y, area, steelTag)      # add to OpenSees model

    # -------------------- Plot (optional) -----------------------
    if plot:
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.set_xlabel('Width (mm)')
        ax.set_ylabel('Height (mm)')
        ax.set_title('Dumbbell-shaped Shear-wall \nI Section Different flange with Rebars (numbers shown)')
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

        max_dim = max(bfT, bfB, h) + 50
        ax.set_xlim(-max_dim/2, max_dim/2)
        ax.set_ylim(-50.0, max_dim)
        ax.set_aspect('equal')
        plt.show()

    ELE_MASS = CONCRETE_DENSITY * AREA   # kg/mm

    return h, ELE_MASS
