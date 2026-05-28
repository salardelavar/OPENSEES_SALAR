def CONCRETE_SHEAR_WALL_SECTION_FUN(secTag, STEEL_TYPE, fc, Kfc,
                           Bsec, Hsec, cover,
                           nFib, CONCRETE_DENSITY,
                           plot=True):
    """
    Create a rectangular confined‑concrete fiber section (OpenSees).

    Parameters
    ----------
    secTag, matTag   : int   – identifiers (only secTag is used in this routine)
    STEEL_TYPE       : str   – 'ELASTIC' or 'INELASTIC'
    fc               : float – unconfined concrete compressive strength (MPa, positive)
    Kc               : float – ratio of confined to unconfined concrete strength
    Bsec, Hsec       : float – width and height of the rectangle (mm)
    nFib             : int   – number of fibers per direction for each patch
    CONCRETE_DENSITY : float – concrete density in kg/mm³  (≈ 2.5e‑9 kg/mm³)
    cover            : float – concrete cover for rebars (mm)
    plot             : bool  – draw a sketch if True

    Returns
    -------
    AREA, ELE_MASS   : float – concrete area (mm²) and mass per unit length (kg/mm)
    
    THIS PYTHON SCRIPT WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
    """
    import numpy as np
    import openseespy.opensees as ops
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    
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
    
    import openseespy.opensees as ops
    
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

    # -------------------- Section area -------------------------
    AREA = Bsec * Hsec
    print(f"Total Section Area = {AREA:.2f} mm²")

    # -------------------- Fiber section -------------------------
    ops.section('Fiber', secTag)
    nFibZ = 10
    def add_rect(mat, y_bot, y_top, x_left, x_right):
        ops.patch('rect', mat, nFibZ, nFib,
                  x_left, y_bot, x_right, y_top)

    # Concrete patch – whole rectangle (cover concrete)
    add_rect(coverTag,
             -Hsec/2,        # y‑bottom
              Hsec/2,        # y‑top
             -Bsec/2,        # x‑left
              Bsec/2)        # x‑right

    # -------------------- Rebars -------------------------------
    #   (diameter [mm], y‑coord [mm], x‑coord [mm])
    #   The coordinates below are the *centre* of each bar.
    #   You can modify them or the cover value as needed.
    REBAR_DIA = 16.0 # [mm] Rebar Diameter
    DIST = 220.0     # [mm] Rebar Distance
    rebars = [
        (REBAR_DIA, -Hsec/2 + cover,  Bsec/2 - cover),           # 1
        (REBAR_DIA,  Hsec/2 - cover,  Bsec/2 - cover),           # 2
        (REBAR_DIA, -Hsec/2 + cover + 1*DIST,  Bsec/2 - cover),   # 3
        (REBAR_DIA,  Hsec/2 - cover - 1*DIST,  Bsec/2 - cover),   # 4
        (REBAR_DIA, -Hsec/2 + cover + 2*DIST,  Bsec/2 - cover),   # 5
        (REBAR_DIA,  Hsec/2 - cover - 2*DIST,  Bsec/2 - cover),   # 6
        (REBAR_DIA, -Hsec/2 + cover + 3*DIST,  Bsec/2 - cover),   # 7
        (REBAR_DIA,  Hsec/2 - cover - 3*DIST,  Bsec/2 - cover),   # 8
        (REBAR_DIA, -Hsec/2 + cover + 4*DIST,  Bsec/2 - cover),   # 9
        (REBAR_DIA,  Hsec/2 - cover - 4*DIST,  Bsec/2 - cover),   # 10
        (REBAR_DIA, -Hsec/2 + cover + 5*DIST,  Bsec/2 - cover),   # 11
        (REBAR_DIA,  Hsec/2 - cover - 5*DIST,  Bsec/2 - cover),   # 12
        (REBAR_DIA, -Hsec/2 + cover + 6*DIST,  Bsec/2 - cover),   # 13
        (REBAR_DIA,  Hsec/2 - cover - 6*DIST,  Bsec/2 - cover),   # 14
        (REBAR_DIA, -Hsec/2 + cover + 7*DIST,  Bsec/2 - cover),   # 15
        (REBAR_DIA,  Hsec/2 - cover - 7*DIST,  Bsec/2 - cover),   # 16
        (REBAR_DIA, -Hsec/2 + cover + 8*DIST,  Bsec/2 - cover),   # 17
        (REBAR_DIA,  Hsec/2 - cover - 8*DIST,  Bsec/2 - cover),   # 18
        (REBAR_DIA, -Hsec/2 + cover + 9*DIST,  Bsec/2 - cover),   # 19
        (REBAR_DIA,  Hsec/2 - cover - 9*DIST,  Bsec/2 - cover),   # 20
        (REBAR_DIA, -Hsec/2 + cover + 10*DIST,  Bsec/2 - cover),   # 21
        (REBAR_DIA,  Hsec/2 - cover - 10*DIST,  Bsec/2 - cover),   # 22
        (REBAR_DIA, -Hsec/2 + cover + 11*DIST,  Bsec/2 - cover),   # 23
        (REBAR_DIA,  Hsec/2 - cover - 11*DIST,  Bsec/2 - cover),   # 24
        (REBAR_DIA, -Hsec/2 + cover + 12*DIST,  Bsec/2 - cover),   # 25
        (REBAR_DIA,  Hsec/2 - cover - 12*DIST,  Bsec/2 - cover),   # 26
        (REBAR_DIA,  0.0,  Bsec/2 - cover),                       # 27 - x=0.0
        
        (REBAR_DIA, -Hsec/2 + cover, 0.0),           # 28
        (REBAR_DIA,  Hsec/2 - cover, 0.0),           # 29
        (REBAR_DIA, -Hsec/2 + cover + 1*DIST, 0.0),   # 30
        (REBAR_DIA,  Hsec/2 - cover - 1*DIST, 0.0),   # 31
        (REBAR_DIA, -Hsec/2 + cover + 2*DIST, 0.0),   # 32
        (REBAR_DIA,  Hsec/2 - cover - 2*DIST, 0.0),   # 33
        (REBAR_DIA, -Hsec/2 + cover + 3*DIST, 0.0),   # 34
        (REBAR_DIA,  Hsec/2 - cover - 3*DIST, 0.0),   # 35
        (REBAR_DIA, -Hsec/2 + cover + 4*DIST, 0.0),   # 36
        (REBAR_DIA,  Hsec/2 - cover - 4*DIST, 0.0),   # 37
        (REBAR_DIA, -Hsec/2 + cover + 5*DIST, 0.0),   # 38
        (REBAR_DIA,  Hsec/2 - cover - 5*DIST, 0.0),   # 39
        (REBAR_DIA, -Hsec/2 + cover + 6*DIST, 0.0),   # 40
        (REBAR_DIA,  Hsec/2 - cover - 6*DIST, 0.0),   # 41
        (REBAR_DIA, -Hsec/2 + cover + 7*DIST, 0.0),   # 42
        (REBAR_DIA,  Hsec/2 - cover - 7*DIST, 0.0),   # 43
        (REBAR_DIA, -Hsec/2 + cover + 8*DIST, 0.0),   # 44
        (REBAR_DIA,  Hsec/2 - cover - 8*DIST, 0.0),   # 45
        (REBAR_DIA, -Hsec/2 + cover + 9*DIST, 0.0),   # 46
        (REBAR_DIA,  Hsec/2 - cover - 9*DIST, 0.0),   # 47
        (REBAR_DIA, -Hsec/2 + cover + 10*DIST, 0.0),   # 48
        (REBAR_DIA,  Hsec/2 - cover - 10*DIST, 0.0),   # 49
        (REBAR_DIA, -Hsec/2 + cover + 11*DIST, 0.0),   # 50
        (REBAR_DIA,  Hsec/2 - cover - 11*DIST, 0.0),   # 51
        (REBAR_DIA, -Hsec/2 + cover + 12*DIST, 0.0),   # 52
        (REBAR_DIA,  Hsec/2 - cover - 12*DIST, 0.0),   # 53
        (REBAR_DIA,  0.0, 0.0),                      # 54 - x=0.0        
        
        (REBAR_DIA, -Hsec/2 + cover,  -Bsec/2 + cover),           # 55
        (REBAR_DIA,    Hsec/2 - cover,  -Bsec/2 + cover),           # 56
        (REBAR_DIA,   -Hsec/2 + cover + 1*DIST,  -Bsec/2 + cover),   # 57
        (REBAR_DIA,    Hsec/2 - cover - 1*DIST,  -Bsec/2 + cover),   # 58
        (REBAR_DIA,   -Hsec/2 + cover + 2*DIST,  -Bsec/2 + cover),   # 59
        (REBAR_DIA,    Hsec/2 - cover - 2*DIST,  -Bsec/2 + cover),   # 60
        (REBAR_DIA,   -Hsec/2 + cover + 3*DIST,  -Bsec/2 + cover),   # 61
        (REBAR_DIA,    Hsec/2 - cover - 3*DIST,  -Bsec/2 + cover),   # 62
        (REBAR_DIA,   -Hsec/2 + cover + 4*DIST,  -Bsec/2 + cover),   # 63
        (REBAR_DIA,    Hsec/2 - cover - 4*DIST,  -Bsec/2 + cover),   # 64
        (REBAR_DIA,   -Hsec/2 + cover + 5*DIST,  -Bsec/2 + cover),   # 65
        (REBAR_DIA,    Hsec/2 - cover - 5*DIST,  -Bsec/2 + cover),   # 66
        (REBAR_DIA,   -Hsec/2 + cover + 6*DIST,  -Bsec/2 + cover),   # 67
        (REBAR_DIA,    Hsec/2 - cover - 6*DIST,  -Bsec/2 + cover),   # 68
        (REBAR_DIA,   -Hsec/2 + cover + 7*DIST,  -Bsec/2 + cover),   # 69
        (REBAR_DIA,    Hsec/2 - cover - 7*DIST,  -Bsec/2 + cover),   # 70
        (REBAR_DIA,   -Hsec/2 + cover + 8*DIST,  -Bsec/2 + cover),   # 71
        (REBAR_DIA,    Hsec/2 - cover - 8*DIST,  -Bsec/2 + cover),   # 72
        (REBAR_DIA,   -Hsec/2 + cover + 9*DIST,  -Bsec/2 + cover),   # 73
        (REBAR_DIA,    Hsec/2 - cover - 9*DIST,  -Bsec/2 + cover),   # 74
        (REBAR_DIA,   -Hsec/2 + cover + 10*DIST,  -Bsec/2 + cover),   # 75
        (REBAR_DIA,    Hsec/2 - cover - 10*DIST,  -Bsec/2 + cover),   # 76
        (REBAR_DIA,   -Hsec/2 + cover + 11*DIST,  -Bsec/2 + cover),   # 77
        (REBAR_DIA,    Hsec/2 - cover - 11*DIST,  -Bsec/2 + cover),   # 78
        (REBAR_DIA,   -Hsec/2 + cover + 12*DIST,  -Bsec/2 + cover),   # 79
        (REBAR_DIA,    Hsec/2 - cover - 12*DIST,  -Bsec/2 + cover),   # 80
        (REBAR_DIA,    0.0,  -Bsec/2 + cover),                      # 81 - x=0.0
    ]

    for dia, y, x in rebars:
        area = np.pi * dia**2 / 4.0            # mm²
        ops.fiber(x, y, area, steelTag)        # (area, material, x, y)

    # -------------------- Plot (optional) -----------------------
    if plot:
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.set_xlabel('Width (mm)')
        ax.set_ylabel('Height (mm)')
        ax.set_title('Rectangular Shear Wall Section with Rebars (numbers shown)')
        ax.grid(True, ls='--', alpha=0.5)

        # Concrete block (light gray)
        rect = patches.Rectangle((-Bsec/2, -Hsec/2),
                                 Bsec, Hsec,
                                 linewidth=1.5,
                                 edgecolor='black',
                                 facecolor='lightgray')
        ax.add_patch(rect)

        # Rebars – red circles + numbers
        for i, (dia, y, x) in enumerate(rebars, start=1):
            circ = patches.Circle((x, y),
                                 radius=dia/2,
                                 edgecolor='red',
                                 facecolor='red',
                                 linewidth=1.5)
            ax.add_patch(circ)

            # label a little above the bar
            ax.text(x, y + dia/2 + 4,
                    f'{i}',
                    color='purple',
                    fontsize=8,
                    ha='center',
                    va='bottom',
                    fontweight='bold')

        max_dim = max(Bsec, Hsec) + 50
        ax.set_xlim(-max_dim/2, max_dim/2)
        ax.set_ylim(-max_dim/2, max_dim/2)
        ax.set_aspect('equal')
        plt.show()

    ELE_MASS = CONCRETE_DENSITY * AREA   # kg/mm
    return Hsec, ELE_MASS

