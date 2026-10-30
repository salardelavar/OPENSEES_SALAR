def CONFINED_CONCRETE_SECTION(secTag, HSec, BSec, cover,
                              STEEL_TYPE, fc, Kc,
                              CONCRETE_DENSITY, PLOT):
 
    """
    Create a rectangular confined‑concrete fiber section (OpenSees) and,
    if PLOT=True, draw a simple 2‑D sketch of the section.

    Parameters
    ----------
    secTag          : int   – section identifier
    HSec, BSec     : float – height and width  (mm)
    cover           : float – concrete cover (mm)
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
    Kfc = 1.3;			      # ratio of confined to unconfined concrete strength
    fcC = Kfc*fcU             # [N/mm²] Concrete Compressive Strength
    Ec = 4700 * np.sqrt(-fcC) # [N/mm^2] Concrete Elastic Modulus
    ec0C = 2*fcC/Ec           # [mm/mm] Concrete Compressive Strain
    fcUC = 0.65*fcC           # [N/mm²] Concrete Compressive Ultimate Strength
    ecuC = 15*ec0C            # [mm/mm] Concrete Compressive Ultimate Strain
    LambdaC = 0.1;	          # ratio between unloading slope
    # tensile-strength properties
    ftC = 0.7 * np.sqrt(-fcC)  # [N/mm²] tensile strength +tension
    ftU = 0.7 * np.sqrt(-fcU)  # [N/mm²] tensile strength +tension
    EtsC = ftC/np.abs(ec0C);   # [N/mm²] tension softening stiffness
    EtsU = ftU/np.abs(ec0U);   # [N/mm²] tension softening stiffness
    
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
        
    #ops.uniaxialMaterial('Concrete01', coreTag, fcC, ec0C, fcUC, ecuC)                     # Core concrete (confined)
    #ops.uniaxialMaterial('Concrete01', coverTag, fcU, ec0U, fcUU, ecuU)                    # Cover concrete (unconfined)
    ops.uniaxialMaterial('Concrete02', coreTag, fcC, ec0C, fcUC, ecuC, LambdaC, ftC, EtsC)  # build core concrete (confined)
    ops.uniaxialMaterial('Concrete02', coverTag, fcU, ec0U, fcUU, ecuU, LambdaU, ftU, EtsU) # build cover concrete (unconfined)
    
    # Some variables derived from the parameters
    y1 = HSec / 2.0
    z1 = BSec / 2.0
    NUMFIBERS = 20  # Number of layers for each fiber
    
    ops.section('Fiber', secTag)
    # Create the concrete core fibers
    ops.patch('rect', coreTag, NUMFIBERS, NUMFIBERS, cover - y1, cover - z1, y1 - cover, z1 - cover)
    
    # Create the concrete cover fibers (top, bottom, left, right)
    ops.patch('rect', coverTag, NUMFIBERS, 10, -y1, z1 - cover, y1, z1)
    ops.patch('rect', coverTag, NUMFIBERS, 10, -y1, -z1, y1, cover - z1)
    ops.patch('rect', coverTag, NUMFIBERS, 10, -y1, cover - z1, cover - y1, z1 - cover)
    ops.patch('rect', coverTag, NUMFIBERS, 10, y1 - cover, cover - z1, y1, z1 - cover)
    REBAR = 25.0 # [mm] Rebar Diameter
    # -------------------- Rebars -------------------------------
    # (diameter [mm], y‑coord [mm], x‑coord [mm])
    rebars = [
        (REBAR,  HSec/2 - cover, -BSec/2 + cover),    # 1
        (REBAR,  HSec/2 - cover, BSec/2 - cover),     # 2
        (REBAR,  HSec/2 - cover, 0.0),                # 3
        (REBAR,  -HSec/2 + cover, -BSec/2 + cover),   # 4
        (REBAR,  -HSec/2 + cover, 0.0),               # 5
        (REBAR,  -HSec/2 + cover, +BSec/2 - cover),   # 6
        (16.0, HSec/4 - cover, -BSec/2 + cover),     # 7
        (16.0, HSec/4 - cover, +BSec/2 - cover),     # 8
        (16.0,  -HSec/4 + cover, -BSec/2 + cover),   # 9
        (16.0,  -HSec/4 + cover, +BSec/2 - cover),   # 10
        (14.0,  0.0, -BSec/2 + cover),               # 11
        (14.0,  0.0, +BSec/2 - cover),               # 12
    ]

    for dia, y, x in rebars:
        area = np.pi * dia**2 / 4.0          # mm²
        ops.fiber(x, y, area, steelTag)      # add to OpenSees model

    
    if PLOT:
        fig, ax = plt.subplots(figsize=(6, 6))
        # concrete block
        rect = plt.Rectangle((-z1, -y1), BSec, HSec,
                             facecolor='lightgray', edgecolor='k', lw=1.5, label='Unconfined')
        ax.add_patch(rect)

        # cover outline (just for visual aid)
        inner = plt.Rectangle((-z1 + cover, -y1 + cover),
                              BSec - 2 * cover, HSec - 2 * cover,
                              facecolor='gray', label='Confined')#, edgecolor='k', ls='--'
        ax.add_patch(inner)

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

        ax.set_aspect('equal')
        ax.set_xlabel('Width  (mm)')
        ax.set_ylabel('Height (mm)')
        ax.set_title(f'Confined Concrete Section – tag {secTag}')
 
        max_dim = max(BSec, HSec) + 50
        ax.set_xlim(-max_dim/2, max_dim/2)
        ax.set_ylim(-max_dim/2, max_dim/2)
        ax.grid(True, ls=':', alpha=0.5)
        ax.legend()
        plt.show()
        
    # Calculate mass
    MASS = HSec * BSec * CONCRETE_DENSITY     

    return HSec, MASS
