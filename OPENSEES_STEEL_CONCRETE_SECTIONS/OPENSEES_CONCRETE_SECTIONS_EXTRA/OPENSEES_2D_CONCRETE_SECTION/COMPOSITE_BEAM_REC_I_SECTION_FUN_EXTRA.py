def COMPOSITE_BEAM_REC_I_SECTION_FUN_EXTRA(secTag, Hsec, Bsec,
                              STEEL_TYPE, nFib, fc, Kfc,
                              bf, tf, h, tw,
                              STEEL_DENSITY,
                              CONCRETE_DENSITY, PLOT):
 
    """
    Create a Composite Beam-Slab fiber section (OpenSees) and,
    if PLOT=True, draw a simple 2‑D sketch of the section.

    Parameters
    ----------
    secTag          : int   – section identifier
    HSec, BSec      : float – height and width  (mm)
    bf              : float – Flange width
    tf              : float – Flange thickness
    h               : float – total I height
    tw              : float – Web thickness
    STEEL_TYPE      : str   – 'ELASTIC' or 'INELASTIC'
    fc              : float – unconfined concrete compressive strength (MPa, negative)
    Kfc             : float – ratio of confined to unconfined concrete strength
    CONCRETE_DENSITY: float – ρc (kg/m³) – will be converted to N·s²/mm³
    STEEL_DENSITY   : float – ρc (kg/m³) – will be converted to N·s²/mm³
    PLOT            : bool  – draw the section if True
    
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
        ops.uniaxialMaterial('Elastic', steelTag, Es)              # REBAR
        ops.uniaxialMaterial('Elastic', steelITag, EsI)            # I SECION
    if STEEL_TYPE == 'INELASTIC':   
        pinchX = 0.8   # Pinching factor in X direction
        pinchY = 0.5   # Pinching factor in Y direction
        damage1 = 0.0  # Damage due to ductility
        damage2 = 0.0  # Damage due to energy
        beta = 0.1     # Stiffness degradation parameter
        # REBAR
        ops.uniaxialMaterial('Hysteretic', steelTag,
                             fy, ey,
                             fu, esu,
                             0.2*fu, 1.1*esu,
                             -fy, -ey,
                             -fu, -esu,
                             -0.2*fu, -1.1*esu,
                             pinchX, pinchY,
                             damage1, damage2, beta)
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
    # --------- Define Fibers for Flanges and Web ---------
    def add_rectangle_patch(Mat_tag, y_bottom, y_top, x_left, x_right, num_y_divisions=nFib, num_x_divisions=nFib):
        """Helper function to add rectangular patches of fibers."""
        ops.patch('rect', Mat_tag, num_y_divisions, num_x_divisions, x_left, y_bottom, x_right, y_top)
    
    # SLAB SECTION
    add_rectangle_patch(coreTag, 0.0, Hsec, -Bsec/2, Bsec/2)
    #%% I SECTION
    # Top flange
    add_rectangle_patch(steelITag, -tf, 0.0, -bf/2, bf/2)
    
    # Bottom flange
    add_rectangle_patch(steelITag, -h, -h + tf, -bf/2, bf/2)
    
    # Web
    add_rectangle_patch(steelITag, -h + tf, -tf, -tw/2, tw/2)
    
    # Calculate mass
    # Area of main I-beam components
    area_top_flange = bf * tf
    area_bottom_flange = bf * tf
    area_web = tw * (h - 2 * tf) # Height of web is total height minus thicknesses of top and bottom flanges
    AREA = area_top_flange + area_bottom_flange + area_web
    print(f"Area of Top Flange: {area_top_flange:.2f} mm^2")
    print(f"Area of Bottom Flange: {area_bottom_flange:.2f} mm^2")
    print(f"Area of Web: {area_web:.2f} mm^2")
    print(f"Total I Section Area: {AREA:.2f} mm^2")

    RD = 16.0      # [mm] Rebar Diameter
    DIST = 150.0   # [mm] Rebar Distance
    
    # -------------------- Rebars -------------------------------
    # (diameter [mm], y‑coord [mm], x‑coord [mm])
    rebars = [
        (RD,  Hsec/4, -Bsec/2 + 1*DIST),     # 1
        (RD,  Hsec/4, +Bsec/2 - 1*DIST),     # 2
        (RD,  Hsec/4, -Bsec/2 + 2*DIST),     # 3
        (RD,  Hsec/4, +Bsec/2 - 2*DIST),     # 4
        (RD,  Hsec/4, -Bsec/2 + 3*DIST),     # 5
        (RD,  Hsec/4, +Bsec/2 - 3*DIST),     # 6
        (RD,  Hsec/4, -Bsec/2 + 4*DIST),     # 7
        (RD,  Hsec/4, +Bsec/2 - 4*DIST),     # 8
        (RD,  Hsec/4, 0.0),                  # 9
        
        (RD,  3*Hsec/4, -Bsec/2 + 1*DIST),   # 10
        (RD,  3*Hsec/4, +Bsec/2 - 1*DIST),   # 11
        (RD,  3*Hsec/4, -Bsec/2 + 2*DIST),   # 12
        (RD,  3*Hsec/4, +Bsec/2 - 2*DIST),   # 13
        (RD,  3*Hsec/4, -Bsec/2 + 3*DIST),   # 14
        (RD,  3*Hsec/4, +Bsec/2 - 3*DIST),   # 15
        (RD,  3*Hsec/4, -Bsec/2 + 4*DIST),   # 16
        (RD,  3*Hsec/4, +Bsec/2 - 4*DIST),   # 17
        (RD,  3*Hsec/4, 0.0),                # 18

    ]

    for dia, y, x in rebars:
        area = np.pi * dia**2 / 4.0          # mm²
        ops.fiber(x, y, area, steelTag)      # add to OpenSees model

    # -------------------- Plot (optional) -----------------------
    if PLOT:
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.set_xlabel('Width (mm)')
        ax.set_ylabel('Height (mm)')
        ax.set_title('Composite Beam Section with Rebars (numbers shown)')
        ax.grid(True, ls='--', alpha=0.5)

        # Concrete geometry (light gray)
        geometry = [
            {"x": -Bsec/2, "y": 0.0, "w": Bsec, "h": Hsec},           # Slab
            {"x": -bf/2, "y": -tf, "w": bf, "h": tf},                 # Top flange
            {"x": -bf/2, "y": -h, "w": bf, "h": tf},                  # Bottom flange
            {"x": -tw/2, "y": -h+tf,    "w": tw, "h": h-2*tf}         # Web
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

        max_dim = max(bf, h + Hsec, Bsec) + 50
        ax.set_xlim(-max_dim/2, max_dim/2)
        ax.set_ylim(-max_dim/2, max_dim/2)
        ax.set_aspect('equal')
        plt.show()

    MASS = Hsec * Bsec * CONCRETE_DENSITY + AREA * STEEL_DENSITY # kg/mm
    
    return Hsec+h, MASS
