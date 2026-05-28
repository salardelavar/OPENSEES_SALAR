def CONCRETE_RECTANGULAR_HOLLOW_SECTION_FUN(secTag, STEEL_TYPE, fc, Kfc,
                           b_o, h_o, b_i, h_i,
                           nFib, CONCRETE_DENSITY,
                           plot=True):
    """
    Create a Hollow Box confined‑concrete fiber section (OpenSees) and,
    if PLOT=True, draw a simple 2‑D sketch of the section.

    Parameters
    ----------
    secTag          : int   – section identifier
    b_o, h_o       : float – outer width and height  (mm)
    b_i, h_i       : float – inner width and height  (mm)
    cover           : float – concrete cover (mm)
    STEEL_TYPE      : str   – 'ELASTIC' or 'INELASTIC'
    fc              : float – unconfined concrete compressive strength (MPa, negative)
    Kfc             : float – ratio of confined to unconfined concrete strength
    CONCRETE_DENSITY: float – ρc (kg/m³) – will be converted to N·s²/mm³
    PLOT            : bool  – draw the section if True
    
    THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
    """
    import openseespy.opensees as ops
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    import numpy as np
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
    ops.uniaxialMaterial('Concrete02', coreTag, fcC, ec0C, fcUC, ecuC, Lambda, ftC, EtsC) # build core concrete (confined)
    ops.uniaxialMaterial('Concrete02', coverTag, fcU, ec0U, fcUU, ecuU, Lambda, ftU, EtsU) # build cover concrete (unconfined)

    # -------------------- Section area -------------------------
    AREA = (b_o * h_o) - (b_i * h_i)
    print(f"Total Section Area = {AREA:.2f} mm²")

    # -------------------- Fiber section -------------------------
    ops.section('Fiber', secTag)

    def add_rect(mat, y_bot, y_top, x_left, x_right):
        ops.patch('rect', mat, nFib, nFib,
                  x_left, y_bot, x_right, y_top)

    T_h = h_o - h_i
    T_b = b_o - b_i
    # ===== Outer cover =====
    # Top‑bottom edges
    add_rect(coverTag,  h_o/2-T_h,  h_o/2, -b_o/2,  b_o/2)   # top edge
    add_rect(coverTag, -h_o/2+T_h, -h_o/2, -b_o/2,  b_o/2)   # bottom edge
    # Left‑right edges
    add_rect(coverTag, -h_o/2+T_h,  h_o/2-T_h, -b_o/2, -b_o/2+T_b)   # left edge
    add_rect(coverTag, -h_o/2+T_h,  h_o/2-T_h,  b_o/2,  b_o/2-T_b)   # right edge

    # -------------------- Rebars -------------------------------
    # (diameter [mm], y‑coord [mm], x‑coord [mm])
    rebars = [
        # (diameter[mm], y[mm], x[mm])
        (32.0,  h_o/2-0.25*T_h,  0.0),      # 1
        (32.0, -h_o/2+0.25*T_h,  0.0),      # 2
        (25.0,  h_o/2-0.25*T_h, -b_o/4),    # 3
        (25.0,  h_o/2-0.25*T_h,  b_o/4),    # 4
        (25.0, -h_o/2+0.25*T_h, -b_o/4),    # 5
        (25.0, -h_o/2+0.25*T_h,  b_o/4),    # 6
        (32.0,  h_i/2+0.25*T_h,  -b_o/2 + 0.25*T_b), # 7
        (32.0, -h_i/2-0.25*T_h,  -b_o/2 + 0.25*T_b), # 8
        (25.0,  h_i/4, -b_o/2 + 0.25*T_b),  # 9
        (25.0,  -h_i/4, -b_o/2 + 0.25*T_b), # 10
        (25.0, h_i/8, -b_o/2 + 0.25*T_b),   # 11
        (25.0, -h_i/8,  -b_o/2 + 0.25*T_b), # 12
        
        (32.0,  h_i/2+0.25*T_h,  +b_o/2 - 0.25*T_b), # 13
        (32.0, -h_i/2-0.25*T_h,  +b_o/2 - 0.25*T_b), # 14
        (25.0,  h_i/4, +b_o/2 - 0.25*T_b),  # 15
        (25.0,  -h_i/4, +b_o/2 - 0.25*T_b), # 16
        (25.0, h_i/8, +b_o/2 - 0.25*T_b),   # 17
        (25.0, -h_i/8,  +b_o/2 - 0.25*T_b), # 18
        (25.0,  h_o/2-0.25*T_h, -b_o/3),    # 19
        (25.0,  h_o/2-0.25*T_h,  b_o/3),    # 20
        (25.0,  h_o/2-0.25*T_h, -b_o/10),   # 21
        (25.0,  h_o/2-0.25*T_h,  b_o/10),   # 22
        (25.0, -h_o/2+0.25*T_h, -b_o/3),    # 23
        (25.0, -h_o/2+0.25*T_h,  b_o/3),    # 24
        (25.0, -h_o/2+0.25*T_h, -b_o/10),   # 25
        (25.0, -h_o/2+0.25*T_h,  b_o/10),   # 26
        (25.0,  h_i/3, +b_o/2 - 0.25*T_b),  # 27
        (25.0,  -h_i/3, +b_o/2 - 0.25*T_b), # 28
        (32.0, 0.0, -b_o/2 + 0.25*T_b),     # 29
        (32.0, 0.0,  +b_o/2 - 0.25*T_b),    # 30
        (25.0,  h_i/2 - 50, -b_o/2 + 0.25*T_b),  # 31
        (25.0,  -h_i/2 +50, -b_o/2 + 0.25*T_b),  # 32
        (25.0,  h_i/2 - 50, +b_o/2 - 0.25*T_b),  # 33
        (25.0,  -h_i/2 + 50, +b_o/2 - 0.25*T_b), # 34
    ]
    
    for dia, y, x in rebars:
        area = np.pi * dia**2 / 4.0          # mm²
        ops.fiber(x, y, area, steelTag)     # add to OpenSees model
    
    # -------------------- Plot (optional) -----------------------
    if plot:
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.set_xlabel('Width (mm)')
        ax.set_ylabel('Height (mm)')
        ax.set_title('Hollow rectangular section with rebars')
        ax.grid(True, ls='--', alpha=0.4)
    
        # Outer cover – light gray
        outer = patches.Rectangle((-b_o/2, -h_o/2), b_o, h_o,
                                  linewidth=1.5, edgecolor='black',
                                  facecolor='lightgray')
        ax.add_patch(outer)

        # Inner hole – white
        inner = patches.Rectangle((-b_i/2, -h_i/2), b_i, h_i,
                                  linewidth=1.5, edgecolor='black',
                                  facecolor='white')
        ax.add_patch(inner)

        # Rebars – red circles + numbers
        for idx, (dia, y, x) in enumerate(rebars, start=1):
            circ = patches.Circle((x, y), dia/2,
                                  edgecolor='red', facecolor='red')
            ax.add_patch(circ)
            ax.text(x, y + dia/2 + 2, f'{idx}',
                    color='purple', fontsize=6,
                    ha='center', va='bottom',
                    fontweight='bold')
    
        # Set limits
        lim = max(b_o, h_o) * 0.6 + 20
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_aspect('equal')
        plt.show()
    
    ELE_MASS = CONCRETE_DENSITY * AREA # kg/mm
    
    return h_o, ELE_MASS
