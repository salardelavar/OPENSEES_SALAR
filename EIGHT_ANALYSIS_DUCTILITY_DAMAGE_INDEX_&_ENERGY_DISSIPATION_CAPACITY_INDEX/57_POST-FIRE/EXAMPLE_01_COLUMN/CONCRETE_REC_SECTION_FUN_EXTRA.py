"""
import openseespy.opensees as ops
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# --------- Section Properties (in mm) ---------
# Rectabgular-beam dimensions
Bsec = 400.0       # width
Hsec  = 600.0      # Theight
cover = 50

fc = 30.0               # [N/mm²] Concrete Compressive Strength
Kc = 1.24               # ratio of confined to unconfined concrete strength
STEEL_TYPE = 'INELASTIC'
# --------- Initialize OpenSees ---------
ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 3)
ops.uniaxialMaterial('Elastic', 1, 100000.0)

# Create a section tag for the fiber section
section_tag = 1
"""
def CONCRETE_REC_SECTION_FUN(secTag, STEEL_TYPE, fc, Kfc,
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
    
    THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
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

    coreTag, coverTag, steelTag = secTag + 100, secTag + 200, secTag + 300
    if STEEL_TYPE == 'ELASTIC':
        ops.uniaxialMaterial('Elastic', steelTag, Es)
    if STEEL_TYPE == 'INELASTIC':
        pinchX = 0.8  # Pinching factor in X direction
        pinchY = 0.5  # Pinching factor in Y direction
        damage1 = 0.0  # Damage due to ductility
        damage2 = 0.0  # Damage due to energy
        beta = 0.1  # Stiffness degradation parameter
        ops.uniaxialMaterial('Hysteretic', steelTag,
                             fy, ey,                     # YIELD STRENGTH
                             fu, esu,                    # UTIMATE STRENGTH
                             0.2 * fu, 1.1 * esu,        # SOFTENING STRENGTH
                             -fy, -ey,
                             -fu, -esu,
                             -0.2 * fu, -1.1 * esu,
                             pinchX, pinchY,
                             damage1, damage2, beta)
        
    #ops.uniaxialMaterial('Concrete01', coreTag, fcC, ec0C, fcUC, ecuC)  # Core concrete (confined)
    #ops.uniaxialMaterial('Concrete01', coverTag, fcU, ec0U, fcUU, ecuU) # Cover concrete (unconfined)
    ops.uniaxialMaterial('Concrete02', coreTag, fcC, ec0C, fcUC, ecuC, LambdaC, ftC, EtsC) # build core concrete (confined)
    ops.uniaxialMaterial('Concrete02', coverTag, fcU, ec0U, fcUU, ecuU, LambdaU, ftU, EtsU) # build cover concrete (unconfined)

    # -------------------- Section area -------------------------
    AREA = Bsec * Hsec
    print(f"Total Section Area = {AREA:.2f} mm²")

    # -------------------- Fiber section -------------------------
    ops.section('Fiber', secTag)

    def add_rect(mat, y_bot, y_top, x_left, x_right):
        ops.patch('rect', mat, nFib, nFib,
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
    rebars = [
        (25.0,  Hsec/2 - cover, -Bsec/2 + cover),   # 1
        (25.0,  Hsec/2 - cover,  Bsec/2 - cover),   # 2
        (18.0,  Hsec/2 - cover,  0.0),              # 3
        (18.0,  0.0,  -Bsec/2 + cover),             # 4
        (18.0,  0.0,   Bsec/2 - cover),             # 5
        (25.0, -Hsec/2 + cover, -Bsec/2 + cover),   # 6
        (25.0, -Hsec/2 + cover,  Bsec/2 - cover),   # 7
        (18.0, -Hsec/2 + cover,  0.0)               # 8
    ]

    for dia, y, x in rebars:
        area = np.pi * dia**2 / 4.0            # mm²
        ops.fiber(x, y, area, steelTag)        # (area, material, x, y)

    # -------------------- Plot (optional) -----------------------
    if plot:
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.set_xlabel('Width (mm)')
        ax.set_ylabel('Height (mm)')
        ax.set_title('Rectangular Section with Rebars (numbers shown)')
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

