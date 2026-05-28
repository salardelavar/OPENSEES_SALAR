"""
import openseespy.opensees as ops
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# --------- Section Properties (in mm) ---------
# I-beam dimensions
bf = 1000.0     # Flange width
tf = 100.0      # Flange thickness
h  = 600.0      # Total web height
tw = 200.0      # Web thickness

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
def CONCRETE_I_SECTION_FUN(secTag, matTag, STEEL_TYPE, fc, Kfc,
                           bf, tf, h, tw,
                           nFib, CONCRETE_DENSITY,
                           plot=True):
    """
    Create a I confined‑concrete fiber section (OpenSees) and,
    if PLOT=True, draw a simple 2‑D sketch of the section.

    Parameters
    ----------
    secTag          : int   – section identifier
    h, bf           : float – height and width  (mm)
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
        ops.uniaxialMaterial('Hysteretic', steelTag, fy, ey, fu, esu, 0.2*fu, 1.1*esu, -fy, -ey, -fu, -esu, -0.2*fu, -1.1*esu, pinchX, pinchY, damage1, damage2, beta)
        # INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Hysteretic_Material
        
    #ops.uniaxialMaterial('Concrete01', coreTag, fcC, ec0C, fcUC, ecuC)  # Core concrete (confined)
    #ops.uniaxialMaterial('Concrete01', coverTag, fcU, ec0U, fcUU, ecuU) # Cover concrete (unconfined)
    ops.uniaxialMaterial('Concrete02', coreTag, fcC, ec0C, fcUC, ecuC, Lambda, ftC, EtsC) # build core concrete (confined)
    ops.uniaxialMaterial('Concrete02', coverTag, fcU, ec0U, fcUU, ecuU, Lambda, ftU, EtsU) # build cover concrete (unconfined)

    # -------------------- Section area -------------------------
    area_top_flange = bf * tf
    area_bottom_flange = bf * tf
    area_web        = tw * (h-2*tf)
    AREA = area_top_flange + area_bottom_flange + area_web
    print(f"Total Section Area = {AREA:.2f} mm²")

    # -------------------- Fiber section -------------------------
    ops.section('Fiber', secTag)

    def add_rect(mat, y_bot, y_top, x_left, x_right):
        ops.patch('rect', mat, nFib, nFib,
                  x_left, y_bot, x_right, y_top)

    # Concrete patches
    add_rect(coverTag,  h/2 - tf,  h/2, -bf/2,  bf/2)     # Top flange
    add_rect(coverTag, -h/2 + tf, -h/2, -bf/2,  bf/2)     # Bottom flange
    add_rect(coverTag, -h/2 + tf,  h/2 - tf, -tw/2, tw/2) # Web
    # -------------------- Rebars -------------------------------
    # (diameter [mm], y‑coord [mm], x‑coord [mm])
    rebars = [
        (25.0,  h/2 - tf/2, -bf/3),   # 1
        (25.0,  h/2 - tf/2,  bf/3),   # 2
        (25.0,  h/2 - tf/2, -bf/6),   # 3
        (25.0,  h/2 - tf/2, bf/6),    # 4
        (18.0, h/4,      -tw/4),      # 5
        (18.0, h/4,       tw/4),      # 6
        (18.0,  -h/4,      -tw/4),    # 7
        (18.0,  -h/4,       tw/4),    # 8
        (18.0, 0.0, -bf/20),          # 9
        (18.0, 0.0,  bf/20),          # 10
        (25.0,  -h/2 + tf/2, -bf/3),  # 11
        (25.0,  -h/2 + tf/2,  bf/3),  # 12
        (25.0,  -h/2 + tf/2, -bf/6),  # 13
        (25.0,  -h/2 + tf/2, bf/6),   # 14
        (25.0,  h/2 - tf/2, 0.0),     # 15
        (25.0,  -h/2 + tf/2, 0.0),    # 16

    ]

    for dia, y, x in rebars:
        area = np.pi * dia**2 / 4.0          # mm²
        ops.fiber(x, y, area, steelTag)      # add to OpenSees model

    # -------------------- Plot (optional) -----------------------
    if plot:
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.set_xlabel('Width (mm)')
        ax.set_ylabel('Height (mm)')
        ax.set_title('I Section with Rebars (numbers shown)')
        ax.grid(True, ls='--', alpha=0.5)

        # Concrete geometry (light gray)
        geometry = [
            {"x": -bf/2, "y": h/2-tf, "w": bf, "h": tf},              # Top flange
            {"x": -bf/2, "y": -h/2+tf, "w": bf, "h": -tf},            # Bottom flange
            {"x": -tw/2, "y": -h/2+tf,    "w": tw, "h": h-2*tf}       # Web
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

        max_dim = max(bf, h) + 50
        ax.set_xlim(-max_dim/2, max_dim/2)
        ax.set_ylim(-max_dim/2, max_dim/2)
        ax.set_aspect('equal')
        plt.show()

    ELE_MASS = CONCRETE_DENSITY * AREA   # kg/mm
    
    return h, ELE_MASS
