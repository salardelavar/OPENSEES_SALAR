def CONCRETE_CIRCULALR_HOLLOW_PIPE_SECTION_FUN(secTag, STEEL_TYPE, Do, t,
                                     numSubdivCirc, numSubdivRad,
                                     fc, Kfc,
                                     CONCRETE_DENSITY, plot=True):
    """
    Create a Hollow Pipe confined‑concrete fiber section (OpenSees) and,
    if PLOT=True, draw a simple 2‑D sketch of the section.
    (Circular Hollow Section)

    Parameters
    ----------
    secTag          : int   – section identifier
    Do              : float – outer width (along y and z) (mm)
    t               : float – wall thickness (mm)
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
    fcU = -fc  # [N/mm²] Concrete Compressive Strength
    Ec = 4700 * np.sqrt(-fcU)  # [N/mm^2] Concrete Elastic Modulus
    ec0U = 2 * fcU / Ec        # [mm/mm] Concrete Compressive Strain
    fcUU = 0.2 * fcU           # [N/mm²] Concrete Compressive Ultimate Strength
    ecuU = 5 * ec0U            # [mm/mm] Concrete Compressive Ultimate Strain

    # Core concrete (confined)
    # Kfc = 1.3;			      # ratio of confined to unconfined concrete strength
    fcC = Kfc * fcU            # [N/mm²] Concrete Compressive Strength
    Ec = 4700 * np.sqrt(-fcC)  # [N/mm^2] Concrete Elastic Modulus
    ec0C = 2 * fcC / Ec        # [mm/mm] Concrete Compressive Strain
    fcUC = 0.65 * fcC          # [N/mm²] Concrete Compressive Ultimate Strength
    ecuC = 15 * ec0C           # [mm/mm] Concrete Compressive Ultimate Strain
    Lambda = 0.1;              # ratio between unloading slope
    # tensile-strength properties
    ftC = 0.7 * np.sqrt(-fcC)  # [N/mm²] tensile strength +tension
    ftU = 0.7 * np.sqrt(-fcU)  # [N/mm²] tensile strength +tension
    EtsC = ftC / 0.002;        # [N/mm²] tension softening stiffness
    EtsU = ftU / 0.002;        # [N/mm²] tension softening stiffness

    # STEEL
    # Reinforcing steel
    fy = 400                   # [N/mm²] Steel Rebar Yield Strength
    Es = 2e5                   # [N/mm²] Modulus of Elasticity
    ey = fy / Es               # [mm/mm] Steel Rebar Yield Strain
    fu = 1.1818 * fy           # [N/mm²] Steel Rebar Ultimate Strength
    esu = 0.09                 # [mm/mm] Steel Rebar Ultimate Strain
    Esh = (fu - fy) / (esu - ey)
    Bs = Esh / Es

    coreTag, coverTag, steelTag = secTag + 100, secTag + 200, secTag + 300
    if STEEL_TYPE == 'ELASTIC':
        #ops.uniaxialMaterial('Steel01', steelTag, fy, Es, Bs)
        ops.uniaxialMaterial('Elastic', steelTag, Es)
    if STEEL_TYPE == 'INELASTIC':
        pinchX = 0.8  # Pinching factor in X direction
        pinchY = 0.5  # Pinching factor in Y direction
        damage1 = 0.0  # Damage due to ductility
        damage2 = 0.0  # Damage due to energy
        beta = 0.1  # Stiffness degradation parameter
        ops.uniaxialMaterial('Hysteretic', steelTag, fy, ey, fu, esu, 0.2 * fu, 1.1 * esu,
                             -fy, -ey, -fu, -esu, -0.2 * fu, -1.1 * esu, pinchX, pinchY,
                             damage1, damage2, beta)
        # INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Hysteretic_Material

    ops.uniaxialMaterial('Concrete01', coreTag, fcC, ec0C, fcUC, ecuC)  # Core concrete (confined)
    ops.uniaxialMaterial('Concrete01', coverTag, fcU, ec0U, fcUU, ecuU)  # Cover concrete (unconfined)
    # ops.uniaxialMaterial('Concrete02', coreTag, fcC, ec0C, fcUC, ecuC, Lambda, ftC, EtsC) # build core concrete (confined)
    # ops.uniaxialMaterial('Concrete02', coverTag, fcU, ec0U, fcUU, ecuU, Lambda, ftU, EtsU) # build cover concrete (unconfined)

    import openseespy.opensees as ops
    import numpy as np
    # ------------------------------
    # Pipe Geometry
    # ------------------------------
    # Do = 260.0           # Outer diameter (mm)
    # t = 10.0             # Wall thickness (mm)
    Di = Do - 2 * t        # Inner diameter (mm)
    r_out = Do / 2
    r_in = Di / 2

    h_o = Do / 2
    T_h = t
    b_o = Do / 2
    T_b = t
    h_i = Di / 2

    # ------------------------------
    # Fiber Section Definition
    # ------------------------------
    # numSubdivCirc = 32   # circumferential divisions
    # numSubdivRad = 4     # radial divisions

    ops.section('Fiber', secTag)

    # Fiber ring of steel material
    ops.patch('circ', coreTag, numSubdivRad, numSubdivCirc, 0.0, 0.0, r_in, r_out)

    # ------------------------------
    # Basic check (optional geometry)
    # ------------------------------
    AREA = np.pi * (r_out ** 2 - r_in ** 2)
    Iz = (np.pi / 4) * (r_out ** 4 - r_in ** 4)
    print(f"Section area A = {AREA:.6f} mm²")
    print(f"Moment of inertia Iz = {Iz:.8f} mm⁴")
    # -------------------- Rebars -------------------------------
    # (diameter [mm], y‑coord [mm], x‑coord [mm])
    rebars = [
        # (diameter[mm], y[mm], x[mm])
        (25.0, h_o - 0.50 * T_h, 0.0),   # 1
        (25.0, -h_o + 0.50 * T_h, 0.0),  # 2
        (25.0, 0.0, h_o - 0.50 * T_h),   # 3
        (25.0, 0.0, -h_o + 0.50 * T_h),  # 4
        (25.0, 0.7071*(h_o - 0.50 * T_h), 0.7071*(h_o - 0.50 * T_h)),   # 5
        (25.0, 0.7071*(-h_o + 0.50 * T_h), 0.7071*(-h_o + 0.50 * T_h)), # 6
        (25.0, 0.7071*(h_o - 0.50 * T_h), -0.7071*(h_o - 0.50 * T_h)),   # 7
        (25.0, 0.7071*(-h_o + 0.50 * T_h),-0.7071*(-h_o + 0.50 * T_h)), # 8
    ]

    # -------------------------
    # Plot section geometry
    # -------------------------
    if plot:
        import matplotlib.pyplot as plt
        import matplotlib.patches as patches
        import numpy as np
        
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.set_xlabel('Width (mm)')
        ax.set_ylabel('Height (mm)')
        ax.set_title('Concrete Pipe Hollow Section with Rebars (numbers shown)')
        ax.grid(True, ls='--', alpha=0.5)

        # Create angles for plotting circles
        theta = np.linspace(0, 2 * np.pi, 100)

        # Calculate x and y coordinates for the outer circle
        x_out = r_out * np.cos(theta)
        y_out = r_out * np.sin(theta)

        # Calculate x and y coordinates for the inner circle
        x_in = r_in * np.cos(theta)
        y_in = r_in * np.sin(theta)
        
        # Plotting Pipe Section
        ax.plot(x_out, y_out, color='black')
        ax.plot(x_in, y_in, color='black')

        # Fill the hollow part with white to make it look hollow
        ax.fill_between(x_out, y_out, color='lightgray', alpha=0.5)
        ax.fill_between(x_in, y_in, color='white')  # Fill the inner circle with white
        
        # Plotting Rebars
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

        max_dim = Do + 50
        ax.set_xlim(-max_dim/2, max_dim/2)
        ax.set_ylim(-max_dim/2, max_dim/2)
        ax.set_aspect('equal')
        plt.show()

    ELE_MASS = CONCRETE_DENSITY * AREA   # kg/mm
    
    return Do, ELE_MASS
