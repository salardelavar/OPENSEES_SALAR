def SHEAR_WALL_ELEMENT_MVLEM_EXTRA(ELE_TAG, WALL_DENSITY,
                                   NODEi, NODEj,
                                   Bsec, Hsec, plot=True):
    """
    Create a rectangular confined‑concrete MVLEM element (OpenSees).

    Parameters
    ----------
    secTag           : int   – identifiers (only secTag is used in this routine)
    Bsec, Hsec       : float – width and height of the rectangle (mm)
    WALL_DENSITY     : float – concrete density in kg/mm³  (≈ 2.5e‑9 kg/mm³)
    
    THIS PYTHON SCRIPT WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
    
    PAPER: Comparative Study of State-of-the-Art Macroscopic Models for Planar Reinforced Concrete Walls - ACI STRUCTURAL JOURNAL 
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
    fcU = -25.0                 # [N/mm²] Concrete Compressive Strength
    EcU = 4700 * np.sqrt(-fcU)  # [N/mm^2] Concrete Elastic Modulus
    ec0U = 2*fcU/EcU            # [mm/mm] Concrete Compressive Strain
    fcUU = 0.2*fcU              # [N/mm²] Concrete Compressive Ultimate Strength
    ecuU = 5*ec0U               # [mm/mm] Concrete Compressive Ultimate Strain
    rtU = 1.2;			        # shape parameter - tension
    rcU = 7.0;			        # shape parameter - compression
    # Core concrete (confined)
    Kfc = 1.3;			       # ratio of confined to unconfined concrete strength
    fcC = Kfc*fcU              # [N/mm²] Concrete Compressive Strength
    EcC = 4700 * np.sqrt(-fcC) # [N/mm^2] Concrete Elastic Modulus
    ec0C = 2*fcC/EcC           # [mm/mm] Concrete Compressive Strain
    fcUC = 0.65*fcC            # [N/mm²] Concrete Compressive Ultimate Strength
    ecuC = 15*ec0C             # [mm/mm] Concrete Compressive Ultimate Strain
    Lambda = 0.1;	           # ratio between unloading slope
    rtC = 1.2;			       # shape parameter - tension
    rcC = 7.42;			       # shape parameter - compression    
    # tensile-strength properties
    ftC = 0.7 * np.sqrt(-fcC)  # [N/mm²] tensile strength +tension
    ftU = 0.7 * np.sqrt(-fcU)  # [N/mm²] tensile strength +tension
    EtsC = ftC/np.abs(ec0C)    # [N/mm²] tension softening stiffness
    EtsU = ftU/np.abs(ec0U)	   # [N/mm²] tension softening stiffness
    etUC = -2*ec0C             # [mm/mm] cracking strain - tension	
    etUU = -2*ec0U             # [mm/mm] cracking strain - tension	
    
    # STEEL
    # Reinforcing steel
    fy = 400          # [N/mm²] Steel Rebar Yield Strength   
    Es = 2e5          # [N/mm²] Modulus of Elasticity
    ey = fy/Es        # [mm/mm] Steel Rebar Yield Strain
    fu = 1.1818*fy    # [N/mm²] Steel Rebar Ultimate Strength
    esu = 0.09        # [mm/mm] Steel Rebar Ultimate Strain
    Esh = (fu - fy)/(esu - ey)
    Bs = Esh / Es
    R0 = 20.0         # initial value of curvature parameter
    a1 = 0.925        # curvature degradation parameter
    a2 = 0.15         # curvature degradation parameter
    
    NN = Hsec / 10
    T1 = NN       # Thickness of depth 01
    T2 = NN       # Thickness of depth 02
    T3 = NN       # Thickness of depth 03
    T4 = NN       # Thickness of depth 04
    T5 = NN       # Thickness of depth 05
    T6 = NN       # Thickness of depth 06
    T7 = NN       # Thickness of depth 07
    T8 = NN       # Thickness of depth 08
    T9 = NN       # Thickness of depth 09
    T10 = NN       # Thickness of depth 10
    B1 = 500.0    # Width of depth 01
    B2 = 500.0    # Width of depth 02
    B3 = 500.0    # Width of depth 03
    B4 = 500.0    # Width of depth 04
    B5 = 500.0    # Width of depth 05
    B6 = 500.0    # Width of depth 06
    B7 = 500.0    # Width of depth 07
    B8 = 500.0    # Width of depth 08
    B9 = 500.0    # Width of depth 09
    B10 = 500.0   # Width of depth 10
    R1 = 0.01     # Reinforcing ratios of depth 01
    R2 = 0.01     # Reinforcing ratios of depth 02
    R3 = 0.01     # Reinforcing ratios of depth 03
    R4 = 0.01     # Reinforcing ratios of depth 04
    R5 = 0.01     # Reinforcing ratios of depth 05
    R6 = 0.01     # Reinforcing ratios of depth 06
    R7 = 0.01     # Reinforcing ratios of depth 07
    R8 = 0.01     # Reinforcing ratios of depth 08
    R9 = 0.01     # Reinforcing ratios of depth 09
    R10 = 0.01    # Reinforcing ratios of depth 10
    
    m = 10
    #print(m)
    c  = 0.4     # location of center of rotation (0.4 recommended)

    # Thickness of each macro‑fiber
    thicks  = [T1] + [T2] + [T3] + [T4] + [T5] + [T6] + [T7] + [T8] + [T9] + [T10]
    #print(thick)
    # Width of each macro‑fiber
    widths = [B1] + [B2] + [B3] + [B4] + [B5] + [B6] + [B7] + [B8] + [B9] + [B10]

    # Density (not required, set to 0)
    rho = [R1] + [R2] + [R3] + [R4] + [R5] + [R6] + [R7] + [R8] + [R9] + [R10]

    # Concrete material
    matConcreteTags = []
    for i in range(m):
        tag = int(ELE_TAG*1000 + i + NODEj * m*100)
        #print(tag)
        matConcreteTags.append(tag)
        ops.uniaxialMaterial('ConcreteCM', tag, fcC, ec0C, EcC, rcC, ecuC, ftC, -ec0C, rtC, etUC)  # confined concrete
        #ops.uniaxialMaterial('ConcreteCM', tag, fcU, ec0U, EcU, rcU, ecuU, ftU, -ec0U, rtU, etUU, '-GapClose', 1)  # unconfined concrete


    # Steel material
    matSteelTags = []
    for i in range(m):
        tag = int(ELE_TAG*2000 + i + NODEj * m*200)
        #print(tag)
        matSteelTags.append(tag) 
        ops.uniaxialMaterial('SteelMPF', tag, fy, fy, Es, Bs, Bs, R0, a1, a2)

    # Shear material
    matShearTag = ELE_TAG*3
    # 1. Area = Bsec * Hsec  
    # Note: In OpenSees, units must be consistent. 
    # If using N and mm: E is in MPa (N/mm^2)
    AREA = Bsec * Hsec
    # 2. Iz = (t * B^3) / 12  (B is the depth of the section in bending)
    Iz = (Bsec * (Hsec**3)) / 12.0
    # 3. Avy = 5/6 * Area (standard for rectangular section)
    Avy = (5.0 / 6.0) * AREA
    # 4. E_mod (Modulus of Elasticity)
    # Formula: E = 4700 * sqrt(fc') in MPa
    E_mod = 4700 * np.sqrt(-fcC)

    # 5. G_mod (Shear Modulus)
    # Formula: G = E / (2 * (1 + nu))
    # Poisson's ratio for concrete (nu) is typically 0.2
    nu = 0.2
    G_mod = E_mod / (2.0 * (1.0 + nu)) 
    E_SHEAR = AREA * G_mod # Shear Stiffness
    # uniaxialMaterial Elastic $matTag $E <$eta> <$Eneg>
    # NOTE: large shear stiffness assigned since only flexural response
    ops.uniaxialMaterial('Elastic', matShearTag, E_SHEAR)

    # MVLEM Element
    ops.element('MVLEM', ELE_TAG, WALL_DENSITY,
            NODEi, NODEj,
            m, c,
            '-thick',  *thicks,
            '-width',  *widths,
            '-rho',    *rho,
            '-matConcrete', *matConcreteTags,
            '-matSteel',    *matSteelTags,
            '-matShear',    matShearTag)

    print(f"Element {ELE_TAG} is Defined - MVLEM defined successfully.")
    SUM_RO = (R1 + R2 + R3 + R4 + R5 + R6 + R7 + R8 + R9+ R10) * 100 
    if plot:
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.set_xlabel('Width (mm)')
        ax.set_ylabel('Height (mm)')
        ax.set_title(f'Rectangular Shear Wall MVLEM Section with Rebars \nSection Rabar Ratio: {SUM_RO} (%)')
        ax.grid(True, ls='--', alpha=0.5)

        # Concrete block (light gray)
        rect = patches.Rectangle((-Bsec/2, -Hsec/2),
                                 Bsec, Hsec,
                                 linewidth=1.5,
                                 edgecolor='black',
                                 facecolor='lightgray')
        ax.add_patch(rect)
        
        max_dim = max(Bsec, Hsec) + 50
        ax.set_xlim(-max_dim/2, max_dim/2)
        ax.set_ylim(-max_dim/2, max_dim/2)
        ax.set_aspect('equal')
        plt.show()
    
    MASS =  AREA * WALL_DENSITY