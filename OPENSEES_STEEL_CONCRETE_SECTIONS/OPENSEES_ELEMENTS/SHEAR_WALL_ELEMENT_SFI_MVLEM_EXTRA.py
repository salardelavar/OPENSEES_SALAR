def SHEAR_WALL_ELEMENT_SFI_MVLEM_EXTRA(ELE_TAG, WALL_DENSITY,
                                   NODEi, NODEj,
                                   Bsec, Hsec, plot=True):
    """
    Create a rectangular confined‑concrete SFI-MVLEM element (OpenSees).

    Parameters
    ----------
    secTag           : int   – identifiers (only secTag is used in this routine)
    Bsec, Hsec       : float – width and height of the rectangle (mm)
    WALL_DENSITY     : float – concrete density in kg/mm³  (≈ 2.5e‑9 kg/mm³)
    
    THIS PYTHON SCRIPT WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
    
    PAPER: Comparative Study of State-of-the-Art Macroscopic Models for Planar Reinforced Concrete Walls - ACI STRUCTURAL JOURNAL 
    INFO LINK: https://opensees.berkeley.edu/wiki/index.php/User_talk:Kkolozvari
    PEER REPORT : https://peer.berkeley.edu/sites/default/files/web_peer607_k._orakcal_l._massone_j._wallace_.pdf
    
    [1] Kolozvari K., Orakcal K., and Wallace J. W. (2015). “Shear-Flexure Interaction Modeling of reinforced Con
    crete Structural Walls and Columns under Reversed Cyclic Loading”, Pacific Earthquake Engineering Research
    Center, University of California, Berkeley, PEER Report No. 2015/12
    [2] Kolozvari K. (2013). “Analytical Modeling of Cyclic Shear-Flexure Interaction in Reinforced Concrete Struc
    tural Walls”, PhD Dissertation, University of California, Los Angeles.
    [3] Orakcal K., Massone L.M., and Ulugtekin D. (2012). “Constitutive Modeling of Reinforced Concrete Panel
    Behavior under Cyclic Loading”, Proceedings, 15th World Conference on Earthquake Engineering, Lisbon,
    Portugal.
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
    # Reinforcing steel in X
    fyX = 400            # [N/mm²] Steel Rebar Yield Strength   
    EsX = 2e5            # [N/mm²] Modulus of Elasticity
    eyX = fyX/EsX        # [mm/mm] Steel Rebar Yield Strain
    fuX = 1.1818*fyX     # [N/mm²] Steel Rebar Ultimate Strength
    esuX = 0.09          # [mm/mm] Steel Rebar Ultimate Strain
    EshX = (fuX - fyX)/(esuX - eyX)
    BsX = EshX / EsX
    R0X = 20.0         # initial value of curvature parameter
    a1X = 0.925        # curvature degradation parameter
    a2X = 0.15         # curvature degradation parameter
    # Reinforcing steel in Y
    fyY = 400            # [N/mm²] Steel Rebar Yield Strength   
    EsY = 2e5            # [N/mm²] Modulus of Elasticity
    eyY = fyY/EsY        # [mm/mm] Steel Rebar Yield Strain
    fuY = 1.1818*fyY     # [N/mm²] Steel Rebar Ultimate Strength
    esuY = 0.09          # [mm/mm] Steel Rebar Ultimate Strain
    EshY = (fuY - fyY)/(esuY - eyY)
    BsY = EshY / EsY
    R0Y = 20.0         # initial value of curvature parameter
    a1Y = 0.925        # curvature degradation parameter
    a2Y = 0.15         # curvature degradation parameter
    
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
    T10 = NN      # Thickness of depth 10
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
    
    m = 10
    print(m)
    c  = 0.4     # location of center of rotation (0.4 recommended)

    # Thickness of each macro‑fiber
    thicks  = [T1] + [T2] + [T3] + [T4] + [T5] + [T6] + [T7] + [T8] + [T9] + [T10]
    #print(thick)
    # Width of each macro‑fiber
    widths = [B1] + [B2] + [B3] + [B4] + [B5] + [B6] + [B7] + [B8] + [B9] + [B10]

    rho = WALL_DENSITY # Material density 
    rouX = 0.01        # Reinforcing ratio in horizontal (x) direction
    rouY = 0.0025      # Reinforcing ratio in vertical (x) direction
    # Shear resisting mechanism parameters
    nu = 1.10          # Concrete friction coefficient (0.0 < ν < 1.5)
    alfadow = 0.02     # Stiffness coefficient of reinforcement dowel action (0.0 < alfadow < 0.05)
    
    # Concrete material
    matTags = []
    for i in range(m):
        tagC = int(ELE_TAG*1000 + i + NODEj * m*1000)
        print(tagC)
        ops.uniaxialMaterial('ConcreteCM', tagC, fcC, ec0C, EcC, rcC, ecuC, ftC, -ec0C, rtC, etUC)  # confined concrete
        #ops.uniaxialMaterial('ConcreteCM', tag, fcU, ec0U, EcU, rcU, ecuU, ftU, -ec0U, rtU, etUU)  # unconfined concrete
        
        # Tag of uniaxialMaterial simulating horizontal (x) reinforcement
        tagRx = int(ELE_TAG*2000 + i + NODEj * m*2000)
        print(tagRx)
        ops.uniaxialMaterial('SteelMPF', tagRx, fyX, fyX, EsX, BsX, BsX, R0X, a1X, a2X)
        
        # Tag of uniaxialMaterial simulating horizontal (y) reinforcement
        tagRy = int(ELE_TAG*3000 + i + NODEj * m*3000)
        print(tagRy)
        ops.uniaxialMaterial('SteelMPF', tagRy, fyY, fyY, EsY, BsY, BsY, R0Y, a1Y, a2Y)
        tag = int(ELE_TAG*4000 + i + NODEj * m*4000)
        
        matTags.append(tag)
        # nDMaterial(’FSAM’, matTag, rho, sXTag, sYTag, concTag, rouX, rouY, nu, alfadow)
        ops.nDMaterial('FSAM', tag, WALL_DENSITY, tagRx, tagRy, tagC, rouX, rouY, nu, alfadow)  # Unconfined Concrete

    # MVLEM Element
    # element('SFI_MVLEM', eleTag, iNode, jNode, m, c, '-thick', *fiberThick, '-width', *fiberWidth, '-mat', *matTags)
    ops.element('SFI_MVLEM', ELE_TAG,
            NODEi, NODEj,
            m, c,
            '-thick',  *thicks,
            '-width',  *widths,
            '-mat',    *matTags)

    print(f"Element {ELE_TAG} is Defined - SFI-MVLEM defined successfully.")

    if plot:
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.set_xlabel('Width (mm)')
        ax.set_ylabel('Height (mm)')
        ax.set_title(f'Rectangular Shear Wall SFI-MVLEM Section with Rebars \nSection Rabar ratio in X: {100*rouX} (%) - Section Rabar ratio in Y: {100*rouY} (%)')
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
    
    AREA =  Bsec * Hsec
    
    MASS =  AREA * WALL_DENSITY