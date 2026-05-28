#----------------------------------------------------------------------    
# CONFINED CONCRETE FIBER RECTANGULAR WITH PLATE SECTION - QUAD FIBERS
#----------------------------------------------------------------------
def CONCRETE_CONFINED_REC_PLATE_SECTION_QUAD_EXTRA(secTag, MAT_TYPE, HSec, BSec, coverH, coverB,
                                      plateThickness=10.0, PLOT=True, CONCRETE_DENSITY=2500/1e9, STEEL_DENSITY=7850/1e9):
    """
    Creates a rectangular reinforced concrete fiber section with surrounding steel plates in OpenSees and plots it using matplotlib.
    Uses plt.Rectangle() for concrete and steel plates, and plt.scatter() for reinforcement bars.

    Parameters:
    secTag (int): Section tag identifier
    HSec (float): Height of the section
    BSec (float): Width of the section
    coverH (float): Cover thickness in the height direction
    coverB (float): Cover thickness in the width direction
    coreID (int): Material tag for confined core concrete
    coverID (int): Material tag for unconfined cover concrete
    steelID (int): Material tag for reinforcing steel
    steelPlateID (int): Material tag for the steel plates
    numBarsTop (int): Number of top reinforcing bars
    barAreaTop (float): Cross-sectional area of each top reinforcing bar
    numBarsBot (int): Number of bottom reinforcing bars
    barAreaBot (float): Cross-sectional area of each bottom reinforcing bar
    numBarsIntTot (int): Total number of intermediate reinforcing bars
    barAreaInt (float): Cross-sectional area of each intermediate reinforcing bar
    plateThickness (float): Thickness of the steel plates (default=10.0)
    PLOT (bool): Whether to plot the section (default=True)
    CONCRETE_DENSITY (float): Density of concrete (default=2500 kg/m³ converted to consistent units)
    STEEL_DENSITY (float): Density of steel (default=7850 kg/m³ converted to consistent units)
    
    THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
    """
    #import openseespy.opensees as ops
    import matplotlib.pyplot as plt
    import numpy as np
    import openseespy.opensees as ops
    from matplotlib.patches import Rectangle
    import matplotlib.patches as patches
    
    # Define materials for nonlinear elements
    # Define parameters (units: mm, N)
    # ------------------------------------------
    # CONCRETE                  tag   f'c        ec0   f'cu        ecu
    # Cover concrete (unconfined)
    fcU = -18                 # [N/mm²] Concrete Compressive Strength
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
    ftC = 0.7 * np.sqrt(-fcC)  # tensile strength +tension
    ftU = 0.7 * np.sqrt(-fcU)  # tensile strength +tension
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

    # Steel Plate
    fyP = 240          # [N/mm²] Steel I Section Yield Strength   
    EsP = 2e5          # [N/mm²] I Section Modulus of Elasticity
    eyP = fyP/EsP      # [mm/mm] Steel I Section Yield Strain
    fuP = 1.1818*fyP   # [N/mm²] Steel I Section Ultimate Strength
    esuP = 0.25        # [mm/mm] Steel I Section Ultimate Strain
    EshP = (fuP - fyP)/(esuP - eyP)
    BsP = EshP / EsP
    
    coreTag, coverTag, steelTag, steelPlateTag = secTag + 100, secTag + 200, secTag + 300, secTag + 400
    numBarsTop, barAreaTop = 5, np.pi *(18**2)/4
    numBarsBot, barAreaBot = 5, np.pi *(20**2)/4
    numBarsIntTot, barAreaInt = 4, np.pi *(5**2)/4
    
    if MAT_TYPE == 'ELASTIC':
        ops.uniaxialMaterial('Steel01', steelTag, fy, Es, Bs)         # REBAR
        ops.uniaxialMaterial('Steel01', steelPlateTag, fyP, EsP, BsP) # PLATE 
    if MAT_TYPE == 'INELASTIC':  
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
                             damage1, damage2, beta) # REBAR
        
        ops.uniaxialMaterial('Hysteretic', steelPlateTag,
                             fyP, eyP,
                             fuP, esuP,
                             0.2*fuP, 1.1*esuP,
                             -fyP, -eyP,
                             -fuP, -esuP,
                             -0.2*fuP, -1.1*esuP,
                             pinchX, pinchY,
                             damage1, damage2, beta) # PLATE
        # INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Hysteretic_Material


    #ops.uniaxialMaterial('Concrete01', coreTag, fcC, ec0C, fcUC, ecuC)  # Core concrete (confined)
    #ops.uniaxialMaterial('Concrete01', coverTag, fcU, ec0U, fcUU, ecuU) # Cover concrete (unconfined)
    ops.uniaxialMaterial('Concrete02', coreTag, fcC, ec0C, fcUC, ecuC, LambdaC, ftC, EtsC) # build core concrete (confined)
    ops.uniaxialMaterial('Concrete02', coverTag, fcU, ec0U, fcUU, ecuU, LambdaU, ftU, EtsU) # build cover concrete (unconfined)
    
    # Number of fibers
    nfCoreY = 10  # Number of fibers in core (y-direction)
    nfCoreZ = 5   # Number of fibers in core (z-direction)
    nfCoverY = 7  # Number of fibers in cover (y-direction)
    nfCoverZ = 7  # Number of fibers in cover (z-direction)
    
    # Calculate mass
    MASS = HSec * BSec * CONCRETE_DENSITY + ((plateThickness*(HSec+2*plateThickness)) + (plateThickness*(BSec-2*plateThickness))) * STEEL_DENSITY
    
    # Calculate dimensions
    coverY = HSec / 2.0  # Distance from the section z-axis to the edge of the cover concrete
    coverZ = BSec / 2.0  # Distance from the section y-axis to the edge of the cover concrete
    coreY = coverY - coverH  # Distance from the section z-axis to the edge of the core concrete
    coreZ = coverZ - coverB  # Distance from the section y-axis to the edge of the core concrete
    numBarsInt = numBarsIntTot // 2  # Number of intermediate bars per side
    DD = (HSec - 2 * coverH) / (numBarsIntTot - 1)  # Spacing between intermediate bars

    # Define the fiber section
    ops.section('Fiber', secTag)

    # Define the core patch
    ops.patch('quad', coreTag, nfCoreZ, nfCoreY, -coreY, coreZ, -coreY, -coreZ, coreY, -coreZ, coreY, coreZ)

    # Define the four cover patches
    ops.patch('quad', coverTag, 2, nfCoverY, -coverY, coverZ, -coreY, coreZ, coreY, coreZ, coverY, coverZ)  # Top cover
    ops.patch('quad', coverTag, 2, nfCoverY, -coreY, -coreZ, -coverY, -coverZ, coverY, -coverZ, coreY, -coreZ)  # Bottom cover
    ops.patch('quad', coverTag, nfCoverZ, 2, -coverY, coverZ, -coverY, -coverZ, -coreY, -coreZ, -coreY, coreZ)  # Left cover
    ops.patch('quad', coverTag, nfCoverZ, 2, coreY, coreZ, coreY, -coreZ, coverY, -coverZ, coverY, coverZ)  # Right cover

    # Define reinforcing layers
    ops.layer('straight', steelTag, numBarsInt, barAreaInt, -coreY+DD, coreZ, coreY-DD, coreZ)  # Intermediate skin reinforcement (+z)
    ops.layer('straight', steelTag, numBarsInt, barAreaInt, -coreY+DD, -coreZ, coreY-DD, -coreZ)  # Intermediate skin reinforcement (-z)
    ops.layer('straight', steelTag, numBarsTop, barAreaTop, coreY, coreZ, coreY, -coreZ)  # Top layer reinforcement
    ops.layer('straight', steelTag, numBarsBot, barAreaBot, -coreY, coreZ, -coreY, -coreZ)  # Bottom layer reinforcement

    # Define steel plates
    plateY = coverY + plateThickness  # Outer edge of the steel plates in the y-direction
    plateZ = coverZ + plateThickness  # Outer edge of the steel plates in the z-direction
    ops.patch('quad', steelPlateTag, 2, 2, -plateY, plateZ, -coverY, coverZ, coverY, coverZ, plateY, plateZ)  # Top plate
    ops.patch('quad', steelPlateTag, 2, 2, -coverY, -coverZ, -plateY, -plateZ, plateY, -plateZ, coverY, -coverZ)  # Bottom plate
    ops.patch('quad', steelPlateTag, 10, 2, -plateY, coverZ, -plateY, -coverZ, -coverY, -coverZ, -coverY, coverZ)  # Left plate
    ops.patch('quad', steelPlateTag, 10, 2, coverY, coverZ, coverY, -coverZ, plateY, -coverZ, plateY, coverZ)  # Right plate
    
    REBAR = 25.0 # [mm] Rebar Diameter
    # -------------------- Rebars -------------------------------
    # (diameter [mm], y‑coord [mm], x‑coord [mm])
    rebars = [
        (REBAR,  HSec/2 - coverB, -BSec/2 + coverB),    # 1
        (REBAR,  HSec/2 - coverB, BSec/2 - coverB),     # 2
        (REBAR,  HSec/2 - coverB, 0.0),                 # 3
        (REBAR,  -HSec/2 + coverB, -BSec/2 + coverB),   # 4
        (REBAR,  -HSec/2 + coverB, +BSec/2 - coverB),   # 5
        (REBAR,  -HSec/2 + coverB, 0.0),                # 6
        (16.0, HSec/4 - coverB, -BSec/2 + coverB),      # 7
        (16.0, HSec/4 - coverB, +BSec/2 - coverB),      # 8
        (16.0,  -HSec/4 + coverB, -BSec/2 + coverB),    # 9
        (16.0,  -HSec/4 + coverB, +BSec/2 - coverB),    # 10
        (14.0,  0.0, -BSec/2 + coverB),                 # 11
        (14.0,  0.0, +BSec/2 - coverB),                 # 12
    ]

    for dia, y, x in rebars:
        area = np.pi * dia**2 / 4.0          # mm²
        ops.fiber(x, y, area, steelTag)      # add to OpenSees model
    # Plotting the section
    if PLOT:
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Core concrete (grey)
        core_width = 2 * coreZ
        core_height = 2 * coreY
        ax.add_patch(Rectangle((-coreZ, -coreY), core_width, core_height, 
                    facecolor='grey', edgecolor='black', label='Core Concrete'))
        
        # Cover concrete (lightgrey)
        # Top cover
        ax.add_patch(Rectangle((-coverZ, coreY), 2 * coverZ, coverH, 
                    facecolor='lightgrey', edgecolor='black', label='Cover Concrete'))
        # Bottom cover
        ax.add_patch(Rectangle((-coverZ, -coreY - coverH), 2 * coverZ, coverH, 
                    facecolor='lightgrey', edgecolor='black'))
        # Left cover
        ax.add_patch(Rectangle((-coverZ, -coreY), coverB, 2 * coreY, 
                    facecolor='lightgrey', edgecolor='black'))
        # Right cover
        ax.add_patch(Rectangle((coreZ, -coreY), coverB, 2 * coreY, 
                    facecolor='lightgrey', edgecolor='black'))

        # Steel plates (blue)
        # Top plate
        ax.add_patch(Rectangle((-plateZ, coverY), 2 * plateZ, plateThickness, 
                    facecolor='blue', edgecolor='black', alpha=0.5, label='Steel Plate'))
        # Bottom plate
        ax.add_patch(Rectangle((-plateZ, -coverY - plateThickness), 2 * plateZ, plateThickness, 
                    facecolor='blue', edgecolor='black', alpha=0.5))
        # Left plate
        ax.add_patch(Rectangle((-plateZ, -coverY), plateThickness, 2 * coverY, 
                    facecolor='blue', edgecolor='black', alpha=0.5))
        # Right plate
        ax.add_patch(Rectangle((coverZ, -coverY), plateThickness, 2 * coverY, 
                    facecolor='blue', edgecolor='black', alpha=0.5))


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
        ax.set_title(f'RC Rectangular Section with Steel Plates (secTag: {secTag})')
 
        max_dim = max(BSec, HSec) + 50
        ax.set_xlim(-max_dim/2, max_dim/2)
        ax.set_ylim(-max_dim/2, max_dim/2)
        ax.grid(True, ls=':', alpha=0.5)
        ax.legend()
        plt.show()

    return HSec, MASS
