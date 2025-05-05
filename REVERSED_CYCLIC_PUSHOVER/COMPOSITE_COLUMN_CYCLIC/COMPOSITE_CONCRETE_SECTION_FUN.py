#----------------------------------------------------------------------    
# CONFINED CONCRETE FIBER RECTANGULAR WITH PLATE SECTION - QUAD FIBERS
#----------------------------------------------------------------------
import openseespy.opensees as ops
def CONCRETE_CONFINED_REC_PLATE_SECTION_QUAD(id, HSec, BSec, coverH, coverB, coreID, coverID, steelID, steelPlateID,
                                      numBarsTop, barAreaTop, numBarsBot, barAreaBot, numBarsIntTot, barAreaInt,
                                      plateThickness=10.0, PLOT=True, CONCRETE_DENSITY=2500/1e9, STEEL_DENSITY=7850/1e9):
    """
    Creates a rectangular reinforced concrete fiber section with surrounding steel plates in OpenSees and plots it using matplotlib.
    Uses plt.Rectangle() for concrete and steel plates, and plt.scatter() for reinforcement bars.

    Parameters:
    id (int): Section tag identifier
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
    """
    #import openseespy.opensees as ops
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.patches import Rectangle
    # Number of fibers
    nfCoreY = 10  # Number of fibers in core (y-direction)
    nfCoreZ = 5   # Number of fibers in core (z-direction)
    nfCoverY = 2  # Number of fibers in cover (y-direction)
    nfCoverZ = 2  # Number of fibers in cover (z-direction)
    
    # Calculate mass
    MASS = HSec * BSec * CONCRETE_DENSITY + (numBarsTop * barAreaTop + numBarsBot * barAreaBot + numBarsIntTot * barAreaInt) * STEEL_DENSITY 
    
    # Calculate dimensions
    coverY = HSec / 2.0  # Distance from the section z-axis to the edge of the cover concrete
    coverZ = BSec / 2.0  # Distance from the section y-axis to the edge of the cover concrete
    coreY = coverY - coverH  # Distance from the section z-axis to the edge of the core concrete
    coreZ = coverZ - coverB  # Distance from the section y-axis to the edge of the core concrete
    numBarsInt = numBarsIntTot // 2  # Number of intermediate bars per side
    DD = (HSec - 2 * coverH) / (numBarsIntTot - 1)  # Spacing between intermediate bars

    # Define the fiber section
    ops.section('Fiber', id)

    # Define the core patch
    ops.patch('quad', coreID, nfCoreZ, nfCoreY, -coreY, coreZ, -coreY, -coreZ, coreY, -coreZ, coreY, coreZ)

    # Define the four cover patches
    ops.patch('quad', coverID, 2, nfCoverY, -coverY, coverZ, -coreY, coreZ, coreY, coreZ, coverY, coverZ)  # Top cover
    ops.patch('quad', coverID, 2, nfCoverY, -coreY, -coreZ, -coverY, -coverZ, coverY, -coverZ, coreY, -coreZ)  # Bottom cover
    ops.patch('quad', coverID, nfCoverZ, 2, -coverY, coverZ, -coverY, -coverZ, -coreY, -coreZ, -coreY, coreZ)  # Left cover
    ops.patch('quad', coverID, nfCoverZ, 2, coreY, coreZ, coreY, -coreZ, coverY, -coverZ, coverY, coverZ)  # Right cover

    # Define reinforcing layers
    ops.layer('straight', steelID, numBarsInt, barAreaInt, -coreY+DD, coreZ, coreY-DD, coreZ)  # Intermediate skin reinforcement (+z)
    ops.layer('straight', steelID, numBarsInt, barAreaInt, -coreY+DD, -coreZ, coreY-DD, -coreZ)  # Intermediate skin reinforcement (-z)
    ops.layer('straight', steelID, numBarsTop, barAreaTop, coreY, coreZ, coreY, -coreZ)  # Top layer reinforcement
    ops.layer('straight', steelID, numBarsBot, barAreaBot, -coreY, coreZ, -coreY, -coreZ)  # Bottom layer reinforcement

    # Define steel plates
    plateY = coverY + plateThickness  # Outer edge of the steel plates in the y-direction
    plateZ = coverZ + plateThickness  # Outer edge of the steel plates in the z-direction
    ops.patch('quad', steelPlateID, 2, 2, -plateY, plateZ, -coverY, coverZ, coverY, coverZ, plateY, plateZ)  # Top plate
    ops.patch('quad', steelPlateID, 2, 2, -coverY, -coverZ, -plateY, -plateZ, plateY, -plateZ, coverY, -coverZ)  # Bottom plate
    ops.patch('quad', steelPlateID, 2, 2, -plateY, coverZ, -plateY, -coverZ, -coverY, -coverZ, -coverY, coverZ)  # Left plate
    ops.patch('quad', steelPlateID, 2, 2, coverY, coverZ, coverY, -coverZ, plateY, -coverZ, plateY, coverZ)  # Right plate

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

        # Reinforcement visualization --------------------------------------------
        # Top bars
        z_top = np.linspace(-coreZ, coreZ, numBarsTop)
        ax.scatter(z_top, [coreY] * numBarsTop, s=100, c='red', marker='o', label='Top Steel')

        # Bottom bars
        z_bot = np.linspace(-coreZ, coreZ, numBarsBot)
        ax.scatter(z_bot, [-coreY] * numBarsBot, s=100, c='purple', marker='o', label='Bottom Steel')

        # Intermediate bars (sides)
        y_int = np.linspace(-coreY + DD, coreY - DD, numBarsIntTot // 2)
        ax.scatter([coreZ] * len(y_int), y_int, s=60, c='blue', marker='o', label='Skin Steel')
        ax.scatter([-coreZ] * len(y_int), y_int, s=60, c='blue', marker='o')

        # Plot formatting
        ax.set_aspect('equal')
        ax.set_xlabel('Width (Z)')
        ax.set_ylabel('Height (Y)')
        ax.set_title(f'RC Rectangular Section with Steel Plates (ID: {id})')
        ax.grid(True, alpha=0.3)
        ax.axhline(0, color='black', lw=0.5)
        ax.axvline(0, color='black', lw=0.5)
        ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1))
        plt.tight_layout()
        plt.show()

    return HSec, MASS
