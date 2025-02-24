import numpy as np
import openseespy.opensees as ops

#--------------------------------------------------    
# CONFINED CONCRETE FIBER SECTION WITH QUAD FIBERS
#--------------------------------------------------

def CONCRETE_CONFINED_REC_SECTION_QUAD(id, HSec, BSec, coverH, coverB, coreID, coverID, steelID, 
                            numBarsTop, barAreaTop, numBarsBot, barAreaBot, numBarsIntTot, 
                            barAreaInt, PLOT=True, CONCRETE_DENSITY=2500/1e9, STEEL_DENSITY=7850/1e9):
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
    DD = (HSec - 2 * coverH) / (numBarsIntTot-1)
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
    
    # Plotting the section
    if PLOT:
        import matplotlib.pyplot as plt
        from matplotlib.patches import Rectangle
        
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Core concrete (grey)
        core_width = 2 * coreZ
        core_height = 2 * coreY
        ax.add_patch(Rectangle((-coreZ, -coreY), core_width, core_height, 
                    facecolor='grey', edgecolor='black', label='Core Concrete'))
        
        # Cover concrete (lightgrey)
        # Top cover
        ax.add_patch(Rectangle((-coverZ, coreY), 2*coverZ, coverH, 
                    facecolor='lightgrey', edgecolor='black', label='Cover Concrete'))
        # Bottom cover
        ax.add_patch(Rectangle((-coverZ, -coreY-coverH), 2*coverZ, coverH, 
                    facecolor='lightgrey', edgecolor='black'))
        # Left cover
        ax.add_patch(Rectangle((-coverZ, -coreY), coverB, 2*coreY, 
                    facecolor='lightgrey', edgecolor='black'))
        # Right cover
        ax.add_patch(Rectangle((coreZ, -coreY), coverB, 2*coreY, 
                    facecolor='lightgrey', edgecolor='black'))

        # Reinforcement visualization --------------------------------------------
        # Top bars
        z_top = np.linspace(-coreZ, coreZ, numBarsTop)
        ax.scatter(z_top, [coreY]*numBarsTop, s=50, c='red', marker='o', label='Top Steel')

        # Bottom bars
        z_bot = np.linspace(-coreZ, coreZ, numBarsBot)
        ax.scatter(z_bot, [-coreY]*numBarsBot, s=50, c='red', marker='o', label='Bottom Steel')

        # Intermediate bars (sides)
        y_int = np.linspace(-coreY+DD, coreY-DD, numBarsIntTot//2)
        ax.scatter([coreZ]*len(y_int), y_int, s=30, c='blue', marker='s', label='Skin Steel')
        ax.scatter([-coreZ]*len(y_int), y_int, s=30, c='blue', marker='s')

        # Plot formatting
        ax.set_aspect('equal')
        ax.set_xlabel('Width (Z)')
        ax.set_ylabel('Height (Y)')
        ax.set_title(f'RC Rectangular Section (ID: {id})')
        ax.grid(True, alpha=0.3)
        ax.axhline(0, color='black', lw=0.5)
        ax.axvline(0, color='black', lw=0.5)
        ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1))
        plt.tight_layout()
        plt.show()

    return HSec, MASS

#----------------------------------------------------------------------------   

#-----------------------------------------------------    
# UNCONFINED CONCRETE FIBER SECTION WITH QUAD FIBERS
#-----------------------------------------------------

def CONCRETE_UNCONFINED_REC_SECTION_QUAD(id, HSec, BSec, coverH, coverB, coreID, coverID, steelID, 
                            numBarsTop, barAreaTop, numBarsBot, barAreaBot, numBarsIntTot, 
                            barAreaInt, PLOT=True, CONCRETE_DENSITY=2500/1e9, STEEL_DENSITY=7850/1e9):
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
    coreY = coverY - 0  # Distance from the section z-axis to the edge of the core concrete
    coreZ = coverZ - 0  # Distance from the section y-axis to the edge of the core concrete
    numBarsInt = numBarsIntTot // 2  # Number of intermediate bars per side
    DD = (HSec - 2 * coverH) / (numBarsIntTot-1)
    # Define the fiber section
    ops.section('Fiber', id)

    # Define the core patch
    ops.patch('quad', coreID, nfCoreZ, nfCoreY, -coreY, coreZ, -coreY, -coreZ, coreY, -coreZ, coreY, coreZ)

    # Define reinforcing layers
    ops.layer('straight', steelID, numBarsInt, barAreaInt, -coreY+coverH+DD, coreZ-coverB, coreY-coverH-DD, coreZ-coverB)  # Intermediate skin reinforcement (+z)
    ops.layer('straight', steelID, numBarsInt, barAreaInt, -coreY+coverH+DD, -coreZ+coverB, coreY-coverH-DD, -coreZ+coverB)  # Intermediate skin reinforcement (-z)
    ops.layer('straight', steelID, numBarsTop, barAreaTop, coreY-coverH, coreZ-coverB, coreY-coverH, -coreZ+coverB)  # Top layer reinforcement
    ops.layer('straight', steelID, numBarsBot, barAreaBot, -coreY+coverH, coreZ-coverB, -coreY+coverH, -coreZ+coverB)  # Bottom layer reinforcement
    
    # Plotting the section
    if PLOT:
        import matplotlib.pyplot as plt
        from matplotlib.patches import Rectangle
        
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Core concrete (grey)
        core_width = 2 * coreZ
        core_height = 2 * coreY
        ax.add_patch(Rectangle((-coreZ, -coreY), core_width, core_height, 
                    facecolor='lightgrey', edgecolor='black', label='Core Concrete'))
        
        # Reinforcement visualization --------------------------------------------
        # Top bars
        z_top = np.linspace(-coreZ+coverB, coreZ-coverB, numBarsTop)
        ax.scatter(z_top, [coreY-coverH]*numBarsTop, s=50, c='red', marker='o', label='Top Steel')

        # Bottom bars
        z_bot = np.linspace(-coreZ+coverB, coreZ-coverB, numBarsBot)
        ax.scatter(z_bot, [-coreY+coverH]*numBarsBot, s=50, c='red', marker='o', label='Bottom Steel')

        # Intermediate bars (sides)
        y_int = np.linspace(-coreY+coverH+DD, coreY-coverH-DD, numBarsInt)
        ax.scatter([coreZ-coverB]*len(y_int), y_int, s=30, c='blue', marker='s', label='Skin Steel')
        ax.scatter([-coreZ+coverB]*len(y_int), y_int, s=30, c='blue', marker='s')

        # Plot formatting
        ax.set_aspect('equal')
        ax.set_xlabel('Width (Z)')
        ax.set_ylabel('Height (Y)')
        ax.set_title(f'RC Rectangular Section (ID: {id})')
        ax.grid(True, alpha=0.3)
        ax.axhline(0, color='black', lw=0.5)
        ax.axvline(0, color='black', lw=0.5)
        ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1))
        plt.tight_layout()
        plt.show()

    return HSec, MASS
    
#----------------------------------------------------------------------------    
