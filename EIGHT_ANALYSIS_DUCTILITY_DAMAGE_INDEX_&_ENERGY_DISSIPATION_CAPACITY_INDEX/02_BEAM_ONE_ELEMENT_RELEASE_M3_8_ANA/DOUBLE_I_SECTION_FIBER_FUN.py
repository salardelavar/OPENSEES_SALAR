import numpy as np
import openseespy.opensees as ops
#----------------------------
# I SECTION
#----------------------------
def I_SECTION(secTag, matTag, PLOT=True, DENSITY=7850/1e9):
    #secTag = 1
    # Define geometric properties of the steel I section
    ops.section('Fiber', secTag)
    # Define section (FiberThermal)
    bf = 300               # [mm] Flange width
    tf = 20                # [mm] Flange thickness
    tw = 10                # [mm] Web thickness
    hw = 400               # [mm] Web height
    NUM = 100              # Number of fibers for web
    d = 2 * tf + hw
    # Mass per Length
    AREA = 2 * bf * tf + tw * hw
    MASS = DENSITY * AREA
    
    # Top flange fibers
    ops.fiber(bf / 2, hw / 2 + tf / 2, bf * tf / 4, matTag)
    ops.fiber(-bf / 2, hw / 2 + tf / 2, bf * tf / 4, matTag)
    # Bottom flange fibers
    ops.fiber(bf / 2, -hw / 2 - tf / 2, bf * tf / 4, matTag)
    ops.fiber(-bf / 2, -hw / 2 - tf / 2, bf * tf / 4, matTag)
    # Web fibers
    for i in range(NUM):
        yLoc = hw / 2 - i * (hw / NUM)
        ops.fiber(0.0, yLoc, tw * hw / NUM, matTag) 
    
    #-------------------
    # PLOT THE SECTION
    #-------------------
    if PLOT == True:
        import matplotlib.pyplot as plt
        import numpy as np
        fig, ax = plt.subplots(figsize=(8, 8))

        # Plot the top flange fibers
        flange_area = bf * tf / 4
        ax.add_patch(plt.Rectangle((-bf / 2, hw / 2), bf, tf, color='grey', alpha=0.7))

        # Plot the bottom flange fibers
        ax.add_patch(plt.Rectangle((-bf / 2, -hw / 2 - tf), bf, tf, color='grey', alpha=0.7))

        # Plot the web fibers
        fiber_height = hw / NUM
        for i in range(NUM):
            yLoc = hw / 2 - i * fiber_height
            ax.add_patch(plt.Rectangle((-tw / 2, yLoc - fiber_height / 2), tw, fiber_height, color='grey', alpha=0.7))

        # Add section outline
        #ax.add_patch(plt.Rectangle((-bf / 2, -hw / 2 - tf), bf, hw + 2 * tf, edgecolor='black', fill=False, linewidth=1.5))

        # Set plot limits
        ax.set_xlim([-bf / 2 - 20, bf / 2 + 20])
        ax.set_ylim([-hw / 2 - tf - 20, hw / 2 + tf + 20])
        ax.set_aspect('equal', adjustable='box')

        # Labels and title
        plt.xlabel('Width (bf)', fontsize=12)
        plt.ylabel('Height (hw + 2 * tf)', fontsize=12)
        plt.title('I-Section with Fibers', fontsize=14)

        # Show plot
        plt.grid(True)
        plt.show()
        
        
    
    return d, MASS # Return Section Height

#----------------------------------------------------------------------------

#------------------------------------------
# DOUBLE I SECTIONS WITH AND WITHOUT PLATE
#------------------------------------------

def DOUBLE_I_SECTION_FIBER(secTag, matTag, PLOT=True, DENSITY=7850/1e9):
    # Define double I-section with centered top/bottom plates
    # You can omit Plate by plate_t=0 and plate_l=0
    ops.section('Fiber', secTag)
    
    # Section dimensions
    bf = 300               # [mm] Flange width per I-section
    tf = 20                # [mm] Flange thickness
    tw = 10                # [mm] Web thickness
    hw = 400               # [mm] Web height
    spacing = 500          # [mm] Distance between I-section centers
    plate_t = 10           # [mm] Plate thickness
    plate_l = 500          # [mm] Plate length
    NUM = 100              # Web fibers per I-section
    
    # Total section depth
    d = 2 * (tf + plate_t) + hw
    # Mass per Length
    AREA = 2*(2 * bf * tf + tw * hw) + 2 *  (plate_t * plate_l)
    MASS = DENSITY * AREA  
    
    # Create two I-sections
    for y_shift in [-spacing/2, spacing/2]:
        # Flanges
        ops.fiber(y_shift - bf/2, hw/2 + tf/2, bf*tf/2, matTag)
        ops.fiber(y_shift + bf/2, hw/2 + tf/2, bf*tf/2, matTag)
        ops.fiber(y_shift - bf/2, -hw/2 - tf/2, bf*tf/2, matTag)
        ops.fiber(y_shift + bf/2, -hw/2 - tf/2, bf*tf/2, matTag)
        
        # Web fibers
        for i in range(NUM):
            yLoc = hw/2 - i*(hw/NUM)
            ops.fiber(y_shift, yLoc, tw*hw/NUM, matTag)
    
    # Add plates (centered between I-sections)
    if plate_t > 0 and plate_l > 0:
        # Top plate
        ops.fiber(0, hw/2 + tf + plate_t/2, plate_l*plate_t, matTag)
        # Bottom plate
        ops.fiber(0, -hw/2 - tf - plate_t/2, plate_l*plate_t, matTag)
    
    # Plotting
    if PLOT:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Draw I-sections
        for y_shift in [-spacing/2, spacing/2]:
            # Flanges
            ax.add_patch(plt.Rectangle((y_shift - bf/2, hw/2), bf, tf, 
                                     color='grey', alpha=0.7))
            ax.add_patch(plt.Rectangle((y_shift - bf/2, -hw/2 - tf), bf, tf, 
                                     color='grey', alpha=0.7))
            # Web
            ax.add_patch(plt.Rectangle((y_shift - tw/2, -hw/2), tw, hw, 
                                     color='grey', alpha=0.7))
        
        # Draw plates
        if plate_t > 0 and plate_l > 0:
            ax.add_patch(plt.Rectangle((-plate_l/2, hw/2 + tf), plate_l, plate_t,
                                     color='blue', alpha=0.5))
            ax.add_patch(plt.Rectangle((-plate_l/2, -hw/2 - tf - plate_t), plate_l, plate_t,
                                     color='blue', alpha=0.5))
        
        # Plot settings
        ax.set_xlim(-spacing - bf/2 - 20, spacing + bf/2 + 20)
        ax.set_ylim(-d/2 - 20, d/2 + 20)
        ax.set_aspect('equal')
        plt.xlabel('Width (mm)')
        plt.ylabel('Height (mm)')
        plt.title(f'Double I-Section with Plates ({plate_t}mm×{plate_l}mm)')
        plt.grid(True)
        plt.show()
    
    return d, MASS  # Return total section height and mass per unit length (omitting plate_t and plate_l)
