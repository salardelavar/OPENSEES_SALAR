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

def DOUBLE_I_SECTION(secTag, matTag, PLOT=True, DENSITY=7850/1e9):
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


#----------------------------------------------------------------------------

#-----------------------------------------------------
# DOUBLE I SECTIONS WITH AND WITHOUT PLATE WITH QUAD
#-----------------------------------------------------

def DOUBLE_I_SECTION_QUAD(secTag, matTag, PLOT=True, DENSITY=7850/1e9):
    """
    Define double I-section with centered top/bottom plates using quadrilateral fibers.
    
    Parameters:
    -----------
    secTag : int
        Section tag
    matTag : int
        Material tag
    PLOT : bool
        If True, plot the section geometry
    DENSITY : float
        Material density in kg/mm³ (default: 7850 kg/m³ converted to kg/mm³)
    
    Returns:
    --------
    tuple : (d, MASS)
        d : total section depth (mm)
        MASS : mass per unit length (kg/mm)
    """
    import openseespy.opensees as ops
    
    ops.section('Fiber', secTag)
    
    # Section dimensions
    bf = 300               # [mm] Flange width per I-section
    tf = 20                # [mm] Flange thickness
    tw = 10                # [mm] Web thickness
    hw = 400               # [mm] Web height
    spacing = 500          # [mm] Distance between I-section centers
    plate_t = 10           # [mm] Plate thickness
    plate_l = 500          # [mm] Plate length
    
    # Number of fibers for quad elements (should be even for symmetric meshing)
    nf_flange = 14          # Flange fibers in width direction
    nf_web = 80             # Web fibers in height direction
    nf_plate = 14           # Plate fibers in length direction
    
    # Total section depth
    d = 2 * (tf + plate_t) + hw
    
    # Calculate area for mass calculation
    # Area of one I-section: 2 flanges + web
    area_one_I = 2 * bf * tf + tw * hw
    # Total area of both I-sections + plates
    AREA = 2 * area_one_I + 2 * (plate_t * plate_l)
    MASS = DENSITY * AREA
    
    # Helper function to create quad fiber
    def add_quad_fiber(yCoords, zCoords, matTag):
        """Add quadrilateral fiber given y and z coordinates (4 points, counter-clockwise)"""
        ops.patch('quad', matTag, 1, 1,
                  *yCoords, *zCoords)
    
    # Create two I-sections
    for y_shift in [-spacing/2, spacing/2]:
        # Top flange (divided into nf_flange fibers along width)
        flange_fibers = nf_flange
        flange_y_inc = bf / flange_fibers
        for i in range(flange_fibers):
            y1 = y_shift - bf/2 + i * flange_y_inc
            y2 = y1 + flange_y_inc
            # Top flange fiber coordinates (counter-clockwise)
            yCoords = [y1, y2, y2, y1]
            zCoords = [hw/2 + tf, hw/2 + tf, hw/2, hw/2]
            add_quad_fiber(yCoords, zCoords, matTag)
        
        # Bottom flange
        for i in range(flange_fibers):
            y1 = y_shift - bf/2 + i * flange_y_inc
            y2 = y1 + flange_y_inc
            # Bottom flange fiber coordinates (counter-clockwise)
            yCoords = [y1, y2, y2, y1]
            zCoords = [-hw/2, -hw/2, -hw/2 - tf, -hw/2 - tf]
            add_quad_fiber(yCoords, zCoords, matTag)
        
        # Web (divided into nf_web fibers along height)
        web_fibers = nf_web
        web_z_inc = hw / web_fibers
        for i in range(web_fibers):
            z1 = hw/2 - i * web_z_inc
            z2 = z1 - web_z_inc
            # Web fiber coordinates (counter-clockwise)
            yCoords = [y_shift - tw/2, y_shift + tw/2, y_shift + tw/2, y_shift - tw/2]
            zCoords = [z1, z1, z2, z2]
            add_quad_fiber(yCoords, zCoords, matTag)
    
    # Add plates (centered between I-sections)
    if plate_t > 0 and plate_l > 0:
        # Top plate (divided into nf_plate fibers along length)
        plate_fibers = nf_plate
        plate_y_inc = plate_l / plate_fibers
        for i in range(plate_fibers):
            y1 = -plate_l/2 + i * plate_y_inc
            y2 = y1 + plate_y_inc
            # Top plate fiber coordinates (counter-clockwise)
            yCoords = [y1, y2, y2, y1]
            zCoords = [hw/2 + tf + plate_t, hw/2 + tf + plate_t, 
                      hw/2 + tf, hw/2 + tf]
            add_quad_fiber(yCoords, zCoords, matTag)
        
        # Bottom plate
        for i in range(plate_fibers):
            y1 = -plate_l/2 + i * plate_y_inc
            y2 = y1 + plate_y_inc
            # Bottom plate fiber coordinates (counter-clockwise)
            yCoords = [y1, y2, y2, y1]
            zCoords = [-hw/2 - tf, -hw/2 - tf,
                      -hw/2 - tf - plate_t, -hw/2 - tf - plate_t]
            add_quad_fiber(yCoords, zCoords, matTag)
    
    # Plotting
    if PLOT:
        import matplotlib.pyplot as plt
        import matplotlib.patches as patches
        
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Draw I-sections
        for y_shift in [-spacing/2, spacing/2]:
            # Top flange
            ax.add_patch(patches.Rectangle(
                (y_shift - bf/2, hw/2), bf, tf,
                facecolor='lightgrey', edgecolor='black', alpha=0.7
            ))
            # Bottom flange
            ax.add_patch(patches.Rectangle(
                (y_shift - bf/2, -hw/2 - tf), bf, tf,
                facecolor='lightgrey', edgecolor='black', alpha=0.7
            ))
            # Web
            ax.add_patch(patches.Rectangle(
                (y_shift - tw/2, -hw/2), tw, hw,
                facecolor='silver', edgecolor='black', alpha=0.7
            ))
        
        # Draw plates
        if plate_t > 0 and plate_l > 0:
            # Top plate
            ax.add_patch(patches.Rectangle(
                (-plate_l/2, hw/2 + tf), plate_l, plate_t,
                facecolor='lightblue', edgecolor='blue', alpha=0.5
            ))
            # Bottom plate
            ax.add_patch(patches.Rectangle(
                (-plate_l/2, -hw/2 - tf - plate_t), plate_l, plate_t,
                facecolor='lightblue', edgecolor='blue', alpha=0.5
            ))
        
        # Add fiber mesh visualization (optional)
        # Plot mesh lines for better visualization
        for y_shift in [-spacing/2, spacing/2]:
            # Flange mesh lines
            for i in range(nf_flange + 1):
                x = y_shift - bf/2 + i * (bf/nf_flange)
                # Top flange
                ax.plot([x, x], [hw/2, hw/2 + tf], 'k-', lw=0.5, alpha=0.3)
                # Bottom flange
                ax.plot([x, x], [-hw/2 - tf, -hw/2], 'k-', lw=0.5, alpha=0.3)
            
            # Web mesh lines
            for i in range(nf_web + 1):
                z = hw/2 - i * (hw/nf_web)
                ax.plot([y_shift - tw/2, y_shift + tw/2], [z, z], 
                       'k-', lw=0.5, alpha=0.3)
        
        # Plate mesh lines
        if plate_t > 0 and plate_l > 0:
            for i in range(nf_plate + 1):
                x = -plate_l/2 + i * (plate_l/nf_plate)
                # Top plate
                ax.plot([x, x], [hw/2 + tf, hw/2 + tf + plate_t], 
                       'b-', lw=0.5, alpha=0.3)
                # Bottom plate
                ax.plot([x, x], [-hw/2 - tf - plate_t, -hw/2 - tf], 
                       'b-', lw=0.5, alpha=0.3)
        
        # Plot settings
        ax.set_xlim(-spacing - bf/2 - 20, spacing + bf/2 + 20)
        ax.set_ylim(-d/2 - 20, d/2 + 20)
        ax.set_aspect('equal')
        ax.set_xlabel('Width (mm)')
        ax.set_ylabel('Height (mm)')
        ax.set_title(f'Double I-Section with Plates ({plate_t}mm×{plate_l}mm)\n'
                    f'Using {nf_flange}x{nf_web} quad fibers per I-section')
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()
    
    return d, MASS