import openseespy.opensees as ops

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
#----------------------------------------    
# STEEL FIBER SECTION WITH QUAD FIBERS
#----------------------------------------

def STEEL_I_SECTION_QUAD(secID, matID,PLOT=True, DENSITY=7850/1e9):
    # input parameters
    # secID - section ID number
    # matID - material ID number 
    # d  = nominal depth
    # tw = web thickness
    # bf = flange width
    # tf = flange thickness
    # PLOT - whether to plot the section (default: True)
    # DENSITY - material density (default: 7850 kg/m^3 converted to N/mm^3)
    bf = 300               # [mm] Flange width
    tf = 20                # [mm] Flange thickness
    tw = 10                # [mm] Web thickness
    hw = 400               # [mm] Web height
    NUM = 100              # Number of fibers for web
    d = 2 * tf + hw 
    
    hw = d - 2 * tf  # Web height
    d = 2 * tf + hw  # Total depth (redundant, but ensures consistency)
    
    # Mass per Length
    AREA = 2 * bf * tf + tw * hw  # Total cross-sectional area
    MASS = DENSITY * AREA  # Mass per unit length
    
    # Number of fibers
    nfdw = 10  # Number of fibers along web depth
    nftw = 2   # Number of fibers along web thickness
    nfbf = 8   # Number of fibers along flange width
    nftf = 2   # Number of fibers along flange thickness
    
    # Calculate dimensions
    y1 = -d / 2  # Bottom of the section
    y2 = -hw / 2  # Bottom of the web
    y3 = hw / 2   # Top of the web
    y4 = d / 2    # Top of the section
    
    z1 = -bf / 2  # Left edge of the flange
    z2 = -tw / 2  # Left edge of the web
    z3 = tw / 2   # Right edge of the web
    z4 = bf / 2   # Right edge of the flange
    
    # Create the fiber section
    ops.section('Fiber', secID)
    
    # Define the patches for the flanges and web
    ops.patch('quad', matID, nfbf, nftf, y1, z4, y1, z1, y2, z1, y2, z4)  # Bottom flange
    ops.patch('quad', matID, nftw, nfdw, y2, z3, y2, z2, y3, z2, y3, z3)  # Web
    ops.patch('quad', matID, nfbf, nftf, y3, z4, y3, z1, y4, z1, y4, z4)  # Top flange
    
    # Plotting the section
    if PLOT:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 6))
        
        # Bottom flange
        bottom_flange = plt.Rectangle((z1, y1), bf, tf, color='grey', label='Bottom Flange')
        ax.add_patch(bottom_flange)
        
        # Web
        web = plt.Rectangle((z2, y2), tw, hw, color='grey', label='Web')
        ax.add_patch(web)
        
        # Top flange
        top_flange = plt.Rectangle((z1, y3), bf, tf, color='grey', label='Top Flange')
        ax.add_patch(top_flange)
        
        # Set plot limits
        ax.set_xlim(z1 - 1, z4 + 1)
        ax.set_ylim(y1 - 1, y4 + 1)
        
        # Labels and title
        ax.set_xlabel('z (width)')
        ax.set_ylabel('y (height)')
        ax.set_title(f'I-Section (ID: {secID})')
        ax.grid(True)
        ax.axhline(0, color='black', linewidth=0.5)
        ax.axvline(0, color='black', linewidth=0.5)
        ax.set_aspect('equal', adjustable='box')
        ax.legend()
        plt.show()
    
    return d, MASS

#---------------------------------------------------------------------------- 
