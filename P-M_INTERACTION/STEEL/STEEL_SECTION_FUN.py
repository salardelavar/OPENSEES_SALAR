#------------------------------------------
# DOUBLE I SECTIONS WITH AND WITHOUT PLATE
#------------------------------------------

def DOUBLE_I_SECTION(secTag, matTag, PLOT=True, DENSITY=7850/1e9):
    import openseespy.opensees as ops
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
