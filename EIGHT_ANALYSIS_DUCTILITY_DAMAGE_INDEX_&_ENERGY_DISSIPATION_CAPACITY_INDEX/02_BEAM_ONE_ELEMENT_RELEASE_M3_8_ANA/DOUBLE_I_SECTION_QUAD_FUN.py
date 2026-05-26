import numpy as np
import openseespy.opensees as ops

#-----------------------------------------------------
# DOUBLE I SECTIONS WITH AND WITHOUT PLATE WITH QUAD
#-----------------------------------------------------

def DOUBLE_I_SECTION_QUAD(secTag, matTag, PLOT=True, DENSITY=7850/1e9):
    """
    Create a Double I-Section with optional top/bottom plates.
    Uses quad patches instead of individual fibers.
    Returns: (section_depth, mass_per_length)
    """
    ops.section('Fiber', secTag)

    # Section geometry
    bf = 300        # [mm] flange width
    tf = 20         # [mm] flange thickness
    tw = 10         # [mm] web thickness
    hw = 400        # [mm] clear web height
    spacing = 500   # [mm] distance between two I-sections (center to center)
    plate_t = 10    # [mm] Plate thickness
    plate_l = 500   # [mm] Plate length
    nf_web = 30     # mesh divisions
    nf_flange = 12

    # Total section depth
    d = hw + 2*(tf + plate_t)

    # Area for mass
    AREA = 2*(2 * bf * tf + tw * hw) + 2*(plate_l * plate_t)
    MASS = AREA * DENSITY

    # -----------------------
    #  Create two I-sections
    # -----------------------
    for y0 in [-spacing/2, spacing/2]:

        # --- Web patch ---
        ops.patch('quad', matTag, nf_web, nf_web,
            y0 - tw/2, -hw/2,      # lower-left
            y0 + tw/2, -hw/2,      # lower-right
            y0 + tw/2,  hw/2,      # upper-right
            y0 - tw/2,  hw/2       # upper-left
        )

        # --- Top flange patch ---
        ops.patch('quad', matTag, nf_flange, nf_flange,
            y0 - bf/2,  hw/2, 
            y0 + bf/2,  hw/2,
            y0 + bf/2,  hw/2 + tf,
            y0 - bf/2,  hw/2 + tf
        )

        # --- Bottom flange patch ---
        ops.patch('quad', matTag, nf_flange, nf_flange,
            y0 - bf/2, -hw/2 - tf,
            y0 + bf/2, -hw/2 - tf,
            y0 + bf/2, -hw/2,
            y0 - bf/2, -hw/2
        )

    # -----------------------
    #  Add plates if needed
    # -----------------------
    if plate_t > 0 and plate_l > 0:

        # Top plate
        ops.patch('quad', matTag, nf_flange, nf_flange,
            -plate_l/2,  hw/2 + tf,
             plate_l/2,  hw/2 + tf,
             plate_l/2,  hw/2 + tf + plate_t,
            -plate_l/2,  hw/2 + tf + plate_t
        )

        # Bottom plate
        ops.patch('quad', matTag, nf_flange, nf_flange,
            -plate_l/2, -hw/2 - tf - plate_t,
             plate_l/2, -hw/2 - tf - plate_t,
             plate_l/2, -hw/2 - tf,
            -plate_l/2, -hw/2 - tf
        )

    # -----------------------
    # Plot section
    # -----------------------
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

    return d, MASS
