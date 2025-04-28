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

#-----------------------------------
# GREEK CROSS SECTION WITH FLANGES
#-----------------------------------

def GREEK_CROSS_SECTION_WITH_FLANGES(secTag, matTag, PLOT=True, DENSITY=7850/1e9):
    # Define geometric properties of the steel Greek section with flanges
    L = 300               # [mm] Length of each arm (tip-to-tip)
    t = 20                # [mm] Thickness of arms
    bf = 100              # [mm] Flange width
    tf = 15               # [mm] Flange thickness
    
    # Define geometric properties of the Greek cross section
    ops.section('Fiber', secTag)
    
    # Parameters
    half_L = L / 2       # Half length of each arm (total length = L)
    half_t = t / 2       # Half thickness of arms
    length_horizontal_part = half_L - half_t  # Length of horizontal arm segments
    length_vertical_part = half_L - half_t    # Length of vertical arm segments

    # Calculate area and mass
    AREA = 2 * L * t - t**2 + 4 * bf * tf  # Add area of flanges
    MASS = DENSITY * AREA

    # Number of fibers per arm segment
    NUM = 100

    # Horizontal arm fibers (left and right)
    if length_horizontal_part > 0:
        area_per_fiber_h = (length_horizontal_part * t) / NUM
        for i in range(NUM):
            # Left segment
            x = -half_L + (i + 0.5) * (length_horizontal_part / NUM)
            ops.fiber(x, 0.0, area_per_fiber_h, matTag)
            # Right segment
            x = half_t + (i + 0.5) * (length_horizontal_part / NUM)
            ops.fiber(x, 0.0, area_per_fiber_h, matTag)

    # Vertical arm fibers (upper and lower)
    if length_vertical_part > 0:
        area_per_fiber_v = (length_vertical_part * t) / NUM
        for i in range(NUM):
            # Lower segment
            y = -half_L + (i + 0.5) * (length_vertical_part / NUM)
            ops.fiber(0.0, y, area_per_fiber_v, matTag)
            # Upper segment
            y = half_t + (i + 0.5) * (length_vertical_part / NUM)
            ops.fiber(0.0, y, area_per_fiber_v, matTag)

    # Central square
    if t > 0:
        ops.fiber(0.0, 0.0, t**2, matTag)

    # Flanges
    # Top flange
    ops.fiber(0.0, half_L + tf / 2, bf * tf, matTag)
    # Bottom flange
    ops.fiber(0.0, -half_L - tf / 2, bf * tf, matTag)
    # Left flange
    ops.fiber(-half_L - tf / 2, 0.0, bf * tf, matTag)
    # Right flange
    ops.fiber(half_L + tf / 2, 0.0, bf * tf, matTag)

    # Plotting
    if PLOT:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 8))

        # Horizontal arms
        ax.add_patch(plt.Rectangle((-half_L, -half_t), length_horizontal_part, t, color='gray', alpha=0.7))
        ax.add_patch(plt.Rectangle((half_t, -half_t), length_horizontal_part, t, color='gray', alpha=0.7))

        # Vertical arms
        ax.add_patch(plt.Rectangle((-half_t, -half_L), t, length_vertical_part, color='gray', alpha=0.7))
        ax.add_patch(plt.Rectangle((-half_t, half_t), t, length_vertical_part, color='gray', alpha=0.7))

        # Central square
        ax.add_patch(plt.Rectangle((-half_t, -half_t), t, t, color='gray', alpha=0.7))

        # Flanges
        ax.add_patch(plt.Rectangle((-bf / 2, half_L), bf, tf, color='blue', alpha=0.5))  # Top flange
        ax.add_patch(plt.Rectangle((-bf / 2, -half_L - tf), bf, tf, color='blue', alpha=0.5))  # Bottom flange
        ax.add_patch(plt.Rectangle((-half_L - tf, -bf / 2), tf, bf, color='blue', alpha=0.5))  # Left flange
        ax.add_patch(plt.Rectangle((half_L, -bf / 2), tf, bf, color='blue', alpha=0.5))  # Right flange

        # Plot settings
        ax.set_xlim([-half_L - bf - 20, half_L + bf + 20])
        ax.set_ylim([-half_L - bf - 20, half_L + bf + 20])
        ax.set_aspect('equal')
        plt.title(f'Greek Cross Section with Flanges (L={L}mm, t={t}mm, bf={bf}mm, tf={tf}mm)')
        plt.xlabel('Width (mm)')
        plt.ylabel('Height (mm)')
        plt.grid(True)
        plt.show()

    return L + 2 * tf, MASS  # Return total section depth (including flanges) and mass per unit length

#----------------------------------------------------------------------------
    
#---------------------------
# 2 UNP SECTION WITH PLATES
#---------------------------

def TWO_UNP_SECTION_WITH_TOP_BOTTOM_PLATES(secTag, matTag, PLOT=True, DENSITY=7850/1e9, include_plates=True):

    #----------------------------
    # Geometry
    #----------------------------
    # UNP (U–profile) geometry for each channel
    H_unp    = 200.0    # [mm] Overall depth of each UNP shape
    L_flange = 50.0     # [mm] Length (horizontal extension) of each flange
    t_flange = 12.0     # [mm] Flange thickness
    t_web    = 8.0      # [mm] Web thickness
    spacing = 200       # [mm] Distance between UNP-section centers

    # Top and bottom plates (horizontal) that connect the flange regions:
    t_plate = 5.0           # [mm] Plate thickness
    plate_extension = -20.0 # [mm] Extra extension beyond the flange length on each side
    # The combined flange region of the two UNP’s spans from x = –L_flange to x = L_flange (i.e., 2*L_flange).
    # The plate will extend an extra plate_extension on each side.
    b_plate = 2 * L_flange + 2 * plate_extension  # overall plate length in x direction

    #----------------------------
    # Discretization Parameters
    #----------------------------
    NUM_web    = 20    # subdivisions along the web height
    NUM_flange = 10    # subdivisions along each flange (horizontal)
    NUM_plate  = 10    # subdivisions for the top and bottom plates (horizontal)

    #----------------------------
    # Compute Areas and Mass
    #----------------------------
    # For one UNP shape:
    AREA_web    = H_unp * t_web
    AREA_flange = L_flange * t_flange  # (each for top and bottom)
    AREA_unp    = AREA_web + 2 * AREA_flange

    # Two UNP shapes:
    AREA_total_unp = 2 * AREA_unp

    # Top and bottom plates:
    if include_plates:
        AREA_plate = b_plate * t_plate
        AREA_total_plate = 2 * AREA_plate
    else:
        AREA_total_plate = 0.0  # No plates

    AREA = AREA_total_unp + AREA_total_plate
    MASS = DENSITY * AREA

    #----------------------------
    # Define the section fibers using OpenSeesPy ops
    #----------------------------
    ops.section('Fiber', secTag)

    # --- Left UNP Shape ---
    # Left UNP Web: rectangle from x = –(L_flange+t_web) to x = –L_flange, y from –H_unp/2 to H_unp/2
    dx_web = t_web  # horizontal thickness of web
    dy_web = H_unp / NUM_web  # vertical subdivision
    x_web_left = -(L_flange + t_web/2)  # center of left web
    for i in range(NUM_web):
        y = -H_unp/2 + (i + 0.5) * dy_web
        area_fiber = dx_web * dy_web
        ops.fiber(x_web_left, y, area_fiber, matTag)

    # Left UNP Top Flange: rectangle from x = –L_flange to x = 0, y from H_unp/2 - t_flange to H_unp/2
    dx_flange = L_flange / NUM_flange
    x0_left_flange = -L_flange  # left edge of left flange
    y_top_flange = H_unp/2 - t_flange/2  # center in the vertical (flange is horizontal)
    for i in range(NUM_flange):
        x = x0_left_flange + (i + 0.5) * dx_flange
        area_fiber = dx_flange * t_flange
        ops.fiber(x, y_top_flange, area_fiber, matTag)

    # Left UNP Bottom Flange: rectangle from x = –L_flange to x = 0, y from –H_unp/2 to –H_unp/2 + t_flange
    y_bot_flange = -H_unp/2 + t_flange/2
    for i in range(NUM_flange):
        x = x0_left_flange + (i + 0.5) * dx_flange
        area_fiber = dx_flange * t_flange
        ops.fiber(x, y_bot_flange, area_fiber, matTag)

    # --- Right UNP Shape ---
    # Right UNP Web: rectangle from x = L_flange + spacing to x = L_flange + spacing + t_web, y from –H_unp/2 to H_unp/2
    x_web_right = L_flange + spacing + t_web/2  # center of right web (shifted by spacing)
    for i in range(NUM_web):
        y = -H_unp/2 + (i + 0.5) * dy_web
        area_fiber = dx_web * dy_web
        ops.fiber(x_web_right, y, area_fiber, matTag)

    # Right UNP Top Flange: rectangle from x = spacing to x = L_flange + spacing, y from H_unp/2 - t_flange to H_unp/2
    dx_flange = L_flange / NUM_flange
    x0_right_flange = spacing  # left edge of right flange (shifted by spacing)
    for i in range(NUM_flange):
        x = x0_right_flange + (i + 0.5) * dx_flange
        area_fiber = dx_flange * t_flange
        ops.fiber(x, y_top_flange, area_fiber, matTag)

    # Right UNP Bottom Flange: rectangle from x = spacing to x = L_flange + spacing, y from –H_unp/2 to –H_unp/2 + t_flange
    for i in range(NUM_flange):
        x = x0_right_flange + (i + 0.5) * dx_flange
        area_fiber = dx_flange * t_flange
        ops.fiber(x, y_bot_flange, area_fiber, matTag)

    # --- Top Plate ---
    if include_plates:
        # The top plate spans from x = –(L_flange + plate_extension) to x = L_flange + spacing + plate_extension.
        x_plate_left = -L_flange - plate_extension
        dx_plate = (b_plate + spacing) / NUM_plate  # adjusted for spacing
        y_plate_top = H_unp/2 + t_plate/2  # plate center in y (above the UNP top flanges)
        for i in range(NUM_plate):
            x = x_plate_left + (i + 0.5) * dx_plate
            area_fiber = dx_plate * t_plate
            ops.fiber(x, y_plate_top, area_fiber, matTag)

    # --- Bottom Plate ---
    if include_plates:
        # Similarly, the bottom plate spans the same horizontal extent.
        y_plate_bot = -H_unp/2 - t_plate/2
        for i in range(NUM_plate):
            x = x_plate_left + (i + 0.5) * dx_plate
            area_fiber = dx_plate * t_plate
            ops.fiber(x, y_plate_bot, area_fiber, matTag)

    #----------------------------
    # Plot the section for visualization
    #----------------------------
    if PLOT:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(10, 6))

        # --- Plot Left UNP ---
        # Left web
        rect_left_web = plt.Rectangle((-L_flange - t_web, -H_unp/2),
                                      t_web, H_unp,
                                      edgecolor='black', facecolor='gray', alpha=0.7)
        ax.add_patch(rect_left_web)
        # Left top flange
        rect_left_top = plt.Rectangle((-L_flange, H_unp/2 - t_flange),
                                      L_flange, t_flange,
                                      edgecolor='black', facecolor='gray', alpha=0.7)
        ax.add_patch(rect_left_top)
        # Left bottom flange
        rect_left_bot = plt.Rectangle((-L_flange, -H_unp/2),
                                      L_flange, t_flange,
                                      edgecolor='black', facecolor='gray', alpha=0.7)
        ax.add_patch(rect_left_bot)

        # --- Plot Right UNP ---
        # Right web
        rect_right_web = plt.Rectangle((L_flange + spacing, -H_unp/2),
                                       t_web, H_unp,
                                       edgecolor='black', facecolor='gray', alpha=0.7)
        ax.add_patch(rect_right_web)
        # Right top flange
        rect_right_top = plt.Rectangle((spacing, H_unp/2 - t_flange),
                                       L_flange, t_flange,
                                       edgecolor='black', facecolor='gray', alpha=0.7)
        ax.add_patch(rect_right_top)
        # Right bottom flange
        rect_right_bot = plt.Rectangle((spacing, -H_unp/2),
                                       L_flange, t_flange,
                                       edgecolor='black', facecolor='gray', alpha=0.7)
        ax.add_patch(rect_right_bot)

        # --- Plot Top Plate ---
        if include_plates:
            rect_top_plate = plt.Rectangle((x_plate_left, H_unp/2),
                                           b_plate + spacing, t_plate,
                                           edgecolor='black', facecolor='blue', alpha=0.7)
            ax.add_patch(rect_top_plate)
        # --- Plot Bottom Plate ---
        if include_plates:
            rect_bot_plate = plt.Rectangle((x_plate_left, -H_unp/2 - t_plate),
                                           b_plate + spacing, t_plate,
                                           edgecolor='black', facecolor='blue', alpha=0.7)
            ax.add_patch(rect_bot_plate)

        # Draw an outer dashed outline of the composite section:
        outer_left = - (L_flange + t_web)
        outer_right = L_flange + spacing + t_web
        outer_bottom = -H_unp/2 - (t_plate if include_plates else 0)
        outer_top = H_unp/2 + (t_plate if include_plates else 0)
        outer_width = outer_right - outer_left
        outer_height = outer_top - outer_bottom
        outer_outline = plt.Rectangle((outer_left, outer_bottom),
                                      outer_width, outer_height,
                                      edgecolor='black', facecolor='none', linestyle='--', linewidth=1.5)
        ax.add_patch(outer_outline)

        ax.set_xlim(outer_left - 20, outer_right + 20)
        ax.set_ylim(outer_bottom - 20, outer_top + 20)
        ax.set_aspect('equal', adjustable='box')
        plt.xlabel('Width (mm)', fontsize=12)
        plt.ylabel('Height (mm)', fontsize=12)
        plt.title('Composite Section: 2 UNP Shapes' + (' with Top & Bottom Plates' if include_plates else ''), fontsize=14)
        plt.grid(True)
        plt.show()

    # Overall section height (including plates if applicable)
    H_total = H_unp + (2 * t_plate if include_plates else 0)
    return H_total, MASS
    
#----------------------------------------------------------------------------   
