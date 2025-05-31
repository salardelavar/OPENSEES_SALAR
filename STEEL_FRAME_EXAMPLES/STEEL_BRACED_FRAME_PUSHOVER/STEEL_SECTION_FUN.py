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
        ax.add_patch(plt.Rectangle((-bf / 2, hw / 2), bf, tf, color='blue', alpha=0.7))

        # Plot the bottom flange fibers
        ax.add_patch(plt.Rectangle((-bf / 2, -hw / 2 - tf), bf, tf, color='blue', alpha=0.7))

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
                                     color='blue', alpha=0.7))
            ax.add_patch(plt.Rectangle((y_shift - bf/2, -hw/2 - tf), bf, tf, 
                                     color='blue', alpha=0.7))
            # Web
            ax.add_patch(plt.Rectangle((y_shift - tw/2, -hw/2), tw, hw, 
                                     color='grey', alpha=0.7))
        
        # Draw plates
        if plate_t > 0 and plate_l > 0:
            ax.add_patch(plt.Rectangle((-plate_l/2, hw/2 + tf), plate_l, plate_t,
                                     color='cyan', alpha=0.5))
            ax.add_patch(plt.Rectangle((-plate_l/2, -hw/2 - tf - plate_t), plate_l, plate_t,
                                     color='cyan', alpha=0.5))
        
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
    
#----------------------------
# RECTANGULAR TUBE SECTION 
#----------------------------
def R_TUBE_SECTION(secTag, matTag, PLOT=True, DENSITY=7850/1e9):
    # Define section tag for rectangular tube steel section
    #secTag = 1
    # Define geometric properties of the rectangular tube
    B = 300  # [mm] Outer width of the tube
    H = 500  # [mm] Outer height of the tube
    t = 20   # [mm] Wall thickness of the tube

    # Number of fibers along each wall direction
    NUM_B = 10  # Number of fibers along the width
    NUM_H = 20  # Number of fibers along the height
    
    # Mass per Length
    AREA = 2 * ((B - 2*t) * t) + 2 *  (H * t)
    MASS = DENSITY * AREA  
    
    # Define the rectangular tube section using fibers
    ops.section('Fiber', secTag)
    # Outer top wall fibers
    for i in range(NUM_B):
        x_loc = -B / 2 + (i + 0.5) * (B / NUM_B)
        ops.fiber(x_loc, H / 2 - t / 2, (B / NUM_B) * t, matTag)
    # Outer bottom wall fibers
    for i in range(NUM_B):
        x_loc = -B / 2 + (i + 0.5) * (B / NUM_B)
        ops.fiber(x_loc, -H / 2 + t / 2, (B / NUM_B) * t, matTag)
    # Outer left wall fibers
    for i in range(NUM_H):
        y_loc = -H / 2 + (i + 0.5) * (H / NUM_H)
        ops.fiber(-B / 2 + t / 2, y_loc, t * (H / NUM_H), matTag)
    # Outer right wall fibers
    for i in range(NUM_H):
        y_loc = -H / 2 + (i + 0.5) * (H / NUM_H)
        ops.fiber(B / 2 - t / 2, y_loc, t * (H / NUM_H), matTag)
        
    
    #-------------------
    # PLOT THE SECTION
    #-------------------
    if PLOT == True:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 8))

        # Plot outer top wall fibers
        fiber_width = B / NUM_B
        for i in range(NUM_B):
            x_loc = -B / 2 + i * fiber_width
            ax.add_patch(plt.Rectangle((x_loc, H / 2 - t), fiber_width, t, color='grey', alpha=0.7))

        # Plot outer bottom wall fibers
        for i in range(NUM_B):
            x_loc = -B / 2 + i * fiber_width
            ax.add_patch(plt.Rectangle((x_loc, -H / 2), fiber_width, t, color='grey', alpha=0.7))

        # Plot outer left wall fibers
        fiber_height = H / NUM_H
        for i in range(NUM_H):
            y_loc = -H / 2 + i * fiber_height
            ax.add_patch(plt.Rectangle((-B / 2, y_loc), t, fiber_height, color='grey', alpha=0.7))

        # Plot outer right wall fibers
        for i in range(NUM_H):
            y_loc = -H / 2 + i * fiber_height
            ax.add_patch(plt.Rectangle((B / 2 - t, y_loc), t, fiber_height, color='grey', alpha=0.7))

        # Add section outline
        outer_outline = plt.Rectangle((-B / 2, -H / 2), B, H, edgecolor='black', fill=False, linewidth=1.5)
        inner_outline = plt.Rectangle((-B / 2 + t, -H / 2 + t), B - 2 * t, H - 2 * t, edgecolor='black', fill=False, linewidth=1.5)
        ax.add_patch(outer_outline)
        ax.add_patch(inner_outline)

        # Set plot limits
        ax.set_xlim([-B / 2 - 20, B / 2 + 20])
        ax.set_ylim([-H / 2 - 20, H / 2 + 20])
        ax.set_aspect('equal', adjustable='box')

        # Labels and title
        plt.xlabel('Width (B)', fontsize=12)
        plt.ylabel('Height (H)', fontsize=12)
        plt.title('Rectangular Tube Section with Fibers', fontsize=14)

        # Show plot
        plt.grid(True)
        plt.show()
    
    return H, MASS# Return Section Height

#----------------------------------------------------------------------------
    
#----------------------------    
# CIRCULAR TUBE SECTION  
#---------------------------- 
def C_TUBE_SECTION(secTag, matTag, PLOT=True, DENSITY=7850/1e9):
    import numpy as np
    # Define geometric properties of the steel tube
    #secTag = 1
    D = 500          # [mm] Outer diameter of the tube
    t = 20           # [mm] Wall thickness of the tube
    r = D / 2        # [mm] Outer radius of the tube (m)
    r_inner = r - t  # [mm] Inner radius of the tube

    # Define material tag
    matTag = 1

    # Number of fibers along the radial and circumferential directions
    NUM_R = 10  # Number of fibers radially
    NUM_C = 20  # Number of fibers circumferentially
    
    # Mass per Length
    AREA = (np.pi * r**2) - (np.pi * r_inner**2)
    MASS = DENSITY * AREA 

    # Define the steel tube section using fibers
    ops.section('Fiber', secTag)

    # Loop over radial and circumferential directions to define fibers
    import numpy as np
    for i in range(NUM_R):
        # Compute the radial location of the fiber
        r_loc = r_inner + (i + 0.5) * (t / NUM_R)
        for j in range(NUM_C):
            # Compute the angular location of the fiber
            theta = j * (2 * np.pi / NUM_C)
            x_loc = r_loc * np.cos(theta)
            y_loc = r_loc * np.sin(theta)
            # Compute the area of the fiber
            fiber_area = (2 * np.pi * r_loc / NUM_C) * (t / NUM_R)
            # Add the fiber to the section
            ops.fiber(x_loc, y_loc, fiber_area, matTag)
    
    #-------------------
    # PLOT THE SECTION
    #-------------------
    if PLOT == True:
        import matplotlib.pyplot as plt
        # Calculate radii
        r = D / 2  # Outer radius
        r_inner = r - t  # Inner radius

        # Initialize the plot
        fig, ax = plt.subplots(figsize=(8, 8))

        # Plot fibers
        for i in range(NUM_R):
            # Radial location of fibers
            r_loc = r_inner + (i + 0.5) * (t / NUM_R)
            for j in range(NUM_C):
                # Angular location of fibers
                theta = j * (2 * np.pi / NUM_C)
                x_loc = r_loc * np.cos(theta)
                y_loc = r_loc * np.sin(theta)
                ax.add_patch(plt.Circle((x_loc, y_loc), 2, color='grey', alpha=0.7))  # Fibers as small circles

        # Add the outer and inner circles for reference
        outer_circle = plt.Circle((0, 0), r, color='black', fill=False, linewidth=1.5)
        inner_circle = plt.Circle((0, 0), r_inner, color='black', fill=False, linewidth=1.5)
        ax.add_patch(outer_circle)
        ax.add_patch(inner_circle)

        # Set aspect ratio and limits
        ax.set_aspect('equal', adjustable='box')
        ax.set_xlim([-r - 20, r + 20])
        ax.set_ylim([-r - 20, r + 20])

        # Labels and title
        plt.xlabel('X (mm)', fontsize=12)
        plt.ylabel('Y (mm)', fontsize=12)
        plt.title('Circular Tube Section with Fibers', fontsize=14)

        # Show the plot
        plt.grid(True)
        plt.show()
    
    return D, MASS # Return Section Height

#----------------------------------------------------------------------------

#---------------------------
# RECTANGULAR SECTION 
#---------------------------
def R_RECTANGULAR_STEEL_SECTION(secTag, matTag, PLOT=True, DENSITY=7850/1e9):
    # Define a rectangular steel section using fibers.
    B = 400       # [mm] Width of the rectangular section
    H = 500       # [mm] Height of the rectangular section
    NUM_B = 10   # Number of fibers along the width of the section
    NUM_H = 10   # Number of fibers along the height of the section
    
    # Mass per Length
    AREA = B * H
    MASS = DENSITY * AREA  
    
    secTag = 1    # Section tag identifier
    matTag = 1    # Material tag for the concrete
    ops.section('Fiber', secTag)
    # Create fibers for the entire section
    for i in range(NUM_B):
        for j in range(NUM_H):
            # Calculate the x and y location of the fiber centroid
            x_loc = -B / 2 + (i + 0.5) * (B / NUM_B)
            y_loc = -H / 2 + (j + 0.5) * (H / NUM_H)
            # Define the fiber area
            fiber_area = (B / NUM_B) * (H / NUM_H)
            # Add the fiber to the section
            ops.fiber(x_loc, y_loc, fiber_area, matTag)
            
    #-------------------
    # PLOT THE SECTION
    #------------------- 
    if PLOT == True:
        import matplotlib.pyplot as plt
        fiber_width = B / NUM_B
        fiber_height = H / NUM_H
        
        # Prepare the plot
        fig, ax = plt.subplots()
        ax.set_aspect('equal', 'box')
        ax.set_xlim(-B / 2, B / 2)
        ax.set_ylim(-H / 2, H / 2)
        
        # Plot the fibers
        for i in range(NUM_B):
            for j in range(NUM_H):
                # Calculate the x and y location of the fiber centroid
                x_loc = -B / 2 + i * fiber_width
                y_loc = -H / 2 + j * fiber_height
                # Draw the rectangle representing the fiber
                rect = plt.Rectangle((x_loc, y_loc), fiber_width, fiber_height, color='grey', edgecolor='black', alpha=0.7)
                ax.add_patch(rect)
        
        # Add labels and grid
        ax.set_xlabel("Width [mm]")
        ax.set_ylabel("Height [mm]")
        ax.grid(True)
        ax.set_title("Rectangular Steel Section with Fibers")
        
        # Show the plot
        plt.show()           

    return H, MASS # Return Section Height

#----------------------------------------------------------------------------

#-------------------------
# BOX SECTION WITH ANGLES
#-------------------------
def BOX_SECTION_WITH_JUST_ANGLES(secTag, matTag, PLOT=True, DENSITY=7850/1e9):
    # ----------------------------
    # Geometry
    # ----------------------------
    B = 300         # [mm] Outer width of the box
    H = 500         # [mm] Outer height of the box
    t_angle = 10    # [mm] Thickness of the angle sections
    L_angle = 50    # [mm] Length of the angle legs

    # Number of subdivisions along each leg for fiber discretization
    NUM_leg = 10  # for both horizontal and vertical legs

    # Compute total cross-sectional area (each corner is L-shaped)
    # Area of one corner: horizontal leg (L_angle*t_angle) + vertical leg ((L_angle-t_angle)*t_angle)
    AREA_corner = 2 * L_angle * t_angle - t_angle**2
    AREA = 4 * AREA_corner
    MASS = DENSITY * AREA

    # ----------------------------
    # Define section fibers (using OpenSeesPy ops module)
    # ----------------------------
    ops.section('Fiber', secTag)
    
    # Bottom-left corner:
    # Horizontal leg: from x = -B/2 to -B/2 + L_angle, at constant y = -H/2 + t_angle/2
    dx = L_angle / NUM_leg
    for i in range(NUM_leg):
        x = -B/2 + (i + 0.5) * dx
        y = -H/2 + t_angle/2
        area_fiber = dx * t_angle
        ops.fiber(x, y, area_fiber, matTag)
    # Vertical leg: from y = -H/2 + t_angle to -H/2 + L_angle, at constant x = -B/2 + t_angle/2
    dy = (L_angle - t_angle) / NUM_leg
    for i in range(NUM_leg):
        y = -H/2 + t_angle + (i + 0.5) * dy
        x = -B/2 + t_angle/2
        area_fiber = dy * t_angle
        ops.fiber(x, y, area_fiber, matTag)
        
    # Bottom-right corner:
    # Horizontal leg: from x = B/2 - L_angle to B/2, at constant y = -H/2 + t_angle/2
    for i in range(NUM_leg):
        x = B/2 - L_angle + (i + 0.5) * dx
        y = -H/2 + t_angle/2
        area_fiber = dx * t_angle
        ops.fiber(x, y, area_fiber, matTag)
    # Vertical leg: from y = -H/2 + t_angle to -H/2 + L_angle, at constant x = B/2 - t_angle/2
    for i in range(NUM_leg):
        y = -H/2 + t_angle + (i + 0.5) * dy
        x = B/2 - t_angle/2
        area_fiber = dy * t_angle
        ops.fiber(x, y, area_fiber, matTag)
        
    # Top-left corner:
    # Horizontal leg: from x = -B/2 to -B/2 + L_angle, at constant y = H/2 - t_angle/2
    for i in range(NUM_leg):
        x = -B/2 + (i + 0.5) * dx
        y = H/2 - t_angle/2
        area_fiber = dx * t_angle
        ops.fiber(x, y, area_fiber, matTag)
    # Vertical leg: from y = H/2 - L_angle to H/2 - t_angle, at constant x = -B/2 + t_angle/2
    for i in range(NUM_leg):
        y = H/2 - L_angle + (i + 0.5) * dy
        x = -B/2 + t_angle/2
        area_fiber = dy * t_angle
        ops.fiber(x, y, area_fiber, matTag)
        
    # Top-right corner:
    # Horizontal leg: from x = B/2 - L_angle to B/2, at constant y = H/2 - t_angle/2
    for i in range(NUM_leg):
        x = B/2 - L_angle + (i + 0.5) * dx
        y = H/2 - t_angle/2
        area_fiber = dx * t_angle
        ops.fiber(x, y, area_fiber, matTag)
    # Vertical leg: from y = H/2 - L_angle to H/2 - t_angle, at constant x = B/2 - t_angle/2
    for i in range(NUM_leg):
        y = H/2 - L_angle + (i + 0.5) * dy
        x = B/2 - t_angle/2
        area_fiber = dy * t_angle
        ops.fiber(x, y, area_fiber, matTag)

    # ----------------------------
    # Plot the section for visualization
    # ----------------------------
    if PLOT:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 8))
        
        # Draw each corner's two legs as patches.
        # Bottom-left corner:
        # Horizontal leg:
        rect_bl_h = plt.Rectangle((-B/2, -H/2), L_angle, t_angle, color='grey', alpha=0.7)
        # Vertical leg:
        rect_bl_v = plt.Rectangle((-B/2, -H/2 + t_angle), t_angle, L_angle - t_angle, color='grey', alpha=0.7)
        ax.add_patch(rect_bl_h)
        ax.add_patch(rect_bl_v)
        
        # Bottom-right corner:
        rect_br_h = plt.Rectangle((B/2 - L_angle, -H/2), L_angle, t_angle, color='grey', alpha=0.7)
        rect_br_v = plt.Rectangle((B/2 - t_angle, -H/2 + t_angle), t_angle, L_angle - t_angle, color='grey', alpha=0.7)
        ax.add_patch(rect_br_h)
        ax.add_patch(rect_br_v)
        
        # Top-left corner:
        rect_tl_h = plt.Rectangle((-B/2, H/2 - t_angle), L_angle, t_angle, color='grey', alpha=0.7)
        rect_tl_v = plt.Rectangle((-B/2, H/2 - L_angle), t_angle, L_angle - t_angle, color='grey', alpha=0.7)
        ax.add_patch(rect_tl_h)
        ax.add_patch(rect_tl_v)
        
        # Top-right corner:
        rect_tr_h = plt.Rectangle((B/2 - L_angle, H/2 - t_angle), L_angle, t_angle, color='grey', alpha=0.7)
        rect_tr_v = plt.Rectangle((B/2 - t_angle, H/2 - L_angle), t_angle, L_angle - t_angle, color='grey', alpha=0.7)
        ax.add_patch(rect_tr_h)
        ax.add_patch(rect_tr_v)
        
        # Draw the outer boundary of the box section
        outer_outline = plt.Rectangle((-B/2, -H/2), B, H, edgecolor='black', fill=False, linewidth=1.5)
        ax.add_patch(outer_outline)
        
        ax.set_xlim([-B/2 - 20, B/2 + 20])
        ax.set_ylim([-H/2 - 20, H/2 + 20])
        ax.set_aspect('equal', adjustable='box')
        plt.xlabel('Width (mm)', fontsize=12)
        plt.ylabel('Height (mm)', fontsize=12)
        plt.title('Box Section with Four Angle (L-Shaped) Sections', fontsize=14)
        plt.grid(True)
        plt.show()
    
    return H, MASS

#----------------------------------------------------------------------------

#-------------------------------------
# BOX SECTION WITH ANGLES AND PLATES
#-------------------------------------    
def BOX_SECTION_WITH_ANGLES_PLATES(secTag, matTag, PLOT=True, DENSITY=7850/1e9):
    #----------------------------
    # Geometry
    #----------------------------
    # Overall (dummy) outer dimensions of the box section
    B = 300  # [mm] Outer width of the box
    H = 500  # [mm] Outer height of the box

    # Angle (corner) parameters
    t_angle = 10   # [mm] Thickness of the angle sections
    L_angle = 50   # [mm] Length of the angle legs

    # Plate parameters (clear lengths and thickness)
    t_plate = 5                # [mm] Plate thickness
    b_plate = B-L_angle*2 + 20 # [mm] Clear horizontal plate length (top and bottom)
    h_plate = H-L_angle*2 + 20 # [mm] Clear vertical plate length (left and right)

    #----------------------------
    # Discretization Parameters
    #----------------------------
    NUM_leg = 10      # subdivisions for each leg of the angle sections
    NUM_plate_x = 10  # subdivisions for horizontal plates (top & bottom)
    NUM_plate_y = 20  # subdivisions for vertical plates (left & right)

    #----------------------------
    # Compute Areas
    #----------------------------
    # Each angle (corner) is L-shaped:
    #   Horizontal leg area: L_angle * t_angle
    #   Vertical leg area: (L_angle - t_angle) * t_angle
    #   Combined area per corner (avoiding double-counting the overlap):
    AREA_corner = L_angle * t_angle + (L_angle - t_angle) * t_angle
    # Total area for four corners:
    AREA_angles = 4 * AREA_corner

    # Plate areas:
    AREA_plate_horizontal = b_plate * t_plate  # top or bottom
    AREA_plate_vertical   = h_plate * t_plate   # left or right

    # Total area (angles + plates)
    AREA = AREA_angles + 2 * AREA_plate_horizontal + 2 * AREA_plate_vertical
    MASS = DENSITY * AREA

    #----------------------------
    # Define the section fibers using OpenSeesPy ops
    #----------------------------
    ops.section('Fiber', secTag)

    # 1. Define fibers for the four angle (corner) sections.
    # Discretization for the angle legs:
    dx_angle = L_angle / NUM_leg         # for horizontal leg
    dy_angle = (L_angle - t_angle) / NUM_leg  # for vertical leg

    # --- Bottom-Left Corner ---
    # Horizontal leg: from x = -B/2 to -B/2 + L_angle, at constant y = -H/2 + t_angle/2
    for i in range(NUM_leg):
        x = -B/2 + (i + 0.5) * dx_angle
        y = -H/2 + t_angle/2
        area_fiber = dx_angle * t_angle
        ops.fiber(x, y, area_fiber, matTag)
    # Vertical leg: from y = -H/2 + t_angle to -H/2 + L_angle, at constant x = -B/2 + t_angle/2
    for i in range(NUM_leg):
        y = -H/2 + t_angle + (i + 0.5) * dy_angle
        x = -B/2 + t_angle/2
        area_fiber = dy_angle * t_angle
        ops.fiber(x, y, area_fiber, matTag)

    # --- Bottom-Right Corner ---
    # Horizontal leg: from x = B/2 - L_angle to B/2, at constant y = -H/2 + t_angle/2
    for i in range(NUM_leg):
        x = B/2 - L_angle + (i + 0.5) * dx_angle
        y = -H/2 + t_angle/2
        area_fiber = dx_angle * t_angle
        ops.fiber(x, y, area_fiber, matTag)
    # Vertical leg: from y = -H/2 + t_angle to -H/2 + L_angle, at constant x = B/2 - t_angle/2
    for i in range(NUM_leg):
        y = -H/2 + t_angle + (i + 0.5) * dy_angle
        x = B/2 - t_angle/2
        area_fiber = dy_angle * t_angle
        ops.fiber(x, y, area_fiber, matTag)

    # --- Top-Left Corner ---
    # Horizontal leg: from x = -B/2 to -B/2 + L_angle, at constant y = H/2 - t_angle/2
    for i in range(NUM_leg):
        x = -B/2 + (i + 0.5) * dx_angle
        y = H/2 - t_angle/2
        area_fiber = dx_angle * t_angle
        ops.fiber(x, y, area_fiber, matTag)
    # Vertical leg: from y = H/2 - L_angle to H/2 - t_angle, at constant x = -B/2 + t_angle/2
    for i in range(NUM_leg):
        y = H/2 - L_angle + (i + 0.5) * dy_angle
        x = -B/2 + t_angle/2
        area_fiber = dy_angle * t_angle
        ops.fiber(x, y, area_fiber, matTag)

    # --- Top-Right Corner ---
    # Horizontal leg: from x = B/2 - L_angle to B/2, at constant y = H/2 - t_angle/2
    for i in range(NUM_leg):
        x = B/2 - L_angle + (i + 0.5) * dx_angle
        y = H/2 - t_angle/2
        area_fiber = dx_angle * t_angle
        ops.fiber(x, y, area_fiber, matTag)
    # Vertical leg: from y = H/2 - L_angle to H/2 - t_angle, at constant x = B/2 - t_angle/2
    for i in range(NUM_leg):
        y = H/2 - L_angle + (i + 0.5) * dy_angle
        x = B/2 - t_angle/2
        area_fiber = dy_angle * t_angle
        ops.fiber(x, y, area_fiber, matTag)

    # 2. Define fibers for the plate regions.
    # -- Horizontal Plates (Top and Bottom) --
    dx_plate = b_plate / NUM_plate_x
    # Top plate: centered horizontally on the top face.
    for i in range(NUM_plate_x):
        x = -b_plate/2 + (i + 0.5) * dx_plate
        y = H/2 + t_plate/2
        area_fiber = dx_plate * t_plate
        ops.fiber(x, y, area_fiber, matTag)
    # Bottom plate:
    for i in range(NUM_plate_x):
        x = -b_plate/2 + (i + 0.5) * dx_plate
        y = -H/2 - t_plate/2
        area_fiber = dx_plate * t_plate
        ops.fiber(x, y, area_fiber, matTag)

    # -- Vertical Plates (Left and Right) --
    dy_plate = h_plate / NUM_plate_y
    # Left plate: centered vertically on the left face.
    for i in range(NUM_plate_y):
        y = -h_plate/2 + (i + 0.5) * dy_plate
        x = -B/2 - t_plate/2
        area_fiber = dy_plate * t_plate
        ops.fiber(x, y, area_fiber, matTag)
    # Right plate:
    for i in range(NUM_plate_y):
        y = -h_plate/2 + (i + 0.5) * dy_plate
        x = B/2 + t_plate/2
        area_fiber = dy_plate * t_plate
        ops.fiber(x, y, area_fiber, matTag)

    #----------------------------
    # Plot the section for visualization
    #----------------------------
    if PLOT:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 8))
        
        # Plot Angle (Corner) Sections as L-shaped patches.
        # Bottom-Left Corner:
        rect_bl_h = plt.Rectangle((-B/2, -H/2), L_angle, t_angle, color='grey', alpha=0.7)
        rect_bl_v = plt.Rectangle((-B/2, -H/2 + t_angle), t_angle, L_angle - t_angle, color='grey', alpha=0.7)
        ax.add_patch(rect_bl_h)
        ax.add_patch(rect_bl_v)
        # Bottom-Right Corner:
        rect_br_h = plt.Rectangle((B/2 - L_angle, -H/2), L_angle, t_angle, color='grey', alpha=0.7)
        rect_br_v = plt.Rectangle((B/2 - t_angle, -H/2 + t_angle), t_angle, L_angle - t_angle, color='grey', alpha=0.7)
        ax.add_patch(rect_br_h)
        ax.add_patch(rect_br_v)
        # Top-Left Corner:
        rect_tl_h = plt.Rectangle((-B/2, H/2 - t_angle), L_angle, t_angle, color='grey', alpha=0.7)
        rect_tl_v = plt.Rectangle((-B/2, H/2 - L_angle), t_angle, L_angle - t_angle, color='grey', alpha=0.7)
        ax.add_patch(rect_tl_h)
        ax.add_patch(rect_tl_v)
        # Top-Right Corner:
        rect_tr_h = plt.Rectangle((B/2 - L_angle, H/2 - t_angle), L_angle, t_angle, color='grey', alpha=0.7)
        rect_tr_v = plt.Rectangle((B/2 - t_angle, H/2 - L_angle), t_angle, L_angle - t_angle, color='grey', alpha=0.7)
        ax.add_patch(rect_tr_h)
        ax.add_patch(rect_tr_v)
        
        # Plot Plate Regions:
        # Top Plate:
        rect_top = plt.Rectangle((-b_plate/2, H/2), b_plate, t_plate, color='blue', alpha=0.7)
        ax.add_patch(rect_top)
        # Bottom Plate:
        rect_bot = plt.Rectangle((-b_plate/2, -H/2 - t_plate), b_plate, t_plate, color='blue', alpha=0.7)
        ax.add_patch(rect_bot)
        # Left Plate:
        rect_left = plt.Rectangle((-B/2-t_plate, -h_plate/2), t_plate, h_plate, color='blue', alpha=0.7)
        ax.add_patch(rect_left)
        # Right Plate:
        rect_right = plt.Rectangle((B/2, -h_plate/2), t_plate, h_plate, color='blue', alpha=0.7)
        ax.add_patch(rect_right)
        
        # Outer boundary of the box section:
        outer_outline = plt.Rectangle((-B/2, -H/2), B, H, edgecolor='black', fill=False, linewidth=1.5)
        ax.add_patch(outer_outline)
        
        ax.set_xlim([-B/2 - 20, B/2 + 20])
        ax.set_ylim([-H/2 - 20, H/2 + 20])
        ax.set_aspect('equal', adjustable='box')
        plt.xlabel('Width (mm)', fontsize=12)
        plt.ylabel('Height (mm)', fontsize=12)
        plt.title('Box Section with Four L-Shaped Angles and Plates', fontsize=14)
        plt.grid(True)
        plt.show()
    
    return H, MASS

#----------------------------------------------------------------------------

#-----------------------
# GREEK CROSS SECTION
#-----------------------

def GREEK_CROSS_SECTION(secTag, matTag, PLOT=True, DENSITY=7850/1e9):
    # Define geometric properties of the steel Greek section
    L = 300               # [mm] Flange width
    t = 20                # [mm] Flange thickness
    
    # Define geometric properties of the Greek cross section
    ops.section('Fiber', secTag)
    
    # Parameters
    half_L = L / 2       # Half length of each arm (total length = L)
    half_t = t / 2       # Half thickness of arms
    length_horizontal_part = half_L - half_t  # Length of horizontal arm segments
    length_vertical_part = half_L - half_t    # Length of vertical arm segments

    # Calculate area and mass
    AREA = 2 * L * t - t**2
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

        # Plot settings
        ax.set_xlim([-half_L-20, half_L+20])
        ax.set_ylim([-half_L-20, half_L+20])
        ax.set_aspect('equal')
        plt.title(f'Greek Cross Section (L={L}mm, t={t}mm)')
        plt.xlabel('Width (mm)')
        plt.ylabel('Height (mm)')
        plt.grid(True)
        plt.show()

    return L, MASS  # Return total section depth and mass per unit length

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
    H_unp    = 100.0    # [mm] Overall depth of each UNP shape
    L_flange = 50.0     # [mm] Length (horizontal extension) of each flange
    t_flange = 6.0      # [mm] Flange thickness
    t_web    = 7.0      # [mm] Web thickness
    spacing = 50        # [mm] Distance between UNP-section centers

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

#----------------------------------------    
# STEEL FIBER SECTION WITH QUAD FIBERS
#----------------------------------------

def STEEL_I_SECTION_QUAD(secID, matID, PLOT=True, DENSITY=7850/1e9):
    # input parameters
    # secID - section ID number
    # matID - material ID number 
    # d  = nominal depth
    # tw = web thickness
    # bf = flange width
    # tf = flange thickness
    # PLOT - whether to plot the section (default: True)
    # DENSITY - material density (default: 7850 kg/m^3 converted to N/mm^3)
    bf = 300               # [mm] Flange width per I-section
    tf = 20                # [mm] Flange thickness
    tw = 10                # [mm] Web thickness
    hw = 400               # [mm] Web height
    d = hw + 2 * tf
    
    #hw = d - 2 * tf  # Web height
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
        bottom_flange = plt.Rectangle((z1, y1), bf, tf, color='blue', label='Bottom Flange')
        ax.add_patch(bottom_flange)
        
        # Web
        web = plt.Rectangle((z2, y2), tw, hw, color='grey', label='Web')
        ax.add_patch(web)
        
        # Top flange
        top_flange = plt.Rectangle((z1, y3), bf, tf, color='blue', label='Top Flange')
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

#--------------------------------
# RECTANGULAR TUBE SECTION QUAD
#--------------------------------

def R_TUBE_SECTION_QUAD(secTag, matTag, PLOT=True, DENSITY=7850/1e9):
    # Define geometric properties of the rectangular tube
    B = 100  # [mm] Outer width of the tube (along z-axis)
    H = 100  # [mm] Outer height of the tube (along y-axis)
    t = 6    # [mm] Wall thickness of the tube
    # Number of fibers along each wall direction
    NUM_B = 10  # Number of fibers along the width (z-axis)
    NUM_H = 20  # Number of fibers along the height (y-axis)
    
    # Mass per Length
    AREA = 2 * ((B - 2*t) * t + H * t)  # Correct area calculation
    MASS = DENSITY * AREA  
    
    # Define the rectangular tube section using quad fibers
    ops.section('Fiber', secTag)
    # Top wall (y from H/2 - t to H/2, z from -B/2 to B/2)
    ops.patch('rect', matTag, NUM_B, 10, 
              H/2 - t, -B/2,
              H/2 - t, B/2,
              H/2, B/2,
              H/2, -B/2)
    # Bottom wall (y from -H/2 to -H/2 + t, z from -B/2 to B/2)
    ops.patch('rect', matTag, NUM_B, 10,
              -H/2, -B/2,
              -H/2, B/2,
              -H/2 + t, B/2,
              -H/2 + t, -B/2)
    # Left wall (z from -B/2 to -B/2 + t, y from -H/2 to H/2)
    ops.patch('rect', matTag, NUM_H, 10,
              -H/2, -B/2,
              -H/2, -B/2 + t,
              H/2, -B/2 + t,
              H/2, -B/2)
    # Right wall (z from B/2 - t to B/2, y from -H/2 to H/2)
    ops.patch('rect', matTag, NUM_H, 10,
              -H/2, B/2 - t,
              -H/2, B/2,
              H/2, B/2,
              H/2, B/2 - t)
        
    #-------------------
    # PLOT THE SECTION
    #-------------------
    if PLOT:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 8))

        # Plot outer top wall fibers
        fiber_width = B / NUM_B
        for i in range(NUM_B):
            x_loc = -B / 2 + i * fiber_width
            ax.add_patch(plt.Rectangle((x_loc, H / 2 - t), fiber_width, t, color='grey', alpha=0.7))

        # Plot outer bottom wall fibers
        for i in range(NUM_B):
            x_loc = -B / 2 + i * fiber_width
            ax.add_patch(plt.Rectangle((x_loc, -H / 2), fiber_width, t, color='grey', alpha=0.7))

        # Plot outer left wall fibers
        fiber_height = H / NUM_H
        for i in range(NUM_H):
            y_loc = -H / 2 + i * fiber_height
            ax.add_patch(plt.Rectangle((-B / 2, y_loc), t, fiber_height, color='grey', alpha=0.7))

        # Plot outer right wall fibers
        for i in range(NUM_H):
            y_loc = -H / 2 + i * fiber_height
            ax.add_patch(plt.Rectangle((B / 2 - t, y_loc), t, fiber_height, color='grey', alpha=0.7))

        # Add section outline
        outer_outline = plt.Rectangle((-B / 2, -H / 2), B, H, edgecolor='black', fill=False, linewidth=1.5)
        inner_outline = plt.Rectangle((-B / 2 + t, -H / 2 + t), B - 2 * t, H - 2 * t, edgecolor='black', fill=False, linewidth=1.5)
        ax.add_patch(outer_outline)
        ax.add_patch(inner_outline)

        # Set plot limits
        ax.set_xlim([-B / 2 - 20, B / 2 + 20])
        ax.set_ylim([-H / 2 - 20, H / 2 + 20])
        ax.set_aspect('equal', adjustable='box')

        # Labels and title
        plt.xlabel('Width (B)', fontsize=12)
        plt.ylabel('Height (H)', fontsize=12)
        plt.title('Rectangular Tube Section with Fibers', fontsize=14)

        # Show plot
        plt.grid(True)
        plt.show()
    
    return H, MASS  # Return Section Height and Mass per unit length
#----------------------------------------------------------------------------


 






 



