import numpy as np
import openseespy.opensees as ops
#----------------------------------------
# Concrete Fiber Section without Rebars 
#----------------------------------------
def R_RECTANGULAR_CONCRETE_SECTION(secTag, B, H, NUM_B, NUM_H, matTag, DENSITY=2500/1e9):
    #matTag = 1    # Material tag for the concrete
    
    # Mass per Length
    AREA = B * H
    MASS = DENSITY * AREA
    
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

    return H, MASS

#--------------------------------------------------------------------------
    
#-------------------------------------------------------   
# Concrete Fiber Beam Section with Rebars (NUM_LAYERS)
#-------------------------------------------------------
def R_RECTANGULAR_CONCRETE_SECTION_REBAR_B(secTag, B, H, COVER, RD, NUM_B, NUM_H, NUM_LAYERS, matTag_concrete, matTag_steel, PLOT=True, DENSITY=2500/1e9):
    
    # Mass per Length
    AREA = B * H
    MASS = DENSITY * AREA
    
    # Define the concrete section using fibers
    ops.section('Fiber', secTag)

    # Create concrete fibers for the rectangular section
    for i in range(NUM_B):
        for j in range(NUM_H):
            # Calculate the x and y location of the fiber centroid
            x_loc = -B / 2 + (i + 0.5) * (B / NUM_B)
            y_loc = -H / 2 + (j + 0.5) * (H / NUM_H)
            # Define the fiber area
            fiber_area = (B / NUM_B) * (H / NUM_H)
            # Add the concrete fiber to the section
            ops.fiber(x_loc, y_loc, fiber_area, matTag_concrete)

    # Add steel rebars
    rebar_area = np.pi * (RD / 2) ** 2  # Calculate rebar cross-sectional area

    # Calculate rebar layer positions, starting from COVER
    effective_height = H - 2 * COVER  # Total height available for rebar layers
    layer_spacing = effective_height / (NUM_LAYERS - 1) if NUM_LAYERS > 1 else 0
    rebar_positions = [-H / 2 + COVER + i * layer_spacing for i in range(NUM_LAYERS)]

    for y_loc in rebar_positions:
        x_loc_1 = -B / 2 + COVER  # First rebar in the layer
        x_loc_2 = B / 2 - COVER   # Second rebar in the layer

        # Add first rebar in the layer
        ops.fiber(x_loc_1, y_loc, rebar_area, matTag_steel)
        # Add second rebar in the layer
        ops.fiber(x_loc_2, y_loc, rebar_area, matTag_steel)

    #-------------------
    # PLOT THE SECTION
    #-------------------
    if PLOT == True:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 8))

        # Plot the concrete fibers
        fiber_width = B / NUM_B
        fiber_height = H / NUM_H

        for i in range(NUM_B):
            for j in range(NUM_H):
                x_loc = -B / 2 + i * fiber_width
                y_loc = -H / 2 + j * fiber_height
                rect = plt.Rectangle((x_loc, y_loc), fiber_width, fiber_height, color='lightgrey', alpha=0.7)
                ax.add_patch(rect)

        # Plot the rebars
        rebar_radius = RD / 2

        for y_loc in rebar_positions:
            x_loc_1 = -B / 2 + COVER  # First rebar in the layer
            x_loc_2 = B / 2 - COVER   # Second rebar in the layer

            # Plot first rebar
            circle1 = plt.Circle((x_loc_1, y_loc), rebar_radius, color='red')
            ax.add_patch(circle1)

            # Plot second rebar
            circle2 = plt.Circle((x_loc_2, y_loc), rebar_radius, color='red')
            ax.add_patch(circle2)

        # Add section outline
        section_outline = plt.Rectangle((-B / 2, -H / 2), B, H, edgecolor='black', fill=False, linewidth=1.5)
        ax.add_patch(section_outline)

        # Set plot limits
        ax.set_xlim([-B / 2 - RD, B / 2 + RD])
        ax.set_ylim([-H / 2 - RD, H / 2 + RD])
        ax.set_aspect('equal', adjustable='box')

        # Labels and title
        plt.xlabel('Width (B)', fontsize=12)
        plt.ylabel('Height (H)', fontsize=12)
        plt.title('Rectangular Concrete Section with Steel Rebars', fontsize=14)

        # Show plot
        plt.grid(True)
        plt.show()

    return H, MASS # Return Section Height

#--------------------------------------------------------------------------
 
#---------------------------------------------------------   
# Concrete Fiber Column Section with Rebars (NUM_LAYERS)
#--------------------------------------------------------- 
def R_RECTANGULAR_CONCRETE_SECTION_REBAR_C(secTag, B, H, COVER, RD, NUM_B, NUM_H, NUM_LAYERS, NUM_FIRST, matTag_concrete, matTag_steel, PLOT=True, DENSITY=2500/1e9):
    
    # Mass per Length
    AREA = B * H
    MASS = DENSITY * AREA
    
    # Define the concrete section using fibers
    ops.section('Fiber', secTag)

    # Create concrete fibers for the rectangular section
    for i in range(NUM_B):
        for j in range(NUM_H):
            # Calculate the x and y location of the fiber centroid
            x_loc = -B / 2 + (i + 0.5) * (B / NUM_B)
            y_loc = -H / 2 + (j + 0.5) * (H / NUM_H)
            # Define the fiber area
            fiber_area = (B / NUM_B) * (H / NUM_H)
            # Add the concrete fiber to the section
            ops.fiber(x_loc, y_loc, fiber_area, matTag_concrete)

    # Add steel rebars
    rebar_area = np.pi * (RD / 2) ** 2  # Calculate rebar cross-sectional area

    # Calculate rebar layer positions, starting from COVER
    effective_height = H - 2 * COVER  # Total height available for rebar layers
    layer_spacing = effective_height / (NUM_LAYERS - 1) if NUM_LAYERS > 1 else 0
    rebar_positions = [-H / 2 + COVER + i * layer_spacing for i in range(NUM_LAYERS)]

    for idx, y_loc in enumerate(rebar_positions):
        if idx == 0 or idx == (len(rebar_positions) - 1):
            # First or last layer: 5 rebars spaced across the width
            start_x = -B / 2 + COVER
            end_x = B / 2 - COVER
            spacing_x = (end_x - start_x) / (NUM_FIRST - 1) if NUM_FIRST > 1 else 0
            for k in range(NUM_FIRST):
                x_loc = start_x + k * spacing_x
                ops.fiber(x_loc, y_loc, rebar_area, matTag_steel)
        else:
            # Middle layers: 2 rebars at edges
            x_loc_1 = -B / 2 + COVER
            x_loc_2 = B / 2 - COVER
            ops.fiber(x_loc_1, y_loc, rebar_area, matTag_steel)
            ops.fiber(x_loc_2, y_loc, rebar_area, matTag_steel)

    #-------------------
    # PLOT THE SECTION
    #-------------------
    if PLOT:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 8))

        # Plot the concrete fibers
        fiber_width = B / NUM_B
        fiber_height = H / NUM_H

        for i in range(NUM_B):
            for j in range(NUM_H):
                x_loc = -B / 2 + i * fiber_width
                y_loc = -H / 2 + j * fiber_height
                rect = plt.Rectangle((x_loc, y_loc), fiber_width, fiber_height, color='lightgrey', alpha=0.7)
                ax.add_patch(rect)

        # Plot the rebars
        rebar_radius = RD / 2

        for idx, y_loc in enumerate(rebar_positions):
            if idx == 0 or idx == (len(rebar_positions) - 1):
                # First or last layer: 5 rebars
                start_x = -B / 2 + COVER
                end_x = B / 2 - COVER
                spacing_x = (end_x - start_x) / (NUM_FIRST - 1) if NUM_FIRST > 1 else 0
                for k in range(NUM_FIRST):
                    x = start_x + k * spacing_x
                    circle = plt.Circle((x, y_loc), rebar_radius, color='red')
                    ax.add_patch(circle)
            else:
                # Middle layers: 2 rebars
                x1 = -B / 2 + COVER
                x2 = B / 2 - COVER
                circle1 = plt.Circle((x1, y_loc), rebar_radius, color='red')
                ax.add_patch(circle1)
                circle2 = plt.Circle((x2, y_loc), rebar_radius, color='red')
                ax.add_patch(circle2)

        # Add section outline
        section_outline = plt.Rectangle((-B / 2, -H / 2), B, H, edgecolor='black', fill=False, linewidth=1.5)
        ax.add_patch(section_outline)

        # Set plot limits
        ax.set_xlim([-B / 2 - RD, B / 2 + RD])
        ax.set_ylim([-H / 2 - RD, H / 2 + RD])
        ax.set_aspect('equal', adjustable='box')

        # Labels and title
        plt.xlabel('Width (B)', fontsize=12)
        plt.ylabel('Height (H)', fontsize=12)
        plt.title('Rectangular Concrete Section with Steel Rebars', fontsize=14)

        # Show plot
        plt.grid(True)
        plt.show()

    return H, MASS  # Return Section Height as per original function    
    
#--------------------------------------------------------------------------

#---------------------------------------------------  
# Concrete Fiber T Section with Rebars (NUM_LAYERS)
#---------------------------------------------------

def T_SECTION(secTag, Bf, Hf, Bw, Hw, cover, rebarDiam, num_rebars_flange, num_rebars_web, concMatTag, rebarMatTag, PLOT=True, concDensity=2400/1e9, rebarDensity=7850/1e9):
    """
    Defines a T-section with reinforcement bars (rebars) and optionally plots the section.

    Parameters:
        secTag (int): Tag for the section.
        Bf (float): Width of the flange (mm).
        Hf (float): Height of the flange (mm).
        Bw (float): Width of the web (mm).
        Hw (float): Height of the web (mm).
        rebarDiam (float): Diameter of the rebars (mm).
        num_rebars_flange (int): Number of rebars in the flange.
        num_rebars_web (int): Number of rebars in the web.
        concMatTag (int): Material tag for concrete.
        rebarMatTag (int): Material tag for rebars.
        PLOT (bool): If True, plots the section. Default is True.
        concDensity (float): Density of concrete (kg/mm³). Default is 2400 kg/m³ converted to kg/mm³.
        rebarDensity (float): Density of steel rebars (kg/mm³). Default is 7850 kg/m³ converted to kg/mm³.

    Returns:
        total_height (float): Total height of the T-section (mm).
        mass_per_unit_length (float): Mass per unit length of the section (kg/mm).
    """
    import numpy as np

    # Rebar parameters
    #cover = 40  # Cover for rebars (mm) - can be adjusted as needed
    rebar_spacing_flange = (Bf - 2 * cover) / (num_rebars_flange - 1)  # Spacing in flange
    rebar_spacing_web = (Hw - 2 * cover) / (num_rebars_web - 1)  # Spacing in web

    # Rebar area calculation
    A_rebar = np.pi * (rebarDiam / 2) ** 2  # Area of a single rebar

    # Number of fibers for concrete discretization
    NUM_X_FLANGE = 10  # Number of fibers along the flange width
    NUM_Y_FLANGE = 5   # Number of fibers along the flange height
    NUM_X_WEB = 5      # Number of fibers along the web width
    NUM_Y_WEB = 10     # Number of fibers along the web height

    # Create the fiber section
    ops.section('Fiber', secTag)

    # Add concrete fibers for flange
    dx_flange = Bf / NUM_X_FLANGE
    dy_flange = Hf / NUM_Y_FLANGE
    for i in range(NUM_X_FLANGE):
        x = -Bf / 2 + (i + 0.5) * dx_flange
        for j in range(NUM_Y_FLANGE):
            y = Hw / 2 + (j + 0.5) * dy_flange
            fiber_area = dx_flange * dy_flange
            ops.fiber(x, y, fiber_area, concMatTag)

    # Add concrete fibers for web
    dx_web = Bw / NUM_X_WEB
    dy_web = Hw / NUM_Y_WEB
    for i in range(NUM_X_WEB):
        x = -Bw / 2 + (i + 0.5) * dx_web
        for j in range(NUM_Y_WEB):
            y = -Hw / 2 + (j + 0.5) * dy_web
            fiber_area = dx_web * dy_web
            ops.fiber(x, y, fiber_area, concMatTag)

    # Add rebar fibers (flange and web)
    # Flange rebars
    for i in range(num_rebars_flange):
        x = -Bf / 2 + cover + i * rebar_spacing_flange
        y = Hw / 2 + Hf - cover
        ops.fiber(x, y, A_rebar, rebarMatTag)

    # Web rebars
    for j in range(num_rebars_web):
        y = -Hw / 2 + cover + j * rebar_spacing_web
        x = 0  # Centered in the web
        ops.fiber(x, y, A_rebar, rebarMatTag)

    # Calculate mass per unit length
    flange_area = Bf * Hf
    web_area = Bw * Hw
    rebar_area_total = (num_rebars_flange + num_rebars_web) * A_rebar
    mass_per_unit_length = ((flange_area + web_area) * concDensity) + (rebar_area_total * rebarDensity)

    # Plotting
    if PLOT:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 8))

        # Plot concrete fibers for flange
        for i in range(NUM_X_FLANGE):
            x = -Bf / 2 + (i + 0.5) * dx_flange
            for j in range(NUM_Y_FLANGE):
                y = Hw / 2 + (j + 0.5) * dy_flange
                ax.add_patch(plt.Rectangle((x - dx_flange / 2, y - dy_flange / 2), dx_flange, dy_flange, color='grey', alpha=0.5))

        # Plot concrete fibers for web
        for i in range(NUM_X_WEB):
            x = -Bw / 2 + (i + 0.5) * dx_web
            for j in range(NUM_Y_WEB):
                y = -Hw / 2 + (j + 0.5) * dy_web
                ax.add_patch(plt.Rectangle((x - dx_web / 2, y - dy_web / 2), dx_web, dy_web, color='grey', alpha=0.5))

        # Plot rebar fibers
        rebar_radius = max(2, rebarDiam / 2)  # Ensure visibility
        # Flange rebars
        for i in range(num_rebars_flange):
            x = -Bf / 2 + cover + i * rebar_spacing_flange
            y = Hw / 2 + Hf - cover
            ax.add_patch(plt.Circle((x, y), rebar_radius, color='red', alpha=0.8))
        # Web rebars
        for j in range(num_rebars_web):
            y = -Hw / 2 + cover + j * rebar_spacing_web
            x = 0  # Centered in the web
            ax.add_patch(plt.Circle((x, y), rebar_radius, color='red', alpha=0.8))

        # Plot section boundaries
        ax.add_patch(plt.Rectangle((-Bf / 2, Hw / 2), Bf, Hf, color='black', fill=False, linewidth=1.5))  # Flange
        ax.add_patch(plt.Rectangle((-Bw / 2, -Hw / 2), Bw, Hw, color='black', fill=False, linewidth=1.5))  # Web

        ax.set_aspect('equal')
        ax.set_xlim(-Bf / 2 - 20, Bf / 2 + 20)
        ax.set_ylim(-Hw / 2 - 20, (Hw / 2 + Hf) + 20)
        plt.xlabel('X (mm)')
        plt.ylabel('Y (mm)')
        plt.title('T-Section with Rebar')
        plt.grid(True)
        plt.show()

    return Hw + Hf, mass_per_unit_length  # Return total height and mass per unit length
    
#--------------------------------------------------------------------------  
#---------------------------------    
# CIRCULAR TUBE CONCRETE SECTION  
#--------------------------------- 

def C_TUBE_SECTION(secTag, D, t, cover, rebarDiam, num_rebars_per_ring, concMatTag, rebarMatTag, PLOT=True, concDensity=2400/1e9, rebarDensity=7850/1e9):
    # Define geometric properties of the concrete tube
    import numpy as np
    
    r_outer = D / 2          # Outer radius (mm)
    r_inner = r_outer - t    # Inner radius (mm)

    # Rebar parameters
    cover_outer = cover#t / 4      # Cover for outer rebar ring (mm)
    cover_inner = cover#t / 4      # Cover for inner rebar ring (mm)
    r_rebar_outer = r_outer - cover_outer  # Radius of outer rebar ring
    r_rebar_inner = r_inner + cover_inner  # Radius of inner rebar ring

    # Rebar area calculation
    A_rebar = np.pi * (rebarDiam / 2) ** 2  # Area of a single rebar

    # Number of fibers for concrete discretization
    NUM_R_CONC = 10  # Radial fibers for concrete
    NUM_C_CONC = 20  # Circumferential fibers for concrete

    # Create the fiber section
    ops.section('Fiber', secTag)

    # Add concrete fibers
    for i in range(NUM_R_CONC):
        dr = t / NUM_R_CONC  # Radial increment
        r_loc = r_inner + (i + 0.5) * dr  # Radial position of the fiber
        for j in range(NUM_C_CONC):
            theta = j * (2 * np.pi / NUM_C_CONC)  # Angle theta
            x = r_loc * np.cos(theta)
            y = r_loc * np.sin(theta)
            fiber_area = (2 * np.pi * r_loc / NUM_C_CONC) * dr  # Area of the concrete fiber
            ops.fiber(x, y, fiber_area, concMatTag)

    # Add rebar fibers (two rings, 20 rebars each)
    for j in range(num_rebars_per_ring):
        theta = j * (2 * np.pi / num_rebars_per_ring)
        
        # Outer rebar ring
        x_outer = r_rebar_outer * np.cos(theta)
        y_outer = r_rebar_outer * np.sin(theta)
        ops.fiber(x_outer, y_outer, A_rebar, rebarMatTag)
        
        # Inner rebar ring
        x_inner = r_rebar_inner * np.cos(theta)
        y_inner = r_rebar_inner * np.sin(theta)
        ops.fiber(x_inner, y_inner, A_rebar, rebarMatTag)

    # Calculate mass per unit length
    concrete_area = np.pi * (r_outer**2 - r_inner**2)
    rebar_area_total = 2 * num_rebars_per_ring * A_rebar  # Total rebar area (two rings)
    MASS = (concrete_area * concDensity) + (rebar_area_total * rebarDensity)

    # Plotting the section
    if PLOT:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 8))

        # Plot concrete fibers
        for i in range(NUM_R_CONC):
            dr = t / NUM_R_CONC
            r_loc = r_inner + (i + 0.5) * dr
            for j in range(NUM_C_CONC):
                theta = j * (2 * np.pi / NUM_C_CONC)
                x = r_loc * np.cos(theta)
                y = r_loc * np.sin(theta)
                ax.add_patch(plt.Circle((x, y), 2, color='grey', alpha=0.5))

        # Plot rebar fibers
        rebar_radius = max(2, rebarDiam / 2)  # Ensure visibility
        for j in range(num_rebars_per_ring):
            theta = j * (2 * np.pi / num_rebars_per_ring)
            
            # Outer rebar
            x_outer = r_rebar_outer * np.cos(theta)
            y_outer = r_rebar_outer * np.sin(theta)
            ax.add_patch(plt.Circle((x_outer, y_outer), rebar_radius, color='red', alpha=0.8))
            
            # Inner rebar
            x_inner = r_rebar_inner * np.cos(theta)
            y_inner = r_rebar_inner * np.sin(theta)
            ax.add_patch(plt.Circle((x_inner, y_inner), rebar_radius, color='red', alpha=0.8))

        # Plot concrete boundaries
        outer_circle = plt.Circle((0, 0), r_outer, color='black', fill=False, linewidth=1.5)
        inner_circle = plt.Circle((0, 0), r_inner, color='black', fill=False, linewidth=1.5)
        ax.add_patch(outer_circle)
        ax.add_patch(inner_circle)

        ax.set_aspect('equal')
        ax.set_xlim(-r_outer - 20, r_outer + 20)
        ax.set_ylim(-r_outer - 20, r_outer + 20)
        plt.xlabel('X (mm)')
        plt.ylabel('Y (mm)')
        plt.title('Tubular Concrete Section with Rebar Rings')
        plt.grid(True)
        plt.show()

    return D, MASS  # Return outer diameter and mass per unit length
    
#----------------------------------------------------------------------------  
 
#------------------------------------    
# RECTANGULAR TUBE CONCRETE SECTION  
#------------------------------------  

def RECT_BOX_SECTION(secTag, B, H, t, cover, rebarDiam,num_rebars_width, num_rebars_depth, concMatTag, rebarMatTag,  PLOT=True, concDensity=2400/1e9, rebarDensity=7850/1e9):
    # Define geometric properties of the rectangular box section
    import numpy as np
    
    #num_rebars_width = 10
    #num_rebars_depth = 5
    
    # Outer and inner dimensions
    B_outer = B  # Outer width (mm)
    H_outer = H  # Outer height (mm)
    B_inner = B_outer - 2 * t  # Inner width (mm)
    H_inner = H_outer - 2 * t  # Inner height (mm)

    # Rebar parameters
    #cover = t / 4  # Cover for rebars (mm)
    rebar_spacing_width = (B_outer - 2 * cover) / (num_rebars_width - 1)  # Spacing between rebars along width
    rebar_spacing_depth = (H_outer - 2 * cover) / (num_rebars_depth - 1)  # Spacing between rebars along depth

    # Rebar area calculation
    A_rebar = np.pi * (rebarDiam / 2) ** 2  # Area of a single rebar

    # Number of fibers for concrete discretization
    NUM_X_CONC = 10  # Number of fibers along the width
    NUM_Y_CONC = 10  # Number of fibers along the height

    # Create the fiber section
    ops.section('Fiber', secTag)

    # Add concrete fibers
    dx = B_outer / NUM_X_CONC  # Width increment
    dy = H_outer / NUM_Y_CONC  # Height increment
    for i in range(NUM_X_CONC):
        x = -B_outer / 2 + (i + 0.5) * dx  # X position of the fiber
        for j in range(NUM_Y_CONC):
            y = -H_outer / 2 + (j + 0.5) * dy  # Y position of the fiber
            # Check if the fiber is within the inner box (to exclude the hollow part)
            if not (-B_inner / 2 < x < B_inner / 2 and -H_inner / 2 < y < H_inner / 2):
                fiber_area = dx * dy  # Area of the concrete fiber
                ops.fiber(x, y, fiber_area, concMatTag)

    # Add rebar fibers (along the perimeter of the box)
    # Top and bottom sides (width direction)
    for i in range(num_rebars_width):
        x = -B_outer / 2 + cover + i * rebar_spacing_width
        # Top side
        y_top = H_outer / 2 - cover
        ops.fiber(x, y_top, A_rebar, rebarMatTag)
        # Bottom side
        y_bottom = -H_outer / 2 + cover
        ops.fiber(x, y_bottom, A_rebar, rebarMatTag)

    # Left and right sides (depth direction)
    for j in range(1, num_rebars_depth-1):
        y = -H_outer / 2 + cover + j * rebar_spacing_depth
        # Right side
        x_right = B_outer / 2 - cover
        ops.fiber(x_right, y, A_rebar, rebarMatTag)
        # Left side
        x_left = -B_outer / 2 + cover
        ops.fiber(x_left, y, A_rebar, rebarMatTag)

    # Calculate mass per unit length
    concrete_area = B_outer * H_outer - B_inner * H_inner  # Area of concrete
    rebar_area_total = 2 * (num_rebars_width + num_rebars_depth) * A_rebar  # Total rebar area (4 sides)
    MASS = (concrete_area * concDensity) + (rebar_area_total * rebarDensity)

    # Plotting the section
    if PLOT:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 8))

        # Plot concrete fibers
        for i in range(NUM_X_CONC):
            x = -B_outer / 2 + (i + 0.5) * dx
            for j in range(NUM_Y_CONC):
                y = -H_outer / 2 + (j + 0.5) * dy
                if not (-B_inner / 2 < x < B_inner / 2 and -H_inner / 2 < y < H_inner / 2):
                    ax.add_patch(plt.Rectangle((x - dx / 2, y - dy / 2), dx, dy, color='grey', alpha=0.5))

        # Plot rebar fibers
        rebar_radius = max(2, rebarDiam / 2)  # Ensure visibility
        # Top and bottom sides (width direction)
        for i in range(num_rebars_width):
            x = -B_outer / 2 + cover + i * rebar_spacing_width
            # Top side
            y_top = H_outer / 2 - cover
            ax.add_patch(plt.Circle((x, y_top), rebar_radius, color='red', alpha=0.8))
            # Bottom side
            y_bottom = -H_outer / 2 + cover
            ax.add_patch(plt.Circle((x, y_bottom), rebar_radius, color='red', alpha=0.8))
        # Left and right sides (depth direction)
        for j in range(1, num_rebars_depth-1):
            y = -H_outer / 2 + cover + j * rebar_spacing_depth
            # Right side
            x_right = B_outer / 2 - cover
            ax.add_patch(plt.Circle((x_right, y), rebar_radius, color='red', alpha=0.8))
            # Left side
            x_left = -B_outer / 2 + cover
            ax.add_patch(plt.Circle((x_left, y), rebar_radius, color='red', alpha=0.8))

        # Plot concrete boundaries
        outer_rect = plt.Rectangle((-B_outer / 2, -H_outer / 2), B_outer, H_outer, color='black', fill=False, linewidth=1.5)
        inner_rect = plt.Rectangle((-B_inner / 2, -H_inner / 2), B_inner, H_inner, color='black', fill=False, linewidth=1.5)
        ax.add_patch(outer_rect)
        ax.add_patch(inner_rect)

        ax.set_aspect('equal')
        ax.set_xlim(-B_outer / 2 - 20, B_outer / 2 + 20)
        ax.set_ylim(-H_outer / 2 - 20, H_outer / 2 + 20)
        plt.xlabel('X (mm)')
        plt.ylabel('Y (mm)')
        plt.title('Rectangular Box Section with Rebar')
        plt.grid(True)
        plt.show()

    return H, MASS  # Return outer dimensions and mass per unit length

#----------------------------------------------------------------------------

#--------------------    
# I CONCRETE SECTION  
#--------------------

def I_SECTION(secTag, Bf, Hf, Bw, Hw, cover, rebarDiam, num_rebars_flange, num_rebars_web, concMatTag, rebarMatTag, PLOT=True, concDensity=2400/1e9, rebarDensity=7850/1e9):
    """
    Defines an I-section with reinforcement bars (rebars) and optionally plots the section.

    Parameters:
        secTag (int): Tag for the section.
        Bf (float): Width of the flange (mm).
        Hf (float): Height of the flange (mm).
        Bw (float): Width of the web (mm).
        Hw (float): Height of the web (mm).
        rebarDiam (float): Diameter of the rebars (mm).
        num_rebars_flange (int): Number of rebars in each flange.
        num_rebars_web (int): Number of rebars in the web.
        concMatTag (int): Material tag for concrete.
        rebarMatTag (int): Material tag for rebars.
        PLOT (bool): If True, plots the section. Default is True.
        concDensity (float): Density of concrete (kg/mm³). Default is 2400 kg/m³ converted to kg/mm³.
        rebarDensity (float): Density of steel rebars (kg/mm³). Default is 7850 kg/m³ converted to kg/mm³.

    Returns:
        total_height (float): Total height of the I-section (mm).
        mass_per_unit_length (float): Mass per unit length of the section (kg/mm).
    """
    import numpy as np

    # Rebar parameters
    #cover = 40  # Cover for rebars (mm) - can be adjusted as needed
    rebar_spacing_flange = (Bf - 2 * cover) / (num_rebars_flange - 1)  # Spacing in flange
    rebar_spacing_web = (Hw - 2 * cover) / (num_rebars_web - 1)  # Spacing in web

    # Rebar area calculation
    A_rebar = np.pi * (rebarDiam / 2) ** 2  # Area of a single rebar

    # Number of fibers for concrete discretization
    NUM_X_FLANGE = 10  # Number of fibers along the flange width
    NUM_Y_FLANGE = 5   # Number of fibers along the flange height
    NUM_X_WEB = 5      # Number of fibers along the web width
    NUM_Y_WEB = 10     # Number of fibers along the web height

    # Create the fiber section
    ops.section('Fiber', secTag)

    # Add concrete fibers for top flange
    dx_flange = Bf / NUM_X_FLANGE
    dy_flange = Hf / NUM_Y_FLANGE
    for i in range(NUM_X_FLANGE):
        x = -Bf / 2 + (i + 0.5) * dx_flange
        for j in range(NUM_Y_FLANGE):
            y = Hw / 2 + Hf - (j + 0.5) * dy_flange
            fiber_area = dx_flange * dy_flange
            ops.fiber(x, y, fiber_area, concMatTag)

    # Add concrete fibers for bottom flange
    for i in range(NUM_X_FLANGE):
        x = -Bf / 2 + (i + 0.5) * dx_flange
        for j in range(NUM_Y_FLANGE):
            y = -Hw / 2 - Hf + (j + 0.5) * dy_flange
            fiber_area = dx_flange * dy_flange
            ops.fiber(x, y, fiber_area, concMatTag)

    # Add concrete fibers for web
    dx_web = Bw / NUM_X_WEB
    dy_web = Hw / NUM_Y_WEB
    for i in range(NUM_X_WEB):
        x = -Bw / 2 + (i + 0.5) * dx_web
        for j in range(NUM_Y_WEB):
            y = -Hw / 2 + (j + 0.5) * dy_web
            fiber_area = dx_web * dy_web
            ops.fiber(x, y, fiber_area, concMatTag)

    # Add rebar fibers (flanges and web)
    # Top flange rebars
    for i in range(num_rebars_flange):
        x = -Bf / 2 + cover + i * rebar_spacing_flange
        y = Hw / 2 + Hf - cover
        ops.fiber(x, y, A_rebar, rebarMatTag)

    # Bottom flange rebars
    for i in range(num_rebars_flange):
        x = -Bf / 2 + cover + i * rebar_spacing_flange
        y = -Hw / 2 - Hf + cover
        ops.fiber(x, y, A_rebar, rebarMatTag)

    # Web rebars
    for j in range(num_rebars_web):
        y = -Hw / 2 + cover + j * rebar_spacing_web
        x = 0  # Centered in the web
        ops.fiber(x, y, A_rebar, rebarMatTag)

    # Calculate mass per unit length
    flange_area = 2 * Bf * Hf  # Two flanges
    web_area = Bw * Hw
    rebar_area_total = (2 * num_rebars_flange + num_rebars_web) * A_rebar
    mass_per_unit_length = ((flange_area + web_area) * concDensity) + (rebar_area_total * rebarDensity)

    # Plotting
    if PLOT:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 8))

        # Plot concrete fibers for top flange
        for i in range(NUM_X_FLANGE):
            x = -Bf / 2 + (i + 0.5) * dx_flange
            for j in range(NUM_Y_FLANGE):
                y = Hw / 2 + Hf - (j + 0.5) * dy_flange
                ax.add_patch(plt.Rectangle((x - dx_flange / 2, y - dy_flange / 2), dx_flange, dy_flange, color='grey', alpha=0.5))

        # Plot concrete fibers for bottom flange
        for i in range(NUM_X_FLANGE):
            x = -Bf / 2 + (i + 0.5) * dx_flange
            for j in range(NUM_Y_FLANGE):
                y = -Hw / 2 - Hf + (j + 0.5) * dy_flange
                ax.add_patch(plt.Rectangle((x - dx_flange / 2, y - dy_flange / 2), dx_flange, dy_flange, color='grey', alpha=0.5))

        # Plot concrete fibers for web
        for i in range(NUM_X_WEB):
            x = -Bw / 2 + (i + 0.5) * dx_web
            for j in range(NUM_Y_WEB):
                y = -Hw / 2 + (j + 0.5) * dy_web
                ax.add_patch(plt.Rectangle((x - dx_web / 2, y - dy_web / 2), dx_web, dy_web, color='grey', alpha=0.5))

        # Plot rebar fibers
        rebar_radius = max(2, rebarDiam / 2)  # Ensure visibility
        # Top flange rebars
        for i in range(num_rebars_flange):
            x = -Bf / 2 + cover + i * rebar_spacing_flange
            y = Hw / 2 + Hf - cover
            ax.add_patch(plt.Circle((x, y), rebar_radius, color='red', alpha=0.8))
        # Bottom flange rebars
        for i in range(num_rebars_flange):
            x = -Bf / 2 + cover + i * rebar_spacing_flange
            y = -Hw / 2 - Hf + cover
            ax.add_patch(plt.Circle((x, y), rebar_radius, color='red', alpha=0.8))
        # Web rebars
        for j in range(num_rebars_web):
            y = -Hw / 2 + cover + j * rebar_spacing_web
            x = 0  # Centered in the web
            ax.add_patch(plt.Circle((x, y), rebar_radius, color='red', alpha=0.8))

        # Plot section boundaries
        ax.add_patch(plt.Rectangle((-Bf / 2, Hw / 2), Bf, Hf, color='black', fill=False, linewidth=1.5))  # Top flange
        ax.add_patch(plt.Rectangle((-Bf / 2, -Hw / 2 - Hf), Bf, Hf, color='black', fill=False, linewidth=1.5))  # Bottom flange
        ax.add_patch(plt.Rectangle((-Bw / 2, -Hw / 2), Bw, Hw, color='black', fill=False, linewidth=1.5))  # Web

        ax.set_aspect('equal')
        ax.set_xlim(-Bf / 2 - 20, Bf / 2 + 20)
        ax.set_ylim(-Hw / 2 - Hf - 20, Hw / 2 + Hf + 20)
        plt.xlabel('X (mm)')
        plt.ylabel('Y (mm)')
        plt.title('I-Section with Rebar')
        plt.grid(True)
        plt.show()

    return Hw + 2 * Hf, mass_per_unit_length  # Return total height and mass per unit length

#----------------------------------------------------------------------------
#----------------------------------------    
# CONCRETE FIBER SECTION FROM EXCEL DATA
#----------------------------------------

def CONCRETE_SECTION_FROM_EXCEL(secTag, excel_file, B, H, FIB_NUM_X, FIB_NUM_Y, matTag_concrete, matTag_steel, PLOT=True, CONCRETE_DENSITY=2500/1e9, STEEL_DENSITY=7850/1e9):
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    #import opensees.openseespy as ops
    # Read data from Excel
    df = pd.read_excel(excel_file)
    # Compute H and B based on fiber centers
    y_centers = df['Y-center']
    #H = y_centers.max() - y_centers.min()# + (y_centers[0] + y_centers[1]) * 0.5
    x_centers = df['X-center']
    #B = x_centers.max() - x_centers.min()# + (x_centers[0] + x_centers[1]) * 0.5

    # Calculate total mass
    concrete_mask = df['MatTag'] == 1 # 1 of concrete 2 for steel
    steel_mask = df['MatTag'] == 2 # 1 of concrete 2 for steel
    
    concrete_area = df.loc[concrete_mask, 'Fiber_Area'].sum()
    steel_area = df.loc[steel_mask, 'Fiber_Area'].sum()
    mass = (concrete_area * CONCRETE_DENSITY) + (steel_area * STEEL_DENSITY)
    
    # Create the fiber section
    ops.section('Fiber', secTag)
    
    # Add fibers to the section
    for idx, row in df.iterrows():
        mat_tag = row['MatTag']
        x = row['X-center']
        y = row['Y-center']
        area = row['Fiber_Area']
        
        if mat_tag == 1:
            ops.fiber(x, y, area, matTag_concrete)
        elif mat_tag == 2:
            ops.fiber(x, y, area, matTag_steel)
        else:
            print(f"Warning: Unknown material tag {mat_tag} at row {idx}")
    
    # Plotting the section
    if PLOT:
        fig, ax = plt.subplots(figsize=(8, 8))
        min_x, max_x = np.inf, -np.inf
        min_y, max_y = np.inf, -np.inf
        
        for idx, row in df.iterrows():
            mat_tag = row['MatTag']
            x = row['X-center']
            y = row['Y-center']
            area = row['Fiber_Area']
            
            if mat_tag == 1:  # Concrete fiber
                #side = np.sqrt(area)
                sideX = B / FIB_NUM_X
                left = x - sideX/2
                right = x + sideX/2
                sideY = H / FIB_NUM_Y
                bottom = y - sideY/2
                top = y + sideY/2
                #rect = plt.Rectangle((left, bottom), side, side, color='grey', alpha=0.7)
                rect = plt.Rectangle((left, bottom), sideX, sideY, color='grey', alpha=0.7)
                ax.add_patch(rect)
            elif mat_tag == 2:  # Steel rebar
                radius = np.sqrt(area / np.pi)
                left = x - radius
                right = x + radius
                bottom = y - radius
                top = y + radius
                circle = plt.Circle((x, y), radius, color='red')
                ax.add_patch(circle)
            
            # Update bounds
            min_x = min(min_x, left, right)
            max_x = max(max_x, left, right)
            min_y = min(min_y, bottom, top)
            max_y = max(max_y, bottom, top)
            #print(min_x, max_x, min_y, max_y)
        
        # Add section outline
        outline_width = max_x - min_x
        outline_height = max_y - min_y
        #outline = plt.Rectangle((min_x, min_y), outline_width, outline_height, edgecolor='black', fill=False, linewidth=1.5)
        #ax.add_patch(outline)
        
        # Adjust plot settings
        ax.set_xlim(min_x - 0.1*outline_width, max_x + 0.1*outline_width)
        ax.set_ylim(min_y - 0.1*outline_height, max_y + 0.1*outline_height)
        ax.set_aspect('equal', adjustable='box')
        plt.xlabel('Width', fontsize=12)
        plt.ylabel('Height', fontsize=12)
        plt.title('Fiber Section from Excel Data', fontsize=14)
        plt.grid(True)
        plt.show()
        #print(B, H)
    
    return H, mass

#----------------------------------------------------------------------------   

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

#-----------------------------------------------------    
# CONFINED CONCRETE FIBER T SECTION WITH QUAD FIBERS
#-----------------------------------------------------

def CONCRETE_CONFINED_T_SECTION_QUAD(id, HSec, BFlange, TFlange, TWeb, coverFlange, coverWeb, coreID, coverID, steelID,
    numBarsTop, barAreaTop, numBarsBot, barAreaBot, numBarsIntTot, barAreaInt,
    PLOT=True, CONCRETE_DENSITY=2500/1e9, STEEL_DENSITY=7850/1e9):
    #import openseespy.opensees as ops
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.patches import Rectangle

    # Fiber discretization parameters
    nfCoreY_web = 10    # Fibers along web height
    nfCoreZ_web = 5     # Fibers along web thickness
    nfCoreY_flange = 5  # Fibers along flange thickness
    nfCoreZ_flange = 10 # Fibers along flange width
    nfCoverY = 2        # Fibers in cover along height
    nfCoverZ = 2        # Fibers in cover along width

    # Geometric calculations
    HWeb = HSec - TFlange  # Web height for T-section (total height minus top flange thickness)

    # Mass calculation
    web_area = HWeb * TWeb
    flange_area = BFlange * TFlange  # Only one flange
    concrete_mass = (web_area + flange_area) * CONCRETE_DENSITY
    steel_mass = (numBarsTop * barAreaTop + numBarsBot * barAreaBot + numBarsIntTot * barAreaInt) * STEEL_DENSITY
    MASS = concrete_mass + steel_mass

    # Define coordinates
    # Web (origin at top of web, y=0)
    web_bottom = -HWeb
    web_top = 0
    web_left = -TWeb / 2
    web_right = TWeb / 2
    # Web core
    web_core_bottom = web_bottom + coverWeb
    web_core_top = web_top - coverWeb
    web_core_left = web_left + coverWeb
    web_core_right = web_right - coverWeb

    # Top flange
    top_flange_bottom = web_top
    top_flange_top = top_flange_bottom + TFlange
    top_flange_left = -BFlange / 2
    top_flange_right = BFlange / 2
    # Top flange core
    top_flange_core_bottom = top_flange_bottom + coverFlange
    top_flange_core_top = top_flange_top - coverFlange
    top_flange_core_left = top_flange_left + coverFlange
    top_flange_core_right = top_flange_right - coverFlange

    # Define fiber section in OpenSees
    ops.section('Fiber', id)

    # Core patches
    # Web core
    ops.patch('quad', coreID, nfCoreZ_web, nfCoreY_web,
             web_core_bottom, web_core_right,
             web_core_bottom, web_core_left,
             web_core_top, web_core_left,
             web_core_top, web_core_right)

    # Top flange core
    ops.patch('quad', coreID, nfCoreZ_flange, nfCoreY_flange,
             top_flange_core_bottom, top_flange_core_right,
             top_flange_core_bottom, top_flange_core_left,
             top_flange_core_top, top_flange_core_left,
             top_flange_core_top, top_flange_core_right)

    # Cover patches
    # Web left cover
    ops.patch('quad', coverID, nfCoverZ, nfCoreY_web,
             web_bottom, web_left,
             web_bottom, web_core_left,
             web_top, web_core_left,
             web_top, web_left)

    # Web right cover
    ops.patch('quad', coverID, nfCoverZ, nfCoreY_web,
             web_bottom, web_right,
             web_bottom, web_core_right,
             web_top, web_core_right,
             web_top, web_right)

    # Web bottom cover
    ops.patch('quad', coverID, nfCoreZ_web, nfCoverY,
             web_bottom, web_right,
             web_bottom, web_left,
             web_core_bottom, web_left,
             web_core_bottom, web_right)

    # Top flange covers
    # Bottom cover
    ops.patch('quad', coverID, nfCoreZ_flange, nfCoverY,
             top_flange_bottom, top_flange_right,
             top_flange_bottom, top_flange_left,
             top_flange_core_bottom, top_flange_left,
             top_flange_core_bottom, top_flange_right)

    # Top cover
    ops.patch('quad', coverID, nfCoreZ_flange, nfCoverY,
             top_flange_core_top, top_flange_right,
             top_flange_core_top, top_flange_left,
             top_flange_top, top_flange_left,
             top_flange_top, top_flange_right)

    # Left cover
    ops.patch('quad', coverID, nfCoverZ, nfCoreY_flange,
             top_flange_bottom, top_flange_left,
             top_flange_bottom, top_flange_core_left,
             top_flange_top, top_flange_core_left,
             top_flange_top, top_flange_left)

    # Right cover
    ops.patch('quad', coverID, nfCoverZ, nfCoreY_flange,
             top_flange_bottom, top_flange_right,
             top_flange_bottom, top_flange_core_right,
             top_flange_top, top_flange_core_right,
             top_flange_top, top_flange_right)

    # Reinforcement layers
    # Top bars in top flange
    y_top_bars = top_flange_top - coverFlange
    z_left_top = top_flange_left + coverFlange
    z_right_top = top_flange_right - coverFlange
    ops.layer('straight', steelID, numBarsTop, barAreaTop,
             y_top_bars, z_left_top, y_top_bars, z_right_top)

    # Bottom bars at bottom of web
    y_bot_bars = web_bottom + coverWeb
    z_left_bot = web_left + coverWeb
    z_right_bot = web_right - coverWeb
    ops.layer('straight', steelID, numBarsBot, barAreaBot,
             y_bot_bars, z_left_bot, y_bot_bars, z_right_bot)

    # Intermediate bars in web (split evenly on both sides)
    numBarsInt = numBarsIntTot // 2
    DD = (HWeb - 2 * coverWeb) / (numBarsInt - 1) if numBarsInt > 1 else 0
    y_start = web_bottom + coverWeb
    y_end = web_top - coverWeb
    z_left_int = web_left + coverWeb
    z_right_int = web_right - coverWeb
    ops.layer('straight', steelID, numBarsInt, barAreaInt,
             y_start, z_left_int, y_end, z_left_int)  # Left side
    ops.layer('straight', steelID, numBarsInt, barAreaInt,
             y_start, z_right_int, y_end, z_right_int)  # Right side

    # Plotting
    if PLOT:
        fig, ax = plt.subplots(figsize=(10, 6))

        # Plot core regions
        ax.add_patch(Rectangle((web_core_left, web_core_bottom),
                              web_core_right - web_core_left, web_core_top - web_core_bottom,
                              facecolor='grey', edgecolor='black', label='Web Core'))
        ax.add_patch(Rectangle((top_flange_core_left, top_flange_core_bottom),
                              top_flange_core_right - top_flange_core_left, top_flange_core_top - top_flange_core_bottom,
                              facecolor='grey', edgecolor='black', label='Flange Core'))

        # Plot cover regions
        # Web covers
        ax.add_patch(Rectangle((web_left, web_bottom), coverWeb, HWeb,
                              facecolor='lightgrey', edgecolor='black', label='Cover'))
        ax.add_patch(Rectangle((web_core_right, web_bottom), coverWeb, HWeb,
                              facecolor='lightgrey', edgecolor='black'))
        ax.add_patch(Rectangle((web_left, web_bottom), TWeb, coverWeb,
                              facecolor='lightgrey', edgecolor='black'))  # Bottom cover
        # Top flange covers
        ax.add_patch(Rectangle((top_flange_left, top_flange_bottom), BFlange, coverFlange,
                              facecolor='lightgrey', edgecolor='black'))  # Bottom
        ax.add_patch(Rectangle((top_flange_left, top_flange_core_top), BFlange, coverFlange,
                              facecolor='lightgrey', edgecolor='black'))  # Top
        ax.add_patch(Rectangle((top_flange_left, top_flange_bottom), coverFlange, TFlange,
                              facecolor='lightgrey', edgecolor='black'))  # Left
        ax.add_patch(Rectangle((top_flange_core_right, top_flange_bottom), coverFlange, TFlange,
                              facecolor='lightgrey', edgecolor='black'))  # Right

        # Plot reinforcement
        z_top = np.linspace(z_left_top, z_right_top, numBarsTop)
        ax.scatter(z_top, [y_top_bars]*numBarsTop, s=50, c='red', marker='o', label='Top Steel')
        z_bot = np.linspace(z_left_bot, z_right_bot, numBarsBot)
        ax.scatter(z_bot, [y_bot_bars]*numBarsBot, s=50, c='red', marker='o', label='Bottom Steel')
        y_int = np.linspace(y_start, y_end, numBarsInt)
        ax.scatter([z_left_int]*numBarsInt, y_int, s=30, c='blue', marker='s', label='Skin Steel')
        ax.scatter([z_right_int]*numBarsInt, y_int, s=30, c='blue', marker='s')

        # Plot formatting
        ax.set_aspect('equal')
        ax.set_xlabel('Width (Z)')
        ax.set_ylabel('Height (Y)')
        ax.set_title(f'Confined Concrete T-Section (ID: {id})')
        ax.grid(True, alpha=0.3)
        ax.axhline(0, color='black', lw=0.5)
        ax.axvline(0, color='black', lw=0.5)
        ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1))
        plt.tight_layout()
        plt.show()

    return HSec, MASS

#----------------------------------------------------------------------------   

#-----------------------------------------------------------    
# CONFINED CONCRETE FIBER CIRCULAR SECTION WITH QUAD FIBERS
#-----------------------------------------------------------

def CONCRETE_CONFINED_CIRCULAR_SECTION(SecTag, DSec, coverSec, numBarsSec, barAreaSec, IDconcCore, IDconcCover, IDreinf, ri=0.0, PLOT=True, CONCRETE_DENSITY=2500/1e9, STEEL_DENSITY=7850/1e9):
    """
    Creates a circular reinforced concrete fiber section in OpenSees and plots it using matplotlib.
    Uses plt.Circle() for concrete fibers and reinforcement bars.

    Parameters:
    SecTag (int): Section tag identifier
    DSec (float): Diameter of the section
    coverSec (float): Cover thickness from concrete surface to reinforcement
    numBarsSec (int): Number of longitudinal reinforcing bars
    barAreaSec (float): Cross-sectional area of each reinforcing bar
    IDconcCore (int): Material tag for confined core concrete
    IDconcCover (int): Material tag for unconfined cover concrete
    IDreinf (int): Material tag for reinforcing steel
    ri (float): Inner radius (for hollow sections, default=0.0)
    PLOT (bool): Whether to plot the section (default=True)
    CONCRETE_DENSITY (float): Density of concrete (default=2500 kg/m³ converted to consistent units)
    STEEL_DENSITY (float): Density of steel (default=7850 kg/m³ converted to consistent units)
    """
    nfCoreR = 8  # Number of radial divisions in the core
    nfCoreT = 8  # Number of angular divisions in the core
    nfCoverR = 4  # Number of radial divisions in the cover
    nfCoverT = 8  # Number of angular divisions in the cover

    # Calculate geometric parameters
    ro = DSec / 2.0          # Outer radius
    rc = ro - coverSec       # Core radius (to reinforcement center)
    MASS = CONCRETE_DENSITY * (np.pi * DSec**2) / 4 + numBarsSec * barAreaSec * STEEL_DENSITY

    # Verify valid core radius
    if rc <= ri:
        raise ValueError("Invalid core radius - check cover and diameter values")
    
    # Create fiber section
    ops.section('Fiber', SecTag)
    
    # Confined core concrete patch
    ops.patch('circ', IDconcCore, nfCoreT, nfCoreR, 0.0, 0.0, ri, rc, 0.0, 360.0)
    
    # Unconfined cover concrete patch
    ops.patch('circ', IDconcCover, nfCoverT, nfCoverR, 0.0, 0.0, rc, ro, 0.0, 360.0)
    
    # Reinforcing steel layer
    ops.layer('circ', IDreinf, numBarsSec, barAreaSec, 0.0, 0.0, rc, 0.0, 360.0)
    
    # Plotting the section
    import matplotlib.pyplot as plt
    if PLOT:
        fig, ax = plt.subplots(figsize=(8, 8))

        # Plot concrete fibers
        for i in range(nfCoreR + nfCoverR):
            if i < nfCoreR:
                r_loc = ri + (i + 0.5) * (rc - ri) / nfCoreR
                color = 'lightblue'  # Confined core concrete
            else:
                r_loc = rc + (i - nfCoreR + 0.5) * (ro - rc) / nfCoverR
                color = 'gray'  # Unconfined cover concrete

            for j in range(nfCoreT + nfCoverT):
                theta = j * (2 * np.pi / (nfCoreT + nfCoverT))
                x = r_loc * np.cos(theta)
                y = r_loc * np.sin(theta)
                ax.add_patch(plt.Circle((x, y), 2, color=color, alpha=0.5))

        # Plot rebar fibers
        rebar_radius = max(2, np.sqrt(barAreaSec / np.pi))  # Ensure visibility
        for j in range(numBarsSec):
            theta = j * (2 * np.pi / numBarsSec)
            x = rc * np.cos(theta)
            y = rc * np.sin(theta)
            ax.add_patch(plt.Circle((x, y), rebar_radius, color='red', alpha=0.8))

        # Plot concrete boundaries
        outer_circle = plt.Circle((0, 0), ro, color='black', fill=False, linewidth=1.5)
        inner_circle = plt.Circle((0, 0), ri, color='black', fill=False, linewidth=1.5)
        ax.add_patch(outer_circle)
        ax.add_patch(inner_circle)

        ax.set_aspect('equal')
        ax.set_xlim(-ro - 20, ro + 20)
        ax.set_ylim(-ro - 20, ro + 20)
        plt.xlabel('X (mm)')
        plt.ylabel('Y (mm)')
        plt.title(f'Circular Concrete Section (Tag: {SecTag})')
        plt.grid(True)
        plt.show()
    
    return DSec, MASS


#----------------------------------------------------------------------------   

#-------------------------------------------------------------------    
# CONFINED CONCRETE FIBER CIRCULAR WITH PLATE SECTION - QUAD FIBERS
#-------------------------------------------------------------------


def CONCRETE_CONFINED_CIRCULAR_PLATE_SECTION(SecTag, DSec, coverSec, pipe_thickness, numBarsSec, barAreaSec, IDconcCore, IDconcCover, IDreinf, IDsteelPipe, ri=0.0, PLOT=True, CONCRETE_DENSITY=2500/1e9, STEEL_DENSITY=7850/1e9):
    """
    Creates a circular reinforced concrete fiber section with a surrounding steel pipe in OpenSees and plots it using matplotlib.
    Uses plt.Circle() for concrete fibers, reinforcement bars, and steel pipe fibers.

    Parameters:
    SecTag (int): Section tag identifier
    DSec (float): Diameter of the section
    coverSec (float): Cover thickness from concrete surface to reinforcement
    numBarsSec (int): Number of longitudinal reinforcing bars
    barAreaSec (float): Cross-sectional area of each reinforcing bar
    IDconcCore (int): Material tag for confined core concrete
    IDconcCover (int): Material tag for unconfined cover concrete
    IDreinf (int): Material tag for reinforcing steel
    IDsteelPipe (int): Material tag for the steel pipe
    ri (float): Inner radius (for hollow sections, default=0.0)
    PLOT (bool): Whether to plot the section (default=True)
    CONCRETE_DENSITY (float): Density of concrete (default=2500 kg/m³ converted to consistent units)
    STEEL_DENSITY (float): Density of steel (default=7850 kg/m³ converted to consistent units)
    """
    #import openseespy.opensees as ops
    import matplotlib.pyplot as plt
    import numpy as np
    
    nfCoreR = 8  # Number of radial divisions in the core
    nfCoreT = 8  # Number of angular divisions in the core
    nfCoverR = 4  # Number of radial divisions in the cover
    nfCoverT = 8  # Number of angular divisions in the cover

    # Calculate geometric parameters
    ro = DSec / 2.0          # Outer radius of concrete
    rc = ro - coverSec       # Core radius (to reinforcement center)
    pipe_outer_radius = ro + pipe_thickness  # Outer radius of the steel pipe
    MASS = CONCRETE_DENSITY * (np.pi * DSec**2) / 4 + numBarsSec * barAreaSec * STEEL_DENSITY

    # Verify valid core radius
    if rc <= ri:
        raise ValueError("Invalid core radius - check cover and diameter values")
    
    # Create fiber section
    ops.section('Fiber', SecTag)
    
    # Confined core concrete patch
    ops.patch('circ', IDconcCore, nfCoreT, nfCoreR, 0.0, 0.0, ri, rc, 0.0, 360.0)
    
    # Unconfined cover concrete patch
    ops.patch('circ', IDconcCover, nfCoverT, nfCoverR, 0.0, 0.0, rc, ro, 0.0, 360.0)
    
    # Reinforcing steel layer
    ops.layer('circ', IDreinf, numBarsSec, barAreaSec, 0.0, 0.0, rc, 0.0, 360.0)
    
    # Steel pipe layer
    ops.patch('circ', IDsteelPipe, nfCoreT, 2, 0.0, 0.0, ro, pipe_outer_radius, 0.0, 360.0)
    
    # Plotting the section
    if PLOT:
        fig, ax = plt.subplots(figsize=(8, 8))

        # Plot concrete fibers
        for i in range(nfCoreR + nfCoverR):
            if i < nfCoreR:
                r_loc = ri + (i + 0.5) * (rc - ri) / nfCoreR
                color = 'lightblue'  # Confined core concrete
            else:
                r_loc = rc + (i - nfCoreR + 0.5) * (ro - rc) / nfCoverR
                color = 'gray'  # Unconfined cover concrete

            for j in range(nfCoreT + nfCoverT):
                theta = j * (2 * np.pi / (nfCoreT + nfCoverT))
                x = r_loc * np.cos(theta)
                y = r_loc * np.sin(theta)
                ax.add_patch(plt.Circle((x, y), 2, color=color, alpha=0.5))

        # Plot rebar fibers
        rebar_radius = max(2, np.sqrt(barAreaSec / np.pi))  # Ensure visibility
        for j in range(numBarsSec):
            theta = j * (2 * np.pi / numBarsSec)
            x = rc * np.cos(theta)
            y = rc * np.sin(theta)
            ax.add_patch(plt.Circle((x, y), rebar_radius, color='red', alpha=0.8))

        # Plot steel pipe fibers
        for i in range(2):  # Two radial divisions for the pipe
            r_loc = ro + (i + 0.5) * pipe_thickness / 2
            for j in range(nfCoreT):
                theta = j * (2 * np.pi / nfCoreT)
                x = r_loc * np.cos(theta)
                y = r_loc * np.sin(theta)
                ax.add_patch(plt.Circle((x, y), 2, color='blue', alpha=0.5))

        # Plot concrete boundaries
        outer_circle = plt.Circle((0, 0), ro, color='black', fill=False, linewidth=1.5)
        inner_circle = plt.Circle((0, 0), ri, color='black', fill=False, linewidth=1.5)
        ax.add_patch(outer_circle)
        ax.add_patch(inner_circle)

        # Plot steel pipe boundaries
        pipe_outer_circle = plt.Circle((0, 0), pipe_outer_radius, color='blue', fill=False, linewidth=1.5)
        ax.add_patch(pipe_outer_circle)

        ax.set_aspect('equal')
        ax.set_xlim(-pipe_outer_radius - 20, pipe_outer_radius + 20)
        ax.set_ylim(-pipe_outer_radius - 20, pipe_outer_radius + 20)
        plt.xlabel('X (mm)')
        plt.ylabel('Y (mm)')
        plt.title(f'Circular Concrete Section with Steel Pipe (Tag: {SecTag})')
        plt.grid(True)
        plt.show()
    
    return DSec, MASS

#----------------------------------------------------------------------------   

#----------------------------------------------------------------------    
# CONFINED CONCRETE FIBER RECTANGULAR WITH PLATE SECTION - QUAD FIBERS
#----------------------------------------------------------------------

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
    ops.layer('straight', steelID, numBarsInt, barAreaInt, -coreY + DD, coreZ, coreY - DD, coreZ)  # Intermediate skin reinforcement (+z)
    ops.layer('straight', steelID, numBarsInt, barAreaInt, -coreY + DD, -coreZ, coreY - DD, -coreZ)  # Intermediate skin reinforcement (-z)
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
        ax.scatter(z_top, [coreY] * numBarsTop, s=50, c='red', marker='o', label='Top Steel')

        # Bottom bars
        z_bot = np.linspace(-coreZ, coreZ, numBarsBot)
        ax.scatter(z_bot, [-coreY] * numBarsBot, s=50, c='red', marker='o', label='Bottom Steel')

        # Intermediate bars (sides)
        y_int = np.linspace(-coreY + DD, coreY - DD, numBarsIntTot // 2)
        ax.scatter([coreZ] * len(y_int), y_int, s=30, c='blue', marker='s', label='Skin Steel')
        ax.scatter([-coreZ] * len(y_int), y_int, s=30, c='blue', marker='s')

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

#---------------------------------------------------------------------------- 
    
#-----------------------------------------------------    
# CONFINED CONCRETE FIBER I SECTION WITH QUAD FIBERS
#-----------------------------------------------------   
    
def CONCRETE_CONFINED_I_SECTION_QUAD(id, HSec, BFlange, TFlange, TWeb, coverFlange, coverWeb, coreID, coverID, steelID,
                                    numBarsTop, barAreaTop, numBarsBot, barAreaBot, numBarsIntTot, barAreaInt,
                                    PLOT=True, CONCRETE_DENSITY=2500/1e9, STEEL_DENSITY=7850/1e9):
    import openseespy.opensees as ops
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.patches import Rectangle

    # Fiber discretization parameters
    nfCoreY_web = 10    # Fibers along web height
    nfCoreZ_web = 5     # Fibers along web thickness
    nfCoreY_flange = 5  # Fibers along flange thickness
    nfCoreZ_flange = 10 # Fibers along flange width
    nfCoverY = 2        # Fibers in cover along height
    nfCoverZ = 2        # Fibers in cover along width

    # Geometric calculations
    HWeb = HSec - 2 * TFlange  # Web height

    # Mass calculation
    web_area = HWeb * TWeb
    flange_area = BFlange * TFlange * 2  # Two flanges
    concrete_mass = (web_area + flange_area) * CONCRETE_DENSITY
    steel_mass = (numBarsTop * barAreaTop + numBarsBot * barAreaBot + numBarsIntTot * barAreaInt) * STEEL_DENSITY
    MASS = concrete_mass + steel_mass

    # Define coordinates
    # Web
    web_bottom = -HWeb / 2
    web_top = HWeb / 2
    web_left = -TWeb / 2
    web_right = TWeb / 2

    # Web core (no vertical cover)
    web_core_bottom = web_bottom - coverFlange
    web_core_top = web_top + coverFlange
    web_core_left = web_left + coverWeb
    web_core_right = web_right - coverWeb

    # Top flange
    top_flange_bottom = web_top
    top_flange_top = top_flange_bottom + TFlange
    top_flange_left = -BFlange / 2
    top_flange_right = BFlange / 2

    # Top flange core (with cover on all sides)
    top_flange_core_bottom = top_flange_bottom + coverFlange
    top_flange_core_top = top_flange_top - coverFlange
    top_flange_core_left = top_flange_left + coverFlange
    top_flange_core_right = top_flange_right - coverFlange

    # Bottom flange
    bottom_flange_top = web_bottom
    bottom_flange_bottom = bottom_flange_top - TFlange
    bottom_flange_left = -BFlange / 2
    bottom_flange_right = BFlange / 2

    # Bottom flange core (with cover on all sides)
    bottom_flange_core_bottom = bottom_flange_bottom + coverFlange
    bottom_flange_core_top = bottom_flange_top - coverFlange
    bottom_flange_core_left = bottom_flange_left + coverFlange
    bottom_flange_core_right = bottom_flange_right - coverFlange

    # Define fiber section in OpenSees
    ops.section('Fiber', id)

    # Core patches
    # Web core
    ops.patch('quad', coreID, nfCoreZ_web, nfCoreY_web,
              web_core_bottom, web_core_right,
              web_core_bottom, web_core_left,
              web_core_top, web_core_left,
              web_core_top, web_core_right)

    # Top flange core
    ops.patch('quad', coreID, nfCoreZ_flange, nfCoreY_flange,
              top_flange_core_bottom, top_flange_core_right,
              top_flange_core_bottom, top_flange_core_left,
              top_flange_core_top, top_flange_core_left,
              top_flange_core_top, top_flange_core_right)

    # Bottom flange core
    ops.patch('quad', coreID, nfCoreZ_flange, nfCoreY_flange,
              bottom_flange_core_bottom, bottom_flange_core_right,
              bottom_flange_core_bottom, bottom_flange_core_left,
              bottom_flange_core_top, bottom_flange_core_left,
              bottom_flange_core_top, bottom_flange_core_right)

    # Cover patches
    # Web left cover
    ops.patch('quad', coverID, nfCoverZ, nfCoreY_web,
              web_bottom, web_left,
              web_bottom, web_core_left,
              web_top, web_core_left,
              web_top, web_left)

    # Web right cover
    ops.patch('quad', coverID, nfCoverZ, nfCoreY_web,
              web_bottom, web_right,
              web_bottom, web_core_right,
              web_top, web_core_right,
              web_top, web_right)

    # Top flange covers
    # Top cover
    ops.patch('quad', coverID, nfCoreZ_flange, nfCoverY,
              top_flange_core_top, top_flange_right,
              top_flange_core_top, top_flange_left,
              top_flange_top, top_flange_left,
              top_flange_top, top_flange_right)

    # Left cover
    ops.patch('quad', coverID, nfCoverZ, nfCoreY_flange,
              top_flange_bottom, top_flange_left,
              top_flange_bottom, top_flange_core_left,
              top_flange_top, top_flange_core_left,
              top_flange_top, top_flange_left)

    # Right cover
    ops.patch('quad', coverID, nfCoverZ, nfCoreY_flange,
              top_flange_bottom, top_flange_right,
              top_flange_bottom, top_flange_core_right,
              top_flange_top, top_flange_core_right,
              top_flange_top, top_flange_right)

    # Bottom flange covers
    # Bottom cover
    ops.patch('quad', coverID, nfCoreZ_flange, nfCoverY,
              bottom_flange_core_bottom, bottom_flange_right,
              bottom_flange_core_bottom, bottom_flange_left,
              bottom_flange_bottom, bottom_flange_left,
              bottom_flange_bottom, bottom_flange_right)
    # Left cover
    ops.patch('quad', coverID, nfCoverZ, nfCoreY_flange,
              bottom_flange_bottom, bottom_flange_left,
              bottom_flange_bottom, bottom_flange_core_left,
              bottom_flange_top, bottom_flange_core_left,
              bottom_flange_top, bottom_flange_left)

    # Right cover
    ops.patch('quad', coverID, nfCoverZ, nfCoreY_flange,
              bottom_flange_bottom, bottom_flange_right,
              bottom_flange_bottom, bottom_flange_core_right,
              bottom_flange_top, bottom_flange_core_right,
              bottom_flange_top, bottom_flange_right)

    # Reinforcement layers
    y_top_bars = top_flange_top - coverFlange
    z_left_top = top_flange_left + coverFlange
    z_right_top = top_flange_right - coverFlange
    ops.layer('straight', steelID, numBarsTop, barAreaTop,
              y_top_bars, z_left_top, y_top_bars, z_right_top)

    y_bot_bars = bottom_flange_bottom + coverFlange
    z_left_bot = bottom_flange_left + coverFlange
    z_right_bot = bottom_flange_right - coverFlange
    ops.layer('straight', steelID, numBarsBot, barAreaBot,
              y_bot_bars, z_left_bot, y_bot_bars, z_right_bot)

    numBarsInt = numBarsIntTot // 2
    DD = (HWeb - 2 * coverWeb) / (numBarsInt - 1) if numBarsInt > 1 else 0
    y_start = web_bottom + coverWeb
    y_end = web_top - coverWeb
    z_left_int = web_left + coverWeb
    z_right_int = web_right - coverWeb
    ops.layer('straight', steelID, numBarsInt, barAreaInt,
              y_start, z_left_int, y_end, z_left_int)  # Left side
    ops.layer('straight', steelID, numBarsInt, barAreaInt,
              y_start, z_right_int, y_end, z_right_int)  # Right side

    # Plotting
    if PLOT:
        fig, ax = plt.subplots(figsize=(10, 6))

        # Plot core regions
        ax.add_patch(Rectangle((web_core_left, web_core_bottom),
                     web_core_right - web_core_left, web_core_top - web_core_bottom,
                     facecolor='grey', edgecolor='black', label='Web Core'))
        ax.add_patch(Rectangle((top_flange_core_left, top_flange_core_bottom),
                     top_flange_core_right - top_flange_core_left, top_flange_core_top - top_flange_core_bottom,
                     facecolor='grey', edgecolor='black', label='Top Flange Core'))
        ax.add_patch(Rectangle((bottom_flange_core_left, bottom_flange_core_bottom),
                     bottom_flange_core_right - bottom_flange_core_left, bottom_flange_core_top - bottom_flange_core_bottom,
                     facecolor='grey', edgecolor='black', label='Bottom Flange Core'))

        # Plot cover regions
        # Web covers
        ax.add_patch(Rectangle((web_left, web_bottom), coverWeb, HWeb,
                     facecolor='lightgrey', edgecolor='black', label='Cover'))
        ax.add_patch(Rectangle((web_core_right, web_bottom), coverWeb, HWeb,
                     facecolor='lightgrey', edgecolor='black'))

        # Top flange covers
        ax.add_patch(Rectangle((top_flange_left, top_flange_core_top), BFlange, coverFlange,
                     facecolor='lightgrey', edgecolor='black'))  # Top
        ax.add_patch(Rectangle((top_flange_left, top_flange_bottom), coverFlange, TFlange,
                     facecolor='lightgrey', edgecolor='black'))  # Left
        ax.add_patch(Rectangle((top_flange_core_right, top_flange_bottom), coverFlange, TFlange,
                     facecolor='lightgrey', edgecolor='black'))  # Right

        # Bottom flange covers
        ax.add_patch(Rectangle((bottom_flange_left, bottom_flange_bottom), BFlange, coverFlange,
                     facecolor='lightgrey', edgecolor='black'))  # Bottom
        ax.add_patch(Rectangle((bottom_flange_left, bottom_flange_bottom), coverFlange, TFlange,
                     facecolor='lightgrey', edgecolor='black'))  # Left
        ax.add_patch(Rectangle((bottom_flange_core_right, bottom_flange_bottom), coverFlange, TFlange,
                     facecolor='lightgrey', edgecolor='black'))  # Right

        # Plot reinforcement
        z_top = np.linspace(z_left_top, z_right_top, numBarsTop)
        ax.scatter(z_top, [y_top_bars]*numBarsTop, s=50, c='red', marker='o', label='Top Steel')
        z_bot = np.linspace(z_left_bot, z_right_bot, numBarsBot)
        ax.scatter(z_bot, [y_bot_bars]*numBarsBot, s=50, c='red', marker='o', label='Bottom Steel')
        y_int = np.linspace(y_start, y_end, numBarsInt)
        ax.scatter([z_left_int]*numBarsInt, y_int, s=30, c='blue', marker='s', label='Skin Steel')
        ax.scatter([z_right_int]*numBarsInt, y_int, s=30, c='blue', marker='s')

        # Plot formatting
        ax.set_aspect('equal')
        ax.set_xlabel('Width (Z, mm)')
        ax.set_ylabel('Height (Y, mm)')
        ax.set_title(f'Confined Concrete I-Section with Connected Cores (ID: {id})')
        ax.grid(True, alpha=0.3)
        ax.axhline(0, color='black', lw=0.5)
        ax.axvline(0, color='black', lw=0.5)
        ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1))
        plt.tight_layout()
        plt.show()

    return HSec, MASS
    
#----------------------------------------------------------------------------    