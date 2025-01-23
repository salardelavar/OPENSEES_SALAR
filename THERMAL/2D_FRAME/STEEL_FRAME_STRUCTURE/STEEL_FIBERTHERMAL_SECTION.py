import numpy as np
import openseespy.opensees as ops
#----------------------------
# I SECTION
#----------------------------
def I_SECTION(secTag, matTag, PLOT=True):
    #secTag = 1
    # Define geometric properties of the steel I section
    ops.section('FiberThermal', secTag)
    # Define section (Fiber)
    bf = 300               # [mm] Flange width
    tf = 20                # [mm] Flange thickness
    tw = 10                # [mm] Web thickness
    hw = 400               # [mm] Web height
    NUM = 100              # Number of fibers for web
    d = 2 * tf + hw
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
        ax.add_patch(plt.Rectangle((-bf / 2, hw / 2), bf, tf, color='lightgrey', alpha=0.7))

        # Plot the bottom flange fibers
        ax.add_patch(plt.Rectangle((-bf / 2, -hw / 2 - tf), bf, tf, color='lightgrey', alpha=0.7))

        # Plot the web fibers
        fiber_height = hw / NUM
        for i in range(NUM):
            yLoc = hw / 2 - i * fiber_height
            ax.add_patch(plt.Rectangle((-tw / 2, yLoc - fiber_height / 2), tw, fiber_height, color='lightgrey', alpha=0.7))

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
    
    return d # Return Section Height
    
#----------------------------    
# CIRCULAR TUBE SECTION  
#---------------------------- 
def C_TUBE_SECTION(secTag, matTag, PLOT=True):
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

    # Define the steel tube section using fibers
    ops.section('FiberThermal', secTag)

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
                ax.add_patch(plt.Circle((x_loc, y_loc), 2, color='lightgrey', alpha=0.7))  # Fibers as small circles

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
    
    return D # Return Section Height

#----------------------------
# RECTANGULAR TUBE SECTION 
#----------------------------
def R_TUBE_SECTION(secTag, matTag, PLOT=True):
    # Define section tag for rectangular tube steel section
    #secTag = 1
    # Define geometric properties of the rectangular tube
    B = 300  # [mm] Outer width of the tube
    H = 500  # [mm] Outer height of the tube
    t = 20   # [mm] Wall thickness of the tube
    # Define material tag
    matTag = 1
    # Number of fibers along each wall direction
    NUM_B = 10  # Number of fibers along the width
    NUM_H = 20  # Number of fibers along the height
    # Define the rectangular tube section using fibers
    ops.section('FiberThermal', secTag)
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
            ax.add_patch(plt.Rectangle((x_loc, H / 2 - t), fiber_width, t, color='lightgrey', alpha=0.7))

        # Plot outer bottom wall fibers
        for i in range(NUM_B):
            x_loc = -B / 2 + i * fiber_width
            ax.add_patch(plt.Rectangle((x_loc, -H / 2), fiber_width, t, color='lightgrey', alpha=0.7))

        # Plot outer left wall fibers
        fiber_height = H / NUM_H
        for i in range(NUM_H):
            y_loc = -H / 2 + i * fiber_height
            ax.add_patch(plt.Rectangle((-B / 2, y_loc), t, fiber_height, color='lightgrey', alpha=0.7))

        # Plot outer right wall fibers
        for i in range(NUM_H):
            y_loc = -H / 2 + i * fiber_height
            ax.add_patch(plt.Rectangle((B / 2 - t, y_loc), t, fiber_height, color='lightgrey', alpha=0.7))

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
    
    return H # Return Section Height
#---------------------------
# RECTANGULAR SECTION 
#---------------------------
def R_RECTANGULAR_STEEL_SECTION(secTag, matTag, PLOT=True):
    # Define a rectangular steel section using fibers.
    B = 400       # [mm] Width of the rectangular section
    H = 500       # [mm] Height of the rectangular section
    NUM_B = 10   # Number of fibers along the width of the section
    NUM_H = 10   # Number of fibers along the height of the section
    secTag = 1    # Section tag identifier
    matTag = 1    # Material tag for the concrete
    ops.section('FiberThermal', secTag)
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
                rect = plt.Rectangle((x_loc, y_loc), fiber_width, fiber_height, color='lightgrey', edgecolor='black', alpha=0.7)
                ax.add_patch(rect)
        
        # Add labels and grid
        ax.set_xlabel("Width [mm]")
        ax.set_ylabel("Height [mm]")
        ax.grid(True)
        ax.set_title("Rectangular Steel Section with Fibers")
        
        # Show the plot
        plt.show()           

    return H # Return Section Height
    