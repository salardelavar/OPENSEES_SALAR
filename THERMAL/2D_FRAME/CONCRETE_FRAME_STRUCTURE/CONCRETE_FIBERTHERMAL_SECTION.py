import numpy as np
import openseespy.opensees as ops
#----------------------------------------------
# Concrete Thermal Fiber Section without Rebars 
#----------------------------------------------
def R_RECTANGULAR_CONCRETE_SECTION(secTag, B, H, NUM_B, NUM_H, matTag):
    #matTag = 1    # Material tag for the concrete
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

    return H
    
#-----------------------------------------------------------    
# Concrete Thermal Fiber Section with Rebars (NUM_LAYERS)
#-----------------------------------------------------------
def R_RECTANGULAR_CONCRETE_SECTION_REBAR_FF(secTag, B, H, COVER, RD, NUM_B, NUM_H, NUM_LAYERS, matTag_concrete, matTag_steel, PLOT=True):
    # Define the concrete section using fibers
    ops.section('FiberThermal', secTag)

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

    return H # Return Section Height