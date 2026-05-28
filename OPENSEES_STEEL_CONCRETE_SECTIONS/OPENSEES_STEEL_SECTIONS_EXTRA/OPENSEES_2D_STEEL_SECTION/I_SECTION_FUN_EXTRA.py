def I_SECTION_FUN(secTag, matTag, PLOT=True, DENSITY=7850/1e9):
    # THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
    import numpy as np
    import openseespy.opensees as ops
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
    fb = bf * tf / 4 # Fiber Area
    x1 = bf / 4
    y1 = hw / 2 + tf / 4
    x2 = -bf / 4
    y2 = hw / 2 + tf / 2 + tf / 4
    ops.fiber(x1, y1, fb, matTag)    # LAYER 01 RIGHT
    ops.fiber(x1, y2, fb, matTag)    # LAYER 02 RIGHT
    ops.fiber(x2, y1, fb, matTag)    # LAYER 03 LEFT
    ops.fiber(x2, y2, fb, matTag)    # LAYER 04 LEFT
    # Bottom flange fibers
    x1 = bf / 4
    y1 = -hw / 2 - tf / 4
    x2 = -bf / 4
    y2 = -hw / 2 - tf / 2 - tf / 4
    ops.fiber(x1, y1, fb, matTag)    # LAYER 01 RIGHT
    ops.fiber(x1, y2, fb, matTag)    # LAYER 02 RIGHT
    ops.fiber(x2, y1, fb, matTag)    # LAYER 03 LEFT
    ops.fiber(x2, y2, fb, matTag)    # LAYER 04 LEFT
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