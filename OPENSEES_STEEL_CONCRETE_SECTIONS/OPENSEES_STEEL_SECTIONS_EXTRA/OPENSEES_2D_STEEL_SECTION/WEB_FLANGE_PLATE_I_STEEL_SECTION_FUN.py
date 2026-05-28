"""
# --------- Section Properties (in mm) ---------
# I-beam dimensions
bf = 200     # Flange width
tf = 20      # Flange thickness
h  = 400     # Total web height
tw = 12      # Web thickness

# Stiffener plates dimensions on flanges
tsF = 20      # Stiffener thickness
hsF = 120     # Stiffener height on web
bsF = 100     # Stiffener length on flange

# Stiffener plates dimensions on webs
tsW = 10      # Stiffener thickness
hsW = 160     # Stiffener height on web
bsW = 60      # Stiffener length on flange
THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
"""
def WEB_FLANGE_PLATE_I_STEEL_SECTION_FUN(secTag, matTag, bf, tf, h, tw,
                                   tsF, bsF,
                                   tsW, hsW, bsW,
                                   nFib, DENSITY_STEEL,
                                   plot=True):
    
    # THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)      
    import openseespy.opensees as ops
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    # -------------------------
    # Section Area Calculation
    # -------------------------
    # Area of main I-beam components
    area_top_flange = bf * tf
    area_bottom_flange = bf * tf
    area_web = tw * (h - 2 * tf) # Height of web is total height minus thicknesses of top and bottom flanges
    
    # Area of stiffener plates
    area_left_web_stiffener = tsW * hsW
    area_right_web_stiffener = tsW * hsW
    area_top_flange_stiffener = 2 * bsF * tsF
    AREA = (area_top_flange + area_bottom_flange + area_web +
              area_left_web_stiffener + area_right_web_stiffener +
              area_top_flange_stiffener)
    
    print(f"Area of Top Flange: {area_top_flange:.2f} mm^2")
    print(f"Area of Bottom Flange: {area_bottom_flange:.2f} mm^2")
    print(f"Area of Web: {area_web:.2f} mm^2")
    print(f"Area of Left Web Stiffener: {area_left_web_stiffener:.2f} mm^2")
    print(f"Area of Right Web Stiffener: {area_right_web_stiffener:.2f} mm^2")
    print(f"Area of Top Flange Stiffener: {area_top_flange_stiffener:.2f} mm^2")
    print(f"Total Section Area: {AREA:.2f} mm^2")

    # -------------------------
    # OpenSees fiber section
    # -------------------------
    ops.section('Fiber', secTag)
    
    # --------- Define Fibers for Flanges and Web ---------
    def add_rectangle_patch(y_bottom, y_top, x_left, x_right, num_y_divisions=nFib, num_x_divisions=nFib):
        """Helper function to add rectangular patches of fibers."""
        ops.patch('rect', matTag, num_y_divisions, num_x_divisions, x_left, y_bottom, x_right, y_top)
    
    #%% I SECTION
    # Top flange
    add_rectangle_patch(h/2 - tf, h/2, -bf/2, bf/2)
    
    # Bottom flange
    add_rectangle_patch(-h/2, -h/2 + tf, -bf/2, bf/2)
    
    # Web
    add_rectangle_patch(-h/2 + tf, h/2 - tf, -tw/2, tw/2)
    
    # --------- Add Stiffener Plates ---------
    # Stiffeners on the web (left and right sides)
    # Left stiffener
    add_rectangle_patch(-hsW/2, hsW/2, -tw/2 - tsW, -tw/2)
    # Right stiffener
    add_rectangle_patch(-hsW/2, hsW/2, tw/2, tw/2 + tsW)
    # Stiffener on the top flange
    add_rectangle_patch(h/2, h/2 + tsF, -bsF/2, bsF/2)
    add_rectangle_patch(-h/2, -h/2 - tsF, -bsF/2, bsF/2)
    
    
    # -------------------------
    # Plot section geometry
    # -------------------------
    if plot:
        fig, ax = plt.subplots(figsize=(8, 8)) # Increased figure size for better visibility
        ax.set_xlabel("Width (mm)")
        ax.set_ylabel("Height (mm)")
        ax.set_title('I-Section with Stiffener Plates')
        ax.grid(True, linestyle='--', alpha=0.6)
        
        # Define colors
        main_color = '#778899'  # Light Slate Gray for main I-beam parts
        stiffener_color = '#4682B4' # Steel Blue for stiffeners
        outline_color = 'black'
        
        # Define the geometry of all parts for drawing
        # (x_bottom_left, y_bottom_left, width, height)
        geometry = [
            # Main I-beam parts
            {"x": -bf/2, "y": h/2 - tf, "w": bf, "h": tf, "color": main_color, "edgecolor": outline_color, "fill": True}, # Top flange
            {"x": -bf/2, "y": -h/2, "w": bf, "h": tf, "color": main_color, "edgecolor": outline_color, "fill": True},      # Bottom flange
            {"x": -tw/2, "y": -h/2 + tf, "w": tw, "h": h - 2*tf, "color": main_color, "edgecolor": outline_color, "fill": True}, # Web
        
            # Stiffener plates
            {"x": -tw/2 - tsW, "y": -hsW/2, "w": tsW, "h": hsW, "color": stiffener_color, "edgecolor": outline_color, "fill": True}, # Left web stiffener
            {"x": tw/2, "y": -hsW/2, "w": tsW, "h": hsW, "color": stiffener_color, "edgecolor": outline_color, "fill": True},      # Right web stiffener
            {"x": -bsF/2, "y": h/2, "w": bsF, "h": tsF, "color": stiffener_color, "edgecolor": outline_color, "fill": True},       # Top flange stiffener
            {"x": -bsF/2, "y": -h/2, "w": bsF, "h": -tsF, "color": stiffener_color, "edgecolor": outline_color, "fill": True}       # Bottom flange stiffener
        ]
        
        # Draw all parts with specified colors
        for part in geometry:
            rect = patches.Rectangle(
                (part["x"], part["y"]),
                part["w"],
                part["h"],
                linewidth=1.5,
                edgecolor=part["edgecolor"],
                facecolor=part["color"] if part["fill"] else 'none',
                label=part.get("label", "") # Add label if present
            )
            ax.add_patch(rect)
        
        # Set plot limits and aspect ratio
        max_dim = max(bf, h, bsF) + tsF + 5 # Add some padding
        ax.set_xlim(-max_dim / 2, max_dim / 2)
        ax.set_ylim(-max_dim / 2, max_dim / 2)
        ax.set_aspect('equal', adjustable='box')
        
        # Add a legend to identify the parts
        handles = [
            patches.Patch(color=main_color, label='Main I-Beam Section'),
            patches.Patch(color=stiffener_color, label='Stiffener Plates')
        ]
        ax.legend(handles=handles, loc='upper right')
        
        plt.show()
        
    ELE_MASS  = DENSITY_STEEL * AREA   # Mass Per Length
     
    return h, ELE_MASS           