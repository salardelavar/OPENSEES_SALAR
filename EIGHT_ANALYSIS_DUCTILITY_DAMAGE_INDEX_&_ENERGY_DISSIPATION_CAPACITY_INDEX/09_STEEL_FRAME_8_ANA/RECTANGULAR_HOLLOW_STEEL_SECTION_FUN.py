def RECTANGULAR_HOLLOW_STEEL_SECTION_FUN(secTag, matTag, B, H, t, nFib, DENSITY_STEEL, plot=True):
    
    import openseespy.opensees as ops
    import matplotlib.pyplot as plt
    import numpy as np
    # Create patches for the steel area
    from matplotlib.patches import Polygon, Rectangle
    from matplotlib.collections import PatchCollection

    """
    Fiber section definition for a Rectangular Hollow Steel Section (RHS)
    and optional visualization of the section geometry.
    """
    # -------------------------
    # Basic geometric check
    # -------------------------
    if t >= B/2:
        raise ValueError("Thickness t must be smaller than B/2")
    
    # -------------------------
    # Section Area Calculation
    # -------------------------
    AREA = B**2 - (B - 2*t) * (H - 2*t)  # [mm^2] 

    # -------------------------
    # OpenSees fiber section
    # -------------------------
    ops.section('Fiber', secTag) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/fibersection.html

    # Outer coordinates
    y1 = -B / 2.0
    y2 =  B / 2.0
    z1 = -H / 2.0
    z2 =  H / 2.0

    # Inner (hollow) coordinates
    yi1 = y1 + t
    yi2 = y2 - t
    zi1 = z1 + t
    zi2 = z2 - t

    # -------------------------
    # Steel plates
    # -------------------------

    # Top plate
    ops.patch('quad', matTag, nFib, 1,
              y1,  z2,
              y2,  z2,
              yi2, zi2,
              yi1, zi2)

    # Bottom plate
    ops.patch('quad', matTag, nFib, 1,
              yi1, zi1,
              yi2, zi1,
              y2,  z1,
              y1,  z1)

    # Right plate
    ops.patch('quad', matTag, 1, nFib,
              yi2, zi1,
              y2,  z1,
              y2,  z2,
              yi2, zi2)

    # Left plate
    ops.patch('quad', matTag, 1, nFib,
              y1,  z1,
              yi1, zi1,
              yi1, zi2,
              y1,  z2)

    # -------------------------
    # Plot section geometry (improved)
    # -------------------------
    if plot:
        # Create figure with better styling
        fig, ax = plt.subplots(figsize=(8, 8))
                
        # Define the four wall patches
        # Top wall
        top_wall = Polygon([
            [y1, z2], [y2, z2], [yi2, zi2], [yi1, zi2]
        ], closed=True)
        
        # Bottom wall
        bottom_wall = Polygon([
            [yi1, zi1], [yi2, zi1], [y2, z1], [y1, z1]
        ], closed=True)
        
        # Right wall
        right_wall = Polygon([
            [yi2, zi1], [y2, z1], [y2, z2], [yi2, zi2]
        ], closed=True)
        
        # Left wall
        left_wall = Polygon([
            [y1, z1], [yi1, zi1], [yi1, zi2], [y1, z2]
        ], closed=True)
        
        # Create a collection of patches with grey fill
        patches = [top_wall, bottom_wall, right_wall, left_wall]
        collection = PatchCollection(patches, 
                                   facecolor='lightgray', 
                                   edgecolor='darkgray',
                                   linewidth=1,
                                   alpha=0.7)
        
        # Add patches to the plot
        ax.add_collection(collection)
        
        # Add boundary lines for better definition
        # Outer rectangle
        outer_y = [y1, y2, y2, y1, y1]
        outer_z = [z1, z1, z2, z2, z1]
        ax.plot(outer_y, outer_z, 'k-', linewidth=2, label='Outer boundary')
        
        # Inner rectangle (hole)
        inner_y = [yi1, yi2, yi2, yi1, yi1]
        inner_z = [zi1, zi1, zi2, zi2, zi1]
        ax.plot(inner_y, inner_z, 'r--', linewidth=2, label='Inner boundary')
        
        # Add dimensions
        # Center lines
        ax.axhline(y=0, color='gray', linestyle=':', linewidth=0.5, alpha=0.5)
        ax.axvline(x=0, color='gray', linestyle=':', linewidth=0.5, alpha=0.5)
        
        # Add dimension arrows and text
        # Width dimension
        ax.annotate('', xy=(y1, z2+5), xytext=(y2, z2+5),
                   arrowprops=dict(arrowstyle='<->', color='blue'))
        ax.text(0, z2+10, f'B = {B} mm', ha='center', va='bottom', color='blue')
        
        # Height dimension
        ax.annotate('', xy=(y2+5, z1), xytext=(y2+5, z2),
                   arrowprops=dict(arrowstyle='<->', color='blue'))
        ax.text(y2+15, 0, f'H = {H} mm', ha='left', va='center', color='blue', rotation=90)
        
        # Thickness dimension
        # Show thickness at one corner
        corner_x = [y1, yi1]
        corner_y = [z1, zi1]
        ax.plot(corner_x, corner_y, 'g-', linewidth=2)
        ax.annotate('', xy=(y1-5, z1-5), xytext=(yi1-5, zi1-5),
                   arrowprops=dict(arrowstyle='<->', color='green'))
        ax.text(y1-15, z1-20, f't = {t} mm', ha='right', color='green')
        
        # Add corner markers
        ax.plot(y1, z1, 'ko', markersize=4)
        ax.plot(y2, z1, 'ko', markersize=4)
        ax.plot(y2, z2, 'ko', markersize=4)
        ax.plot(y1, z2, 'ko', markersize=4)
        
        # Add centroid marker
        ax.plot(0, 0, 'r+', markersize=10, markeredgewidth=2, label='Centroid')
        
        # Add cross-hatch pattern for steel area (optional)
        # This creates a subtle pattern to distinguish the steel
        for i in range(-int(B), int(B), 10):
            if i > y1 and i < y2:
                ax.axvline(x=i, ymin=(z1+zi1)/2/H+0.5, ymax=(z2+zi2)/2/H+0.5, 
                          color='gray', linewidth=0.3, alpha=0.2)
        
        # Set equal aspect ratio
        ax.set_aspect('equal')
        
        # Set labels and title
        ax.set_xlabel('y (mm)')
        ax.set_ylabel('z (mm)')
        ax.set_title(f'Rectangular Hollow Steel Section\nB={B}mm, H={H}mm, t={t}mm')
        
        # Add grid
        ax.grid(True, alpha=0.3)
        
        # Add legend
        ax.legend(loc='upper right')
        
        # Set margins to show dimensions clearly
        margin = max(B, H) * 0.15
        ax.set_xlim([y1 - margin, y2 + margin])
        ax.set_ylim([z1 - margin, z2 + margin])
        
        # Add section properties box
        area = 2 * t * (B + H - 2 * t)  # Approximate area
        props_text = f'Section Properties:\nArea ≈ {area:.0f} mm²'
        ax.text(0.02, 0.98, props_text, transform=ax.transAxes,
               verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        plt.tight_layout()
        plt.show()

    ELE_MASS  = DENSITY_STEEL * AREA   # Mass Per Length
     
    return AREA, ELE_MASS