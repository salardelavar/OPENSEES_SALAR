def SQUARE_HOLLOW_STEEL_SECTION_FUN(secTag, matTag, B, t, nFib, DENSITY_STEEL, plot=True):  
    """
    Fiber section definition for a Square Hollow Steel Section (SHS)
    and optional visualization of the section geometry.
    B      # [mm] outer width (along y and z)
    t      # [mm] wall thickness
    THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
    """   
    import openseespy.opensees as ops
    import matplotlib.pyplot as plt
    import numpy as np
    # Create patches for the steel area
    from matplotlib.patches import Polygon, Rectangle
    from matplotlib.collections import PatchCollection

    # -------------------------
    # Basic geometric check
    # -------------------------
    if t >= B/2:
        raise ValueError("Thickness t must be smaller than B/2")
    
    # -------------------------
    # Section Area Calculation
    # -------------------------
    AREA = B**2 - (B - 2*t)**2   # [mm^2]  
    
    # -------------------------
    # OpenSees fiber section
    # -------------------------
    ops.section('Fiber', secTag) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/fibersection.html

    # Outer coordinates
    y1 = -B / 2.0
    y2 =  B / 2.0
    z1 = -B / 2.0
    z2 =  B / 2.0

    # Inner (hollow) coordinates
    yi1 = y1 + t
    yi2 = y2 - t
    zi1 = z1 + t
    zi2 = z2 - t

    # Top plate
    ops.patch('quad', matTag, nFib, 5,
              y1,  z2,
              y2,  z2,
              yi2, zi2,
              yi1, zi2)

    # Bottom plate
    ops.patch('quad', matTag, nFib, 5,
              yi1, zi1,
              yi2, zi1,
              y2,  z1,
              y1,  z1)

    # Right plate
    ops.patch('quad', matTag, 5, nFib,
              yi2, zi1,
              y2,  z1,
              y2,  z2,
              yi2, zi2)

    # Left plate
    ops.patch('quad', matTag, 5, nFib,
              y1,  z1,
              yi1, zi1,
              yi1, zi2,
              y1,  z2)

    # -------------------------
    # Plot section geometry
    # -------------------------
    if plot:
        # Create figure with better size
        fig, ax = plt.subplots(figsize=(8, 8))
        
        # Create the steel section as a grey-filled area
        # Method 1: Using patch collections for the four walls
        patches = []
        
        # Top wall
        top_wall = Polygon([(y1, z2), (y2, z2), (yi2, zi2), (yi1, zi2)], 
                          closed=True, facecolor='lightgray', edgecolor='black', 
                          linewidth=1.5, alpha=0.8)
        patches.append(top_wall)
        
        # Bottom wall
        bottom_wall = Polygon([(yi1, zi1), (yi2, zi1), (y2, z1), (y1, z1)], 
                             closed=True, facecolor='lightgray', edgecolor='black', 
                             linewidth=1.5, alpha=0.8)
        patches.append(bottom_wall)
        
        # Right wall
        right_wall = Polygon([(yi2, zi1), (y2, z1), (y2, z2), (yi2, zi2)], 
                            closed=True, facecolor='lightgray', edgecolor='black', 
                            linewidth=1.5, alpha=0.8)
        patches.append(right_wall)
        
        # Left wall
        left_wall = Polygon([(y1, z1), (yi1, zi1), (yi1, zi2), (y1, z2)], 
                           closed=True, facecolor='lightgray', edgecolor='black', 
                           linewidth=1.5, alpha=0.8)
        patches.append(left_wall)
        
        # Add all patches to the axis
        for patch in patches:
            ax.add_patch(patch)
        
        # Mark the centroid
        ax.plot(0, 0, 'ro', markersize=8, label='Centroid', zorder=5)
        
        # Add dimension annotations
        # Width dimension
        ax.annotate('', xy=(y1, z1-0.1*B), xytext=(y2, z1-0.1*B),
                   arrowprops=dict(arrowstyle='<->', color='blue', lw=1.5))
        ax.text(0, z1-0.15*B, f'B = {B} mm', ha='center', va='top', 
               fontsize=10, bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.7))
        
        # Height dimension (optional - same as B for square)
        ax.annotate('', xy=(y1-0.1*B, z1), xytext=(y1-0.1*B, z2),
                   arrowprops=dict(arrowstyle='<->', color='blue', lw=1.5))
        ax.text(y1-0.15*B, 0, f'B = {B} mm', ha='center', va='center', 
               rotation=90, fontsize=10, bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.7))
        
        # Thickness annotation
        mid_y = (y1 + yi1) / 2
        mid_z = (z2 + zi2) / 2
        ax.annotate('', xy=(y1, z2-0.05*B), xytext=(yi1, zi2-0.05*B),
                   arrowprops=dict(arrowstyle='<->', color='green', lw=1.5))
        ax.text(mid_y, mid_z+0.05*B, f't = {t} mm', ha='center', va='bottom',
               fontsize=9, bbox=dict(boxstyle='round,pad=0.2', facecolor='lightgreen', alpha=0.7))
        
        # Set axis properties
        ax.set_aspect('equal')
        ax.set_xlabel('y (mm)', fontsize=12)
        ax.set_ylabel('z (mm)', fontsize=12)
        ax.set_title(f'Square Hollow Steel Section (SHS)\nB = {B} mm, t = {t} mm', 
                    fontsize=14, fontweight='bold')
        
        # Add grid with custom style
        ax.grid(True, linestyle='--', alpha=0.3, color='gray')
        
        # Set axis limits with some padding
        padding = 0.1 * B
        ax.set_xlim(y1 - padding, y2 + padding)
        ax.set_ylim(z1 - padding, z2 + padding)
        
        # Add legend
        ax.plot([], [], 's', color='lightgray', markersize=15, label='Steel section', alpha=0.8)
        ax.plot([], [], 'r-', linewidth=1.5, label='Outer boundary')
        ax.plot([], [], 'k--', linewidth=1.5, label='Inner boundary')
        ax.legend(loc='upper right', fontsize=10, framealpha=0.9)
        
        # Add section properties text box
        props_text = f"Section Properties:\nAREA = {AREA:.0f} mm²\nMass/length = {DENSITY_STEEL * AREA:.2f}"
        ax.text(0.02, 0.98, props_text, transform=ax.transAxes,
               fontsize=9, verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        
        plt.tight_layout()
        plt.show()

        
    ELE_MASS  = DENSITY_STEEL * AREA   # Mass Per Length
     
    return AREA, ELE_MASS   