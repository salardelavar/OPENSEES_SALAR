def RECTANGULAR_HOLLOW_STEEL_SECTION_FUN(secTag, matTag, B, H, t, nFib, DENSITY_STEEL, plot=True):    
    """
    Fiber section definition for a Rectangular Hollow Steel Section (RHS)
    and optional visualization of the section geometry.
    B      # [mm] outer width (along y)
    H      # [mm] outer height (along z)
    t      # [mm] wall thickness
    THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
    """
    import openseespy.opensees as ops
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.patches import Polygon, Rectangle
    from matplotlib.collections import PatchCollection
    
    # -------------------------
    # Basic geometric check
    # -------------------------
    if t >= B/2:
        raise ValueError("Thickness t must be smaller than B/2")
    if t >= H/2:
        raise ValueError("Thickness t must be smaller than H/2")        
    
    # -------------------------
    # Section Area Calculation
    # -------------------------
    AREA = B*H - (B - 2*t)*(H - 2*t)   # [mm^2]  
    
    # -------------------------
    # OpenSees fiber section
    # -------------------------
    ops.section('Fiber', secTag)

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
        # Create figure with better styling
        fig, ax = plt.subplots(figsize=(8, 8))
        
        # Create the hollow section as a grey filled area
        # Method: Create outer rectangle patch and subtract inner rectangle
        outer_rect = Rectangle((y1, z1), B, H, linewidth=2, edgecolor='black', facecolor='none')
        
        # Create the steel area as a grey polygon (ring)
        # We'll create 4 separate patches for each wall
        steel_patches = []
        
        # Top wall
        top_wall = Polygon([
            (y1, z2), (y2, z2), (yi2, zi2), (yi1, zi2)
        ], closed=True, facecolor='lightgray', edgecolor='black', linewidth=1)
        
        # Bottom wall
        bottom_wall = Polygon([
            (yi1, zi1), (yi2, zi1), (y2, z1), (y1, z1)
        ], closed=True, facecolor='lightgray', edgecolor='black', linewidth=1)
        
        # Right wall
        right_wall = Polygon([
            (yi2, zi1), (y2, z1), (y2, z2), (yi2, zi2)
        ], closed=True, facecolor='lightgray', edgecolor='black', linewidth=1)
        
        # Left wall
        left_wall = Polygon([
            (y1, z1), (yi1, zi1), (yi1, zi2), (y1, z2)
        ], closed=True, facecolor='lightgray', edgecolor='black', linewidth=1)
        
        # Add all steel patches
        for patch in [top_wall, bottom_wall, right_wall, left_wall]:
            ax.add_patch(patch)
        
        # Add inner boundary as dashed line for clarity
        inner_rect = Rectangle((yi1, zi1), B-2*t, H-2*t, linewidth=1.5, 
                               edgecolor='darkred', facecolor='none', linestyle='--')
        ax.add_patch(inner_rect)
        
        # Add center marker
        ax.plot(0, 0, 'b+', markersize=10, label='Centroid')
        
        # Add dimensions annotation
        ax.annotate(f'B = {B} mm', xy=(0, z2+5), ha='center', fontsize=9)
        ax.annotate(f'H = {H} mm', xy=(y2+5, 0), va='center', fontsize=9)
        ax.annotate(f't = {t} mm', xy=(B/4, H/4), fontsize=9, 
                   bbox=dict(boxstyle="round,pad=0.3", facecolor="yellow", alpha=0.7))
        
        # Set axis properties
        ax.set_aspect('equal')
        ax.set_xlabel('y (mm)', fontsize=11)
        ax.set_ylabel('z (mm)', fontsize=11)
        ax.set_title(f'Rectangular Hollow Steel Section (RHS)\nB={B}mm, H={H}mm, t={t}mm, Area={AREA:.1f}mm²', 
                    fontsize=12, fontweight='bold')
        
        # Add grid
        ax.grid(True, linestyle=':', alpha=0.6)
        
        # Set axis limits with some padding
        padding = max(B, H) * 0.1
        ax.set_xlim(y1 - padding, y2 + padding)
        ax.set_ylim(z1 - padding, z2 + padding)
        
        # Add legend
        from matplotlib.lines import Line2D
        legend_elements = [
            Line2D([0], [0], color='black', linewidth=2, label='Outer boundary'),
            Line2D([0], [0], color='darkred', linestyle='--', linewidth=1.5, label='Inner boundary'),
            Line2D([0], [0], marker='s', color='w', markerfacecolor='lightgray', markersize=10, label='Steel area'),
            Line2D([0], [0], marker='+', color='b', markersize=10, label='Centroid')
        ]
        ax.legend(handles=legend_elements, loc='upper right', fontsize=9)
        
        # Add dimension lines
        # Horizontal dimension (B)
        ax.plot([y1, y2], [z2+2, z2+2], 'k-', linewidth=1)
        ax.plot([y1, y1], [z2, z2+4], 'k-', linewidth=1)
        ax.plot([y2, y2], [z2, z2+4], 'k-', linewidth=1)
        
        # Vertical dimension (H)
        ax.plot([y2+2, y2+2], [z1, z2], 'k-', linewidth=1)
        ax.plot([y2, y2+4], [z1, z1], 'k-', linewidth=1)
        ax.plot([y2, y2+4], [z2, z2], 'k-', linewidth=1)
        
        plt.tight_layout()
        plt.show()
        
    ELE_MASS  = DENSITY_STEEL * AREA   # Mass Per Length
     
    return AREA, ELE_MASS        