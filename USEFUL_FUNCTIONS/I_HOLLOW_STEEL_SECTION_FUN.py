def I_HOLLOW_STEEL_SECTION_FUN(secTag, matTag, D, B, tf, tw, t_hollow, nFibWeb, nFibFlange, DENSITY_STEEL, plot=True):
    """
    Fiber section definition for an I-shaped Hollow Steel Section 
    (built-up I-section with hollow flanges and web)
    
    Parameters:
    -----------
    secTag : int
        Section tag for OpenSees
    matTag : int
        Material tag for steel
    D : float
        Total depth of section (mm)
    B : float
        Flange width (mm)
    tf : float
        Flange thickness (mm)
    tw : float
        Web thickness (mm)
    t_hollow : float
        Wall thickness of hollow sections (mm)
    nFibWeb : int
        Number of fibers along web depth
    nFibFlange : int
        Number of fibers along flange width
    DENSITY_STEEL : float
        Steel density (kg/mm³ or consistent units)
    plot : bool
        Whether to generate visualization
    """
    
    import openseespy.opensees as ops
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.patches import Polygon, Rectangle
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Patch
    
    # -------------------------
    # Basic geometric checks
    # -------------------------
    if 2*t_hollow >= tf:
        raise ValueError(f"2*t_hollow ({2*t_hollow}) must be less than tf ({tf})")
    if 2*t_hollow >= tw:
        raise ValueError(f"2*t_hollow ({2*t_hollow}) must be less than tw ({tw})")
    
    # -------------------------
    # Calculate section properties
    # -------------------------
    h_web = D - 2*tf  # clear web height
    
    # Gross area (solid I-section)
    gross_area = 2*B*tf + h_web*tw
    
    # Hollow cavities area
    flange_cavity_area = (B - 2*t_hollow) * (tf - 2*t_hollow)
    web_cavity_area = (h_web - 2*t_hollow) * (tw - 2*t_hollow)
    hollow_area = 2 * flange_cavity_area + web_cavity_area
    AREA = gross_area - hollow_area
    
    # -------------------------
    # OpenSees fiber section
    # -------------------------
    ops.section('Fiber', secTag)
    
    # Coordinate definitions
    y_left = -B/2.0
    y_right = B/2.0
    z_bottom = -D/2.0
    z_top = D/2.0
    
    z_flange_bottom = -D/2.0 + tf
    z_flange_top = D/2.0 - tf
    
    y_web_left = -tw/2.0
    y_web_right = tw/2.0
    
    # Cavity boundaries
    # Top flange cavity
    y_flange_cavity_left = y_left + t_hollow
    y_flange_cavity_right = y_right - t_hollow
    z_cavity_top = z_top - t_hollow
    z_flange_cavity_bottom = z_flange_top + t_hollow  # bottom of top cavity
    
    # Bottom flange cavity
    z_cavity_bottom = z_bottom + t_hollow
    z_flange_cavity_top = z_flange_bottom - t_hollow  # top of bottom cavity
    
    # Web cavity
    y_web_cavity_left = y_web_left + t_hollow
    y_web_cavity_right = y_web_right - t_hollow
    z_web_cavity_bottom = z_flange_bottom + t_hollow
    z_web_cavity_top = z_flange_top - t_hollow
    
    # -------------------------
    # Steel patches (all solid areas)
    # -------------------------
    steel_patches = []
    
    # Top flange - left outer part (from left edge to web)
    steel_patches.append(Polygon([
        [y_left, z_flange_top],
        [y_web_left, z_flange_top],
        [y_web_left, z_top],
        [y_left, z_top]
    ], closed=True))
    
    # Top flange - right outer part (from web to right edge)
    steel_patches.append(Polygon([
        [y_web_right, z_flange_top],
        [y_right, z_flange_top],
        [y_right, z_top],
        [y_web_right, z_top]
    ], closed=True))
    
    # Bottom flange - left outer part
    steel_patches.append(Polygon([
        [y_left, z_bottom],
        [y_web_left, z_bottom],
        [y_web_left, z_flange_bottom],
        [y_left, z_flange_bottom]
    ], closed=True))
    
    # Bottom flange - right outer part
    steel_patches.append(Polygon([
        [y_web_right, z_bottom],
        [y_right, z_bottom],
        [y_right, z_flange_bottom],
        [y_web_right, z_flange_bottom]
    ], closed=True))
    
    # Web - left solid part (if web thickness allows)
    if tw > 2*t_hollow:
        steel_patches.append(Polygon([
            [y_web_left, z_flange_bottom],
            [y_web_cavity_left, z_flange_bottom],
            [y_web_cavity_left, z_flange_top],
            [y_web_left, z_flange_top]
        ], closed=True))
        
        # Web - right solid part
        steel_patches.append(Polygon([
            [y_web_cavity_right, z_flange_bottom],
            [y_web_right, z_flange_bottom],
            [y_web_right, z_flange_top],
            [y_web_cavity_right, z_flange_top]
        ], closed=True))
    
    # Top flange - inner left part (between web and cavity)
    if B > 2*t_hollow + tw:
        steel_patches.append(Polygon([
            [y_web_left, z_flange_cavity_bottom],
            [y_flange_cavity_left, z_flange_cavity_bottom],
            [y_flange_cavity_left, z_flange_top],
            [y_web_left, z_flange_top]
        ], closed=True))
        
        # Top flange - inner right part
        steel_patches.append(Polygon([
            [y_flange_cavity_right, z_flange_cavity_bottom],
            [y_web_right, z_flange_cavity_bottom],
            [y_web_right, z_flange_top],
            [y_flange_cavity_right, z_flange_top]
        ], closed=True))
    
    # Bottom flange - inner left part
    if B > 2*t_hollow + tw:
        steel_patches.append(Polygon([
            [y_web_left, z_flange_bottom],
            [y_flange_cavity_left, z_flange_bottom],
            [y_flange_cavity_left, z_flange_cavity_top],
            [y_web_left, z_flange_cavity_top]
        ], closed=True))
        
        # Bottom flange - inner right part
        steel_patches.append(Polygon([
            [y_flange_cavity_right, z_flange_bottom],
            [y_web_right, z_flange_bottom],
            [y_web_right, z_flange_cavity_top],
            [y_flange_cavity_right, z_flange_cavity_top]
        ], closed=True))
    
    # -------------------------
    # Plot section geometry (improved)
    # -------------------------
    if plot:
        fig, ax = plt.subplots(figsize=(12, 10))
        
        # Draw steel patches with a light blue fill and subtle hatch
        steel_collection = PatchCollection(steel_patches, 
                                           facecolor='#add8e6',  # light blue
                                           edgecolor='#2f4f4f',   # dark slate
                                           linewidth=0.8,
                                           hatch='...',           # small dot pattern
                                           alpha=0.9)
        ax.add_collection(steel_collection)
        
        # Draw hollow cavities as white areas with red dashed borders
        # Top flange cavity
        if B > 2*t_hollow and tf > 2*t_hollow:
            top_cavity = Polygon([
                [y_flange_cavity_left, z_flange_cavity_bottom],
                [y_flange_cavity_right, z_flange_cavity_bottom],
                [y_flange_cavity_right, z_cavity_top],
                [y_flange_cavity_left, z_cavity_top]
            ], closed=True, facecolor='white', edgecolor='red', 
               linewidth=2, linestyle='--', alpha=0.8)
            ax.add_patch(top_cavity)
        
        # Bottom flange cavity
        if B > 2*t_hollow and tf > 2*t_hollow:
            bottom_cavity = Polygon([
                [y_flange_cavity_left, z_cavity_bottom],
                [y_flange_cavity_right, z_cavity_bottom],
                [y_flange_cavity_right, z_flange_cavity_top],
                [y_flange_cavity_left, z_flange_cavity_top]
            ], closed=True, facecolor='white', edgecolor='red', 
               linewidth=2, linestyle='--', alpha=0.8)
            ax.add_patch(bottom_cavity)
        
        # Web cavity
        if tw > 2*t_hollow and h_web > 2*t_hollow:
            web_cavity = Polygon([
                [y_web_cavity_left, z_web_cavity_bottom],
                [y_web_cavity_right, z_web_cavity_bottom],
                [y_web_cavity_right, z_web_cavity_top],
                [y_web_cavity_left, z_web_cavity_top]
            ], closed=True, facecolor='white', edgecolor='red', 
               linewidth=2, linestyle='--', alpha=0.8)
            ax.add_patch(web_cavity)
        
        # Centroid marker
        ax.plot(0, 0, 'r+', markersize=15, markeredgewidth=3, 
                label='Centroid', color='red', zorder=5)
        
        # Dimension arrows (cleaner layout)
        # Depth dimension (right side)
        ax.annotate('', xy=(y_right+10, z_bottom), xytext=(y_right+10, z_top),
                   arrowprops=dict(arrowstyle='<->', color='blue', lw=1.5))
        ax.text(y_right+20, 0, f'D = {D} mm', ha='left', va='center', 
               color='blue', fontsize=10, rotation=90,
               bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7))
        
        # Width dimension (top side)
        ax.annotate('', xy=(y_left, z_top+10), xytext=(y_right, z_top+10),
                   arrowprops=dict(arrowstyle='<->', color='blue', lw=1.5))
        ax.text(0, z_top+20, f'B = {B} mm', ha='center', va='bottom', 
               color='blue', fontsize=10,
               bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7))
        
        # Flange thickness (left side)
        ax.annotate('', xy=(y_left-15, z_flange_top), xytext=(y_left-15, z_top),
                   arrowprops=dict(arrowstyle='<->', color='green', lw=1.5))
        ax.text(y_left-30, (z_flange_top+z_top)/2, f'tf = {tf} mm', 
               ha='right', va='center', color='green', fontsize=9,
               bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7))
        
        # Web thickness (bottom side)
        ax.annotate('', xy=(y_web_left, z_bottom-15), xytext=(y_web_right, z_bottom-15),
                   arrowprops=dict(arrowstyle='<->', color='green', lw=1.5))
        ax.text(0, z_bottom-30, f'tw = {tw} mm', ha='center', va='top', 
               color='green', fontsize=9,
               bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7))
        
        # Hollow wall thickness (bottom left corner)
        ax.annotate('', xy=(y_left, z_bottom+tf/2), xytext=(y_left+t_hollow, z_bottom+tf/2),
                   arrowprops=dict(arrowstyle='<->', color='purple', lw=1.5))
        ax.text(y_left+t_hollow/2, z_bottom+tf/2+12, f't_hollow = {t_hollow} mm', 
               ha='center', va='bottom', color='purple', fontsize=9,
               bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7))
        
        # Aesthetics
        ax.set_aspect('equal')
        ax.set_xlabel('y (mm)', fontsize=12)
        ax.set_ylabel('z (mm)', fontsize=12)
        ax.set_title(f'Hollow I-Section\n' +
                    f'D={D}mm, B={B}mm, tf={tf}mm, tw={tw}mm, t_hollow={t_hollow}mm', 
                    fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.2, linestyle=':')
        
        # Custom legend
        legend_elements = [
            Patch(facecolor='#add8e6', edgecolor='#2f4f4f', hatch='...', label='Steel material'),
            Patch(facecolor='white', edgecolor='red', linestyle='--', label='Hollow cavity'),
            plt.Line2D([0], [0], marker='+', color='red', markersize=15, 
                      linewidth=0, label='Centroid')
        ]
        ax.legend(handles=legend_elements, loc='upper right', fontsize=10,
                 framealpha=0.9)
        
        # Set margins
        margin = max(D, B) * 0.15
        ax.set_xlim([y_left - margin, y_right + margin])
        ax.set_ylim([z_bottom - margin, z_top + margin])
        
        # Section properties box
        props_text = (f'Section Properties:\n'
                      f'Area = {AREA:.0f} mm²\n'
                      f'Mass/length = {DENSITY_STEEL * AREA:.3f} kg/mm\n'
                      f'Web height = {h_web:.1f} mm\n'
                      f'Cavity area = {hollow_area:.0f} mm²')
        ax.text(0.02, 0.98, props_text, transform=ax.transAxes,
               verticalalignment='top', fontsize=10,
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.9,
                        edgecolor='orange', linewidth=2))
        
        plt.tight_layout()
        plt.show()
    
    ELE_MASS = DENSITY_STEEL * AREA
    return AREA, ELE_MASS