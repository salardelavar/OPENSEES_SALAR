def CONFINED_CONCRETE_WITH_PLATE_SECTION_FUN(secTag, HSec, BSec, cover, As,
                                             coreTag, coverTag, steelTag, plateTag, T,
                                             PLATE=True, COL=True, PLOT=True,
                                             CONCRETE_DENSITY=2500/1e9, STEEL_DENSITY=7850/1e9):
    import openseespy.opensees as ops

    # Derived parameters
    y1 = HSec / 2.0
    z1 = BSec / 2.0
    NUMFIBERS = 20
    
    # Calculate mass
    MASS = HSec * BSec * CONCRETE_DENSITY + ((T*(HSec+2*T)) + (T*(BSec-2*T))) * STEEL_DENSITY
    
    ops.section('Fiber', secTag)

    # Core concrete
    ops.patch('rect', coreTag, NUMFIBERS, 5, cover - y1, cover - z1, y1 - cover, z1 - cover)

    # Cover concrete
    ops.patch('rect', coverTag, NUMFIBERS, 5, -y1, z1 - cover, y1, z1)
    ops.patch('rect', coverTag, NUMFIBERS, 5, -y1, -z1, y1, cover - z1)
    ops.patch('rect', coverTag, NUMFIBERS, 5, -y1, cover - z1, cover - y1, z1 - cover)
    ops.patch('rect', coverTag, NUMFIBERS, 5, y1 - cover, cover - z1, y1, z1 - cover)

    # Longitudinal steel (3 bars on each face, optional middle)
    ops.layer('straight', steelTag, 3, As, y1 - cover, z1 - cover, y1 - cover, cover - z1)
    if not COL:
        ops.layer('straight', steelTag, 2, As, 0.0, z1 - cover, 0.0, cover - z1)
    ops.layer('straight', steelTag, 3, As, cover - y1, z1 - cover, cover - y1, cover - z1)
    
    if PLATE:
        if COL: # COLUMN SECTION 
            # Steel Plate Jacket (Top, Bottom, Left, Right)
            ops.patch('rect', plateTag, NUMFIBERS, 2, -(y1 + T), z1, y1 + T, z1 + T)    # Top plate
            ops.patch('rect', plateTag, NUMFIBERS, 2, -(y1 + T), -z1 - T, y1 + T, -z1)  # Bottom plate
            ops.patch('rect', plateTag, 2, NUMFIBERS, -(y1 + T), -z1, -y1, z1)          # Left plate
            ops.patch('rect', plateTag, 2, NUMFIBERS, y1, -z1, y1 + T, z1)              # Right plate
        
        if not COL: # BEAM SECTION 
            # Steel Plate Jacket (Top, Bottom, Left, Right)
            ops.patch('rect', plateTag, NUMFIBERS, 2, -(y1 + T), -z1 - T, y1 + T, -z1)  # Bottom plate
   
    if PLOT:
        import matplotlib.pyplot as plt
        import matplotlib.patches as patches
        fig, ax = plt.subplots(figsize=(7, 7))

        # Core concrete
        core = patches.Rectangle(
            (cover - y1, cover - z1),
            (y1 - cover) - (cover - y1),
            (z1 - cover) - (cover - z1),
            facecolor="lightgray"
        )
        ax.add_patch(core)
        
        # Cover concrete (4 sides)
        cover_regions = [
            (-y1, z1 - cover, 2*y1, cover),            # top
            (-y1, -z1, 2*y1, cover),                   # bottom
            (-y1, cover - z1, cover, (z1 - cover)*2),  # left
            (y1 - cover, cover - z1, cover, (z1 - cover)*2)  # right
        ]
        for (x, y, w, h_) in cover_regions:
            ax.add_patch(patches.Rectangle((x, y), w, h_, facecolor="lightyellow"))
        
        # Steel bars (approx only endpoints)
        bar_coords = [
            (y1 - cover, z1 - cover),
            (y1 - cover, cover - z1),
            (y1 - cover, (z1 - cover + cover - z1)/2)
        ]
        for (x, y) in bar_coords:
            circ = patches.Circle((x, y), radius=As**0.5, color="red")
            ax.add_patch(circ)
        
        # Plate jacket
        if PLATE:
            if COL:  # four sides
                plate_boxes = [
                    (-(y1+T),  z1,      (y1+T)-( -y1-T ),   T),      # top
                    (-(y1+T), -z1-T,    (y1+T)-( -y1-T ),   T),      # bottom
                    (-(y1+T), -z1,      T, 2*z1),                      # left
                    (y1,      -z1,      T, 2*z1)                       # right
                ]
            else:  # beam
                plate_boxes = [
                    (-(y1+T), -z1-T, (y1+T)-( -y1-T ), T)  # bottom only
                ]
            for (x, y, w, h_) in plate_boxes:
                ax.add_patch(patches.Rectangle((x, y), w, h_, facecolor="lightblue"))
        
        # Formatting
        ax.set_aspect("equal")
        ax.set_xlabel("y axis")
        ax.set_ylabel("z axis")
        ax.set_title("Fiber Section Visualization")
        ax.grid(True)
        plt.show()
    
    return HSec, MASS
          
