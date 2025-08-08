def CONFINED_CONCRETE_WITH_PLATE_SECTION_FUN(secTag, h, b, cover, As, coreTag, coverTag, steelTag, plateTag, T, PLATE=True, COL=True):
    import openseespy.opensees as ops

    # Derived parameters
    y1 = h / 2.0
    z1 = b / 2.0
    NUMFIBERS = 20

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
