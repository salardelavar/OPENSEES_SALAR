def CONFINED_CONCRETE_SECTION(secTag, h, b, cover, As, coreTag, coverTag, steelTag, COL=True):
    import openseespy.opensees as ops
    # Some variables derived from the parameters
    y1 = h / 2.0
    z1 = b / 2.0
    NUMFIBERS = 20  # Number of layers for each fiber
    
    ops.section('Fiber', secTag)
    # Create the concrete core fibers
    ops.patch('rect', coreTag, NUMFIBERS, 5, cover - y1, cover - z1, y1 - cover, z1 - cover)
    
    # Create the concrete cover fibers (top, bottom, left, right)
    ops.patch('rect', coverTag, NUMFIBERS, 5, -y1, z1 - cover, y1, z1)
    ops.patch('rect', coverTag, NUMFIBERS, 5, -y1, -z1, y1, cover - z1)
    ops.patch('rect', coverTag, NUMFIBERS, 5, -y1, cover - z1, cover - y1, z1 - cover)
    ops.patch('rect', coverTag, NUMFIBERS, 5, y1 - cover, cover - z1, y1, z1 - cover)
    
    # Create the reinforcing fibers (left, middle, right)
    ops.layer('straight', steelTag, 3, As, y1 - cover, z1 - cover, y1 - cover, cover - z1)
    if COL == False:
        ops.layer('straight', steelTag, 2, As, 0.0, z1 - cover, 0.0, cover - z1)
    ops.layer('straight', steelTag, 3, As, cover - y1, z1 - cover, cover - y1, cover - z1)
