def CONFINED_CONCRETE_SECTION_THERMAL(secTag, h, b, cover, As, coreTag, coverTag, steelTag, COL=True):
    import openseespy.opensees as ops

    # For plotting
    fiber_data = []

    NUMFIBERS = 20
    NTHICK = 5

    y_half = h / 2.0
    z_half = b / 2.0

    ops.section('FiberThermal', secTag)

    # Core concrete fibers
    y_core_min = cover - y_half
    y_core_max = y_half - cover
    z_core_min = cover - z_half
    z_core_max = z_half - cover
    dy_core = (y_core_max - y_core_min) / NUMFIBERS
    dz_core = (z_core_max - z_core_min) / NTHICK

    for i in range(NUMFIBERS):
        y_i = y_core_min + (i + 0.5) * dy_core
        for j in range(NTHICK):
            z_j = z_core_min + (j + 0.5) * dz_core
            area_ij = dy_core * dz_core
            ops.fiber(y_i, z_j, area_ij, coreTag)
            fiber_data.append((y_i, z_j, area_ij, coreTag))

    # Cover strips (top, bottom, left, right)
    def add_cover_fibers(y_min, y_max, z_min, z_max):
        dy = (y_max - y_min) / NUMFIBERS
        dz = (z_max - z_min) / NTHICK
        for i in range(NUMFIBERS):
            y_i = y_min + (i + 0.5) * dy
            for j in range(NTHICK):
                z_j = z_min + (j + 0.5) * dz
                area_ij = dy * dz
                ops.fiber(y_i, z_j, area_ij, coverTag)
                fiber_data.append((y_i, z_j, area_ij, coverTag))

    # Rebar layers
    z_extent = z_half - cover
    def _bar_z_positions(nBars):
        if nBars == 1:
            return [0.0]
        dz_bar = 2.0 * z_extent / (nBars - 1)
        return [-z_extent + k * dz_bar for k in range(nBars)]

    for y_bar, nBars in [(+y_half - cover, 3), (-y_half + cover, 3)]:
        for z_bar in _bar_z_positions(nBars):
            ops.fiber(y_bar, z_bar, As, steelTag)
            fiber_data.append((y_bar, z_bar, As, steelTag))

    if not COL:
        for z_bar in _bar_z_positions(2):
            ops.fiber(0.0, z_bar, As, steelTag)
            fiber_data.append((0.0, z_bar, As, steelTag))

    return fiber_data  # return this for plotting
