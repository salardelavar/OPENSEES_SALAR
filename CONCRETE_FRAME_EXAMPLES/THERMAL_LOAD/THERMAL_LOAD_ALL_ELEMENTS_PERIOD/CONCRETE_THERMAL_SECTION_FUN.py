def CONFINED_CONCRETE_SECTION_THERMAL(secTag, fc, Kc, h, b, cover, As, COL=True):
    import openseespy.opensees as ops
    
    # Define parameters (units: mm, N)
    # ------------------------------------------
    # CONCRETE                  tag   f'c        ec0   f'cu        ecu

    # Cover concrete (unconfined)
    fcU = -fc                                # [N/mm²] Concrete Compressive Strength
    ec0U = -0.002 * (abs(fcU)/30)**0.25      # [mm/mm] Concrete Compressive Strain
    fcUU = 0.1 * fcU                         # [N/mm²] Concrete Compressive Ultimate Strength
    ecuU = ec0U * 3.5                        # [mm/mm] Concrete Compressive Ultimate Strain
    lamdaU = 0.1                             # Ratio between unloading slope at epsU and initial slope
    ftU = 0.05 * fcU                         # [N/mm²] Tensile strength
    EtsU = ftU/0.002                         # [N/mm²] Tension softening stiffness 
    
    # Core concrete (confined)
    fcC = fcU * Kc                # [N/mm²] Concrete Compressive Strength
    ec0C = ec0U * 1.8             # [mm/mm] Concrete Compressive Strain
    fcUC = 0.95 * fcC             # [N/mm²] Concrete Compressive Ultimate Strength
    ecuC = 1.475 * ecuU * Kc      # [mm/mm] Concrete Compressive Ultimate Strain
    lamdaC = 0.1                  # Ratio between unloading slope at epsU and initial slope
    ftC = 0.05 * fcC              # [N/mm²] Tensile strength
    EtsC = ftC/0.002              # [N/mm²] Tension softening stiffness

     
    # STEEL
    # Reinforcing steel
    fy = 400            # [N/mm²] Steel Rebar Yield Strength   
    Es = 2e5            # [N/mm²] Modulus of Elasticity
    ey = fy/Es          # [mm/mm] Steel Rebar Yield Strain
    fu = 1.1818*fy      # [N/mm²] Steel Rebar Ultimate Strength
    esu = 0.12          # [mm/mm] Steel Rebar Ultimate Strain
    Esh = (fu - fy)/(esu - ey)
    Bs = Esh / Es
    
    coreTag = 1 + secTag
    coverTag = 2 + secTag
    steelTag = 3 + secTag
    
    # Steel Rebar
    ops.uniaxialMaterial('Steel01Thermal', steelTag, fy, Es, Bs) 
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/steel01thermal.html
    # Concrete 
    
    """
    Concrete02Thermal:
    Concrete02Thermal is created for modelling concrete, which is derived from
    the standard "Concrete02" material and incorprates with temperature dependent
    properties sugggested in Eurocode 2 (EN1992-1-2)
    """
    ops.uniaxialMaterial('Concrete02Thermal', coreTag, fcC, ec0C, fcUC, ecuC, lamdaC, ftC, EtsC)  # Core concrete (confined)
    ops.uniaxialMaterial('Concrete02Thermal', coverTag, fcU, ec0U, fcUU, ecuU, lamdaU, ftU, EtsU) # Cover concrete (unconfined)
    # INFO LINK: https://openseesforfire.github.io/Subpages/MaterialCmds.html
    print(' Thermal Material Done.')

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
