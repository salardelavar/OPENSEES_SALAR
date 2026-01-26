def CONFINED_CONCRETE_SECTION(secTag, h, b, cover, RD, STEEL_KIND, CONCRETE_KIND, fc, Kc, COL=True, DENSITY=2500/1e9):
    import openseespy.opensees as ops
    import numpy as np
    
    # Mass per Length
    AREA = b * h
    MASS = DENSITY * AREA
    
    # Rebar Area
    As = (np.pi * RD**2)/4
    
    # Define parameters (units: mm, N)
    # ------------------------------------------
    # CONCRETE                  tag   f'c        ec0   f'cu        ecu    
    # Parametric definitions for unconfined concrete
    fcU = -fc                                 # [N/mm²] Unconfined concrete compressive strength
    EcU = 4700 * np.sqrt(-fcU)                # [N/mm^2] Unconfined Concrete Elastic Modulus
    ec0U  = 2*fcU/EcU                         # [mm/mm] Initial strain at peak strength (semi-empirical)
    fcUU  = 0.667 * fcU                       # [N/mm²] Ultimate stress (~10% of peak compressive strength)
    ecuU  = ec0U * 4.0                        # [mm/mm] Ultimate strain (e.g., 3.5× ec0U)
    lamdaU = 0.1                              # Ratio between unloading slope at epsU and initial slope
    ftU = 0.7 * np.sqrt(-fcU)                 # [N/mm²] Tensile strength
    EtsU = ftU/np.abs(ec0U)                   # [N/mm²] Tension softening stiffness 
    
    # Parametric definitions for confined core concrete
    fcC   = fcU * Kc                          # [N/mm²] Confined strength (increased by confinement factor Kc)
    EcC = 4700 * np.sqrt(-fcC)                # [N/mm^2] Confined Concrete Elastic Modulus
    ec0C  = 2*fcC/EcC                         # [mm/mm] Peak strain increases with confinement
    fcUC  = 0.6667 * fcC                      # [N/mm²] Ultimate confined stress (~95% of peak)
    ecuC  = 4.0 * ecuU                        # [mm/mm] Ultimate confined strain (larger ductility due to confinement)
    lamdaC = 0.1                              # Ratio between unloading slope at epsU and initial slope
    ftC = 0.7 * np.sqrt(-fcC)                 # [N/mm²] Tensile strength
    EtsC = ftC/np.abs(ec0C)                   # [N/mm²] Tension softening stiffness

    
    # STEEL
    # Reinforcing steel
    fy = 400                                  # [N/mm²] Steel Rebar Yield Strength   
    Es = 2e5                                  # [N/mm²] Modulus of Elasticity
    ey = fy/Es                                # [mm/mm] Steel Rebar Yield Strain
    fu = 405.8                                # [N/mm²] Steel Rebar Ultimate Strength
    esu = 0.09                                # [mm/mm] Steel Rebar Ultimate Strain
    Esh = (fu - fy)/(esu - ey)
    Bs = Esh / Es
    
    
    
    coreTag, coverTag, steelTag = secTag + 100, secTag + 200, secTag + 300
    if STEEL_KIND == 1:# WITHOUT HARDENING AND ULTIMATE STRAIN
        ops.uniaxialMaterial('Steel01', steelTag, fy, Es, Bs) 
    if STEEL_KIND == 2:# WITH HARDENING AND ULTIMATE STRAIN    
        pinchX = 0.8   # Pinching factor in X direction
        pinchY = 0.5   # Pinching factor in Y direction
        damage1 = 0.0  # Damage due to ductility
        damage2 = 0.0  # Damage due to energy
        beta = 0.1     # Stiffness degradation parameter
        ops.uniaxialMaterial('Hysteretic', steelTag, fy, ey, fu, esu, 0.2*fu, 1.1*esu, -fy, -ey, -fu, -esu, -0.2*fu, -1.1*esu, pinchX, pinchY, damage1, damage2, beta)
        # INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Hysteretic_Material
        
    if CONCRETE_KIND == 1:
        ops.uniaxialMaterial('Concrete01', coreTag, fcC, ec0C, fcUC, ecuC)  # Core concrete (confined)
        ops.uniaxialMaterial('Concrete01', coverTag, fcU, ec0U, fcUU, ecuU) # Cover concrete (unconfined)
    if CONCRETE_KIND == 2:
        ops.uniaxialMaterial('Concrete02', coreTag, fcC, ec0C, fcUC, ecuC, lamdaC, ftC, EtsC)  # Core concrete (confined)
        ops.uniaxialMaterial('Concrete02', coverTag, fcU, ec0U, fcUU, ecuU, lamdaU, ftU, EtsU) # Cover concrete (unconfined)
    
    # Some variables derived from the parameters
    y1 = h / 2.0
    z1 = b / 2.0
    NUMFIBERS = 40  # Number of layers for each fiber
    
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
    
    return h, MASS # Return Section Height
