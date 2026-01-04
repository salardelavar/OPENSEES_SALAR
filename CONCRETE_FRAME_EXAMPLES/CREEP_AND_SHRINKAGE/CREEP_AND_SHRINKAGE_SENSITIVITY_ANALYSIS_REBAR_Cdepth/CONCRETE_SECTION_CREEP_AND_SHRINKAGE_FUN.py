def CONFINED_CONCRETE_CREEP_AND_SHRINKAGE_SECTION(secTag, h, b, cover, As, fc,
                                                  Kfc, STEEL_KIND, COL=True):
    # THIS FUNCTION WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
    import openseespy.opensees as ops
    import numpy as np
    # unconfined concrete
    Ec = 4700 * np.sqrt(-fc) # [N/mm^2] Concrete Elastic Modulus
    fc1U = fc;			    # UNCONFINED concrete maximum stress
    eps1U = -0.0025;	    # strain at maximum strength of unconfined concrete
    fc2U = 0.2*fc1U;		# ultimate stress
    eps2U = -0.012;			# strain at ultimate stress
    Lambda = 0.1;		    # ratio between unloading slope at $eps2 and initial slope $Ec
    EcU = 4700 * np.sqrt(-fc1U) # [N/mm^2] Concrete Elastic Modulus
    
    # confined concrete - bottom and top section
    Kfc = 1.2;			    # ratio of confined to unconfined concrete strength
    fc1C = Kfc*fc;		    # CONFINED concrete (mander model), maximum stress
    eps1C = 2*fc1C/Ec;	    # strain at maximum stress 
    fc2C = 0.2*fc1C;		# ultimate stress
    eps2C = 5*eps1C;		# strain at ultimate stress 
    EcC = 4700 * np.sqrt(-fc1C) # [N/mm^2] Concrete Elastic Modulus

    # tensile-strength properties
    ftC = -0.55*fc1C;		# tensile strength +tension
    ftU = -0.55*fc1U;		# tensile strength +tension
    EtC = ftC/0.002;		# tension softening stiffness
    EtU = ftU/0.002;		# tension softening stiffness
     
    # STEEL
    # Reinforcing steel
    fy = 4000           # [N/mm²] Steel Rebar Yield Strength   
    Es = 2e5            # [N/mm²] Modulus of Elasticity
    ey = fy/Es          # [mm/mm] Steel Rebar Yield Strain
    fu = 1.1818*fy      # [N/mm²] Steel Rebar Ultimate Strength
    esu = ey*75.2       # [mm/mm] Steel Rebar Ultimate Strain
    Esh = (fu - fy)/(esu - ey)
    Bs = Esh / Es
    
    # TDConcrete
    # LINK INFO: https://openseespydoc.readthedocs.io/en/latest/src/TDConcrete.html
    BETA = 0.4       # Tension softening parameter (tension softening exponent) (Recommended value)
    tD = 14          # Analysis time at initiation of drying (in days)
    epsshu = 400e-6  # Ultimate shrinkage strain as per ACI 209R-92 (shrinkage is negative)
    psish = 64.174   # Fitting parameter of the shrinkage time evolution function as per ACI 209R-92
    Tcr = 28         # Creep model age (in days)

    phiu = 3.0       # Ultimate creep coefficient as per ACI 209R-92
    psicr1 = 1.0     # Fitting parameter of the creep time evolution function as per ACI 209R-92 (Recommended value) 
    psicr2 = 75.4218 # Fitting parameter of the creep time evolution function as per ACI 209R-92 (Based on section dimensions)
    tcast = 2        # Analysis time corresponding to concrete casting (in days; minimum value 2.0)

    # TDConcreteEXP
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/TDConcreteEXP.html
    epscru = 0.0002  # Ultimate creep strain (e.g., taken from experimental measurements)
    sigCr = -5       # Concrete compressive stress (input as negative) associated with $epscru (e.g., experimentally applied)
    psicr1 = 1.0     # Fitting parameter of the creep time evolution function as per ACI 209R-92 (Recommended value) 
    psicr2 = 75.4218 # Fitting parameter of the creep time evolution function as per ACI 209R-92 (Based on section dimensions)
    tcast = 2        # Analysis time corresponding to concrete casting (in days; minimum value 2.0)
    
    # TDConcreteMC10
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/TDConcreteMC10.html

    # TDConcreteMC10NL
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/TDConcreteMC10NL.html
    
    
    coreTag = secTag+1   # Concrete Core Tag
    coverTag = secTag+2  # Concrete Cover Tag
    steelTag = secTag+3  # Concrete Steel Rebar Tag
    
    if STEEL_KIND == 1:# WITHOUT HARDENING AND ULTIMATE STRAIN
        ops.uniaxialMaterial('Steel01', steelTag, fy, Es, 0.0) 
        
    if STEEL_KIND == 2:# WITH HARDENING AND ULTIMATE STRAIN    
        pinchX = 0.8   # Pinching factor in X direction
        pinchY = 0.5   # Pinching factor in Y direction
        damage1 = 0.0  # Damage due to ductility
        damage2 = 0.0  # Damage due to energy
        beta = 0.1 # Stiffness degradation parameter
        ops.uniaxialMaterial('Hysteretic', steelTag, fy, ey, fu, esu, 0.2*fu, 1.1*esu, -fy, -ey, -fu, -esu, -0.2*fu, -1.1*esu, pinchX, pinchY, damage1, damage2, beta)
        # INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Hysteretic_Material

    # Time-dependent concrete material - Core concrete (confined)
    ops.uniaxialMaterial('TDConcrete', coreTag, fc1C, ftC, EcC, BETA, tD, -epsshu, psish, Tcr, phiu, psicr1, psicr2, tcast) 
    # Time-dependent concrete material - Cover concrete (unconfined)
    ops.uniaxialMaterial('TDConcrete', coverTag, fc1U, ftU, EcU, BETA, tD, -epsshu, psish, Tcr, phiu, psicr1, psicr2, tcast) 
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/TDConcrete.html
    
    #ops.uniaxialMaterial('TDConcreteEXP', coverTag, fc1U, ftU, EcU, BETA, tD, -epsshu, psish, Tcr, epscru, sigCr, psicr1, psicr2, tcast)
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/TDConcreteEXP.html
    
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
    
    return Tcr
