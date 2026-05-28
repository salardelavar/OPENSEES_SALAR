def AGGREGATED_AXIAL_FLEXURAL_SECTION_SHEAR_WALL_ELEMENT_FUN(ELE_TAG,
                                                  NODE_I, NODE_J,
                                                  SECTION_Es, SECTION_Iz, SECTION_Area,
                                                  DENSITY_CONCRETE, STEEL_TYPE):
    # THIS PYTHON SCRIPT WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)    
    import openseespy.opensees as ops
    
    ColMatTagFlex = ELE_TAG + 1000
    ColMatTagAxial = ELE_TAG + 2000
    ColMatTagShear = ELE_TAG + 3000
    ColSecTag = ELE_TAG + 4000

    EI = SECTION_Es * SECTION_Iz
    EA = SECTION_Es * SECTION_Area
    
    if STEEL_TYPE == 'ELASTIC':
        # Define Nonlinear Sectional Moment-Curvature
        FiY = 0.000001474      # [1/mm] Yield Curvature of element
        MY = EI * FiY          # [N.mm] Yield Moment of element
        MU = 1.36 * MY         # [N.mm] Utimate Moment of element
        FiSU = 0.000541        # [1/mm] Utimate Curvature of element
        EIColCrack = MY / FiY  # [N.mm^2] Elasic Flextural Rigidity
        # MOMENT-CURVATURE RELATION FOR SECTION
        #ops.uniaxialMaterial('Elastic', ColMatTagFlex, EIColCrack)             # TESNSION AND COMPRESSION IS SAME VALUES
        ops.uniaxialMaterial('Elastic', ColMatTagFlex, EIColCrack ,0.0, 0.5*EIColCrack) # TESNSION AND COMPRESSION IS NOT SAME VALUES   
        # SHEAR-SHEAR STRAIN RELATION FOR SECTION
        # Define Nonlinear Sectional Shear-shear Strain
        ViY = 0.0009474         # [mm/mm] Yield Shear-strain of element
        VY = EA * ViY           # [N] Yield Shear of element
        VU = 1.1818 * VY        # [N] Utimate Shear of element
        ViSU = 0.0541           # [mm/mm] Utimate Shear-strain of element
        EIColShear = VY / ViY   # [N.mm^2] Elasic Shear Rigidity
        #ops.uniaxialMaterial('Elastic', ColMatTagShear, EIColShear)             # TESNSION AND COMPRESSION IS SAME VALUES
        ops.uniaxialMaterial('Elastic', ColMatTagShear, EIColShear ,0.0, 0.5*EIColShear) # TESNSION AND COMPRESSION IS NOT SAME VALUES
        
    if STEEL_TYPE == 'INELASTIC':
        # Define Nonlinear Sectional Moment-Curvature
        FiY = 0.000001474      # [1/mm] Yield Curvature of element
        MY = EI * FiY          # [N.mm] Yield Moment of element
        MU = 1.36 * MY         # [N.mm] Utimate Moment of element
        FiSU = 0.000541        # [1/mm] Utimate Curvature of element
        EIColCrack = MY / FiY  # [N.mm^2] Elasic Flextural Rigidity
        pinchX = 0.4           # Pinching factor in X direction
        pinchY = 0.2           # Pinching factor in Y direction
        damage1 = 0.0          # Damage due to ductility
        damage2 = 0.0          # Damage due to energy
        beta = 0.1             # Stiffness degradation 
    
        # MOMENT-CURVATURE RELATION FOR SECTION
        ops.uniaxialMaterial('Hysteretic', ColMatTagFlex,
                                            MY, FiY,
                                            MU, FiSU,
                                            0.23*MU, 1.13*FiSU,
                                            -MY, -FiY,
                                            -MU, -FiSU,
                                            -0.23*MU, -1.07*FiSU,
                                            pinchX, pinchY,
                                            damage1, damage2, beta)
        # INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Hysteretic_Material
        
        # Define Nonlinear Sectional Shear-shear Strain
        ViY = 0.0009474         # [mm/mm] Yield Shear-strain of element
        VY = EA * ViY           # [N] Yield Shear of element
        VU = 1.1818 * VY        # [N] Utimate Shear of element
        ViSU = 0.0541           # [mm/mm] Utimate Shear-strain of element
        EIColShear = VY / ViY   # [N.mm^2] Elasic Shear Rigidity
        pinchXs = 0.4           # Pinching factor in X direction
        pinchYs = 0.2           # Pinching factor in Y direction
        damage1s = 0.0          # Damage due to ductility
        damage2s = 0.0          # Damage due to energy
        betaS = 0.1             # Stiffness degradation 
    
        # SHEAR-SHEAR STRAIN RELATION FOR SECTION
        ops.uniaxialMaterial('Hysteretic', ColMatTagShear,
                                            VY, ViY,
                                            VU, ViSU,
                                            0.23*VU, 1.13*ViSU,
                                            -VY, -ViY,
                                            -VU, -ViSU,
                                            -0.23*VU, -1.07*ViSU,
                                            pinchXs, pinchYs,
                                            damage1s, damage2s, betaS)
        # INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Hysteretic_Material
    #------------------------------------------------------------------------------------------------
    ops.uniaxialMaterial("Elastic", ColMatTagAxial, EA)
    #------------------------------------------------------------------------------------------------
    ops.section("Aggregator", ColSecTag,
                ColMatTagAxial, "P",     # Axial force-deformation
                ColMatTagFlex, "Mz",     # Moment-curvature about section local z axis
                ColMatTagShear, "Vy")    # Shear force-deformation along section local y-axis
    """
    ops.section("Aggregator", ColSecTag,
                ColMatTagAxial, "P",     # Axial force-deformation
                ColMatTagFlex, "Mz")     # Moment-curvature about section local z axis
    """            
    #------------------------------------------------------------------------------------------------
    # Define Geometric Transformation
    ColTransfTag = ELE_TAG + 10
    ops.geomTransf("Linear", ColTransfTag)
    #ops.geomTransf('PDelta', ColTransfTag)
    #ops.geomTransf('Corotational', ColTransfTag)

    # Define Element
    ele_tag = ELE_TAG
    transfTag = ELE_TAG
    ops.geomTransf('Linear', transfTag)
    #ops.geomTransf('PDelta', transfTag)
    #ops.geomTransf('Corotational', transfTag)
    numIntgrPts = 5
    ele_Mass = SECTION_Area * DENSITY_CONCRETE
    ops.element("nonlinearBeamColumn", ELE_TAG, NODE_I, NODE_J, numIntgrPts, ColSecTag, ColTransfTag, '-mass', ele_Mass)