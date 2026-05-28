def THERMAL_ELEMENT_FUN(ELE_TAG, NODE_I, NODE_J, SEC_TAG, ELE_TYPE, distributed_load, Max_Thermal, DEPTH):
    # THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
    # CREATE ELEMENTS
    # Define beam-column elements
    # Thermo-mechaical beam-column Elements
    # INFO LINK: https://openseesforfire.github.io/Subpages/Elecmds.html
    import openseespy.opensees as ops
    transfTag = 1000 + ELE_TAG
    #ops.geomTransf('Linear', transfTag)
    #ops.geomTransf('PDelta', transfTag)
    ops.geomTransf('Corotational', transfTag)
    
    # Define beam integration (Lobatto integration)
    numIntegrationPoints = 5
    biTag01 = 2000 + ELE_TAG
    ops.beamIntegration('Lobatto', biTag01, SEC_TAG, numIntegrationPoints)
    #ops.element('forceBeamColumnThermal', ELE_TAG, NODE_I, NODE_J, transfTag, biTag01)
    ops.element('dispBeamColumnThermal', ELE_TAG, NODE_I, NODE_J, transfTag, biTag01)
    # eleLoad -ele $eleTag -type -beamThermal $T1 $y1 $T2 $Y2
    DD = 0.5 * DEPTH
    ops.eleLoad('-ele', ELE_TAG, '-type', '-beamThermal', Max_Thermal, -DD, Max_Thermal, DD) 
    