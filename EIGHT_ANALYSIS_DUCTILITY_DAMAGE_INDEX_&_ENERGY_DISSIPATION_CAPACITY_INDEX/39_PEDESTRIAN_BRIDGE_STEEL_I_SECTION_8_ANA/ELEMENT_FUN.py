def ELEMENT_FUN(ELE_TAG, NODE_I, NODE_J, SEC_TAG, ELE_TYPE, Ele_Mass):
    # THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
    # CREATE ELEMENTS
    import openseespy.opensees as ops
    transfTag = 1000 + ELE_TAG
    #ops.geomTransf('Linear', transfTag)
    #ops.geomTransf('PDelta', transfTag)
    ops.geomTransf('Corotational', transfTag)
    numIntgrPts = 5
    if ELE_TYPE == 'elasticBeamColumn':
        ops.element('elasticBeamColumn', ELE_TAG, NODE_I, NODE_J, numIntgrPts, SEC_TAG, transfTag,'-mass', Ele_Mass)

    if ELE_TYPE == 'nonlinearBeamColumn':
        ops.element('nonlinearBeamColumn', ELE_TAG, NODE_I, NODE_J, numIntgrPts, SEC_TAG, transfTag,'-mass', Ele_Mass) 