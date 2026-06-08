def SECTION_CURVATURE(eleTag, secTag):
    import openseespy.opensees as ops
    # WITH THIS FUNCTION IN EACH ANALYSIS STEPS, WE CAN EVALUATE SECTION CURVATURE
    # WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
    # INFO LINK: https://opensees.berkeley.edu/OpenSees/manuals/usermanual/259.htm
    # recorder Element -file ele1sec1Defo.out –time -ele 1 section 1 deformation
    response = ops.eleResponse(eleTag, 'section', secTag, 'deformation')
    # Get response and separate
    cur = response[0]
    return cur