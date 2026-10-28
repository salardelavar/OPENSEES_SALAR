def SECTION_CURVATURE(eleTag, secTag):
    import openseespy.opensees as ops
    # WITH THIS FUNCTION IN EACH ANALYSIS STEPS, WE CAN EVALUATE SECTION FORCE
    # WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
    # INFO LINK: https://opensees.berkeley.edu/OpenSees/manuals/usermanual/259.htm
    # recorder Element -file ele1sec1Force.out –time -ele 1 section 1 force
    response = ops.eleResponse(eleTag, 'section', secTag, 'force')
    # Get response and separate
    force = response[0]
    return force