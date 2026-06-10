def SECTION_CURVATURE(eleTag, secTag):
    import openseespy.opensees as ops
    # WITH THIS FUNCTION IN EACH ANALYSIS STEPS, WE CAN EVALUATE SECTION STIFFNESS
    # WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
    # INFO LINK: https://opensees.berkeley.edu/OpenSees/manuals/usermanual/259.htm
    # recorder Element -file ele1sec1Stiff.out –time -ele 1 section 1 stiffness
    response = ops.eleResponse(eleTag, 'section', secTag, 'stiffness')
    # Get response and separate
    stiffness = response[0]
    return stiffness