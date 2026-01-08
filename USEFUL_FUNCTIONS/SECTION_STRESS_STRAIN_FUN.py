def SECTION_STRESS_STRAIN(eleTag, secTag, FiberZ, FiberY):
    import openseespy.opensees as ops
    # WITH THIS FUNCTION IN EACH ANALYSIS STEPS, WE CAN EVALUATE SECTION STRESS-STRAIN
    # WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
    # INFO LINK: https://opensees.berkeley.edu/OpenSees/manuals/usermanual/259.htm
    # recorder Element -file ele1sec1StressStrain.out –time -ele 1 section 1 fiber $y $z <$matID> stressStrain
    response = ops.eleResponse(eleTag, 'section', secTag, 'fiber', FiberY, FiberZ, 'stressStrain')
    # Get response and separate
    stress = response[0]
    strain = response[1]
    return stress, strain