def SECTION_CRACK_DEPTH_FUN(eleTag, secTag,
                        STRAINt,          # TOP FIBER MATERIAL ULTIMATE STRAIN FOR CRACK DEPTH EVALUTION
                        STRAINb,          # BOTTOM FIBER MATERIAL ULTIMATE STRAIN FOR CRACK DEPTH EVALUTION
                        FiberYb, FiberZb, # BOTTOM FIBER FROM NEUTRAL AXIS
                        FiberYt, FiberZt, # TOP FIBER FROM NEUTRAL AXIS 
                        ):
    import openseespy.opensees as ops
    import numpy as np
    # WITH THIS FUNCTION IN EACH ANALYSIS STEPS, WE CAN EVALUATE SECTION STRESS-STRAIN
    # WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
    # INFO LINK: https://opensees.berkeley.edu/OpenSees/manuals/usermanual/259.htm
    # recorder Element -file ele1sec1StressStrain.out –time -ele 1 section 1 fiber $y $z <$matID> stressStrain
    responseB = ops.eleResponse(eleTag, 'section', secTag, 'fiber', FiberYb, FiberZb, 'stressStrain')
    # Get response and separate
    stressB = responseB[0]
    strainB = responseB[1]
    responseT = ops.eleResponse(eleTag, 'section', secTag, 'fiber', FiberYt, FiberZt, 'stressStrain')
    # Get response and separate
    stressT = responseT[0]
    strainT = responseT[1]
    h = np.abs(FiberYb) + np.abs(FiberYt)
    X = (np.abs(strainT) * h) / (np.abs(strainT) + np.abs(strainB)) # NEUTRAL AXIS
    
    if strainT <= STRAINt:
        Z = (strainT / STRAINt)
        Ct = (Z * X - X) / Z      # CRACK DEPTH
    else:
        Ct = 0.0
    
    x = h - X    # DISTANCE FROM NEUTRAL AXIS TO BOTTOM FIBER
    
    if strainB >= STRAINb:
        Z = (strainB / STRAINb)
        Cb = (Z * x - x) / Z      # CRACK DEPTH
    else:
        Cb = 0.0        
    #print(X, Ct) 
    
    return X, x, Ct, Cb