def DAMAGE_INDEX_FUN(DISP, DY, DSU):
    # EVALUATION OF DUCTILITY DAMAGE INDEX
    import numpy as np
    DI = 100*(np.abs(DISP)-DY)/(DSU-DY) 
    if DI <= 0.0: # DAMAGE INDEX -> 0.0 <= DI <=100 
        DI = 0.0
    if DI >= 100: 
        DI = 100.0
    return  DI   