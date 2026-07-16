def OPENSEEES_HYSTERETICSM_FORCE_DISP_FUN(matTag, DP, FP, DN, FN,
                                            pinchX, pinchY,
                                            damage1, damage2, beta,
                                            PLOT = True, X_LABEL='DISP', Y_LABEL='FORCE', TITLE='FORCE-DISPLACEMENT CURVE'):
    import openseespy.opensees as ops
    import numpy as np
    import matplotlib.pyplot as plt
    # Positive branch points
    pos_disp = [0, DP[0], DP[1], DP[2], DP[3]]
    pos_force = [0, FP[0], FP[1], FP[2], FP[3]]
    KP = np.array([FP[0], DP[0], 
           FP[1], DP[1], 
           FP[2], DP[2],
           FP[3], DP[3]])
    
    # Negative branch points
    neg_disp = [0, DN[0], DN[1], DN[2], DN[3]]
    neg_force = [0, FN[0], FN[1], FN[2], FN[3]]
    KN = np.array([FN[0], DN[0], 
           FN[1], DN[1], 
           FN[2], DN[2],
           FN[3], DN[3]])
    
    if PLOT == True:
        # Plot
        plt.figure(0, figsize=(12, 10))
        plt.plot(pos_disp, pos_force, marker='o', color='red')
        plt.plot(neg_disp, neg_force, marker='o', color='black')
        
        plt.xlabel(X_LABEL)
        plt.ylabel(Y_LABEL)
        plt.title(TITLE)
        plt.grid(True)
        plt.axhline(0, linewidth=0.5)
        plt.axvline(0, linewidth=0.5)
        plt.show()
        
    ops.uniaxialMaterial('HystereticSM', matTag, *KP.flatten(), *KN.flatten(), pinchX, pinchY, damage1, damage2, beta)
    # INFO LINK: https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/HystereticSM.html