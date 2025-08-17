def ZERO_LENGTH_ELEMENT(ele_tag, NODE_I, NODE_J, KP, KN, M, ELE_KIND, DRD, DOF):
    import openseespy.opensees as ops
    import numpy as np

    #%%
    Kd = KP[0] / KP[1]        # [N/mm] Elastic Stiffness
    omega = np.sqrt(Kd/M)     # [rad/sec] Viscous Damper Natural Frequency
    alpha = 0.7               # velocity exponent (usually 0.3–1.0)
    Cd = 2 * DRD * omega * M  # [N·s/mm] Damping coefficient 
    # Elements Tag
    mat_tag01 = ele_tag + 100; mat_tag02 = ele_tag + 200;
    
    if ELE_KIND == 'ELASTIC':
        ops.uniaxialMaterial('Elastic', mat_tag01, Kd, Cd)
        ops.element("zeroLength", ele_tag, NODE_I, NODE_J, '-mat', mat_tag01, '-dir', DOF)
    elif ELE_KIND == 'INELASTIC':
        ops.uniaxialMaterial('HystereticSM', mat_tag01, '-posEnv', *KP.flatten(), '-negEnv', *KN.flatten(), '-pinch', 1, 1)
        #INFO LINK: https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/HystereticSM.html
        ops.uniaxialMaterial('Viscous', mat_tag02, Cd, alpha)  # Material for C (alpha=1.0 for linear)
        #INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/Fatigue.html
        ops.element('zeroLength', ele_tag, NODE_I, NODE_J, '-mat', mat_tag01, mat_tag02, '-dir', DOF, DOF)