def TRUSS_ELEMENT(ele_tag, NODE_I, NODE_J, KP, KN, AREA, M, ELE_KIND, DRD, PLOT):
    import openseespy.opensees as ops
    import numpy as np

    Kd = KP[0] / KP[1]        # [N/mm] Elastic Stiffness
    omega = np.sqrt(Kd/M)     # [rad/sec] Viscous Damper Natural Frequency
    alpha = 0.7               # velocity exponent (usually 0.3–1.0)
    Cd = 2 * DRD * omega * M  # [N·s/mm] Damping coefficient 
    # Elements Tag
    mat_tag01 = ele_tag + 100; mat_tag02 = ele_tag + 200;
    
    if ELE_KIND == 'ELASTIC':
        ops.uniaxialMaterial('Elastic', mat_tag01, Kd, Cd)
        ops.element('corotTruss', ele_tag, NODE_I, NODE_J, AREA, mat_tag01)
    elif ELE_KIND == 'INELASTIC':
        ops.uniaxialMaterial('HystereticSM', mat_tag01, '-posEnv', *KP.flatten(), '-negEnv', *KN.flatten(), '-pinch', 1, 1)
        #INFO LINK: https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/HystereticSM.html
        ops.uniaxialMaterial('Viscous', mat_tag02, Cd, alpha)  # Material for C (alpha=1.0 for linear)
        #INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/Fatigue.html
        ops.element('corotTruss', ele_tag, NODE_I, NODE_J, AREA, mat_tag01)
        ops.element('corotTruss', ele_tag+1, NODE_I, NODE_J, AREA, mat_tag02)
        
    if PLOT == True:
        import matplotlib.pyplot as plt
        KPT = KP.reshape(-1, 2)  # Reshape into pairs of (stress, strain)
        KNC = KN.reshape(-1, 2)  # Reshape into pairs of (stress, strain)
        # Extract force and displacement values
        strain_tension = KPT[:, 1]
        stress_tension = KPT[:, 0]
        strain_compression = KNC[:, 1]
        stress_compression = KNC[:, 0]

        # Plot force-displacement relation
        plt.figure(0, figsize=(8, 6))
        plt.plot(strain_tension, stress_tension, 'bo-', label='Tension Strength')
        plt.plot(strain_compression, stress_compression, 'ro-', label='Compression Strength')
        plt.axhline(0, color='black', linewidth=1)
        plt.axvline(0, color='black', linewidth=1)
        plt.xlabel("Displacement [mm]")
        plt.ylabel("Force [N]")
        plt.title("Stress-Strain Nonlinear Relation for Viscous Damper")
        plt.legend()
        plt.grid(True)
        plt.show()
        