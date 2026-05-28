def PIPE_HOLLOW_STEEL_SECTION_FUN_3D(SEC_TAG, MAT_TYPE, Do, t, numSubdivCirc, numSubdivRad, DENSITY_STEEL, plot=True):
    """
    Fiber section definition for a Steel Pipe Section (Circular Hollow Section)
    and optional visualization of the section geometry.
    Do     # [mm] outer width (along y and z)
    t      # [mm] wall thickness
    THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
    """  
    import openseespy.opensees as ops
    import numpy as np
    MAT_TAG = SEC_TAG + 1000
    if MAT_TYPE == 'ELASTIC':
        E = 200000.0  # [N/mm²] Modulus of steel
        #ops.uniaxialMaterial('Elastic', MAT_TAG, E)             # TESNSION AND COMPRESSION IS SAME VALUES
        ops.uniaxialMaterial('Elastic', MAT_TAG, E ,0.0, 0.5*E) # TESNSION AND COMPRESSION IS NOT SAME VALUES
        # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/ElasticUni.html
    # Material (inelastic steel)    
    if MAT_TYPE == 'INELASTIC':
        Fy = 240.0			    # [N/mm²] Steel yield stress
        Es = 200000.0		    # [N/mm²] Modulus of steel
        ey = Fy/Es			    # [mm/mm] Steel yield strain
        Fu = 1.1818*Fy          # [N/mm²] Steel Ultimate Strength
        esu = 0.12              # [mm/mm] Steel Ultimate Strain
        Esh = (Fu - Fy)/(esu - ey)
        Bs = Esh / Es           # strain-hardening ratio 
        pinchX = 0.8            # Pinching factor in X direction
        pinchY = 0.5            # Pinching factor in Y direction
        damage1 = 0.0           # Damage due to ductility
        damage2 = 0.0           # Damage due to energy
        beta = 0.1              # Stiffness degradation parameter
        ops.uniaxialMaterial('Hysteretic', MAT_TAG,
                             Fy, ey,
                             Fu, esu,
                             0.2*Fu, 1.1*esu,
                             -Fy, -ey,
                             -Fu, -esu,
                             -0.2*Fu, -1.1*esu,
                             pinchX, pinchY,
                             damage1, damage2, beta)
        # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/Hysteretic.html
    # ------------------------------
    # Pipe Geometry
    # ------------------------------
    #Do = 260.0           # Outer diameter (mm)
    #t = 10.0             # Wall thickness (mm)
    Di = Do - 2*t        # Inner diameter (mm)
    r_out = Do / 2
    r_in = Di / 2
    
    # ------------------------------
    # Fiber Section Definition
    # ------------------------------
    #numSubdivCirc = 32   # circumferential divisions
    #numSubdivRad = 4     # radial divisions
    
    ops.section('Fiber', SEC_TAG, '-GJ', 1.0e6)
    
    # Fiber ring of steel material
    ops.patch('circ', MAT_TAG, numSubdivRad, numSubdivCirc, 0.0, 0.0, r_in, r_out)
    
    # ------------------------------
    # Basic check (optional geometry)
    # ------------------------------
    AREA = np.pi * (r_out**2 - r_in**2)
    Iz = (np.pi / 4) * (r_out**4 - r_in**4)
    print(f"Section area A = {AREA:.6f} mm²")
    print(f"Moment of inertia Iz = {Iz:.8f} mm⁴")
    
    # -------------------------
    # Plot section geometry
    # -------------------------
    if plot:
        import matplotlib.pyplot as plt
        
        # Create angles for plotting circles
        theta = np.linspace(0, 2*np.pi, 100)
        
        # Calculate x and y coordinates for the outer circle
        x_out = r_out * np.cos(theta)
        y_out = r_out * np.sin(theta)
        
        # Calculate x and y coordinates for the inner circle
        x_in = r_in * np.cos(theta)
        y_in = r_in * np.sin(theta)
        
        # Plotting
        plt.figure(figsize=(6, 6))
        plt.plot(x_out, y_out, color='blue', label='Outer Wall')
        plt.plot(x_in, y_in, color='red', label='Inner Wall')
        
        # Fill the hollow part with white to make it look hollow
        plt.fill_between(x_out, y_out, color='grey', alpha=0.5)
        plt.fill_between(x_in, y_in, color='white') # Fill the inner circle with white
        
        plt.gca().set_aspect('equal', adjustable='box')
        plt.title('Steel Pipe Hollow Section')
        plt.xlabel('X-axis (mm)')
        plt.ylabel('Y-axis (mm)')
        plt.legend()
        plt.grid(True)
        plt.show()

    ELE_MASS  = DENSITY_STEEL * AREA   # Mass Per Length
     
    return Do, ELE_MASS   
