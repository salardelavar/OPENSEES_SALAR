def PIPE_HOLLOW_STEEL_SECTION_FUN(secTag, matTag, Do, t, numSubdivCirc, numSubdivRad, DENSITY_STEEL, plot=True):
    """
    Fiber section definition for a Steel Pipe Section (Circular Hollow Section)
    and optional visualization of the section geometry.
    Do     # [mm] outer width (along y and z)
    t      # [mm] wall thickness
    THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
    """  
    import openseespy.opensees as ops
    import numpy as np
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
    
    ops.section('Fiber', secTag)
    
    # Fiber ring of steel material
    ops.patch('circ', matTag, numSubdivRad, numSubdivCirc, 0.0, 0.0, r_in, r_out)
    
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
