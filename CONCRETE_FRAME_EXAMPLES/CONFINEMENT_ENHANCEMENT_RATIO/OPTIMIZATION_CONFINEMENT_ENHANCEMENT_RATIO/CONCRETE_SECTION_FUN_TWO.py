
"""
def CONFINED_CONCRETE_SECTION(
    secTag: int, h: float, b: float, cover: float, As: float,
    coreTag: int, coverTag: int, steelTag: int,
    STEEL_KIND: str, Es: float, fy: float, ey: float, fu: float, esu: float,
    fcC: float, ec0C: float, fcUC: float, ecuC: float,
    fcU: float, ec0U: float, fcUU: float, ecuU: float,
    COL: bool = True
):
"""
def CONFINED_CONCRETE_SECTION(secTag, h, b, cover, As, STEEL_KIND, fc, Kc, COL=True, PLOT_STRESS=True):
    
    
    # Define parameters (units: mm, N)
    # ------------------------------------------
    # CONCRETE                  tag   f'c        ec0   f'cu        ecu    
    # Parametric definitions for unconfined concrete
    fcU = -fc                                 # [N/mm²] Unconfined concrete compressive strength
    ec0U  = -0.002 * (abs(fcU)/30)**0.25      # [mm/mm] Initial strain at peak strength (semi-empirical)
    fcUU  = 0.1 * fcU                         # [N/mm²] Ultimate stress (~10% of peak compressive strength)
    ecuU  = ec0U * 3.5                        # [mm/mm] Ultimate strain (e.g., 3.5× ec0U)
    
    # Parametric definitions for confined core concrete
    fcC   = fcU * Kc                          # [N/mm²] Confined strength (increased by confinement factor Kc)
    ec0C  = ec0U * 1.8                        # [mm/mm] Peak strain increases with confinement
    fcUC  = 0.95 * fcC                        # [N/mm²] Ultimate confined stress (~95% of peak)
    ecuC  = 1.475 * ecuU * Kc                 # [mm/mm] Ultimate confined strain (larger ductility due to confinement)
    """
    # Exponents chosen so that at Kc=2: ec0C/ec0U=1.8 and ecuC/ecuU=1.875
    alpha = math.log(1.8)  / math.log(2)                 # ≃0.847
    beta  = math.log(1.875)/ math.log(2)                 # ≃0.907
    
    ec0C = ec0U * Kc**alpha                              # [mm/mm] confined peak strain
    fcUC = 0.95 * fcC                                    # [N/mm²] confined ultimate stress
    ecuC = ecuU * Kc**beta                               # [mm/mm] confined ultimate strain    
    """
    if PLOT_STRESS == True:
        import numpy as np
        import matplotlib.pyplot as plt
        # Strain axis from zero to max of ultimate strains
        min_strain = min(ecuU, ecuC)
        strain = np.linspace(0, min_strain, 600)
    
        # Kent–Park stress–strain function
        def kent_park(stress_peak, strain_peak, stress_ult, strain_ult, epsilon):
            # ascending branch (parabolic)
            sigma = np.where(
                epsilon <= strain_peak,
                stress_peak * (2 * (epsilon/strain_peak) - (epsilon/strain_peak)**2),
                # descending branch (linear to ultimate)
                stress_peak + (stress_ult - stress_peak) * (epsilon - strain_peak) / (strain_ult - strain_peak)
            )
            return sigma
        
        # Compute for unconfined and confined
        stress_U = kent_park(fcU, ec0U, fcUU, ecuU, strain)
        stress_C = kent_park(fcC, ec0C, fcUC, ecuC, strain)
        
        # Plot
        plt.figure()
        plt.plot(strain, stress_U, label="Unconfined (Kent–Park)")
        plt.plot(strain, stress_C, label="Confined (Kent–Park)")
        plt.xlabel("Strain (mm/mm)")
        plt.ylabel("Stress (N/mm²)")
        plt.title("Stress–Strain Curves (Kent–Park Model)")
        plt.legend()
        plt.grid(True)
        plt.show()
    
    # STEEL
    # Reinforcing steel
    fy = 4000                                 # [N/mm²] Steel Rebar Yield Strength   
    Es = 2e5                                  # [N/mm²] Modulus of Elasticity
    ey = fy/Es                                # [mm/mm] Steel Rebar Yield Strain
    fu = 1.1818*fy                            # [N/mm²] Steel Rebar Ultimate Strength
    esu = ey*75.2                             # [mm/mm] Steel Rebar Ultimate Strain
    Esh = (fu - fy)/(esu - ey)
    Bs = Esh / Es
    
    import openseespy.opensees as ops
    
    coreTag, coverTag, steelTag = secTag + 100, secTag + 200, secTag + 300
    if STEEL_KIND == 1:# WITHOUT HARDENING AND ULTIMATE STRAIN
        ops.uniaxialMaterial('Steel01', steelTag, fy, Es, 0.0) 
    if STEEL_KIND == 2:# WITH HARDENING AND ULTIMATE STRAIN    
        pinchX = 0.8   # Pinching factor in X direction
        pinchY = 0.5   # Pinching factor in Y direction
        damage1 = 0.0  # Damage due to ductility
        damage2 = 0.0  # Damage due to energy
        beta = 0.1     # Stiffness degradation parameter
        ops.uniaxialMaterial('Hysteretic', steelTag, fy, ey, fu, esu, 0.2*fu, 1.1*esu, -fy, -ey, -fu, -esu, -0.2*fu, -1.1*esu, pinchX, pinchY, damage1, damage2, beta)
        # INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Hysteretic_Material
        
    ops.uniaxialMaterial('Concrete01', coreTag, fcC, ec0C, fcUC, ecuC)  # Core concrete (confined)
    ops.uniaxialMaterial('Concrete01', coverTag, fcU, ec0U, fcUU, ecuU) # Cover concrete (unconfined)
    
    # Some variables derived from the parameters
    y1 = h / 2.0
    z1 = b / 2.0
    NUMFIBERS = 40  # Number of layers for each fiber
    
    ops.section('Fiber', secTag)
    # Create the concrete core fibers
    ops.patch('rect', coreTag, NUMFIBERS, 5, cover - y1, cover - z1, y1 - cover, z1 - cover)
    
    # Create the concrete cover fibers (top, bottom, left, right)
    ops.patch('rect', coverTag, NUMFIBERS, 5, -y1, z1 - cover, y1, z1)
    ops.patch('rect', coverTag, NUMFIBERS, 5, -y1, -z1, y1, cover - z1)
    ops.patch('rect', coverTag, NUMFIBERS, 5, -y1, cover - z1, cover - y1, z1 - cover)
    ops.patch('rect', coverTag, NUMFIBERS, 5, y1 - cover, cover - z1, y1, z1 - cover)
    
    # Create the reinforcing fibers (left, middle, right)
    ops.layer('straight', steelTag, 3, As, y1 - cover, z1 - cover, y1 - cover, cover - z1)
    if COL == False:
        ops.layer('straight', steelTag, 2, As, 0.0, z1 - cover, 0.0, cover - z1)
    ops.layer('straight', steelTag, 3, As, cover - y1, z1 - cover, cover - y1, cover - z1)
