def RECONSTRUCTION_COST_FUN(DISP, DY, DSU):
    # EVALUATION OF RECONSTRUCTION COST
    # THIS PYTHON SCRIPT WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
    """
    The `RECONSTRUCTION_COST_FUN` function calculates the damage index (DI) from the absolute
    displacement (`DISP`), yield displacement (`DY`), and ultimate displacement (`DSU`),
    limiting DI to a range of 0–100%. It classifies structural damage into 11 zones,
    from Zone 0 (No Damage) to Zone 10 (Failure Damage), and prints the corresponding damage category.
    Finally, it estimates the reconstruction cost (`rc`) as a percentage of a fixed $1,000,000
    reference cost, based on DI, and returns the cost in dollars.
    """
    import numpy as np
    DI = 100*(np.abs(DISP)-DY)/(DSU-DY) 
    if DI <= 0.0: # DAMAGE INDEX -> 0.0 <= DI <=100 
        DI = 0.0
    if DI >= 100: 
        DI = 100.0
        
    if DI <= 0:
        rc = 0.0
        print("\n\tZONE 0: NO DAMAGE\n")
    elif DI > 0 and DI <= 10:
        rc = int(1_000_000 * DI * 0.01) # [$]
        print("\n\tZONE 1: VERY MINOR DAMAGE\n")
    elif DI > 10 and DI <= 20:
        rc = int(1_000_000 * DI * 0.01)
        print("\n\tZONE 2: MINOR DAMAGE\n")
    elif DI > 20 and DI <= 30:
        rc = int(1_000_000 * DI * 0.01)
        print("\n\tZONE 3: MODERATE–LOW DAMAGE\n")
    elif DI > 30 and DI <= 40:
        rc = int(1_000_000 * DI * 0.01)
        print("\n\tZONE 4: MODERATE DAMAGE\n")
    elif DI > 40 and DI <= 50:
        rc = int(1_000_000 * DI * 0.01)
        print("\n\tZONE 5: MODERATE–HIGH DAMAGE\n")
    elif DI > 50 and DI <= 60:
        rc = int(1_000_000 * DI * 0.01)
        print("\n\tZONE 6: SEVERE–LOW DAMAGE\n")
    elif DI > 60 and DI <= 70:
        rc = int(1_000_000 * DI * 0.01)
        print("\n\tZONE 7: SEVERE–MEDIUM DAMAGE\n")    
    elif DI > 70 and DI <= 80:
        rc = int(1_000_000 * DI * 0.01)
        print("\n\tZONE 8: SEVERE–HIGH DAMAGE\n")
    elif DI > 80 and DI <= 90:
        rc = int(1_000_000 * DI * 0.01)
        print("\n\tZONE 9: VERY SEVERE DAMAGE\n")
    elif DI > 90 and DI <= 100:
        rc = int(1_000_000 * DI * 0.01)
        print("\n\tZONE 10: FAILURE DAMAGE\n")
            
      
    return rc   