import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def PLOT_PUSHOVER(FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y, ROT, KA, KS, KI, STEP):
    plt.figure(1, figsize=(12, 8))
    plt.plot(MOMENT, FORCE_A, color='black')
    #plt.scatter(MOMENT, FORCE_A, color='black', linewidth=2)
    plt.title('P-M Interaction')
    plt.ylabel('Axial Force [N]')
    plt.xlabel('Bending Moment [N.mm]')
    plt.grid()
    plt.savefig(f'PUSHOVER_PLOT_{1}.png')
    plt.show()

    plt.figure(2, figsize=(12, 8))
    plt.plot(DISP_X, FORCE_S, color='green', linewidth=2)
    #plt.scatter(DISP_X, FORCE_S, color='green', linewidth=2)
    plt.title('SHEAR FORCE-DISPLACEMENT DIAGRAM')
    plt.ylabel('Shear Force [N]')
    plt.xlabel('Displacement  in X [mm]')
    plt.grid()
    plt.savefig(f'PUSHOVER_PLOT_{2}.png')
    plt.show()

    plt.figure(3, figsize=(12, 8))
    plt.plot(DISP_Y, FORCE_A, color='purple', linewidth=2)
    #plt.scatter(DISP_Y, FORCE_A, color='purple', linewidth=2)
    plt.title('AXIAL FORCE-DISPLACEMENT DIAGRAM')
    plt.ylabel('Axial Force [N]')
    plt.xlabel('Displacement in Y [mm]')
    plt.savefig(f'PUSHOVER_PLOT_{3}.png')
    plt.grid()
    plt.show()

    plt.figure(4, figsize=(12, 8))
    plt.plot(ROT, MOMENT, color='red', linewidth=2)
    #plt.scatter(ROT, MOMENT, color='red', linewidth=2)
    plt.title('MOMENT-ROTATION DIAGRAM')
    plt.ylabel('Moment [kN.mm]')
    plt.xlabel('Rotation [rad]')
    plt.savefig(f'PUSHOVER_PLOT_{4}.png')
    plt.grid()
    plt.show()

    plt.figure(5, figsize=(12, 8))
    #plt.plot(KI, KS, color='black', linewidth=2)
    plt.scatter(KI, KS, color='black', linewidth=2)
    plt.title('ROTATIONAL STIFFNESS-LATERAL STIFFNESS DIAGRAM')
    plt.ylabel('Rotational Stiffness [N.mm/Rad]')
    plt.xlabel('Lateral Stiffness in X Dir. [N/mm]')
    plt.semilogx()
    plt.semilogy()
    plt.grid()
    plt.savefig(f'PUSHOVER_PLOT_{5}.png')
    plt.show()

    plt.figure(6, figsize=(12, 8))
    #plt.plot(KI, KA, color='black', linewidth=2)
    plt.scatter(KI, KA, color='black', linewidth=2)
    plt.title('ROTATIONAL STIFFNESS-LATERAL STIFFNESS DIAGRAM')
    plt.ylabel('Rotational Stiffness [N.mm/Rad]')
    plt.xlabel('Lateral Stiffness in Y Dir. [N/mm]')
    plt.semilogx()
    plt.semilogy()
    plt.grid()
    plt.savefig(f'PUSHOVER_PLOT_{6}.png')
    plt.show()
    
    
    import BILINEAR_CURVE as BC
    
    # --------------------------------------
    #  Plot BaseShear-Displacement Analysis 
    # --------------------------------------
    XX = np.abs(DISP_X); YY = np.abs(FORCE_S); # ABSOLUTE VALUE
    SLOPE_NODE = 10
    
    DATA = BC.BILNEAR_CURVE(XX, YY, SLOPE_NODE)
    X, Y, Elastic_ST, Plastic_ST, Tangent_ST, Ductility_Rito, Over_Strength_Factor = DATA
    
    XLABEL = 'Displacement in X [mm]'
    YLABEL = 'Base-Shear Reaction [N]'
    LEGEND01 = 'Curve'
    LEGEND02 = 'Bilinear Fitted'
    LEGEND03 = 'Undefined'
    TITLE = f'Last Data of BaseShear-Displacement Analysis - Ductility Ratio: {X[2]/X[1]:.4f} - Over Strength Factor: {Y[2]/Y[1]:.4f}'
    COLOR = 'black'
    BC.PLOT_2D(np.abs(DISP_X), np.abs(FORCE_S), X, Y, X, Y, XLABEL, YLABEL, TITLE, LEGEND01, LEGEND02, LEGEND03, COLOR='black', Z=2) 
    #print(f'\t\t Ductility Ratio: {Y[2]/Y[1]:.4f}')
    
    # Calculate Over Strength Coefficient (Ω0)
    Omega_0 = Y[2] / Y[1]
    # Calculate Displacement Ductility Ratio (μ)
    mu = X[2] / X[1]
    # Calculate Ductility Coefficient (Rμ)
    #R_mu = 1
    #R_mu = (2 * mu - 1) ** 0.5
    R_mu = mu
    # Calculate Structural Behavior Coefficient (R)
    R = Omega_0 * R_mu
    print(f'Over Strength Coefficient (Ω0):      {Omega_0:.4f}')
    print(f'Displacement Ductility Ratio (μ):    {mu:.4f}')
    print(f'Ductility Coefficient (Rμ):          {R_mu:.4f}')
    print(f'Structural Behavior Coefficient (R): {R:.4f}')
    # ---------------------------------------
    #  Plot BaseAxial-Displacement Analysis
    # ---------------------------------------
    XX = np.abs(DISP_Y); YY = np.abs(FORCE_A); # ABSOLUTE VALUE
    SLOPE_NODE = 10
    
    DATA = BC.BILNEAR_CURVE(XX, YY, SLOPE_NODE)
    X, Y, Elastic_ST, Plastic_ST, Tangent_ST, Ductility_Rito, Over_Strength_Factor = DATA
    
    XLABEL = 'Displacement in Y [mm]'
    YLABEL = 'Base-Axial Reaction [N]'
    LEGEND01 = 'Curve'
    LEGEND02 = 'Bilinear Fitted'
    LEGEND03 = 'Undefined'
    TITLE = f'Last Data of BaseAxial-Displacement Analysis - Ductility Ratio: {X[2]/X[1]:.4f} - Over Strength Factor: {Y[2]/Y[1]:.4f}'
    COLOR = 'black'
    BC.PLOT_2D(np.abs(DISP_Y), np.abs(FORCE_A), X, Y, X, Y, XLABEL, YLABEL, TITLE, LEGEND01, LEGEND02, LEGEND03, COLOR='black', Z=2) 
    #print(f'\t\t Ductility Ratio: {YY[2]/YY[1]:.4f}')
    
    # Calculate Over Strength Coefficient (Ω0)
    Omega_0 = Y[2] / Y[1]
    # Calculate Displacement Ductility Ratio (μ)
    mu = X[2] / X[1]
    # Calculate Ductility Coefficient (Rμ)
    #R_mu = 1
    #R_mu = (2 * mu - 1) ** 0.5
    R_mu = mu
    # Calculate Structural Behavior Coefficient (R)
    R = Omega_0 * R_mu
    print(f'Over Strength Coefficient (Ω0):      {Omega_0:.4f}')
    print(f'Displacement Ductility Ratio (μ):    {mu:.4f}')
    print(f'Ductility Coefficient (Rμ):          {R_mu:.4f}')
    print(f'Structural Behavior Coefficient (R): {R:.4f}')
     
#%%------------------------------------------------------------------------------
def PLOT_DYNAMIC(time, DISP_X, DISP_Y, velocity_X, velocity_Y, acceleration_X, acceleration_Y, FORCE_S, FORCE_A, MOMENT, ROT, delta):
    # Compute the Cumulative Maximum Absolute Value of Last Analysis Data
    def MAX_ABS(X):
        import numpy as np
        X = np.asarray(X)  # Convert input to a numpy array for faster operations
        X_MAX = np.zeros_like(X)  # Initialize an array to store cumulative max values
        X_MAX[0] = np.abs(X[0])  # Set the first value

        # Compute cumulative maximum absolute values
        for i in range(1, len(X)):
            X_MAX[i] = max(X_MAX[i-1], np.abs(X[i]))
        
        return X_MAX 
    #------------------------------------------------------------------
    # EXACT SOLUTION:
    from scipy.optimize import fsolve

    # Define the equation for natural logarithm of this ratio, called the logarithmic decrement, we denote by δ
    def EQUATION(x, delta):
        import numpy as np
        if np.any(x == 0):  # Avoid division by zero
            return np.inf  
            
        # Calculate the value of the equation
        A = x**2 - 1 + ((2 * np.pi * x) / np.mean(delta)) ** 2
        #print(f"x: {x}, A: {A}")  # Debugging output
        # Return the difference (for root finding)
        return A
          
    # Initial guess for root(s)
    x0 = 1  # Intial Guess for Damping Ratio
    # Solve for x
    solution = fsolve(EQUATION, x0, args=(delta))
    print(f"Exact Damping Ratio: {solution[0]:.8e}")
    #------------------------------------------------------------------
    DISP_ZX = MAX_ABS(DISP_X)  
    DISP_ZY = MAX_ABS(DISP_Y) 
    VELO_ZX = MAX_ABS(velocity_X) 
    VELO_ZY = MAX_ABS(velocity_Y)
    ACCE_ZX = MAX_ABS(acceleration_X) 
    ACCE_ZY = MAX_ABS(acceleration_Y)
    BASE_ZS = MAX_ABS(FORCE_S)
    BASE_ZA = MAX_ABS(FORCE_A)
    BASE_ZM = MAX_ABS(MOMENT)
    ROT_Z = MAX_ABS(ROT)
    
    plt.figure(1, figsize=(8, 6))
    plt.plot(time, DISP_X, color='blue', linewidth=2)
    plt.plot(time, DISP_ZX, color='red', linewidth=2)
    plt.xlabel('Time [s]')
    plt.ylabel('Displacement in X [mm]')
    plt.title(f'Time vs Displacement - MAX. ABS: {DISP_ZX[-1]} | ξ (Calculated): {100*solution[0]:.5e} %')
    plt.grid()
    plt.savefig(f'DYNAMIC_PLOT_{1}.png')
    plt.show()
    
    plt.figure(2, figsize=(8, 6))
    plt.plot(time, DISP_Y, color='blue', linewidth=2)
    plt.plot(time, DISP_ZY, color='red', linewidth=2)
    plt.xlabel('Time [s]')
    plt.ylabel('Displacement in Y [mm]')
    plt.title(f'Time vs Displacement - MAX. ABS: {DISP_ZY[-1]}')
    plt.grid()
    plt.savefig(f'DYNAMIC_PLOT_{2}.png')
    plt.show()
    
    plt.figure(3, figsize=(8, 6))
    plt.plot(time, velocity_X, color='blue', linewidth=2)
    plt.plot(time, VELO_ZX, color='red', linewidth=2)
    plt.xlabel('Time [s]')
    plt.ylabel('Velocity in X [mm/s]')
    plt.title(f'Time vs Velocity - MAX. ABS: {VELO_ZX[-1]}')
    plt.grid()
    plt.savefig(f'DYNAMIC_PLOT_{3}.png')
    plt.show()
    
    plt.figure(4, figsize=(8, 6))
    plt.plot(time, velocity_Y, color='blue', linewidth=2)
    plt.plot(time, VELO_ZY, color='red', linewidth=2)
    plt.xlabel('Time [s]')
    plt.ylabel('Velocity in Y [mm/s]')
    plt.title(f'Time vs Velocity - MAX. ABS: {VELO_ZY[-1]}')
    plt.grid()
    plt.savefig(f'DYNAMIC_PLOT_{4}.png')
    plt.show()
    
    plt.figure(5, figsize=(8, 6))
    plt.plot(time, acceleration_X, color='blue', linewidth=2)
    plt.plot(time, ACCE_ZX, color='red', linewidth=2)
    plt.xlabel('Time [s]')
    plt.ylabel('Acceleration in X [mm/s^2]')
    plt.title(f'Time vs Acceleration - MAX. ABS: {ACCE_ZX[-1]}')
    plt.grid()
    plt.savefig(f'DYNAMIC_PLOT_{5}.png')
    plt.show()
    
    plt.figure(6, figsize=(8, 6))
    plt.plot(time, acceleration_Y, color='blue', linewidth=2)
    plt.plot(time, ACCE_ZY, color='red', linewidth=2)
    plt.xlabel('Time [s]')
    plt.ylabel('Acceleration in Y [mm/s^2]')
    plt.title(f'Time vs Acceleration - MAX. ABS: {ACCE_ZY[-1]}')
    plt.grid()
    plt.savefig(f'DYNAMIC_PLOT_{6}.png')
    plt.show()
    
    plt.figure(7, figsize=(8, 6))
    plt.plot(time, FORCE_S, color='blue', linewidth=2)
    plt.plot(time, BASE_ZS, color='red', linewidth=2)
    plt.xlabel('Time [s]')
    plt.ylabel('Shear Base-reaction [N]')
    plt.title(f'Time vs Shear Base-reaction - MAX. ABS: {BASE_ZS[-1]}')
    plt.grid()
    plt.savefig(f'DYNAMIC_PLOT_{7}.png')
    plt.show()  
    
    plt.figure(8, figsize=(8, 6))
    plt.plot(time, FORCE_A, color='blue', linewidth=2)
    plt.plot(time, BASE_ZA, color='red', linewidth=2)
    plt.xlabel('Time [s]')
    plt.ylabel('Axial Base-reaction [N]')
    plt.title(f'Time vs Axial Base-reaction - MAX. ABS: {BASE_ZA[-1]}')
    plt.grid()
    plt.savefig(f'DYNAMIC_PLOT_{8}.png')
    plt.show() 
    
    plt.figure(9, figsize=(8, 6))
    plt.plot(time, MOMENT, color='blue', linewidth=2)
    plt.plot(time, BASE_ZM, color='red', linewidth=2)
    plt.title(f'Time vs Moment Base-reaction - MAX. ABS: {BASE_ZM[-1]}')
    plt.xlabel('Time [s]')
    plt.ylabel('Moment Base-reaction [N.mm]')
    plt.grid()
    plt.savefig(f'DYNAMIC_PLOT_{9}.png')
    plt.show() 
    
    plt.figure(10, figsize=(8, 6))
    plt.plot(time, ROT, color='blue', linewidth=2)
    plt.plot(time, ROT_Z, color='red', linewidth=2)
    plt.xlabel('Time [s]')
    plt.ylabel('Rotation in Z [rad]')
    plt.title(f'Time vs Rotation in Z - MAX. ABS: {ROT_Z[-1]}')
    plt.grid()
    plt.savefig(f'DYNAMIC_PLOT_{10}.png')
    plt.show()
    