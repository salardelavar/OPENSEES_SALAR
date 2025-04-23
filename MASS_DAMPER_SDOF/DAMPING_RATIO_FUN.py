# %% DAMPING RATIO EXACT SOLUTION:
def DAMPING_RATIO(DISP):
    import numpy as np
    from scipy.optimize import fsolve
    # Calculating Damping Ratio Using Logarithmic Decrement Analysis    
    peaks = np.array([DISP[i] for i in range(1, len(DISP)-1) if DISP[i] > DISP[i-1] and DISP[i] > DISP[i+1]])
    #print(peaks)
    # Natural logarithm
    delta = np.log(peaks[:-1] / peaks[1:])
    #print(delta)
    # Define the equation for natural logarithm of this ratio, called the logarithmic decrement, we denote by δ
    def EQUATION(x, delta):
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
    return solution