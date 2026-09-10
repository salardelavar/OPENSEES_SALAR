def DAMPING_RATIO(disp):
    from scipy.optimize import fsolve
    import numpy as np
    
    # Calculating Damping Ratio Using Logarithmic Decrement Analysis 
    displacement = np.array(disp)
    peaks = np.array([displacement[i] for i in range(1, len(displacement)-1) if displacement[i] > displacement[i-1] and displacement[i] > displacement[i+1]])
    
    # Natural logarithm
    delta = np.log(peaks[:-1] / peaks[1:]) 
    
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
    if solution[0] != 1: 
        print(f"Exact Damping Ratio: {solution[0]:.8e}")
        return 100*solution[0]
    else:
        return 0.0