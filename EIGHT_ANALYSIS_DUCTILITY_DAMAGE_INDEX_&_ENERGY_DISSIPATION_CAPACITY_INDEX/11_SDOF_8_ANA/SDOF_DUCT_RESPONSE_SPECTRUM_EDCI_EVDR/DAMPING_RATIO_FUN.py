def DAMPING_RATIO(disp):
    from scipy.optimize import fsolve
    import numpy as np
    
    # Convert to numpy array
    displacement = np.array(disp, dtype=float)
    
    # Find peaks
    peaks = np.array([
        displacement[i] for i in range(1, len(displacement) - 1)
        if displacement[i] > displacement[i - 1] and displacement[i] > displacement[i + 1]
    ])
    
    # Avoid division by zero or negative peaks
    if len(peaks) < 2 or np.any(peaks <= 0):
        return 0.0
    
    # Calculate logarithmic decrement
    delta = np.log(peaks[:-1] / peaks[1:])
    
    # Remove NaNs or infinite values
    delta = delta[np.isfinite(delta)]
    if len(delta) == 0:
        return 0.0
    
    # Define equation for root finding
    def EQUATION(x, delta):
        if np.any(x == 0):  # avoid divide by zero
            return np.inf
        return x**2 - 1 + ((2 * np.pi * x) / np.mean(delta)) ** 2
          
    # Solve for damping ratio
    x0 = 1.0
    solution = fsolve(EQUATION, x0, args=(delta))
    
    if np.isfinite(solution[0]) and solution[0] != 1:
        print(f"Exact Damping Ratio: {solution[0]:.8e}")
        return 100 * solution[0]
    else:
        return 0.0