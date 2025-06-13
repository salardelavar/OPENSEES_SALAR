# EXACT SOLUTION:
from scipy.optimize import fsolve
import numpy as np

# Define the equation for natural logarithm of this ratio, called the logarithmic decrement, we denote by δ
def EQUATION(x, delta):
    if np.any(x == 0):  # Avoid division by zero
        return np.inf  
        
    # Calculate the value of the equation
    A = x**2 - 1 + ((2 * np.pi * x) / np.mean(delta)) ** 2
    #print(f"x: {x}, A: {A}")  # Debugging output
    # Return the difference (for root finding)
    return A

def DAMPING_RATIO(delta):
    # Initial guess for root(s)
    x0 = 1  # Intial Guess for Damping Ratio
    # Solve for x
    solution = fsolve(EQUATION, x0, args=(delta))
    print(f"Exact Damping Ratio: {solution[0]:.8e}")
    return solution[0]