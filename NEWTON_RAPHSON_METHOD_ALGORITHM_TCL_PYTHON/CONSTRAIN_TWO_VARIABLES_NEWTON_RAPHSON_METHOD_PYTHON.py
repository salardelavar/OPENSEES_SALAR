# ############################################################################
#     >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<       #
#                       NEWTON-RAPHSON ALGORITHM METHOD                      #
#                 WITH CONSTRAINED OPTIMIZATION (UB, LB)                     #
#     THIS PYTHON SCRIPT WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)       #
# ############################################################################
import time as TI
import numpy as np

# Initial Parameters
X01 = 100         # Initial Guess for X1
X02 = 100         # Initial Guess for X2
ESP = 1e-3        # Finite difference derivative Convergence Tolerance
TOLERANCE = 1e-6  # Convergence Tolerance
RESIDUAL = 100    # Convergence Residual 
IT = 0            # Initial Iteration
ITMAX = 100000    # Max. Iteration
DEMAND = 0.0      # Target Value 

# BOUNDS FOR CONSTRAINED OPTIMIZATION
# Upper bounds (Ub) and Lower bounds (Lb) for [X1, X2]
LB = np.array([-10.0, -10.0])   # Lower bounds for X1 and X2
UB = np.array([10.0, 10.0])     # Upper bounds for X1 and X2

# Analysis Durations:
starttime = TI.process_time()

# ---------------------------------------------------------------------------
# FIND THE OPTIMUM VALUE (NEWTON-RAPHSON SOLVER FOR OPTIMAL X1 and X2)
# WITH CONSTRAINT PROJECTION
# ---------------------------------------------------------------------------
def FUNCTION(X1, X2):
    A01 = (X1**2 + X2**2 - 25)           # FUNCTION 01
    A02 = (X1*X2 - 12)                   # FUNCTION 02
    return A01, A02

def project_to_bounds(x, lb, ub):
    """Project the solution vector onto the feasible bounds"""
    return np.clip(x, lb, ub)

def check_bounds_violation(x, lb, ub, tolerance=1e-10):
    """Check if the current solution violates bounds"""
    violation_lower = np.any(x < lb - tolerance)
    violation_upper = np.any(x > ub + tolerance)
    return violation_lower or violation_upper

while (RESIDUAL > TOLERANCE) and (IT < ITMAX):
    # Check and project initial bounds
    X01 = np.clip(X01, LB[0], UB[0])
    X02 = np.clip(X02, LB[1], UB[1])
    
    # X -------------------
    SUPPLY01, SUPPLY02 = FUNCTION(X01, X02)
    F01 = SUPPLY01 - DEMAND
    F02 = SUPPLY02 - DEMAND
    print(f'F:    {F01:.6f}, {F02:.6f}')
    
    # XMIN -------------------
    # Evaluate at Xmin and Fmin
    Xmin01 = np.clip(X01 - ESP, LB[0], UB[0])
    Xmin02 = np.clip(X02 - ESP, LB[1], UB[1])
    print(f'Xmin: {Xmin01:.6f}, {Xmin02:.6f}')
    
    SUPPLYmin01, SUPPLYmin02 = FUNCTION(Xmin01, Xmin02)
    Fmin01 = SUPPLYmin01 - DEMAND
    Fmin02 = SUPPLYmin02 - DEMAND
    print(f'Fmin: {Fmin01:.6f}, {Fmin02:.6f}')
    
    # XMAX -------------------
    # Evaluate at Xmax and Fmax
    Xmax01 = np.clip(X01 + ESP, LB[0], UB[0])
    Xmax02 = np.clip(X02 + ESP, LB[1], UB[1])
    print(f'Xmax: {Xmax01:.6f}, {Xmax02:.6f}')
    
    SUPPLYmax01, SUPPLYmax02 = FUNCTION(Xmax01, Xmax02)
    Fmax01 = SUPPLYmax01 - DEMAND
    Fmax02 = SUPPLYmax02 - DEMAND
    print(f'Fmax: {Fmax01:.6f}, {Fmax02:.6f}')
    
    # DF -------------------
    # Calculate the Finite difference derivative with safeguards
    denom1 = 2 * ESP if (Xmax01 - Xmin01) > 1e-15 else 1.0
    denom2 = 2 * ESP if (Xmax02 - Xmin02) > 1e-15 else 1.0
    
    DF01 = (Fmax01 - Fmin01) / denom1  # Calculate the Finite difference derivative of F1
    DF02 = (Fmax02 - Fmin02) / denom2  # Calculate the Finite difference derivative of F2
    
    # Avoid division by zero or very small derivatives
    if abs(DF01) < 1e-15:
        DF01 = 1e-15 * np.sign(DF01) if DF01 != 0 else 1e-15
    if abs(DF02) < 1e-15:
        DF02 = 1e-15 * np.sign(DF02) if DF02 != 0 else 1e-15
    
    print(f'DF:   {DF01:.6f}, {DF02:.6f}')
    
    # DX -------------------
    DX01 = F01 / DF01  # Calculate dx
    DX02 = F02 / DF02  # Calculate dx
    
    # Apply damping factor for better convergence near bounds
    damping_factor = 0.7  # Helps prevent oscillation near bounds
    DX01 *= damping_factor
    DX02 *= damping_factor
    
    print(f'DX:   {DX01:.6f}, {DX02:.6f}')
    
    # RESIDUAL -------------------
    RESIDUAL = np.abs(max(DX01, DX02))  # Calculate residual
    print(f'IT: {IT + 1} - RESIDUAL: {RESIDUAL:.10e}\n')
    print(f'X01: {X01:.6f} - X02: {X02:.6f}\n')
    
    # Update with projection to bounds
    X01_new = X01 - DX01
    X02_new = X02 - DX02
    
    # Project to bounds
    X01 = np.clip(X01_new, LB[0], UB[0])
    X02 = np.clip(X02_new, LB[1], UB[1])
    
    # Check if we hit bounds and adjust
    if X01 == LB[0] or X01 == UB[0]:
        print(f"\t\t X1 hit bound: {X01:.6f}")
    if X02 == LB[1] or X02 == UB[1]:
        print(f"\t\t X2 hit bound: {X02:.6f}")
    
    IT += 1  # update iteration
    
    # CONTROLLING -------------------
    if IT == ITMAX:
        print("\t\t Iteration reached to Max. Iteration")
        print("\t\t Change ESP and TOLERANCE for better Convergence")
        break
    
    if RESIDUAL < TOLERANCE:
        print(f'\t\t Optimum X1:                      {X01:.6f}')
        print(f'\t\t Optimum X2:                      {X02:.6f}')
        print(f'\t\t Iteration Counts:                {IT}')
        print(f'\t\t Convergence Residual:            {RESIDUAL:.10e}')
        
        # Check final constraint satisfaction
        F1_final, F2_final = FUNCTION(X01, X02)
        constraint_violation = np.abs(F1_final) + np.abs(F2_final)
        print(f'\t\t Constraint Violation:            {constraint_violation:.10e}')
        
        # Check if solution is at bounds
        at_lower_bound = [X01 <= LB[0] + TOLERANCE, X02 <= LB[1] + TOLERANCE]
        at_upper_bound = [X01 >= UB[0] - TOLERANCE, X02 >= UB[1] - TOLERANCE]
        
        for i, (at_lb, at_ub) in enumerate(zip(at_lower_bound, at_upper_bound)):
            if at_lb:
                print(f'\t\t X{i+1} at LOWER bound: {LB[i]:.6f}')
            if at_ub:
                print(f'\t\t X{i+1} at UPPER bound: {UB[i]:.6f}')

# Final summary
totaltime = TI.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f}')

# Print final bounds information
print("\n" + "="*60)
print("CONSTRAINT BOUNDS SUMMARY:")
print(f"Lower Bounds (LB): X1 ≥ {LB[0]:.4f}, X2 ≥ {LB[1]:.4f}")
print(f"Upper Bounds (UB): X1 ≤ {UB[0]:.4f}, X2 ≤ {UB[1]:.4f}")
print(f"Final Solution:    X1 = {X01:.6f}, X2 = {X02:.6f}")
print("="*60 + "\n")