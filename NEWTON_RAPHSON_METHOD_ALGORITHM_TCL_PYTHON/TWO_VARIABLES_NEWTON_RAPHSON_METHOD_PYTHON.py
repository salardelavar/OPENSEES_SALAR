# ############################################################################
#     >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<       #
#                       NEWTON-RAPHSON ALGORITHM METHOD                      #
#     THIS PYTHON SCRIPT WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)       #
# ############################################################################
import time as TI
import numpy as np
 
X01 = 100         # Intial Guess for X1
X02 = 100         # Intial Guess for X2
ESP = 1e-3        # Finite difference derivative Convergence Tolerance
TOLERANCE = 1e-6  # Convergence Tolerance
RESIDUAL = 100    # Convergence Residual 
IT = 0            # Intial Iteration
ITMAX = 100000    # Max. Iteration
DEMAND = 0.0      # Target Value 

# Analysis Durations:
starttime = TI.process_time()

# ---------------------------------------------------------------------------
# FIND THE OPTIMUM VALUE (NEWTON-RAPHSON SOLVER FOR OPTIMAL X1 and X2)
# ---------------------------------------------------------------------------
def FUNCTION(X1, X2):
    A01 = (X1**2 + X2**2 - 25)           # FUNCTION 01
    A02 = ( X1*X2 - 12)                  # FUNCTION 02
    return A01, A02
while (RESIDUAL > TOLERANCE):
    # X -------------------
    SUPPLY01, SUPPLY02 = FUNCTION(X01, X02)
    F01 = SUPPLY01 - DEMAND
    F02 = SUPPLY02 - DEMAND
    print('F:    ', F01, F02)
    # XMIN -------------------
    # Evaluate at Xmin and Fmin
    Xmin01 = X01 - ESP
    Xmin02 = X02 - ESP
    print('Xmin:    ', Xmin01, Xmin02)
    SUPPLYmin01, SUPPLYmin02 = FUNCTION(Xmin01, Xmin02)
    Fmin01 = SUPPLYmin01 - DEMAND
    Fmin02 = SUPPLYmin02 - DEMAND
    print('Fmin: ', Fmin01, Fmin02)
    # XMAX -------------------
    # Evaluate at Xmax and Fmax
    Xmax01 = X01 + ESP
    Xmax02 = X02 + ESP
    print('Xmax:    ', Xmax01, Xmax02)
    SUPPLYmax01, SUPPLYmax02 = FUNCTION(Xmax01, Xmax02)
    Fmax01 = SUPPLYmax01 - DEMAND
    Fmax02 = SUPPLYmax02 - DEMAND
    print('Fmax: ', Fmax01, Fmax02)
    # DF -------------------
    DF01 = (Fmax01 - Fmin01)/(2 * ESP);# Calculate the Finite difference derivative of F1
    DF02 = (Fmax02 - Fmin02)/(2 * ESP);# Calculate the Finite difference derivative of F2
    print('DF:   ', DF01, DF02)
    # DX -------------------
    DX01 = F01 / DF01; # Calculate dx
    DX02 = F02 / DF02; # Calculate dx
    print('DX:   ', DX01, DX02)
    # RESIDUAL -------------------
    RESIDUAL = np.abs(max(DX01, DX02)); # Calculate residual
    print('IT: ', IT + 1, ' - RESIDUAL: ', RESIDUAL,'\n\nX01: ', X01,' - X02: ', X02,'\n')
    X01 -= DX01; # update X1
    X02 -= DX02; # update X2
    IT += 1; # update iteration
    # CONTROLLING -------------------
    if IT == ITMAX:
        print("\t\t Iteration reached to Max. Iteration")
        print("\t\t Change ESP and TOLERANCE for better Convergence")
        X01 =- DX01 # update X1
        X02 =- DX02 # update X2
        break;
    if RESIDUAL < TOLERANCE:
        print(f'\t\t Optimum X1:                      {X01:.4f}')
        print(f'\t\t Optimum X2:                      {X02:.4f}')
        print(f'\t\t Iteration Counts:                {IT}')
        print(f'\t\t Convergence Residual:            {RESIDUAL:.10e}')
    #print(X)

totaltime = TI.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n')