# ############################################################################
#     >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<       #
#                       NEWTON-RAPHSON ALGORITHM METHOD                      #
#     THIS PYTHON SCRIPT WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)       #
# ############################################################################
import time as TI
import numpy as np

X = 100           # Intial Guess
ESP = 1e-3        # Finite difference derivative Convergence Tolerance
TOLERANCE = 1e-6  # Convergence Tolerance
RESIDUAL = 100    # Convergence Residual 
IT = 0            # Intial Iteration
ITMAX = 100000    # Max. Iteration
DEMAND = 60.0       # Target Value 

# Analysis Durations:
starttime = TI.process_time()

# ---------------------------------------------------------------------------
# FIND THE OPTIMUM VALUE (NEWTON-RAPHSON SOLVER FOR OPTIMAL X)
# ---------------------------------------------------------------------------
def FUNCTION(X):
    return (X**3 - 4 *X**2 + X**1 -250.0)
while (RESIDUAL > TOLERANCE):
    # X -------------------
    answer = FUNCTION(X)
    SUPPLY = answer
    F = SUPPLY - DEMAND
    print('F:    ', F)
    # XMIN -------------------
    # Evaluate at Xmin and Fmin
    Xmin = X - ESP
    print('Xmin:    ', Xmin)
    answer = FUNCTION(Xmin)
    SUPPLYmin = answer
    Fmin = SUPPLYmin - DEMAND
    print('Fmin: ', Fmin)
    # XMAX -------------------
    # Evaluate at Xmax and Fmax
    Xmax = X + ESP
    print('Xmmax:    ', Xmax)
    answer = FUNCTION(Xmax)
    SUPPLYmax = answer
    Fmax = SUPPLYmax - DEMAND
    print('Fmax: ', Fmax)
    # DF -------------------
    DF = (Fmax - Fmin)/(2 * ESP);# Calculate the Finite difference derivative of F
    print('DF:   ', DF)
    # DX -------------------
    DX = F / DF; # Calculate dx
    print('DX:   ', DX)
    # RESIDUAL -------------------
    RESIDUAL = np.abs(DX); # Calculate residual
    print('IT: ', IT + 1, ' - RESIDUAL: ', RESIDUAL,' - X: ', X,'\n')
    X -= DX; # update X
    IT += 1; # update iteration
    # CONTROLLING -------------------
    if IT == ITMAX:
        print("\t\t Iteration reached to Max. Iteration")
        print("\t\t Change ESP and TOLERANCE for better Convergence")
        X =- DX # update X
        break;
    if RESIDUAL < TOLERANCE:
        print(f'\t\t Optimum X:  {X:.4f}')
        print(f'\t\t Iteration Counts:                {IT}')
        print(f'\t\t Convergence Residual:            {RESIDUAL:.10e}')
    #print(X)

totaltime = TI.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n')