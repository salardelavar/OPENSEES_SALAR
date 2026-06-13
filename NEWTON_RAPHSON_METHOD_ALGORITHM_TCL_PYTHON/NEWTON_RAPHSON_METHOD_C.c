/*
############################################################################
    >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<
                      NEWTON-RAPHSON ALGORITHM METHOD
    THIS C PROGRAM CONVERTED FROM PYTHON SCRIPT BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
############################################################################
*/
#include <stdio.h>
#include <math.h>
#include <time.h>

#define DEMAND      60.0      // Target Value
#define ESP         1e-3      // Finite difference derivative Convergence Tolerance
#define TOLERANCE   1e-6      // Convergence Tolerance
#define ITMAX       100000    // Max. Iteration

/* The function to be evaluated: SUPPLY = f(X) */
double FUNCTION(double X) {
    return (X*X*X - 4.0*X*X + X - 250.0);
}

int main() {
    double X = 100.0;          // Initial Guess
    double RESIDUAL = 100.0;   // Convergence Residual
    int IT = 0;                // Initial Iteration

    clock_t start_time, end_time;
    start_time = clock();

    while (RESIDUAL > TOLERANCE) {
        double SUPPLY = FUNCTION(X);
        double F = SUPPLY - DEMAND;
        printf("F:    %f\n", F);

        // Evaluate at Xmin
        double Xmin = X - ESP;
        printf("Xmin:    %f\n", Xmin);
        double SUPPLYmin = FUNCTION(Xmin);
        double Fmin = SUPPLYmin - DEMAND;
        printf("Fmin: %f\n", Fmin);

        // Evaluate at Xmax
        double Xmax = X + ESP;
        printf("Xmmax:    %f\n", Xmax);
        double SUPPLYmax = FUNCTION(Xmax);
        double Fmax = SUPPLYmax - DEMAND;
        printf("Fmax: %f\n", Fmax);

        // Central finite difference derivative
        double DF = (Fmax - Fmin) / (2.0 * ESP);
        printf("DF:   %f\n", DF);

        // Newton step
        double DX = F / DF;
        printf("DX:   %f\n", DX);

        RESIDUAL = fabs(DX);
        printf("IT: %d - RESIDUAL: %.10e - X: %f\n\n", IT + 1, RESIDUAL, X);

        X -= DX;   // Update X
        IT++;      // Update iteration

        if (IT == ITMAX) {
            printf("\t\t Iteration reached to Max. Iteration\n");
            printf("\t\t Change ESP and TOLERANCE for better Convergence\n");
            X = -DX;   // As in original code
            break;
        }
    }

    if (RESIDUAL < TOLERANCE) {
        printf("\t\t Optimum X:                     %.4f\n", X);
        printf("\t\t Iteration Counts:                %d\n", IT);
        printf("\t\t Convergence Residual:            %.10e\n", RESIDUAL);
    }

    end_time = clock();
    double total_time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
    printf("\nTotal time (s): %.4f \n\n", total_time);

    return 0;
}