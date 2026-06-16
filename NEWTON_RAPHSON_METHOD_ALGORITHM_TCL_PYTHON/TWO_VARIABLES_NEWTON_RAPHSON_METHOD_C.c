/******************************************************************************
 *     >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<       *
 *                       NEWTON-RAPHSON ALGORITHM METHOD                      *
 *     THIS C SCRIPT TRANSLATED FROM PYTHON BY SALAR DELAVAR GHASHGHAEI       *
 ******************************************************************************/

#include <stdio.h>
#include <math.h>
#include <time.h>

/* Function to evaluate the system of equations */
void FUNCTION(double X1, double X2, double *A01, double *A02) {
    *A01 = (X1 * X1 + X2 * X2 - 25.0);    /* FUNCTION 01 */
    *A02 = (X1 * X2 - 12.0);              /* FUNCTION 02 */
}

int main(void) {
    double X01 = 100.0;        /* Initial Guess for X1 */
    double X02 = 100.0;        /* Initial Guess for X2 */
    const double ESP = 1e-3;   /* Finite difference derivative Convergence Tolerance */
    const double TOLERANCE = 1e-6;  /* Convergence Tolerance */
    double RESIDUAL = 100.0;   /* Convergence Residual */
    int IT = 0;                /* Initial Iteration */
    const int ITMAX = 100000;  /* Max. Iteration */
    const double DEMAND = 0.0; /* Target Value */

    double SUPPLY01, SUPPLY02;
    double F01, F02;
    double Xmin01, Xmin02, SUPPLYmin01, SUPPLYmin02, Fmin01, Fmin02;
    double Xmax01, Xmax02, SUPPLYmax01, SUPPLYmax02, Fmax01, Fmax02;
    double DF01, DF02;
    double DX01, DX02;

    /* Analysis Duration start */
    clock_t starttime = clock();

    /* --------------------------------------------------------------------- */
    /* FIND THE OPTIMUM VALUE (NEWTON-RAPHSON SOLVER FOR OPTIMAL X1 and X2)  */
    /* --------------------------------------------------------------------- */
    while (RESIDUAL > TOLERANCE) {
        /* X ------------------- */
        FUNCTION(X01, X02, &SUPPLY01, &SUPPLY02);
        F01 = SUPPLY01 - DEMAND;
        F02 = SUPPLY02 - DEMAND;
        printf("F:     %g %g\n", F01, F02);

        /* XMIN ------------------- */
        /* Evaluate at Xmin and Fmin */
        Xmin01 = X01 - ESP;
        Xmin02 = X02 - ESP;
        printf("Xmin:  %g %g\n", Xmin01, Xmin02);
        FUNCTION(Xmin01, Xmin02, &SUPPLYmin01, &SUPPLYmin02);
        Fmin01 = SUPPLYmin01 - DEMAND;
        Fmin02 = SUPPLYmin02 - DEMAND;
        printf("Fmin:  %g %g\n", Fmin01, Fmin02);

        /* XMAX ------------------- */
        /* Evaluate at Xmax and Fmax */
        Xmax01 = X01 + ESP;
        Xmax02 = X02 + ESP;
        printf("Xmax:  %g %g\n", Xmax01, Xmax02);
        FUNCTION(Xmax01, Xmax02, &SUPPLYmax01, &SUPPLYmax02);
        Fmax01 = SUPPLYmax01 - DEMAND;
        Fmax02 = SUPPLYmax02 - DEMAND;
        printf("Fmax:  %g %g\n", Fmax01, Fmax02);

        /* DF ------------------- */
        DF01 = (Fmax01 - Fmin01) / (2.0 * ESP);  /* Calculate the Finite difference derivative of F1 */
        DF02 = (Fmax02 - Fmin02) / (2.0 * ESP);  /* Calculate the Finite difference derivative of F2 */
        printf("DF:    %g %g\n", DF01, DF02);

        /* DX ------------------- */
        DX01 = F01 / DF01;  /* Calculate dx */
        DX02 = F02 / DF02;  /* Calculate dx */
        printf("DX:    %g %g\n", DX01, DX02);

        /* RESIDUAL ------------------- */
        RESIDUAL = fmax(fabs(DX01), fabs(DX02));  /* Calculate residual */
        printf("IT: %d - RESIDUAL: %g - X01: %g - X02: %g\n\n", IT + 1, RESIDUAL, X01, X02);

        X01 -= DX01;  /* update X1 */
        X02 -= DX02;  /* update X2 */
        IT++;          /* update iteration */

        /* CONTROLLING ------------------- */
        if (IT == ITMAX) {
            printf("\t\t Iteration reached to Max. Iteration\n");
            printf("\t\t Change ESP and TOLERANCE for better Convergence\n");
            X01 -= DX01;  /* update X1 */
            X02 -= DX02;  /* update X2 */
            break;
        }
        if (RESIDUAL < TOLERANCE) {
            printf("\t\t Optimum X1:                      %.4f\n", X01);
            printf("\t\t Optimum X2:                      %.4f\n", X02);
            printf("\t\t Iteration Counts:                %d\n", IT);
            printf("\t\t Convergence Residual:            %.10e\n", RESIDUAL);
        }
    }

    clock_t endtime = clock();
    double totaltime = (double)(endtime - starttime) / CLOCKS_PER_SEC;
    printf("\nTotal time (s): %.4f \n\n", totaltime);

    return 0;
}