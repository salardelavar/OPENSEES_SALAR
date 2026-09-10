# DETERMINATION OF OPTIMUM STRUCTURAL INITIAL ELASTIC STIFFNESS FROM STRUCTURAL PERIOD USING FINITE‑DIFFERENCE NEWTON ITERATION AND OPENSEES WITH AND WITHOUT PARALLEL COMPUTING

![alt text](https://github.com/salardelavar/OPENSEES_SALAR/blob/main/EIGHT_ANALYSIS_DUCTILITY_DAMAGE_INDEX_%26_ENERGY_DISSIPATION_CAPACITY_INDEX/11_SDOF_8_ANA/OPTIMIZATION_PERIOD_SDOF/KE_OPTIMIZATION_PERIOD_SDOF/COVER_KE.png) 

 # PERIOD ANALYSIS WITH NEWTON-RAPHSON OPTIMIZATION
--------------------------------------------------------------------------------
1. OBJECTIVE
--------------------------------------------------------------------------------
    Perform a Period analysis of a single-degree-of-freedom
    (SDOF) system and iteratively determine the optimum initial elastic
    stiffness required to achieve a target structural period.

    The material behaviour may be either elastic or inelastic, selected via
    the MAT_TYPE flag ('ELASTIC' or 'INELASTIC').

    In this configuration:
        MAT_TYPE = 'INELASTIC'

    The analysis is driven by a Newton-Raphson root-finding scheme that
    repeatedly calls the nonlinear SDOF solver with ANAL_TYPE = 'PERIOD'
    until the structural period converges to the prescribed demand.

2. OPTIMIZATION STRATEGY
--------------------------------------------------------------------------------
    A Newton-Raphson root-finding scheme is used to solve for the stiffness X
    such that the maximum structural period matches a target demand.

    Target demand:
        DEMAND = 0.10   [in the chosen period units]

    Residual function:
        F(X) = PERIOD_MAX - DEMAND

        where PERIOD_MAX is returned by the SDOF solver when called with
        ANAL_TYPE = 'PERIOD'. The residual measures how far the maximum
        structural period is from the target. A value of zero indicates
        exact convergence.

    Derivative approximation:
        Central finite differences with step ESP = 1e-3

            dF/dX ~ [F(X + ESP) - F(X - ESP)] / (2 * ESP)

        Central differences are preferred over forward or backward
        differences because they cancel the leading-order truncation error
        and yield second-order accuracy.

    Update rule:
        X <- X - F / (dF/dX)

        The stiffness is corrected along the tangent direction of the
        residual function until the root is reached.

    Convergence criteria:
        - The iteration stops when |DX| < TOLERANCE = 1e-6
        - Or when ITMAX = 100000 iterations are reached (safety limit)

    Each evaluation of F requires a full nonlinear analysis via the SDOF
    routine, which returns:
        (PERIOD_MIN, PERIOD_MAX)

    The solver evaluates F at three points per iteration:
        - X          (current estimate)
        - X - ESP    (lower finite-difference point)
        - X + ESP    (upper finite-difference point)

    This triple evaluation is what makes each Newton step computationally
    equivalent to three full SDOF period analyses.

3. WORKFLOW SUMMARY
--------------------------------------------------------------------------------
    Step 1 : Initialize X (initial stiffness guess) and solver parameters.
    Step 2 : Evaluate F(X)         via SDOF(ANAL_TYPE='PERIOD').
    Step 3 : Evaluate F(X - ESP)   via SDOF(ANAL_TYPE='PERIOD').
    Step 4 : Evaluate F(X + ESP)   via SDOF(ANAL_TYPE='PERIOD').
    Step 5 : Compute the central finite-difference derivative dF/dX.
    Step 6 : Compute the Newton update DX = F / (dF/dX).
    Step 7 : Compute the residual RESIDUAL = |DX|.
    Step 8 : Update X <- X - DX and increment the iteration counter.
    Step 9 : Check termination (TOLERANCE or ITMAX).
    Step 10: Report the optimum stiffness, iteration count, residual, and
             total elapsed CPU time.

4. INPUT PARAMETERS
--------------------------------------------------------------------------------
    MAT_TYPE  : Material behaviour                         ('INELASTIC')
    X         : Initial guess for elastic stiffness        (4500000.0 N/m)
    ESP       : Finite-difference step                     (1e-3)
    TOLERANCE : Convergence tolerance on DX                (1e-6)
    RESIDUAL  : Initial residual (any value > TOLERANCE)   (100)
    IT        : Iteration counter                          (0)
    ITMAX     : Maximum allowed iterations                 (100000)
    DEMAND    : Target structural period                   (0.10)
    ANAL_TYPE : Analysis type passed to SDOF               ('PERIOD')

    Note: TOTAL_MASS must be defined in the calling scope before the loop.

5. OUTPUT
--------------------------------------------------------------------------------
    - Optimum initial elastic stiffness X                [N/m]
    - Iteration count at convergence
    - Final convergence residual
    - Total elapsed CPU time in seconds

6. NOTES ON NUMERICAL ROBUSTNESS
--------------------------------------------------------------------------------
    - The choice of ESP balances truncation error (too large) against
      subtractive cancellation and round-off (too small).
    - The central-difference scheme is second-order accurate in ESP.
    - If ITMAX is reached without convergence, the user is advised to
      adjust ESP and TOLERANCE to improve conditioning.
    - Because the SDOF solver may involve non-smooth hysteretic behaviour,
      convergence is not guaranteed for every configuration; the
      tolerance and finite-difference step may require tuning.
    - The initial stiffness guess X = 4.5e6 N/m should be chosen close
      enough to the expected solution to keep the Newton iteration stable.

7. VARIABLE DICTIONARY
--------------------------------------------------------------------------------
    X          : Current elastic stiffness estimate              [N/m]
    Xmin       : Lower finite-difference point   (X - ESP)       [N/m]
    Xmax       : Upper finite-difference point   (X + ESP)       [N/m]
    F          : Residual at X                   (PERIOD_MAX - DEMAND)
    Fmin       : Residual at Xmin
    Fmax       : Residual at Xmax
    DF         : Central finite-difference derivative dF/dX
    DX         : Newton update step               (F / DF)
    RESIDUAL   : Absolute value of DX             (|DX|)
    IT         : Iteration counter
    ITMAX      : Maximum allowed iterations
    SUPPLY     : PERIOD_MAX at X
    SUPPLYmin  : PERIOD_MAX at Xmin
    SUPPLYmax  : PERIOD_MAX at Xmax
    DEMAND     : Target structural period
    TOLERANCE  : Convergence tolerance on DX
    ESP        : Finite-difference step
THIS PYTHON SCRIPT IS WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)


