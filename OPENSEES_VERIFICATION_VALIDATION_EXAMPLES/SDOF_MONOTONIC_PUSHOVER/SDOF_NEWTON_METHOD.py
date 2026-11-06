###################################################################################################
#                >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                 #
#      MONOTONIC PUSHOVER ANALYSIS OF A SINGLE-DEGREE-OF-FREEDOM (SDOF) NONLINEAR SPRING WITH     #
#                                         DISPLACEMENT CONTROL                                    #
#-------------------------------------------------------------------------------------------------#
#                  FINITE ELEMENT METHOD IN PYTHON AND OPENSEESPY VERIFICATIONS                   #
#-------------------------------------------------------------------------------------------------#
#                   THIS PROGRAM IS WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                 # 
#                               E-mail:salar.d.ghashghaei@gmail.com                               #
###################################################################################################
"""
DISPLACEMENT-CONTROLLED NEWTON-RAPHSON MONOTONIC PUSHOVER OF PARALLEL NONLINEAR SPRINGS:
    
 This Python script implements and verifies a displacement-controlled Newton-Raphson monotonic
 pushover analysis of a single-degree-of-freedom system composed of parallel nonlinear springs.
 It defines a custom numerical solver (`SDOF_NUMERICAL_PUSHOVER_NEWTON`) that applies incremental
 displacements, predicts the load factor from the current tangent stiffness, updates internal forces
 from a multilinear backbone curve scaled by coefficients `Coff`, and applies a single-step Newton
 corrector (which simplifies to adjusting the load factor by the residual) until equilibrium is
 satisfied within a tight tolerance. The same problem is then solved with OpenSeesPy using zeroLength
 elements and HystereticSM (or Elastic) materials under DisplacementControl, after which the resulting
 force-displacement curves from both approaches are compared for accuracy, residual statistics are
 reported, and the curves are plotted.

----------------------------------------------------------------------
Unknowns          : displacement u of the control DOF and load factor lambda
Equilibrium       : R(u, lambda) = lambda * P_ref - F_int(u)          (unbalanced force)
Internal force    : F_int(u) = sum_j  Coff_j * sign(u) * Fb(|u|)       (Fb = backbone)
Tangent stiffness : K_T(u)   = sum_j  Coff_j * dFb/d|u|                (segment slope)
Constraint        : u_(n+1) = u_n + dU                                 (imposed increment)

Per increment
  1. PREDICTOR : d_lambda = K_T(u_n) * dU / P_ref   -> u = u_n + dU,  lambda = lambda_n + d_lambda
  2. Update spring states at the trial u, compute  R = lambda*P_ref - F_int(u)
  3. CORRECTOR (Newton) : solve  K_T*du_R = R  and  K_T*du_P = P_ref, then choose
        d_lambda = -du_R / du_P   so that the control DOF does not move  (du_R + d_lambda*du_P = 0)
     K_T cancels in the ratio, so d_lambda = -R / P_ref (valid even when K_T <= 0, i.e. softening)
  4. Check  |R| / F_ref <= TOL ; repeat 3-4 otherwise

Because there is only ONE controlled DOF, the corrector removes the unbalance in a single
iteration. A residual appears only where the predictor tangent differs from the true
secant path, i.e. when a step crosses a backbone kink (change of segment).
"""
#%% ---------------------------------------
import numpy as np

def SDOF_NUMERICAL_PUSHOVER_NEWTON(Displacement, Force, Coff, ANA_STEPS, MAX_ITERATIONS, TOL, Dmax=None):
    """
    Returns
    -------
    U, F        : displacement and base reaction histories (including the origin)
    R_init      : normalised unbalance right after the predictor, per increment
    R_final     : normalised unbalance after the corrector, per increment
    ITER        : number of corrector iterations per increment
    CONVERGED   : bool array, one flag per increment
    """
    D_bb = np.concatenate(([0.0], np.asarray(Displacement, dtype=float)))
    F_bb = np.concatenate(([0.0], np.asarray(Force, dtype=float)))
    Coff = np.asarray(Coff, dtype=float)

    slopes = np.diff(F_bb) / np.diff(D_bb)          # TANGENT of each segment (not secant)
    Dmax = D_bb[-1] if Dmax is None else float(Dmax)
    u_inc = Dmax / ANA_STEPS
    P_ref = 1.0
    F_ref = Coff.sum() * F_bb.max()                  # force scale for the normalised residual

    def spring_state(u):
        """Force vector and tangent vector of all springs at displacement u."""
        a = abs(u)
        f = np.sign(u) * Coff * np.interp(a, D_bb, F_bb)       # flat beyond the last node
        seg = np.searchsorted(D_bb, a, side='right') - 1       # segment ahead of |u|
        k = Coff * (slopes[seg] if seg < len(slopes) else 0.0)
        return f, k

    # committed state
    u, lam = 0.0, 0.0
    f, k = spring_state(u)
    Kt = k.sum()

    U, F = [0.0], [0.0]
    R_init, R_final, ITER, CONVERGED = [], [], [], []

    for i in range(1, ANA_STEPS + 1):
        u_target = Dmax if i == ANA_STEPS else u_inc * i     # no floating-point overshoot
        dU = u_target - u

        # ---- 1. predictor with the committed tangent ----
        lam_t = lam + Kt * dU / P_ref
        u_t = u + dU

        # ---- 2. state update and unbalanced force ----
        f, k = spring_state(u_t)
        R = lam_t * P_ref - f.sum()
        R_init.append(abs(R) / F_ref)

        # ---- 3. Newton corrector at fixed control displacement ----
        it = 0
        while abs(R) / F_ref > TOL and it < MAX_ITERATIONS:
            d_lam = -R / P_ref                    # K_T cancels (see header)
            lam_t += d_lam
            f, k = spring_state(u_t)              # u_t is unchanged by the constraint
            R = lam_t * P_ref - f.sum()
            it += 1

        ok = abs(R) / F_ref <= TOL
        R_final.append(abs(R) / F_ref)
        ITER.append(it)
        CONVERGED.append(ok)

        if not ok:
            print(f'(-) Increment {i}: NOT converged, |R|/Fref = {abs(R)/F_ref:.3e}. '
                  f'Analysis stopped, state not committed.')
            break

        # ---- commit ----
        u, lam, Kt = u_t, lam_t, k.sum()
        U.append(u)
        F.append(lam * P_ref)                     # base reaction = lambda * P_ref

    return (np.array(U), np.array(F), np.array(R_init), np.array(R_final),
            np.array(ITER), np.array(CONVERGED))

#%% ---------------------------------------
 # OPENSEESPY 

def SDOF_OPENSEES(MAT_TYPE):
    import openseespy.opensees as ops
    import numpy as np
    import PLOT_FORCE_DISP_FUN as S01
    import ANALYSIS_FUNCTION as S02
    
    # FREE UNITS
    FY, DY = 2000.0, 0.01                 # Yield force & displacement
    FU, DU = 1.2 * FY, DY * 10            # Ultimate (peak) force & displacement

    FR1, DR1 = FY / 4,  DY * 12           # 1st residual / post-peak point
    FR2, DR2 = FY / 10, DY * 15           # 2nd residual / post-peak point
    
    Force = np.array([FY, FU, FR1, FR2])
    Displacement = np.array([DY, DU, DR1, DR2])
    
    S01.PLOT_FORCE_DISP_FUN(Displacement, Force)
    
    ANA_STEPS    = 500 # Number of Displacement Increments for OpenSeespy Coding
    Dmax  = float(np.max(Displacement))
    u_inc = Dmax / ANA_STEPS
    nSprings = len(Coff)                # Number of Elements
    
    MAX_ITERATIONS = 500     # Convergence iteration for test
    MAX_TOLERANCE = 1.0e-6   # Convergence tolerance for test
    #%% ---------------- Segment slopes ----------------
    #%% ---------------- Model ----------------
    # Reset model
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)

    ops.node(1, 0.0)      # fixed support
    ops.node(2, 0.0)      # loaded node
    ops.fix(1, 1)

    # Each spring j = parallel zeroLength element with MultiLinear backbone
    for j in range(nSprings):
        matTag = 1 + j          # keep separate from element tags
        eleTag = matTag
        if MAT_TYPE == 'ELASTIC':
            kke = Coff[j] * (Force[0] / Displacement[0])   # initial stiffness
            ops.uniaxialMaterial('Elastic', matTag, kke)
    
        elif MAT_TYPE == 'INELASTIC':
            pos_pts, neg_pts = [], []
            for f, d in zip(Force, Displacement):
                pos_pts.extend([Coff[j] * f,  d])
                neg_pts.extend([-Coff[j] * f, -d])
    
            ops.uniaxialMaterial('HystereticSM', matTag,
                                 *pos_pts, *neg_pts,
                                 1.0, 1.0, 0.0, 0.0, 0.0)   # pinchX,pinchY,damage1,damage2,beta
    
        ops.element('zeroLength', eleTag, 1, 2, '-mat', matTag, '-dir', 1)

    #%% ---------------- Loading & Analysis ----------------
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(2, 1.0)      # reference load (direction of imposed disp.)

    ops.system('BandGeneral')
    ops.numberer('Plain')
    ops.constraints('Plain')
    ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/normDispIncr.html
    ops.integrator('DisplacementControl', 2, 1, u_inc)
    if MAT_TYPE == 'INELASTIC':
        ops.algorithm('Newton')  # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/algorithm.html
    if MAT_TYPE == 'ELASTIC':
        ops.algorithm('Linear')    
    ops.analysis('Static')

    #%% ---------------- Pushover Loop ----------------
    U_os = [0.0]
    F_os = [0.0]

    print('#################################################')
    print('#  OpenSeesPy : Pushover (Displacement Control) #')
    print('#################################################')

    for i in range(1, ANA_STEPS + 1):
        OK = ops.analyze(1)
        S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS)

        u = ops.nodeDisp(2, 1)
        ops.reactions()
        Ftot = ops.nodeReaction(1, 1)       

        U_os.append(u)                      # DISPLACEMENT 
        F_os.append(Ftot)                   # REACTION

        print(f'(+) Step {i:>6d}  u = {u:8.5f}  F_total = {Ftot:10.5f}')
        
    return U_os, F_os

#%% ----------------------------------------------------------------------------
# MAIN

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    
    # FREE UNITS
    FY, DY = 2000.0, 0.01                 # Yield force & displacement
    FU, DU = 1.2 * FY, DY * 10            # Ultimate (peak) force & displacement

    FR1, DR1 = FY / 4,  DY * 12           # 1st residual / post-peak point
    FR2, DR2 = FY / 10, DY * 15           # 2nd residual / post-peak point
    
    Force = np.array([FY, FU, FR1, FR2])
    Displacement = np.array([DY, DU, DR1, DR2])
    
    """
    Coff is an array of stiffness scaling coefficients (one per parallel spring).
     Each value multiplies the shared multilinear force–displacement backbone,
     so the springs can have different initial stiffnesses and force capacities
     while keeping the same shape of the nonlinear curve. In the analytical solution
     the coefficients scale the segment stiffnesses; in OpenSeesPy they scale the
     material points of each zeroLength spring. When all entries equal 1.0 the springs
     are identical; different values produce a set of parallel springs with proportional
     strength and stiffness.
    """
    Coff = np.array([1.0, 1.0, 1.0, 1.0, 1.0])   # STIFFNESS COEFFICIENTS

    ANA_STEPS, MAX_ITERATIONS, TOL = 1000, 50, 1.0e-10

    #%% RUN NUMERICAL INELASTIC ANALYSIS
    D1, F1, R0, R1, IT1, OK1 = SDOF_NUMERICAL_PUSHOVER_NEWTON(Displacement, Force, Coff,
                                                ANA_STEPS, MAX_ITERATIONS, TOL)
    
    #%% RUN OPENSEES INELASTIC ANALYSIS
    MAT_TYPE = 'INELASTIC'    
    U_os, F_os = SDOF_OPENSEES(MAT_TYPE) # RUN THE FUNCTION

    U_os = np.array(U_os)
    F_os = np.array(F_os)

    print(f'## Analysis finished at u = {U_os[-1]:.5f}, F = {F_os[-1]:.5f}')

    #%% RUN OPENSEES ELASTIC ANALYSIS
    MAT_TYPE = 'ELASTIC'    
    U_el, F_el = SDOF_OPENSEES(MAT_TYPE) # RUN THE FUNCTION
    
    U_el = np.array(U_el)
    F_el = np.array(F_el)
    
    print(f'## Analysis finished at u = {U_el[-1]:.5f}, F = {F_el[-1]:.5f}')
    
    F_an_i = np.interp(U_os, D1, F1)
    abs_err = np.abs(F_os - F_an_i)
    rel_err = abs_err[1:] / np.maximum(np.abs(F_an_i[1:]), 1e-12)
    print(f"max abs err = {abs_err.max():.3e}, max rel err = {rel_err.max():.3e}")
    
    # --- elastic comparison: same solver, single straight backbone ---
    Dmax = Displacement.max()
    Ke = Force[0] / Displacement[0]
    D2, F2, R0e, R1e, IT2, OK2 = SDOF_NUMERICAL_PUSHOVER_NEWTON([10 * Dmax], [Ke * 10 * Dmax], Coff,
                                                  ANA_STEPS, MAX_ITERATIONS, TOL, Dmax=Dmax)
    plt.figure(1, figsize=(12, 8))
    plt.plot(abs(D1), abs(F1), color='black', linewidth=4, label='Analytical - Inelastic')
    #plt.plot(abs(D2), abs(F2), color='blue', linestyle='--', linewidth=2, label='Analytical - Elastic')
    plt.plot(abs(U_os), abs(F_os),  color='red',  linestyle='--', linewidth=4, label='OpenSeesPy - Inelastic')
    #plt.plot(abs(U_el), abs(F_el),  color='purple', linestyle='--', linewidth=2, label='OpenSeesPy - Elastic')
    plt.grid(True)
    plt.xlabel('Displacement')
    plt.ylabel('Base Reaction')
    plt.title('Base Reaction–Displacement : Analytical vs OpenSeesPy', color='b')
    plt.legend()
    plt.tight_layout()
    plt.show()