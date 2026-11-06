###################################################################################################
#                >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                 #
#   DISPLACEMENT-CONTROLLED NEWTON-RAPHSON MONOTONIC PUSHOVER OF PARALLEL NONLINEAR SPRINGS       #
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
In the function  SDOF_NUMERICAL_PUSHOVER_NEWTON():   
This function pushes a spring system sideways in small, equal displacement steps and
 records the force needed at each step. It first builds the backbone curve (adding the origin)
 and combines all parallel springs into one by summing `Coff`. Then, for each step, it 'predicts'
 the new force by moving along the current tangent stiffness, 'checks the residual'
 (the gap between that predicted force and the true backbone force at the new displacement),
 and 'corrects' the force by subtracting the residual until the gap falls below the tolerance.
 If a step converges, it 'commits' the new displacement, force, and tangent as the starting 
 point for the next step; if not, the analysis stops. The result is the force-displacement history,
 plus diagnostics (residual before and after correction, iteration counts, and convergence flags)
 for each increment. The residual is nonzero only when a step crosses a kink in the backbone,
 because the predictor then uses the old segment's slope for part of a step that belongs to
 the next segment.
"""
#%% ---------------------------------------
import matplotlib.pyplot as plt
import numpy as np

def SDOF_NUMERICAL_PUSHOVER_NEWTON(Displacement, Force, Coff, ANA_STEPS,
                                   MAX_ITERATIONS, TOL, Dmax=None):
    # Backbone including the origin
    D = np.concatenate(([0.0], Displacement))
    F = np.concatenate(([0.0], Force))
    slopes = np.diff(F) / np.diff(D)          # tangent of each segment

    Csum  = np.sum(Coff)                      # parallel springs -> one scaled spring
    Dmax  = D[-1] if Dmax is None else Dmax
    dU    = Dmax / ANA_STEPS
    F_ref = Csum * F.max()                    # force scale for the residual

    def INTERNAL_FORCE(u):
        """
        INTERNAL_FORCE(u): Returns the total spring force at displacement u by
        linearly interpolating the backbone curve (D, F) and multiplying by Csum
        to account for all parallel springs, holding the force constant beyond the
        last backbone point.
        """
        return Csum * np.interp(u, D, F)      # flat after the last point
    def TANGENT(u):
        """
        TANGENT(u): Returns the current stiffness (slope) at u by finding which backbone
        segment u lies in (searchsorted), taking that segment’s slope, and scaling it by
        Csum, or returning 0 if u is past the last point where the curve is flat.
        """
        if u >= D[-1]:
            return 0.0
        seg = np.searchsorted(D, u, side='right') - 1
        return Csum * slopes[seg]

    # Initial state
    u, lam = 0.0, 0.0
    Kt = TANGENT(u)
    U, Fb = [0.0], [0.0]
    R_init, R_final, ITER, CONVERGED = [], [], [], []

    for i in range(1, ANA_STEPS + 1):
        u_new = dU * i

        # 1) Predictor: committed tangent
        lam_new = lam + Kt * (u_new - u)

        # 2) Residual (unbalanced force)
        R = lam_new - INTERNAL_FORCE(u_new)
        R_init.append(abs(R) / F_ref)

        # 3) Newton corrector at fixed displacement
        it = 0
        while abs(R) / F_ref > TOL and it < MAX_ITERATIONS:
            lam_new -= R                       # K_T cancels: d_lam = -R / P_ref
            R = lam_new - INTERNAL_FORCE(u_new)
            it += 1

        ok = abs(R) / F_ref <= TOL
        R_final.append(abs(R) / F_ref)
        ITER.append(it)
        CONVERGED.append(ok)

        if not ok:
            print(f'(-) Increment {i}: NOT converged, |R|/Fref = {abs(R)/F_ref:.3e}')
            break

        # 4) Commit
        u, lam, Kt = u_new, lam_new, TANGENT(u_new)
        U.append(u)
        Fb.append(lam)

    return (np.array(U), np.array(Fb), np.array(R_init),
            np.array(R_final), np.array(ITER), np.array(CONVERGED))

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
