###################################################################################################
#                >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                 #
#          NONLINEAR DYNAMIC (TIME-HISTORY) ANALYSIS OF PARALLEL NONLINEAR SDOF SPRINGS           #
#-------------------------------------------------------------------------------------------------#
#                  FINITE ELEMENT METHOD IN PYTHON AND OPENSEESPY VERIFICATIONS                   #
#   NEWMARK-BETA + NEWTON-RAPHSON  vs.  OPENSEESPY (HystereticSM, Transient analysis)             #
#-------------------------------------------------------------------------------------------------#
#                THIS PYTHON SCRIPT IS WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)              # 
#                               E-mail:salar.d.ghashghaei@gmail.com                               #
###################################################################################################
"""
 This script runs a nonlinear time-history analysis of a single-degree-of-freedom system,
 made of five parallel hysteretic springs that share one multilinear backbone, and checks
 a hand-written solver against OpenSeesPy. The system obeys m·u'' + c·u' + Fs(u) = −m·ag(t),
 with mass chosen for a 0.5 s initial period, 5% mass-proportional damping, and the scaled
 ground acceleration record as excitation. In the numerical solver,
 SDOF_NUMERICAL_DYNAMIC_NEWMARK marches through the record with Newmark-β
 (average acceleration). At each step it guesses the displacement, computes the implied
 acceleration and velocity, and calls SPRING to get the spring force. It then forms the
 dynamic residual and corrects the displacement with Newton-Raphson, using the effective
 tangent k + γc/(βΔt) + m/(βΔt²), until the normalized residual drops below the tolerance.
 SPRING is a pure function of the committed state and the trial displacement, so iterating
 never corrupts the history. It stores the last reversal point, the largest peaks on each
 side and the start of the last unloading, mirrors the motion so it always looks like it
 goes toward +x, and picks a branch. The branches are the backbone beyond the previous peak,
 unloading along k₀, retracing k₀ after a partial unloading, or reloading toward the
 previous peak. Once a step converges, the state is committed. The solver returns time,
 displacement, velocity, acceleration, force with the sign flipped to match OpenSees'
 reaction convention, iteration counts and residuals. SDOF_OPENSEES_DYNAMIC builds the
 same system with HystereticSM zeroLength springs, a nodal mass, rayleigh damping equal
 to 2ξω₀m, a Path time series with UniformExcitation, and a Newton/Newmark Transient
 analysis. It falls back on the ANALYSIS helper's alternative tests and algorithms if
 a step fails. The main block defines the backbone and the record, runs both solvers,
 prints the maximum absolute and relative errors in displacement and force, and plots
 the ground motion, displacement history, hysteresis loops and absolute acceleration.
 The two solvers agree to roundoff, which confirms that the hysteresis rule and the
 Newmark-Newton scheme reproduce OpenSees.
 
 Equation of motion (relative displacement u, ground acceleration ag):

        m*u'' + c*u' + Fs(u, history) = -m*ag(t)

 NUMERICAL: Newmark-beta (default: average acceleration, gamma=0.5, beta=0.25), solved with a full
 Newton-Raphson iteration at every time step. The hysteretic spring (SPRING) is exactly the same
 peak-oriented rule used in the cyclic pushover script; it is a pure function of (committed state,
 trial displacement), so it can be re-evaluated during iterations without corrupting the history.

 OPENSEES: zeroLength springs with HystereticSM, nodal mass, UniformExcitation with a Path
 time series, mass-proportional damping (c = 2*xi*omega0*m, identical to the numerical model),
 Newmark integrator, Transient analysis.

 Units are free but must be consistent. The record file is assumed to be in units of g
 (set GM_UNIT_G = 1.0 if it is already in your acceleration unit).
"""
# BOOK: Dynamics of Structures in SI Units -  Anil Kumar Chopra
'https://share.google/TDN5O4eWmw5pH8zUH'
# BOOK: Differential Equations for Engineers-Wei-Chau Xie-CAMBRIDGE-2010

import numpy as np
import matplotlib.pyplot as plt
import PLOT_FORCE_DISP_FUN as S01

#%% ---------------------------------------------------------------
# NUMERICAL: Newmark-beta + Newton-Raphson, peak-oriented hysteresis
def SDOF_NUMERICAL_DYNAMIC_NEWMARK(Displacement, Force, Coff, ag, dt, mass, xi,
                                   MAX_ITERATIONS, MAX_TOLERANCE):
                                       
    """
    # Properties
    K0 = Σ Coff · (FY / DY)            initial stiffness
    ω0 = √(K0 / m)                     natural frequency
    T0 = 2π / ω0   →   m = K0·(T0 / 2π)²
    c  = 2·ξ·ω0·m                      mass-proportional damping (OpenSees: alphaM = 2ξω0)
    """    
    # Newmark-beta method parameters
    gamma, beta = 0.5, 0.25
    
    D = np.concatenate(([0.0], Displacement))
    F = np.concatenate(([0.0], Force))
    slopes = np.diff(F) / np.diff(D)
    k0 = slopes[0]                                # initial (unloading) stiffness
    Csum = np.sum(Coff)                           # parallel springs -> one scaled spring
    F_ref = Csum * F.max()

    K0 = Csum * k0
    omega0 = np.sqrt(K0 / mass)
    c = 2.0 * xi * omega0 * mass                  # constant (mass-proportional) damping

    def BACKBONE(x):
        """
        BACKBONE(x):
        For a non-negative displacement $x$, returns the force and tangent stiffness
        from the multilinear backbone curve (interpolating between points or returning
        the last force with zero stiffness if x exceeds the ultimate displacement).
        """
        if x >= D[-1]:
            return F[-1], 0.0
        seg = np.searchsorted(D, x, side='right') - 1
        return np.interp(x, D, F), slopes[seg]

    def SPRING(u, st):
        # Pure function: committed state + trial u -> (force, tangent, trial state)
        """
        SPRING(u, st):
        Given a trial displacement u and the committed hysteresis state st,
        detects reversals, mirrors the problem so motion is always treated as
        positive, selects the appropriate branch (backbone, unloading on k_0,
        or reloading toward the previous peak), updates peaks and state, and
        returns the scaled force, tangent, and new trial state.
        """
        # Newmark-beta method parameters
        gamma, beta = 0.5, 0.25
        
        st = dict(st)
        du = u - st['u']
        if du != 0 and np.sign(du) != st['d']:    # reversal
            st['u_a'], st['f_a'] = st['u_rev'], st['f_rev']      # point where the last unloading began
            st['u_rev'], st['f_rev'], st['d'] = st['u'], st['f'], int(np.sign(du))

        d = st['d']
        x, xr, gr = d * u, d * st['u_rev'], d * st['f_rev']
        xt, gt = (st['umax'], st['fmax']) if d == 1 else (-st['umin'], -st['fmin'])
        xz = xr - gr / k0

        if x >= xt:                               # backbone
            g, k = BACKBONE(x)
        elif x <= xz:                             # unloading on k0
            g, k = gr + k0 * (x - xr), k0
        elif gr > 0:                              # reversal during partial unloading (F still > 0)
            xa, ga = d * st['u_a'], d * st['f_a']
            if ga > gr and x <= xa:               # retrace the unloading line back to where it began
                g, k = gr + k0 * (x - xr), k0
            elif ga > gr:                         # then head for the previous peak
                k = (gt - ga) / (xt - xa)
                g = ga + k * (x - xa)
            else:                                 # fallback: reversal point -> peak
                k = (gt - gr) / (xt - xr)
                g = gr + k * (x - xr)
        else:                                     # reloading from F=0 crossing toward previous peak
            k = gt / (xt - xz)
            g = k * (x - xz)

        f = d * g
        st['u'], st['f'] = u, f
        if u > st['umax']: st['umax'], st['fmax'] = u, f
        if u < st['umin']: st['umin'], st['fmin'] = u, f
        return Csum * f, Csum * k, st

    state = dict(u=0.0, f=0.0, d=1, u_rev=0.0, f_rev=0.0, u_a=0.0, f_a=0.0,
                 umax=D[1], fmax=F[1], umin=-D[1], fmin=-F[1])

    # initial conditions (system at rest)
    u, v = 0.0, 0.0
    fs = 0.0
    a = 0.0   # OpenSees Transient also starts with a(0)=0 (no initial-acceleration solve)

    n = len(ag)
    U, V, A, FS = np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n)
    A[0] = a
    ITER, R_final = np.zeros(n, dtype=int), np.zeros(n)

    for i in range(1, n):
        """
        For each time step, it guesses the new displacement as the previous one, computes Newmark acceleration 
        and velocity from it, and evaluates the dynamic residual R = m·a + c·v + Fs(u) − (−m·ag).
        While R is above tolerance, it corrects u by −R/K_eff (a Newton-Raphson update with the effective
                                                               tangent k + γc/(βΔt) + m/(βΔt²)).

        If the residual never gets below TOL within MAX_ITERATIONS, it prints a “NOT converged”
        message and stops the analysis at that step.
        """
        p_new = -mass * ag[i]                       # effective ground-motion load
        u_new = u                                   # initial guess: previous displacement

        it = 0
        while True:
            # Newmark kinematics as functions of the trial displacement
            a_new = (u_new - u) / (beta * dt**2) - v / (beta * dt) - (1.0 / (2 * beta) - 1.0) * a # ACCELERATION
            v_new = v + dt * ((1 - gamma) * a + gamma * a_new)                                    # VELOCITY

            fs_new, k, trial = SPRING(u_new, state)
            R = mass * a_new + c * v_new + fs_new - p_new                                         # DYNAMIC RESIDUAL
            DX = abs(R) / F_ref 
            if DX <= MAX_TOLERANCE or it >= MAX_ITERATIONS:
                break
            K_eff = k + gamma / (beta * dt) * c + mass / (beta * dt**2)                           # STRUCTRE STIFFNESS
            u_new -= R / K_eff
            it += 1

        R_final[i], ITER[i] = DX, it
                
        if DX > MAX_TOLERANCE:
            print(f'(-) Step {i}: NOT converged, |R|/Fref = {abs(R)/F_ref:.3e}')
            n = i
            break

        # commit
        state, u, v, a, fs = trial, u_new, v_new, a_new, fs_new
        U[i], V[i], A[i], FS[i] = u, v, a, fs

    DATA = ((np.arange(n) * dt, U[:n], V[:n], A[:n], -FS[:n], ITER[:n], R_final[:n]))
    
    return DATA

#%% ---------------------------------------------------------------
# OPENSEES: transient analysis with uniform base excitation
def SDOF_OPENSEES_DYNAMIC(MAT_TYPE, Displacement, Force, Coff, ag, dt, mass, xi,
                          GM_FACTOR, MAX_ITERATIONS, MAX_TOLERANCE):
    
    # ag is the raw record; GM_FACTOR scales it to the model's acceleration unit
    import openseespy.opensees as ops
    import ANALYSIS_FUNCTION as S02
    
    # Reset model
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)
    
    # Nodes
    ops.node(1, 0.0)
    ops.node(2, 0.0)
    
    # Constrains
    ops.fix(1, 1)
    
    # Mass
    ops.mass(2, mass)

    for j in range(len(Coff)):
        tag = 1 + j
        if MAT_TYPE == 'ELASTIC':
            ops.uniaxialMaterial('Elastic', tag, Coff[j] * Force[0] / Displacement[0])
        if MAT_TYPE == 'INELASTIC':
            pos, neg = [], []
            for f, d in zip(Force, Displacement):
                pos.extend([Coff[j] * f, d]); neg.extend([-Coff[j] * f, -d])
            ops.uniaxialMaterial('HystereticSM', tag, *pos, *neg, 1.0, 1.0, 0.0, 0.0, 0.0)
            # INFO LINK: https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/HystereticSM.html
            
        ops.element('zeroLength', tag, 1, 2, '-mat', tag, '-dir', 1)
        # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/ZeroLength.html

    # mass-proportional damping: c = alphaM * m = 2*xi*omega0*m
    K0 = np.sum(Coff) * Force[0] / Displacement[0]
    omega0 = np.sqrt(K0 / mass)
    ops.rayleigh(2.0 * xi * omega0, 0.0, 0.0, 0.0)
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/reyleigh.html

    # ground motion
    ops.timeSeries('Path', 1, '-dt', dt, '-values', *ag, '-factor', GM_FACTOR)
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/pathTs.html
    ops.pattern('UniformExcitation', 1, 1, '-accel', 1)
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/uniformExcitation.html

    ops.system('BandGeneral')
    ops.numberer('Plain')
    ops.constraints('Plain')
    ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS)
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/normDispIncr.html

    if MAT_TYPE == 'INELASTIC':
        ops.algorithm('Newton')  # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/algorithm.html
    if MAT_TYPE == 'ELASTIC':
        ops.algorithm('Linear')
        
    ops.integrator('Newmark', 0.5, 0.25) # INFO LINK: https://opensees.github.io/OpenSeesDocumentation/user/manual/analysis/integrator/Newmark.html
    ops.analysis('Transient')            # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/analysis.html

    t_os, U_os, F_os = [0.0], [0.0], [0.0]
    for i in range(1, len(ag)):
        OK = ops.analyze(1, dt)
        if OK != 0:
            S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS)
        ops.reactions()
        t_os.append(ops.getTime())
        U_os.append(ops.nodeDisp(2, 1))                        # DISPLACEMENT
        F_os.append(ops.nodeReaction(1, 1))                    # REACTION

    DATA = (np.array(t_os), np.array(U_os), np.array(F_os))
    
    return DATA

#%% ---------------------------------------------------------------
# MAIN
if __name__ == '__main__':

    # ----- Backbone (units: N, m, kg, s) -----
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
    
    # Plot the Force-Displacemnt for Spring
    S01.PLOT_FORCE_DISP_FUN(Displacement, Force)

    # Ground motion
    GM_FILE = 'Ground_Acceleration_X.txt'
    DT = 0.01                    # <-- time step of the record (CHECK THIS!)
    GM_UNIT_G = 9.81             # record in g -> m/s^2 (use 1.0 if already in m/s^2)
    GM_SCALE = 10.0              # scale factor so the structure actually yields
    ag_raw = np.loadtxt(GM_FILE)
    GM_FACTOR = GM_UNIT_G * GM_SCALE
    ag = ag_raw * GM_FACTOR

    # Dynamic properties
    T0 = 0.5                                                   # target initial period (s)
    K0 = np.sum(Coff) * Force[0] / Displacement[0]
    MASS = K0 * (T0 / (2 * np.pi))**2                          # m = K0 / omega0^2
    XI = 0.05                                                  # damping ratio

    print(f'K0 = {K0:.3e},  m = {MASS:.3f},  T0 = {T0:.3f} s,  N = {len(ag)},  dt = {DT}')

    #%% DEFINE ANALYSIS PROPERTIES
    MAX_ITERATIONS = 50    # Convergence iteration for test
    MAX_TOLERANCE = 1.0e-6   # Convergence tolerance for test
    
    #%% RUN NUMERICAL METHOD
    t1, U1, V1, A1, F1, IT, RES = SDOF_NUMERICAL_DYNAMIC_NEWMARK(
        Displacement, Force, Coff, ag, DT, MASS, XI, MAX_ITERATIONS=50, TOL=1.0e-10)

    #%% RUN OPENSEES
    t_os, U_os, F_os = SDOF_OPENSEES_DYNAMIC('INELASTIC', Displacement, Force, Coff,
                                             ag_raw, DT, MASS, XI, GM_FACTOR)

    n = min(len(U1), len(U_os))
    errU = np.abs(U_os[:n] - U1[:n]); errF = np.abs(F_os[:n] - F1[:n])
    print(f'max |U| = {np.abs(U1).max():.4e}   max |Fs| = {np.abs(F1).max():.4e}')
    print(f'Disp  : max abs err = {errU.max():.3e}, max rel err = {errU.max()/np.abs(U1).max():.3e}')
    print(f'Force : max abs err = {errF.max():.3e}, max rel err = {errF.max()/np.abs(F1).max():.3e}')
    print(f'Max Newton iterations per step = {IT.max()}')

    # ----- Plots -----
    fig, ax = plt.subplots(2, 2, figsize=(15, 9))
    ax[0, 0].plot(np.arange(len(ag)) * DT, ag, 'g')
    ax[0, 0].set_title('Ground acceleration')
    ax[0, 0].set_xlabel('Time [s]')
    ax[0, 0].set_ylabel('ag [m/s^2]')
    ax[0, 1].plot(t1, U1, 'k', lw=3, label='Numerical')
    ax[0, 1].plot(t_os, U_os, 'r--', lw=2, label='OpenSeesPy')
    ax[0, 1].set_title('Displacement time history')
    ax[0, 1].set_xlabel('Time [s]')
    ax[0, 1].set_ylabel('Displacement [m]')
    ax[0, 1].legend()
    ax[1, 0].plot(U1, F1, 'k', lw=3, label='Numerical')
    ax[1, 0].plot(U_os, F_os, 'r--', lw=2, label='OpenSeesPy')
    ax[1, 0].set_title('Hysteresis loops')
    ax[1, 0].set_xlabel('Displacement [m]')
    ax[1, 0].set_ylabel('Base Reaction [N]')
    ax[1, 0].legend()
    ax[1, 1].plot(t1, A1 + ag[:len(A1)], 'b')
    ax[1, 1].set_title('Absolute acceleration (numerical) [m/s^2]')
    ax[1, 1].set_xlabel('Time [s]')
    for a_ in ax.ravel(): a_.grid(True)
    plt.tight_layout()
    plt.savefig('SDOF_DYNAMIC_TIME_HISTORY.png', dpi=150)
    plt.show()