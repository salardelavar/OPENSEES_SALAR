######################################################################################################
#                >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                    #
#   CYCLIC DISPLACEMENT-CONTROLLED NEWTON-RAPHSON PUSHOVER OF A MULTI-STORY NONLINEAR SHEAR BUILDING #
#----------------------------------------------------------------------------------------------------#
#                  FINITE ELEMENT METHOD IN PYTHON AND OPENSEESPY VERIFICATIONS                      #
#----------------------------------------------------------------------------------------------------#
#                   THIS PROGRAM IS WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                    # 
#                               E-mail:salar.d.ghashghaei@gmail.com                                  #
######################################################################################################
"""
 This python script drives the same parallel-spring MDOF system through a back-and-forth displacement
 history (0 → +a → −a → 0 at growing amplitudes) using two independent approaches and compares them.
 DISP_HISTORY builds the list of target displacements step by step. The numerical solver
 MDOF_NUMERICAL_CYCLIC_NEWTON works like the monotonic version (predict with the last tangent,
 compute the residual, correct the force until balanced), but because force now depends on loading history,
 the SPRING function remembers the last reversal point and the largest peaks on each side.
 It then decides the branch: follow the backbone when pushing past the previous peak, unload along 
 the initial stiffness after a reversal, or reload in a straight line toward the previous peak on the
 opposite side. MDOF_OPENSEES_CYCLIC builds the equivalent model with HystereticSM springs and re-issues
 the DisplacementControl increment at every step so the node lands exactly on each target displacement.
 Finally, the main block runs both, prints the maximum absolute and relative error, and plots the two hysteresis
 loops together with the displacement history, so you can check that the hand-written hysteresis rule reproduces OpenSees.
"""
import numpy as np
import matplotlib.pyplot as plt
import PLOT_FORCE_DISP_FUN as S01

#%% ---------------------------------------------------------------
# CYCLIC DISPLACEMENT HISTORY (0 -> +a -> -a -> 0 for every amplitude)
def DISP_HISTORY(amplitudes, n_leg=50):
    legs = [0.0]
    for a in amplitudes:
        for p in (a, -a, 0.0):
            legs.extend(np.linspace(legs[-1], p, n_leg + 1)[1:])
    return np.array(legs)

#%% ---------------------------------------------------------------
# LATERAL LOAD PATTERN (normalised so that sum(P_ref) = 1)
def LOAD_PATTERN(N, kind='TRIANGULAR', heights=None, masses=None):
    h = np.arange(1, N + 1, dtype=float) if heights is None else np.asarray(heights, float)
    m = np.ones(N) if masses is None else np.asarray(masses, float)
    if kind == 'TRIANGULAR':
        P = m * h
    elif kind == 'UNIFORM':
        P = m.copy()
    else:
        raise ValueError("kind must be 'TRIANGULAR' or 'UNIFORM'")
    return P / P.sum()

#%% ---------------------------------------------------------------
# ONE STORY SPRING: peak-oriented hysteresis (same rule as the SDOF code, one instance per story)
def MAKE_HYSTERESIS(Displacement, Force):
    D = np.concatenate(([0.0], Displacement))
    F = np.concatenate(([0.0], Force))
    slopes = np.diff(F) / np.diff(D)
    k0 = slopes[0]                                # unloading stiffness (beta = 0)

    def BACKBONE(x):                              # x >= 0, flat after last point
        if x >= D[-1]:
            return F[-1], 0.0
        seg = np.searchsorted(D, x, side='right') - 1
        return np.interp(x, D, F), slopes[seg]

    def INITIAL_STATE():
        # initial peaks = yield point, as in OpenSees Hysteretic
        return dict(u=0.0, f=0.0, d=1, u_rev=0.0, f_rev=0.0,
                    umax=D[1], fmax=F[1], umin=-D[1], fmin=-F[1])

    def SPRING(u, st):
        # Pure function: committed state + trial drift u -> (unit force, unit tangent, trial state)
        st = dict(st)
        du = u - st['u']
        if du != 0 and np.sign(du) != st['d']:    # reversal: remember where it happened
            st['u_rev'], st['f_rev'], st['d'] = st['u'], st['f'], int(np.sign(du))

        d = st['d']                               # mirror so we always move toward +x
        x, xr, gr = d * u, d * st['u_rev'], d * st['f_rev']
        xt, gt = (st['umax'], st['fmax']) if d == 1 else (-st['umin'], -st['fmin'])
        xz = xr - gr / k0                         # where the unloading line hits F = 0

        if x >= xt:                               # backbone
            g, k = BACKBONE(x)
        elif x <= xz:                             # unloading, slope k0
            g, k = gr + k0 * (x - xr), k0
        else:                                     # reloading toward the previous peak
            k = gt / (xt - xz)
            g = k * (x - xz)

        f = d * g
        st['u'], st['f'] = u, f
        if u > st['umax']: st['umax'], st['fmax'] = u, f
        if u < st['umin']: st['umin'], st['fmin'] = u, f
        return f, k, st

    return SPRING, INITIAL_STATE, k0

#%% ---------------------------------------------------------------
# NUMERICAL: displacement-controlled Newton-Raphson on the roof, N-story shear building
def MDOF_NUMERICAL_CYCLIC_NEWTON(Displacement, Force, STORY_COFF, P_ref, U_target,
                                 MAX_ITERATIONS=100, TOL=1.0e-8):
    """
    STORY_COFF[i] : scale of the shared multilinear backbone for story i (stiffness AND strength)
                    (= sum of the parallel-spring coefficients of that story)
    P_ref         : reference lateral load vector (one value per floor), load = lam * P_ref
    U_target      : roof displacement history (control DOF = top floor)
    """
    scale = np.asarray(STORY_COFF, float)
    N = len(scale)
    P_ref = np.asarray(P_ref, float)
    c = N - 1                                           # control DOF = roof
    SPRING, INITIAL_STATE, k0 = MAKE_HYSTERESIS(Displacement, Force)
    F_ref = scale.max() * np.max(Force)                 # normalisation of the residual norm
    k_floor = 1.0e-8 * k0 * scale.min()                 # keeps K_t non-singular on flat branches

    def ASSEMBLE(k_story):                              # tridiagonal tangent stiffness matrix
        """
        ASSEMBLE(k_story) builds the N×N tridiagonal tangent stiffness matrix of the shear
        building by adding each story’s tangent k_i to its two connected floor DOFs (+k_i on
        the diagonal of floors i−1 and i, −k_i on the off-diagonals between them),
        after flooring near-zero tangents so the matrix stays solvable.
        """
        k_story = np.where(np.abs(k_story) < k_floor, k_floor, k_story)
        K = np.zeros((N, N))
        for i in range(N):
            K[i, i] += k_story[i]
            if i + 1 < N:
                K[i, i] += k_story[i + 1]
                K[i, i + 1] -= k_story[i + 1]
                K[i + 1, i] -= k_story[i + 1]
        return K

    def EVALUATE(u_vec, states):                        # story drifts -> shears, tangents, trial states
        """
        EVALUATE(u_vec, states) takes a trial set of floor displacements and returns what
        the structure’s springs would do at that displacement, without committing anything.
        """
        drift = np.diff(np.concatenate(([0.0], u_vec)))
        """
        Story drifts: it prepends the ground (0) to u_vec and takes differences,
        so drift[i] = u_i − u_{i−1}
        """
        V, k, trial = np.zeros(N), np.zeros(N), []
        for i in range(N):
            f, ki, ti = SPRING(drift[i], states[i])
            """
            Spring call per story: each story’s drift goes into SPRING(drift[i], states[i]),
            which uses that story’s committed hysteresis state to decide the branch (backbone,
            unloading, or reloading). It returns a unit force f, unit tangent ki, and a trial state ti
            """
            V[i], k[i] = scale[i] * f, scale[i] * ki
            """
            Scaling: the force and tangent are multiplied by scale[i] (that story’s STORY_COFF),
            giving the story shear V[i] and tangent stiffness k[i].
            """
            trial.append(ti)
        f_int = V - np.concatenate((V[1:], [0.0]))      # internal force at each floor
        """
        Internal force at floors: f_int = V − [V[1:], 0] is the vertical equilibrium of
        each floor: floor i feels the shear of the story below pushing it back (V[i]) minus
        the shear of the story above pulling it forward (V[i+1]). The roof has no story above,
        so its second term is 0.
        """
        return V, k, f_int, trial

    # committed state
    states = [INITIAL_STATE() for _ in range(N)]
    u, lam = np.zeros(N), 0.0
    kt = scale * k0

    U_floor, Lam, V_story = [u.copy()], [0.0], [np.zeros(N)]
    R_init, R_final, ITER, CONVERGED = [], [], [], []

    for i in range(1, len(U_target)):
        # 1) predictor: move the roof to the target with the last committed tangent
        K = ASSEMBLE(kt)
        dUp = np.linalg.solve(K, P_ref)
        dlam = (U_target[i] - u[c]) / dUp[c]
        u_new, lam_new = u + dlam * dUp, lam + dlam

        # 2) residual
        V, k, f_int, trial = EVALUATE(u_new, states)
        R = lam_new * P_ref - f_int
        R_init.append(np.linalg.norm(R) / F_ref)

        # 3) Newton corrector (roof displacement stays fixed at the target)
        it = 0
        while np.linalg.norm(R) / F_ref > TOL and it < MAX_ITERATIONS:
            K = ASSEMBLE(k)
            dUp, dUr = np.linalg.solve(K, P_ref), np.linalg.solve(K, R)
            dlam = -dUr[c] / dUp[c]
            u_new, lam_new = u_new + dlam * dUp + dUr, lam_new + dlam
            V, k, f_int, trial = EVALUATE(u_new, states)
            R = lam_new * P_ref - f_int
            it += 1

        ok = np.linalg.norm(R) / F_ref <= TOL
        R_final.append(np.linalg.norm(R) / F_ref)
        ITER.append(it)
        CONVERGED.append(ok)

        if not ok:
            print(f'(-) Step {i}: NOT converged, |R|/Fref = {np.linalg.norm(R) / F_ref:.3e}')
            break

        # 4) commit
        states, u, lam, kt = trial, u_new, lam_new, k
        U_floor.append(u.copy()); Lam.append(lam); V_story.append(V.copy())

    U_floor = np.array(U_floor)                         # (steps, N) floor displacements
    drift = np.diff(np.hstack((np.zeros((len(U_floor), 1)), U_floor)), axis=1)
    DATA = dict(U_floor=U_floor, drift=drift,
                base_reaction=-np.array(Lam) * P_ref.sum(),   # same sign as OpenSees base reaction
                V_story=np.array(V_story),
                R_init=np.array(R_init), R_final=np.array(R_final),
                ITER=np.array(ITER), CONVERGED=np.array(CONVERGED))
    
    return DATA

#%% ---------------------------------------------------------------
# OPENSEES: N-story shear building, DisplacementControl on the roof re-issued every step
def MDOF_OPENSEES_CYCLIC(MAT_TYPE, Displacement, Force, STORY_COFF, P_ref, U_target,
                         MAX_ITERATIONS=500, MAX_TOLERANCE=1.0e-6):
    import openseespy.opensees as ops
    import ANALYSIS_FUNCTION as S02

    N = len(STORY_COFF)
    # Reset model
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)

    # Nodes: 0 = ground, 1..N = floors
    for n in range(N + 1):
        ops.node(n, 0.0)
    ops.fix(0, 1)

    # Story i = zeroLength element between node i-1 and node i
    for i in range(1, N + 1):
        s = STORY_COFF[i - 1]
        if MAT_TYPE == 'ELASTIC':
            ops.uniaxialMaterial('Elastic', i, s * Force[0] / Displacement[0])
        if MAT_TYPE == 'INELASTIC':
            pos, neg = [], []
            for f, d in zip(Force, Displacement):
                pos.extend([s * f, d]); neg.extend([-s * f, -d])
            ops.uniaxialMaterial('HystereticSM', i, *pos, *neg, 1.0, 1.0, 0.0, 0.0, 0.0)
        ops.element('zeroLength', i, i - 1, i, '-mat', i, '-dir', 1)

    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    for i in range(1, N + 1):
        ops.load(i, float(P_ref[i - 1]))

    ops.system('BandGeneral')
    ops.numberer('Plain')
    ops.constraints('Plain')
    ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/normDispIncr.html
    if MAT_TYPE == 'INELASTIC':
        ops.algorithm('Newton')  # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/algorithm.html
    if MAT_TYPE == 'ELASTIC':
        ops.algorithm('Linear')
    ops.analysis('Static')

    U_os, F_os = [np.zeros(N)], [0.0]
    for i in range(1, len(U_target)):
        du = U_target[i] - ops.nodeDisp(N, 1)
        ops.integrator('DisplacementControl', N, 1, du)
        OK = ops.analyze(1)
        if OK != 0:
            S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS)
        ops.reactions()
        U_os.append([ops.nodeDisp(n, 1) for n in range(1, N + 1)])   # FLOOR DISPLACEMENTS
        F_os.append(ops.nodeReaction(0, 1))                          # BASE REACTION

    U_floor = np.array(U_os)
    drift = np.diff(np.hstack((np.zeros((len(U_floor), 1)), U_floor)), axis=1)
    F_os = np.array(F_os)
    # proportional loading -> story shear_i = lam * sum(P_ref[i:]), lam = -F_base / sum(P_ref)
    V_story = (-F_os[:, None] / np.sum(P_ref)) * np.cumsum(P_ref[::-1])[::-1][None, :]
    return dict(U_floor=U_floor, drift=drift, base_reaction=F_os, V_story=V_story)

#%% ---------------------------------------------------------------
# MAIN
if __name__ == '__main__':

    # FREE UNITS
    # BACKBONE (per unit story, in story-drift terms)
    FY, DY = 2000.0, 0.01                 # Yield force & drift
    FU, DU = 1.2 * FY, DY * 10            # Ultimate (peak) force & drift
    FR1, DR1 = FY / 4,  DY * 12           # 1st residual / post-peak point
    FR2, DR2 = FY / 10, DY * 15           # 2nd residual / post-peak point

    Force = np.array([FY, FU, FR1, FR2])
    Displacement = np.array([DY, DU, DR1, DR2])

    # BUILDING
    N_STORY = 10
    # Story scale factors (ground story strongest). Several parallel springs per story could be
    # lumped here as sum(Coff_story).
    STORY_COFF = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
    P_ref = LOAD_PATTERN(N_STORY, 'TRIANGULAR')   # or 'UNIFORM'; add heights=..., masses=...

    S01.PLOT_FORCE_DISP_FUN(Displacement, Force)

    # ROOF displacement history
    U_target = DISP_HISTORY(DY * np.array([1, 2, 4, 6, 8, 12]), n_leg=50)

    # RUN NUMERICAL METHOD
    NUM = MDOF_NUMERICAL_CYCLIC_NEWTON(Displacement, Force, STORY_COFF, P_ref, U_target,
                                       MAX_ITERATIONS=100, TOL=1.0e-8)
    print(f'Max Newton iterations in one step: {NUM["ITER"].max()}')

    # RUN OPENSEES
    OS = MDOF_OPENSEES_CYCLIC('INELASTIC', Displacement, Force, STORY_COFF, P_ref, U_target)

    n = min(len(NUM['base_reaction']), len(OS['base_reaction']))
    err = np.abs(OS['base_reaction'][:n] - NUM['base_reaction'][:n])
    print(f'Base shear: max abs err = {err.max():.3e}, '
          f'max rel err = {err.max() / np.abs(NUM["base_reaction"]).max():.3e}')
    err_d = np.abs(OS['drift'][:n] - NUM['drift'][:n]).max()
    print(f'Story drifts: max abs err = {err_d:.3e}')

    # PLOTS
    fig, ax = plt.subplots(2, 2, figsize=(15, 10))
    roof_num, roof_os = NUM['U_floor'][:, -1], OS['U_floor'][:, -1]
    ax[0, 0].plot(roof_num, NUM['base_reaction'], 'k', lw=4, label='Numerical')
    ax[0, 0].plot(roof_os, OS['base_reaction'], 'r--', lw=3, label='OpenSeesPy')
    ax[0, 0].set_xlabel('Roof displacement'); ax[0, 0].set_ylabel('Base reaction')
    ax[0, 0].set_title('Global hysteresis loops'); ax[0, 0].grid(True); ax[0, 0].legend()

    ax[0, 1].plot(U_target, 'b'); ax[0, 1].set_xlabel('Step'); ax[0, 1].set_ylabel('Roof displacement')
    ax[0, 1].set_title('Cyclic roof displacement history'); ax[0, 1].grid(True)

    colors = plt.cm.viridis(np.linspace(0, 0.9, N_STORY))
    for i in range(N_STORY):
        ax[1, 0].plot(NUM['drift'][:, i], NUM['V_story'][:, i], color=colors[i], lw=3,
                      label=f'Story {i + 1}')
        ax[1, 0].plot(OS['drift'][:n, i], OS['V_story'][:n, i], 'r--', lw=1.2)
    ax[1, 0].set_xlabel('Story drift'); ax[1, 0].set_ylabel('Story shear')
    ax[1, 0].set_title('Story hysteresis (solid: numerical, red dashed: OpenSeesPy)')
    ax[1, 0].grid(True); ax[1, 0].legend()

    # envelope of peak story drifts (damage / soft-story check)
    ax[1, 1].plot(np.abs(NUM['drift']).max(axis=0) / DY, np.arange(1, N_STORY + 1), 'ko-', lw=2, ms=8,
                  label='Numerical')
    ax[1, 1].plot(np.abs(OS['drift']).max(axis=0) / DY, np.arange(1, N_STORY + 1), 'r--s', lw=2,
                  label='OpenSeesPy')
    ax[1, 1].set_xlabel('Peak story drift / yield drift'); ax[1, 1].set_ylabel('Story')
    ax[1, 1].set_title('Peak drift profile'); ax[1, 1].grid(True); ax[1, 1].legend()
    plt.tight_layout()
    plt.show()
