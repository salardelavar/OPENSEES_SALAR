###################################################################################################
#                >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                 #
#     NONLINEAR DYNAMIC (TIME-HISTORY) ANALYSIS OF A MULTI-STORY NONLINEAR SHEAR BUILDING         #
#-------------------------------------------------------------------------------------------------#
#                  FINITE ELEMENT METHOD IN PYTHON AND OPENSEESPY VERIFICATIONS                   #
#   NEWMARK-BETA + NEWTON-RAPHSON  vs.  OPENSEESPY (HystereticSM, Transient analysis)             #
#-------------------------------------------------------------------------------------------------#
#             THIS PYTHON SCRIPT IS WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                 #
#                               E-mail:salar.d.ghashghaei@gmail.com                               #
###################################################################################################
"""
 This Python script performs nonlinear dynamic time-history analysis of an N-story shear
 building whose stories are modeled as peak-oriented hysteretic springs acting on
 inter-story drifts. It assembles the initial elastic stiffness and mass matrices,
 computes modal properties and Rayleigh damping coefficients calibrated to a target
 damping ratio in two selected modes, and integrates the equations of motion with the
 average-acceleration Newmark-β method combined with full Newton-Raphson iteration at
 every time step; an equivalent OpenSeesPy model that uses zero-length elements with
 HystereticSM materials, the same Rayleigh damping, UniformExcitation, and Newmark
 integrator is run for verification. After loading a ground-acceleration record, scaling
 the floor masses to achieve a prescribed fundamental period, and executing both solvers,
 the script reports peak responses, displacement and force errors between the two
 implementations, and generates diagnostic plots of ground motion, roof displacement
 histories, global and story-level hysteresis loops, peak-drift profiles, and absolute roof
 acceleration.
 
 N-story shear building: floor i has one lateral DOF u_i (relative to the ground), story i is a
 nonlinear hysteretic spring acting on the story drift d_i = u_i - u_(i-1)   (u_0 = 0).

 Equation of motion (relative displacements, uniform ground acceleration ag):

        M*u'' + C*u' + f_int(u, history) = -M*r*ag(t)          r = {1,...,1}

 NUMERICAL: Newmark-beta (default: average acceleration, gamma = 0.5, beta = 0.25) with a full
 Newton-Raphson iteration at every step. The effective tangent is the assembled tridiagonal matrix

        K_eff = K_t + gamma/(beta*dt) * C + 1/(beta*dt^2) * M

 Every story owns its own committed hysteresis state (peak-oriented rule with partial-unloading
 handling, identical to the SDOF time-history script). SPRING is a pure function of
 (committed state, trial drift), so it can be re-evaluated during iterations without corrupting history.

 DAMPING: Rayleigh, C = a0*M + a1*K0 with K0 = INITIAL stiffness (constant in time), calibrated so that
 the damping ratio equals xi in two chosen modes. OpenSees equivalent: ops.rayleigh(a0, 0, a1, 0)
 (third argument = betaKinit), so both models use exactly the same constant C matrix.

 OPENSEES: zeroLength story springs with HystereticSM, floor masses, UniformExcitation with a Path
 time series, Newmark integrator, Transient analysis.

 Units are free but must be consistent. The record is assumed to be in g (GM_UNIT_G converts to the
 model acceleration unit; use 1.0 if it is already in your unit).
"""
# BOOK: Dynamics of Structures in SI Units -  Anil Kumar Chopra
'https://share.google/TDN5O4eWmw5pH8zUH'
# BOOK: Differential Equations for Engineers-Wei-Chau Xie-CAMBRIDGE-2010

import numpy as np
import matplotlib.pyplot as plt
import PLOT_FORCE_DISP_FUN as S01

#%% ---------------------------------------------------------------
# ONE STORY SPRING: peak-oriented hysteresis (same rule as the SDOF dynamic code, one instance per story)
def MAKE_HYSTERESIS(Displacement, Force):
    """
    MAKE_HYSTERESIS: builds a peak-oriented hysteretic story spring from the multilinear backbone
    and returns the SPRING, INITIAL_STATE and initial stiffness k0 (the unloading stiffness).
    """
    D = np.concatenate(([0.0], Displacement))
    F = np.concatenate(([0.0], Force))
    slopes = np.diff(F) / np.diff(D)
    k0 = slopes[0]                                # initial (unloading) stiffness

    def BACKBONE(x):                              # x >= 0, flat after last point
        """
        BACKBONE (inner): returns the backbone force and tangent stiffness at a positive displacement,
        holding the force constant after the last point.
        """
        if x >= D[-1]:
            return F[-1], 0.0
        seg = np.searchsorted(D, x, side='right') - 1
        return np.interp(x, D, F), slopes[seg]

    def INITIAL_STATE():
        """
        INITIAL_STATE (inner): creates the starting hysteresis state (zero displacement and force,
        with the initial peaks set at the yield point, as in OpenSees Hysteretic).
        """
        # initial peaks = yield point, as in OpenSees Hysteretic
        return dict(u=0.0, f=0.0, d=1, u_rev=0.0, f_rev=0.0, u_a=0.0, f_a=0.0,
                    umax=D[1], fmax=F[1], umin=-D[1], fmin=-F[1])

    def SPRING(u, st):
        """
        SPRING (inner): a pure function that takes the committed state and a trial drift,
        and returns force, tangent and trial state. It covers reversals, unloading on k0,
        partial unloading and reloading toward the previous peak, so it can be re-evaluated
        during Newton iterations without corrupting history.
        """
        # Pure function: committed state + trial drift u -> (unit force, unit tangent, trial state)
        st = dict(st)
        du = u - st['u']
        if du != 0 and np.sign(du) != st['d']:    # reversal
            st['u_a'], st['f_a'] = st['u_rev'], st['f_rev']      # point where the last unloading began
            st['u_rev'], st['f_rev'], st['d'] = st['u'], st['f'], int(np.sign(du))

        d = st['d']                               # mirror so we always move toward +x
        x, xr, gr = d * u, d * st['u_rev'], d * st['f_rev']
        xt, gt = (st['umax'], st['fmax']) if d == 1 else (-st['umin'], -st['fmin'])
        xz = xr - gr / k0                         # where the unloading line hits F = 0

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
        else:                                     # reloading from F = 0 crossing toward previous peak
            k = gt / (xt - xz)
            g = k * (x - xz)

        f = d * g
        st['u'], st['f'] = u, f
        if u > st['umax']: st['umax'], st['fmax'] = u, f
        if u < st['umin']: st['umin'], st['fmin'] = u, f
        return f, k, st

    return SPRING, INITIAL_STATE, k0

#%% ---------------------------------------------------------------
# STRUCTURAL MATRICES AND DYNAMIC PROPERTIES
def ASSEMBLE_K(k_story):
    # Tridiagonal stiffness matrix of a shear building from the story stiffnesses
    """
    ASSEMBLE_K: assembles the tridiagonal shear-building stiffness matrix from the story stiffnesses.
    """
    N = len(k_story)
    K = np.zeros((N, N))
    for i in range(N):
        K[i, i] += k_story[i]
        if i + 1 < N:
            K[i, i] += k_story[i + 1]
            K[i, i + 1] -= k_story[i + 1]
            K[i + 1, i] -= k_story[i + 1]
    return K

def MODAL_PROPERTIES(Displacement, Force, STORY_COFF, MASS):
    # Natural circular frequencies of the INITIAL elastic structure (K0, M) and the mode shapes
    """
    MODAL_PROPERTIES: solves the eigenproblem of the initial elastic structure (K0, M) 
    to get natural circular frequencies, mode shapes and the K0 and M matrices.
    """
    k0 = Force[0] / Displacement[0]
    K0 = ASSEMBLE_K(np.asarray(STORY_COFF, float) * k0)
    M = np.diag(MASS)
    Mih = np.diag(1.0 / np.sqrt(MASS))
    w2, phi = np.linalg.eigh(Mih @ K0 @ Mih)
    return np.sqrt(w2), Mih @ phi, K0, M

def RAYLEIGH_COEFFICIENTS(omega, xi, modes=(1, 2)):
    # C = a0*M + a1*K0 gives damping ratio xi in the two chosen modes (1-based)
    """
    RAYLEIGH_COEFFICIENTS: computes a0 and a1 so that C = a0·M + a1·K0 gives the damping ratio ξ
    in the two chosen modes.
    """
    wi, wj = omega[modes[0] - 1], omega[modes[1] - 1]
    a0 = 2.0 * xi * wi * wj / (wi + wj)
    a1 = 2.0 * xi / (wi + wj)
    return a0, a1

def SCALE_MASS_TO_PERIOD(Displacement, Force, STORY_COFF, MASS_SHAPE, T1_target):
    # Scales the floor masses (keeping their ratios) so that the first period equals T1_target
    """
    SCALE_MASS_TO_PERIOD: scales the floor masses (keeping their ratios) so the first-mode
    period equals the target T1.
    """
    shape = np.asarray(MASS_SHAPE, float)
    omega, _, _, _ = MODAL_PROPERTIES(Displacement, Force, STORY_COFF, shape)
    return shape * (omega[0] * T1_target / (2 * np.pi)) ** 2

#%% ---------------------------------------------------------------
# NUMERICAL: Newmark-beta + Newton-Raphson, N-story shear building
def MDOF_NUMERICAL_DYNAMIC_NEWMARK(Displacement, Force, STORY_COFF, MASS, ag, dt, xi,
                                   MAX_ITERATIONS, MAX_TOLERANCE, MODES=(1, 2)):
    """
    MDOF_NUMERICAL_DYNAMIC_NEWMARK: runs the nonlinear time-history analysis with Newmark-beta and
    a full Newton-Raphson iteration at each step, and returns floor.
    
    STORY_COFF[i] : scale of the shared multilinear backbone for story i (stiffness AND strength)
    MASS[i]       : lumped mass of floor i
    ag            : ground acceleration history (already in model units), sampled every dt
    """
    # Newmark-beta method parameters
    gamma, beta= 0.5, 0.25
    
    scale = np.asarray(STORY_COFF, float)
    MASS = np.asarray(MASS, float)
    N = len(scale)
    SPRING, INITIAL_STATE, k0 = MAKE_HYSTERESIS(Displacement, Force)
    F_ref = scale.max() * np.max(Force)                 # normalisation of the residual norm
    k_floor = 1.0e-8 * k0 * scale.min()                 # keeps K_t non-singular on flat branches

    omega, phi, K0, M = MODAL_PROPERTIES(Displacement, Force, STORY_COFF, MASS)
    a0, a1 = RAYLEIGH_COEFFICIENTS(omega, xi, MODES)
    C = a0 * M + a1 * K0                                # constant Rayleigh damping matrix
    r = np.ones(N)                                      # influence vector

    def ASSEMBLE(k_story):
        """
        ASSEMBLE (inner): assembles the tangent stiffness matrix, replacing near-zero story tangents
        with a tiny floor value so the matrix stays non-singular on flat branches.
        """
        k_story = np.where(np.abs(k_story) < k_floor, k_floor, k_story)
        
        return ASSEMBLE_K(k_story)

    def EVALUATE(u_vec, states):                        # story drifts -> shears, tangents, trial states
        """
        EVALUATE (inner): converts trial floor displacements into story drifts, then into story shears,
        tangents, floor internal forces and trial hysteresis states.
        """
        drift = np.diff(np.concatenate(([0.0], u_vec)))
        V, k, trial = np.zeros(N), np.zeros(N), []
        for i in range(N):
            f, ki, ti = SPRING(drift[i], states[i])
            V[i], k[i] = scale[i] * f, scale[i] * ki
            trial.append(ti)
        f_int = V - np.concatenate((V[1:], [0.0]))      # internal force at each floor
        
        return V, k, f_int, trial

    # committed state: system at rest (OpenSees Transient also starts with a(0) = 0)
    states = [INITIAL_STATE() for _ in range(N)]
    u, v, a = np.zeros(N), np.zeros(N), np.zeros(N)

    n = len(ag)
    U, VEL, ACC = np.zeros((n, N)), np.zeros((n, N)), np.zeros((n, N))
    VS = np.zeros((n, N))                               # story shears
    ITER, R_final = np.zeros(n, dtype=int), np.zeros(n)

    for i in range(1, n):
        p_new = -MASS * r * ag[i]                       # effective ground-motion load
        u_new = u.copy()                                # initial guess: previous displacement

        it = 0
        while True:
            # Newmark kinematics as functions of the trial displacement vector
            a_new = (u_new - u) / (beta * dt**2) - v / (beta * dt) - (1.0 / (2 * beta) - 1.0) * a
            v_new = v + dt * ((1 - gamma) * a + gamma * a_new)

            V, k, f_int, trial = EVALUATE(u_new, states)
            R = M @ a_new + C @ v_new + f_int - p_new   # DYNAMIC RESIDUAL (vector)
            DX = np.linalg.norm(R) / F_ref
            if DX <= MAX_TOLERANCE or it >= MAX_ITERATIONS:
                break
            K_eff = ASSEMBLE(k) + gamma / (beta * dt) * C + M / (beta * dt**2)
            u_new = u_new - np.linalg.solve(K_eff, R)
            it += 1

        R_final[i], ITER[i] = DX, it

        if DX > MAX_TOLERANCE:
            print(f'(-) Step {i}: NOT converged, |R|/Fref = {DX:.3e}')
            n = i
            break

        # commit
        states, u, v, a = trial, u_new, v_new, a_new
        U[i], VEL[i], ACC[i], VS[i] = u, v, a, V

    U, VEL, ACC, VS = U[:n], VEL[:n], ACC[:n], VS[:n]
    drift = np.diff(np.hstack((np.zeros((n, 1)), U)), axis=1)
    
    DATA = dict(t=np.arange(n) * dt, U_floor=U, V_floor=VEL, A_floor=ACC, drift=drift,
                V_story=VS, base_reaction=-VS[:, 0],          # same sign as OpenSees base reaction
                ITER=ITER[:n], R_final=R_final[:n], omega=omega, phi=phi, a0=a0, a1=a1)
    return DATA

#%% ---------------------------------------------------------------
# OPENSEES: transient analysis with uniform base excitation
def MDOF_OPENSEES_DYNAMIC(MAT_TYPE, Displacement, Force, STORY_COFF, MASS, ag_raw, dt, xi,
                          GM_FACTOR,  MAX_ITERATIONS, MAX_TOLERANCE, MODES=(1, 2)):
    # ag_raw is the raw record; GM_FACTOR scales it to the model's acceleration unit
    """
    MDOF_OPENSEES_DYNAMIC: builds the same shear building in OpenSeesPy (zeroLength HystereticSM springs,
    lumped masses, Rayleigh damping, UniformExcitation), runs the transient analysis, and returns results
    for comparison.
    """
    import openseespy.opensees as ops
    import ANALYSIS_FUNCTION as S02
    
    # Newmark-beta method parameters
    gamma, beta = 0.5, 0.25
    
    N = len(STORY_COFF)
    
    # Reset model
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)

    # Nodes: 0 = ground, 1..N = floors
    for n_ in range(N + 1):
        ops.node(n_, 0.0)
        
    # Constrains
    ops.fix(0, 1)
    
    # Masses
    for i in range(1, N + 1):
        ops.mass(i, float(MASS[i - 1]))

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
            # INFO LINK: https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/HystereticSM.html
            
        # '-doRayleigh', 1 is REQUIRED: zeroLength ignores the stiffness-proportional Rayleigh term by default
        ops.element('zeroLength', i, i - 1, i, '-mat', i, '-dir', 1, '-doRayleigh', 1)
        # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/ZeroLength.html

    # Rayleigh damping on the INITIAL stiffness: C = a0*M + a1*K0  (same constant C as the numerical model)
    omega, _, _, _ = MODAL_PROPERTIES(Displacement, Force, STORY_COFF, np.asarray(MASS, float))
    a0, a1 = RAYLEIGH_COEFFICIENTS(omega, xi, MODES)
    ops.rayleigh(a0, 0.0, a1, 0.0)
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/reyleigh.html

    # ground motion
    ops.timeSeries('Path', 1, '-dt', dt, '-values', *ag_raw, '-factor', GM_FACTOR)
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
        
    ops.integrator('Newmark', gamma, beta) # INFO LINK: https://opensees.github.io/OpenSeesDocumentation/user/manual/analysis/integrator/Newmark.html
    ops.analysis('Transient')              # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/analysis.html

    t_os, U_os, F_os, V_os = [0.0], [np.zeros(N)], [0.0], [np.zeros(N)]
    for i in range(1, len(ag_raw)):
        OK = ops.analyze(1, dt)
        if OK != 0:
            S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS)
        ops.reactions()
        t_os.append(ops.getTime())
        U_os.append([ops.nodeDisp(n_, 1) for n_ in range(1, N + 1)])   # FLOOR DISPLACEMENTS
        F_os.append(ops.nodeReaction(0, 1))                            # BASE REACTION
        V_os.append([-ops.eleResponse(j, 'force')[0] for j in range(1, N + 1)])  # STORY SHEARS

    U_floor = np.array(U_os)
    drift = np.diff(np.hstack((np.zeros((len(U_floor), 1)), U_floor)), axis=1)
    
    DATA = dict(t=np.array(t_os), U_floor=U_floor, drift=drift, base_reaction=np.array(F_os), V_story=np.array(V_os))
    
    return DATA

#%% ---------------------------------------------------------------
# MAIN
if __name__ == '__main__':
    """
    defines the backbone, building and ground motion, runs both solvers, prints the error comparison
    (displacement, base shear, story shear), and plots the time histories, hysteresis loops and drift profile.
    """
    # ---- BACKBONE (per unit story, in story-drift terms; units: N, m, kg, s) ----
    FY, DY = 2000.0, 0.01                 # Yield force & drift
    FU, DU = 1.2 * FY, DY * 10            # Ultimate (peak) force & drift
    FR1, DR1 = FY / 4,  DY * 12           # 1st residual / post-peak point
    FR2, DR2 = FY / 10, DY * 15           # 2nd residual / post-peak point

    Force = np.array([FY, FU, FR1, FR2])
    Displacement = np.array([DY, DU, DR1, DR2])
    
    # BUILDING
    N_STORY = 10
    # Story scale factors (stiffness AND strength). Several parallel springs per story can be
    # lumped here as sum(Coff_story).
    STORY_COFF = np.ones(N_STORY)
    MASS_SHAPE = np.ones(N_STORY)         # floor-mass ratios (scaled below to reach T1)
    T1 = 1.0                              # target first-mode period (s)
    XI = 0.05                             # damping ratio in the two Rayleigh modes
    MODES = (1, 2)

    MASS = SCALE_MASS_TO_PERIOD(Displacement, Force, STORY_COFF, MASS_SHAPE, T1)
    
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

    omega, phi, _, _ = MODAL_PROPERTIES(Displacement, Force, STORY_COFF, MASS)
    print('Periods of the first 3 modes [s]:', np.round(2 * np.pi / omega[:3], 4))
    print(f'Floor mass = {MASS[0]:.3f},  N = {len(ag)},  dt = {DT}')
    
    #%% DEFINE ANALYSIS PROPERTIES
    MAX_ITERATIONS = 50    # Convergence iteration for test
    MAX_TOLERANCE = 1.0e-6   # Convergence tolerance for test
    
    #%% RUN NUMERICAL METHOD
    NUM = MDOF_NUMERICAL_DYNAMIC_NEWMARK(Displacement, Force, STORY_COFF, MASS, ag, DT, XI,
                                         MAX_ITERATIONS, MAX_TOLERANCE, MODES=MODES)
    print(f'Rayleigh: a0 = {NUM["a0"]:.5f},  a1 = {NUM["a1"]:.5f}')

    #%% RUN OPENSEES
    OS = MDOF_OPENSEES_DYNAMIC('INELASTIC', Displacement, Force, STORY_COFF, MASS, ag_raw, DT, XI,
                               GM_FACTOR, MAX_ITERATIONS, MAX_TOLERANCE, MODES=MODES)

    n = min(len(NUM['t']), len(OS['t']))
    errU = np.abs(OS['U_floor'][:n] - NUM['U_floor'][:n])
    errF = np.abs(OS['base_reaction'][:n] - NUM['base_reaction'][:n])
    errV = np.abs(OS['V_story'][:n] - NUM['V_story'][:n])
    print(f'max |roof U| = {np.abs(NUM["U_floor"][:, -1]).max():.4e}   '
          f'max |base shear| = {np.abs(NUM["base_reaction"]).max():.4e}')
    print(f'Floor disp  : max abs err = {errU.max():.3e}, max rel err = {errU.max() / np.abs(NUM["U_floor"]).max():.3e}')
    print(f'Base shear  : max abs err = {errF.max():.3e}, max rel err = {errF.max() / np.abs(NUM["base_reaction"]).max():.3e}')
    print(f'Story shear : max abs err = {errV.max():.3e}')
    print(f'Max Newton iterations per step = {NUM["ITER"].max()}')
    print('Peak drift / yield drift per story:', np.round(np.abs(NUM['drift']).max(axis=0) / DY, 2))

    # ---- Plots ----
    t1 = NUM['t']
    roof = NUM['U_floor'][:, -1]
    roof_os = OS['U_floor'][:, -1]
    fig, ax = plt.subplots(2, 3, figsize=(19, 10))

    ax[0, 0].plot(np.arange(len(ag)) * DT, ag, 'g')
    ax[0, 0].set_title('Ground acceleration'); ax[0, 0].set_xlabel('Time [s]'); ax[0, 0].set_ylabel('ag [m/s^2]')

    ax[0, 1].plot(t1, roof, 'k', lw=3, label='Numerical')
    ax[0, 1].plot(OS['t'], roof_os, 'r--', lw=2, label='OpenSeesPy')
    ax[0, 1].set_title('Roof displacement time history'); ax[0, 1].set_xlabel('Time [s]')
    ax[0, 1].set_ylabel('Displacement [m]'); ax[0, 1].legend()

    ax[0, 2].plot(roof, NUM['base_reaction'], 'k', lw=3, label='Numerical')
    ax[0, 2].plot(roof_os, OS['base_reaction'], 'r--', lw=2, label='OpenSeesPy')
    ax[0, 2].set_title('Global hysteresis loops'); ax[0, 2].set_xlabel('Roof displacement [m]')
    ax[0, 2].set_ylabel('Base reaction [N]'); ax[0, 2].legend()

    colors = plt.cm.viridis(np.linspace(0, 0.9, N_STORY))
    for i in range(N_STORY):
        ax[1, 0].plot(NUM['drift'][:, i], NUM['V_story'][:, i], color=colors[i], lw=2, label=f'Story {i + 1}')
        ax[1, 0].plot(OS['drift'][:n, i], OS['V_story'][:n, i], 'r--', lw=0.8)
    ax[1, 0].set_title('Story hysteresis (solid: numerical, red dashed: OpenSeesPy)')
    ax[1, 0].set_xlabel('Story drift [m]'); ax[1, 0].set_ylabel('Story shear [N]')
    ax[1, 0].legend(fontsize=7, ncol=2)

    ax[1, 1].plot(np.abs(NUM['drift']).max(axis=0) / DY, np.arange(1, N_STORY + 1), 'ko-', lw=2, ms=8, label='Numerical')
    ax[1, 1].plot(np.abs(OS['drift']).max(axis=0) / DY, np.arange(1, N_STORY + 1), 'r--s', lw=2, label='OpenSeesPy')
    ax[1, 1].set_title('Peak story drift profile'); ax[1, 1].set_xlabel('Peak drift / yield drift')
    ax[1, 1].set_ylabel('Story'); ax[1, 1].legend()

    ax[1, 2].plot(t1, NUM['A_floor'][:, -1] + ag[:len(t1)], 'b')
    ax[1, 2].set_title('Roof absolute acceleration (numerical) [m/s^2]'); ax[1, 2].set_xlabel('Time [s]')

    for a_ in ax.ravel(): a_.grid(True)
    plt.tight_layout()
    plt.savefig('MDOF_DYNAMIC_TIME_HISTORY.png', dpi=150)
    plt.show()
