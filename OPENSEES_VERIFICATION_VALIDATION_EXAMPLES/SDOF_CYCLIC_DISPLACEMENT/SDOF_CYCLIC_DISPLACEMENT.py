###################################################################################################
#                >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                 #
#   CYCLIC DISPLACEMENT-CONTROLLED NEWTON-RAPHSON PUSHOVER OF PARALLEL NONLINEAR SDOF SPRINGS     #
#-------------------------------------------------------------------------------------------------#
#                  FINITE ELEMENT METHOD IN PYTHON AND OPENSEESPY VERIFICATIONS                   #
#-------------------------------------------------------------------------------------------------#
#                   THIS PROGRAM IS WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                 # 
#                               E-mail:salar.d.ghashghaei@gmail.com                               #
###################################################################################################
"""
 This python script drives the same parallel-spring SDOF system through a back-and-forth displacement
 history (0 → +a → −a → 0 at growing amplitudes) using two independent approaches and compares them.
 DISP_HISTORY builds the list of target displacements step by step. The numerical solver
 SDOF_NUMERICAL_CYCLIC_NEWTON works like the monotonic version (predict with the last tangent,
 compute the residual, correct the force until balanced), but because force now depends on loading history,
 the SPRING function remembers the last reversal point and the largest peaks on each side.
 It then decides the branch: follow the backbone when pushing past the previous peak, unload along 
 the initial stiffness after a reversal, or reload in a straight line toward the previous peak on the
 opposite side. SDOF_OPENSEES_CYCLIC builds the equivalent model with HystereticSM springs and re-issues
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
# NUMERICAL: displacement-controlled Newton, peak-oriented hysteresis
"""
This function runs the parallel-spring SDOF system through a prescribed displacement history (U_target)
 and records the force at every step. It first builds the backbone curve with the origin added, takes
 the first slope k0 as the unloading stiffness, and lumps all parallel springs into one spring scaled by Csum.
 At each step it predicts the force from the last committed tangent, calls SPRING to get the true force at the
 new displacement, computes the residual (predicted minus true), and corrects the force by subtracting the residual
 until it falls below TOL. The SPRING function holds the hysteresis logic: it detects a load reversal, stores the
 reversal point, and mirrors the coordinates so motion always looks like it goes toward +x. It then chooses one of
 three branches: the backbone when the displacement goes beyond the previous peak, unloading along k0 right after
 a reversal, or reloading along a straight line from the zero-force crossing toward the previous peak.
 It also updates the largest positive and negative excursions. If a step converges, the new state, displacement,
 and tangent are committed, and the displacement and force are stored. If not, the analysis stops. The only change
 from the earlier version is that the stored force is now -lam, which flips the sign so the numerical reaction matches
 the sign convention of the OpenSees base reaction (the reaction at the fixed support is opposite to the applied load).

  One thing to check: Fb.append(-lam) flips only the stored output, while the solver itself still works with lam as
 the internal force. That is fine for plotting and comparison, but if you reuse Fb for anything inside the solver
 (such as R_init or energy calculations), remember it has the opposite sign to lam.
"""
def SDOF_NUMERICAL_CYCLIC_NEWTON(Displacement, Force, Coff, U_target, MAX_ITERATIONS, TOL):
    D = np.concatenate(([0.0], Displacement))
    F = np.concatenate(([0.0], Force))
    slopes = np.diff(F) / np.diff(D)
    k0 = slopes[0]                                # unloading stiffness (beta = 0)
    Csum = np.sum(Coff)                           # parallel springs -> one scaled spring
    F_ref = Csum * F.max()                        # reference force used only for normalizing residuals during the Newton–Raphson iterations

    def BACKBONE(x):                              # x >= 0, flat after last point
        """
        BACKBONE(x):
        For a non-negative displacement $x$, returns the force and tangent stiffness
        from the multilinear backbone curve (interpolating between points or returning
        the last force with zero stiffness if x exceeds the ultimate displacement).
        """
        if x >= D[-1]:
            return F[-1], 0.0
        seg = np.searchsorted(D, x, side='right') - 1 # finds the index of the segment of the multilinear backbone that contains the displacement x
        return np.interp(x, D, F), slopes[seg]        # this pair is later used by the Newton–Raphson solver as the current force and the current tangent stiffness of the spring when it is on the backbone branch

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
        return Csum * f, Csum * k, st

    # committed state (initial peaks = yield point, as in OpenSees Hysteretic)
    state = dict(u=0.0, f=0.0, d=1, u_rev=0.0, f_rev=0.0,
                 umax=D[1], fmax=F[1], umin=-D[1], fmin=-F[1])
    
    u, lam, Kt = 0.0, 0.0, Csum * k0
    U, Fb = [0.0], [0.0]
    R_init, R_final, ITER, CONVERGED = [], [], [], []

    for i in range(1, len(U_target)):
        u_new = U_target[i]

        lam_new = lam + Kt * (u_new - u)                    # 1) predictor
        f, k, trial = SPRING(u_new, state)
        R = lam_new - f                                     # 2) residual
        R_init.append(abs(R) / F_ref)

        it = 0
        # 3) corrector
        while abs(R) / F_ref > TOL and it < MAX_ITERATIONS:
            # This is the iterative correction step that drives the residual to zero
            lam_new -= R                        # Correct the predicted force by subtracting the residual
            f, k, trial = SPRING(u_new, state)  # Re-evaluate the spring at the fixed target displacement u_new with the (still uncommitted) state, obtaining the true force f, tangent k, and trial state.
            R = lam_new - f                     # Recompute the residual
            it += 1                             # Increment the iteration counter

        ok = abs(R) / F_ref <= TOL              # Sets a Boolean flag that is True only when the final normalized residual is smaller than (or equal to) the prescribed tolerance, i.e. the step has converged.
        R_final.append(abs(R) / F_ref)          # Stores the magnitude of that final residual (normalized by the reference force) so you can later inspect how well each displacement step satisfied equilibrium.
        ITER.append(it)
        CONVERGED.append(ok)
        
        if not ok:
            print(f'(-) Step {i}: NOT converged, |R|/Fref = {abs(R)/F_ref:.3e}')
            break
        
        # 4) commit
        # After a successful convergence, this single line permanently updates the solution state for the next displacement increment:
        state, u, lam, Kt = trial, u_new, lam_new, k   
        
        U.append(u)                                      # DISPLACEMENT
        Fb.append(-lam)                                  # REACTION
        
        DATA = (np.array(U), np.array(Fb), np.array(R_init),
                np.array(R_final), np.array(ITER), np.array(CONVERGED))
        
    return DATA

#%% ---------------------------------------------------------------
# OPENSEES: Varying displacement increment each step
def SDOF_OPENSEES_CYCLIC(MAT_TYPE, Displacement, Force, Coff, U_target,
                         MAX_ITERATIONS=500, MAX_TOLERANCE=1.0e-6):
    import openseespy.opensees as ops
    import ANALYSIS_FUNCTION as S02
    
    # Reset model
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)
    
    # Nodes
    ops.node(1, 0.0)
    ops.node(2, 0.0)
    
    ops.fix(1, 1)

    # Each spring j = parallel zeroLength element with MultiLinear backbone
    for j in range(len(Coff)):
        tag = 1 + j
        if MAT_TYPE == 'ELASTIC':
            ops.uniaxialMaterial('Elastic', tag, Coff[j] * Force[0] / Displacement[0])
        if MAT_TYPE == 'INELASTIC':
            pos, neg = [], []
            for f, d in zip(Force, Displacement):
                pos.extend([Coff[j] * f, d]); neg.extend([-Coff[j] * f, -d])
            ops.uniaxialMaterial('HystereticSM', tag, *pos, *neg, 1.0, 1.0, 0.0, 0.0, 0.0)
        ops.element('zeroLength', tag, 1, 2, '-mat', tag, '-dir', 1)

    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(2, 1.0)

    
    ops.system('BandGeneral')
    ops.numberer('Plain')
    ops.constraints('Plain')
    ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/normDispIncr.html

    if MAT_TYPE == 'INELASTIC':
        ops.algorithm('Newton')  # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/algorithm.html
    if MAT_TYPE == 'ELASTIC':
        ops.algorithm('Linear') 
        
    ops.analysis('Static')

    U_os, F_os = [0.0], [0.0]
    for i in range(1, len(U_target)):
        du = U_target[i] - ops.nodeDisp(2, 1)                  # increment to reach the target
        ops.integrator('DisplacementControl', 2, 1, du)        # re-issue every step
        OK = ops.analyze(1)
        S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS)
        ops.reactions()
        U_os.append(ops.nodeDisp(2, 1))                        # DISPLACEMENT
        F_os.append(ops.nodeReaction(1, 1))                    # REACTION
        
    DATA = (np.array(U_os), np.array(F_os))  
    
    return DATA

#%% ---------------------------------------------------------------
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
    
    S01.PLOT_FORCE_DISP_FUN(Displacement, Force)

    U_target = DISP_HISTORY(DY * np.array([1, 2, 4, 6, 8, 10, 12, 15]), n_leg=50)
    # RUN NUMERICAL METHOD
    D1, F1, R0, R1, IT, OK = SDOF_NUMERICAL_CYCLIC_NEWTON(Displacement, Force, Coff,
                                                          U_target, 50, 1.0e-10)
    # RUN OPENSEES
    U_os, F_os = SDOF_OPENSEES_CYCLIC('INELASTIC', Displacement, Force, Coff, U_target)

    n = min(len(F1), len(F_os))
    err = np.abs(F_os[:n] - (-F1[:n]))
    print(f'max abs err = {err.max():.3e}, max rel err = {err.max() / np.abs(F1).max():.3e}')

    fig, ax = plt.subplots(1, 2, figsize=(14, 6))
    ax[0].plot(D1, F1, 'k', lw=4, label='Numerical')
    ax[0].plot(U_os, F_os, 'r--', lw=3, label='OpenSeesPy')
    ax[0].set_xlabel('Displacement')
    ax[0].set_ylabel('Base Reaction')
    ax[0].set_title('Hysteresis loops')
    ax[0].grid(True)
    ax[0].legend()
    ax[1].plot(U_target, 'b'); ax[1].set_xlabel('Step')
    ax[1].set_ylabel('Displacement')
    ax[1].set_title('Cyclic displacement history')
    ax[1].grid(True)
    plt.tight_layout()
    plt.show()
