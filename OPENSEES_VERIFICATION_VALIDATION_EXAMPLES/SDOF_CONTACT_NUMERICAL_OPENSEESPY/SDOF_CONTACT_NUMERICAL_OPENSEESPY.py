######################################################################################################################
#                    >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                                #
#       NUMERICAL NEWTON-RAPHSON VERIFICATION OF THE OPENSEES CONTACT-DRIVEN SDOF MONOTINIC PUSHOVER MODEL           #
#--------------------------------------------------------------------------------------------------------------------#
#                               THIS PROGRAM IS WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                        # 
#                                               E-mail:salar.d.ghashghaei@gmail.com                                  #
######################################################################################################################
"""
 This Python script implements a Numerical Newton-Raphson solver to independently verify
 an OpenSees contact-driven single-degree-of-freedom (SDOF) monotonic pushover analysis.
 It defines a symmetric piecewise-linear force-displacement backbone, assembles internal
 forces and tangent stiffness for two parallel springs (one always active, the second
 activated after a gap closes), and performs displacement-controlled pushover analysis under
 two interpretations of the contact spring: an “OPENSEES” mode that matches the original model’s
 strain definition and a more physical “gap” mode. Force-controlled Newton-Raphson is used as an
 inverse check, periods are computed from the tangent stiffness, and results (pushover curves and
 period evolution) are compared with an optional OpenSeesPy run that dynamically adds the second
 zeroLength element once the gap is closed.
"""

import numpy as np
import matplotlib.pyplot as plt

#%%------------------------------------------------------------------
# Structural properties (N, m)
FY = 85000.0                    # [N]   Yield force
FU = 1.2 * FY                   # [N]   Ultimate force
KE = 4500000.0                  # [N/m] Elastic stiffness
DY = FY / KE                    # [m]   Yield displacement
DSU = 0.36                      # [m]   Ultimate displacement
M = 150000.0                    # [kg]  Mass

DMAX = 0.35                     # [m] Max. Pushover Incremental Displacement
DINCR = 0.001                   # [m] Pushover Increment
GAP_DIST = 0.25                 # [m] Gap Distance (Conact Distance)

# Positive branch points
pos_disp = [0, DY, DSU, 1.1*DSU, 1.25*DSU]
pos_force = [0, FY, FU, 0.2*FU, 0.1*FU]
KP = np.array([FY, DY, FU, DSU, 0.2*FU, 1.1*DSU, 0.1*FU, 1.25*DSU])

# Negative branch points
neg_disp = [0, -DY, -DSU, -1.1*DSU, -1.25*DSU]
neg_force = [0, -FY, -FU, -0.2*FU, -0.1*FU]
KN = np.array([-FY, -DY, -FU, -DSU, -0.2*FU, -1.1*DSU, -0.1*FU, -1.25*DSU])

# Plot
plt.plot(pos_disp, pos_force, marker='o', color='red')
plt.plot(neg_disp, neg_force, marker='o', color='black')

plt.xlabel("Displacement [m]")
plt.ylabel("Force [N]")
plt.title("Force–Displacement Curve")
plt.grid(True)
plt.axhline(0, linewidth=0.5)
plt.axvline(0, linewidth=0.5)
plt.show()

#%% DEFINE ANALYSIS PROPERTIES
MAX_ITERATIONS = 500     # Convergence iteration for test
MAX_TOLERANCE = 1.0e-6   # Convergence tolerance for test

#%%------------------------------------------------------------------
# 1) MATERIAL: symmetric piecewise-linear backbone (monotonic loading only)

def backbone(u):
    # Return (force, tangent stiffness) of the backbone at strain u
    """
    The backbone() function evaluates the piecewise-linear force-displacement envelope of
    the hysteretic material at a given strain u. It first determines the sign s of the strain
    and works with its absolute value a; if a exceeds the last defined displacement point it
    returns a constant residual force with zero tangent stiffness. Otherwise it locates the
    linear segment that contains a by binary search on the positive displacement breakpoints,
    computes the slope k of that segment, interpolates the corresponding force magnitude, and finally
    returns the signed force together with the (unsigned) tangent stiffness.
    """
    s = 1.0 if u >= 0.0 else -1.0
    a = abs(u)
    if a >= pos_disp[-1]:                                   # beyond last envelope point -> flat residual
        return s * pos_force[-1], 0.0
    i = int(np.searchsorted(pos_disp, a, side='right') - 1)  # segment to the right of the point
    k = (pos_force[i + 1] - pos_force[i]) / (pos_disp[i + 1] - pos_disp[i])
    f = pos_force[i] + k * (a - pos_disp[i])
    
    return s * f, k


#%%------------------------------------------------------------------
# 2) STRUCTURE: assemble internal force and tangent stiffness at node 2

def assemble(u, spring2_active, mode):
    """
    The assemble() function evaluates the internal force and tangent stiffness of 
    the SDOF system at a given displacement u. It always queries the piecewise-linear
    backbone for spring 1 (which is permanently active) to obtain its force f1 and stiffness k1.
    If the second spring is flagged as active, it computes that spring’s strain according to the
    chosen contact interpretation—simply u in 'OPENSEES' mode or u - GAP_DIST in the physically
    correct gap mode—then obtains the corresponding force f2 and stiffness k2 from the same backbone.
    Finally it returns the individual spring forces together with the assembled tangent stiffness (k1 + k2)
    and total restoring force (f1 + f2) that are required by the subsequent Newton-Raphson iterations.
    """
    f1, k1 = backbone(u)                                 # Spring 1 (always active)
    f2, k2 = 0.0, 0.0
    if spring2_active:
        strain2 = u if mode == 'OPENSEES' else u - GAP_DIST
        f2, k2 = backbone(strain2)                       # Spring 2 (parallel, node 3 fixed)
    
    return f1, f2, k1 + k2, f1 + f2


#%%------------------------------------------------------------------
# 3) DISPLACEMENT-CONTROL NEWTON-RAPHSON (OpenSees DisplacementControl logic, reference load P = 1)

def newton_disp_control(u_n, lam_n, du_target, spring2_active, mode, P=1.0):
    """
    The newton_disp_control() function implements a displacement-controlled Newton-Raphson
    algorithm that advances the system by a prescribed displacement increment `du_target` while
    solving for the corresponding load factor `lam`. Starting from the previous converged state `(u_n, lam_n)`,
    it first computes a predictor step using the current tangent stiffness, then iteratively corrects the residual
    unbalanced force `R = lam*P - F_int(u)` by updating both displacement and load factor so that the controlled
    degree-of-freedom always advances exactly by `du_target`; the loop terminates when the residual falls below
    a relative tolerance or the maximum iteration count is reached, returning the updated state, the number of
    iterations, and a convergence flag.
    """
    u, lam = u_n, lam_n
    Kt = assemble(u, spring2_active, mode)[2]
    dUp = P / Kt
    dlam = du_target / dUp                                # predictor: control DOF moves by du_target
    lam += dlam
    u += dlam * dUp
    for it in range(1, MAX_ITERATIONS + 1):
        fint = assemble(u, spring2_active, mode)[3]
        R = lam * P - fint                                # unbalanced load
        if abs(R) <= MAX_TOLERANCE * max(1.0, abs(lam * P)):
            return u, lam, it - 1, True
        Kt = assemble(u, spring2_active, mode)[2]
        dUp, dUr = P / Kt, R / Kt
        dlam = -(dUr - 0.0) / dUp                          # keep control-DOF increment = du_target
        u += dUr + dlam * dUp
        lam += dlam
        
    return u, lam, MAX_ITERATIONS, False


#%%------------------------------------------------------------------
# 4) FORCE-CONTROL NEWTON-RAPHSON (inverse check):  find u such that  F_int(u) = F_target

def newton_force_control(F_target, u0, spring2_active, mode):
    """
    The newton_force_control function performs a force-controlled Newton-Raphson iteration
    to solve the inverse problem of finding the displacement `u` that produces a prescribed
    internal force `F_target`. Starting from an initial guess `u0`, it repeatedly assembles
    the current tangent stiffness `Kt` and internal force `fint` via the `assemble` routine
    (which accounts for whether the second spring is active and which contact mode is used),
    computes the residual force imbalance `R = F_target - fint`, and updates the displacement
    by the classical Newton step `u += R / Kt` until the residual falls below a relative tolerance
    or the maximum number of iterations is reached, finally returning the converged displacement
    together with the number of iterations performed.
    """
    u = u0
    for it in range(1, MAX_ITERATIONS + 1):
        _, _, Kt, fint = assemble(u, spring2_active, mode)
        R = F_target - fint
        if abs(R) <= MAX_TOLERANCE * max(1.0, abs(F_target)):
            return u, it - 1
        u += R / Kt
    
    return u, MAX_ITERATIONS


#%%------------------------------------------------------------------
# 5) NUMERICAL MONOTONIC PUSHOVER

def NUMERICAL_MONOTONIC_PUSHOVER(mode):
    """
    The function `NUMERICAL_MONOTONIC_PUSHOVER(mode)` drives a displacement-controlled monotonic pushover
    analysis by advancing the system in fixed increments of `DINCR` up to `DMAX`. At each step it calls
    the displacement-control Newton–Raphson solver (`newton_disp_control`) to obtain a converged displacement
    `u` and load factor `lam`, assembles the spring forces and tangent stiffness, and then checks whether the
    gap has closed; if so it activates the second spring (matching the order used in the OpenSees script) and
    updates the reported tangent. An independent force-controlled Newton–Raphson solve is performed as an inverse 
    verification, after which the function records displacement, individual and total forces, secant and tangent
    stiffnesses, period, iteration count and inverse-error into a dictionary that is finally returned as NumPy arrays.
    """
    
    n_steps = int(round(abs(DMAX / DINCR)))
    
    u, lam, active = 0.0, 0.0, False
    out = {k: [] for k in ['u', 'f1', 'f2', 'V', 'K', 'Kt', 'T', 'iters', 'inv_err']}
    for step in range(n_steps):
        u, lam, nit, ok = newton_disp_control(u, lam, DINCR, active, mode)
        if not ok:
            raise RuntimeError(f'Newton did not converge at step {step + 1}')
        f1, f2, Kt, V = assemble(u, active, mode)

        # Activation check performed AFTER the converged step (same order as the OpenSees script)
        Kt_report = Kt
        if (u >= GAP_DIST) and not active:
            active = True
            f2, V = 0.0, f1                           # element just added: no force recorded for it this step
            Kt_report = assemble(u, True, mode)[2]    # eigen analysis in the script sees the new tangent

        # Inverse check with force-controlled NR, started from the previous converged displacement
        Vcheck = assemble(u, active, mode)[3]
        u_inv, _ = newton_force_control(Vcheck, u - DINCR, active, mode)

        out['u'].append(u);  out['f1'].append(f1);  out['f2'].append(f2)
        out['V'].append(abs(V)); out['K'].append(abs(V) / abs(u))
        out['Kt'].append(Kt_report)
        out['T'].append(2.0 * np.pi * np.sqrt(M / Kt_report))
        out['iters'].append(nit)
        out['inv_err'].append(abs(u_inv - u))
    
    return {k: np.array(v) for k, v in out.items()}


#%%------------------------------------------------------------------
# 6) OPENSEESPY
def OPENSEES_PUSHOVER_CONTACT():
    import openseespy.opensees as ops
    import ANALYSIS_FUNCTION as S02

    Ke = FY / DY                          # [N/m] Spring Elastic Stiffness
    KP = np.array([FY, DY,
                   FU, DSU,
                   0.2 * FU, 1.1 * DSU,
                   0.1 * FU, 1.25 * DSU]) # [N/m] Spring Plastic Stiffness
    KN = -KP
        
    #%% Wipe and create model
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)
    
    #%% Create nodes (both at x=0 for zeroLength)
    ops.node(1, 0.0)
    ops.node(2, 0.0)
    ops.node(3, 0.0)
    
    #%% Define Constrain
    ops.fix(1, 1) 
    ops.fix(3, 1) 
    
    #%% Define Mass
    ops.mass(2, M)
    
    #%% Define uniaxial materials for each spring
    C = 2 * 0.05 * np.sqrt(Ke / M) * M
    ops.uniaxialMaterial('HystereticSM', 1,
                         '-posEnv', *KP,
                         '-negEnv', *KN,
                         '-pinch', 1, 1)
    
    ops.uniaxialMaterial('Viscous', 2, C, 1.0)
    ops.element('zeroLength', 1, 1, 2, '-mat', 1, 2, '-dir', 1, 1)
    #%% Create a dummy load pattern for displacement control
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(2, 1.0)
    
    ops.integrator('DisplacementControl', 2, 1, DINCR) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/integrator.html
    ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/test.html
    ops.algorithm('Newton')   # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/algorithm.html
    ops.system('BandGeneral') # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/system.html
    ops.numberer('Plain')     # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/numberer.html
    ops.constraints('Plain')  # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/constraints.html
    ops.analysis('Static')    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/analysis.html
    
    u_l, V_l, T_l, added = [], [], [], False
    for _ in range(int(round(DMAX / DINCR))):
        OK = ops.analyze(1)
        if OK != 0:
            S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS)
            
        u = ops.nodeDisp(2, 1)
        f1 = ops.eleResponse(1, 'force')[1]
        if u >= GAP_DIST and not added:
            ops.element('zeroLength', 2, 3, 2, '-mat', 1, 2, '-dir', 1, 1)
            ops.domainChange(); added = True
        f2 = ops.eleResponse(2, 'force')[1] if added else 0.0
        lam = ops.eigen('-fullGenLapack', 1)[0]
        u_l.append(u); V_l.append(abs(f1 + f2)); T_l.append(2 * np.pi / np.sqrt(lam))
    
    DATA = (np.array(u_l), np.array(V_l), np.array(T_l))
    return DATA


#%%------------------------------------------------------------------
if __name__ == '__main__':
    #%% RUN NUMERICAL INELASTIC ANALYSIS
    R_os = NUMERICAL_MONOTONIC_PUSHOVER('OPENSEES')
    R_ph = NUMERICAL_MONOTONIC_PUSHOVER('NUMERICAL')

    print('\n================  NEWTON-RAPHSON SUMMARY  ================')
    for name, R in [('OPENSEES-equivalent', R_os), ('NUMERICAL', R_ph)]:
        print(f'\n[{name}]')
        print(f'  Max NR iterations / step        : {R["iters"].max()}')
        print(f'  Max |u_inverse - u| (force NR)  : {R["inv_err"].max():.3e} m')
        print(f'  Reaction at u=GAP (before jump) : {R["V"][int(round(GAP_DIST / DINCR)) - 1]:.1f} N')
        print(f'  Final reaction                  : {R["V"][-1]:.1f} N')
        print(f'  Period range                    : {R["T"].min():.4f} - {R["T"].max():.4f} s')

    Ksh = (FU - FY) / (DSU - DY)
    print('\n--- Hand checks ---')
    print(f'  T elastic            = {2*np.pi*np.sqrt(M/KE):.4f} s')
    print(f'  T hardening (1 spr.) = {2*np.pi*np.sqrt(M/Ksh):.4f} s   (Ksh = {Ksh:.1f} N/m)')
    print(f'  T hardening (2 spr.) = {2*np.pi*np.sqrt(M/(2*Ksh)):.4f} s')
    f_gap = FY + Ksh * (GAP_DIST - DY)
    print(f'  Spring force at gap  = {f_gap:.1f} N  -> 2 springs in parallel = {2*f_gap:.1f} N')

    #%% RUN OPENSEES INELASTIC ANALYSIS
    OS = OPENSEES_PUSHOVER_CONTACT()
    if OS is not None:
        u_o, V_o, T_o = OS
        print('\n--- OpenSees vs Newton-Raphson (OPENSEES-equivalent) ---')
        print(f'  max |dV|/Vmax = {np.max(np.abs(V_o - R_os["V"])) / V_o.max():.3e}')
        print(f'  max |dT|/Tmax = {np.max(np.abs(T_o - R_os["T"])) / T_o.max():.3e}')

    #%% Plots
    fig, ax = plt.subplots(1, 2, figsize=(14, 6))
    #ax[0].plot(R_os['u'], R_os['V'], 'k', lw=3, label='NR - as coded in OpenSees (strain2 = u)')
    ax[0].plot(R_ph['u'], R_ph['V'], color='black', lw=3, label='Numerical')
    if OS is not None:
        ax[0].plot(OS[0], OS[1], 'r--', lw=3, label='OpenSees')
    ax[0].axvline(GAP_DIST, color='purple', ls='--')
    ax[0].set_xlabel('Displacement [m]')
    ax[0].set_ylabel('Total reaction [N]')
    ax[0].set_title('Pushover curve')
    ax[0].grid(True)
    ax[0].legend()

    #ax[1].plot(R_os['u'], R_os['T'], 'k', lw=3, label='NR - as coded')
    ax[1].plot(R_ph['u'], R_ph['T'], color='black', lw=3, label='Numerical')
    if OS is not None:
        ax[1].plot(OS[0], OS[2], 'r--', lw=3, label='OpenSees')
    ax[1].axvline(GAP_DIST, color='purple', ls='--')
    ax[1].set_xlabel('Displacement [m]')
    ax[1].set_ylabel('Period [s]')
    ax[1].set_title('Period shift'); ax[1].grid(True)
    ax[1].legend()
    plt.tight_layout() 
    plt.savefig('SDOF_CONTACT_NUMERICAL_OPENSEESPY.png', dpi=150); plt.show()
