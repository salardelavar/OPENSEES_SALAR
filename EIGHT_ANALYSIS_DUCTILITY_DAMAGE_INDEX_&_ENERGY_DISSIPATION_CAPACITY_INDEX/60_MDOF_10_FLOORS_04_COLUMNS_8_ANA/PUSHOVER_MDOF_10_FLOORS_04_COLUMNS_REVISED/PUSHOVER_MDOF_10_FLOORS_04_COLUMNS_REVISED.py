###########################################################################################################
#                   >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                      #
#   NONLINEAR PUSHOVER AND TIME-HISTORY VALIDATION OF A 10-STORY SHEAR BUILDING: MDOF-TO-SDOF REDUCTION   #
#                      USING N2 AND DISPLACEMENT-BASED METHODS IN OPENSEESPY                              #
#---------------------------------------------------------------------------------------------------------#
#                  THIS PYTHON SCRIPT IS WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                    #
#                                   EMAIL: salar.d.ghashghaei@gmail.com                                   #
###########################################################################################################
"""
 This script builds a 10-story shear-building model in OpenSeesPy (one degree of freedom per
 floor, with four parallel columns per story, each modeled as a zero-length element with a
 HystereticSM material that has yield, peak, and two-stage strength-degrading branches)
 and uses it to compare two ways of reducing a multi-degree-of-freedom (MDOF) structure to
 an equivalent single-degree-of-freedom (SDOF) system. It first runs an eigenvalue analysis
 to get periods and the first mode shape, then performs a displacement-controlled pushover
 using a modal load pattern, tracking story-1 drift until the target is reached or
 convergence fails (with fallback algorithms at each step). From the resulting capacity
 curve, it applies the N2/ATC-40 transformation, which uses a fixed first-mode shape to 
 compute the participation factor Γ, equivalent mass m*, and SDOF curve. It also applies
 a displacement-based design (DBD) approach, which uses the actual displaced shape at each
 step to compute effective displacement, mass, stiffness, and period. The SDOF curve is
 idealized as an equal-energy elastic-perfectly-plastic bilinear curve, truncated at 20%
 strength loss after the peak, from which yield point, ductility, and overstrength are
 extracted. A ductility-based damage index per story is computed from interstory drift.
 Finally, the script validates the SDOF idealization by running nonlinear time-history
 analyses on both the full MDOF model (with Rayleigh damping calibrated to modes 1 and 3)
 and the equivalent SDOF system under a synthetic band-limited ground motion, then compares
 the peak roof displacements (scaling the SDOF result by Γ) and plots the pushover curves,
 effective periods, displaced shapes, damage indices, story shears, and time histories.

Units: N, m, kg, s
"""
import numpy as np
import matplotlib.pyplot as plt
import openseespy.opensees as ops

#%% ----------------------------------------------
# 1) INPUT DATA

NUM_FLOOR = 10
NUM_COL = 4
TOTAL_MASS = 500000.0                                   # [kg]
G = 9.81                                                # [m/s2]

FYi = np.array([85000.0, 95000.0, 75000.0, 65000.0])    # [N] Column yield force
FUi = np.array([1.18, 1.12, 1.20, 1.10]) * FYi          # [N] Column ultimate force
Kei = np.array([4.5e6, 4.1e6, 4.6e6, 4.3e6])            # [N/m] Column elastic stiffness
DSUi = np.array([0.32, 0.34, 0.30, 0.36])               # [m] Column ultimate drift
DYi = FYi / Kei                                         # [m] Column yield drift
MASS = np.full(NUM_FLOOR, TOTAL_MASS / NUM_FLOOR)       # [kg] Floor masses

RESID_1, RESID_2 = 0.20, 0.10                           # residual strength ratios after DSU
DAMPING = 0.05                                          # Rayleigh damping ratio

CTRL_NODE = 2                                           # pushover control node (story-1 drift)
DU_STEP = 0.01                                          # [m] control increment
D_TARGET = 1.3 * DSUi.max()                             # [m] target story-1 drift
LOSS = 0.20                                             # ultimate point = 20 % strength loss after the peak

#%% ----------------------------------------------
# 2) MODEL

def envelope(FY, DY, FU, DSU):
    """HystereticSM envelope points (force, disp) for one column, positive and negative branch."""
    kp = [FY, DY, FU, DSU, RESID_1 * FU, 1.1 * DSU, RESID_2 * FU, 1.25 * DSU]
    kn = [-k for k in kp]
    return kp, kn


def build_mdof():
    """10-floor shear building: node 1 fixed, floors = nodes 2..11, 4 parallel columns per story."""
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)
    ops.node(1, 0.0)
    ops.fix(1, 1)
    for k in range(NUM_FLOOR):
        ops.node(k + 2, 0.0)
        ops.mass(k + 2, MASS[k])                       # FIX: mass on the floor node, not on node k+1
    ele_by_story = []
    tag = 0
    for ii in range(NUM_FLOOR):
        tags = []
        for jj in range(NUM_COL):
            tag += 1
            kp, kn = envelope(FYi[jj], DYi[jj], FUi[jj], DSUi[jj])
            ops.uniaxialMaterial('HystereticSM', 1000 + tag, '-posEnv', *kp, '-negEnv', *kn, '-pinch', 1, 1)
            ops.element('zeroLength', tag, ii + 1, ii + 2, '-mat', 1000 + tag, '-dir', 1)
            tags.append(tag)
        ele_by_story.append(tags)
    return ele_by_story


def first_modes(n=3):
    """Eigen-analysis of the initial (elastic) structure. Returns omega [rad/s] and mode shapes (roof = 1)."""
    lam = np.array(ops.eigen('-fullGenLapack', n))
    omega = np.sqrt(lam)
    phi = np.zeros((n, NUM_FLOOR))
    for m in range(n):
        v = np.array([ops.nodeEigenvector(k + 2, m + 1, 1) for k in range(NUM_FLOOR)])
        phi[m] = v / v[-1]
    return omega, phi


def solve_step():
    """One analysis step with algorithm fall-back. Returns True if converged."""
    for alg in (('Newton',), ('ModifiedNewton', '-initial'), ('KrylovNewton',), ('NewtonLineSearch',)):
        ops.algorithm(*alg)
        if ops.analyze(1) == 0:
            ops.algorithm('Newton')
            return True
    ops.algorithm('Newton')
    return False


#%% ----------------------------------------------
# 3) PUSHOVER  (modal load pattern, control on story-1 drift)

def pushover(phi1):
    ele_by_story = build_mdof()
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    for k in range(NUM_FLOOR):
        ops.load(k + 2, MASS[k] * phi1[k] / MASS.sum())           # P_i = m_i * phi_i (normalised)
    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.system('BandGeneral')
    ops.test('NormDispIncr', 1.0e-9, 200)
    ops.integrator('DisplacementControl', CTRL_NODE, 1, DU_STEP)
    ops.analysis('Static')

    U = [np.zeros(NUM_FLOOR)]                                     # floor displacements per step
    VB = [0.0]
    VS = [np.zeros(NUM_FLOOR)]                                    # story shears
    n_steps = int(D_TARGET / DU_STEP)
    for s in range(n_steps):
        if not solve_step():
            print(f'Pushover stopped at step {s}, story-1 drift = {ops.nodeDisp(CTRL_NODE, 1):.4f} m')
            break
        ops.reactions()
        U.append(np.array([ops.nodeDisp(k + 2, 1) for k in range(NUM_FLOOR)]))
        VB.append(-ops.nodeReaction(1, 1))
        # story shear = sum of the 4 column forces of the SAME story
        VS.append(np.array([-sum(ops.eleResponse(t, 'force')[0] for t in tags) for tags in ele_by_story]))
    U, VB, VS = np.array(U), np.array(VB), np.array(VS)
    print(f'Pushover done: {len(VB)-1} steps, peak base shear = {VB.max()/1e3:.1f} kN')
    return U, VB, VS


#%% ----------------------------------------------
# 4) MDOF -> SDOF TRANSFORMATIONS

def n2_transformation(U, VB, phi1):
    """N2 / ATC-40 / ASCE-41: fixed shape phi_1.  Returns D*, F*, Gamma, m*."""
    m_star = np.sum(MASS * phi1)
    gamma = m_star / np.sum(MASS * phi1 ** 2)
    D_star = U[:, -1] / gamma                                     # roof displacement / Gamma
    F_star = VB / gamma
    return D_star, F_star, gamma, m_star


def dbd_properties(U, VB):
    """Displacement-based design: actual displaced shape at every step.
       Delta_eff = sum(m u^2)/sum(m u);  m_eff = sum(m u)/Delta_eff;  K_eff = Vb/Delta_eff."""
    mu = U @ MASS
    mu2 = (U ** 2) @ MASS
    with np.errstate(divide='ignore', invalid='ignore'):
        d_eff = np.where(mu > 0, mu2 / mu, np.nan)
        m_eff = np.where(d_eff > 0, mu / d_eff, np.nan)
        k_eff = np.where(d_eff > 0, VB / d_eff, np.nan)
        t_eff = 2 * np.pi * np.sqrt(m_eff / k_eff)
    return d_eff, m_eff, k_eff, t_eff


#%% ----------------------------------------------
# 5) EQUAL-ENERGY BILINEAR FIT (ultimate point at 20 % strength loss after the peak)

def bilinear_fit(D, F, k0):
    """Equal-energy elastic-perfectly-plastic idealisation (ASCE 41 / ASTM E2126 type).
    The curve is truncated at the ultimate point = first step after the peak with LOSS (20 %) strength drop,
    k0 = initial stiffness, Fy follows from equal areas:  Fy = k0*(Du - sqrt(Du^2 - 2*A/k0))."""
    i_peak = int(np.argmax(F))
    below = np.where(F[i_peak:] <= (1.0 - LOSS) * F[i_peak])[0]
    i_u = i_peak + (int(below[0]) if below.size else len(F) - 1 - i_peak)
    d, f = D[:i_u + 1], F[:i_u + 1]
    area = np.sum(0.5 * (f[1:] + f[:-1]) * np.diff(d))
    Du = d[-1]
    Fy = k0 * (Du - np.sqrt(Du ** 2 - 2.0 * area / k0))
    Dy = Fy / k0
    X, Y = np.array([0.0, Dy, Du]), np.array([0.0, Fy, Fy])
    out = dict(Dy=Dy, Fy=Fy, Du=Du, Fmax=F[i_peak], Ke=k0, mu=Du / Dy, omega_s=F[i_peak] / Fy, i_u=i_u)
    return X, Y, out


#%% ----------------------------------------------
# 6) DAMAGE INDEX PER STORY (ductility based, from drift)

def damage_index(U):
    drift = np.diff(np.hstack([np.zeros((U.shape[0], 1)), U]), axis=1)
    dy, du = DYi.mean(), DSUi.mean()
    di = np.clip(100.0 * (np.abs(drift) - dy) / (du - dy), 0.0, 100.0)
    return drift, di


#%% ----------------------------------------------
# 7) VALIDATION: NONLINEAR TIME HISTORY  MDOF  vs  EQUIVALENT SDOF

def synthetic_record(pga_g=0.40, duration=20.0, dt=0.01, seed=7):
    """Band-limited (0.3-3 Hz) filtered noise with an envelope. Only used if no record file is available."""
    from scipy.signal import butter, sosfiltfilt
    rng = np.random.default_rng(seed)
    n = int(duration / dt)
    t = np.arange(n) * dt
    sos = butter(2, [0.3, 3.0], btype='bandpass', fs=1.0 / dt, output='sos')
    ag = sosfiltfilt(sos, rng.standard_normal(n))
    ag *= np.minimum(t / 2.0, 1.0) * np.exp(-np.maximum(t - 10.0, 0.0) / 4.0)
    DATA = ag / np.max(np.abs(ag)) * pga_g * G, dt
    return DATA


def run_transient(ag, dt, roof_node, kind):
    """Run a transient analysis on the model already built. Returns time, roof disp (None for failure)."""
    ops.timeSeries('Path', 2, '-dt', dt, '-values', *ag)
    ops.pattern('UniformExcitation', 2, 1, '-accel', 2)
    ops.constraints('Plain') 
    ops.numberer('Plain') 
    ops.system('BandGeneral')
    ops.test('NormDispIncr', 1.0e-8, 200)
    ops.integrator('Newmark', 0.5, 0.25)
    ops.analysis('Transient')
    t_hist, u_hist = [0.0], [0.0]
    for _ in range(len(ag)):
        ok = 0
        for alg in (('Newton',), ('ModifiedNewton', '-initial'), ('KrylovNewton',)):
            ops.algorithm(*alg)
            ok = ops.analyze(1, dt)
            if ok == 0:
                break
        if ok != 0:
            print(f'{kind}: analysis failed at t = {ops.getTime():.2f} s (collapse / non-convergence)')
            break
        t_hist.append(ops.getTime()); u_hist.append(ops.nodeDisp(roof_node, 1))
    
    DATA = np.array(t_hist), np.array(u_hist)
    return DATA


def mdof_time_history(ag, dt, omega):
    build_mdof()
    w1, w3 = omega[0], omega[2]                                   # Rayleigh on INITIAL stiffness (modes 1 and 3)
    a0 = 2 * DAMPING * w1 * w3 / (w1 + w3)
    a1 = 2 * DAMPING / (w1 + w3)
    ops.rayleigh(a0, 0.0, a1, 0.0)
    
    DATA = run_transient(ag, dt, NUM_FLOOR + 1, 'MDOF')
    return DATA

def sdof_time_history(ag, dt, m_star, Dy, Fy, Du, gamma):
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)
    ops.node(1, 0.0); ops.node(2, 0.0); ops.fix(1, 1); ops.mass(2, m_star)
    kp = [Fy, Dy, Fy, Du, RESID_1 * Fy, 1.1 * Du, RESID_2 * Fy, 1.25 * Du]   # plateau to Du, then strength loss
    kn = [-k for k in kp]
    ops.uniaxialMaterial('HystereticSM', 1, '-posEnv', *kp, '-negEnv', *kn, '-pinch', 1, 1)
    ops.element('zeroLength', 1, 1, 2, '-mat', 1, '-dir', 1)
    w_star = np.sqrt((Fy / Dy) / m_star)
    ops.rayleigh(0.0, 0.0, 2 * DAMPING / w_star, 0.0)             # c = 2*zeta*m*omega on initial stiffness
    t, d = run_transient(ag, dt, 2, 'SDOF')
    return t, d * gamma                                           # back to roof displacement


#%% ----------------------------------------------
# 8) MAIN

if __name__ == '__main__':
    # --- modal properties of the initial structure
    build_mdof()
    omega, phi = first_modes(3)
    phi1 = phi[0]
    print('Periods [s]:', np.round(2 * np.pi / omega, 4))
    print('Mode-1 shape (floor 1..10):', np.round(phi1, 3))

    # --- pushover and transformation
    U, VB, VS = pushover(phi1)
    D_star, F_star, gamma, m_star = n2_transformation(U, VB, phi1)
    # initial SDOF stiffness from the first steps of the curve (elastic range)
    k0 = F_star[3] / D_star[3]
    X, Y, bl = bilinear_fit(D_star, F_star, k0)
    T_star = 2 * np.pi * np.sqrt(m_star / bl['Ke'])
    d_eff, m_eff, k_eff, t_eff = dbd_properties(U, VB)
    i_y = int(np.argmin(np.abs(D_star - bl['Dy'])))
    i_u = bl['i_u']

    print('\n================ N2 (fixed shape) ================')
    print(f"Gamma = {gamma:.4f}   m* = {m_star:,.0f} kg   T* = {T_star:.3f} s")
    print(f"Dy* = {bl['Dy']*1e3:.2f} mm  Fy* = {bl['Fy']/1e3:.1f} kN  Du* = {bl['Du']*1e3:.2f} mm  Fmax* = {bl['Fmax']/1e3:.1f} kN")
    print(f"Ductility = {bl['mu']:.2f}   Over-strength (Fmax/Fy) = {bl['omega_s']:.2f}")
    print('\n================ DBD (actual displaced shape) ================')
    for name, i in (('yield', i_y), ('ultimate', i_u)):
        print(f"{name:9s}: Delta_eff = {d_eff[i]*1e3:7.2f} mm   m_eff = {m_eff[i]:,.0f} kg   "
              f"K_eff = {k_eff[i]/1e6:6.2f} MN/m   T_eff = {t_eff[i]:.3f} s")

    drift, DI = damage_index(U)
    print(f"\nMax story damage index at the end of the pushover: {DI[i_u].max():.1f} % (story {int(DI[i_u].argmax())+1})")

    # --- plots: pushover
    fig, ax = plt.subplots(2, 2, figsize=(14, 10))
    ax[0, 0].plot(U[:, -1], VB / 1e3, 'k', lw=2)
    ax[0, 0].set(xlabel='Roof displacement [m]', ylabel='Base shear [kN]', title='MDOF pushover (modal pattern)')
    ax[0, 1].plot(D_star, F_star / 1e3, 'k', lw=2, label='SDOF capacity (N2)')
    ax[0, 1].plot(X, Y / 1e3, 'r--o', lw=2, label='Equal-energy EPP idealisation')
    ax[0, 1].set(xlabel="D* [m]", ylabel="F* [kN]", title='Equivalent SDOF'); ax[0, 1].legend()
    ax[1, 0].plot(U[:, -1], t_eff, lw=2, label='T_eff (DBD)')
    ax[1, 0].axhline(T_star, color='r', ls='--', label='T* (N2, elastic)')
    ax[1, 0].set(xlabel='Roof displacement [m]', ylabel='Period [s]', title='Effective period', ylim=(0, 2 * T_star)); ax[1, 0].legend()
    for i, lab in ((i_y, 'yield'), (i_u, 'ultimate')):
        ax[1, 1].plot(np.hstack([0, U[i]]), np.arange(NUM_FLOOR + 1), 'o-', label=f'{lab}')
    ax[1, 1].plot(np.hstack([0, phi1 * U[i_y, -1]]), np.arange(NUM_FLOOR + 1), 'k--', label='phi_1 scaled to yield')
    ax[1, 1].set(xlabel='Displacement [m]', ylabel='Floor', title='Displaced shapes'); ax[1, 1].legend()
    for a in ax.flat:
        a.grid(alpha=0.3)
    fig.tight_layout()

    fig2, ax2 = plt.subplots(1, 2, figsize=(14, 5))
    ax2[0].plot(U[:, -1], DI.max(axis=1), 'k', lw=2)
    ax2[0].set(xlabel='Roof displacement [m]', ylabel='Max story damage index [%]', title='Damage index')
    for k in range(NUM_FLOOR):
        ax2[1].plot(U[:, -1], VS[:, k] / 1e3, lw=1, label=f'story {k+1}')
    ax2[1].set(xlabel='Roof displacement [m]', ylabel='story shear [kN]', title='story shears'); ax2[1].legend(ncol=2, fontsize=7)
    for a in ax2:
        a.grid(alpha=0.3)
    fig2.tight_layout()

    # --- validation with a nonlinear time history (replace with your record: np.loadtxt('Ground_Acceleration_X.txt') * G)
    ag, dt = synthetic_record(pga_g=0.40)
    t1, u1 = mdof_time_history(ag, dt, omega)
    t2, u2 = sdof_time_history(ag, dt, m_star, bl['Dy'], bl['Fy'], bl['Du'], gamma)
    p1, p2 = np.max(np.abs(u1)), np.max(np.abs(u2))
    print('\n================ VALIDATION (NLTHA) ================')
    print(f"Peak roof displacement  MDOF = {p1*1e3:.1f} mm   SDOF = {p2*1e3:.1f} mm   SDOF/MDOF = {p2/p1:.2f}")

    fig3, ax3 = plt.subplots(figsize=(14, 5))
    ax3.plot(t1, u1 * 1e3, 'k', lw=1.5, label='MDOF roof')
    ax3.plot(t2, u2 * 1e3, 'r--', lw=1.5, label='Equivalent SDOF x Gamma')
    ax3.set(xlabel='Time [s]', ylabel='Roof displacement [mm]', title='Time-history validation'); ax3.legend(); ax3.grid(alpha=0.3)
    fig3.tight_layout()
    plt.show()
