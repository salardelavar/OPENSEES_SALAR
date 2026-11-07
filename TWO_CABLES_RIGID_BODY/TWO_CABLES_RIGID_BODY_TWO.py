###################################################################
# >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<#
# PUSHOVER ANALYSIS OF CABLE SUBJECTED TO LATERAL DISPLACEMENT,   #
# MATERIAL AND GEOMETRIC NONLINEARITY EFFECTS, AND DISPLACEMENT   #
#           CONTROL VERIFICATION USING OPENSEES                   #
#-----------------------------------------------------------------#
#  THIS PROGRAM IS WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)  #
#              E-MAIL: SALAR.D.GHASHGHAEI@GMAIL.COM               #
###################################################################
"""
# OpenSees & Python, Analysis of Cable Structures Incorporating Material and Geometric Nonlinearities

This Python code performs a pushover analysis of a cable subjected to lateral displacement
 using both custom calculations and the OpenSees library.
 
Key steps include:
1. Parameter Initialization: Sets geometry, material properties, and simulation parameters.
2. Custom Calculation: Uses Newton-Raphson iterations to calculate displacement and force responses for truss members based on incremental displacement.
3. OpenSees Modeling: Creates a finite element model, applies boundary conditions, loads, and performs displacement-controlled analysis.
4. Post-processing: Reads displacements and forces from results files.
5. Plotting Results: Compares Python and OpenSees results to validate the analysis.


Compares
  (1) NUMERICAL : exact kinematics of a rigid bar rotating about pin A (cables -> strain -> stress),
  (2) OPENSEES  : corotTruss cables (HystereticSM envelope) + bar modelled with elastic beam elements.

BEAM_TYPE = 'RIGID'     -> very stiff elements + Corotational transformation  (same assumptions as the
                           numerical model; this is the correct verification case)
BEAM_TYPE = 'NON-RIGID' -> 50x50 mm elastic beams, Linear transformation (your original OpenSees model;
                           the bar bends and node 2 is NOT forced onto the rigid-rotation circle)

Both runs use identical step sequence: u_k = k * D6, k = 0..N  (k = 0 is the undeformed state).
Requires: numpy, matplotlib, openseespy.  (PLOT_2D_TRUSS.py optional for the deformed-shape plot.)
"""
import numpy as np
import time as ti
import openseespy.opensees as ops
import ANALYSIS_FUNCTION as S01
import PLOT_2D_TRUSS as S02
import matplotlib.pyplot as plt

trapz = getattr(np, 'trapezoid', None) or np.trapz     # NumPy 1.x / 2.x compatibility
#%% ============================ INPUT ============================
BEAM_TYPE = 'RIGID'        # 'RIGID' or 'NON-RIGID'
D6 = -0.01                 # [mm] displacement increment at node C (negative = downward)
D6max = 600.0              # [mm] maximum displacement
SHOW_DEFORMED = True       # True -> use PLOT_2D_TRUSS.PLOT_2D_FRAME at the end of the OpenSees run
TOL, MAXIT = 1.0e-6, 100

# Geometry [mm]
A_pt = np.array([0.0, 0.0])      # A  pinned
D_pt = np.array([500.0, 500.0])  # D  fixed (top of cable 1)
E_pt = np.array([1000.0, 500.0]) # E  fixed (top of cable 2)
L_total, L1_0, L2_0 = 1000.0, 500.0, 500.0
A1 = A2 = np.pi * 50.0**2 / 4.0  # [mm^2]

# Material envelope (same data given to HystereticSM)  [kN, mm]
FY = 0.240                       # [kN/mm^2] Yield strength of steel section
FU = 1.1818 * FY                 # [kN/mm^2] Ultimate strength of steel section
Ke = 200.0                       # [kN/mm^2] Modulus of elasticity of steel section
DY = FY / Ke                     # [mm/mm] Yield steel strain
DSU = 0.36                       # [mm/mm] Ultimate steel strain
pos_eps = np.array([0.0, DY, DSU, 1.1 * DSU, 1.25 * DSU])
pos_sig = np.array([0.0, FY, FU, 0.2 * FU, 0.1 * FU])
all_eps = np.concatenate((-pos_eps[:0:-1], pos_eps))
all_sig = np.concatenate((-pos_sig[:0:-1], pos_sig))

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

plt.xlabel("Strain [mm/mm]")
plt.ylabel("Stress [kN]")
plt.title("Stress-Strain Curve")
plt.grid(True)
plt.axhline(0, linewidth=0.5)
plt.axvline(0, linewidth=0.5)
plt.show()

N = int(round(abs(D6max / D6)))
up = D6 * np.arange(0, N + 1)    # imposed vertical displacement history (k = 0 .. N)

#%% ======================= 1. NUMERICAL (vectorised) =======================
def NUMERICAL(up):
    theta = np.arcsin(up / L_total)
    c, s = np.cos(theta), np.sin(theta)
    B = np.column_stack((L1_0 * c, L1_0 * s))
    C = np.column_stack((L_total * c, L_total * s))
    L1 = np.linalg.norm(B - D_pt, axis=1)
    L2 = np.linalg.norm(C - E_pt, axis=1)
    e1, e2 = (L1 - L1_0) / L1_0, (L2 - L2_0) / L2_0
    N1 = np.interp(e1, all_eps, all_sig) * A1       # tension positive [kN]
    N2 = np.interp(e2, all_eps, all_sig) * A2
    Ry = -N1 * (B[:, 1] - D_pt[1]) / L1 - N2 * (C[:, 1] - E_pt[1]) / L2   # support reactions
    return dict(R=Ry, N1=N1, N2=N2, e1=e1, e2=e2, ux=C[:, 0] - L_total)

num = NUMERICAL(up)

#%% ======================= 2. OPENSEESPY MODEL =======================
#%% Model setup
# Reset model
ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 3)

# Nodes
for tag, xy in zip(range(1, 6), [(0, 0), (500, 0), (1000, 0), (500, 500), (1000, 500)]):
    ops.node(tag, *map(float, xy))

ops.fix(1, 1, 1, 0)
ops.fix(4, 1, 1, 1)
ops.fix(5, 1, 1, 1)

KP = [FY, DY, FU, DSU, 0.2 * FU, 1.1 * DSU, 0.1 * FU, 1.25 * DSU]
KN = [-v for v in KP]
ops.uniaxialMaterial('HystereticSM', 1, '-posEnv', *KP, '-negEnv', *KN, '-pinch', 1, 1)
ops.element('corotTruss', 1, 2, 4, A1, 1)    # cable 1
ops.element('corotTruss', 2, 3, 5, A2, 1)    # cable 2

if BEAM_TYPE == 'RIGID':
    Eb, Ab, Ib = 2.0e6, 5.0e4, 5.0e9                  # >> cable stiffness (EA/L = 785 kN/mm)
    ops.geomTransf('Corotational', 1)                 # large rotations (up to ~37 deg)
else:
    Eb, Ab, Ib = 200.0, 50.0 * 50.0, 50.0 * 50.0**3 / 12.0
    ops.geomTransf('Linear', 1)
ops.element('elasticBeamColumn', 3, 1, 2, Ab, Eb, Ib, 1)
ops.element('elasticBeamColumn', 4, 2, 3, Ab, Eb, Ib, 1)

ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)
ops.load(2, 0.0, -1.0, 0.0)          # reference load shape (as in the original script)
ops.load(3, 0.0, -1.0, 0.0)

ops.constraints('Plain')
ops.numberer('Plain')
ops.system('BandGeneral')
ops.test('NormDispIncr', TOL, MAXIT, 0)
ops.algorithm('Newton')
ops.integrator('DisplacementControl', 3, 2, D6)
ops.analysis('Static')

def solve_step(dU, depth=0):
    """One displacement increment with fall-back algorithms and automatic sub-stepping."""
    ops.integrator('DisplacementControl', 3, 2, dU)
    for alg in (('Newton',), ('ModifiedNewton',), ('KrylovNewton',), ('NewtonLineSearch',)):
        ops.algorithm(*alg)
        if ops.analyze(1) == 0:
            ops.algorithm('Newton')
            return 0
    ops.algorithm('Newton')
    if depth >= 6:
        return -1
    for _ in range(2):                                # split increment in two halves
        if solve_step(dU / 2.0, depth + 1) != 0:
            return -1
    ops.integrator('DisplacementControl', 3, 2, D6)
    return 0

def snapshot():
    ops.reactions()
    return (ops.nodeDisp(3, 2), ops.nodeDisp(3, 1),
            ops.nodeReaction(4, 2) + ops.nodeReaction(5, 2),
            ops.eleResponse(1, 'axialForce')[0], ops.eleResponse(2, 'axialForce')[0])

hist = [snapshot()]                                   # k = 0 (undeformed)
for k in range(1, N + 1):
    if solve_step(D6) != 0:
        print(f'*** OpenSees did not converge at step {k} (u = {D6 * k:.3f} mm) - comparison truncated')
        break
    ops.integrator('DisplacementControl', 3, 2, D6)
    hist.append(snapshot())

hist = np.array(hist)
m = len(hist)
u_ops, ux_ops, R_ops, N1_ops, N2_ops = hist.T

if SHOW_DEFORMED:
    import PLOT_2D_TRUSS as S02
    S02.PLOT_2D_FRAME(deformed_scale=1)

#%% ======================= 3. COMPARISON METRICS =======================
R_num, N1_num, N2_num = num['R'][:m], num['N1'][:m], num['N2'][:m]
Rpk = np.max(np.abs(R_num))

def metrics(a, b):                        # a = numerical (reference), b = OpenSees
    err = b - a
    ref = np.max(np.abs(a))
    ss_res, ss_tot = np.sum(err**2), np.sum((a - a.mean())**2)
    return dict(max_abs=np.max(np.abs(err)), rmse=np.sqrt(np.mean(err**2)),
                max_rel_pct=100 * np.max(np.abs(err)) / ref, r2=1 - ss_res / ss_tot)

mR, m1, m2 = metrics(R_num, R_ops), metrics(N1_num, N1_ops), metrics(N2_num, N2_ops)
E_num = trapz(np.abs(R_num), np.abs(up[:m]))
E_ops = trapz(np.abs(R_ops), np.abs(u_ops))
k0_num = abs(R_num[1] / up[1])
k0_ops = abs(R_ops[1] / u_ops[1])

print('\n' + '=' * 72)
print(f'NUMERICAL vs OPENSEES  |  BEAM_TYPE = {BEAM_TYPE}  |  steps compared: {m - 1} of {N}')
print('=' * 72)
print(f'{"":28s}{"Numerical":>14s}{"OpenSees":>14s}{"Diff [%]":>12s}')
print(f'{"Initial stiffness [kN/mm]":28s}{k0_num:14.4f}{k0_ops:14.4f}{100*(k0_ops-k0_num)/abs(k0_num):12.3f}')
print(f'{"Peak reaction [kN]":28s}{np.max(np.abs(R_num)):14.3f}{np.max(np.abs(R_ops)):14.3f}'
      f'{100*(np.max(np.abs(R_ops))-np.max(np.abs(R_num)))/np.max(np.abs(R_num)):12.3f}')
print(f'{"Disp. at peak [mm]":28s}{up[:m][np.argmax(np.abs(R_num))]:14.3f}{u_ops[np.argmax(np.abs(R_ops))]:14.3f}')
print(f'{"Area under curve [kN.mm]":28s}{E_num:14.2f}{E_ops:14.2f}{100*(E_ops-E_num)/E_num:12.3f}')
print('-' * 72)
for name, mm in (('Base reaction', mR), ('Cable 1 force', m1), ('Cable 2 force', m2)):
    print(f'{name:14s} max|err| = {mm["max_abs"]:.4e} kN | RMSE = {mm["rmse"]:.4e} kN | '
          f'max err = {mm["max_rel_pct"]:.4f} % of peak | R2 = {mm["r2"]:.8f}')
print('=' * 72)

#%% ======================= 4. PLOTS =======================
fig, ax = plt.subplots(2, 2, figsize=(16, 10))

ax[0, 0].plot(up[:m], R_num, 'k-', lw=4, label='Numerical')
ax[0, 0].plot(u_ops, R_ops, 'r--', lw=2.5, label='OpenSees')
ax[0, 0].set(title='Base Reaction vs. Displacement', xlabel='Y-Displacement at C (mm)', ylabel='Base Reaction (kN)')

ax[0, 1].plot(up[:m], N1_num, 'k-', lw=4, label='Cable 1 - Numerical')
ax[0, 1].plot(up[:m], N2_num, color='purple', lw=4, label='Cable 2 - Numerical')
ax[0, 1].plot(u_ops, N1_ops, 'c--', lw=2.5, label='Cable 1 - OpenSees')
ax[0, 1].plot(u_ops, N2_ops, 'r--', lw=2.5, label='Cable 2 - OpenSees')
ax[0, 1].set(title='Cable Axial Forces', xlabel='Y-Displacement at C (mm)', ylabel='Axial Force (kN)')

ax[1, 0].plot(up[:m], R_ops - R_num, 'r', lw=1.5, label='Reaction error (OpenSees - Numerical)')
ax[1, 0].plot(up[:m], N1_ops - N1_num, 'c', lw=1.0, label='Cable 1 force error')
ax[1, 0].plot(up[:m], N2_ops - N2_num, 'm', lw=1.0, label='Cable 2 force error')
ax[1, 0].set(title='Absolute Error', xlabel='Y-Displacement at C (mm)', ylabel='Error (kN)')

ax[1, 1].plot(up[:m], num['ux'][:m], 'k-', lw=4, label='Numerical')
ax[1, 1].plot(u_ops, ux_ops, 'r--', lw=2.5, label='OpenSees')
ax[1, 1].set(title='Horizontal Displacement of C (arc of rigid bar)', xlabel='Y-Displacement at C (mm)',
             ylabel='X-Displacement (mm)')

for a in ax.ravel():
    a.grid(True, linestyle='--', alpha=0.6)
    a.legend()
fig.suptitle(f'Numerical vs OpenSeesPy - {BEAM_TYPE} bar', fontsize=14)
plt.tight_layout()
plt.show()
