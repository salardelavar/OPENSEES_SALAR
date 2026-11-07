###########################################################################################################
#                   >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                      #
#   MONOTONIC PUSHOVER ANALYSIS OF STRUCTURE WITH FORCE ANALOGY METHOD (FAM) AND OPENSEESPY VERIFICATION  #
#---------------------------------------------------------------------------------------------------------#
#                 THIS PYTHON SCRIPT IS WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                     #
#                                   EMAIL: salar.d.ghashghaei@gmail.com                                   #
###########################################################################################################
"""
This python script performs a displacement-controlled pushover analysis of a fixed-guided column 
using the Force Analogy Method (FAM). 
1. The global stiffness of the structure remains constant at its initial elastic state (K_elastic).
2. Material nonlinearity (hinge formation) is captured by solving for plastic rotations (theta_p)
   at the base and top, which act as analog strain-like parameters reducing internal forces.
3. The moment and shear equations are decomposed into elastic and plastic components.
4. A return-mapping algorithm tracks the onset of yielding and computes the plastic rotations.
5. Three subplots are generated: Pushover Curve, Moment Evolution, and Plastic Rotation Evolution.
"""
# BOOK: Theory of Nonlinear Structural Analysis: The Force Analogy Method for Earthquake Engineering
# Gang Li, Kevin K. F. Wong (2014) - Wiley
'https://onlinelibrary.wiley.com/doi/book/10.1002/9781118718070'
# YOUTUBE: Pushover Analysis Force Analogy Method with Force Control Based on Timoshenko Beam Theory in C
'https://www.youtube.com/watch?v=u9M6klJO6uU&t=58s'


import numpy as np
import matplotlib.pyplot as plt
import openseespy.opensees as ops
import STEEL_ELASTIC_SECTION_ANALYSIS_FUN as S01
import ANALYSIS_FUNCTION as S02
import CONCRETE_ELASTIC_SECTION_ANALYSIS_FUN as S03
import COMPOSITE_ELASTIC_SECTION_ANALYSIS_FUN as S04
#%% Input Parameters (Units: N, mm)
# STEEL SECTION
#x_c, y_c, A_total, Ix_total, Iy_total = S01.STEEL_ELASTIC_SECTION_ANALYSIS_FUN()
# CONCRETE SECTION
#x_c, y_c, A_total, Ix_total, Iy_total = S02.CONCRETE_ELASTIC_SECTION_ANALYSIS_FUN()
# COMPOSITE SECTION
x_el, y_el, A_trans_total, Ix_trans_total, Iy_trans_total, E_ref = S04.COMPOSITE_ELASTIC_SECTION_ANALYSIS_FUN()

E = E_ref                 # [N/mm^2]
I = Ix_trans_total        # [mm^4]
A = A_trans_total         # [mm^2] (axial is irrelevant here)
L = 3000.0                # [mm]
Mp_base = 150e6           # [N.mm]
Mp_top = 250e6            # [N.mm]
D_target = 5.0            # [mm]
N_POINTS = 50             # FAM Numerical Increments

#%% FAM reference solution
def FAM(d):
    """
    FAM with proper hinge-sequence logic: candidate yield sets (none / base / top / both) are tested
    against the CORRECTED moments (after redistribution), not against the elastic trial moments.
    A set is admissible if the non-yielded moment stays <= its capacity and plastic rotations >= 0.
    """
    K_el = 12 * E * I / L**3
    K_d = np.array([6 * E * I / L**2] * 2)
    K_R = np.array([[4 * E * I / L, 2 * E * I / L], [2 * E * I / L, 4 * E * I / L]])
    K_p = np.array([6 * E * I / L**2] * 2)
    Mp = np.array([Mp_base, Mp_top])
    tol = 1e-9
    for active in ([], [0], [1], [0, 1]):
        th = np.zeros(2)
        if active:
            idx = np.array(active)
            th[idx] = np.linalg.solve(K_R[np.ix_(idx, idx)], K_d[idx] * d - Mp[idx])   # M_i = Mp_i
        M = K_d * d - K_R @ th
        free = [i for i in range(2) if i not in active]
        if all(abs(M[i]) <= Mp[i] * (1 + tol) for i in free) and np.all(th >= -tol):
            break
    return K_el * d - K_p @ th, M[0], M[1], th[0], th[1]



# RUN FORCE ANALOGY METHOD
disp_steps = np.linspace(0, D_target, N_POINTS)
fam = np.array([FAM(d) for d in disp_steps])    # np.array() is used for running faster      
#%% OpenSeesPy model
"""
Model (2D, ndm=2, ndf=3):
  - Fixed base, top node: rotation fixed + vertical fixed, horizontal displacement-controlled
    (fixed-guided column -> antisymmetric double curvature, K = 12EI/L^3).
  - Elastic beam-column + two concentrated zero-length rotational hinges (elastic-perfectly plastic,
    Steel01 with b ~ 0), using the Ibarra-Medina-Krawinkler (IMK) stiffness modification so the
    ELASTIC response is EXACTLY that of the original member (no artificial added flexibility):
        K_spring = (n+1) * 6EI/L      I_beam = I * (n+1)/n
  - Plastic hinge rotation  theta_p = |theta_spring| - |M| / K_spring   (comparable with FAM theta_p).

"""
MAX_ITERATIONS = 500     # Convergence iteration for test
MAX_TOLERANCE = 1.0e-6   # Convergence tolerance for test


# Reset model
ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 3)

n = 100.0                                  # IMK factor
K_spring = (n + 1.0) * 6.0 * E * I / L     # rotational spring stiffness [N.mm/rad]
I_beam = I * (n + 1.0) / n
b_hard = 1e-9                              # tiny hardening to avoid singular tangent

# Nodes: 1 (fixed base) | 2 (column bottom) | 3 (column top, loaded) | 4 (top rotation-fixed ground)
ops.node(1, 0.0, 0.0)
ops.node(2, 0.0, 0.0)
ops.node(3, 0.0, L)
ops.node(4, 0.0, L)

ops.fix(1, 1, 1, 1)        # fixed base
ops.fix(3, 0, 1, 0)        # beam top: vertical fixed, horizontal free (loaded node), rotation via spring
ops.fix(4, 0, 0, 1)        # spring ground side at top: rotation fixed (guided)
ops.equalDOF(1, 2, 1, 2)   # base hinge: translations tied, rotation free
ops.equalDOF(3, 4, 1, 2)   # top hinge: translations tied, rotation free

# Hinge materials (elastic - perfectly plastic moment-rotation)
ops.uniaxialMaterial('Steel01', 1, Mp_base, K_spring, b_hard)
ops.uniaxialMaterial('Steel01', 2, Mp_top, K_spring, b_hard)

ops.element('zeroLength', 1, 1, 2, '-mat', 1, '-dir', 3)   # base hinge (rot. dof)
ops.element('zeroLength', 2, 3, 4, '-mat', 2, '-dir', 3)   # top hinge

ops.geomTransf('Linear', 1)
ops.element('elasticBeamColumn', 3, 2, 3, A, E, I_beam, 1)

# Load pattern (unit reference load at top, displacement controlled)
ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)
ops.load(3, 1.0, 0.0, 0.0)

# Analysis
dU = D_target / (N_POINTS - 1)

ops.constraints('Transformation')
ops.numberer('RCM')
ops.system('BandGeneral')
ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS, 0) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/normDispIncr.html
ops.algorithm('Newton')                                    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/algorithm.html
ops.integrator('DisplacementControl', 3, 1, dU)

ops.analysis('Static')

def record():
    # base shear from the elastic column element (local y = shear for a vertical member)
    V = abs(ops.eleResponse(3, 'localForce')[1])
    # zeroLength 'force' = spring moment, 'deformation' = spring rotation (2D: dir 3 = rotation)
    Mb = ops.eleResponse(1, 'force')[5]
    Mt = ops.eleResponse(2, 'force')[5]
    rb = ops.eleResponse(1, 'deformation')[0]
    rt = ops.eleResponse(2, 'deformation')[0]
    Mb, Mt = abs(Mb), abs(Mt)
    tpb = max(abs(rb) - Mb / K_spring, 0.0)
    tpt = max(abs(rt) - Mt / K_spring, 0.0)
    return V, Mb, Mt, tpb, tpt, ops.nodeDisp(3, 1)

ops_res = [(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)]

for step in range(1, N_POINTS):
    OK = ops.analyze(1)
    S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS)
    ops_res.append(record())
    print(f'STEP: {step+1} - DISP.: {ops_res[-1][5]:.3f} mm - REACTION: {ops_res[-1][0]:.3f}')

ops_res = np.array(ops_res)
m = len(ops_res)
d_ops = ops_res[:, 5]

#%% Error summary
V_fam = fam[:m, 0]
err = np.abs(ops_res[:, 0] - V_fam)
rel = err[1:] / np.maximum(np.abs(V_fam[1:]), 1e-12) * 100
print('\n' + '=' * 60)
print('FORCE ANALOGY METHOD vs OPENSEES')
print('=' * 60)
print(f'Elastic stiffness  FAM : {12*E*I/L**3:.4f} N/mm')
print(f'Elastic stiffness  OPS : {ops_res[1,0]/d_ops[1]:.4f} N/mm')
print(f'Max base-shear error   : {err.max():.4e} N  ({rel.max():.4e} %)')
print(f'Final shear  FAM / OPS : {V_fam[-1]:.3f} / {ops_res[-1,0]:.3f} N')
print(f'Final theta_p base  FAM / OPS: {fam[m-1,3]:.6e} / {ops_res[-1,3]:.6e} rad')
print(f'Final theta_p top   FAM / OPS: {fam[m-1,4]:.6e} / {ops_res[-1,4]:.6e} rad')
print('=' * 60)

#%% Plots the Results
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5.5))

ax1.plot(disp_steps, fam[:, 0], color='navy', lw=3, label='FAM (corrected)')
ax1.plot(d_ops, ops_res[:, 0], 'o', color='crimson', ms=4, mfc='none', label='OpenSeesPy')
ax1.set_title('Pushover Curve (Base Shear vs. Disp.)')
ax1.set_xlabel('Lateral Displacement (mm)');
ax1.set_ylabel('Base Shear (N)')

ax2.plot(disp_steps, fam[:, 1], color='crimson', lw=3, label='FAM $M_{base}$')
ax2.plot(disp_steps, fam[:, 2], color='forestgreen', lw=3, label='FAM $M_{top}$')
ax2.plot(d_ops, ops_res[:, 1], 'k--', lw=3, label='OpenSees $M_{base}$')
ax2.plot(d_ops, ops_res[:, 2], 'b--', lw=3, label='OpenSees $M_{top}$')
ax2.set_title('Bending Moment Evolution')
ax2.set_xlabel('Lateral Displacement (mm)')
ax2.set_ylabel('Bending Moment (N.mm)')

ax3.plot(disp_steps, fam[:, 3], color='crimson', lw=3, label='FAM base')
ax3.plot(disp_steps, fam[:, 4], color='forestgreen', lw=3, label='FAM top')
ax3.plot(d_ops, ops_res[:, 3], 'k--', lw=3, label='OpenSees base')
ax3.plot(d_ops, ops_res[:, 4], 'b--', lw=3, label='OpenSees top')
ax3.set_title('Plastic Rotation Evolution')
ax3.set_xlabel('Lateral Displacement (mm)')
ax3.set_ylabel('Plastic Rotation (rad)')

for ax in (ax1, ax2, ax3):
    ax.grid(True, linestyle='--', alpha=0.6)
    ax.legend(fontsize=9)
    ax.set_xlim(0, D_target)

plt.tight_layout()
plt.show()