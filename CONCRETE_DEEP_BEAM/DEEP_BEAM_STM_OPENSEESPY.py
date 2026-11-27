#######################################################################################
#          >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<           #
#        NONLINEAR STRUT-AND-TIE MODEL OF A CONCRETE DEEP BEAM IN OPENSEESPY          #
#-------------------------------------------------------------------------------------#
#   Simply supported deep beam, midspan point load, displacement-controlled           #
#   Struts: Concrete01 (compression only, softened)   Tie: Steel02 (bottom bars)      #
#   Units: N, mm, MPa                                                                 #
#-------------------------------------------------------------------------------------#
#       THIS PYTHON SCRIPT IS WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)           #
#######################################################################################
"""
 The script models a concrete deep beam as a three-node strut-and-tie truss in OpenSeesPy, in N, mm and MPa.
 The inputs at the top set the span, lever arm, beam thickness, effective strut width, concrete strength with
 a softening factor of 0.6, and the tie steel (4Ø20 bars). It then builds a 2D, 2-DOF-per-node model. Node 1 is
 a pin and node 2 a roller at the supports, and node 3 is the load point at top midspan. The two inclined
 `corotTruss` elements from the supports to the load point carry the concrete arch action, using `Concrete01`,
 which has no tension capacity. The bottom `corotTruss` between the supports carries the tension, using `Steel02`
 with strain hardening. A unit downward load is applied at node 3, and the analysis uses `DisplacementControl` on
 that node's vertical DOF, so the load rises until the tie yields and the full load-deflection curve is traced.
 Newton iteration is used, with a fallback to `ModifiedNewton` if a step fails. At each step the code records the midspan
 deflection, the applied load (taken from the sum of the support reactions) and the axial forces in the strut and tie.
 It then compares the result with hand calculations: the tie-yield load of 2·As·fy·tanθ (about 844 kN) and the strut-crushing
 load of 2·β·fc·A·sinθ (about 1124 kN). Because the tie governs, the OpenSees peak of about 848 kN, slightly above the tie-yield
 load because of steel hardening, shows a ductile flexural-type failure. Finally it plots three panels: load against deflection
 with the two capacity lines, member forces against deflection, and the original and deformed truss geometry.
"""
import numpy as np
import matplotlib.pyplot as plt
import openseespy.opensees as ops

#%% ---------------- GEOMETRY & MATERIAL INPUTS ----------------
L      = 2000.0      # [mm] span between supports
H_ARM  = 800.0       # [mm] internal lever arm z: tie centroid -> top node
B_W    = 200.0       # [mm] beam thickness
W_STR  = 250.0       # [mm] effective strut width  -> A_strut = W_STR * B_W

FC     = 30.0        # [MPa] concrete compressive strength
BETA_S = 0.60        # strut efficiency factor (cracked, bottle-shaped strut)
EPS_C0 = 0.002       # strain at peak
EPS_CU = 0.0035      # ultimate strain
FY, ES = 420.0, 200000.0   # steel yield (MPa) and modulus
AS_TIE = 4 * np.pi * 20**2 / 4   # 4 bars of 20 mm  (mm^2)
B_HARD = 0.01        # steel strain-hardening ratio

U_MAX, N_STEP = 12.0, 1200   # [mm] max midspan deflection, number of steps
dU = U_MAX / N_STEP

#%% ---------------- BUILD MODEL ----------------
ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 2)

# Nodes: 1 = left support, 2 = right support, 3 = load node (top, midspan)
ops.node(1, 0.0, 0.0)
ops.node(2, L,   0.0)
ops.node(3, L / 2, H_ARM)

ops.fix(1, 1, 1)     # pin
ops.fix(2, 0, 1)     # roller

# Materials
fpc = -BETA_S * FC
ops.uniaxialMaterial('Concrete01', 1, fpc, -EPS_C0, 0.2 * fpc, -EPS_CU)   # struts
ops.uniaxialMaterial('Steel02', 2, FY, ES, B_HARD, 18.5, 0.925, 0.15)     # tie

# Elements: truss(tag, iNode, jNode, A, matTag)
A_STRUT = W_STR * B_W
ops.element('corotTruss', 1, 1, 3, A_STRUT, 1)    # left  diagonal strut
ops.element('corotTruss', 2, 2, 3, A_STRUT, 1)    # right diagonal strut
ops.element('corotTruss', 3, 1, 2, AS_TIE, 2)     # bottom tie

# Load pattern (unit load, scaled by displacement control)
ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)
ops.load(3, 0.0, -1.0)

#%% ---------------- ANALYSIS ----------------
ops.system('BandGeneral')
ops.numberer('RCM')
ops.constraints('Plain')
ops.test('NormDispIncr', 1e-6, 50)
ops.algorithm('Newton')
ops.integrator('DisplacementControl', 3, 2, -dU)
ops.analysis('Static')

disp, load = [0.0], [0.0]
F_strut, F_tie = [0.0], [0.0]
for i in range(N_STEP):
    ok = ops.analyze(1)
    if ok != 0:                         # fall back to a stronger algorithm
        ops.algorithm('ModifiedNewton', '-initial')
        ok = ops.analyze(1)
        ops.algorithm('Newton')
    if ok != 0:
        print(f'Analysis stopped at step {i+1}, u = {-ops.nodeDisp(3, 2):.3f} mm')
        break
    ops.reactions()
    disp.append(-ops.nodeDisp(3, 2))
    load.append(-(ops.nodeReaction(1, 2) + ops.nodeReaction(2, 2)) * -1.0)   # total vertical reaction = applied P
    F_strut.append(ops.basicForce(1)[0])    # axial force (negative = compression)
    F_tie.append(ops.basicForce(3)[0])      # axial force (positive = tension)

disp, load = np.array(disp), np.array(load)
F_strut, F_tie = np.array(F_strut), np.array(F_tie)

#%% ---------------- HAND-CALC CHECK ----------------
theta = np.arctan(H_ARM / (L / 2))
P_tie   = 2 * AS_TIE * FY * np.tan(theta)                     # tie yields: T = P/(2 tan(theta))
P_strut = 2 * BETA_S * FC * A_STRUT * np.sin(theta)           # strut crushes: C = P/(2 sin(theta))
print(f'Strut angle      = {np.degrees(theta):.1f} deg')
print(f'P (tie yield)    = {P_tie/1e3:.1f} kN')
print(f'P (strut crush)  = {P_strut/1e3:.1f} kN')
print(f'OpenSees peak P  = {load.max()/1e3:.1f} kN  at u = {disp[load.argmax()]:.2f} mm')

#%% ---------------- PLOTS ----------------
fig, ax = plt.subplots(1, 3, figsize=(18, 5.5))

ax[0].plot(disp, load / 1e3, 'k', lw=2.5, label='OpenSeesPy STM')
ax[0].axhline(P_tie / 1e3, color='g', ls='--', label=f'Tie yield = {P_tie/1e3:.0f} kN')
ax[0].axhline(P_strut / 1e3, color='r', ls='--', label=f'Strut crushing = {P_strut/1e3:.0f} kN')
ax[0].set_xlabel('Midspan deflection [mm]'); ax[0].set_ylabel('Applied load P [kN]')
ax[0].set_title('Load - deflection'); ax[0].legend()

ax[1].plot(disp, F_tie / 1e3, 'g', lw=2, label='Tie (tension)')
ax[1].plot(disp, -F_strut / 1e3, 'r', lw=2, label='Strut (compression)')
ax[1].set_xlabel('Midspan deflection [mm]'); ax[1].set_ylabel('Member force [kN]')
ax[1].set_title('Member forces'); ax[1].legend()

# Truss geometry with final deformed shape (scaled)
nodes = {1: (0, 0), 2: (L, 0), 3: (L / 2, H_ARM)}
elems = [(1, 3, 'r'), (2, 3, 'r'), (1, 2, 'g')]
sc = 20.0
for i, j, c in elems:
    xi, yi = nodes[i]; xj, yj = nodes[j]
    ax[2].plot([xi, xj], [yi, yj], c + '--', lw=1, alpha=0.5)
    xi2 = xi + sc * ops.nodeDisp(i, 1); yi2 = yi + sc * ops.nodeDisp(i, 2)
    xj2 = xj + sc * ops.nodeDisp(j, 1); yj2 = yj + sc * ops.nodeDisp(j, 2)
    ax[2].plot([xi2, xj2], [yi2, yj2], c, lw=3)
ax[2].set_aspect('equal'); ax[2].set_title(f'STM geometry (dashed) & deformed shape (x{sc:.0f})')
ax[2].set_xlabel('x [mm]'); ax[2].set_ylabel('y [mm]')

for a in ax: a.grid(True, alpha=0.4)
plt.tight_layout()
plt.savefig('DEEP_BEAM_STM.png', dpi=150)
plt.show()
