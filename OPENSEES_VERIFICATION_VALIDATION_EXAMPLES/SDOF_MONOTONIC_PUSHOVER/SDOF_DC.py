###################################################################################################
#                >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                 #
#      MONOTONIC PUSHOVER ANALYSIS OF A SINGLE-DEGREE-OF-FREEDOM (SDOF) NONLINEAR SPRING WITH     #
#                                         DISPLACEMENT CONTROL                                    #
#-------------------------------------------------------------------------------------------------#
#                  FINITE ELEMENT METHOD IN PYTHON AND OPENSEESPY VERIFICATIONS                   #
#-------------------------------------------------------------------------------------------------#
#                   THIS PROGRAM IS WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                 # 
#                               E-mail:salar.d.ghashghaei@gmail.com                               #
###################################################################################################
"""
This script performs a displacement-controlled monotonic pushover analysis of a
 single-degree-of-freedom (SDOF) system made of parallel nonlinear springs, first 
 analytically in Python and then using OpenSeesPy for verification. It defines a
 multilinear force–displacement backbone (`Force` and `Displacement`) and per-spring
 coefficients (`Coff`), computes segment stiffnesses, and steps the imposed displacement up
 to the maximum specified displacement; at each step it calculates each spring’s tangent
 stiffness and force, sums them to obtain the total base reaction, and also constructs a
 linear-elastic comparison curve. It then plots residual–increment, iteration–increment,
 and base reaction–displacement diagrams. After that, it builds the same SDOF model in
 OpenSeesPy using `zeroLength` spring elements with either `HystereticSM` (inelastic multilinear)
 or `Elastic` materials, applies `DisplacementControl`, runs the pushover, and compares OpenSeesPy
 results with the analytical solution, reporting maximum absolute and relative errors. Finally,
 it post-processes the OpenSeesPy pushover curve by fitting area-preserving bilinear, trilinear,
 quadrilinear, and pentalinear approximations and reports elastic/plastic/tangent stiffness,
 ductility ratio, and overstrength factor.
"""
import numpy as np
import matplotlib.pyplot as plt
import PLOT_FORCE_DISP_FUN as S01
import ANALYSIS_FUNCTION as S02
#%% ---------------- Parameters ----------------
ANA_STEPS = 10000        # number of displacement increments

MAX_ITERATIONS = 500     # Convergence iteration for test
MAX_TOLERANCE = 1.0e-6   # Convergence tolerance for test
u = 0.0

#%% ---------------- Spring properties ----------------
# UNIT: Free
Force        = np.array([2000.0, 3000.0, 500.0, 200.0])
Displacement = np.array([0.01, 0.1, 0.12, 0.15])
Coff         = np.array([1.0, 1.0, 1.0, 1.0, 1.0])

S01.PLOT_FORCE_DISP_FUN(Displacement, Force)

Dmax  = float(np.max(Displacement))
u_inc = Dmax / ANA_STEPS            # imposed displacement increment

#%% ---------------- Segment slopes ----------------
Rk1 = (Force[0] - 0.0)       / (Displacement[0] - 0.0)
Rk2 = (Force[1] - Force[0])  / (Displacement[1] - Displacement[0])
Rk3 = (Force[2] - Force[1])  / (Displacement[2] - Displacement[1])
Rk4 = (Force[3] - Force[2])  / (Displacement[3] - Displacement[2])

#%% ---------------- Normalized force on multilinear curve ----------------
def F_normalized(u_val):
    """Force on the normalized (unit-coefficient) multilinear curve at u_val."""
    a = abs(u_val)
    s = np.sign(u_val) if u_val != 0.0 else 1.0
    if a <= Displacement[0]:
        return Rk1 * u_val
    elif a <= Displacement[1]:
        return s * (Force[0] + Rk2 * (a - Displacement[0]))
    elif a <= Displacement[2]:
        return s * (Force[1] + Rk3 * (a - Displacement[1]))
    elif a <= Displacement[3]:
        return s * (Force[2] + Rk4 * (a - Displacement[2]))
    else:
        return 0.0

#%% ---------------- Header for my own coding ----------------
print('#################################################')
print('#  Pushover Analysis of Nonlinear Springs       #')
print('#       (Displacement Control)                  #')
print('#################################################')

#%% ---------------- Nonlinear springs analysis ----------------
F1i, U1, DU1, I1, IT1 = [], [], [], [], []
f = np.zeros(5)
last_i = 0

for i in range(1, ANA_STEPS + 1):
    u = u_inc * i                      # ← imposed displacement

    it = 0
    residual = 100.0

    # For displacement control, the spring state is uniquely
    # determined by u; the loop below is kept for structural
    # similarity and converges in one iteration.
    while residual > MAX_TOLERANCE and it < MAX_ITERATIONS:
        Kt = np.zeros(5)
        a  = abs(u)
        for j in range(5):
            if a <= Displacement[0]:
                Kt[j] = Coff[j] * Rk1
            elif a <= Displacement[1]:
                Kt[j] = Coff[j] * (Force[0] + Rk2 * (a - Displacement[0])) / a
            elif a <= Displacement[2]:
                Kt[j] = Coff[j] * (Force[1] + Rk3 * (a - Displacement[1])) / a
            elif a <= Displacement[3]:
                Kt[j] = Coff[j] * (Force[2] + Rk4 * (a - Displacement[2])) / a
            else:
                Kt[j] = 0.0

        f_new = Kt * u
        residual = float(np.max(np.abs(f_new - f)))
        f = f_new
        it += 1

    if it >= MAX_ITERATIONS:
        print(f'(-)Increment {i:.0f}: reached max iterations {it}')
    else:
        print(f'(+)Increment {i:.0f}: converged in {it} iteration(s)')

    F_total = float(np.sum(f))         # ← total reaction force

    F1i.append(F_total)
    U1.append(u)
    DU1.append(residual)
    I1.append(i)
    IT1.append(it)
    last_i = i

    if abs(u) >= Dmax:
        print('  ## Displacement reached to ultimate displacement ##')
        break

D1 = np.concatenate(([0.0], np.array(U1)))
F1 = np.concatenate(([0.0], np.array(F1i)))

#%% ---------------- Linear springs (comparison) ----------------
print('#################################################')
print('#  Pushover Analysis of Linear Springs           #')
print('#################################################')

F2i, U2, DU2, I2, IT2 = [], [], [], [], []
K_linear = float(np.sum(Coff * Rk1))

for i in range(1, last_i + 1):
    u = u_inc * i
    F_total = K_linear * u

    F2i.append(F_total)
    U2.append(u)
    DU2.append(0.0)
    I2.append(i)
    IT2.append(1)

D2 = np.concatenate(([0.0], np.array(U2)))
F2 = np.concatenate(([0.0], np.array(F2i)))

#%% ---------------- Figure 1: Residual–Increment ----------------
plt.figure(2)
plt.plot(I1, DU1, color='black', linewidth=2, label='INELASTIC')
plt.plot(I2, DU2, color='green', linestyle='--', linewidth=2, label='ELASTIC')
plt.grid(True)
plt.xlabel('increment')
plt.ylabel('Residual')
plt.legend(loc='upper left', bbox_to_anchor=(1.02, 1))
plt.title('Residual–Increment diagram', color='b')

#%% ---------------- Figure 2: Iteration–Increment ----------------
plt.figure(3)
plt.plot(I1, IT1, color='black', linewidth=2, label='INELASTIC')
plt.plot(I2, IT2, color='green', linestyle='--', linewidth=2, label='ELASTIC')
plt.grid(True)
plt.xlabel('increment')
plt.ylabel('Iteration')
plt.legend(loc='upper left', bbox_to_anchor=(1.02, 1))
plt.title('Iteration–Increment diagram', color='b')

#%% ---------------- Figure 3: Force–Displacement ----------------
plt.figure(4)
plt.plot(D1, F1, color='black', linewidth=2, label='INELASTIC')
plt.plot(D2, F2, color='red', linestyle='--', linewidth=2, label='ELASTIC')
plt.grid(True)
plt.legend(loc='upper left', bbox_to_anchor=(1.02, 1))
plt.xlabel('Displacement')
plt.ylabel('Base Reaction')
plt.title('Base Reaction–Displacement Diagram (Displacement Control)', color='b')

plt.tight_layout()
plt.show()

#%% ------------------------------------------------------
"""
#*******************************************************************#
#  OpenSeesPy version of the Pushover Analysis of Nonlinear Springs #
#                     (Displacement Control)                        #
#                      OPENSEESPY VERSION                           #
#*******************************************************************#
"""
import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt

#%% ---------------- Parameters ----------------

"""
Force        = np.array([2000.0, 3000.0, 500.0, 200.0])
Displacement = np.array([0.01, 0.1, 0.12, 0.15])
Coff         = np.array([0.5, 0.8, 1.0, 0.8, 0.5])
"""

ANA_STEPS    = 500
Dmax  = float(np.max(Displacement))
u_inc = Dmax / ANA_STEPS
nSprings = len(Coff)

def SDOF_OPENSEES(ANA_STEPS, MAT_TYPE):
    #%% ---------------- Model ----------------
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
        Ftot = ops.nodeReaction(1, 1)       # reaction at fixed node

        U_os.append(u)
        F_os.append(Ftot)

        print(f'(+) Step {i:>6d}  u = {u:8.5f}  F_total = {Ftot:10.5f}')
        
    return U_os, F_os

# RUN OPENSEES INELASTIC ANALYSIS
MAT_TYPE = 'INELASTIC'    
U_os, F_os = SDOF_OPENSEES(ANA_STEPS, MAT_TYPE) # RUN THE FUNCTION

U_os = np.array(U_os)
F_os = np.array(F_os)

print(f'## Analysis finished at u = {U_os[-1]:.5f}, F = {F_os[-1]:.5f}')

# RUN OPENSEES ELASTIC ANALYSIS
MAT_TYPE = 'ELASTIC'    
U_el, F_el = SDOF_OPENSEES(ANA_STEPS, MAT_TYPE) # RUN THE FUNCTION

U_el = np.array(U_el)
F_el = np.array(F_el)




U_an = np.linspace(0.0, Dmax, ANA_STEPS + 1)
F_an = np.array([np.sum(Coff * F_normalized(u)) for u in U_an])

K_lin = float(np.sum(Coff * Rk1))
F_lin = K_lin * U_an

# Interpolate OpenSees results onto the analytical grid
F_os_i = np.interp(U_an, U_os, F_os)
err_abs = np.abs(F_os_i - F_an)
mask    = np.abs(F_an) > 1e-6
err_rel = np.max(err_abs[mask] / np.abs(F_an[mask]))

print('----------------------------------------------')
print(f'Max absolute error  : {err_abs.max():.6e}')
print(f'Max relative error  : {err_rel:.6e}')
print('----------------------------------------------')

#%% ---------------------------------------
plt.figure(1, figsize=(12, 8))
plt.plot(abs(U_os), abs(F_os), color='red',  linestyle='--', linewidth=2, label='OpenSeesPy Inelastic')
plt.plot(abs(U_el), abs(F_el), color='black', linestyle=':',linewidth=2.5, label='OpenSeesPy Elastic')
plt.grid(True)
plt.xlabel('Displacement')
plt.ylabel('Base Reaction')
plt.title('Base Reaction–Displacement — OpenSeesPy Pushover', color='b')
plt.legend(); plt.tight_layout(); plt.show()

plt.figure(2, figsize=(12, 8))
plt.plot(abs(D1), abs(F1), color='black', linewidth=2, label='Analytical - Inelastic')
plt.plot(abs(U_os), abs(F_os),  color='red',  linewidth=2, linestyle='--', label='OpenSeesPy - Inelastic')
plt.grid(True)
plt.xlabel('Displacement')
plt.ylabel('Base Reaction')
plt.title('Base Reaction–Displacement — Inelastic Pushover', color='b')
plt.legend(); plt.tight_layout(); plt.show()

plt.figure(3, figsize=(12, 8))
plt.plot(abs(D1), abs(F1), color='black', linewidth=2, label='Analytical - Inelastic')
plt.plot(abs(D2), abs(F2), color='blue', linestyle='--', linewidth=2, label='Analytical - Elastic')
plt.plot(abs(U_os), abs(F_os),  color='purple',  linewidth=2, label='OpenSeesPy - Inelastic')
plt.plot(abs(U_el), abs(F_el),  color='red', linestyle='--', linewidth=2, label='OpenSeesPy - Elastic')
plt.grid(True)
plt.xlabel('Displacement')
plt.ylabel('Base Reaction')
plt.title('Base Reaction–Displacement : Analytical vs OpenSeesPy', color='b')
plt.legend()
plt.tight_layout()
plt.show()

#%% ---------------------------------------
# ---------------------------------------------
#  Plot Base Reaction - Displacement Analysis 
# ---------------------------------------------

#%% Core fitting functions
def _area_under_curve(Cur, Mom):
    Cur = np.asarray(Cur, dtype=float)
    Mom = np.asarray(Mom, dtype=float)
    return np.sum((Mom[:-1] + Mom[1:]) * 0.5 * np.diff(Cur))


def _piecewise_linear_fit(Cur, Mom, fixed_nodes, SLOPE_NODE):
    """
    Fit a piecewise linear curve while preserving the area under the original curve.

    fixed_nodes : indices of intermediate breakpoints after yield.
                  Example: [] for bilinear, [8] for trilinear,
                           [7, 11] for quadrilinear, etc.
    SLOPE_NODE  : index used to compute the initial elastic slope k0.
    """
    Cur = np.asarray(Cur, dtype=float)
    Mom = np.asarray(Mom, dtype=float)

    if len(Cur) != len(Mom):
        raise ValueError("Cur and Mom must have the same length")
    if not (0 <= SLOPE_NODE < len(Cur)):
        raise ValueError("SLOPE_NODE out of range")
    if Cur[SLOPE_NODE] == 0:
        raise ValueError("Cur[SLOPE_NODE] must be non-zero")

    k0 = Mom[SLOPE_NODE] / Cur[SLOPE_NODE]
    Area = _area_under_curve(Cur, Mom)

    Cu = Cur[-1]
    Mu = Mom[-1]

    fixed_nodes = list(fixed_nodes)
    for n in fixed_nodes:
        if not (0 <= n < len(Cur) - 1):
            raise ValueError("fixed node must be a valid index before the last point")
    if any(fixed_nodes[i] >= fixed_nodes[i + 1] for i in range(len(fixed_nodes) - 1)):
        raise ValueError("fixed_nodes must be strictly increasing")

    X_fixed = [Cur[i] for i in fixed_nodes] + [Cu]
    Y_fixed = [Mom[i] for i in fixed_nodes] + [Mu]

    x2, y2 = X_fixed[0], Y_fixed[0]

    # Area from the first fixed breakpoint to the end
    A_rest = 0.0
    for j in range(len(X_fixed) - 1):
        A_rest += 0.5 * (Y_fixed[j] + Y_fixed[j + 1]) * (X_fixed[j + 1] - X_fixed[j])

    denom = k0 * x2 - y2
    if abs(denom) < 1e-12:
        raise ValueError("Cannot solve for first yield point: denominator near zero")

    x1 = (2.0 * Area - 2.0 * A_rest - y2 * x2) / denom
    y1 = k0 * x1

    if not (0.0 < x1 < x2):
        print(f"Warning: computed x1 = {x1:.4f} is not between 0 and x2 = {x2:.4f}")

    X = np.concatenate(([0.0, x1], X_fixed))
    Y = np.concatenate(([0.0, y1], Y_fixed))
    return X, Y, k0


def _postprocess(X, Y):
    Elastic_ST = Y[1] / X[1]
    Plastic_ST = Y[-1] / X[-1]
    Tangent_ST = np.diff(Y) / np.diff(X)
    Ductility_Rito = X[-1] / X[1]
    Over_Strength_Factor = Y[-1] / Y[1]

    print('+==========================+')
    print('=   Analysis curve fitted =')
    print('     Disp        Base Shear')
    print('----------------------------')
    print(np.column_stack((X, Y)))
    print('+==========================+')
    print('+----------------------------------------------------+')
    print(f' Structure Elastic Stiffness :     {Elastic_ST:.2f}')
    print(f' Structure Plastic Stiffness :     {Plastic_ST:.2f}')
    print(f' Structure Tangent Stiffness :     {np.array2string(Tangent_ST, precision=2)}')
    print(f' Structure Ductility Ratio :       {Ductility_Rito:.2f}')
    print(f' Structure Over Strength Factor:   {Over_Strength_Factor:.2f}')
    print('+----------------------------------------------------+')

    return X, Y, Elastic_ST, Plastic_ST, Tangent_ST, Ductility_Rito, Over_Strength_Factor


def MULTILINEAR_CURVE(Cur, Mom, SLOPE_NODE, FIXED_NODES):
    X, Y, k0 = _piecewise_linear_fit(Cur, Mom, FIXED_NODES, SLOPE_NODE)
    return _postprocess(X, Y)


def BILINEAR_CURVE(Cur, Mom, SLOPE_NODE):
    return MULTILINEAR_CURVE(Cur, Mom, SLOPE_NODE, [])


def TRILINEAR_CURVE(Cur, Mom, SLOPE_NODE, NODE2):
    return MULTILINEAR_CURVE(Cur, Mom, SLOPE_NODE, [NODE2])


def QUADRILINEAR_CURVE(Cur, Mom, SLOPE_NODE, NODE2, NODE3):
    return MULTILINEAR_CURVE(Cur, Mom, SLOPE_NODE, [NODE2, NODE3])


def PENTALINEAR_CURVE(Cur, Mom, SLOPE_NODE, NODE2, NODE3, NODE4):
    return MULTILINEAR_CURVE(Cur, Mom, SLOPE_NODE, [NODE2, NODE3, NODE4])


SLOPE_NODE = 1

# Choose breakpoint indices from the original curve
# %%
TRI_NODE = int(0.65 * len(U_os))

QUAD_NODES = (int(0.65 * (len(U_os))),
              int(0.70 * (len(U_os))))

PENTA_NODES = (int(0.65 * (len(U_os))),
               int(0.70 * (len(U_os))),
               int(0.80 * (len(U_os))))

#%% Fit curves
U_os, F_os = np.abs(U_os), np.abs(F_os)
X_bi, Y_bi, *_ = BILINEAR_CURVE(U_os, F_os, SLOPE_NODE)
X_tri, Y_tri, *_ = TRILINEAR_CURVE(U_os, F_os, SLOPE_NODE, TRI_NODE)
X_quad, Y_quad, *_ = QUADRILINEAR_CURVE(U_os, F_os, SLOPE_NODE, *QUAD_NODES)
X_penta, Y_penta, *_ = PENTALINEAR_CURVE(U_os, F_os, SLOPE_NODE, *PENTA_NODES)

#%% Plot
plt.figure(figsize=(10, 6))

plt.plot(U_os, F_os, 'ko-', linewidth=2, markersize=5, label='Original curve')

plt.plot(X_bi, Y_bi, 'b--', linewidth=2, marker='s', markersize=5,
         label='Bilinear')

plt.plot(X_tri, Y_tri, 'g-.', linewidth=2, marker='^', markersize=5,
         label='Trilinear')

plt.plot(X_quad, Y_quad, 'm:', linewidth=2, marker='D', markersize=5,
         label='Quadrilinear')

plt.plot(X_penta, Y_penta, 'r-', linewidth=2, marker='o', markersize=5,
         label='Pentalinear')

plt.xlabel('Displacement')
plt.ylabel('Base Reaction')
plt.title('Multilinear Fitting of Pushover Curve')
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()