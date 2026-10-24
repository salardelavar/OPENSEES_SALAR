###########################################################################################################
#                   >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                      #
#            MULTILINEAR CURVE FITTING AND PLOTTING FOR PUSHOVER AND MOMENT-CURVATURE ANALYSIS            #
#---------------------------------------------------------------------------------------------------------#
#                  THIS PYTHON SCRIPT IS WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                    #
#                                   EMAIL: salar.d.ghashghaei@gmail.com                                   #
###########################################################################################################
"""
This code generalizes a multilinear framework that fits bilinear, trilinear, quadrilinear,
 and pentalinear idealizations to a nonlinear pushover or moment-curvature curve by preserving
 the area under the original curve, solving for the first yield point through area equivalence,
 keeping user-selected intermediate breakpoints from the original data, computing elastic, plastic,
 and tangent stiffnesses as well as ductility ratio and over-strength factor, and finally plotting
 the original curve together with all fitted multilinear curves for visual comparison.
"""
import numpy as np
import matplotlib.pyplot as plt


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


#%% Example data: replace with your own Cur and Mom
Cur = np.array([
    0.0, 0.2, 0.5, 1.0, 1.8, 2.8, 4.0, 5.5, 7.0,
    8.5, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 25.0
])

Mom = np.array([
    0.0, 35.0, 80.0, 140.0, 200.0, 245.0, 275.0, 295.0, 308.0,
    316.0, 321.0, 325.0, 327.0, 270.0, 210.5, 150.8, 130.0, 110.2
])

# If you have a CSV file, use:
# data = np.loadtxt("data.csv", delimiter=",", skiprows=1)
# Cur = data[:, 0]
# Mom = data[:, 1]

SLOPE_NODE = 1

# Choose breakpoint indices from the original curve
TRI_NODE = 8
QUAD_NODES = (7, 11)
PENTA_NODES = (6, 9, 13)

#%% Fit curves
X_bi, Y_bi, *_ = BILINEAR_CURVE(Cur, Mom, SLOPE_NODE)
X_tri, Y_tri, *_ = TRILINEAR_CURVE(Cur, Mom, SLOPE_NODE, TRI_NODE)
X_quad, Y_quad, *_ = QUADRILINEAR_CURVE(Cur, Mom, SLOPE_NODE, *QUAD_NODES)
X_penta, Y_penta, *_ = PENTALINEAR_CURVE(Cur, Mom, SLOPE_NODE, *PENTA_NODES)

#%% Plot
plt.figure(figsize=(10, 6))

plt.plot(Cur, Mom, 'ko-', linewidth=2, markersize=5, label='Original curve')

plt.plot(X_bi, Y_bi, 'b--', linewidth=2, marker='s', markersize=5,
         label='Bilinear')
plt.plot(X_tri, Y_tri, 'g-.', linewidth=2, marker='^', markersize=5,
         label='Trilinear')
plt.plot(X_quad, Y_quad, 'm:', linewidth=2, marker='D', markersize=5,
         label='Quadrilinear')
plt.plot(X_penta, Y_penta, 'r-', linewidth=2, marker='o', markersize=5,
         label='Pentalinear')

plt.xlabel('Displacement / Curvature')
plt.ylabel('Base Shear / Moment')
plt.title('Multilinear Fitting of Pushover / Moment-Curvature Curve')
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()
