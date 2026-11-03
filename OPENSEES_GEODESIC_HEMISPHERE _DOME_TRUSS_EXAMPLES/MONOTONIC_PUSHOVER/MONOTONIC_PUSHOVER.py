######################################################################################################################
#                            >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                        #
# MONOTONIC PUSHOVER ANALYSIS OF GEODESIC HEMISPHERE DOME TRUSS WITH GEOMETRIC AND MATERIAL NONLINEARITY IN OPENSEES #
#--------------------------------------------------------------------------------------------------------------------#
#                           THIS PROGRAM IS WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                            #
#                                       EMAIL: salar.d.ghashghaei@gmail.com                                          #
######################################################################################################################
"""
This OpenSeesPy script performs a monotonic displacement-controlled pushover analysis of a geodesic hemisphere dome
 truss that accounts for both geometric and material nonlinearity. It first generates the geometry by creating a
 frequency-subdivided icosahedron, extracting the upper hemisphere nodes and truss edges, then builds a 3D OpenSees
 model in which base nodes are fixed, steel members are assigned either elastic or hysteretic inelastic material,
 and either ordinary Truss or CorotTruss (large-displacement) elements are used. A linearly varying vertical (or horizontal)
 reference load is applied to all non-base nodes, after which a DisplacementControl integrator drives the apex node to a
 prescribed target displacement while recording the corresponding base reaction; finally the script plots the original versus
 deformed geometry together with the resulting pushover curve.
"""
# WIKIPEDIA: Geodesic dome
'https://en.wikipedia.org/wiki/Geodesic_dome'
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import openseespy.opensees as ops
import ANALYSIS_FUNCTION as S02

#%% -------------------------------------------
# 1. Geometry Generation

def normalize(v):
    n = np.linalg.norm(v)
    return v / n if n > 1e-12 else v

def create_icosahedron():
    phi = (1.0 + np.sqrt(5.0)) / 2.0
    verts = np.array([
        [-1,  phi, 0], [ 1,  phi, 0], [-1, -phi, 0], [ 1, -phi, 0],
        [0, -1,  phi], [0,  1,  phi], [0, -1, -phi], [0,  1, -phi],
        [ phi, 0, -1], [ phi, 0,  1], [-phi, 0, -1], [-phi, 0,  1]
    ], dtype=float)
    verts = np.array([normalize(v) for v in verts])
    faces = np.array([
        [0,11,5], [0,5,1], [0,1,7], [0,7,10], [0,10,11],
        [1,5,9], [5,11,4], [11,10,2], [10,7,6], [7,1,8],
        [3,9,4], [3,4,2], [3,2,6], [3,6,8], [3,8,9],
        [4,9,5], [2,4,11], [6,2,10], [8,6,7], [9,8,1]
    ])
    return verts, faces

def subdivide_face(v1, v2, v3, freq):
    points = []
    for i in range(freq + 1):
        for j in range(freq + 1 - i):
            a = i / freq
            b = j / freq
            c = 1.0 - a - b
            p = a*v1 + b*v2 + c*v3
            points.append(normalize(p))
    points = np.array(points)
    faces = []
    def get_idx(row, col):
        return row*(freq+1) - (row*(row-1))//2 + col
    for i in range(freq):
        for j in range(freq - i):
            p1 = get_idx(i, j)
            p2 = get_idx(i, j+1)
            p3 = get_idx(i+1, j)
            faces.append([p1, p2, p3])
            if j < freq - i - 1:
                p4 = get_idx(i+1, j+1)
                faces.append([p2, p4, p3])
    return points, faces

def generate_geodesic_hemisphere(frequency=4, radius=1.0, z_threshold=0.0):
    base_verts, base_faces = create_icosahedron()
    all_points = []
    point_map = {}
    def add_point(p):
        key = tuple(np.round(p, decimals=12))
        if key not in point_map:
            point_map[key] = len(all_points)
            all_points.append(p.copy())
        return point_map[key]
    all_faces = []
    for face in base_faces:
        v1, v2, v3 = base_verts[face]
        local_pts, local_faces = subdivide_face(v1, v2, v3, frequency)
        local_to_global = [add_point(p) for p in local_pts]
        for lf in local_faces:
            all_faces.append([local_to_global[k] for k in lf])
    vertices = np.array(all_points) * radius
    mask = vertices[:, 2] >= z_threshold - 1e-8
    hemi_nodes = vertices[mask]
    old_to_new = {old: new for new, old in enumerate(np.where(mask)[0])}
    hemi_faces = []
    for f in all_faces:
        if all(mask[i] for i in f):
            hemi_faces.append([old_to_new[i] for i in f])
    edge_set = set()
    for f in hemi_faces:
        for a,b in [(f[0],f[1]),(f[1],f[2]),(f[2],f[0])]:
            edge_set.add(tuple(sorted((a,b))))
    edges = list(edge_set)
    return hemi_nodes, np.array(hemi_faces), edges

#%% -------------------------------------------
# 2. OpenSeesPy Model

def PUSHOVER_ANALYSIS(nodes, edges, radius, maxU, dU, SECTION_AREA, ELE_MAT, ELE_TYPE, DIREACTION):
    
    MAX_ITERATIONS = 50      # Convergence iteration for test
    MAX_TOLERANCE = 1.0e-6   # Convergence tolerance for test
        
    # Reset model
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    
    # Nodes
    for i, (x,y,z) in enumerate(nodes):
        ops.node(i+1, float(x), float(y), float(z))
    
    # Fix base nodes
    base_nodes = []
    for i, (x,y,z) in enumerate(nodes):
        if abs(z) < 1e-5:
            ops.fix(i+1, 1, 1, 1)
            base_nodes.append(i+1)
    print(f"Fixed base nodes: {len(base_nodes)}")
    
    # ---------- Steel Material ----------
    fy   = 240.0        # [N/mm²] Section Yield Stress
    E0   = 2.1e5        # [N/mm²] Section Elastic Modulus

    steelTag = 1
    
    if ELE_MAT == 'ELASTIC':
        ops.uniaxialMaterial('Elastic', steelTag, E0)
        # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/ElasticUni.html
    if ELE_MAT == 'INELASTIC': 
        """
        b    = 0.01         # Section Strain Hardening Ratio
        R0   = 18.0
        cR1  = 0.925
        cR2  = 0.15
        ops.uniaxialMaterial('Steel02', steelTag, Fy, E0, b, R0, cR1, cR2)
        # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/steel02.html
        """
        pinchX = 0.8   # Pinching factor in X direction
        pinchY = 0.5   # Pinching factor in Y direction
        damage1 = 0.0  # Damage due to ductility
        damage2 = 0.0  # Damage due to energy
        beta = 0.1     # Stiffness degradation parameter
        ey = fy/E0     # [mm/mm] Yield Strain
        fu = 1.18 * fy # [N/mm²] Ultimate Stress
        esu = ey * 50  # [mm/mm] Ultimate strain (at peak stress)

        ops.uniaxialMaterial('Hysteretic', steelTag,
                             fy, ey,
                             fu, esu,
                             0.2*fu, 1.1*esu,
                             -fy, -ey,
                             -fu, -esu,
                             -0.2*fu, -1.1*esu,
                             pinchX, pinchY,
                             damage1, damage2, beta)
        # INFO LINK: https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/Hysteretic.html
    
    # Truss elements
    for eid, (n1, n2) in enumerate(edges, start=1):
        if ELE_TYPE == 'TRUSS':
            ops.element('Truss', eid, n1+1, n2+1, SECTION_AREA, steelTag)
            # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/trussEle.html
        if ELE_TYPE == 'COROTTRUSS':
            ops.element('CorotTruss', eid, n1+1, n2+1, SECTION_AREA, steelTag)
            # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/corotTruss.html
    
    # Apex (control node)
    apex_id = int(np.argmax(nodes[:, 2])) + 1
    print(f"Control (Apex) node: {apex_id}")
    
    # Reference load pattern
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    for i, (x, y, z) in enumerate(nodes):
        if abs(z) > 1e-5:
            if DIREACTION == 'X':
                ops.load(i+1, z * 1.0, 0.0, 0.0)
            if DIREACTION == 'Y':
                ops.load(i+1, 0.0, z * 1.0, 0.0)        
            if DIREACTION == 'Z':
                ops.load(i+1, 0.0, 0.0, z * 1.0)   
                
    # Analysis settings
    ops.system('BandSPD')     # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/BandSPD.html
    #ops.system('BandGen')     # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/BandGen.html
    ops.numberer('RCM')       # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/RCM.html
    #ops.numberer('Plain')       # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/PlainNumberer.html
    ops.constraints('Plain')  # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/constraints.html
    ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS, 0) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/normDispIncr.html
    ops.algorithm('Newton') # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/algorithm.html
    
    # ---------- DisplacementControl ----------
    nsteps = int(abs(maxU / dU))
    
    if DIREACTION == 'X':
        DD = 1
    if DIREACTION == 'Y':
        DD = 2
    if DIREACTION == 'Z':
        DD = 3
        
    ops.integrator('DisplacementControl', apex_id, DD, dU)  # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/displacementControl.html
    ops.analysis('Static')   # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/analysis.html
    
    # Pushover loop
    disp_list = []
    reaction_list = []
    
    print(f"\nRunning Displacement-Controlled Pushover (target = {maxU:.0f} mm)...")
    print(f"{'Step':>6} {'Displacement (mm)':>15} {'Base Reaction (N)':>18}")
    print("-" * 60)
    
    for step in range(nsteps):
        OK = ops.analyze(1)
        S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS)

        ops.reactions()
        ux = ops.nodeDisp(apex_id, DD)
        rx = sum(ops.nodeReaction(bn, DD) for bn in base_nodes)
            
        disp_list.append(ux)
        reaction_list.append(rx)
        
        print(f"{step+1:6d} {ux:12.2f} mm,  {rx:18.2f} N")
    
    print("\nMonotonic Pushover finished.\n")
    
    # Deformed geometry
    deformed = np.zeros_like(nodes)
    for i in range(len(nodes)):
        deformed[i, 0] = nodes[i, 0] + ops.nodeDisp(i+1, 1)
        deformed[i, 1] = nodes[i, 1] + ops.nodeDisp(i+1, 2)
        deformed[i, 2] = nodes[i, 2] + ops.nodeDisp(i+1, 3)
    
    DATA = (np.array(disp_list), np.array(reaction_list), deformed, apex_id)
    
    return DATA

#%% -------------------------------------------
# 3. Plot Results

def PLOT_RESULTS(nodes, edges, deformed, disp, base_reaction, frequency, radius, DIREACTION):
    fig = plt.figure(figsize=(15, 10))
    
    max_disp = np.max(np.abs(deformed - nodes))
    scale = max(5.0, min(40.0, 0.8 / (max_disp + 1e-9)))
    def_vis = nodes + scale * (deformed - nodes)
    
    # 3D
    ax1 = fig.add_subplot(221, projection='3d')
    for a, b in edges:
        p1, p2 = nodes[a], nodes[b]
        ax1.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]], 'b-', lw=0.35, alpha=0.5)
        p1d, p2d = def_vis[a], def_vis[b]
        ax1.plot([p1d[0], p2d[0]], [p1d[1], p2d[1]], [p1d[2], p2d[2]], 'r-', lw=0.65)
    ax1.set_title(f'3D: Original (blue) / Deformed x{scale:.0f} (red)')
    ax1.view_init(elev=22, azim=40)
    ax1.set_xlabel('X'); ax1.set_ylabel('Y'); ax1.set_zlabel('Z')
    
    # Side view
    ax2 = fig.add_subplot(222)
    ax2.set_aspect('equal')
    for a, b in edges:
        p1, p2 = nodes[a], nodes[b]
        ax2.plot([p1[0], p2[0]], [p1[2], p2[2]], 'b-', lw=0.45, alpha=0.55)
        p1d, p2d = def_vis[a], def_vis[b]
        ax2.plot([p1d[0], p2d[0]], [p1d[2], p2d[2]], 'r-', lw=0.7)
    ax2.set_title('Side View (XZ)')
    ax2.set_xlabel('X (mm)'); ax2.set_ylabel('Z (mm)')
    ax2.grid(True, alpha=0.3)
    ax2.axhline(0, color='k', lw=0.8)
    
    # Pushover curve
    ax3 = fig.add_subplot(223)
    ax3.plot(disp, base_reaction, 'b-o', ms=3.5, lw=1.5)
    ax3.set_xlabel(f'Apex Displacement U{DIREACTION} (mm)')
    ax3.set_ylabel(f'Apex Base Reaction R{DIREACTION} (N)')
    ax3.set_title('Pushover Curve with Displacement Control')
    ax3.grid(True, alpha=0.4)
    
    # Top view
    ax4 = fig.add_subplot(224)
    ax4.set_aspect('equal')
    for a, b in edges:
        p1d, p2d = def_vis[a], def_vis[b]
        ax4.plot([p1d[0], p2d[0]], [p1d[1], p2d[1]], 'r-', lw=0.55)
    ax4.scatter(def_vis[:, 0], def_vis[:, 1], c='darkred', s=8, zorder=5)
    ax4.set_title('Top View - Deformed')
    ax4.set_xlabel('X (mm)'); ax4.set_ylabel('Y (mm)')
    ax4.grid(True, alpha=0.3)
    
    plt.suptitle(f'Geodesic {frequency}v Hemisphere | R={radius} mm | Displacement Control', fontsize=13)
    plt.tight_layout()
    plt.savefig(f'geodesic_{frequency}v_pushover_DC.png', dpi=160, bbox_inches='tight', facecolor='white')
    print(f"Plot saved -> geodesic_{frequency}v_pushover_DC.png")
    plt.show()

#%% -------------------------------------------
# 1. Geometry Generation
SHAPE_FREQUENCY = 4
"""
NOTICE: SHAPE_FREQUENCY Controls the subdivision frequency of the base icosahedron faces
 when generating the geodesic dome geometry (higher values produce finer meshes with more
 nodes and truss elements).
"""
RADIUS    = 10000.0   # [mm] Structure Radius
    
print("=" * 80)
print("  Geodesic Hemisphere Truss")
print("  OpenSeesPy  |  Steel Material  |  Displacement Control")
print("=" * 80)
    
nodes, faces, edges = generate_geodesic_hemisphere(frequency=SHAPE_FREQUENCY, radius=RADIUS)
print(f"\nGeometry: {len(nodes)} nodes | {len(edges)} truss elements | {SHAPE_FREQUENCY}v | R={RADIUS} mm")
    
np.savetxt("Geodesic_Nodes.txt", nodes, header="X Y Z (mm)", fmt="%.6f")

# 2. OpenSeesPy Model    
dU   = -0.01              # [mm] Incremental Displacement
maxU = -170.0             # [mm] Target Displacement
    
A = 1600.0                 # [mm²] Section Area
ELE_MAT = 'INELASTIC'      # SELECT MATERIAL  ->  'ELASTIC' OR 'INELASTIC'     
ELE_TYPE = 'COROTTRUSS'    # SELECT GEOMETRIC -> 'TRUSS' OR 'COROTTRUSS'

DIREACTION = 'X'           # SELECT MONOTONIC PUSHOVER DIREACTION 'X' OR 'Y' OR 'Z'
  
disp, reaction, deformed, apex = PUSHOVER_ANALYSIS(nodes, edges, RADIUS, maxU, dU, A, ELE_MAT, ELE_TYPE, DIREACTION)
   
# 3. Plot Results 
PLOT_RESULTS(nodes, edges, deformed, disp, reaction, SHAPE_FREQUENCY, RADIUS, DIREACTION)


# 4. Multilinear Curve Fitting
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
TRI_NODE = int(0.6 * len(disp))

QUAD_NODES = (int(0.6 * (len(disp))),
              int(0.65 * (len(disp))))

PENTA_NODES = (int(0.6 * (len(disp))),
               int(0.65 * (len(disp))),
               int(0.9 * (len(disp))))

#%% Fit curves
disp, reaction = np.abs(disp), np.abs(reaction)
X_bi, Y_bi, *_ = BILINEAR_CURVE(disp, reaction, SLOPE_NODE)
X_tri, Y_tri, *_ = TRILINEAR_CURVE(disp, reaction, SLOPE_NODE, TRI_NODE)
X_quad, Y_quad, *_ = QUADRILINEAR_CURVE(disp, reaction, SLOPE_NODE, *QUAD_NODES)
X_penta, Y_penta, *_ = PENTALINEAR_CURVE(disp, reaction, SLOPE_NODE, *PENTA_NODES)

#%% Plot
plt.figure(figsize=(10, 6))

plt.plot(disp, reaction, 'ko-', linewidth=2, markersize=5, label='Original curve')

plt.plot(X_bi, Y_bi, 'b--', linewidth=2, marker='s', markersize=5,
         label='Bilinear')

plt.plot(X_tri, Y_tri, 'g-.', linewidth=2, marker='^', markersize=5,
         label='Trilinear')

plt.plot(X_quad, Y_quad, 'm:', linewidth=2, marker='D', markersize=5,
         label='Quadrilinear')

plt.plot(X_penta, Y_penta, 'r-', linewidth=2, marker='o', markersize=5,
         label='Pentalinear')

plt.xlabel('Displacement [mm]')
plt.ylabel('Base Reaction [N]')
plt.title('Multilinear Fitting of Monotonic Pushover Curve')
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(f'Multilinear_Fitting_of_Monotonic_Pushover_Curve.png', dpi=160, bbox_inches='tight', facecolor='white')
print(f"Plot saved -> Multilinear_Fitting_of_Monotonic_Pushover_Curve.png")
plt.show()
    
print("\nAll done.")