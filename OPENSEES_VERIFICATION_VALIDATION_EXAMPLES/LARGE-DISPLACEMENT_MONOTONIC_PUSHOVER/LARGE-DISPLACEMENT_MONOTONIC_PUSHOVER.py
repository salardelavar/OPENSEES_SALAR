###################################################################################################
#                >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                 #
#       LARGE-DISPLACEMENT MONOTONIC PUSHOVER OF A 2-BAR TRUSS: NUMERICAL NEWTON-RAPHSON METHOD   #
#                                    WITH OPENSEESPY VERIFICATION                                 #
#-------------------------------------------------------------------------------------------------#
#                     THIS PYTHON SCRIPT IS WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)         #
###################################################################################################
"""
 This script verifies a large-displacement two-bar truss pushover by solving
 it two independent ways and comparing them. The model has two pinned supports
 (nodes 1 and 2) and one free node 3 at (500, 250), connected by a flexible
 diagonal bar (element 1) and a very stiff vertical bar (element 2), with parameters
 taken from your original example. TRUSS_INTERNAL takes node 3’s displacement and computes,
 for each bar, its current length, unit direction and axial force N = EA·(L−L0)/L0.
 It sums these into the internal force vector and builds the consistent tangent matrix:
 the material part (EA/L0)·e·eᵀ plus the geometric part (N/L)·(I − e·eᵀ), which captures
 how bar rotation changes stiffness. TRUSS_NUMERICAL uses that function for a
 displacement-controlled Newton–Raphson: at each 1 mm increment it imposes the displacement
 on one DOF (Y for case ‘Y’, X for case ‘X’), then iterates on the other, free DOF until its
 force residual is zero, since no external load acts there. It stores the node displacements,
 the base reaction (−f in the imposed direction) and the two axial forces. TRUSS_OPENSEES builds 
 the same structure with two corotTruss elements, an Elastic material, a unit reference load along
 the imposed DOF, and the DisplacementControl integrator with Newton iterations. It falls back to
 your ANALYSIS() helper if an increment fails, then records the same five quantities from node displacements,
 reactions and basicForce. The main block runs both methods for each case, prints the maximum relative error
 for every quantity, and plots them in a 2×4 grid, with the numerical result as a thick black line and OpenSees
 as a red dashed line, so overlapping curves confirm the two methods agree.
"""
import numpy as np
import matplotlib.pyplot as plt

#%% ---------------------------------------------------------------
# PARAMETERS UNITS (mm, kN)
NODE_01 = [0.0, 0.0]                               # Node 1 X, Y Coordinate
NODE_02 = [500.0, 0.0]                             # Node 2 X, Y Coordinate
NODE_03 = [500.0, 250.0]                           # Node 3 X, Y Coordinate
SD01 = 50.0                                        # [mm] Section 01 Diameter
SD02 = 1000.0                                      # [mm] Section 02 Diameter

A1 = np.pi * SD01**2 / 4.0                         # [mm^2] element 1: node 1 -> node 3 (diagonal)
A2 = np.pi * SD02**2 / 4.0                         # [mm^2] element 2: node 2 -> node 3 (rigid vertical bar)
E = 200.0                                          # [kN/mm^2]  Elastic Modulus
XY_FIX = np.array([NODE_01, NODE_02])              # pinned supports (node 1, node 2)
XY3_0 = np.array(NODE_03)                          # free node 3 (initial position)
EA = np.array([E * A1, E * A2])
L0 = np.array([np.linalg.norm(XY3_0 - XY_FIX[0]), np.linalg.norm(XY3_0 - XY_FIX[1])])

#%% DEFINE ANALYSIS PROPERTIES
MAX_ITERATIONS = 50      # Convergence iteration for test
MAX_TOLERANCE = 1.0e-6   # Convergence tolerance for test
#%% ---------------------------------------------------------------
# NUMERICAL: large-displacement truss, Newton-Raphson, displacement control
def TRUSS_INTERNAL(u):
    """
    Computes the internal force vector, consistent tangent stiffness (material + geometric),
    and axial forces of both bars at the current displaced position of the free node.
    """
    # internal force vector and consistent tangent at node 3 for u = [ux, uy]
    f, K, N = np.zeros(2), np.zeros((2, 2)), np.zeros(2)
    p = XY3_0 + u                                  # current position of node 3
    for j in range(2):
        d = p - XY_FIX[j]
        L = np.linalg.norm(d)                      # current length
        e = d / L                                  # current unit direction (support -> node 3)
        N[j] = EA[j] * (L - L0[j]) / L0[j]         # axial force (tension positive)
        f += N[j] * e
        # material stiffness + geometric stiffness
        K += EA[j] / L0[j] * np.outer(e, e) + N[j] / L * (np.eye(2) - np.outer(e, e))
    return f, K, N

def TRUSS_NUMERICAL(case, n_inc=200, inc=0.1, MAX_ITERATIONS=50, MAX_TOLERANCE=1.0e-6):
    """
    TRUSS_NUMERICAL: Performs displacement-controlled Newton-Raphson iteration on the free DOF of
    the 2-bar truss (while imposing incremental displacement on the controlled DOF)
    and returns the history of node displacements, base reaction, and axial forces.
    """
    c = 1 if case == 'Y' else 0                    # imposed DOF (0 = X, 1 = Y)
    fr = 1 - c                                     # free DOF
    u = np.zeros(2)
    out = []
    for i in range(1, n_inc + 1):
        u[c] = i * inc                             # impose the displacement
        it = 0
        while True:
            f, K, N = TRUSS_INTERNAL(u)
            R = f[fr]                              # residual of the free DOF (no external load)
            if abs(R) / EA[0] <= MAX_TOLERANCE or it >= MAX_ITERATIONS:
                break
            u[fr] -= R / K[fr, fr]                 # Newton update of the free DOF
            it += 1
        if abs(R) / EA[0] > MAX_TOLERANCE:
            print(f'(-) Increment {i}: NOT converged')
            break
        out.append((u[0], u[1], -f[c], N[0], N[1]))
    
    DATA = np.array(out) # columns: ux, uy, base reaction, N1, N2
    return DATA    

#%% ---------------------------------------------------------------
# OPENSEESPY: corotTruss + DisplacementControl
def TRUSS_OPENSEES(case, n_inc=200, inc=0.1, MAX_ITERATIONS=200, MAX_TOLERANCE=1.0e-12):
    """
    TRUSS_OPENSEES: Builds the identical 2-bar truss in OpenSeesPy using `corotTruss` elements
    and `DisplacementControl`, then returns the same history of node displacements,
    base reaction and axial forces for verification against the numerical method.
    """
    import openseespy.opensees as ops
    import ANALYSIS_FUNCTION as S02
    
    # Reset model
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 2)
    ops.node(1, *XY_FIX[0])
    ops.node(2, *XY_FIX[1])
    ops.node(3, *XY3_0)
    
    ops.fix(1, 1, 1)
    ops.fix(2, 1, 1)
    
    ops.uniaxialMaterial('Elastic', 1, E)
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/ElasticUni.html
    
    ops.element('corotTruss', 1, 1, 3, A1, 1)      # large-displacement truss
    ops.element('corotTruss', 2, 2, 3, A2, 1)
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/corotTruss.html
              
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    # imposed DOF of node 3
    if case == 'X':
        dof = 1
        ops.load(3, 1.0, 0.0)
    if case == 'Y':
        dof = 2
        ops.load(3, 0.0, 1.0)   # reference load along imposed DOF only

    ops.constraints('Plain')                                # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/constraints.html
    ops.numberer('Plain')                                   # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/numberer.html
    ops.system('BandGeneral')                               # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/BandGen.html
    ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/normDispIncr.html
    ops.algorithm('Newton')                                 # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/algorithm.html
    ops.integrator('DisplacementControl', 3, dof, inc)      # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/displacementControl.html
    ops.analysis('Static')                                  # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/analysis.html

    out = []
    for i in range(1, n_inc + 1):
        OK = ops.analyze(1)
        if OK != 0:
            S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS)
            print(f'(-) Increment {i}: OpenSees did not converge'); break
        
        ops.reactions()
        out.append((ops.nodeDisp(3, 1), ops.nodeDisp(3, 2),
                    ops.nodeReaction(1, dof) + ops.nodeReaction(2, dof),   # base reaction
                    ops.eleResponse(1, 'basicForce')[0],                   # axial force ele 1
                    ops.eleResponse(2, 'basicForce')[0]))                  # axial force ele 2
    
    DATA = np.array(out)
    return DATA

#%% ---------------------------------------------------------------
# MAIN: compare both methods for the two pushover cases
if __name__ == '__main__':
    
    names = ['ux', 'uy', 'Base reaction', 'N1', 'N2']
    fig, ax = plt.subplots(2, 4, figsize=(20, 9))
    for r, case in enumerate(('Y', 'X')):          # 'Y': impose DOF 6, 'X': impose DOF 5
        
        #%% RUN NUMERICAL METHOD
        NUM = TRUSS_NUMERICAL(case)
        
        #%% RUN OPENSEES
        OS = TRUSS_OPENSEES(case)
        n = min(len(NUM), len(OS))
        err = [np.abs(NUM[:n, k] - OS[:n, k]).max() / np.abs(OS[:, k]).max() for k in range(5)]
        print(f'Case {case}: max relative error  ' +
              '  '.join(f'{nm} = {e:.2e}' for nm, e in zip(names, err)))

        imp = 1 if case == 'Y' else 0              # column of the imposed displacement
        
        #%% PLOT
        for k, col in enumerate([1 - imp, 2, 3, 4]):
            ax[r, k].plot(NUM[:, imp], NUM[:, col], 'k', lw=4, label='Numerical')
            ax[r, k].plot(OS[:, imp], OS[:, col], 'r--', lw=2.5, label='OpenSeesPy')
            ax[r, k].set_title(['Free-DOF displacement', 'Base reaction', 'Axial force ele. 1', 'Axial force ele. 2'][k])
            ax[r, k].set_xlabel(f'Imposed {"Y" if case == "Y" else "X"} displacement [mm]')
            ax[r, k].grid(True)
        ax[r, 0].legend()
    plt.tight_layout()
    plt.show()