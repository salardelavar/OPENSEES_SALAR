###########################################################################################################
#                          >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<               #
#                                     LINEAR ELASTIC BEAM ANALYSIS USING OPENSEES                         #
#---------------------------------------------------------------------------------------------------------#
#                          THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                     #
#                                   EMAIL: salar.d.ghashghaei@gmail.com                                   #
###########################################################################################################
"""
1. Model setup: Creates a 2D linear elastic beam of length 6 m divided into 699 elements, with specified material (E, A, Iz) and a uniform distributed load `w`.
2. Nodes and supports: Nodes are generated along the beam, with one end hinged and the other on a roller to create a simply supported boundary condition.
3. Elements and loading: Each segment is defined as an `elasticBeamColumn` element, and the uniform load is applied to all elements.
4. Static analysis: Performs a linear static analysis using `LoadControl` to calculate nodal displacements.
5. Post-processing: The function `PLOT_ELEMENT_FORCES_AND_DEFORMED` extracts element end forces (axial, shear, moment), interpolates transverse displacements with cubic beam shape functions, and plots both internal force diagrams and the deformed shape of the beam.
"""
import openseespy.opensees as ops
#import numpy as np
#import matplotlib.pyplot as plt

L  = 6.0          # [m] Beam length
E  = 2.1e11       # [Pa] Young's modulus
A  = 0.02         # [m^2] Cross-sectional area
Iz = 8e-4         # [m^4] Moment of inertia about the z-axis
w  = -10e3        # [N/m] Uniformly distributed load
NIDE_I, NIDE_J = 1, 700 # START , END NODE
nEle = NIDE_J-NIDE_I    # Number of elements 
npts = 30               # Number of element division      

#%%------------------------------------------------------------------------------------------------
# clear data
ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 3)
#%%------------------------------------------------------------------------------------------------
# Nodes
dx = L / nEle
for i in range(nEle + 1):
    ops.node(i+1, i*dx, 0.0)
#%%------------------------------------------------------------------------------------------------
# Boundary Conditions
ops.fix(1,        1, 1, 1)   # fixed
ops.fix(nEle+1,   0, 1, 0)   # roller
#%%------------------------------------------------------------------------------------------------
# Geometry Transformation
ops.geomTransf('Linear', 1)
#ops.geomTransf('PDelta', 1)
#ops.geomTransf('Corotational', 1)
#%%------------------------------------------------------------------------------------------------
# Elements
for i in range(nEle):
    ops.element(
        'elasticBeamColumn',
        i+1,
        i+1,
        i+2,
        A, E, Iz,
        1
    )

ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)
WA = []
for i in range(1, nEle+1):
    WA.append(w + w*i*2)
    #print(WA[-1])
    ops.eleLoad('-ele', i, '-type', '-beamUniform', w + w*i)

ops.system('BandGeneral')
ops.numberer('Plain')
ops.constraints('Plain')
ops.integrator('LoadControl', 1.0)
ops.algorithm('Linear')
ops.analysis('Static')

ops.analyze(1)
#%%------------------------------------------------------------------------------------------------
def PLOT_ELEMENT_FORCES_AND_DEFORMED(NIDE_I, NIDE_J, npts=30, scale=1.0, nEle=1):
    #npts: each element nodes
    import openseespy.opensees as ops
    import numpy as np
    import matplotlib.pyplot as plt
    x1, y1 = ops.nodeCoord(NIDE_I)
    x2, y2 = ops.nodeCoord(NIDE_J)

    v1  = ops.nodeDisp(NIDE_I,   2)
    th1 = ops.nodeDisp(NIDE_I,   3)
    v2  = ops.nodeDisp(NIDE_J, 2)
    th2 = ops.nodeDisp(NIDE_J, 3)
    
    L = np.sqrt((x2-x1)**2 + (y2-y1)**2)
    x = np.linspace(0, L, nEle + 1)
    
    N = np.zeros(nEle + 1)
    V = np.zeros(nEle + 1)
    M = np.zeros(nEle + 1)
    
    # Force at the left end of the element
    for i in range(NIDE_I, NIDE_J):
        f = ops.eleForce(i)
        N[i-1] = f[0]
        V[i-1] = f[1]
        M[i-1] = f[2]
    
    # Force at the Right end of the element
    f_last = ops.eleForce(nEle)
    N[-1] = f_last[3]
    V[-1] = f_last[4]
    M[-1] = f_last[5]
    
    plt.figure(figsize=(10,8))
    
    # ------------------
    # Distributed Load
    # ------------------
    plt.subplot(4,1,1)
    plt.plot(x[1:], WA, 'black', linewidth=2)
    plt.ylabel('DISTRIBUTED LOAD')
    plt.grid()
    
    # ------------------
    # Axial
    # ------------------
    plt.subplot(4,1,2)
    plt.plot(x, N, 'r', linewidth=2)
    plt.ylabel('INTERNAL AXIAL FORCE (N)')
    plt.grid()
    
    # ------------------
    # Shear
    # ------------------
    plt.subplot(4,1,3)
    plt.plot(x, V, 'b', linewidth=2)
    plt.ylabel('SHEAR FORCE (V)')
    plt.grid()
    
    # ------------------
    # Moment
    # ------------------
    plt.subplot(4,1,4)
    plt.plot(x, M, 'g', linewidth=2)
    plt.ylabel('MOMENT FORCE (M)')
    plt.xlabel('Element Length')
    plt.grid()
    
    plt.suptitle('Element Intenal Forces')
    plt.tight_layout()
    plt.show()    

    # Shape functions
    xi = np.linspace(0, 1, npts)
    N1 = 1 - 3*xi**2 + 2*xi**3
    N2 = xi - 2*xi**2 + xi**3
    N3 = 3*xi**2 - 2*xi**3
    N4 = -xi**2 + xi**3

    # Transverse displacement (local)
    v = (N1*v1 +
         N2*th1*L +
         N3*v2 +
         N4*th2*L)

    # Undeformed coordinates
    x_und = np.linspace(x1, x2, npts)
    y_und = np.linspace(y1, y2, npts)

    # Local normal
    nx =  (y2 - y1) / L
    ny = (x2 - x1) / L

    # Deformed coordinates
    x_def = x_und + scale * v * nx
    y_def = y_und + scale * v * ny
    
    # DEFORMED SHAPE
    plt.figure(figsize=(12,8))
    plt.plot(x_def, y_def, color='r', linewidth=2)
    # INTIAL SHAPE
    plt.plot([0, L], [0, 0], 'k--', linewidth=1)
    plt.title('Deformed Shape of Element')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.axis('equal')
    plt.grid()
    plt.show()



PLOT_ELEMENT_FORCES_AND_DEFORMED(NIDE_I, NIDE_J, npts=npts, scale=10.0, nEle=nEle)
#%%------------------------------------------------------------------------------------------------