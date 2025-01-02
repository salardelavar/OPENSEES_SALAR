#*****************************************************************%
#                   >> IN THE NAME OF ALLAH <<                    %
# Pushover Analysis of truss  subjected to lateral displacement   %
# Checking the analysis by Displacement - Large displacement      %
#          Verified by ABAQUS, SEISMSOSTRUCT, OpenSees            %
#-----------------------------------------------------------------%
#  This Program is Written by Salar Delavar Ghashghaei (Qashqai)  %
#              E-mail:salar.d.ghashghaei@gmail.com                %
#*****************************************************************%
import numpy as np

# Define parameters (units: mm, kN)
P5 = 0   # [kN] External Load [DOF (5)]
P6 = 0   # [kN] External Load [DOF (6)]
D5 = -1  # [mm] Initial displacement [DOF (5)] Incremental Displacement
D5max = 1000  # [mm] Maximum displacement [DOF (5)]
XY1i = np.array([0, 0])  # [x y] Point 1 Coordinate
XY2i = np.array([0, 500])  # [x y] Point 2 Coordinate
XY3i = np.array([500, 250])  # [x y] Point 3 Coordinate
A1 = np.pi * (50**2) / 4  # [mm^2]
A2 = np.pi * (50**2) / 4  # [mm^2]
E = 200  # [N/mm^2] Modulus of elasticity 
Nsteps =  int(np.abs(D5max / D5))# Number of steps for calculation

MAX_ITERATIONS = 5000 # convergence iteration for test
TOLERANCE = 1.0e-10   # convergence tolerance for test

L1i = np.sqrt((XY3i[0] - XY1i[0])**2 + (XY3i[1] - XY1i[1])**2)
L2i = np.sqrt((XY3i[0] - XY2i[0])**2 + (XY3i[1] - XY2i[1])**2)
EA1 = E * A1  # [kN]
EA2 = E * A2  # [kN]
u = np.zeros((1,))  # Initial guess value

# Initialize storage for results
lanX1L1, lanY1L1, lanX2L1, lanY2L1 = [], [], [], []
DU2, I2, IT2 = [], [], []
XY3ox2, XY3oy2 = [], []
ess1o2, ess2o2, LL1, LL2 = [], [], [], []
U1L, U2L, INT_L_f1, INT_L_f2 = [], [], [], []
INT_L_f1y, INT_L_f2y, TBSL1 = [], [], []

print("#####################################################")
print("# Python Code: Large Displacement Analysis [DOF(5)] #")
print("#####################################################")

# Gradually increase the applied displacement
for i in range(Nsteps):
    up = D5 * i
    XY1 = XY1i.copy()
    XY2 = XY2i.copy()
    XY3 = XY3i + np.array([up, u[0]])

    # Compute lengths and direction cosines
    L1 = np.sqrt((XY3[0] - XY1[0])**2 + (XY3[1] - XY1[1])**2)
    L2 = np.sqrt((XY3[0] - XY2[0])**2 + (XY3[1] - XY2[1])**2)
    lanX1, lanY1 = (XY3[0] - XY1[0]) / L1, (XY3[1] - XY1[1]) / L1
    lanX2, lanY2 = (XY3[0] - XY2[0]) / L2, (XY3[1] - XY2[1]) / L2

    # Element stiffness
    G1 = EA1 / L1
    G2 = EA2 / L2

    # Assemble global stiffness matrix
    Kp = np.array([
        [G1 * lanX1**2 + G2 * lanX2**2, G1 * lanX1 * lanY1 + G2 * lanX2 * lanY2],
        [G1 * lanX1 * lanY1 + G2 * lanX2 * lanY2, G1 * lanY1**2 + G2 * lanY2**2]
    ])
    Fii = Kp[:, 0] * up
    Kini = G1 * lanY1**2 + G2 * lanY2**2

    # Define the applied load
    Fi = np.array([P5, P6])
    F = Fi - Fii
    F = np.array([F[1]])

    # Iterative Newton-Raphson method
    it = 0
    residual = 100
    while residual > TOLERANCE:
        # Update lengths and direction cosines
        L1 = np.sqrt((XY3[0] - XY1[0])**2 + (XY3[1] - XY1[1])**2)
        L2 = np.sqrt((XY3[0] - XY2[0])**2 + (XY3[1] - XY2[1])**2)
        lanX1, lanY1 = (XY3[0] - XY1[0]) / L1, (XY3[1] - XY1[1]) / L1
        lanX2, lanY2 = (XY3[0] - XY2[0]) / L2, (XY3[1] - XY2[1]) / L2

        G1 = EA1 / L1
        G2 = EA2 / L2
        K = G1 * lanY1**2 + G2 * lanY2**2
        f = K * u - F

        # Calculate du and update residual
        du = -f / Kini
        residual = np.max(np.abs(du))
        it += 1

        if it == MAX_ITERATIONS:
            print(f"(-) For increment {i}, trial iterations reached to Ultimate {it}")
            print("    ## The solution for this step is not converged ##")
            break

        u += du

    if it < MAX_ITERATIONS:
        print(f"(+) Increment {i} : It is converged in {it} iterations")

    # Store results
    lanX1L1.append(lanX1)
    lanY1L1.append(lanY1)
    lanX2L1.append(lanX2)
    lanY2L1.append(lanY2)
    DU2.append(residual)
    I2.append(i)
    IT2.append(it)
    XY3ox2.append(XY3[0])
    XY3oy2.append(XY3[1])

    es1 = (L1 - L1i) / L1i
    es2 = (L2 - L2i) / L2i
    ess1o2.append(es1)
    ess2o2.append(es2)
    LL1.append(L1)
    LL2.append(L2)

    # Internal element forces
    INT_L_f1.append(-EA1 * es1 * lanX1)
    INT_L_f2.append(-EA2 * es2 * lanX2)
    INT_L_f1y.append(-EA1 * es1 * lanY1)
    INT_L_f2y.append(-EA2 * es2 * lanY2)
    TBSL1.append(INT_L_f1[-1] + INT_L_f2[-1])

    U1L.append(up)
    U2L.append(u[0])

    if abs(up) >= D5max:
        print("      ## Displacement at [DOF (5)] reached to Ultimate Displacement ##")
        break

# Convert results to arrays
D1 = np.array(U1L)# X-displacement [DOF 5]
F1 = np.array(TBSL1)

print("########################################################")
print("# OpensSees Code: Large Displacement Analysis [DOF(5)] #")
print("########################################################")

import openseespy.opensees as ops
from Analysis_Function import ANALYSIS

XY1 = [0, 0]  # [x y] Point 1
XY2 = [0, 500]  # [x y] Point 2
XY3 = [500, 250]  # [x y] Point 3

# Model setup
ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 2)

# Nodes
ops.node(1, *XY1)
ops.node(2, *XY2)
ops.node(3, *XY3)

# Boundary Conditions
ops.fix(1, 1, 1)
ops.fix(2, 1, 1)

# Materials
ops.uniaxialMaterial('Elastic', 1, EA1)
ops.uniaxialMaterial('Elastic', 2, EA2)

# Elements
ops.element('truss', 1, 1, 3, A1, 1)
ops.element('truss', 2, 2, 3, A2, 2)
#ops.element('corotTruss', 1, 1, 3, A1, 1)
#ops.element('corotTruss', 2, 2, 3, A2, 2)

# Time series and pattern
ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)
# Apply external force 
ops.load(3, P5, P6)

print('Model Built')

NstepGravity = 10
DGravity = 1 / NstepGravity
ops.integrator('LoadControl', DGravity) # determine the next time step for an analysis
ops.numberer('Plain') # renumber dof's to minimize band-width (optimization), if you want to
ops.system('BandGeneral') # how to store and solve the system of equations in the analysis
ops.constraints('Plain') # how it handles boundary conditions
ops.test('NormDispIncr', TOLERANCE, MAX_ITERATIONS) # determine if convergence has been achieved at the end of an iteration step
ops.algorithm('Newton') # use Newton's solution algorithm: updates tangent stiffness at every iteration
ops.analysis('Static') # define type of analysis static or transient
ops.analyze(NstepGravity) # apply gravity

ops.loadConst('-time', 0.0) #maintain constant gravity loads and reset time to zero

ops.timeSeries('Linear', 2)
ops.pattern('Plain', 200, 2)
ops.load(3, -1, 0.0)
    
# Analysis setup
ops.wipeAnalysis()
ops.system('BandGeneral')
ops.numberer('Plain')
ops.constraints('Plain')
ops.integrator('DisplacementControl', 3, 1, D5) # Apply displacement control at DOF 5 (node 3, X)
ops.algorithm('Newton')
ops.analysis('Static')

# Results storage
displacements = []
forces = []

# Perform analysis
displacement_X = 0
step = 0
while abs(displacement_X) < abs(D5max) and step < Nsteps:
    print(f'STEP {step+1}')
    stable = ops.analyze(1)
    #ANALYSIS(stable, 1, TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS
    # Record results
    displacement_X = ops.nodeDisp(3, 1) # Displacement in X Direction
    displacement_Y = ops.nodeDisp(3, 2) # Displacement in Y Direction
    #reaction_force = ops.nodeResponse(1, 1, 6) + ops.nodeResponse(2, 1, 6) # Reaction node 1 & 2 force in x-direction 
    displacements.append(displacement_X)
    #forces.append(reaction_force)
    step += 1
    
    # Update geometry
    Xnew = XY3[0] + displacement_X
    Ynew = XY3[1] + displacement_Y
    #print('Geometry:', Xnew)
    L1 = np.sqrt((Xnew - XY1[0])**2 + (Ynew - XY1[1])**2)
    L2 = np.sqrt((Xnew - XY2[0])**2 + (Ynew - XY2[1])**2)
    lanX1, lanY1 = (Xnew - XY1[0]) / L1, (Ynew - XY1[1]) / L1
    lanX2, lanY2 = (Xnew - XY2[0]) / L2, (Ynew - XY2[1]) / L2
    new_x3 = L1 * lanX1
    new_y3 = L1 * lanY1 
    es1 = (L1 - L1i) / L1i
    es2 = (L2 - L2i) / L2i
    reaction_force = (-EA1 * es1 * lanX1) + (-EA2 * es2 * lanX2)
    forces.append(reaction_force)
    
    """
    # Remove and redefine node and elements
    ops.remove('node', 3)
    ops.remove('element', 1)
    ops.remove('element', 2)
    
    ops.node(3, new_x3, new_y3)
    ops.element('truss', 1, 1, 3, A1, 1)
    ops.element('truss', 2, 2, 3, A2, 2)
    """
# Post-processing
displacements = np.array(displacements)
forces = np.array(forces)


### Plot the Results from the Output Data:
import matplotlib.pyplot as plt
# Plot Y-displacement (D6) vs Base Reaction (kN)
plt.figure(figsize=(10, 5))
plt.plot(D1, F1, label='Python')
plt.plot(displacements, forces, label='OpenSees', color='r', linestyle='--')
plt.xlabel('X-displacement (DOF 5) [mm]')
plt.ylabel('Base Reaction [kN]')
plt.title('X-displacement (DOF 5) vs Base Reaction')
plt.grid(True)
plt.legend()
plt.show()




