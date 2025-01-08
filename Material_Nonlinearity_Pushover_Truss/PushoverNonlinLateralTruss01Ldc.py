###################################################################
#                   >> IN THE NAME OF ALLAH <<                    #
#  PUSHOVER ANALYSIS OF TRUSS SUBJECTED TO LATERAL DISPLACEMENT   #
#           MATERIAL AND GEOMETRIC NONLINEARITIES EFFECT          #
#          CHECKING THE ANALYSIS BY DISPLACEMENT CONTROL          #
#                       VERIFIED BY OPENSEES                      #
#-----------------------------------------------------------------#
#  THIS PROGRAM IS WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)  #
#              E-MAIL: SALAR.D.GHASHGHAEI@GMAIL.COM               #
###################################################################

import numpy as np
import time as ti

# Define parameters (units: mm, kN)
P5 = 0                       # [kN] External Load [DOF (5)]
P6 = 0                       # [kN] External Load [DOF (6)]
D5 = -0.0005                 # [mm] Initial Displacement [DOF (5)] Incremental Displacement
D5max = 20                   # [mm] Maximum displacement [DOF (5)]
XY1i = np.array([0, 0])      # [x y] Point 1 Coordinate
XY2i = np.array([0, 500])    # [x y] Point 2 Coordinate
XY3i = np.array([500, 250])  # [x y] Point 3 Coordinate
A1 = np.pi * (50**2) / 4     # [mm^2] Cross Section Area for Element 1
A2 = np.pi * (50**2) / 4     # [mm^2] Cross Section Area for Element 2

Nsteps =  int(np.abs(D5max / D5))# Number of steps for calculation

# Steel Section Properties
fy = 0.24       # [kN/mm^2] Yield strength of steel section
Es = 200        # [kN/mm^2] Modulus of elasticity of steel section
fu = 1.5 * fy   # Ultimate steel stress
ey = fy / Es    # Yield steel strain
esu = 0.35      # Ultimate steel strain

Esh = (fu - fy) / (esu - ey)  # Strain hardening modulus
b = Esh / Es                  # Strain hardening ratio


MAX_ITERATIONS = 10000 # convergence iteration for test
TOLERANCE = 1.0e-12    # convergence tolerance for test

# Truss Length
L1i = np.sqrt((XY3i[0] - XY1i[0])**2 + (XY3i[1] - XY1i[1])**2)
L2i = np.sqrt((XY3i[0] - XY2i[0])**2 + (XY3i[1] - XY2i[1])**2)

u = np.array([0.0])  # Initial guess value

# Analysis Durations:
starttime = ti.process_time()

print("###############")
print("# Python Code #")
print("###############")
print("--------------------------------------------------")
print("            Increments               Iterations   ")
print("--------------------------------------------------")

# Initialize storage for results
u = np.zeros(1)  # initial guess value
lanX1L1, lanY1L1, lanX2L1, lanY2L1 = [], [], [], []
DU2, I2, IT2 = [], [], []
ess1o2, ess2o2 = [], []
XY3ox2, XY3oy2 = [], []
ess1o2, ess2o2, LL1, LL2 = [], [], [], []
FS1o2, FSTAN1o2, FS2o2, FSTAN2o2 = [], [], [], []
U1L, U2L = [], []
Null2 = []
U5N = []
TBS5N = []

# Gradually increase the applied displacement
for i in range(Nsteps):
    up = D5 * i  # Incremental displacement for each step
    XY1 = XY1i.copy()  # Copy initial coordinates of point 1
    XY2 = XY2i.copy()  # Copy initial coordinates of point 2
    XY3 = XY3i + np.array([up, u[0]])  # Update coordinates of point 3
    
    # Compute lengths and direction cosines
    L1 = np.sqrt((XY3[0] - XY1[0])**2 + (XY3[1] - XY1[1])**2)  # Length of element 1
    L2 = np.sqrt((XY3[0] - XY2[0])**2 + (XY3[1] - XY2[1])**2)  # Length of element 2
    lanX1, lanY1 = (XY3[0] - XY1[0]) / L1, (XY3[1] - XY1[1]) / L1  # Direction cosines for element 1
    lanX2, lanY2 = (XY3[0] - XY2[0]) / L2, (XY3[1] - XY2[1]) / L2  # Direction cosines for element 2
    
    # Element strain
    es1 = (L1 - L1i) / L1i # Strain of element 1
    es2 = (L2 - L2i) / L2i # Strain of element 2
     
    # Stress-strain relationship for element 1
    if 0 <= es1 <= ey:
        fs1 = Es * es1
        fstan1 = Es
    elif -ey <= es1 < 0:
        fs1 = Es * es1
        fstan1 = Es
    elif ey < es1 <= esu:
        fs1 = fy + Esh * (abs(es1) - ey)
        fstan1 = (fy + Esh * (abs(es1) - ey)) / abs(es1)
    elif -esu <= es1 < -ey:
        fs1 = -fy - Esh * (abs(es1) - ey)
        fstan1 = (fy + Esh * (abs(es1) - ey)) / abs(es1)
    else:
        fs1, fstan1 = 0, 0

    # Stress-strain relationship for element 2
    if 0 <= es2 <= ey:
        fs2 = Es * es2
        fstan2 = Es
    elif -ey <= es2 < 0:
        fs2 = Es * es2
        fstan2 = Es
    elif ey < es2 <= esu:
        fs2 = fy + Esh * (abs(es2) - ey)
        fstan2 = (fy + Esh * (abs(es2) - ey)) / abs(es2)
    elif -esu <= es2 < -ey:
        fs2 = -fy - Esh * (abs(es2) - ey)
        fstan2 = (fy + Esh * (abs(es2) - ey)) / abs(es2)
    else:
        fs2, fstan2 = 0, 0
   
    # Element stiffness
    G1 = fstan1 * A1 / L1 # Stiffness of element 1
    G2 = fstan2 * A2 / L2 # Stiffness of element 2

    # Assemble global stiffness matrix
    Kp = np.array([
        [G1 * lanX1 ** 2 + G2 * lanX2 ** 2, G1 * lanX1 * lanY1 + G2 * lanX2 * lanY2],
        [G1 * lanX1 * lanY1 + G2 * lanX2 * lanY2, G1 * lanY1 ** 2 + G2 * lanY2 ** 2]
    ])
    Fii = Kp[:, 0] * up # Internal force vector
    Kini = np.array([G1 * lanY1 ** 2 + G2 * lanY2 ** 2]) # Initial stiffness

    # Define applied load
    Fi = np.array([P5, P6]) # External load vector
    F = Fi - Fii            # Net force vector
    F = np.array([F[1]])    # Select relevant component of force

    # Iterative Newton-Raphson method
    it, residual = 0, 100
    while residual > TOLERANCE:
        # Update lengths and direction cosines
        L1 = np.sqrt((XY3[0] - XY1[0])**2 + (XY3[1] - XY1[1])**2)  # Length of element 1
        L2 = np.sqrt((XY3[0] - XY2[0])**2 + (XY3[1] - XY2[1])**2)  # Length of element 2
        lanX1, lanY1 = (XY3[0] - XY1[0]) / L1, (XY3[1] - XY1[1]) / L1  # Direction cosines for element 1
        lanX2, lanY2 = (XY3[0] - XY2[0]) / L2, (XY3[1] - XY2[1]) / L2  # Direction cosines for element 2
        
        # Element strain
        es1 = (L1 - L1i) / L1i # Strain of element 1
        es2 = (L2 - L2i) / L2i # Strain of element 2
        
        # Stress-strain relationship for element 1
        if 0 <= es1 <= ey or -ey <= es1 < 0:
            fs1 = Es * es1
            fstan1 = Es
        elif ey < abs(es1) <= esu:
            fs1 = fy + Esh * (abs(es1) - ey)
            fstan1 = (fy + Esh * (abs(es1) - ey)) / abs(es1)
        else:
            fs1, fstan1 = 0, 0

        # Stress-strain relationship for element 2
        if 0 <= es2 <= ey or -ey <= es2 < 0:
            fs2 = Es * es2
            fstan2 = Es
        elif ey < abs(es2) <= esu:
            fs2 = fy + Esh * (abs(es2) - ey)
            fstan2 = (fy + Esh * (abs(es2) - ey)) / abs(es2)
        else:
            fs2, fstan2 = 0, 0

        # Element stiffness
        G1 = fstan1 * A1 / L1 # Stiffness of element 1
        G2 = fstan2 * A2 / L2 # Stiffness of element 2
        
        # Updated stiffness    
        K = np.array([G1 * lanY1 ** 2 + G2 * lanY2 ** 2]) 
        
        # Force imbalance
        f = K * u - F 
        # Calculate du and update residual
        du = -f / Kini # Displacement increment
        residual = np.max(np.abs(du)) # Update residual
        it += 1 # Increment iteration counter

        if it == MAX_ITERATIONS:
            print(f"             {i+1}                      {it}")
            print("    ## Trial iterations reached to Ultimate .The solution for this step is not converged ##")
            break

        u += du # Update displacement

    if it <= MAX_ITERATIONS:
        print(f"             {i+1}                      {it}")


    # Store results
    lanX1L1.append(lanX1)  # Store direction cosine for element 1
    lanY1L1.append(lanY1)  # Store direction cosine for element 1
    lanX2L1.append(lanX2)  # Store direction cosine for element 2
    lanY2L1.append(lanY2)  # Store direction cosine for element 2
    DU2.append(residual)   # Store residuals
    I2.append(i)           # Store step index
    IT2.append(it)         # Store iteration count
    XY3ox2.append(XY3[0])  # Store X-coordinate of point 3
    XY3oy2.append(XY3[1])  # Store Y-coordinate of point 3
    es1 = (L1 - L1i) / L1i  # Strain in element 1
    es2 = (L2 - L2i) / L2i  # Strain in element 2
    ess1o2.append(es1)      # Store strain for element 1
    ess2o2.append(es2)      # Store strain for element 2
    LL1.append(L1)          # Store length of element 1
    LL2.append(L2)          # Store length of element 2

    TBS5N.append(-1 * (fstan1 * A1 * es1 * lanX1 + fstan2 * A2 * es2 * lanX2))  # Total base reaction

    U1L.append(up)  
    U2L.append(u[0])  # Store displacement
    U5N.append(up)    # Store displacement increment


    # Convergence Conditions
    if abs(es1) >= esu and abs(es2) >= esu:
        print('      ## Strain in element(1) and element(2) reached to Ultimate Strain ##')
        break
    if abs(es1) >= esu:
        print('      ## Strain in element(1) reached to Ultimate Strain ##')
        break
    if abs(es2) >= esu:
        print('      ## Strain in element(2) reached to Ultimate Strain ##')
        break
    if abs(up) >= D5max:
        print('      ## Displacement at [DOF (5)] reached to Ultimate Displacement ##')
        break


# Convert results to arrays
U5N = np.array(U5N)     # X-displacement [DOF 5]
TBS5N = np.array(TBS5N) # Forces in Base-Reaction

#-----------------------------------------------------


print("##################")
print("# OpensSees Code #")
print("##################")

import openseespy.opensees as ops
from Analysis_Function import ANALYSIS

XY1 = [0, 0]      # [x y] Point 1
XY2 = [0, 500]    # [x y] Point 2
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

# Materials (MATERIAL NONLINEARITY)
ops.uniaxialMaterial('Steel01', 1, fy, Es, b)  # Steel with bilinear kinematic hardening Material

# Elements (GEOMETRIC NONLINEARITY)
ops.element('corotTruss', 1, 1, 3, A1, 1)
ops.element('corotTruss', 2, 2, 3, A2, 1)

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

# Output Data
ops.recorder('Node', '-file', f"DTH_PUSH.txt",'-time', '-node', 3, '-dof', 1,2, 'disp')        # Displacement Time History Node 3
ops.recorder('Node', '-file', f"BTH_PUSH_01.txt",'-time', '-node', 1, '-dof', 1,2, 'reaction') # Base Shear Time History Node 1
ops.recorder('Node', '-file', f"BTH_PUSH_02.txt",'-time', '-node', 2, '-dof', 1,2, 'reaction') # Base Shear Time History Node 2

# Results storage
displacements = []
forces = []

# Perform analysis
displacement_X = 0
step = 0
while abs(displacement_X) < abs(D5max) and step < Nsteps:
    print(f'STEP {step+1}')
    stable = ops.analyze(1)
    ANALYSIS(stable, 1, TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS
    # Record results
    displacement_X = ops.nodeDisp(3, 1) # Displacement in X Direction
    displacement_Y = ops.nodeDisp(3, 2) # Displacement in Y Direction
    displacements.append(displacement_X)
    step += 1
    """
    # Update geometry
    Xnew = XY3[0] + displacement_X
    Ynew = XY3[1] + displacement_Y
    #print('Geometry:', Xnew)
    L1 = np.sqrt((Xnew - XY1[0])**2 + (Ynew - XY1[1])**2)
    L2 = np.sqrt((Xnew - XY2[0])**2 + (Ynew - XY2[1])**2)
    lanX1, lanY1 = (Xnew - XY1[0]) / L1, (Ynew - XY1[1]) / L1
    lanX2, lanY2 = (Xnew - XY2[0]) / L2, (Ynew - XY2[1]) / L2
    
    # Reaction Forces from Elements: 
    # You are also calculating reaction forces using element strains in the following section of your code
    new_x3 = L1 * lanX1
    new_y3 = L1 * lanY1 
    es1 = (L1 - L1i) / L1i
    es2 = (L2 - L2i) / L2i
    reaction_force = (-EA1 * es1 * lanX1) + (-EA2 * es2 * lanX2)
    forces.append(reaction_force)
    """

#displacements = np.array(displacements)
#forces = np.array(forces)

totaltime = ti.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n')

#-----------------------------------------------------
# Post-processing
#-----------------------------------------------------
def OUTPUT_SECOND_COLUMN(X, COLUMN):
    import numpy as np
    # Time History
    filename = f"{X}.txt"
    data_collected = np.loadtxt(filename)
    X = data_collected[:, COLUMN]   
    return X 
    
dispP = OUTPUT_SECOND_COLUMN('DTH_PUSH', 1) # Reading Disp from Text file - NODE 1
base01 = OUTPUT_SECOND_COLUMN('BTH_PUSH_01', 1) # Reading base reaction from Text file - NODE 1
base02 = OUTPUT_SECOND_COLUMN('BTH_PUSH_02', 1) # Reading base reaction from Text file - NODE 2
forces = base01 + base02

#-----------------------------------------------------
### Plot the Results from the Output Data:
import matplotlib.pyplot as plt

# Plot Y-displacement (D5) vs Base Reaction (kN)
plt.figure(figsize=(10, 5))
plt.plot(U5N, TBS5N, label='Python', color='black')
plt.plot(dispP, forces, label='OpenSees', color='r', linestyle='--')
plt.xlabel('X-Displacement (mm) [DOF 5]')
plt.ylabel('Base Reaction (kN) [DOF 1] + [DOF 3]')
plt.title('X-Displacement vs X-Base Reaction')
plt.grid(True)
plt.legend()
plt.show()
