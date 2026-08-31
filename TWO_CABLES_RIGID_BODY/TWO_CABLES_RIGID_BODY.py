###################################################################
#                   >> IN THE NAME OF ALLAH <<                    #
# PUSHOVER ANALYSIS OF CABLE SUBJECTED TO LATERAL DISPLACEMENT,   #
# MATERIAL AND GEOMETRIC NONLINEARITY EFFECTS, AND DISPLACEMENT   #
#           CONTROL VERIFICATION USING OPENSEES                   #
#-----------------------------------------------------------------#
#  THIS PROGRAM IS WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)  #
#              E-MAIL: SALAR.D.GHASHGHAEI@GMAIL.COM               #
###################################################################
"""
This Python code performs a pushover analysis of a cable subjected to lateral displacement
 using both custom calculations and the OpenSees library.
 
Key steps include:
1. Parameter Initialization: Sets geometry, material properties, and simulation parameters.
2. Custom Calculation: Uses Newton-Raphson iterations to calculate displacement and force responses for truss members based on incremental displacement.
3. OpenSees Modeling: Creates a finite element model, applies boundary conditions, loads, and performs displacement-controlled analysis.
4. Post-processing: Reads displacements and forces from results files.
5. Plotting Results: Compares Python and OpenSees results to validate the analysis.
"""

import numpy as np
import time as ti
import openseespy.opensees as ops
import ANALYSIS_FUNCTION as S01
import PLOT_2D_TRUSS as S02
import matplotlib.pyplot as plt

# Define parameters (units: mm, kN)
P5 = 0                       # [kN] External Load [DOF (5)]
P6 = 0                       # [kN] External Load [DOF (6)]
D6 = -0.01                   # [mm] Initial Displacement [DOF (6)] Incremental Displacement
D6max = 600.0                # [mm] Maximum displacement [DOF (6)]
XY1i = [0.0, 0.0]            # [x y] Point 1
XY2i = [500.0, 0.0]          # [x y] Point 2
XY3i = [1000.0, 0.0]         # [x y] Point 3

XY4i = [500.0, 500.0]        # [x y] Point 4
XY5i = [1000.0, 500.0]       # [x y] Point 5

A1 = np.pi * (50**2) / 4     # [mm^2] Cross Section Area for Element 1
A2 = np.pi * (50**2) / 4     # [mm^2] Cross Section Area for Element 2

Nsteps =  int(np.abs(D6max / D6))# Number of steps for calculation

# Steel Section Properties
fy = 0.240        # [kN/mm^2] Yield strength of steel section
fu = 1.1818 * fy  # [kN/mm^2] Ultimate strength of steel section
Es = 200.0        # [kN/mm^2] Modulus of elasticity of steel section
ey = fy / Es      # [mm/mm] Yield steel strain
esu = 0.35        # [mm/mm] Ultimate steel strain

Esh = (fu - fy) / (esu - ey)  # Strain hardening modulus
b = Esh / Es                  # Strain hardening ratio


MAX_ITERATIONS = 10000     # convergence iteration for test
MAX_TOLERANCE = 1.0e-6     # convergence tolerance for test

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
U6N = []      # imposed Y displacement at node C
TBS6N = []    # total base reaction (vertical)
XY3ox2 = []   # X coordinate of node C
XY3oy2 = []   # Y coordinate of node C
ess1o2 = []   # strain cable 1
ess2o2 = []   # strain cable 2
LL1 = []      # length cable 1
LL2 = []      # length cable 2
lanX1L1, lanY1L1, lanX2L1, lanY2L1 = [], [], [], []
DU2, I2, IT2 = [], [], []
U1L, U2L = [], []
Null2 = []

# Geometry (mm)
A = np.array([0.0, 0.0])          # Node 1, pinned
B0 = np.array([500.0, 0.0])       # Node 2 initial
C0 = np.array([1000.0, 0.0])      # Node 3 initial
D = np.array([500.0, 500.0])      # Node 4, fixed
E = np.array([1000.0, 500.0])     # Node 5, fixed

L_total = 1000.0                  # length of rigid bar (A to C)
L1_0 = 500.0                      # initial length cable 1 (B-D)
L2_0 = 500.0                      # initial length cable 2 (C-E)

# Material parameters (same as OpenSees HystereticSM envelope)
FY = 0.240
FU = 1.1818 * FY
Ke = 200.0
DY = FY / Ke
DSU = 0.36
Ksh = (FU - FY) / (DSU - DY)
b = Ksh / Ke

# Positive envelope points (strain, stress)
pos_strains = np.array([0.0, DY, DSU, 1.1*DSU, 1.25*DSU])
pos_stresses = np.array([0.0, FY, FU, 0.2*FU, 0.1*FU])
# Negative envelope points
neg_strains = np.array([-DY, -DSU, -1.1*DSU, -1.25*DSU])
neg_stresses = np.array([-FY, -FU, -0.2*FU, -0.1*FU])

# Combine for interpolation (include zero)
all_strains = np.concatenate(([neg_strains[-1]], neg_strains, [0.0], pos_strains))
all_stresses = np.concatenate(([neg_stresses[-1]], neg_stresses, [0.0], pos_stresses))
# Sort by strain
idx = np.argsort(all_strains)
all_strains = all_strains[idx]
all_stresses = all_stresses[idx]

def material_response(eps):
    """Return stress for given strain using piecewise linear interpolation."""
    return np.interp(eps, all_strains, all_stresses)

# Gradually increase the applied displacement (vertical at node C)
for i in range(Nsteps):
    up = D6 * i  # vertical displacement at node C (negative downward)
    
    # Compute rotation angle of rigid bar
    sin_theta = up / L_total
    theta = np.arcsin(sin_theta)   # negative for downward
    cos_theta = np.cos(theta)
    
    # New positions of B and C due to rigid rotation about A
    B_new = np.array([L1_0 * cos_theta, L1_0 * sin_theta])  # since AB=500
    C_new = np.array([L_total * cos_theta, L_total * sin_theta])
    
    # Cable lengths
    L1 = np.linalg.norm(B_new - D)
    L2 = np.linalg.norm(C_new - E)
    
    # Strains
    es1 = (L1 - L1_0) / L1_0
    es2 = (L2 - L2_0) / L2_0
    
    # Stresses from material model
    fs1 = material_response(es1)
    fs2 = material_response(es2)
    
    # Axial forces (tension positive)
    N1 = fs1 * A1
    N2 = fs2 * A2
    
    # Total vertical reaction at supports (D and E)
    # Vertical component of force on support from cable
    # For cable 1: force on D is N1 * (B_new - D)/L1, vertical component = N1 * (B_new[1]-D[1])/L1
    # Reaction at D is opposite sign
    R_y1 = -N1 * (B_new[1] - D[1]) / L1
    R_y2 = -N2 * (C_new[1] - E[1]) / L2
    total_reaction = R_y1 + R_y2
    
    # Store results
    U6N.append(up)
    TBS6N.append(total_reaction)
    XY3ox2.append(C_new[0])
    XY3oy2.append(C_new[1])
    ess1o2.append(es1)
    ess2o2.append(es2)
    LL1.append(L1)
    LL2.append(L2)
    
    # Direction cosines (for compatibility with later code, not strictly needed)
    lanX1 = (B_new[0] - D[0]) / L1
    lanY1 = (B_new[1] - D[1]) / L1
    lanX2 = (C_new[0] - E[0]) / L2
    lanY2 = (C_new[1] - E[1]) / L2
    lanX1L1.append(lanX1)
    lanY1L1.append(lanY1)
    lanX2L1.append(lanX2)
    lanY2L1.append(lanY2)
    
    U1L.append(up)                  # same as U6N
    U2L.append(C_new[0] - XY5i[0])  # horizontal displacement at C
    DU2.append(0.0)                 # no residual (direct solution)
    I2.append(i)
    IT2.append(1)                   # one iteration (direct)
    """
    # Convergence conditions (ultimate strain)
    if abs(es1) >= DSU and abs(es2) >= DSU:
        print('      ## Strain in element(1) and element(2) reached to Ultimate Strain ##')
        break
    if abs(es1) >= DSU:
        print('      ## Strain in element(1) reached to Ultimate Strain ##')
        break
    if abs(es2) >= DSU:
        print('      ## Strain in element(2) reached to Ultimate Strain ##')
        break
    if abs(up) >= D6max:
        print('      ## Displacement at [DOF (6)] reached to Ultimate Displacement ##')
        break
    """
    print(f"             {i+1}                      1")

# Convert results to arrays
U5N = np.array(U6N)     # Y-displacement [DOF 6]
TBSYN = np.array(TBS6N) # Forces in Y-Base-Reaction


#%% -----------------------------------------------------


print("##################")
print("# OpensSees Code #")
print("##################")



XY1 = [0.0, 0.0]      # [x y] Point 1
XY2 = [500.0, 0.0]    # [x y] Point 2
XY3 = [1000.0, 0.0]   # [x y] Point 3

XY4 = [500.0, 500.0]  # [x y] Point 4
XY5 = [1000.0, 500.0] # [x y] Point 5

FY = 0.240                                     # [kN] Yield Strength of Element
FU = 1.1818 * FY                               # [kN] Ultimate Srength of Element
Ke = 200.0                                     # [kN/mm] Elastic Module
DY = FY / Ke                                   # [mm] Yield Strain
DSU = 0.36                                     # [mm] Ultimate Srain
Ksh = (FU - FY) / (DSU - DY)                   # [kN/mm] Strain Hardening Modulus
Kp = FU / DSU                                  # [kN/mm] Plastic Module
b = Ksh / Ke                                   # Strain Hardening Ratio

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
#%% Model setup
ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 3)

#%% Nodes
ops.node(1, *XY1)
ops.node(2, *XY2)
ops.node(3, *XY3)
ops.node(4, *XY4)
ops.node(5, *XY5)

#%% Boundary Conditions
ops.fix(1, 1, 1, 0)
ops.fix(4, 1, 1, 1)
ops.fix(5, 1, 1, 1)

#%% Materials (MATERIAL NONLINEARITY)
MatTag = 1
#ops.uniaxialMaterial('Steel01', MatTag, fy, Es, b)  # Steel with bilinear kinematic hardening Material
ops.uniaxialMaterial('HystereticSM', MatTag, '-posEnv', *KP.flatten(), '-negEnv', *KN.flatten(), '-pinch', 1, 1)
#%% Elements (GEOMETRIC NONLINEARITY)
ops.element('corotTruss', 1, 2, 4, A1, MatTag) # CABLE 01
ops.element('corotTruss', 2, 3, 5, A2, MatTag) # CABLE 02

BEAM_TYPE = 'NON-RIGID'

if BEAM_TYPE == 'RIGID':
    ops.rigidLink('-beam', 1, 2) # RIGID ELEMENT 01
    ops.rigidLink('-beam', 2, 3) # RIGID ELEMENT 02
    
if BEAM_TYPE == 'NON-RIGID':
    Ec = 200.0                   # [kN/mm^2] Young's Modulus
    Bc = 50.0                    # [mm] Beam Section Width
    Hc = 50.0                    # [mm] Beam Section Height
    AREAc = Bc*Hc
    IZc = Bc*(Hc**3) / 12
    # Geometric transformation
    Tratag = 1
    ops.geomTransf('Linear', Tratag)
    #ops.geomTransf('Corotational', Tratag)
    ops.element('elasticBeamColumn', 3, 1, 2, AREAc, Ec, IZc, Tratag) # BEAM ELEMENT 01 
    ops.element('elasticBeamColumn', 4, 2, 3, AREAc, Ec, IZc, Tratag) # BEAM ELEMENT 02
    

#%% Time series and pattern
ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)

#%% Apply external force 
ops.load(3, P5, P6, 0.0)

print('Model Built')

NstepGravity = 10
DGravity = 1 / NstepGravity
ops.integrator('LoadControl', DGravity) # determine the next time step for an analysis
ops.numberer('Plain')                   # renumber dof's to minimize band-width (optimization), if you want to
ops.system('BandGeneral')               # how to store and solve the system of equations in the analysis
ops.constraints('Plain')                # how it handles boundary conditions
ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS) # determine if convergence has been achieved at the end of an iteration step
ops.algorithm('Newton')                 # use Newton's solution algorithm: updates tangent stiffness at every iteration
ops.analysis('Static')                  # define type of analysis static or transient
ops.analyze(NstepGravity)               # apply gravity

ops.loadConst('-time', 0.0)             # maintain constant gravity loads and reset time to zero

ops.timeSeries('Linear', 2)
ops.pattern('Plain', 200, 2)
ops.load(2, 0.0, -1.0, 0.0)
ops.load(3, 0.0, -1.0, 0.0)
    
#%% Analysis setup
ops.wipeAnalysis()
ops.system('BandGeneral')
ops.numberer('Plain')
ops.constraints('Plain')
#ops.integrator('DisplacementControl', 2, 2, 0.5*D6) # Apply displacement control at DOF 4 (node 2, Y)

ops.integrator('DisplacementControl', 3, 2, D6) # Apply displacement control at DOF 6 (node 3, Y)
ops.algorithm('Newton')
ops.analysis('Static')

#%% Output Data
ops.recorder('Node', '-file', "DTH_PUSH.txt",'-time', '-node', 3, '-dof', 1, 2, 'disp')        # Displacement Time History Node 3
ops.recorder('Node', '-file', "BTH_PUSH_01.txt",'-time', '-node', 4, '-dof', 1, 2, 'reaction') # Base Reaction Time History Node 4
ops.recorder('Node', '-file', "BTH_PUSH_02.txt",'-time', '-node', 5, '-dof', 1, 2, 'reaction') # Base Reaction Time History Node 5

#%% Results storage
displacement_X, displacement_Y = [], []
force1 = []; force2 = []; reaction = [];

#%% Perform analysis

step = 0
while step < Nsteps:
    print(f'STEP {step+1}')
    stable = ops.analyze(1)
    S01.ANALYSIS(stable, 1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS
    # Record results
    displacement_X.append(ops.nodeDisp(3, 1)) # Displacement in X Direction
    displacement_Y.append(ops.nodeDisp(3, 2)) # Displacement in Y Direction
    f1 = ops.eleResponse(1, 'force')[4]       # INTERNAL FORCE FOR ELEMENT 01
    f2 = ops.eleResponse(2, 'force')[4]       # INTERNAL FORCE FOR ELEMENT 02
    force1.append(f1)
    force2.append(f2)
    ops.reactions()
    reaction.append(ops.nodeReaction(4, 2) + ops.nodeReaction(5, 2)) # Total Reaction Force
    #reaction.append(np.abs((f1+f2)))    # Total Reaction Force
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

# %% Plot 2D Frame Shapes for Nonlinear Static Analysis
S02.PLOT_2D_FRAME(deformed_scale=1)  # Adjust scale factor as needed

#-------------------
# Post-processing
#-------------------

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


# Plot Y-displacement (D6) vs Eleemnts Force (kN)
plt.figure(0, figsize=(10, 5))
plt.plot(displacement_Y, force1, label='Ele. Force 1', color='black', linestyle='--', linewidth=4)
plt.plot(displacement_Y, force2, label='Ele. Force 2', color='purple', linewidth=4)
plt.xlabel('Y-Displacement (mm) [DOF 6]')
plt.ylabel('Elements Force')
plt.title('Y-Displacement vs Elements Force')
plt.grid(True)
plt.legend()
plt.show()

# Plot Y-displacement (D6) vs Base Reaction (kN)
plt.figure(1, figsize=(10, 5))
plt.plot(U5N, TBSYN, label='Python', color='black', linewidth=4)
plt.plot(displacement_Y, reaction, label='OpenSees', color='r', linestyle='--', linewidth=4)
plt.xlabel('Y-Displacement (mm) [DOF 6]')
plt.ylabel('Base Reaction (kN) [DOF 8] + [DOF 10]')
plt.title('Y-Displacement vs Y-Base Reaction')
plt.grid(True)
plt.legend()
plt.show()
