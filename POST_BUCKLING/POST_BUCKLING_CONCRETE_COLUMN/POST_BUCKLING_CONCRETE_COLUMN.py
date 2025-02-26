###########################################################################################################
#                                         IN THE NAME OF ALLAH                                            #
#  SIMULATION OF THE POST-BUCKLING BEHAVIOR OF A SIMPLY SUPPORTED CONFINED CONCRETE COLUMN USING OPENSEES #
#                        CONSIDERING THE GEOMETRIC AND MATERIAL PROPERTIES NONLINEARITY                   #
#---------------------------------------------------------------------------------------------------------#
# IT MODELS A 2D ELASTIC BEAM-COLUMN WITH AN INITIAL IMPERFECTION AND APPLIES AN AXIAL COMPRESSIVE        #
# LOAD TO ANALYZE LARGE DISPLACEMENTS.                                                                    #
# 1. MODEL SETUP: A COLUMN OF LENGTH L IS DEFINED WITH N ELEMENTS, INCORPORATING A SMALL INITIAL          #
# IMPERFECTION (HALF-SINE WAVE).                                                                          #
# 2. NODES & BOUNDARY CONDITIONS: NODES ARE CREATED, WITH THE BOTTOM FIXED IN X, Y AND THE TOP FIXED IN   #
# X BUT FREE IN Y AND ROTATION.                                                                           #
# 3. ELEMENT DEFINITION: THE COLUMN IS MODELED USING ELASTIC BEAM-COLUMN ELEMENTS WITH COROTATIONAL       #
# TRANSFORMATION FOR GEOMETRIC NONLINEARITY.                                                              #
# 4. LOAD APPLICATION: A STATIC AXIAL FORCE P IS APPLIED AT THE TOP NODE.                                 #
# 5. ANALYSIS SETUP: A DISPLACEMENTCONTROL INTEGRATOR IS USED TO INCREMENTALLY PUSH THE COLUMN DOWNWARDS. #  
# 6. NONLINEAR SOLVER: THE NEWTON METHOD IS USED WITH A NORMDISPINCR TEST FOR CONVERGENCE.                #
# 7. ANALYSIS EXECUTION: THE LOOP PERFORMS INCREMENTAL LOADING STEPS, RECORDING AXIAL DISPLACEMENT,       #
# LATERAL DISPLACEMENT, AND AXIAL FORCE.                                                                  #
# 8. BUCKLING BEHAVIOR CAPTURE: LATERAL DISPLACEMENTS AT THE MID-HEIGHT NODE INDICATE POST-BUCKLING       #
# DEFORMATION.                                                                                            #
# 9. RESULTS EXTRACTION: REACTION FORCES AT THE BASE NODE PROVIDE THE AXIAL COMPRESSIVE LOAD.             #
# 10. PLOTTING: THE SCRIPT VISUALIZES AXIAL FORCE VS. LATERAL DISPLACEMENT, SHOWING THE POST-BUCKLING     #
# RESPONSE OF THE COLUMN.                                                                                 #
#---------------------------------------------------------------------------------------------------------#
#                          THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                     #
#                                   EMAIL: salar.d.ghashghaei@gmail.com                                   #
###########################################################################################################
import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
import CONCRETE_FIBER_SECTION as S01
import ANALYSIS_FUNCTION as S02

# Define parameters (units: m, N)
# Define model parameters
P = 400        # [N] Axial Load
L = 3000       # [mm] Column length
N = 100        # Number of elements
dy = L / N     # Element length
ε = 0.001 * L  # [m] Imperfection amplitude

#--------------------------------------------------------------------
# Define  Steel Rebar Material Properties (Steel01)
fy = 4000              # [N/mm²] Yield strength of steel rebar
fu = 1.5 * fy          # [N/mm²] Ultimate strength of steel rebar
Es = 200000            # [N/mm²] Modulus of elasticity of steel rebar
ey = fy / Es           # [mm/mm] Yield steel strain
esu = 0.35             # [mm/mm] Ultimate steel strain
Esh = (fu - fy) / (esu - ey)  # [N/mm^2] Strain hardening modulus
b = Esh / Es                  # Strain hardening ratio
#--------------------------------------------------------------------
# Define material properties for (Concrete02)
fcp = -30              # [N/mm²] Compressive strength
epsc0 = -0.0025        # [mm/mm] Strain at peak compressive strength
E0 = fcp / epsc0
fpcu = fcp * 0.05      # [N/mm²] Crushing strength
epsU = -0.004          # [mm/mm] Strain at ultimate stress
lamda = 0.1            # Ratio between unloading slope at epsU and initial slope
ft = 5                 # [N/mm²] Tensile strength
Ets = ft/0.002         # [N/mm²] Tension softening stiffness 
#--------------------------------------------------------------------
# Define a rectangular concrete section using fibers.
secTag = 1
HSec = 500        # [mm] Section depth
BSec = 400        # [mm] Section width
coverH = 50       # [mm] Cover depth
coverB = 50       # [mm] Cover width
coreID = 2        # Core concrete material ID
coverID = 2       # Cover concrete material ID
steelID = 1       # Steel material ID
numBarsTop = 3    # Number of top bars
numBarsBot = 3    # Number of bottom bars
numBarsIntTot = 8 # Total number of intermediate bars (must be even)
barAreaTop = 200  # [mm^2] Area of top bars
barAreaBot = 150  # [mm^2] Area of intermediate bars
barAreaInt = 100  # [mm^2] Area of intermediate bars
#------------------------------------------------------------------------------------------------
DENSITY_CONCRETE = 2500/1e9      # [kg/mm^3] Concrete Density
#------------------------------------------------------------------------------------------------
# Define Analysis Properties
MAX_ITERATIONS = 1000      # Convergence iteration for test
MAX_TOLERANCE = 1.0e-10    # Convergence tolerance for test
#------------------------------------------------------------------------------------------------
# Create model (2D, 3 DOF per node: u_x, u_y, theta)
ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 3)

# Define nodes with initial imperfection (half-sine wave)
for i in range(N+1):
    y = i * dy
    x = ε * np.sin(np.pi * y / L)  # Imperfection in x-direction
    ops.node(i+1, x, y)
#------------------------------------------------------------------------------------------------
# Define boundary conditions for simply supported column
ops.fix(1, 1, 1, 0)    # Bottom node: fix u_x, u_y, free theta
ops.fix(N+1, 1, 0, 0)  # Top node: fix u_x, free u_y, free theta
#------------------------------------------------------------------------------------------------
# Materials (STEEL MATERIAL NONLINEARITY)
matTag01 = 1
ops.uniaxialMaterial('Steel01', matTag01, fy, Es, b) # Steel with bilinear kinematic hardening Material
#------------------------------------------------------------------------------------------------
# Materials (CONCRETE MATERIAL NONLINEARITY)
matTag02 = 2
ops.uniaxialMaterial('Concrete02', matTag02, fcp, epsc0, fpcu, epsU, lamda, ft, Ets)
#ops.uniaxialMaterial('Concrete01', matTag02, fcp, epsc0, fpcu, epsU)

# Create the RC rectangular section
Depth, Ele_Mass = S01.CONCRETE_CONFINED_REC_SECTION_QUAD(secTag, HSec, BSec, coverH, coverB, coreID, coverID, steelID, numBarsTop, barAreaTop, numBarsBot, barAreaBot, numBarsIntTot, barAreaInt, PLOT=True, CONCRETE_DENSITY=2500/1e9, STEEL_DENSITY=7850/1e9)
#------------------------------------------------------------------------------------------------
# Define a plain load pattern (for static analysis)
ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)
ops.load(N+1, 0.0, P, 0.0)
#------------------------------------------------------------------------------------------------
# Define geometric transformation (corotational for large displacements)
transfTag = 1
ops.geomTransf('Corotational', transfTag)
numIntgrPts = 5
# Define elastic beam elements
for i in range(1, N+1):
    ops.element('nonlinearBeamColumn', i, i, i+1, numIntgrPts, secTag, transfTag,'-mass', Ele_Mass)
#------------------------------------------------------------------------------------------------
# Analysis settings
ops.constraints('Plain')
ops.numberer('Plain')
ops.system('BandGeneral')  # System of equations
ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS, 2)  # Convergence test
ops.algorithm('Newton')  # Nonlinear solver
node_tag = N+1  # Top node for displacement control
dof = 2  # y-direction (axial)
disp_inc = -0.01  # Displacement increment (m, negative for compression)
ops.integrator('DisplacementControl', node_tag, dof, disp_inc)
ops.analysis('Static')  # Static analysis
#------------------------------------------------------------------------------------------------
# Run analysis and record results
mid_node = (N//2) + 1  # Mid-height node for lateral displacement
axial_disp = []  # Axial displacement at top node
lateral_disp = []  # Lateral displacement at mid-height
axial_load = []  # Axial compressive load (reaction force)
#total_steps = int(np.abs(1/disp_inc))  # Number of analysis steps
total_steps = 690

for step in range(total_steps):
    OK = ops.analyze(1)  # Perform one analysis step
    S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS)  # CHECK THE ANALYSIS
    axial_disp.append(ops.nodeDisp(node_tag, dof))      # Record axial displacement
    lateral_disp.append(ops.nodeDisp(mid_node, 1))      # Record lateral displacement
    axial_load.append(-ops.eleResponse(1, 'force')[1])  # Record reaction force at bottom node (compressive load)
    print(f'STEP {step+1} DONE')
else:
    print('Analysis completed successfully')
#------------------------------------------------------------------------------------------------    
# Plot axial load vs. axial displacement
I = (BSec * HSec**3)/12
Ec = 4700 * np.sqrt(abs(fcp))
PCR = (np.pi**2 * Ec * I)/ L**2
plt.figure(figsize=(8, 6))
plt.plot(axial_disp, np.abs(axial_load), 'black', label=f"Post Buckling Load: {np.abs(axial_load[-1]) : .3f}", linewidth=4)
# ADD CRITICAL BUCKLING LOAD AS A HORIZONTAL LINE
#plt.axhline(y=PCR, color='r', linestyle='--', label=f"Critical Load (Pcr): {PCR : .3f}", linewidth=4) 
plt.xlabel('Axial displacement at top-height (mm)')
plt.ylabel('Axial compressive load (N)')
plt.title('Post-buckling behavior of simply supported column')
plt.grid(True)
plt.legend()
plt.show()    
#------------------------------------------------------------------------------------------------
# Plot axial load vs. lateral displacement
plt.figure(figsize=(8, 6))
plt.plot(lateral_disp, np.abs(axial_load), 'black', linewidth=4)
plt.xlabel('Lateral displacement at mid-height (mm)')
plt.ylabel('Axial compressive load (N)')
plt.title('Post-buckling behavior of simply supported column')
plt.grid(True)
plt.show()
#------------------------------------------------------------------------------------------------
# Define a function to plot the frame shapes
def PLOT_FRAME(deformed_scale=1.0):
    fig, ax = plt.subplots(1, figsize=(20, 16))

    # Extract node coordinates
    nodes = ops.getNodeTags()
    node_coords = {node: ops.nodeCoord(node) for node in nodes}

    # Plot undeformed shape
    for ele in ops.getEleTags():
        node1, node2 = ops.eleNodes(ele)
        x1, y1 = node_coords[node1]
        x2, y2 = node_coords[node2]
        ax.plot([x1, x2], [y1, y2], 'k-', label='Undeformed' if ele == 1 else "")  # Black line for undeformed

    # Plot deformed shape
    for ele in ops.getEleTags():
        node1, node2 = ops.eleNodes(ele)
        x1, y1 = node_coords[node1]
        x2, y2 = node_coords[node2]

        ux1, uy1, _ = ops.nodeDisp(node1)  # Displacement at node1
        ux2, uy2, _ = ops.nodeDisp(node2)  # Displacement at node2

        ax.plot([x1 + deformed_scale * ux1, x2 + deformed_scale * ux2],
                [y1 + deformed_scale * uy1, y2 + deformed_scale * uy2],
                'r--', label='Deformed' if ele == 1 else "")  # Red dashed line for deformed
                
    # Annotate nodes with their tags
    for node, (x, y) in node_coords.items():
        ux, uy, _ = ops.nodeDisp(node)  # Displacement at node
        ax.text(x, y, f"{node}", color='blue', fontsize=12, ha='center', label='Node Tags' if node == 1 else "")  # Undeformed
        ax.text(x + deformed_scale * ux, y + deformed_scale * uy, f"{node}", color='purple', fontsize=12, ha='center')  # Deformed            

    #ax.set_aspect('equal', 'box')
    ax.set_xlabel('X [mm]')
    ax.set_ylabel('Y [mm]')
    ax.set_title(f'Undeformed and Deformed Shapes - SCALE: {deformed_scale:.3f}')
    ax.legend()
    ax.grid()
    plt.show()
    
# Plot frame shapes
PLOT_FRAME(deformed_scale=100)  # Adjust scale factor as needed

#------------------------------------------------------------------------------------------------
