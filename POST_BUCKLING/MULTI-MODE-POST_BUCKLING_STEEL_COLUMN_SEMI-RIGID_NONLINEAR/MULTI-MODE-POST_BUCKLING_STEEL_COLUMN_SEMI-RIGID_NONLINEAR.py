###########################################################################################################
#                     >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                    #
#          INVESTIGATION OF MULTI-MODE POST-BUCKLING PHENOMENA IN SEMI-RIGID STEEL COLUMNS                #
#           USING OPENSEES CONSIDERING THE GEOMETRIC AND MATERIAL PROPERTIES NONLINEARITY                 #
#---------------------------------------------------------------------------------------------------------#
# IT MODELS A 2D INELASTIC BEAM-COLUMN WITH AN INITIAL IMPERFECTION (FOUR DIFFRENET SHAPES)               #
#  AND APPLIES AN AXIAL COMPRESSIVE LOAD TO ANALYZE LARGE DISPLACEMENTS.                                  #
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
import STEEL_FIBER_SECTION as S01
import ANALYSIS_FUNCTION as S02
import OPENSEEES_HYSTERETIC_FORCE_DISP_FUN as S03
import time as TI
#%%------------------------------------------------------------------------------------------------
# Define parameters (units: mm, N)
# Define model parameters
P = 1          # [N] Axial Load
M = 1          # [N.mm] Moment Load
L = 3000       # [mm] Column length
N = 100        # Number of elements
dy = L / N     # Element length
ε = 0.001 * L  # [mm] Imperfection amplitude
#MODE_NUM = 2   # Mode Shape Number
#%%------------------------------------------------------------------------------------------------
total_steps = 9200
rotat_inc = 0.00001  # Rotation increment (rad)
disp_inc = -0.01     # Displacement increment (mm, negative for compression)
#%%------------------------------------------------------------------------------------------------
# Define  Steel Rebar Material Properties (Steel01)
secTag = 1
fy = 240.0                    # [N/mm²] Yield strength of steel section
fu = 1.5 * fy                 # [N/mm²] Ultimate strength of steel section
Es = 200000.0                 # [N/mm²] Modulus of elasticity of steel section
ey = fy / Es                  # [mm/mm] Yield steel strain
esu = 0.35                    # [mm/mm] Ultimate steel strain
Esh = (fu - fy) / (esu - ey)  # [N/mm^2] Strain hardening modulus
b = Esh / Es                  # Strain hardening ratio
#%%------------------------------------------------------------------------------------------------
DENSITY_STEEL = 7850/1e9      # [kg/m^3] -> [kg/mm^3] Steel Material Density
#%%------------------------------------------------------------------------------------------------
# Define Analysis Properties
MAX_ITERATIONS = 5000         # Convergence iteration for test
MAX_TOLERANCE = 1.0e-8        # Convergence tolerance for test
#%%------------------------------------------------------------------------------------------------
# Define Post-Buckling Funtion
def POST_BUCKLING(MODE_NUM):
    # Create model (2D, 3 DOF per node: u_x, u_y, theta)
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    
    # Define nodes with initial imperfection ((1/MODE_NUM) sine wave)
    for i in range(N+1):
        y = i * dy
        x = ε * np.sin(MODE_NUM * np.pi * y / L)  # Imperfection in x-direction
        ops.node(i+1, x, y)
    
    #------------------------------------------------------------------------------------------------
    # Define Semi-rigid connection with rotational spring
    
    ops.node(N+2, 0, 0) # Rotational spring in the bottom of column   
    ops.fix(N+2, 1, 1, 1)
    ops.node(N+3, 0, L) # Rotational spring in the top of column   
    ops.fix(N+3, 1, 1, 1)
    #------------------------------------------------------------------------------------------------
    # Define boundary conditions for simply supported column
    ops.fix(1, 1, 1, 0)    # Bottom node: fix u_x, u_y, free theta
    ops.fix(N+1, 1, 0, 0)  # Top node: fix u_x, free u_y, free theta
    #------------------------------------------------------------------------------------------------
    # Define Nonlinear Connection in Bottom and Top Rotational Springs
    matTag_S = 500
    MY = 12430.0           # [N.mm] Yield Moment of springs
    TY = 0.001474          # [Rad] Yield Rotation of springs
    MU = 1.36 * MY         # [N.mm] Utimate Moment of springs
    TSU = 0.0541           # [Rad] Utimate Rotation of springs
    #ROTATIONAl_STIFFNESS = 1e14 # [N/mm] Rotational stiffness of springs
    pinchX = 0.4           # Pinching factor in X direction
    pinchY = 0.2           # Pinching factor in Y direction
    damage1 = 0.0          # Damage due to ductility
    damage2 = 0.0          # Damage due to energy
    beta = 0.1             # Stiffness degradation parameter
    # MOMENT-ROTATION RELATIONSHIP OF ROTATIONAL SPRING AND PLOT 
    DP = [0, 0, 0]
    FP = [0, 0, 0]
    DN = [0, 0, 0]
    FN = [0, 0, 0]
    DP[0], FP[0] = TY, MY
    DP[1], FP[1] = TSU, MU 
    DP[2], FP[2] = 1.13*TSU, 0.43*MU
    DN[0], FN[0] = -TY, -MY 
    DN[1], FN[1] = -TSU, -MU
    DN[2], FN[2] = -1.07*TSU, -0.43*MU   
    S03.OPENSEEES_HYSTERETIC_FORCE_DISP_FUN(matTag_S, DP, FP, DN, FN, 
                                        pinchX, pinchY,
                                        damage1, damage2, beta,
                                        PLOT = True, X_LABEL='Rotation (rad)', Y_LABEL='Moment [N.mm]', TITLE='MOMENT-ROTATION CURVE')
    
    #ops.uniaxialMaterial('Elastic', matTag_S, ROTATIONAl_STIFFNESS)
    ops.element('zeroLength', N+200, N+2, 1, '-mat', matTag_S, '-dir', 3)   # ROTATIONAL SPRING AT BOTTOM OF COLUMN
    ops.element('zeroLength', N+300, N+1, N+3, '-mat', matTag_S, '-dir', 3) # ROTATIONAL SPRING AT TOP OF COLUMN 
    #------------------------------------------------------------------------------------------------
    # Materials (STEEL MATERIAL NONLINEARITY)
    """
    matTag = 1
    ops.uniaxialMaterial('Steel01', matTag, fy, Es, b) # Steel with bilinear kinematic hardening Material
    """   
    matTag = 1
    pinchX = 0.8           # Pinching factor in X direction
    pinchY = 0.5           # Pinching factor in Y direction
    damage1 = 0.0          # Damage due to ductility
    damage2 = 0.0          # Damage due to energy
    beta = 0.1             # Stiffness degradation parameter
    ops.uniaxialMaterial('Hysteretic', matTag, fy, ey, fu, esu, 0.2*fu, 1.1*esu, -fy, -ey, -fu, -esu, -0.2*fu, -1.1*esu, pinchX, pinchY, damage1, damage2, beta)
    # INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Hysteretic_Material
    #------------------------------------------------------------------------------------------------
    # Define a plain load pattern (for static analysis)
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(N+1, 0.0, P, M)
    #------------------------------------------------------------------------------------------------
    # Create the Double I section section
    matTag = 1
    #Depth, Ele_Mass = S01.DOUBLE_I_SECTION_FIBER(secTag, matTag, PLOT=True, DENSITY=DENSITY_STEEL)
    Depth, Ele_Mass = S01.DOUBLE_I_SECTION_QUAD(secTag, matTag, PLOT=True, DENSITY=DENSITY_STEEL)
    AREA = Ele_Mass / DENSITY_STEEL
    PY = AREA * fy
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
    
    # Rotational Increments
    ops.integrator('DisplacementControl', node_tag, 3, rotat_inc)
    # Lateral Increments
    ops.integrator('DisplacementControl', node_tag, 2, disp_inc) # y-direction (axial)
    ops.analysis('Static')  # Static analysis
    #------------------------------------------------------------------------------------------------
    # Run analysis and record results
    mid_node = (N//2) + 1  # Mid-height node for lateral displacement
    axial_disp = []        # Axial displacement at top node
    rotation = []          # Rotation of column
    lateral_disp = []      # Lateral displacement at mid-height
    axial_force_BOT = []        # BOTTOM NODE - Axial compressive load (reaction force)
    shear_force_BOT = []        # BOTTOM NODE - Shear load (reaction force)
    moment_force_BOT = []       # BOTTOM NODE - Moment load (reaction force)
    axial_force_TOP = []        # TOP NODE - Axial compressive load (reaction force)
    shear_force_TOP = []        # TOP NODE - Shear load (reaction force)
    moment_force_TOP = []       # TOP NODE - Moment load (reaction force)
    ELE_D_FORCE, SEC_D_FORCE, SEC_D_STIFFNESS, ELE_D_FLEX = [], [], [], []
    #total_steps = int(np.abs(1/disp_inc))  # Number of analysis steps
    
    for step in range(total_steps):
        OK = ops.analyze(1)  # Perform one analysis step
        S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS)  # CHECK THE ANALYSIS
        axial_disp.append(ops.nodeDisp(node_tag, 2))         # Record axial displacement
        rotation.append(ops.nodeDisp(node_tag, 3))           # Record rotation
        lateral_disp.append(ops.nodeDisp(mid_node, 1))       # Record lateral displacement
        axial_force_BOT.append(-ops.eleResponse(1, 'force')[1])   # Record reaction force at bottom node (Axial compressive load)
        shear_force_BOT.append(-ops.eleResponse(1, 'force')[0])   # Record reaction force at bottom node (shear load)
        moment_force_BOT.append(-ops.eleResponse(1, 'force')[2])  # Record reaction force at bottom node (moment load)
        # IN EACH STEP, OUTPUT MIDDLE ELEMENT AND SECTION DATA
        # ELEMENT FORCE
        ELE_D_FORCE.append(ops.eleForce(mid_node-1)) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/eleForce.html
        #print(ELE_D_FORCE[-1])
        # SECTION FORCE
        SEC_D_FORCE.append(ops.sectionForce(mid_node-1, 1)) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/sectionForce.html
        #print(SEC_D_FORCE[-1])
        # SECTION SIFFNESS
        SEC_D_STIFFNESS.append(ops.sectionStiffness(mid_node-1, 1)) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/sectionStiff.html
        #print(SEC_D_STIFFNESS[-1])
        # SECTION FLEXIBILIY
        ELE_D_FLEX.append(ops.sectionFlexibility(mid_node-1, 1)) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/sectionFlexibility.html
        #print(ELE_D_FLEX[-1])
        print(f'STEP {step+1} DONE')
    else:
        print('Analysis completed successfully')
        
    DATA = (axial_disp, rotation, lateral_disp, axial_force_BOT,
            shear_force_BOT, moment_force_BOT, PY, ELE_D_FORCE,
            SEC_D_FORCE, SEC_D_STIFFNESS, ELE_D_FLEX)
    return DATA
        
#%%------------------------------------------------------------------------------------------------
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
    

#%%------------------------------------------------------------------------------------------------
# Analysis Durations:
current_time = TI.strftime("%H:%M:%S", TI.localtime())
print("Start Time:", current_time)

# MODE 1
MODE_NUM = 1
DATA = POST_BUCKLING(MODE_NUM)
(axial_disp_01, rotation_01, lateral_disp_01,
 axial_load_01, shear_load_01, moment_load_01, PY,
 ELE_D_FORCE_01, SEC_D_FORCE_01, SEC_D_STIFFNESS_01, ELE_D_FLEX_01) = DATA 

# Plot frame shapes
PLOT_FRAME(deformed_scale=1)  # Adjust scale factor as needed

# MODE 2
MODE_NUM = 2
DATA = POST_BUCKLING(MODE_NUM)
(axial_disp_02, rotation_02, lateral_disp_02,
 axial_load_02, shear_load_02, moment_load_02, PY,
 ELE_D_FORCE_02, SEC_D_FORCE_02, SEC_D_STIFFNESS_02, ELE_D_FLEX_02) = DATA 

# Plot frame shapes
PLOT_FRAME(deformed_scale=1)  # Adjust scale factor as needed

# MODE 3
MODE_NUM = 3
DATA = POST_BUCKLING(MODE_NUM)
(axial_disp_03, rotation_03, lateral_disp_03,
 axial_load_03, shear_load_03, moment_load_03, PY,
 ELE_D_FORCE_03, SEC_D_FORCE_03, SEC_D_STIFFNESS_03, ELE_D_FLEX_03) = DATA 

# Plot frame shapes
PLOT_FRAME(deformed_scale=1)  # Adjust scale factor as needed

# MODE 4
MODE_NUM = 4
DATA = POST_BUCKLING(MODE_NUM)
(axial_disp_04, rotation_04, lateral_disp_04,
 axial_load_04, shear_load_04, moment_load_04, PY,
 ELE_D_FORCE_04, SEC_D_FORCE_04, SEC_D_STIFFNESS_04, ELE_D_FLEX_04) = DATA 

# Plot frame shapes
PLOT_FRAME(deformed_scale=1)  # Adjust scale factor as needed
    
current_time = TI.strftime("%H:%M:%S", TI.localtime())
print("Finish Time:", current_time)
#%%------------------------------------------------------------------------------------------------    
# Plot axial load vs. axial displacement
plt.figure(1, figsize=(8, 6))
plt.plot(axial_disp_01, np.abs(axial_load_01)/PY, label=f"Mode 01 - Post Buckling Load: {np.abs(axial_load_01[-1]) : .3f}", linewidth=4)
plt.plot(axial_disp_02, np.abs(axial_load_02)/PY, label=f"Mode 02 - Post Buckling Load: {np.abs(axial_load_02[-1]) : .3f}", linewidth=4)
plt.plot(axial_disp_03, np.abs(axial_load_03)/PY, label=f"Mode 03 - Post Buckling Load: {np.abs(axial_load_03[-1]) : .3f}", linewidth=4)
plt.plot(axial_disp_04, np.abs(axial_load_04)/PY, label=f"Mode 04 - Post Buckling Load: {np.abs(axial_load_04[-1]) : .3f}", linewidth=4)
# ADD CRITICAL BUCKLING LOAD AS A HORIZONTAL LINE
plt.xlabel('Axial displacement at top-height (mm)')
#plt.ylabel('Axial compressive load (N)')
plt.ylabel('Axial compressive load / Axial Yield Load Capacity')
plt.title('Post-buckling behavior of column')
plt.grid(True)
#plt.semilogy()
plt.legend()
plt.show()    
#%%------------------------------------------------------------------------------------------------
# Plot axial load vs. lateral displacement
plt.figure(2, figsize=(8, 6))
plt.plot(lateral_disp_01, np.abs(axial_load_01)/PY, linewidth=4, label=f"Mode 01 - Lateral Disp: {np.abs(lateral_disp_01[-1]) : .3f}")
plt.plot(lateral_disp_02, np.abs(axial_load_02)/PY, linewidth=4, label=f"Mode 02 - Lateral Disp: {np.abs(lateral_disp_02[-1]) : .3f}")
plt.plot(lateral_disp_03, np.abs(axial_load_03)/PY, linewidth=4, label=f"Mode 03 - Lateral Disp: {np.abs(lateral_disp_03[-1]) : .3f}")
plt.plot(lateral_disp_04, np.abs(axial_load_04)/PY, linewidth=4, label=f"Mode 04 - Lateral Disp: {np.abs(lateral_disp_04[-1]) : .3f}")
plt.xlabel('Lateral displacement at mid-height (mm)')
#plt.ylabel('Axial compressive load (N)')
plt.ylabel('Axial compressive load / Axial Yield Load Capacity')
plt.title('Post-buckling behavior of column')
plt.grid(True)
plt.legend()
plt.show()
#%%------------------------------------------------------------------------------------------------
# Plot moment vs, rotation
plt.figure(3, figsize=(8, 6))
plt.plot(rotation_01, np.abs(moment_load_01), linewidth=4, label=f"Mode 01 - Moment: {np.abs(moment_load_01[-1]) : .3f}")
plt.plot(rotation_02, np.abs(moment_load_02), linewidth=4, label=f"Mode 02 - Moment: {np.abs(moment_load_02[-1]) : .3f}")
plt.plot(rotation_03, np.abs(moment_load_03), linewidth=4, label=f"Mode 03 - Moment: {np.abs(moment_load_03[-1]) : .3f}")
plt.plot(rotation_04, np.abs(moment_load_04), linewidth=4, label=f"Mode 04 - Moment: {np.abs(moment_load_04[-1]) : .3f}")
plt.xlabel('Rotation at top-height (rad)')
plt.ylabel('Moment (N.mm)')
plt.title('Post-buckling behavior of column')
plt.grid(True)
plt.legend()
plt.show()
#%%------------------------------------------------------------------------------------------------
# EXCEL OUTPUT
import pandas as pd

# Create DataFrame function
def create_df(axial, rotation, lateral, axial_load, shear_load, moment_load):
    df = pd.DataFrame({
        "Axial_Displacement": axial,
        "Rotation": rotation,
        "Lateral_Displacement": lateral,
        "Axial_Load": axial_load,
        "Shear_Load": shear_load,
        "Moment_Load": moment_load
    })
    return df


# Save to Excel
with pd.ExcelWriter("MULTI-MODE-POST_BUCKLING_STEEL_COLUMN_SEMI-RIGID_NONLINEAR_OUTPUT.xlsx", engine='openpyxl') as writer:
    
    # MODE 1
    df1 = create_df(axial_disp_01, rotation_01, lateral_disp_01,
                    axial_load_01, shear_load_01, moment_load_01)
    df1.to_excel(writer, sheet_name="Mode_1", index=False)
    
    # MODE 2
    df2 = create_df(axial_disp_02, rotation_02, lateral_disp_02,
                    axial_load_02, shear_load_02, moment_load_02)
    df2.to_excel(writer, sheet_name="Mode_2", index=False)

    # MODE 3
    df3 = create_df(axial_disp_03, rotation_03, lateral_disp_03,
                    axial_load_03, shear_load_03, moment_load_03)
    df3.to_excel(writer, sheet_name="Mode_3", index=False)

    # MODE 4
    df4 = create_df(axial_disp_04, rotation_04, lateral_disp_04,
                    axial_load_04, shear_load_04, moment_load_04)
    df4.to_excel(writer, sheet_name="Mode_4", index=False)
#%%------------------------------------------------------------------------------------------------