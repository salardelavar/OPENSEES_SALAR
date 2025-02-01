###########################################################################################################
#                                         IN THE NAME OF ALLAH                                            #
#                THERMAL LOADING STRUCTURAL ANALYSIS OF A 2D CONCRETE FRAME USING OPENSEES                #
#                          WITH UNCERTAINTY USING MONTE CARLO SIMULATION:                                 #
#                 INCORPORATING BETA PROBABILITY DISTRIBUTION FOR STOCHASTIC PARAMETERS                   #
#---------------------------------------------------------------------------------------------------------#
# This program uses OpenSees to model a 2D concrete frame subjected to thermal and distributed loads      #
# with Monte Carlo simulations to incorporate uncertainty using Beta probability distributions for        #
# material properties and geometrical parameters.                                                         #
#---------------------------------------------------------------------------------------------------------#
#                   THIS PROGRAM IS WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                         #
#                                EMAIL: SALAR.D.GHASHGHAEI@GMAIL.COM                                      #
###########################################################################################################

""" 
 Key Features:
 [1] 2D Frame Model Definition: 
 The structure has a specified number of stories and bays, with fixed supports at the base.
 Each story and bay's height and width are modeled probabilistically with a Beta distribution to simulate uncertainty.
 
 [2] Material Properties with Thermal Effects: 
 Steel and concrete material properties are defined with Beta distribution to account for variability
 in yield strength, modulus of elasticity, and other characteristics. The concrete model includes thermal effects,
 especially for temperature-dependent properties.
 
 [3] Element Setup: 
 The frame elements are modeled using beam-column elements with corotational transformations and Lobatto beam
 integration for accurate bending response. Distributed and thermal loads are applied to the beams, with thermal
 gradients modeled specifically in the first story.
 
 [4] Thermal Loading: 
 A linear thermal gradient is applied to simulate the effect of temperature variation within the structural frame, with
 the maximum thermal load defined using the Beta distribution. These effects are incorporated into the analysis using
 the specialized material models 'Steel01Thermal' and 'Concrete02Thermal'.
 
 [5] Analysis Procedure: 
 A static analysis is performed using load control with thermal increments. The Newton-Raphson method is used to solve the system,
 and convergence is monitored with specified tolerances and iterations.
 
 [6] Monte Carlo Simulation: 
 Multiple simulations (NUM_SIM) are run with varying material properties and geometry, using Monte Carlo methods to study
 the probabilistic behavior of the frame under thermal loading.
 
 [7] Data Collection & Visualization: 
 The program records displacements, reactions, and stress-strain data for each simulation. These are output into text files for
 post-processing, visualization, and analysis of the structural response.
"""

import time as TI
import numpy as np
import matplotlib.pyplot as plt
import openseespy.opensees as ops
from Analysis_Function import ANALYSIS
import SALAR_MATH as S01
import MARKOV_CHAIN as S03
from CONCRETE_FIBERTHERMAL_SECTION import R_RECTANGULAR_CONCRETE_SECTION_REBAR_B, R_RECTANGULAR_CONCRETE_SECTION_REBAR_C

#--------------------------------------------------------------------
# Define parameters (units: mm, kN)
num_stories = 3        # Number of Stories
num_bays = 4           # Number of Bays

NUM_SIM = 6000                                          # NUMBER OF SIMULATIONS
story_height = S01.BETA_PDF(2950, 3150, 2, 1, NUM_SIM)  # [mm] Height of each story
bay_width = S01.BETA_PDF(6950, 7150, 2, 1, NUM_SIM)     # [mm] Width of each bay

S01.HISROGRAM_BOXPLOT(story_height, HISTO_COLOR='pink', LABEL='Height of each story [mm]')
S01.HISROGRAM_BOXPLOT(bay_width, HISTO_COLOR='lightblue', LABEL='Width of each bay [mm]')
#--------------------------------------------------------------------
# Define  Steel Rebar Material Properties (Steel01Thermal)
fy = S01.BETA_PDF(0.38, 0.41, 1, 2, NUM_SIM)               # [kN/mm²] Yield strength of steel rebar
fu = 1.5 * fy                                              # [kN/mm²] Ultimate strength of steel rebar
Es = S01.BETA_PDF(190, 210, 1, 2, NUM_SIM)                 # [kN/mm²] Modulus of elasticity of steel rebar
ey = fy / Es                                               # [mm/mm] Yield steel strain
esu = S01.BETA_PDF(0.32, 0.36, 1, 2, NUM_SIM)              # [mm/mm] Ultimate steel strain
Esh = (fu - fy) / (esu - ey)                               # [kN/mm^2] Strain hardening modulus
b = Esh / Es                                               # Strain hardening ratio

S01.HISROGRAM_BOXPLOT(fy, HISTO_COLOR='cyan', LABEL='Yield strength of steel rebar [kN/mm²]')
S01.HISROGRAM_BOXPLOT(fu, HISTO_COLOR='orange', LABEL='Ultimate strength of steel rebar [kN/mm²]')
S01.HISROGRAM_BOXPLOT(Es, HISTO_COLOR='lightgreen', LABEL='Modulus of elasticity of steel rebar [kN/mm²]')
S01.HISROGRAM_BOXPLOT(esu, HISTO_COLOR='grey', LABEL='Ultimate steel strain [mm/mm]')
#--------------------------------------------------------------------
# Define material properties for (Concrete02Thermal)
fcp = S01.BETA_PDF(-0.032, -0.028, 1, 2, NUM_SIM)            # [kN/mm²] Compressive strength
epsc0 = S01.BETA_PDF(-0.0026, -0.0022, 1, 2, NUM_SIM)        # [mm/mm] Strain at peak compressive strength
E0 = fcp / epsc0
fpcu = fcp * 0.05                                            # [kN/mm²] Crushing strength
epsU = S01.BETA_PDF(-0.0043, -0.0035, 1, 2, NUM_SIM)         # [mm/mm] Strain at ultimate stress
lamda = S01.BETA_PDF(0.09, 0.11, 1, 2, NUM_SIM)              # Ratio between unloading slope at epsU and initial slope
ft = S01.BETA_PDF(0.004, 0.005, 1, 2, NUM_SIM)               # [kN/mm²] Tensile strength
Ets = ft/S01.BETA_PDF(0.0018, 0.0022, 1, 2, NUM_SIM)         # [kN/mm²] Tension softening stiffness 

S01.HISROGRAM_BOXPLOT(fcp, HISTO_COLOR='lime', LABEL='Compressive strength [kN/mm²]')
S01.HISROGRAM_BOXPLOT(epsc0, HISTO_COLOR='yellow', LABEL='Strain at peak compressive strength [mm/mm]')
S01.HISROGRAM_BOXPLOT(epsU, HISTO_COLOR='brown', LABEL='Strain at ultimate stress [mm/mm]')
#--------------------------------------------------------------------
# Define Thermal and Distributed load
Max_Thermal = S01.BETA_PDF(700, 800, 1, 1, NUM_SIM)               # [°C] Temperature
distributed_load = S01.BETA_PDF(-0.0030, -0.0025, 1, 1, NUM_SIM)  # [kN/mm] Uniform Distributed Loads

S01.HISROGRAM_BOXPLOT(Max_Thermal, HISTO_COLOR='gold', LABEL='Temperature [°C]')
S01.HISROGRAM_BOXPLOT(distributed_load, HISTO_COLOR='lightgreen', LABEL='Distributed load [kN/mm]')
#--------------------------------------------------------------------
# Define Analysis Properties
MAX_ITERATIONS = 1000      # Convergence iteration for test
MAX_TOLERANCE = 1.0e-10    # Convergence tolerance for test

Nstep = 200           # Number of incremental steps
Incr_Temp = 1/Nstep   # Incremental temperature step
#--------------------------------------------------------------------
# Define a rectangular concrete section using fibers.
B = S01.BETA_PDF(395, 405, 1, 2, NUM_SIM)       # [mm] Width of the rectangular section
H = S01.BETA_PDF(495, 505, 1, 2, NUM_SIM)       # [mm] Height of the rectangular section
COVER = S01.BETA_PDF(20, 40, 1, 2, NUM_SIM)     # [mm] Cover of the rectangular section
RD = S01.BETA_PDF(22, 25, 1, 2, NUM_SIM)        # [mm] Diameter of the rebar
NUM_B = 3     # Number of fibers along the width of the section
NUM_H = 100   # Number of fibers along the height of the section
NUM_LAYERS = 2 # Number of Rebar fiber layers

S01.HISROGRAM_BOXPLOT(B, HISTO_COLOR='brown', LABEL='Width of the rectangular section [mm]')
S01.HISROGRAM_BOXPLOT(H, HISTO_COLOR='lime', LABEL='Height of the rectangular section [mm]')
S01.HISROGRAM_BOXPLOT(COVER, HISTO_COLOR='blue', LABEL='Cover of the rectangular section [mm]')
S01.HISROGRAM_BOXPLOT(RD, HISTO_COLOR='purple', LABEL='Diameter of the rebar [mm]')
#--------------------------------------------------------------------
# Thermal Loading Analysis
def TEMP_ANAL(I):
    # Clear pervious model
    ops.wipe()
    # Define model
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    # Node coordinates
    node_id = 1
    for i in range(num_stories + 1):  # Including ground level
        for j in range(num_bays + 1):  # Including leftmost column
            ops.node(node_id, j * bay_width[I], i * story_height[I])
            if i == 0:  # Fix supports at ground level
                ops.fix(node_id, 1, 1, 1)
            else:  # Free to move for higher levels
                ops.fix(node_id, 0, 0, 0)
            node_id += 1

    #print(' Structure Coordinates Done.') 
    #--------------------------------------------------------------------
    # Materials (STEEL MATERIAL NONLINEARITY)
    matTag01 = 1
    ops.uniaxialMaterial('Steel01Thermal', matTag01, fy[I], Es[I], b[I]) # Steel with bilinear kinematic hardening Material with thermaal effect
    #--------------------------------------------------------------------
    # Materials (CONCRETE MATERIAL NONLINEARITY)
    matTag02 = 2
    """
    Concrete02Thermal:
    Concrete02Thermal is created for modelling concrete, which is derived from
     the standard "Concrete02" material and incorprates with temperature dependent
     properties sugggested in Eurocode 2 (EN1992-1-2)
    """
    ops.uniaxialMaterial('Concrete02Thermal', matTag02, fcp[I], epsc0[I], fpcu[I], epsU[I], lamda[I], ft[I], Ets[I])
    #ops.uniaxialMaterial('Concrete01Thermal', matTag02, fcp[I], epsc0[I], fpcu[I], epsU[I])
    """
    ConcreteECThermal:
    ConcreteECThermal is derived by modification of the existing concrete material
     class "Concrete02" and "Concrete02Thermal" to include the temperature-dependent
     properties according to EN 1992-1-2 about concrete with siliceous aggregates at elevated temperature.
    """
    #ops.uniaxialMaterial('ConcreteECThermal', matTag02, fcp, epsc0[I], fpcu[I], epsU[I], lamda, ft[I], Ets[I]) # Concrete hardening Material with thermal effect
    #--------------------------------------------------------------------
    # BEAMS SECTION
    secTag01 = 1    # Section beam tag identifier
    # Concrete Fiber Section without Rebars
    #Depth = R_RECTANGULAR_CONCRETE_SECTION(secTag01, B, H, NUM_B, NUM_H, matTag01)
    # Concrete Fiber Section with Rebars 
    NUM_LAYERS = 2  # Number of rebar layers
    Depth01 = R_RECTANGULAR_CONCRETE_SECTION_REBAR_B(secTag01, B[I], H[I], COVER[I], RD[I], NUM_B, NUM_H, NUM_LAYERS, matTag02, matTag01, PLOT=False)
    # COLUMNS SECTION
    secTag02 = 2    # Section column tag identifier
    NUM_LAYERS = 4  # Number of rebar layers
    NUM_FIRST = 4   # Number of rebars in first and last layer
    Depth02 = R_RECTANGULAR_CONCRETE_SECTION_REBAR_C(secTag02, B[I], B[I], COVER[I], RD[I], NUM_B, NUM_H, NUM_LAYERS, NUM_FIRST, matTag02, matTag01, PLOT=False)

    #print(' Thermal Section Done.')
    #--------------------------------------------------------------------
    # Define geometric transformation
    transfTag = 1
    ops.geomTransf('Corotational', transfTag)
    #--------------------------------------------------------------------
    # Define beam integration (Lobatto integration)
    numIntegrationPoints = 5
    biTag01 = 1
    ops.beamIntegration('Lobatto', biTag01, secTag01, numIntegrationPoints)

    # Define column integration (Lobatto integration)
    numIntegrationPoints = 5
    biTag02 = 2
    ops.beamIntegration('Lobatto', biTag02, secTag02, numIntegrationPoints)
    #--------------------------------------------------------------------
    # Define beam-column elements
    # Thermo-mechaical beam-column Elements
    # INFO LINK: https://openseesforfire.github.io/Subpages/Elecmds.html
    element_id = 1
    beam_id = [] # Initialize the beam ID list
    for i in range(num_stories):
        for j in range(num_bays + 1):  # Vertical elements (columns)
            node1 = i * (num_bays + 1) + j + 1
            node2 = (i + 1) * (num_bays + 1) + j + 1
            # element dispBeamColumnThermal $eleID $NodeTag0 $NodeTag1 $NumInts $secTag $TransfTag
            ops.element('dispBeamColumnThermal', element_id, node1, node2, transfTag, biTag02)
            # element forceBeamColumnThermal $eleID $NodeTag0 $NodeTag1 $NumInts $secTag $TransfTag
            #ops.element('forceBeamColumnThermal', element_id, node1, node2, secTag02, transfTag)
            element_id += 1

        for j in range(num_bays):  # Horizontal elements (beams)
            node1 = (i + 1) * (num_bays + 1) + j + 1
            node2 = (i + 1) * (num_bays + 1) + j + 2
            ops.element('dispBeamColumnThermal', element_id, node1, node2, transfTag, biTag01)
            if i == 0:  # Only store beams from the first story
                beam_id.append(element_id)
            element_id += 1

    #print(' Thermal Element Done.')        
    #--------------------------------------------------------------------
    # Define time series
    Tag = 1
    ops.timeSeries('Linear', Tag)
    # Define load pattern
    patternTag = 1
    ops.pattern('Plain', patternTag, Tag)
    #--------------------------------------------------------------------
    # Apply distributed load to all beams in the structure
    for ele_id in beam_id:
        # eleLoad -ele $eleTag1 <$eleTag2 ....> -type -beamUniform $Wz <$Wx>
        ops.eleLoad('-ele', ele_id, '-type', '-beamUniform', distributed_load[I])
    #--------------------------------------------------------------------    
    # A linear thermal gradient is applied to the elements in the first story.
    # The temperature varies from -DD to DD across the section height.
    # INFO LINK: https://openseesforfire.github.io/Subpages/ThermalActionCmds.html

    # Apply thermal load to all elements
    #for ele_id in range(1, element_id):
    #    ops.eleLoad('-ele', ele_id, '-type', '-beamThermal', Max_Thermal[I], -DD, Max_Thermal[I], DD)

    # Apply thermal load only to the identified beam elements
    DD = 0.5 * Depth01
    for ele_id in beam_id:
        if ele_id == 6 or ele_id == 7 or ele_id == 8 or ele_id == 9: # Apply Thermal Loads only to the Beams of Story 1
            # eleLoad -ele $eleTag -type -beamThermal $T1 $y1 $T2 $Y2
            ops.eleLoad('-ele', ele_id, '-type', '-beamThermal', Max_Thermal[I], -DD, Max_Thermal[I], DD)  

    #print(' Thermal Forces Done.')     
    #--------------------------------------------------------------------
    # Define analysis
    ops.system('BandGeneral')
    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/normDispIncr.html
    ops.algorithm('Newton')
    ops.integrator('LoadControl', Incr_Temp)
    ops.analysis('Static')

    #print(' Thermal Analysis Properties Done.')
    #--------------------------------------------------------------------
    # Output Data
    ops.recorder('Node', '-file', f"DTH_PUSH_{I}.txt",'-time', '-node', 6, '-dof', 1,2,3, 'disp')        # Displacement Time History Node 6
    ops.recorder('Node', '-file', f"BTH_PUSH_01_{I}.txt",'-time', '-node', 1, '-dof', 1,2,3, 'reaction') # Base Reaction Time History Node 1
    ops.recorder('Node', '-file', f"BTH_PUSH_02_{I}.txt",'-time', '-node', 2, '-dof', 1,2,3, 'reaction') # Base Reaction Time History Node 2
    ops.recorder('Node', '-file', f"BTH_PUSH_03_{I}.txt",'-time', '-node', 3, '-dof', 1,2,3, 'reaction') # Base Reaction Time History Node 3
    ops.recorder('Node', '-file', f"BTH_PUSH_04_{I}.txt",'-time', '-node', 4, '-dof', 1,2,3, 'reaction') # Base Reaction Time History Node 4
    ops.recorder('Node', '-file', f"BTH_PUSH_05_{I}.txt",'-time', '-node', 5, '-dof', 1,2,3, 'reaction') # Base Reaction Time History Node 5
    #ops.recorder('Element', '-file', f'STRESS_STRAIN_BOT_{I}.txt', '-time', '-ele', 6, 'section', secTag01, 'fiber', -DD, 0, 'stressStrainTangent')
    #ops.recorder('Element', '-file', f'STRESS_STRAIN_TOP_{I}.txt', '-time', '-ele', 6, 'section', secTag01, 'fiber', +DD, 0, 'stressStrainTangent')
    ops.recorder('Element', '-file', f'STRESS_STRAIN_BOT_{I}.txt', '-time', '-ele', 6, 'section', secTag01, 'fiber', -DD, 0, 'stressStrain')
    ops.recorder('Element', '-file', f'STRESS_STRAIN_TOP_{I}.txt', '-time', '-ele', 6, 'section', secTag01, 'fiber', +DD, 0, 'stressStrain')
    #--------------------------------------------------------------------
    # Perform analysis
    temp = []
    dispX = []
    dispY = []

    mid_node = (num_stories * (num_bays + 1) + num_bays // 2 + 1)  # Middle node at top

    for i in range(Nstep):
        OK = ops.analyze(1)
        #print('STEP: ', i + 1)
        ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS

        temp.append(ops.getLoadFactor(patternTag) * Max_Thermal)
        dispX.append(ops.nodeDisp(mid_node, 1)) # X Displacement
        dispY.append(ops.nodeDisp(mid_node, 2)) # Y Displacement
        
    #ops.wipe()   
    
    return temp, dispX, dispY, mid_node

#--------------------------------------------------------------------
# Analysis Durations:
starttime = TI.process_time()

# Initialize lists to store max values
max_displacement_X = []
max_displacement_Y = []
max_temp = []
max_base_reaction = []
# NUM_SIM is the number of simulations
for i in range(NUM_SIM):
    temp, dispX, dispY, mid_node = TEMP_ANAL(i)
    # Calculate and store the max absolute values
    max_displacement_X.append(np.max(np.abs(dispX)))
    max_displacement_Y.append(np.max(np.abs(dispY)))
    max_temp.append(np.max(np.abs(temp)))
    print(f'STEP {i + 1} DONE')    

totaltime = TI.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n') 
#--------------------------------------------------------------------

# Define a function to plot the frame shapes
def PLOT_FRAME(deformed_scale=1.0):
    fig, ax = plt.subplots(figsize=(20, 16))

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
    ax.set_title('Undeformed and Deformed Shapes')
    ax.legend()
    ax.grid()
    plt.show()
    
# Plot frame shapes
PLOT_FRAME(deformed_scale=10)  # Adjust scale factor as needed

#--------------------------------------------------------------------

def OUTPUT_SECOND_COLUMN(X, COLUMN):
    import numpy as np
    # Time History
    filename = f"{X}.txt"
    data_collected = np.loadtxt(filename)
    X = data_collected[:, COLUMN]   
    return X 
    
DISP_X, DISP_Y, DISP_Z = [], [], []
BASES_AXIAL, BASES_SHEAR, BASES_MOMENT = [], [], []
strain_B, strain_T, stress_B, stress_T = [], [], [], []

# Analysis Durations:
starttime = TI.process_time()

for i in range (NUM_SIM):   
    DISP_X.append(np.max(np.abs(OUTPUT_SECOND_COLUMN(f'DTH_PUSH_{i}', 1)))) # Reading Disp from Text file - X Direaction
    DISP_Y.append(np.max(np.abs(OUTPUT_SECOND_COLUMN(f'DTH_PUSH_{i}', 2)))) # Reading Disp from Text file - Y Direaction
    DISP_Z.append(np.max(np.abs(OUTPUT_SECOND_COLUMN(f'DTH_PUSH_{i}', 3)))) # Reading Disp from Text file - Z Direaction (Rotation)
    # AXIAL
    base01_Y = np.max(np.abs(OUTPUT_SECOND_COLUMN(f'BTH_PUSH_01_{i}', 1))) # Reading base reaction from Text file - NODE 1
    base02_Y = np.max(np.abs(OUTPUT_SECOND_COLUMN(f'BTH_PUSH_02_{i}', 1))) # Reading base reaction from Text file - NODE 2
    base03_Y = np.max(np.abs(OUTPUT_SECOND_COLUMN(f'BTH_PUSH_03_{i}', 1))) # Reading base reaction from Text file - NODE 3
    base04_Y = np.max(np.abs(OUTPUT_SECOND_COLUMN(f'BTH_PUSH_04_{i}', 1))) # Reading base reaction from Text file - NODE 4
    base05_Y = np.max(np.abs(OUTPUT_SECOND_COLUMN(f'BTH_PUSH_05_{i}', 1))) # Reading base reaction from Text file - NODE 5
    BASES_AXIAL.append(base01_Y + base02_Y + base03_Y + base04_Y + base05_Y)
    # SHEAR
    base01_X = np.max(np.abs(OUTPUT_SECOND_COLUMN(f'BTH_PUSH_01_{i}', 2))) # Reading base reaction from Text file - NODE 1
    base02_X = np.max(np.abs(OUTPUT_SECOND_COLUMN(f'BTH_PUSH_02_{i}', 2))) # Reading base reaction from Text file - NODE 2
    base03_X = np.max(np.abs(OUTPUT_SECOND_COLUMN(f'BTH_PUSH_03_{i}', 2))) # Reading base reaction from Text file - NODE 3
    base04_X = np.max(np.abs(OUTPUT_SECOND_COLUMN(f'BTH_PUSH_04_{i}', 2))) # Reading base reaction from Text file - NODE 4
    base05_X = np.max(np.abs(OUTPUT_SECOND_COLUMN(f'BTH_PUSH_05_{i}', 2))) # Reading base reaction from Text file - NODE 5
    BASES_SHEAR.append(base01_X + base02_X + base03_X + base04_X + base05_X)
    # MOMENT
    base01_Z = np.max(np.abs(OUTPUT_SECOND_COLUMN(f'BTH_PUSH_01_{i}', 3))) # Reading base reaction from Text file - NODE 1
    base02_Z = np.max(np.abs(OUTPUT_SECOND_COLUMN(f'BTH_PUSH_02_{i}', 3))) # Reading base reaction from Text file - NODE 2
    base03_Z = np.max(np.abs(OUTPUT_SECOND_COLUMN(f'BTH_PUSH_03_{i}', 3))) # Reading base reaction from Text file - NODE 3
    base04_Z = np.max(np.abs(OUTPUT_SECOND_COLUMN(f'BTH_PUSH_04_{i}', 3))) # Reading base reaction from Text file - NODE 4
    base05_Z = np.max(np.abs(OUTPUT_SECOND_COLUMN(f'BTH_PUSH_05_{i}', 3))) # Reading base reaction from Text file - NODE 5
    BASES_MOMENT.append(base01_Z + base02_Z + base03_Z + base04_Z + base05_Z)
    # STRESS AND STRAIN OF ELEMENT 6
    strain_B.append(np.max(np.abs(OUTPUT_SECOND_COLUMN(f'STRESS_STRAIN_BOT_{i}', 2)))) # Reading bottom strain from Text file
    stress_B.append(np.max(np.abs(OUTPUT_SECOND_COLUMN(f'STRESS_STRAIN_BOT_{i}', 1)))) # Reading bottom stress from Text file
    strain_T.append(np.max(np.abs(OUTPUT_SECOND_COLUMN(f'STRESS_STRAIN_TOP_{i}', 2)))) # Reading top strain from Text file
    stress_T.append(np.max(np.abs(OUTPUT_SECOND_COLUMN(f'STRESS_STRAIN_TOP_{i}', 1)))) # Reading top stress from Text file
    print(f'STEP {i + 1} DATA LOAD DONE')    

totaltime = TI.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n') 
#------------------------------------------------------------------------------------------------
S01.HISROGRAM_BOXPLOT(DISP_X, HISTO_COLOR='purple', LABEL='Displacement X [mm]')
S01.HISROGRAM_BOXPLOT(DISP_Y, HISTO_COLOR='blue', LABEL='Displacement Y [mm]')
S01.HISROGRAM_BOXPLOT(DISP_Z, HISTO_COLOR='green', LABEL='Rotation Z [rad]')

S01.HISROGRAM_BOXPLOT(BASES_AXIAL, HISTO_COLOR='pink', LABEL='BASE-AXIAL [kN]')
S01.HISROGRAM_BOXPLOT(BASES_SHEAR, HISTO_COLOR='cyan', LABEL='BASE-SHEAR [kN]')
S01.HISROGRAM_BOXPLOT(BASES_MOMENT, HISTO_COLOR='brown', LABEL='BASE-MOMENT [kN.mm]')
#------------------------------------------------------------------------------------------------
XLABEL = 'Displacement X [mm]'
YLABEL = 'Temperature [°C]'
TITLE = f'{YLABEL} and {XLABEL} scatter chart'
COLOR = 'orange'
X = DISP_X
Y = Max_Thermal
S01.PLOT_SCATTER(X, Y , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

# CLUSTER DATA
S01.CLUSTER_DATA(X, Y, XLABEL, YLABEL, MAX_CLUSTERS=9)
#------------------------------------------------------------------------------------------------
XLABEL = 'Displacement Y [mm]'
YLABEL = 'Temperature [°C]'
TITLE = f'{YLABEL} and {XLABEL} scatter chart'
COLOR = 'red'
X = DISP_Y
Y = Max_Thermal
S01.PLOT_SCATTER(X, Y , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

# CLUSTER DATA
S01.CLUSTER_DATA(X, Y, XLABEL, YLABEL, MAX_CLUSTERS=9)
#------------------------------------------------------------------------------------------------
XLABEL = 'Displacement X [mm]'
YLABEL = 'BASE-AXIAL [kN]'
TITLE = f'{YLABEL} and {XLABEL} scatter chart'
COLOR = 'purple'
X = DISP_X
Y = BASES_SHEAR
S01.PLOT_SCATTER(X, Y , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

# CLUSTER DATA
S01.CLUSTER_DATA(X, Y, XLABEL, YLABEL, MAX_CLUSTERS=9)
#------------------------------------------------------------------------------------------------
XLABEL = 'Displacement Y [mm]'
YLABEL = 'BASE-SHEAR [kN]'
TITLE = f'{YLABEL} and {XLABEL} scatter chart'
COLOR = 'pink'
X = DISP_Y
Y = BASES_AXIAL
S01.PLOT_SCATTER(X, Y , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

# CLUSTER DATA
S01.CLUSTER_DATA(X, Y, XLABEL, YLABEL, MAX_CLUSTERS=9)
#------------------------------------------------------------------------------------------------
XLABEL = 'Rotation Z [rad]'
YLABEL = 'BASE-MOMENT [kN.mm]'
TITLE = f'{YLABEL} and {XLABEL} scatter chart'
COLOR = 'cyan'
X = DISP_Z
Y = BASES_MOMENT
S01.PLOT_SCATTER(X, Y , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

# CLUSTER DATA
S01.CLUSTER_DATA(X, Y, XLABEL, YLABEL, MAX_CLUSTERS=9)
#------------------------------------------------------------------------------------------------
# PLOT THE TIME-HISTORY
#S01.PLOT_TIME_HISTORY(time, displacement, velocity, acceleration, base_reaction)
#------------------------------------------------------------------------------------------------

import pandas as pd

data = {
    'displacement_x': DISP_X,
    'displacement_y': DISP_Y,
    'rotation_z': DISP_Z,
    'base_axial': BASES_AXIAL,
    'base_shear': BASES_SHEAR,
    'base_moment': BASES_MOMENT,
}

# Convert to DataFrame
df = pd.DataFrame(data)
#------------------------------------------------------------------------------------------------
# PLOT HEATMAP FOR CORRELATION 
S01.PLOT_HEATMAP(df)
#------------------------------------------------------------------------------------------------
# RANDOM FOREST ANALYSIS
"""
This code predicts the seismic safety of a structure using simulation data by training a Random Forest Classifier to
 classify whether the system is "safe" or "unsafe" based on features like maximum displacement and base reaction. 
 
 A regression model is also trained to estimate safety likelihood. It evaluates model performance using
 metrics like classification accuracy, mean squared error, and R² score. Additionally, it identifies key features influencing
 safety through feature importance analysis. The tool aids in seismic risk assessment, structural optimization, and understanding
 critical safety parameters.
"""

#print(df)
threshold_displacement = 50  # [mm] target displacement threshold
S01.RANDOM_FOREST(df, threshold_displacement)
#------------------------------------------------------------------------------------------------
# MULTIPLE REGRESSION MODEL
S01.MULTIPLE_REGRESSION(df) 
#------------------------------------------------------------------------------------------------
# PERFORM RELIABILITY ANALYSIS FOR BASE REACTION AND ELEMENT CAPACITY
mean_capacity = 250       # Mean Base-Shear Capacity
std_dev_capacity = 20     # Std Base-Shear Capacity
num_sim = NUM_SIM
TITLE = 'base shear'
S01.RELIABILITY_ANALYSIS(BASES_SHEAR, num_sim, mean_capacity, std_dev_capacity, TITLE)
#------------------------------------------------------------------------------------------------
# MARKOV CHAIN MODEl (structural damage analysis by evaluating displacement)
FILE_TF = False         # Indicate whether to read data from a file or use provided data
file_path = None        # Not used when 'file_tf' is False
DATA = DISP_X # If not using a file, replace None with a NumPy array of data

S03.MARKOV_CHAIN(FILE_TF, file_path, DATA)
#------------------------------------------------------------------------------------------------
# MACHINE LEARNING: LONG SHORT-TREM MEMERY (LSTM) METHOD
x = DISP_X 
y = BASES_SHEAR 
look_back = int(NUM_SIM * 0.5)
ITERATION = 10
XLABEL = 'Displacement_X'
YLABEL = 'base shear'
S01.PREDICT_LSTM(x, y, look_back, ITERATION, XLABEL, YLABEL)
#------------------------------------------------------------------------------------------------
"""
# Plot results
plt.figure(1, figsize=(8, 6))
plt.plot(disp_Y, BASES_AXIAL, color='blue', linewidth=2)
plt.xlabel('Node 6 Displacement Y (mm)')
plt.ylabel('Base Axial Reaction (kN)')
plt.title('Base Axial Reaction vs Displacement Y')
plt.grid()
plt.show()

plt.figure(2, figsize=(8, 6))
plt.plot(disp_X, BASES_SHEAR, color='red', linewidth=2)
plt.xlabel('Node 6 Displacement X (mm)')
plt.ylabel('Base Shear Reaction (kN)')
plt.title('Base Shear Reaction vs Displacement X')
plt.grid()
plt.show()

plt.figure(3, figsize=(8, 6))
plt.plot(disp_Z, BASES_MOMENT, color='green', linewidth=2)
plt.xlabel('Node 6 Rotation (rad)')
plt.ylabel('Base Moment Reaction (kN.mm)')
plt.title('Base Moment Reaction vs Rotation')
plt.grid()
plt.show()

plt.figure(4, figsize=(8, 6))
plt.plot(temp, dispX, color='purple', linewidth=2)
plt.xlabel('Temperature (°C)')
plt.ylabel(f'Middle Node {mid_node} Displacement X (mm)')
plt.title('Temperature vs Displacement X')
plt.grid()
plt.show()

plt.figure(5, figsize=(8, 6))
plt.plot(temp, dispY, color='black', linewidth=2)
plt.xlabel('Temperature (°C)')
plt.ylabel(f'Middle Node {mid_node} Displacement Y (mm)')
plt.title('Temperature vs Displacement Y')
plt.grid()
plt.show()

plt.figure(6, figsize=(8, 6))
plt.plot(strain_B, stress_B, color='blue', label='Bottom Fiber', linewidth=2)
plt.plot(strain_T, stress_T, color='red', label='Top Fiber', linewidth=2)
plt.xlabel('Strain (mm/mm)')
plt.ylabel('Stress (kN/mm^2)')
plt.title(f'Stress-Strain Relation of Element {6} Top & Bottom Fibers')
plt.grid()
plt.legend()
plt.show()
"""
#--------------------------------------------------------------------
