####################################################################################
#                                 IN THE NAME OF ALLAH                             #
#           CONCRETE SECTION REBAR OPTIMIZATION BASED ON DEMAND DISPLACEMENT       #
#     THERMAL LOADING STRUCTURAL ANALYSIS OF A 2D CONCRETE FRAME USING OPENSEES    #
#----------------------------------------------------------------------------------#
#        OPTIMIZATION ALOGORITHM: NEWTON-RAPHSON METHOD PARALLEL COMPUTING         #
#----------------------------------------------------------------------------------#
# In all the beams of the floors, the extended load due to the dead load has been  #
# applied, and only in the beams of the first floor, the load due to the heat of   #
# the fire has been applied.                                                       #
#----------------------------------------------------------------------------------#
#       THIS PROGRAM IS WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)              #
#                      EMAIL: SALAR.D.GHASHGHAEI@GMAIL.COM                         #
####################################################################################

"""
 Models and Analyzes a 2D Concrete Frame subjected to Thermal and Distributed Loads using OpenSees and Find optimum Steel Rebar Diameter
 with Newthon-raphson Method 
Key points:  
1. Model Definition: The 2D frame has specified node coordinates for stories and bays, with fixed supports at the base.
 Material properties for concrete (with thermal effects) and concrete rectangular section geometries are defined using fiber elements.  
2. Element and Load Setup: Beam-column elements with corotational geometric transformation and Lobatto beam integration
 are created. Distributed loads are applied to beams, and a thermal gradient is applied to the first story beams.  
3. Analysis Setup: The analysis uses static load control with thermal increments, and the Newton-Raphson algorithm ensures convergence.
 Convergence tolerances and maximum iterations are defined.  
4. Output and Post-processing: Displacements, reactions, and deformations are recorded during the analysis.
 Data is extracted from output files for plotting base reactions (axial, shear, moment) and node displacements
 against temperature or applied load.  
5. Visualization: The frame's undeformed and deformed shapes are plotted, and results like temperature-displacement
 relationships and base reactions are visualized.
"""

import time as TI
import numpy as np
import matplotlib.pyplot as plt
import openseespy.opensees as ops
from Analysis_Function import ANALYSIS
from CONCRETE_FIBERTHERMAL_SECTION import R_RECTANGULAR_CONCRETE_SECTION_REBAR_B, R_RECTANGULAR_CONCRETE_SECTION_REBAR_C

#--------------------------------------------------------------------
# Define parameters (units: mm, kN)
num_stories = 3        # Number of Stories
num_bays = 4           # Number of Bays
story_height = 3000.0  # [mm] Height of each story
bay_width = 7000.0     # [mm] Width of each bay
#--------------------------------------------------------------------
# Define  Steel Rebar Material Properties (Steel01Thermal)
fy = 0.4               # [kN/mm²] Yield strength of steel rebar
fu = 1.5 * fy          # [kN/mm²] Ultimate strength of steel rebar
Es = 200               # [kN/mm²] Modulus of elasticity of steel rebar
ey = fy / Es           # [mm/mm] Yield steel strain
esu = 0.35             # [mm/mm] Ultimate steel strain
Esh = (fu - fy) / (esu - ey)  # [kN/mm^2] Strain hardening modulus
b = Esh / Es                  # Strain hardening ratio
#--------------------------------------------------------------------
# Define material properties for (Concrete02Thermal)
fcp = -0.03            # [kN/mm²] Compressive strength
epsc0 = -0.0025        # [mm/mm] Strain at peak compressive strength
E0 = fcp / epsc0
fpcu = fcp * 0.05      # [kN/mm²] Crushing strength
epsU = -0.004          # [mm/mm] Strain at ultimate stress
lamda = 0.1            # Ratio between unloading slope at epsU and initial slope
ft = 0.005             # [kN/mm²] Tensile strength
Ets = ft/0.002         # [kN/mm²] Tension softening stiffness 
#--------------------------------------------------------------------
# Define Thermal and Distributed load
Max_Thermal = 800.0        # [°C] Temperature
distributed_load = -0.003  # [kN/mm] Uniform Distributed Loads
#--------------------------------------------------------------------
# Define Analysis Properties
MAX_ITERATIONS = 100       # Convergence iteration for test
MAX_TOLERANCE = 1.0e-10    # Convergence tolerance for test

Nstep = 400                # Number of incremental steps
Incr_Temp = 1/Nstep        # Incremental temperature step
#--------------------------------------------------------------------
# Define a rectangular concrete section using fibers.
B = 400       # [mm] Width of the rectangular section
H = 500       # [mm] Height of the rectangular section
COVER = 30    # [mm] Cover of the rectangular section
RD = 25       # [mm] Diameter of the rebar
#--------------------------------------------------------------------
# Thermal Loading Analysis
def TEMP_ANAL(B, H, COVER, RD, DEMAND_NODE):
    # Clear pervious model
    ops.wipe()
    # Define model
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    # Node coordinates
    node_id = 1
    for i in range(num_stories + 1):  # Including ground level
        for j in range(num_bays + 1):  # Including leftmost column
            ops.node(node_id, j * bay_width, i * story_height)
            if i == 0:  # Fix supports at ground level
                ops.fix(node_id, 1, 1, 1)
            else:  # Free to move for higher levels
                ops.fix(node_id, 0, 0, 0)
            node_id += 1

    #print(' Structure Coordinates Done.') 
    #--------------------------------------------------------------------
    # Materials (STEEL MATERIAL NONLINEARITY)
    matTag01 = 1
    ops.uniaxialMaterial('Steel01Thermal', matTag01, fy, Es, b) # Steel with bilinear kinematic hardening Material with thermaal effect
    #--------------------------------------------------------------------
    # Materials (CONCRETE MATERIAL NONLINEARITY)
    matTag02 = 2
    """
    Concrete02Thermal:
    Concrete02Thermal is created for modelling concrete, which is derived from
     the standard "Concrete02" material and incorprates with temperature dependent
     properties sugggested in Eurocode 2 (EN1992-1-2)
    """
    ops.uniaxialMaterial('Concrete02Thermal', matTag02, fcp, epsc0, fpcu, epsU, lamda, ft, Ets)
    #ops.uniaxialMaterial('Concrete01Thermal', matTag02, fcp, epsc0, fpcu, epsU)
    """
    ConcreteECThermal:
    ConcreteECThermal is derived by modification of the existing concrete material
     class "Concrete02" and "Concrete02Thermal" to include the temperature-dependent
     properties according to EN 1992-1-2 about concrete with siliceous aggregates at elevated temperature.
    """
    #ops.uniaxialMaterial('ConcreteECThermal', matTag02, fcp, epsc0, fpcu, epsU, lamda, ft, Ets) # Concrete hardening Material with thermal effect
    #--------------------------------------------------------------------
    # Define a rectangular concrete section using fibers.
    B = 400       # [mm] Width of the rectangular section
    H = 500       # [mm] Height of the rectangular section
    COVER = 30    # [mm] Cover of the rectangular section
    #RD = 25       # [mm] Diameter of the rebar
    NUM_B = 3     # Number of fibers along the width of the section
    NUM_H = 100   # Number of fibers along the height of the section
    #secTag01 = 1    # Section tag identifier
    # Concrete Thermal Fiber Section without Rebars
    #Depth = R_RECTANGULAR_CONCRETE_SECTION(secTag01, B, H, NUM_B, NUM_H, matTag01)
    # Concrete Thermal Fiber Section with Rebars 

    # BEAMS SECTION
    secTag01 = 1    # Section beam tag identifier
    # Concrete Fiber Section without Rebars
    #Depth = R_RECTANGULAR_CONCRETE_SECTION(secTag01, B, H, NUM_B, NUM_H, matTag01)
    # Concrete Fiber Section with Rebars 
    NUM_LAYERS = 2  # Number of rebar layers
    Depth01 = R_RECTANGULAR_CONCRETE_SECTION_REBAR_B(secTag01, B, H, COVER, RD, NUM_B, NUM_H, NUM_LAYERS, matTag02, matTag01, PLOT=False)
    # COLUMNS SECTION
    secTag02 = 2    # Section column tag identifier
    NUM_LAYERS = 4  # Number of rebar layers
    NUM_FIRST = 4   # Number of rebars in first and last layer
    Depth02 = R_RECTANGULAR_CONCRETE_SECTION_REBAR_C(secTag02, B, B, COVER, RD, NUM_B, NUM_H, NUM_LAYERS, NUM_FIRST, matTag02, matTag01, PLOT=False)

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
        ops.eleLoad('-ele', ele_id, '-type', '-beamUniform', distributed_load)
    #--------------------------------------------------------------------    
    # A linear thermal gradient is applied to the elements in the first story.
    # The temperature varies from -DD to DD across the section height.
    # INFO LINK: https://openseesforfire.github.io/Subpages/ThermalActionCmds.html

    # Apply thermal load to all elements
    #for ele_id in range(1, element_id):
    #    ops.eleLoad('-ele', ele_id, '-type', '-beamThermal', Max_Thermal, -DD, Max_Thermal, DD)

    # Apply thermal load only to the identified beam elements
    DD = 0.5 * Depth01
    for ele_id in beam_id:
        if ele_id == 6 or ele_id == 7 or ele_id == 8 or ele_id == 9: # Apply Thermal Loads only to the Beams of Story 1
            # eleLoad -ele $eleTag -type -beamThermal $T1 $y1 $T2 $Y2
            ops.eleLoad('-ele', ele_id, '-type', '-beamThermal', Max_Thermal, -DD, Max_Thermal, DD)  

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
    """
    ops.recorder('Node', '-file', "DTH_PUSH.txt",'-time', '-node', 6, '-dof', 1,2,3, 'disp')        # Displacement Time History Node 6
    ops.recorder('Node', '-file', "BTH_PUSH_01.txt",'-time', '-node', 1, '-dof', 1,2,3, 'reaction') # Base Reaction Time History Node 1
    ops.recorder('Node', '-file', "BTH_PUSH_02.txt",'-time', '-node', 2, '-dof', 1,2,3, 'reaction') # Base Reaction Time History Node 2
    ops.recorder('Node', '-file', "BTH_PUSH_03.txt",'-time', '-node', 3, '-dof', 1,2,3, 'reaction') # Base Reaction Time History Node 3
    ops.recorder('Node', '-file', "BTH_PUSH_04.txt",'-time', '-node', 4, '-dof', 1,2,3, 'reaction') # Base Reaction Time History Node 4
    ops.recorder('Node', '-file', "BTH_PUSH_05.txt",'-time', '-node', 5, '-dof', 1,2,3, 'reaction') # Base Reaction Time History Node 5
    ops.recorder('Element', '-file', 'STRESS_STRAIN_BOT.txt', '-time', '-ele', 6, 'section', secTag01, 'fiber', -DD, 0, 'stressStrainTangent')
    ops.recorder('Element', '-file', 'STRESS_STRAIN_TOP.txt', '-time', '-ele', 6, 'section', secTag01, 'fiber', +DD, 0, 'stressStrainTangent')
    ops.recorder('Element', '-file', 'STRESS_STRAIN_BOT.txt', '-time', '-ele', 6, 'section', secTag01, 'fiber', -DD, 0, 'stressStrain')
    ops.recorder('Element', '-file', 'STRESS_STRAIN_TOP.txt', '-time', '-ele', 6, 'section', secTag01, 'fiber', +DD, 0, 'stressStrain')
    """
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
        dispX.append(ops.nodeDisp(DEMAND_NODE, 1)) # X Displacement
        dispY.append(ops.nodeDisp(DEMAND_NODE, 2)) # Y Displacement
        
    #ops.wipe()   
    
    return temp, dispX, dispY, mid_node

#--------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------
# FIND THE OPTIMUM VALUE (NEWTON-RAPHSON SOLVER FOR OPTIMAL REBAR DIAMETER) WITH PARALLEL COMPUTING
# ----------------------------------------------------------------------------------------------------

import concurrent.futures
from multiprocessing import freeze_support

ESP = 1e-3        # Finite difference derivative Convergence Tolerance
TOLERANCE = 1e-6  # Convergence Tolerance
RESIDUAL = 100    # Convergence Residual 
IT = 0            # Intial Iteration
ITMAX = 100000    # Max. Iteration
DEMAND = 60       # [mm] Demand Lateral Displacement
DEMAND_NODE = 6   # Demand Displacement node

# Helper function must be defined outside main guard
def Run_Analysis(X):
    """Helper function for parallel execution"""
    DATA = TEMP_ANAL(B, H, COVER, X, DEMAND_NODE)
    temp, dispX, dispY, mid_node = DATA
    SUPPLY = np.max(np.abs(dispX))
    print(f'SUPPLY: {SUPPLY:.5f}')
    return SUPPLY  # Return SUPPLY value
    
def Optimize_Rebar_Diameter():
    X = RD            # [mm] Intial Guess for Rebar Diameter
    ESP = 1e-3        # Finite difference derivative Convergence Tolerance
    TOLERANCE = 1e-6  # Convergence Tolerance
    RESIDUAL = 100    # Convergence Residual 
    IT = 0            # Intial Iteration
    ITMAX = 100000    # Max. Iteration
    #DEMAND = 60       # [mm] Demand Lateral Displacement
    #DEMAND_NODE = 6   # Demand Displacement node
    

    with concurrent.futures.ProcessPoolExecutor() as executor:
        while RESIDUAL > TOLERANCE and IT < ITMAX:
            params = [X, X - ESP, X + ESP]
            
            futures = [executor.submit(Run_Analysis, param) for param in params]
            SUPPLY, SUPPLYmin, SUPPLYmax = [f.result() for f in futures]

            F = SUPPLY - DEMAND
            Fmin = SUPPLYmin - DEMAND
            Fmax = SUPPLYmax - DEMAND

            DF = (Fmax - Fmin) / (2 * ESP)
            DX = F / DF if abs(DF) > 1e-12 else 0
            RESIDUAL = abs(DX)
            X -= DX
            IT += 1

            print(f'IT: {IT} - RESIDUAL: {RESIDUAL:.6e} - SECTION REBAR DIAMETER: {X:.6e}')
            
            if RESIDUAL < TOLERANCE:
                print(f'\t\t Optimum Section Rebar Diameter:   {X:.4f}')
                print(f'\t\t Iteration Counts:                 {IT}')
                print(f'\t\t Convergence Residual:             {RESIDUAL:.10e}')
                break

if __name__ == '__main__':
    current_time = TI.strftime("%H:%M:%S", TI.localtime())
    print(f"Current time (HH:MM:SS): {current_time}\n\n")
    

    freeze_support()  # Required for Windows support
    Optimize_Rebar_Diameter()

    
    current_time = TI.strftime("%H:%M:%S", TI.localtime())
    print(f"Current time (HH:MM:SS): {current_time}\n\n")    

#%%--------------------------------------------------------------------


