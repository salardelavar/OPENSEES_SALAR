######################################################################################################################
#                            >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                        #
#                        IMPLEMENTATION OF RESPONSE SPECTRUM ANALYSIS WITH CQC MODAL COMBINATION                     #
#                                       FOR 2D REINFORCED CONCRETE FRAMES USING OPENSEES                             #
#--------------------------------------------------------------------------------------------------------------------#
#                              THIS PROGRAM CHANGED BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                            #
#                                       EMAIL: salar.d.ghashghaei@gmail.com                                          #
######################################################################################################################
"""
1. Response Spectrum Analysis Implementation: This code performs comprehensive response spectrum analysis using OpenSees for seismic evaluation of structures.
2. Modal Analysis Foundation: It begins with eigenvalue analysis to extract natural vibration modes and periods essential for dynamic characterization.
3 CQC Modal Combination: Implements the Complete Quadratic Combination method for statistically combining modal responses, accounting for closely-spaced modes.
4. 2D Concrete Frame Modeling: Creates a simplified yet realistic 1-bay 1-story reinforced concrete moment frame with appropriate material properties.
5. Multiple Analysis Approaches: Demonstrates three different methodologies for response spectrum analysis - direct input, time series, and mode-by-mode approaches.
6. Realistic Structural Properties: Uses code-compliant concrete modulus (Ec = 4700√fc) and practical member dimensions for engineering accuracy.
7. Damping Considerations: Incorporates 5% modal damping ratio typical for concrete structures in seismic analysis.
8. Force Response Extraction: Focuses on critical section forces, particularly bending moments at column bases for seismic design.
9. Spectral Matching: Utilizes given response spectrum data to ensure analysis matches specific seismic hazard characteristics.
10. Verification Framework: Provides comparative analysis between different implementation methods for result validation and quality assurance.

This represents advanced earthquake engineering practice combining theoretical dynamics with
 practical implementation for seismic performance assessment.
"""
#%%------------------------------------------------------------------------------
# INFO LINK:
'https://openseespydoc.readthedocs.io/en/latest/src/responseSpectrumAnalysis.html'
# YOUTUBE VIDEOS:
'https://www.youtube.com/watch?v=I7aG25MasNM' 
'https://www.youtube.com/watch?v=48DNGzmznZ0'   
#%%------------------------------------------------------------------------------
import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
import PLOT_2D as S01


# Given data
Tn = [0.0, 0.06, 0.1, 0.12, 0.18, 0.24, 0.3, 0.36, 0.4, 0.42, 
      0.48, 0.54, 0.6, 0.66, 0.72, 0.78, 0.84, 0.9, 0.96, 1.02, 
      1.08, 1.14, 1.2, 1.26, 1.32, 1.38, 1.44, 1.5, 1.56, 1.62, 
      1.68, 1.74, 1.8, 1.86, 1.92, 1.98, 2.04, 2.1, 2.16, 2.22, 
      2.28, 2.34, 2.4, 2.46, 2.52, 2.58, 2.64, 2.7, 2.76, 2.82, 
      2.88, 2.94, 3.0, 3.06, 3.12, 3.18, 3.24, 3.3, 3.36, 3.42, 
      3.48, 3.54, 3.6, 3.66, 3.72, 3.78, 3.84, 3.9, 3.96, 4.02, 
      4.08, 4.14, 4.2, 4.26, 4.32, 4.38, 4.44, 4.5, 4.56, 4.62, 
      4.68, 4.74, 4.8, 4.86, 4.92, 4.98, 5.04, 5.1, 5.16, 5.22, 
      5.28, 5.34, 5.4, 5.46, 5.52, 5.58, 5.64, 5.7, 5.76, 5.82, 
      5.88, 5.94, 6.0]

Sa = [1.9612, 3.72628, 4.903, 4.903, 4.903, 4.903, 4.903, 4.903, 4.903, 4.6696172, 
      4.0861602, 3.6321424, 3.2683398, 2.971218, 2.7241068, 2.5142584, 2.3348086, 2.1788932, 2.0425898, 1.9229566, 
      1.8160712, 1.7199724, 1.6346602, 1.5562122, 1.485609, 1.4208894, 1.3620534, 1.3071398, 1.2571292, 1.211041, 
      1.166914, 1.1267094, 1.0894466, 1.054145, 1.0217852, 0.990406, 0.960988, 0.9335312, 0.9080356, 0.8835206, 
      0.8599862, 0.838413, 0.8168398, 0.7972278, 0.7785964, 0.759965, 0.7432948, 0.7266246, 0.710935, 0.6952454, 
      0.6805364, 0.666808, 0.6540602, 0.6285646, 0.6040496, 0.5814958, 0.5609032, 0.5403106, 0.5206986, 0.5030478, 
      0.485397, 0.4697074, 0.4540178, 0.4393088, 0.4255804, 0.411852, 0.3991042, 0.3863564, 0.3755698, 0.3638026, 
      0.353016, 0.34321, 0.333404, 0.3245786, 0.3157532, 0.3069278, 0.2981024, 0.2902576, 0.2833934, 0.2755486, 
      0.2686844, 0.2618202, 0.254956, 0.2490724, 0.2431888, 0.2373052, 0.2314216, 0.2265186, 0.220635, 0.215732, 
      0.210829, 0.205926, 0.2020036, 0.1971006, 0.1931782, 0.1892558, 0.1853334, 0.181411, 0.1774886, 0.1735662, 
      0.1706244, 0.166702, 0.1637602]

# Plotting the response spectrum
plt.figure(figsize=(10, 6))
plt.plot(Tn, Sa, marker='o', color='b', label='Response Spectrum')
plt.title('Response Spectrum')
plt.xlabel('Natural Period (Tn) [s]')
plt.ylabel('Spectral Acceleration (Sa) [g]')
#plt.semilogx();plt.semilogy();
plt.grid(True)
plt.legend()
plt.show()

#%%-------------------------------------------------------------------
# Define concrete material properties
# Concrete material - using Elastic for simplicity (can be changed to Concrete01/Concrete02 for nonlinear)

# Define parameters (units: mm, N)
# ------------------------------------------
mass = 500000.0          # [kg] (lumped mass)
fc = 30                  # [N/mm²] Concrete compressive strength
Ec = 4700 * np.sqrt(fc)  # [N/mm²] Concrete elastic modulus
# Define steel material for reinforcement (if needed)
fy = 420                 # [N/mm²] Steel Yield Strength
Es = 200e3               # [N/mm²] Steel Elastic Modulus
# Define geometric dimensions
bay_width = 5000.0       # [mm] Beam Length
story_height = 3000.0    # [mm] Column Length
Bb = 300.0               # [mm] Beam Setion Width
Hb = 400.0               # [mm] Beam Setion Height
Bc = 400.0               # [mm] Column Setion Width
Hc = 400.0               # [mm] Column Setion Height

NUM_MODES = 2           # Number of Modes
DAMPIMND_RATIO = 0.05   # Damping Ratio

# DEFINE ANALYSIS PROPERTIES
MAX_ITERATIONS = 5000    # Convergence iteration for test
MAX_TOLERANCE = 1.0e-6   # Convergence tolerance for test
#%%-------------------------------------------------------------------
# Reset model
ops.wipe()
# define a 2D model
ops.model("basic","-ndm",2,"-ndf",3)

# Define time series for response spectrum analysis
# Tn and Sa should be defined with your specific values
# Example: Tn = [0.1, 0.2, 0.3, ...], Sa = [0.5, 0.4, 0.3, ...]
ops.timeSeries("Path", 1, "-time", *Tn, "-values", *Sa)

ops.uniaxialMaterial("Elastic", 1, Ec)
#ops.uniaxialMaterial("Elastic", 2, Es)

# Define section properties for columns and beams
# Column section: 400mm x 400mm
A_col = Bc * Hc          # [mm²] Section Area
I_col = Bc * Hc**3 / 12  # [mm⁴]

# Beam section: 300mm x 400mm  
A_beam = Bb * Hb          # [mm²] Section Area
I_beam = Bb * Hb**3 / 12  # [mm⁴]

# Define sections
ops.section("Elastic", 1, Ec, A_col, I_col)  # Column section
ops.section("Elastic", 2, Ec, A_beam, I_beam)  # Beam section

# Define nodes
ops.node(1, 0.0, 0.0)  # Base left
ops.node(2, bay_width, 0.0)  # Base right  
ops.node(3, 0.0, story_height)  # Top left
ops.node(4, bay_width, story_height)  # Top right

# Define masses (concentrated at top nodes)
ops.mass(3, mass, mass, 0.0)
ops.mass(4, mass, mass, 0.0)

# Define fixities (fixed base)
ops.fix(1, 1, 1, 1)  # Fixed left base
ops.fix(2, 1, 1, 1)  # Fixed right base

# Define geometric transformation
ops.geomTransf("Linear", 1)  # 2D linear transformation

# Define beam integration
ops.beamIntegration("Lobatto", 1, 1, 5)  # For columns
ops.beamIntegration("Lobatto", 2, 2, 5)  # For beam

# Create elements
# Columns
ops.element("forceBeamColumn", 1, 1, 3, 1, 1)  # Left column
ops.element("forceBeamColumn", 2, 2, 4, 1, 1)  # Right column

# Beam
ops.element("forceBeamColumn", 3, 3, 4, 1, 2)  # Beam

# Define analysis settings
ops.constraints("Transformation")
ops.numberer("RCM")
ops.system("UmfPack")
ops.test("NormUnbalance", MAX_TOLERANCE, MAX_ITERATIONS)
ops.algorithm("Linear")
ops.integrator("LoadControl", 0.0)
ops.analysis("Static")

# Run eigenvalue analysis with 2 modes
eigs = ops.eigen("-genBandArpack", NUM_MODES)
print(eigs)

# Compute modal properties
ops.modalProperties("-print", "-file", "SALAR_ModalReport.txt", "-unorm")

# Define recorder for element forces
# COLUMN
filename = 'SALAR_ele_1_sec_1.txt'
ops.recorder('Element', '-file', filename, '-closeOnWrite', '-precision', 16, '-ele', 1, 'section', '1', 'force')
filename = 'SALAR_ele_2_sec_1.txt'
ops.recorder('Element', '-file', filename, '-closeOnWrite', '-precision', 16, '-ele', 2, 'section', '1', 'force')
# BEAM
filename = 'SALAR_ele_3_sec_1.txt'
ops.recorder('Element', '-file', filename, '-closeOnWrite', '-precision', 16, '-ele', 3, 'section', '2', 'force')
# Response spectrum analysis settings
tsTag = 1
direction = 1  # excited DOF = Ux

# Damping for each mode
dmp = [DAMPIMND_RATIO] * len(eigs)
scalf = [1.0] * len(eigs)

# CQC function 
def CQC(mu, lambdas, dmp, scalf):
    u = 0.0
    ne = len(lambdas)
    for i in range(ne):
        for j in range(ne):
            di = dmp[i]
            dj = dmp[j]
            bij = lambdas[i] / lambdas[j]
            rho = ((8.0 * np.sqrt(di * dj) * (di + bij * dj) * (bij ** (3.0 / 2.0))) /
                   ((1.0 - bij ** 2.0) ** 2.0 + 4.0 * di * dj * bij * (1.0 + bij ** 2.0) + 
                    4.0 * (di ** 2.0 + dj ** 2.0) * bij ** 2.0))
            u += scalf[i] * mu[i] * scalf[j] * mu[j] * rho
    return np.sqrt(u)

#%%-------------------------------------------------------------------
# TEST 00 - Response Spectrum Analysis with Tn and Sa lists

ops.responseSpectrumAnalysis(direction, '-Tn', *Tn, '-Sa', *Sa)

# Read results
My = []
with open(filename, 'r') as f:
    lines = f.read().split('\n')
    for line in lines:
        if len(line) > 0:
            tokens = line.split(' ')
            My.append(float(tokens[NUM_MODES-1]))

# CQC modal combination
MyCQC = CQC(My, eigs, dmp, scalf)

print('\n\nTEST 00:\nRun a Response Spectrum Analysis for all modes.')
print('Do CQC combination in post processing.')
print('Use Tn and Sa lists.\n')
print('{0: >10}{1: >15}'.format('Mode', 'My'))
for i in range(len(eigs)):
    print('{0: >10}{1: >15f}'.format(i + 1, My[i]))
print('{0: >10}{1: >15f}'.format('CQC ', MyCQC))

#%% ========================================================================
# TEST 01 - Response Spectrum Analysis with Path timeSeries

ops.responseSpectrumAnalysis(tsTag, direction)

# Read results
My = []
with open(filename, 'r') as f:
    lines = f.read().split('\n')
    for line in lines:
        if len(line) > 0:
            tokens = line.split(' ')
            My.append(float(tokens[NUM_MODES-1]))

# CQC modal combination
MyCQC = CQC(My, eigs, dmp, scalf)

print('\n\nTEST 01:\nRun a Response Spectrum Analysis for all modes.')
print('Do CQC combination in post processing.')
print('Use a Path timeSeries to store Tn-Sa pairs.\n')
print('{0: >10}{1: >15}'.format('Mode', 'My'))
for i in range(len(eigs)):
    print('{0: >10}{1: >15f}'.format(i + 1, My[i]))
print('{0: >10}{1: >15f}'.format('CQC ', MyCQC))

#%%-------------------------------------------------------------------
# TEST 02 - Mode-by-mode Response Spectrum Analysis

My = []
for i in range(len(eigs)):
    ops.responseSpectrumAnalysis(direction, '-Tn', *Tn, '-Sa', *Sa, '-mode', i + 1)
    force = ops.eleResponse(1, 'section', '1', 'force')
    My.append(force[NUM_MODES-1])

# CQC modal combination
MyCQC = CQC(My, eigs, dmp, scalf)

print('\n\nTEST 02:\nRun a Response Spectrum Analysis mode-by-mode.')
print('Grab results during the loop and do CQC combination with them.\n')
print('{0: >10}{1: >15}'.format('Mode', 'My'))
for i in range(len(eigs)):
    print('{0: >10}{1: >15f}'.format(i + 1, My[i]))
print('{0: >10}{1: >15f}'.format('CQC ', MyCQC))

# Clean up
#ops.wipe()
#%%-------------------------------------------------------------------
# %% Plot 2D Frame Shapes for Response Spectrum Analysis
S01.PLOT_2D_FRAME(deformed_scale=100.0)
#%%-------------------------------------------------------------------
# Print out the state of nodes 3 and 4
ops.printModel("node",1, 2, 3, 4)
# Print out the state of element 1 , 2 and 3
ops.printModel("ele", 1, 2 , 3)
# Print the Model
#printModel()
ops.printModel("-JSON", "-file", "RESPONSE_SPECTRUM_ANALYSIS.json")
#%%-------------------------------------------------------------------------------