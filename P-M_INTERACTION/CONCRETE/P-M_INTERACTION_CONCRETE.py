###########################################################################################################
#                                             IN THE NAME OF ALLAH                                        #
#             COMPARATIVE ANALYSIS OF AXIAL FORCE-MOMENT (P-M) INTERACTION BEHAVIOR IN CONFINED           #
#  REINFORCED CONCRETE CROSS-SECTIONS: EVALUATING STRAIN HARDENING EFFECTS AND ULTIMATE STRAIN CRITERIA   #
#---------------------------------------------------------------------------------------------------------#
#                          THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                     #
#                                   EMAIL: salar.d.ghashghaei@gmail.com                                   #
###########################################################################################################
import openseespy.opensees as ops
import matplotlib.pyplot as plt
import numpy as np
import ANALYSIS_FUNCTION as S02

# Define materials for nonlinear columns
# Define parameters (units: mm, N)
# ------------------------------------------
# CONCRETE                  tag   f'c        ec0   f'cu        ecu
# Core concrete (confined)
fcC = -27.6         # [N/mm²] Concrete Compressive Strength
ec0C = -0.0045      # [mm/mm] Concrete Compressive Strain
fcUC = -21          # [N/mm²] Concrete Compressive Ultimate Strength
ecuC = -0.015       # [mm/mm] Concrete Compressive Ultimate Strain

# Cover concrete (unconfined)
fcU = -18           # [N/mm²] Concrete Compressive Strength
ec0U = -0.0025      # [mm/mm] Concrete Compressive Strain
fcUU = -2           # [N/mm²] Concrete Compressive Ultimate Strength
ecuU = -0.008       # [mm/mm] Concrete Compressive Ultimate Strain
#ops.uniaxialMaterial('Concrete01', 2, fcU, ec0U, fcUU, ecuU)
 
# STEEL
# Reinforcing steel
fy = 400          # [N/mm²] Steel Rebar Yield Strength   
Es = 2e5          # [N/mm²] Modulus of Elasticity
ey = fy/Es        # [mm/mm] Steel Rebar Yield Strain
fu = 1.1818*fy    # [N/mm²] Steel Rebar Ultimate Strength
esu = ey*75.2     # [mm/mm] Steel Rebar Ultimate Strain
E = fy / ey  # Young's modulus
Esh = (fu - fy)/(esu - ey)
Bs = Esh / E

B = 400                 # [mm] Depth of the Section 
H = 500                 # [mm] Height of the Section  
cover = 50              # [mm] Concrete Section Cover
DIA = 25                # [mm] # Rebar Size Diameter
As = np.pi*(DIA**2)/4   # [mm²] Area of Rebar
# Define Analysis Properties
MAX_ITERATIONS = 1000      # Convergence iteration for test
MAX_TOLERANCE = 1.0e-10    # Convergence tolerance for test
#STEEL_KIND: 1 -> WITHOUT HARDENING AND ULTIMATE STRAIN
#STEEL_KIND: 2 -> WITH HARDENING AND ULTIMATE STRAIN
#%%------------------------------------------------------------------------------
def P_M_INTERACTION(b, h, cover, Es, fy, As, ECU, NUM, STEEL_KIND):
    # Clear previous model
    ops.wipe()
    # Define model builder
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    # Define two nodes at (0,0)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 0.0, 0.0)
    # Fix all degrees of freedom except axial and bending
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 0, 1, 0)
    
    secTag = 10
    coreTag = 1
    coverTag = 2
    steelTag = 3
    
    if STEEL_KIND == 1:# WITHOUT HARDENING AND ULTIMATE STRAIN
        ops.uniaxialMaterial('Steel01', steelTag, fy, Es, 0.0) 
    if STEEL_KIND == 2:# WITH HARDENING AND ULTIMATE STRAIN    
        pinchX = 0.8   # Pinching factor in X direction
        pinchY = 0.5   # Pinching factor in Y direction
        damage1 = 0.0  # Damage due to ductility
        damage2 = 0.0  # Damage due to energy
        beta = 0.1 # Stiffness degradation parameter
        ops.uniaxialMaterial('Hysteretic', steelTag, fy, ey, fu, esu, 0.2*fu, 1.1*esu, -fy, -ey, -fu, -esu, -0.2*fu, -1.1*esu, pinchX, pinchY, damage1, damage2, beta)
        # INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Hysteretic_Material


    ops.uniaxialMaterial('Concrete01', coreTag, fcC, ec0C, fcUC, ecuC)  # Core concrete (confined)
    ops.uniaxialMaterial('Concrete01', coverTag, fcU, ec0U, fcUU, ecuU) # Cover concrete (unconfined)
    # Some variables derived from the parameters
    y1 = h / 2.0
    z1 = b / 2.0
    NUMFIBERS = 20  # Number of layers for each fiber

    ops.section('Fiber', secTag)
    # Create the concrete core fibers
    ops.patch('rect', coreTag, NUMFIBERS, 5, cover - y1, cover - z1, y1 - cover, z1 - cover)

    # Create the concrete cover fibers (top, bottom, left, right)
    ops.patch('rect', coverTag, NUMFIBERS, 5, -y1, z1 - cover, y1, z1)
    ops.patch('rect', coverTag, NUMFIBERS, 5, -y1, -z1, y1, cover - z1)
    ops.patch('rect', coverTag, NUMFIBERS, 5, -y1, cover - z1, cover - y1, z1 - cover)
    ops.patch('rect', coverTag, NUMFIBERS, 5, y1 - cover, cover - z1, y1, z1 - cover)

    # Create the reinforcing fibers (left, middle, right)
    ops.layer('straight', steelTag, 3, As, y1 - cover, z1 - cover, y1 - cover, cover - z1)
    ops.layer('straight', steelTag, 2, As, 0.0, z1 - cover, 0.0, cover - z1)
    ops.layer('straight', steelTag, 3, As, cover - y1, z1 - cover, cover - y1, cover - z1)
    
    # Define element
    ops.element('zeroLengthSection', 1, 1, 2, secTag)
    EY = fy / Es  # Steel yield strain

    ep_list = np.linspace(ECU, 150*EY, NUM)

    FORCE, MOMENT, STRAIN, CUR, X = [], [], [], [], []

    ops.timeSeries('Constant', 2)
    ops.pattern('Plain', 2, 2)
    
    for ep in ep_list:
        ops.constraints('Transformation')
        ops.system('UmfPack')
        ops.algorithm('Newton')
        ops.numberer('Plain')
        ops.integrator('LoadControl', 1)
        ops.analysis('Static')

        Z = (ep - ECU) / h
        #print(ECU + Z * h / 2)
        x = ECU / Z  # Neutral Axis Depth
        ops.sp(2, 1, ECU + Z * h / 2)  # Strain at the centroid
        ops.sp(2, 3, Z)                 # Curvature
        # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/sp.html

        OK = ops.analyze(1)
        S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS

        ops.reactions()
        p = -ops.nodeReaction(2, 1)
        m = ops.nodeReaction(2, 3)
        strain = ops.nodeDisp(2, 1)
        cur = ops.nodeDisp(2, 3)

        FORCE.append(p)
        MOMENT.append(m)
        STRAIN.append(strain)
        CUR.append(cur)
        X.append(x)
        ops.wipeAnalysis()

    return FORCE, MOMENT, STRAIN, CUR, X
#%%------------------------------------------------------------------------------
# WITHOUT HARDENING AND ULTIMATE STRAIN
AXIAL, MOMENT, STRAIN, CUR, X = P_M_INTERACTION(B, H, cover, Es, fy, As, ecuC, NUM=5000, STEEL_KIND=1)
# WITH HARDENING AND ULTIMATE STRAIN
AXIAL02, MOMENT02, STRAIN02, CUR02, X02  = P_M_INTERACTION(B, H, cover, Es, fy, As, ecuC, NUM=5000, STEEL_KIND=2)

plt.figure(1, figsize=(12, 8))
#plt.plot(MOMENT, AXIAL, color='black')
plt.scatter(MOMENT, AXIAL, color='blue', linewidth=2)
plt.scatter(MOMENT02, AXIAL02, color='cyan', linestyle='--', linewidth=2)
plt.title('P-M Interaction')
plt.ylabel('Axial Force [N]')
plt.xlabel('Bending Moment [N.mm]')
plt.legend(['Steel01', 'Hysteretic'])
plt.grid()
plt.show()

plt.figure(2, figsize=(12, 8))
plt.plot(STRAIN, AXIAL, color='purple', linewidth=2)
plt.plot(STRAIN02, AXIAL02, color='magenta', linestyle='--', linewidth=2)
plt.title('AXIAL FORCE-STRAIN DIAGRAM')
plt.ylabel('Axial Force [N]')
plt.xlabel('Strain [mm/mm]')
plt.legend(['Steel01', 'Hysteretic'])
plt.grid()
plt.show()

plt.figure(3, figsize=(12, 8))
plt.plot(CUR, MOMENT, color='red', linewidth=2)
plt.plot(CUR02, MOMENT02, color='orange', linestyle='--', linewidth=2)
plt.title('MOMENT-CURVATURE DIAGRAM')
plt.ylabel('Moment [kN]')
plt.xlabel('Curvature [1/mm]')
plt.legend(['Steel01', 'Hysteretic'])
plt.grid()
plt.show()

plt.figure(4, figsize=(12, 8))
plt.plot(X, CUR, color='green', linewidth=2)
plt.plot(X02, CUR02, color='lime', linestyle='--', linewidth=2)
plt.title('MOMENT-CURVATURE DIAGRAM')
plt.ylabel('Curvature [1/mm]')
plt.xlabel('Neutral Axis Depth [mm]')
plt.legend(['Steel01', 'Hysteretic'])
plt.grid()
plt.show()

#%%------------------------------------------------------------------------------   
# EXPORT DATA TO EXCEL
import pandas as pd
DATA_TOTAL = {
    'NEUTRAL-AXIS_W': X,
    'NEUTRAL-AXIS_WO': X02,
    'STRAIN_W': STRAIN,
    'STRAIN_WO': STRAIN02,
    'CURVATURE_W': CUR,
    'CURVATURE_WO': CUR02,
    'AXIAL_FORCE_W': AXIAL,
    'AXIAL_FORCE_WO': AXIAL02,
    'MOMENT_W': MOMENT,
    'MOMENT_WO': MOMENT02,
}
# Convert to DataFrame
results_df = pd.DataFrame(DATA_TOTAL)
# Export the DataFrame to an Excel file
results_df.to_excel('P-M_INTERACTION_RESULTS.xlsx', index=False) 