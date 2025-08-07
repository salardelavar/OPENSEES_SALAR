###########################################################################################################
#                      >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                   #
#             OPTIMIZATION OF PLATE THICKNESS WITH COMPARATIVE ANALYSIS OF AXIAL FORCE–MOMENT (P–M)       #
#                       INTERACTION IN CONFINED REINFORCED CONCRETE WITH PLATE                            #
#     COMPOSITE SECTIONS: EVALUATING STRAIN HARDENING AND ULTIMATE STRAIN CRITERIA USING OPENSEES         #
#---------------------------------------------------------------------------------------------------------#
#                                OPTIMIZATION ALGORITHM : NEWTON-RAPHSON METHOD                           # 
#---------------------------------------------------------------------------------------------------------#
#                          THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                     #
#                                   EMAIL: salar.d.ghashghaei@gmail.com                                   #
###########################################################################################################
"""
Performs a comparative analysis of reinforced concrete columns' axial force-moment (P-M)
 interaction behavior using OpenSeesPy.

1. Objective: 
    Evaluates how steel plate , rebar reinforcement strain hardening and ultimate strain criteria affect P-M
    interaction capacity, crucial for seismic design where ductility matters.

2. Materials:
   - Concrete: Confined (core) and unconfined (cover) modeled with `Concrete01` material laws.
   - *Steel: Two models - `Steel01` (elastic-perfectly plastic) vs. `Hysteretic` (includes hardening & fracture).

3. Section Modeling: Creates a rectangular RC cross-section with fiber discretization. Core concrete is confined,
 cover concrete unconfined, with rebars placed at edges and mid-depth.

4. Analysis Method: Uses a displacement-controlled approach to simulate increasing curvature/strain. For each strain level:
   - Applies strain compatibility (plane sections remain plane)
   - Computes axial force (P) and moment (M) using nonlinear static analysis.

5. Key Outputs:
   - P-M interaction diagrams
   - Moment-curvature relationships
   - Neutral axis depth trends
   - Flexural rigidity (EI) variations

6. Comparison: Contrasts two steel models:
   - *Without hardening*: Brittle failure, lower ductility
   - *With hardening*: Enhanced ductility, gradual strength degradation

7. Seismic Relevance: Strain hardening improves energy dissipation capacity - critical for
 earthquake-resistant structures. Ultimate strain criteria help predict failure thresholds under cyclic loading.

8. Visualization: Generates 6 plots including P-M curves, moment-curvature, and EI relationships
 to visualize performance differences.

9. Data Export: Results exported to Excel for quantitative comparison (e.g., moment capacity at 2% drift).

10. Application: Helps engineers choose appropriate material models for performance-based seismic design,
 assessing how confinement and steel hardening influence collapse prevention mechanisms.

This analysis directly supports seismic detailing decisions by quantifying how material nonlinearities
 affect global section behavior under combined loading.
"""
import time as TI
import openseespy.opensees as ops
import matplotlib.pyplot as plt
import numpy as np
import ANALYSIS_FUNCTION as S02
import COMPOSITE_CONCRETE_SECTION_FUN as S03

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
 
# STEEL
# Reinforcing steel
fy = 4000         # [N/mm²] Steel Rebar Yield Strength   
Es = 2e5          # [N/mm²] Modulus of Elasticity
ey = fy/Es        # [mm/mm] Steel Rebar Yield Strain
fu = 1.1818*fy    # [N/mm²] Steel Rebar Ultimate Strength
esu = ey*75.2     # [mm/mm] Steel Rebar Ultimate Strain
Esh = (fu - fy)/(esu - ey)
Bs = Esh / Es

# Steel Plate
fyP = 240          # [N/mm²] Steel Rebar Yield Strength   
EsP = 2e5          # [N/mm²] Modulus of Elasticity
eyP = fyP/EsP      # [mm/mm] Steel Rebar Yield Strain
fuP = 1.1818*fyP   # [N/mm²] Steel Rebar Ultimate Strength
esuP = eyP*41.8    # [mm/mm] Steel Rebar Ultimate Strain
EshP = (fuP - fyP)/(esuP - eyP)
BsP = EshP / EsP

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
def P_M_INTERACTION(b, h, cover, plateThickness, Es, fy, As, ECU, NUM, STEEL_KIND):
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
    steelPlateTag = 4
    numBarsTop, barAreaTop = 5, np.pi *(18**2)/4
    numBarsBot, barAreaBot = 5, np.pi *(20**2)/4
    numBarsIntTot, barAreaInt = 4, np.pi *(5**2)/4
    
    if STEEL_KIND == 1:# WITHOUT HARDENING AND ULTIMATE STRAIN
        ops.uniaxialMaterial('Steel01', steelTag, fy, Es, 0.0) 
        ops.uniaxialMaterial('Steel01', steelPlateTag, fyP, EsP, 0.0) 
    if STEEL_KIND == 2:# WITH HARDENING AND ULTIMATE STRAIN    
        pinchX = 0.8   # Pinching factor in X direction
        pinchY = 0.5   # Pinching factor in Y direction
        damage1 = 0.0  # Damage due to ductility
        damage2 = 0.0  # Damage due to energy
        beta = 0.1 # Stiffness degradation parameter
        ops.uniaxialMaterial('Hysteretic', steelTag, fy, ey, fu, esu, 0.2*fu, 1.1*esu, -fy, -ey, -fu, -esu, -0.2*fu, -1.1*esu, pinchX, pinchY, damage1, damage2, beta)
        
        ops.uniaxialMaterial('Hysteretic', steelPlateTag , fyP, eyP, fuP, esuP, 0.2*fuP, 1.1*esuP, -fyP, -eyP, -fuP, -esuP, -0.2*fuP, -1.1*esuP, pinchX, pinchY, damage1, damage2, beta)
        # INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Hysteretic_Material


    ops.uniaxialMaterial('Concrete01', coreTag, fcC, ec0C, fcUC, ecuC)  # Core concrete (confined)
    ops.uniaxialMaterial('Concrete01', coverTag, fcU, ec0U, fcUU, ecuU) # Cover concrete (unconfined)

    d, MASS = S03.CONCRETE_CONFINED_REC_PLATE_SECTION_QUAD(secTag, H, B, cover, cover, coreTag, coverTag, steelTag, steelPlateTag,
                                          numBarsTop, barAreaTop, numBarsBot, barAreaBot, numBarsIntTot, barAreaInt,
                                          plateThickness, PLOT=True, CONCRETE_DENSITY=2500/1e9, STEEL_DENSITY=7850/1e9)

    
    # Define element
    ops.element('zeroLengthSection', 1, 1, 2, secTag)
    EY = fy / Es  # Steel yield strain

    ep_list = np.linspace(ECU, 150*EY, NUM)

    FORCE, MOMENT, STRAIN, CUR, X, EI = [], [], [], [], [], []

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
        EI.append(np.abs(m)/np.abs(cur))
        ops.wipeAnalysis()

    return FORCE, MOMENT, STRAIN, CUR, X, EI
#%%------------------------------------------------------------------------------
# DAMPING RATIO OPTIMIZATION:
plateThickness = 10    
X = plateThickness # Intial Guess Plate Thickness
ESP = 1e-5         # Finite difference derivative Convergence Tolerance
TOLERANCE = 1e-10  # Convergence Tolerance
RESIDUAL = 100     # Convergence Residual 
IT = 0             # Intial Iteration
ITMAX = 100000     # Max. Iteration

AXIAL_TARGET = 15_000_000

# Analysis Durations:
starttime = TI.process_time()

### FIND THE OPTIMUM VALUE 
while (RESIDUAL > TOLERANCE):
    # X -------------------------------------------------------
    DATA = P_M_INTERACTION(B, H, cover, X, Es, fy, As, ecuC, NUM=5000, STEEL_KIND=2)
    AXIAL, MOMENT, STRAIN, CUR, NA, EI = DATA
    print(f'PLATE THICKNESS: {X:.8e}')
    AXIAL_DEMAND = np.max(np.abs(AXIAL))
    F = AXIAL_DEMAND - AXIAL_TARGET
    print('F: ', F)
    # Xmin -------------------------------------------------------
    XMIN = X - ESP 
    DATA = P_M_INTERACTION(B, H, cover, XMIN, Es, fy, As, ecuC, NUM=5000, STEEL_KIND=2)
    AXIAL, MOMENT, STRAIN, CUR, NA, EI = DATA
    AXIAL_DEMAND_MIN = np.max(np.abs(AXIAL))
    Fmin = AXIAL_DEMAND_MIN - AXIAL_TARGET
    print('Fmin: ', Fmin)
    # Xmax -------------------------------------------------------
    XMAX = X + ESP 
    DATA = P_M_INTERACTION(B, H, cover, XMAX, Es, fy, As, ecuC, NUM=5000, STEEL_KIND=2)
    AXIAL, MOMENT, STRAIN, CUR, NA, EI = DATA
    AXIAL_DEMAND_MAX = np.max(np.abs(AXIAL))
    Fmax = AXIAL_DEMAND_MAX - AXIAL_TARGET
    print('Fmax: ', Fmax)
    # DF -------------------------------------------------------
    DF = (Fmax - Fmin)/(2 * ESP);# Calculate the Finite difference derivative of F
    print('DF: ', DF)
    # DX -------------------------------------------------------
    DX = F / DF;        # Calculate dx
    print('DX: ', DX)
    # RESIDUAL -------------------------------------------------
    RESIDUAL = abs(DX); # Calculate residual
    
    print('IT: ', IT+1,' - RESIDUAL: ', RESIDUAL,' - PLATE THICKNESS: ', X,'\n')
    X -= DX;            # update X
    IT += 1;            # update iteration
    if IT == ITMAX:
        print('\t\t Iteration reached to Max. Iteration')
        print('\t\t Change ESP and TOLERANCE for better Convergence')
        X = -1
        break;
    if RESIDUAL < TOLERANCE:
        print(f'\t\t Optimum Plate Thickness:    {X:.6e}')
        print(f'\t\t Iteration Counts:           {IT}')
        print(f'\t\t Convergence Residual:       {RESIDUAL:.10e}')

totaltime = TI.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n')  

#%%------------------------------------------------------------------------------------------------
plateThickness = X
# WITHOUT HARDENING AND ULTIMATE STRAIN
AXIAL, MOMENT, STRAIN, CUR, NA, EI = P_M_INTERACTION(B, H, cover, plateThickness, Es, fy, As, ecuC, NUM=5000, STEEL_KIND=1)
# WITH HARDENING AND ULTIMATE STRAIN
AXIAL02, MOMENT02, STRAIN02, CUR02, NA02, EI02  = P_M_INTERACTION(B, H, cover, plateThickness, Es, fy, As, ecuC, NUM=5000, STEEL_KIND=2)
#%%------------------------------------------------------------------------------------------------  
plt.figure(1, figsize=(12, 8))
#plt.plot(MOMENT, AXIAL, color='black')
plt.scatter(MOMENT, AXIAL, color='blue', linewidth=2)
plt.scatter(MOMENT02, AXIAL02, color='cyan', linestyle='--', linewidth=2)
plt.title('P-M Interaction')
plt.ylabel('Axial Force [N]')
plt.xlabel('Bending Moment [N.mm]')
plt.legend(['Steel01: WITHOUT HARDENING AND ULTIMATE STRAIN', 'Hysteretic: WITH HARDENING AND ULTIMATE STRAIN'])
plt.grid()
plt.show()

plt.figure(2, figsize=(12, 8))
plt.plot(STRAIN, AXIAL, color='purple', linewidth=2)
plt.plot(STRAIN02, AXIAL02, color='magenta', linestyle='--', linewidth=2)
plt.title('AXIAL FORCE-STRAIN DIAGRAM')
plt.ylabel('Axial Force [N]')
plt.xlabel('Strain [mm/mm]')
plt.legend(['Steel01: WITHOUT HARDENING AND ULTIMATE STRAIN', 'Hysteretic: WITH HARDENING AND ULTIMATE STRAIN'])
plt.grid()
plt.show()

plt.figure(3, figsize=(12, 8))
plt.plot(CUR, MOMENT, color='red', linewidth=2)
plt.plot(CUR02, MOMENT02, color='orange', linestyle='--', linewidth=2)
plt.title('MOMENT-CURVATURE DIAGRAM')
plt.ylabel('Moment [kN.mm]')
plt.xlabel('Curvature [1/mm]')
plt.legend(['Steel01: WITHOUT HARDENING AND ULTIMATE STRAIN', 'Hysteretic: WITH HARDENING AND ULTIMATE STRAIN'])
plt.grid()
plt.show()

plt.figure(4, figsize=(12, 8))
plt.plot(NA, MOMENT, color='green', linewidth=2)
plt.plot(NA02, MOMENT02, color='lime', linestyle='--', linewidth=2)
plt.title('MOMENT-NEUTRAL AXIS DIAGRAM')
plt.ylabel('Moment [kN.mm]')
plt.xlabel('Neutral Axis Depth [mm]')
plt.legend(['Steel01: WITHOUT HARDENING AND ULTIMATE STRAIN', 'Hysteretic: WITH HARDENING AND ULTIMATE STRAIN'])
plt.grid()
plt.show()

plt.figure(5, figsize=(12, 8))
plt.plot(NA, EI, color='brown', linewidth=2)
plt.plot(NA02, EI02, color='gold', linestyle='--', linewidth=2)
plt.title('FLEXURAL RIGIDITY-NEUTRAL AXIS DIAGRAM')
plt.ylabel('Flexural Rigidity [N/mm²]')
plt.xlabel('Neutral Axis Depth [mm]')
plt.legend(['Steel01: WITHOUT HARDENING AND ULTIMATE STRAIN', 'Hysteretic: WITH HARDENING AND ULTIMATE STRAIN'])
plt.semilogy()
plt.grid()
plt.show()

plt.figure(6, figsize=(12, 8))
plt.plot(EI, np.abs(AXIAL), color='black', linewidth=2)
plt.plot(EI02, np.abs(AXIAL02), color='grey', linestyle='--', linewidth=2)
plt.title('FLEXURAL RIGIDITY-AXIAL RIGIDITY DIAGRAM')
plt.ylabel('Flexural Rigidity [N/mm²]')
plt.xlabel('Axial Rigidity [N]')
plt.legend(['Steel01: WITHOUT HARDENING AND ULTIMATE STRAIN', 'Hysteretic: WITH HARDENING AND ULTIMATE STRAIN'])
plt.semilogx()
plt.semilogy()
plt.grid()
plt.show()

#%%------------------------------------------------------------------------------   
# EXPORT DATA TO EXCEL
import pandas as pd
DATA_TOTAL = {
    'NEUTRAL-AXIS_WO': NA,  # WITHOUT HARDENING AND ULTIMATE STRAIN
    'NEUTRAL-AXIS_W': NA02, # WITH HARDENING AND ULTIMATE STRAIN
    'STRAIN_WO': STRAIN,
    'STRAIN_W': STRAIN02,
    'CURVATURE_WO': CUR,
    'CURVATURE_W': CUR02,
    'AXIAL_FORCE_WO': AXIAL,
    'AXIAL_FORCE_W': AXIAL02,
    'MOMENT_WO': MOMENT,
    'MOMENT_W': MOMENT02,
    'AXIAL_RIGIDITY_WO': np.abs(AXIAL),
    'AXIAL_RIGIDITY_W': np.abs(AXIAL02),
    'FLEXURAL_RIGIDITY_WO': EI,
    'FLEXURAL_RIGIDITY_W': EI02,
}
# Convert to DataFrame
results_df = pd.DataFrame(DATA_TOTAL)
# Export the DataFrame to an Excel file
results_df.to_excel('P-M_INTERACTION_COMPSITE_CONCRETE_OPTIMIZATION_RESULTS.xlsx', index=False) 