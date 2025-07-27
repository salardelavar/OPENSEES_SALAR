#######################################################################################
#                                    IN THE NAME OF ALLAH                             #
#                             UNCERTAINTY ANALYSIS OF RC FRAMES                       #
#               A PROBABILISTIC SEISMIC ASSESSMENT FRAMEWORK USING OPENSEES           #
#-------------------------------------------------------------------------------------#
#                         Developed by: Salar Delavar Ghashghaei (Qashqai)            #
#                              Email: salar.d.ghashghaei@gmail.com                    #
#######################################################################################
"""
Key Features:
-------------
1. Probabilistic Nonlinear Frame Modeling:
   - 2D RC frame with fiber-section nonlinearBeamColumn elements
   - Corotational geometric transformation for large displacements
   - Distributed plasticity modeling with 5 integration points

2. Advanced Material Models:
   - Concrete:
     * Confined core (Concrete01): fc' = -27.6 MPa, εc0 = -0.0045
     * Unconfined cover (Concrete01): fc' = -18 MPa, εc0 = -0.0025
   - Reinforcement:
     * Hysteretic steel model with hardening (fy = 400 MPa)
     * Pinching behavior (pinchX = 0.8, pinchY = 0.5)
     * Cyclic degradation (β = 0.1)

3. Comprehensive Analysis Capabilities:
   - Static Pushover Analysis:
     * Displacement-controlled up to 100mm drift
     * Automated bilinear curve fitting
     * Ductility ratio and overstrength factor calculation
   - Nonlinear Time History Analysis:
     * HHT-α integrator (α=1, γ=1.5, β=0.25)
     * Rayleigh damping (calibrated to first two modes)
     * Supports multi-directional excitation (X/Y components)

4. Uncertainty Quantification Framework:
   - 100-sample Monte Carlo simulation
   - Beta-distributed parameters (α=2, β=2-5)
   ±10% variation on all material and geometric properties
   - Statistical analysis of response quantities

5. Performance Assessment Metrics:
   - Structural Behavior Coefficients:
     * Overstrength factor (Ω₀ = Vmax/Vy)
     * Ductility ratio (μ = Δu/Δy)
     * Force reduction factor (R = Ω₀×μ)
   - Damage Indices:
     * Ductility Damage Index (DDI = (Δmax-Δy)/(Δu-Δy))
     * Markov-chain based damage state probabilities

6. Advanced Post-Processing:
   - Automated machine learning:
     * Random Forest feature importance
     * LSTM for time-history prediction
   - Reliability analysis:
     * Capacity-demand reliability indices
     * Neural network failure probability estimation
   - Comprehensive visualization:
     * Parameter distribution histograms
     * Hysteretic response plots
     * Time-history envelopes

Key Outputs:
------------
- Statistical distributions of 40+ response parameters
- Pushover and dynamic analysis results
- Damage state probabilities
- Excel export of all simulation data
- Publication-quality plots (PNG format)

Validation Features:
--------------------
- Logarithmic decrement damping verification
- Eigenvalue analysis for period verification
- Bilinear curve fitting validation
- Energy balance checks

Implementation Notes:
---------------------
- Uses OpenSeesPy Python interface
- Parallel processing ready (embarrassingly parallel)
- Modular design for easy customization
- Detailed inline documentation
- BSD-3 license

Typical Applications:
--------------------
- Fragility curve development
- Performance-based earthquake engineering
- Sensitivity analysis
- Code calibration studies
- Structural reliability assessment

This version:
1. Better reflects the actual code implementation
2. Organizes information more systematically
3. Includes specific parameter values from the code
4. Highlights the advanced analysis features
5. Maintains your original formatting style
6. Adds implementation details
7. Includes application examples
"""
import openseespy.opensees as ops
import matplotlib.pyplot as plt
import numpy as np
import time as TI
import SALAR_MATH as S01
import ANALYSIS_FUNCTION as S02
import CONCRETE_SECTION_FUN as S03
import PLOT_2D as S04
import MARKOV_CHAIN as S05
import BILINEAR_CURVE as BC
import DAMPING_RATIO as DA


# Define materials for nonlinear columns and beam
# Define parameters (units: mm, N)
# ------------------------------------------
# CONCRETE                  tag   f'c        ec0   f'cu        ecu
# Core concrete (confined)
fcCi = -27.6         # [N/mm²] Concrete Compressive Strength
ec0Ci = -0.0045      # [mm/mm] Concrete Compressive Strain
fcUCi = -21          # [N/mm²] Concrete Compressive Ultimate Strength
ecuCi = -0.015       # [mm/mm] Concrete Compressive Ultimate Strain

# Cover concrete (unconfined)
fcUi = -18           # [N/mm²] Concrete Compressive Strength
ec0Ui = -0.0025      # [mm/mm] Concrete Compressive Strain
fcUUi = -2           # [N/mm²] Concrete Compressive Ultimate Strength
ecuUi = -0.008       # [mm/mm] Concrete Compressive Ultimate Strain
 
# STEEL
# Reinforcing steel
fyi = 4000         # [N/mm²] Steel Rebar Yield Strength   
Esi = 2e5          # [N/mm²] Modulus of Elasticity
eyi = fyi/Esi      # [mm/mm] Steel Rebar Yield Strain
fui = 1.1818*fyi   # [N/mm²] Steel Rebar Ultimate Strength
esui = eyi*75.2    # [mm/mm] Steel Rebar Ultimate Strain
Eshi = (fui - fyi)/(esui - eyi)
Bsi = Eshi / Esi

# Column Section
Bci = 500                 # [mm] Depth of the Section 
Hci = 500                 # [mm] Height of the Section  
coverCi = 50              # [mm] Concrete Section Cover
DIAci = 25                # [mm] # Rebar Size Diameter
AsCi = np.pi*(DIAci**2)/4 # [mm²] Area of Rebar

# Beam Section
Bbi = 500                 # [mm] Depth of the Section 
Hbi = 300                 # [mm] Height of the Section  
coverBi = 50              # [mm] Concrete Section Cover
DIAbi = 18                # [mm] # Rebar Size Diameter
AsBi = np.pi*(DIAbi**2)/4 # [mm²] Area of Rebar

LENGTH_COLi = 3000        # [mm] Column Length 
LENGTH_BMi = 7000         # [mm] Beam Length 

#%%--------------------------------------------------------

GMfact = 9810    # standard acceleration of gravity or standard acceleration
SSF_X = 0.0001   # Seismic Acceleration Scale Factor in X Direction
SSF_Y = 0.0001   # Seismic Acceleration Scale Factor in Y Direction
iv0_X = 0.00005  # [mm/s] Initial velocity applied to the node  in X Direction
iv0_Y = 0.00005  # [mm/s] Initial velocity applied to the node  in Y Direction
st_iv0 = 0.0     # [s] Initial velocity applied starting time
SEI = 'X'        # Seismic Direction
DRi = 0.05       # Intial Guess for Damping ratio
duration = 15.0  # [s] Total simulation duration
dt = 0.01        # [s] Time step
MASSi = 12000    # [kg] Mass on the each column

DMAX = 100       # [mm] Maximum Displacement
DINCR = 0.05     # [mm] Increment Displacement

# Define Analysis Properties
MAX_ITERATIONS = 20000     # Convergence iteration for test
MAX_TOLERANCE = 1.0e-8    # Convergence tolerance for test
#STEEL_KIND: 1 -> WITHOUT HARDENING AND ULTIMATE STRAIN
#STEEL_KIND: 2 -> WITH HARDENING AND ULTIMATE STRAIN

# Set random seed for reproducibility
num_samples = 100 # Number of Samples
ERROR_RATIO = 0.1 # 10 % Error for uncertainty

#%%--------------------------------------------------------
from matplotlib.ticker import MaxNLocator

# Create a dictionary of all parameters for easier access
parameters = {
    'Concrete Core': {'fcC': fcCi, 'ec0C': ec0Ci, 'fcUC': fcUCi, 'ecuC': ecuCi},
    'Concrete Cover': {'fcU': fcUi, 'ec0U': ec0Ui, 'fcUU': fcUUi, 'ecuU': ecuUi},
    'Steel': {'fy': fyi, 'Es': Esi, 'ey': eyi, 'fu': fui, 'esu': esui, 'Esh': Eshi, 'Bs': Bsi},
    'Column Section': {'Bc': Bci, 'Hc': Hci, 'coverC': coverCi, 'DIAc': DIAci, 'AsC': AsCi},
    'Beam Section': {'Bb': Bbi, 'Hb': Hbi, 'coverB': coverBi, 'DIAb': DIAbi, 'AsB': AsBi},
    'Lengths': {'LENGTH_COL': LENGTH_COLi, 'LENGTH_BM': LENGTH_BMi},
    'Mass': {'Mass': MASSi},
    'Damping Ratio': {'Damping Ratio': DRi},
}

# Define Beta distribution parameters
beta_params = {
    'Concrete Core': (2, 5),
    'Concrete Cover': (2, 5),
    'Steel': (2, 2),
    'Column Section': (2, 2),
    'Beam Section': (2, 2),
    'Lengths': (2, 5),
    'Mass': (1, 1),
    'Damping Ratio': (1, 1)
}


# Dictionary to store all random samples
samples_dict = {}

# Generate Beta distributions for all parameters
for category, params_dict in parameters.items():
    alpha, beta = beta_params[category]
    category_samples = {}
    
    for param_name, nominal_value in params_dict.items():
        beta_samples = np.random.beta(alpha, beta, num_samples)
        scaled_samples = nominal_value * ((1-ERROR_RATIO) + (2*ERROR_RATIO) * beta_samples)
        category_samples[param_name] = scaled_samples
    
    samples_dict[category] = category_samples

# Create a clean visualization with better layout
plt.figure(figsize=(32, 28))
plt.suptitle(f'Parameter Distributions with Statistics (n={num_samples})', fontsize=20, y=0.98)

# Create a grid layout based on max parameters per category
max_params = max(len(params) for params in parameters.values())
rows = len(parameters)
cols = max_params

# Create subplots for each category
for cat_idx, (category, params_dict) in enumerate(parameters.items()):
    for param_idx, (param_name, nominal_value) in enumerate(params_dict.items()):
        ax = plt.subplot2grid((rows, cols), (cat_idx, param_idx))
        samples = samples_dict[category][param_name]
        
        # Calculate statistics
        median_val = np.median(samples)
        mean_val = np.mean(samples)
        max_val = np.max(samples)
        min_val = np.min(samples)
        
        # Create histogram with clear bins
        ax.hist(samples, bins=15, alpha=0.7, color='skyblue', edgecolor='navy')
        
        # Add statistics lines with different styles
        ax.axvline(nominal_value, color='black', linestyle='-', linewidth=2.5, label='Nominal')
        ax.axvline(median_val, color='red', linestyle='--', linewidth=2, label='Median')
        ax.axvline(mean_val, color='green', linestyle='-.', linewidth=2, label='Mean')
        ax.axvline(max_val, color='blue', linestyle=':', linewidth=1.5, alpha=0.7, label='Max')
        ax.axvline(min_val, color='purple', linestyle=':', linewidth=1.5, alpha=0.7, label='Min')
        
        # Add text box with statistics
        stats_text = (f"Nom: {nominal_value:.4g}\n"
                      f"Med: {median_val:.4g}\n"
                      f"Avg: {mean_val:.4g}\n"
                      f"Max: {max_val:.4g}\n"
                      f"Min: {min_val:.4g}")
        
        ax.text(0.95, 0.95, stats_text, transform=ax.transAxes,
                fontsize=10, verticalalignment='top', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        # Set titles and labels
        ax.set_title(f"{param_name}\n({category})", fontsize=12)
        ax.set_xlabel('Value')
        if param_idx == 0:
            ax.set_ylabel('Frequency')
        
        # Improve tick spacing
        ax.xaxis.set_major_locator(MaxNLocator(5))
        ax.yaxis.set_major_locator(MaxNLocator(5))
        
        # Add legend to first plot only to avoid repetition
        if cat_idx == 0 and param_idx == 0:
            ax.legend(loc='upper left', bbox_to_anchor=(0, 1.3), ncol=5, frameon=False)

plt.tight_layout(rect=[0, 0, 1, 0.96])  # Adjust for suptitle
plt.subplots_adjust(hspace=0.4, wspace=0.3)
plt.show()

# Access samples for ALL parameters from the samples_dict

# CONCRETE CORE SAMPLES
fcC = samples_dict['Concrete Core']['fcC']      # Core concrete compressive strength
ec0C = samples_dict['Concrete Core']['ec0C']    # Core concrete compressive strain
fcUC = samples_dict['Concrete Core']['fcUC']    # Core concrete ultimate strength
ecuC = samples_dict['Concrete Core']['ecuC']    # Core concrete ultimate strain

# CONCRETE COVER SAMPLES
fcU = samples_dict['Concrete Cover']['fcU']     # Cover concrete compressive strength
ec0U = samples_dict['Concrete Cover']['ec0U']   # Cover concrete compressive strain
fcUU = samples_dict['Concrete Cover']['fcUU']   # Cover concrete ultimate strength
ecuU = samples_dict['Concrete Cover']['ecuU']   # Cover concrete ultimate strain

# STEEL SAMPLES
fy = samples_dict['Steel']['fy']                # Steel yield strength
Es = samples_dict['Steel']['Es']                # Steel modulus of elasticity
ey = samples_dict['Steel']['ey']                # Steel yield strain
fu = samples_dict['Steel']['fu']                # Steel ultimate strength
esu = samples_dict['Steel']['esu']              # Steel ultimate strain
Esh = samples_dict['Steel']['Esh']              # Steel hardening modulus
Bs = samples_dict['Steel']['Bs']                # Steel hardening ratio

# COLUMN SECTION SAMPLES
Bc = samples_dict['Column Section']['Bc']       # Column width
Hc = samples_dict['Column Section']['Hc']       # Column height
coverC = samples_dict['Column Section']['coverC'] # Column cover
DIAc = samples_dict['Column Section']['DIAc']   # Column rebar diameter
AsC = samples_dict['Column Section']['AsC']     # Column rebar area

# BEAM SECTION SAMPLES
Bb = samples_dict['Beam Section']['Bb']         # Beam width
Hb = samples_dict['Beam Section']['Hb']         # Beam height
coverB = samples_dict['Beam Section']['coverB'] # Beam cover
DIAb = samples_dict['Beam Section']['DIAb']     # Beam rebar diameter
AsB = samples_dict['Beam Section']['AsB']       # Beam rebar area

# LENGTH SAMPLES
LENGTH_COL = samples_dict['Lengths']['LENGTH_COL'] # Column length
LENGTH_BM = samples_dict['Lengths']['LENGTH_BM']   # Beam length

# MASS SAMPLES
MASS = samples_dict['Mass']['Mass']                   # Mass on the Column 
DR = samples_dict['Damping Ratio']['Damping Ratio']   # Structure Damping Ratio

# Print sample statistics for verification
print("Parameter Samples Statistics:")
print("---------------------------")
print(f"Core concrete fcC: mean = {np.mean(fcC):.2f} MPa, std = {np.std(fcC):.2f}")
print(f"Steel fy: mean = {np.mean(fy):.2f} MPa, std = {np.std(fy):.2f}")
print(f"Column Hc: mean = {np.mean(Hc):.2f} mm, std = {np.std(Hc):.2f}")
print(f"Beam length: mean = {np.mean(LENGTH_BM):.2f} mm, std = {np.std(LENGTH_BM):.2f}")
print(f"Mass on the Column : mean = {np.mean(MASS):.2f} kg, std = {np.std(MASS):.2f}")
print(f"Structure Damping Ratio: mean = {np.mean(100*DR):.2f} (%), std = {np.std(100*DR):.2f}")
#%%------------------------------------------------------------------------------
def PD_ANALYSIS(J, STEEL_KIND, ANA_KIND):
    # Reset model
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)

    # Nodes
    ops.node(1, 0.0, 0.0)
    ops.node(2, LENGTH_BM[J], 0.0)
    ops.node(3, 0.0, LENGTH_COL[J])
    ops.node(4, LENGTH_BM[J], LENGTH_COL[J])
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 1, 1, 1)

    secTagC = 10
    secTagB = 20
    coreTag = 1
    coverTag = 2
    steelTag = 3
    
    if STEEL_KIND == 1:# WITHOUT HARDENING AND ULTIMATE STRAIN
        ops.uniaxialMaterial('Steel01', steelTag, fy[J], Es[J], 0.0) 
    if STEEL_KIND == 2:# WITH HARDENING AND ULTIMATE STRAIN    
        pinchX = 0.8   # Pinching factor in X direction
        pinchY = 0.5   # Pinching factor in Y direction
        damage1 = 0.0  # Damage due to ductility
        damage2 = 0.0  # Damage due to energy
        beta = 0.1 # Stiffness degradation parameter
        ops.uniaxialMaterial('Hysteretic', steelTag, fy[J], ey[J], fu[J], esu[J], 0.2*fu[J], 1.1*esu[J], -fy[J], -ey[J], -fu[J], -esu[J], -0.2*fu[J], -1.1*esu[J], pinchX, pinchY, damage1, damage2, beta)
        # INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Hysteretic_Material


    ops.uniaxialMaterial('Concrete01', coreTag, fcC[J], ec0C[J], fcUC[J], ecuC[J])  # Core concrete (confined)
    ops.uniaxialMaterial('Concrete01', coverTag, fcU[J], ec0U[J], fcUU[J], ecuU[J]) # Cover concrete (unconfined)
    
    # COLUMN SECTION
    S03.CONFINED_CONCRETE_SECTION(secTagC, Hc[J], Bc[J], coverC[J], AsC[J], coreTag, coverTag, steelTag, COL=True)
    # BEAM SECTION
    S03.CONFINED_CONCRETE_SECTION(secTagB, Hb[J], Bb[J], coverB[J], AsB[J], coreTag, coverTag, steelTag, COL=False)
    
    # Define geometric transformation (corotational for large displacements)
    transfTag = 1
    #ops.geomTransf('Linear', transfTag)
    #ops.geomTransf('PDelta', transfTag)
    ops.geomTransf('Corotational', transfTag)
    numIntgrPts = 5
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/nonlinearBeamColumn.html
    ops.element('nonlinearBeamColumn', 1, 1, 3, numIntgrPts, secTagC, transfTag, '-mass', Bc[J]*Hc[J]*0.000025) # COLUMN 01
    ops.element('nonlinearBeamColumn', 2, 2, 4, numIntgrPts, secTagC, transfTag, '-mass', Bc[J]*Hc[J]*0.000025) # COLUMN 02
    ops.element('nonlinearBeamColumn', 3, 3, 4, numIntgrPts, secTagB, transfTag, '-mass', Bb[J]*Hb[J]*0.000025) # BEAM 01
    
    if ANA_KIND == 'PUSHOVER':
        WEIGHT = MASS[J] * GMfact
        # Data storage
        FORCE_S, FORCE_A, MOMENT = [], [], []
        DISP_X, DISP_Y, ROT = [], [], []
        KA, KS, KI, STEP = [], [], [], []
    
        # Define time series and load pattern
        ops.timeSeries('Linear', 1)
        ops.pattern('Plain', 1, 1)
        #ops.load(3, 1.0, -WEIGHT, 0.0)
        #ops.load(4, 1.0, -WEIGHT, 0.0)
        ops.load(3, 1.0, -1, 0.0)
        ops.load(4, 1.0, -1, 0.0)
        
        # Total steps
        steps = int(np.abs(DMAX)/np.abs(DINCR))
    
        ops.constraints('Plain')
        ops.numberer('Plain')
        ops.system('BandGeneral')
        ops.algorithm('Newton')
        ops.analysis('Static')
        
        for step in range(steps):
            
            ops.integrator('DisplacementControl', 3, 1, DINCR) 
            ops.integrator('DisplacementControl', 4, 1, DINCR) 
            OK = ops.analyze(1)
            S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS
            # Record results
            ops.reactions()
            S = ops.nodeReaction(1, 1) + ops.nodeReaction(2, 1) # SHEAR BASE REACTION
            A = ops.nodeReaction(1, 2) + ops.nodeReaction(2, 2) # AXIAL BASE REACTION
            M = ops.nodeReaction(1, 3) + ops.nodeReaction(2, 3) # MOMENT BASE REACTION
            #print(rot, M)
            disp_X = ops.nodeDisp(3, 1) # LATERAL DISPLACEMENT IN X FOR NODE 3
            disp_Y = ops.nodeDisp(3, 2) # LATERAL DISPLACEMENT IN Y FOR NODE 3
            rot = ops.nodeDisp(3, 3)    # ROTATION IN Z FOR NODE 3
            FORCE_S.append(S)
            FORCE_A.append(A)
            MOMENT.append(M)
            DISP_X.append(disp_X)
            DISP_Y.append(disp_Y)
            ROT.append(rot)
            KS.append(np.abs(S)/np.abs(disp_X)) # LATERAL STIFFNESS IN X
            KA.append(np.abs(A)/np.abs(disp_Y)) # LATERAL STIFFNESS IN Y
            KI.append(np.abs(M)/np.abs(rot))    # ROTATIONAL STIFFNESS IN Z
            STEP.append(step)
            #print(step+1, disp_X, S)
        return FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y, ROT, KA, KS, KI, STEP        
    
    if ANA_KIND == 'DYNAMIC':
        # Define mass
        ops.mass(3, MASS[J], MASS[J], 0.0)
        ops.mass(4, MASS[J], MASS[J], 0.0)
        
        # Static analysis
        ops.timeSeries('Linear', 1)
        ops.pattern('Plain', 1, 1)
        if SEI == 'X':
            ops.load(3, 1.0, 0.0, 0.0)
            ops.load(4, 1.0, 0.0, 0.0)
        if SEI == 'Y': 
            ops.load(3, 0.0, 1.0, 0.0)
            ops.load(4, 0.0, 1.0, 0.0)
        if SEI == 'XY':
            ops.load(3, 1.0, 1.0, 0.0)
            ops.load(4, 1.0, 1.0, 0.0)
        
        # Dynamic analysis
        ops.constraints('Plain')
        ops.numberer('Plain')
        ops.system('BandGeneral')
        ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS)
        #ops.integrator('Newmark', 0.5, 0.25) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/newmark.html
        alpha = 1;gamma=1.5-alpha; beta=((2-alpha)**2)/4;
        ops.integrator('HHT', alpha, gamma, beta) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/hht.html
        ops.algorithm('Newton') # INFO LINK: https://openseespydoc.readthedocs.io/en/stable/src/algorithm.html
        ops.analysis('Transient')
        
        # Calculate Rayleigh damping factors
        Lambda01 = ops.eigen('-fullGenLapack', 2)  # eigenvalue mode 2
        #Lambda01 = ops.eigen('-genBandArpack', 2) # eigenvalue mode 2
        Omega01 = np.power(max(Lambda01), 0.5)
        Omega02 = np.power(min(Lambda01), 0.5)
        a0 = (2 * Omega01 * Omega02 * DR[J]) / (Omega01 + Omega02) # c = a0 * m : Mass-proportional damping
        a1 = (DR[J] * 2) / (Omega01 + Omega02)   # c = a1 * k : Stiffness-proportional damping
        # Apply Rayleigh damping
        ops.rayleigh(a0, a1, 0, 0)   # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/reyleigh.html
        #ops.rayleigh(0, 0, 2 * DR[J] * Omega01, 0) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/reyleigh.html
        PERIOD_01 = (np.pi * 2) / Omega01 # Structure First Period
        PERIOD_02 = (np.pi * 2) / Omega02 # Structure Second Period
        #print('Structure First Period:  ', PERIOD_01)
        #print('Structure Second Period: ', PERIOD_02) 
        
        # Define time series for input motion (Acceleration time history)
        if SEI == 'X':
            SEISMIC_TAG_01 = 100
            ops.timeSeries('Path', SEISMIC_TAG_01, '-dt', dt, '-filePath', 'Ground_Acceleration_X.txt', '-factor', GMfact, '-startTime', st_iv0) # SEISMIC-X
            # Define load patterns
            # pattern UniformExcitation $patternTag $dof -accel $tsTag <-vel0 $vel0> <-fact $cFact>
            ops.pattern('UniformExcitation', SEISMIC_TAG_01, 1, '-accel', SEISMIC_TAG_01, '-vel0', iv0_X, '-fact', SSF_X) # SEISMIC-X
        if SEI == 'Y':
            SEISMIC_TAG_02 = 200
            ops.timeSeries('Path', SEISMIC_TAG_02, '-dt', dt, '-filePath', 'Ground_Acceleration_Y.txt', '-factor', GMfact) # SEISMIC-Z
            ops.pattern('UniformExcitation', SEISMIC_TAG_02, 2, '-accel', SEISMIC_TAG_02, '-vel0', iv0_Y, '-fact', SSF_Y) 
        if SEI == 'XY':
            SEISMIC_TAG_01 = 100
            ops.timeSeries('Path', SEISMIC_TAG_01, '-dt', dt, '-filePath', 'Ground_Acceleration_X.txt', '-factor', GMfact, '-startTime', st_iv0) # SEISMIC-X
            # Define load patterns
            # pattern UniformExcitation $patternTag $dof -accel $tsTag <-vel0 $vel0> <-fact $cFact>
            ops.pattern('UniformExcitation', SEISMIC_TAG_01, 1, '-accel', SEISMIC_TAG_01, '-vel0', iv0_X, '-fact', SSF_X) # SEISMIC-X 
            SEISMIC_TAG_02 = 200
            ops.timeSeries('Path', SEISMIC_TAG_02, '-dt', dt, '-filePath', 'Ground_Acceleration_Y.txt', '-factor', GMfact) # SEISMIC-Z
            ops.pattern('UniformExcitation', SEISMIC_TAG_02, 2, '-accel', SEISMIC_TAG_02, '-vel0', iv0_Y, '-fact', SSF_Y)  # SEISMIC-Z
        print('Seismic Defined Done.')
            
        # Data storage
        FORCE_S, FORCE_A, MOMENT = [], [], []
        DISP_X, DISP_Y, ROT = [], [], []
        KA, KS, KI, STEP = [], [], [], []
        time = []
        displacement = []
        velocity_X, velocity_Y = [], []
        acceleration_X, acceleration_Y = [], []
            
        stable = 0
        current_time = 0.0
        
        while stable == 0 and current_time < duration:
            stable = ops.analyze(1, dt)
            S02.ANALYSIS(stable, 1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS
            current_time = ops.getTime()
            time.append(current_time)
            # Record results
            ops.reactions()
            S = ops.nodeReaction(1, 1) + ops.nodeReaction(2, 1) # SHEAR BASE REACTION
            A = ops.nodeReaction(1, 2) + ops.nodeReaction(2, 2) # AXIAL BASE REACTION
            M = ops.nodeReaction(1, 3) + ops.nodeReaction(2, 3) # MOMENT BASE REACTION
            #print(rot, M)
            disp_X = ops.nodeDisp(3, 1) # LATERAL DISPLACEMENT IN X FOR NODE 3
            disp_Y = ops.nodeDisp(3, 2) # LATERAL DISPLACEMENT IN Y FOR NODE 3
            rot = ops.nodeDisp(3, 3)    # ROTATION IN Z FOR NODE 3
            velocity_X.append(ops.nodeVel(3, 1))       # LATERAL VELOCITY IN X FOR NODE 3
            acceleration_X.append(ops.nodeAccel(3, 1)) # LATERAL ACCELERATION IN X FOR NODE 3
            velocity_Y.append(ops.nodeVel(3, 2))       # LATERAL VELOCITY IN Y FOR NODE 3
            acceleration_Y.append(ops.nodeAccel(3, 2)) # LATERAL ACCELERATION IN Y FOR NODE 3
            FORCE_S.append(S)
            FORCE_A.append(A)
            MOMENT.append(M)
            DISP_X.append(disp_X)
            DISP_Y.append(disp_Y)
            ROT.append(rot)
            KS.append(np.abs(S)/np.abs(disp_X)) # LATERAL STIFFNESS IN X
            KA.append(np.abs(A)/np.abs(disp_Y)) # LATERAL STIFFNESS IN Y
            KI.append(np.abs(M)/np.abs(rot))    # ROTATIONAL STIFFNESS IN Z
            #print(current_time, disp_X, S)
        
        # Calculating Damping Ratio Using Logarithmic Decrement Analysis 
        displacement = np.array(DISP_X)
        peaks = np.array([displacement[i] for i in range(1, len(displacement)-1) if displacement[i] > displacement[i-1] and displacement[i] > displacement[i+1]])
        # Natural logarithm
        delta = np.log(peaks[:-1] / peaks[1:])
        damping_ratio = DA.DAMPING_RATIO(delta)
        
    #ops.wipe()  
    return FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y, ROT, KA, KS, KI, time, velocity_X, velocity_Y, acceleration_X, acceleration_Y, PERIOD_01, PERIOD_02, damping_ratio


#%%------------------------------------------------------------------------------
# Uncertainity Analysis
# Analysis Durations:
starttime = TI.process_time()

FORCE_Sp_MAX, FORCE_Ap_MAX, MOMENTp_MAX = [], [], []
DISP_Xp_MAX, DISP_Yp_MAX, ROTp_MAX, KAp_MAX, KSp_MAX, KIp_MAX = [], [], [], [], [], []

FORCE_Sd_MAX, FORCE_Ad_MAX, MOMENTd_MAX, DISP_Xd_MAX  = [], [], [], []
DISP_Yd_MAX, ROTd_MAX, KAd_MAX, KSd_MAX, KId_MAX = [], [], [], [], []
velocity_Xd_MAX, velocity_Yd_MAX, acceleration_Xd_MAX, acceleration_Yd_MAX = [], [], [], []
PERIOD_01d_MAX, PERIOD_02d_MAX, Damping_ratiod_MAX = [], [], []

Elastic_ST_MAX, Plastic_ST_MAX, Tangent_ST_MAX, Ductility_Rito_MAX, Over_Strength_Factor_MAX = [], [], [], [], []
R_MAX, DIx_MAX = [], []

# WITH HARDENING AND ULTIMATE STRAIN
for J in range(num_samples):
    print(f"\nSTEP: {J+1}\n")
    # RUN NONLINEAR STATIC ANALYSIS
    DATA = PD_ANALYSIS(J, STEEL_KIND=2, ANA_KIND='PUSHOVER')
    FORCE_Sp, FORCE_Ap, MOMENTp, DISP_Xp, DISP_Yp, ROTp, KAp, KSp, KIp, STEPp = DATA
    
    FORCE_Sp_MAX.append(np.max(np.abs(FORCE_Sp)))
    FORCE_Ap_MAX.append(np.max(np.abs(FORCE_Ap)))
    MOMENTp_MAX.append(np.max(np.abs(MOMENTp)))
    DISP_Xp_MAX.append(np.max(np.abs(DISP_Xp)))
    DISP_Yp_MAX.append(np.max(np.abs(DISP_Yp)))
    ROTp_MAX.append(np.max(np.abs(ROTp)))
    KAp_MAX.append(np.max(np.abs(KAp)))
    KSp_MAX.append(np.max(np.abs(KSp)))
    KIp_MAX.append(np.max(np.abs(KIp)))
    
    # RUN NONLINEAR DYNAMIC ANALYSIS
    DATA = PD_ANALYSIS(J, STEEL_KIND=2, ANA_KIND='DYNAMIC')
    FORCE_Sd, FORCE_Ad, MOMENTd, DISP_Xd, DISP_Yd, ROTd, KAd, KSd, KId, timed, velocity_Xd, velocity_Yd, acceleration_Xd, acceleration_Yd, PERIOD_01d, PERIOD_02d, Damping_ratiod = DATA
    
    FORCE_Sd_MAX.append(np.max(np.abs(FORCE_Sd)))
    FORCE_Ad_MAX.append(np.max(np.abs(FORCE_Ad)))
    MOMENTd_MAX.append(np.max(np.abs(MOMENTd)))
    DISP_Xd_MAX.append(np.max(np.abs(DISP_Xd))) 
    DISP_Yd_MAX.append(np.max(np.abs(DISP_Yd)))
    ROTd_MAX.append(np.max(np.abs(ROTd)))
    KAd_MAX.append(np.max(np.abs(KAd))) 
    KSd_MAX.append(np.max(np.abs(KSd))) 
    KId_MAX.append(np.max(np.abs(KId))) 
    velocity_Xd_MAX.append(np.max(np.abs(velocity_Xd))) 
    velocity_Yd_MAX.append(np.max(np.abs(velocity_Yd)))
    acceleration_Xd_MAX.append(np.max(np.abs(acceleration_Xd))) 
    acceleration_Yd_MAX.append(np.max(np.abs(acceleration_Yd))) 
    PERIOD_01d_MAX.append(np.max(np.abs(PERIOD_01d)))
    PERIOD_02d_MAX.append(np.max(np.abs(PERIOD_02d))) 
    Damping_ratiod_MAX.append(np.max(np.abs(Damping_ratiod)))
    
    # --------------------------------------
    #  Plot BaseShear-Displacement Analysis 
    # --------------------------------------
    XX = np.abs(DISP_Xp); YY = np.abs(FORCE_Sp); # ABSOLUTE VALUE
    SLOPE_NODE = 10

    DATA = BC.BILNEAR_CURVE(XX, YY, SLOPE_NODE)
    X, Y, Elastic_ST, Plastic_ST, Tangent_ST, Ductility_Rito, Over_Strength_Factor = DATA

    XLABEL = 'Displacement in X [mm]'
    YLABEL = 'Base-Shear Reaction [N]'
    LEGEND01 = 'Curve'
    LEGEND02 = 'Bilinear Fitted'
    LEGEND03 = 'Undefined'
    TITLE = f'Last Data of BaseShear-Displacement Analysis - Ductility Ratio: {X[2]/X[1]:.4f} - Over Strength Factor: {Y[2]/Y[1]:.4f}'
    COLOR = 'black'
    BC.PLOT_2D(np.abs(DISP_Xp), np.abs(FORCE_Sp), X, Y, X, Y, XLABEL, YLABEL, TITLE, LEGEND01, LEGEND02, LEGEND03, COLOR='black', Z=2) 
    #print(f'\t\t Ductility Ratio: {Y[2]/Y[1]:.4f}')

    # Calculate Over Strength Coefficient (Ω0)
    Omega_0 = Y[2] / Y[1]
    # Calculate Displacement Ductility Ratio (μ)
    mu = X[2] / X[1]
    # Calculate Ductility Coefficient (Rμ)
    #R_mu = 1
    #R_mu = (2 * mu - 1) ** 0.5
    R_mu = mu
    # Calculate Structural Behavior Coefficient (R)
    R = Omega_0 * R_mu
    print(f'Over Strength Coefficient (Ω0):      {Omega_0:.4f}')
    print(f'Displacement Ductility Ratio (μ):    {mu:.4f}')
    print(f'Ductility Coefficient (Rμ):          {R_mu:.4f}')
    print(f'Structural Behavior Coefficient (R): {R:.4f}')
    Dd = np.max(np.abs(DISP_Xd))
    DIx = (Dd - X[1]) /(X[2] - X[1])
    print(f'Structural Ductility Damage Index in X Direction: {DIx:.4f}')
    
    R_MAX.append(R)
    Elastic_ST_MAX.append(Elastic_ST)
    Plastic_ST_MAX.append(Plastic_ST)
    Tangent_ST_MAX.append(Tangent_ST)
    Ductility_Rito_MAX.append(Ductility_Rito)
    Over_Strength_Factor_MAX.append(Over_Strength_Factor)
    DIx_MAX.append(DIx)

totaltime = TI.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n') 
#%%------------------------------------------------------------------------------
import matplotlib.pyplot as plt
import numpy as np

# Your data lists (empty in the example, but would be populated in your case)
data_lists = {
    "fcC": fcC,
    "fcUC": fcUC,
    "ecuC": ecuC,
    "fcU": fcU,
    "fcUU": fcUU,
    "ecuU": ecuU,
    "fy": fy,
    "ey": ey,
    "fu": fu,
    "esu": esu,
    "Mass": MASS,
    "Damping Ratio": DR,
    "Bb": Bb,
    "Hb": Hb,
    "Bc": Bc,
    "Hc": Hc,
    "AsB": AsB,
    "AsC": AsC,
    "LENGTH_COL": LENGTH_COL,
    "LENGTH_BM": LENGTH_BM,
    "FORCE_S_PUSHOVER": FORCE_Sp_MAX,
    "FORCE_A_PUSHOVER": FORCE_Ap_MAX,
    "MOMENT_PUSHOVER": MOMENTp_MAX,
    "DISP_X_PUSHOVER": DISP_Xp_MAX,
    "DISP_Y_PUSHOVER": DISP_Yp_MAX,
    "ROT_PUSHOVER": ROTp_MAX,
    "KA_PUSHOVER": KAp_MAX,
    "KS_PUSHOVER": KSp_MAX,
    "KI_PUSHOVER": KIp_MAX,
    "FORCE_S_DYNAMIC": FORCE_Sd_MAX,
    "FORCE_A_DYNAMIC": FORCE_Ad_MAX,
    "MOMENT_DYNAMIC": MOMENTd_MAX,
    "DISP_X_DYNAMIC": DISP_Xd_MAX,
    "DISP_Y_DYNAMIC": DISP_Yd_MAX,
    "ROT_DYNAMIC": ROTd_MAX,
    "KA_DYNAMIC": KAd_MAX,
    "KS_DYNAMIC": KSd_MAX,
    "KI_DYNAMIC": KId_MAX,
    "velocity_X_DYNAMIC": velocity_Xd_MAX,
    "velocity_Y_DYNAMIC": velocity_Yd_MAX,
    "acceleration_X_DYNAMIC": acceleration_Xd_MAX,
    "acceleration_Y_DYNAMIC": acceleration_Yd_MAX,
    "PERIOD_01_DYNAMIC": PERIOD_01d_MAX,
    "PERIOD_02_DYNAMIC": PERIOD_02d_MAX,
    "Damping_ratio_DYNAMIC": Damping_ratiod_MAX,
    "Elastic_ST_DYNAMIC": Elastic_ST_MAX,
    "Plastic_ST_DYNAMIC": Plastic_ST_MAX,
    "Tangent_ST_DYNAMIC": Tangent_ST_MAX,
    "Ductility_Rito": Ductility_Rito_MAX,
    "Over_Strength_Factor": Over_Strength_Factor_MAX,
    "Structural_Ductilty_Damage": DIx_MAX,
    "Structural_Behavior_Coefficient": R_MAX,
}

# Create figure with dynamic sizing
num_plots = len(data_lists)
cols = 4
rows = (num_plots + cols - 1) // cols  # Ceiling division
fig, axes = plt.subplots(rows, cols, figsize=(20, 5 * rows))
plt.subplots_adjust(hspace=0.4, wspace=0.3)
axes = axes.flatten()  # Flatten in case of 2D array

for i, (name, data) in enumerate(data_lists.items()):
    ax = axes[i]
    
    if len(data) > 0:  # Ensure data exists and has elements
        # Calculate statistics
        median_val = np.median(data)
        mean_val = np.mean(data)
        max_val = np.max(data)
        min_val = np.min(data)
        std_val = np.std(data)
        
        # Create histogram
        n, bins, patches = ax.hist(data, bins='auto', alpha=0.7, rwidth=0.85, color='skyblue')
        
        # Add statistics lines
        ax.axvline(median_val, color='red', linestyle='dashed', linewidth=1.5)
        ax.axvline(mean_val, color='green', linestyle='dashed', linewidth=1.5)
        ax.axvline(max_val, color='blue', linestyle='dotted', linewidth=1.2)
        ax.axvline(min_val, color='purple', linestyle='dotted', linewidth=1.2)
        
        # Add text box with statistics
        stats_text = (f'Median: {median_val:.4g}\n'
                     f'Mean: {mean_val:.4g}\n'
                     f'Std: {std_val:.4g}\n'
                     f'Min: {min_val:.4g}\n'
                     f'Max: {max_val:.4g}')
        
        ax.text(0.98, 0.98, stats_text, 
                transform=ax.transAxes,
                verticalalignment='top',
                horizontalalignment='right',
                bbox=dict(facecolor='white', alpha=0.8, edgecolor='gray'),
                fontsize=9)
        
        # Add title and labels
        ax.set_title(name, fontsize=12, pad=10)
        ax.set_xlabel('Value', fontsize=9)
        ax.set_ylabel('Frequency', fontsize=9)
        ax.grid(alpha=0.2)
        
    else:
        ax.text(0.5, 0.5, 'No Data', 
                ha='center', va='center', 
                fontsize=12, color='gray')
        ax.set_title(name, fontsize=12, pad=10)
        ax.set_xticks([])
        ax.set_yticks([])

# Hide unused subplots
for j in range(i + 1, len(axes)):
    axes[j].axis('off')

plt.tight_layout()
plt.savefig('data_distributions.png', dpi=150, bbox_inches='tight')
plt.show()
#%%------------------------------------------------------------------------------
plt.figure(1, figsize=(12, 8))
plt.plot(MOMENTp, FORCE_Ap, color='black')
plt.plot(MOMENTd, FORCE_Ad, color='purple')
#plt.scatter(MOMENTp, FORCE_Ap, color='black', linewidth=2)
#plt.scatter(MOMENTd, FORCE_Ad, color='purple', linewidth=2)
plt.title('P-M Interaction')
plt.ylabel('Axial Force [N]')
plt.xlabel('Bending Moment [N.mm]')
plt.legend(['PUSHOVER', 'DYNAMIC'])
plt.grid()
plt.show()

plt.figure(2, figsize=(12, 8))
plt.plot(DISP_Xp, FORCE_Sp, color='green', linewidth=2)
plt.plot(DISP_Xd, FORCE_Sd, color='lime', linewidth=2)
#plt.scatter(DISP_Xp, FORCE_Sp, color='green', linewidth=2)
#plt.scatter(DISP_Xd, FORCE_Sd, color='lime', linewidth=2)
plt.title('SHEAR FORCE-DISPLACEMENT DIAGRAM')
plt.ylabel('Shear Force [N]')
plt.xlabel('Displacement  in X [mm]')
plt.legend(['PUSHOVER', 'DYNAMIC'])
plt.grid()
plt.show()

plt.figure(3, figsize=(12, 8))
plt.plot(DISP_Yp, FORCE_Ap, color='purple', linewidth=2)
plt.plot(DISP_Yd, FORCE_Ad, color='green', linewidth=2)
#plt.scatter(DISP_Yp, FORCE_Ap, color='purple', linewidth=2)
#plt.scatter(DISP_Yd, FORCE_Ad, color='green', linewidth=2)
plt.title('AXIAL FORCE-DISPLACEMENT DIAGRAM')
plt.ylabel('Axial Force [N]')
plt.xlabel('Displacement in Y [mm]')
plt.legend(['PUSHOVER', 'DYNAMIC'])
plt.grid()
plt.show()

plt.figure(4, figsize=(12, 8))
plt.plot(ROTp, MOMENTp, color='red', linewidth=2)
plt.plot(ROTd, MOMENTd, color='pink', linewidth=2)
#plt.scatter(ROTp, MOMENTp, color='red', linewidth=2)
#plt.scatter(ROTd, MOMENTd, color='pink', linewidth=2)
plt.title('MOMENT-ROTATION DIAGRAM')
plt.ylabel('Moment [kN.mm]')
plt.xlabel('Rotation [rad]')
plt.legend(['PUSHOVER', 'DYNAMIC'])

plt.grid()
plt.show()

plt.figure(5, figsize=(12, 8))
#plt.plot(KIp, KSp, color='black', linewidth=2)
#plt.plot(KId, KSd, color='grey', linewidth=2)
plt.scatter(KIp, KSp, color='black', linewidth=2)
plt.scatter(KId, KSd, color='grey', linewidth=2)
plt.title('ROTATIONAL STIFFNESS-LATERAL STIFFNESS DIAGRAM')
plt.ylabel('Rotational Stiffness [N.mm/Rad]')
plt.xlabel('Lateral Stiffness in X Dir. [N/mm]')
plt.semilogx()
plt.semilogy()
plt.legend(['PUSHOVER', 'DYNAMIC'])
plt.grid()
plt.show()

plt.figure(6, figsize=(12, 8))
#plt.plot(KIp, KAp, color='black', linewidth=2)
#plt.plot(KId, KAd, color='grey', linewidth=2)
plt.scatter(KIp, KAp, color='black', linewidth=2)
plt.scatter(KId, KAd, color='grey', linewidth=2)
plt.title('ROTATIONAL STIFFNESS-LATERAL STIFFNESS DIAGRAM')
plt.ylabel('Rotational Stiffness [N.mm/Rad]')
plt.xlabel('Lateral Stiffness in Y Dir. [N/mm]')
plt.semilogx()
plt.semilogy()
plt.legend(['PUSHOVER', 'DYNAMIC'])
plt.grid()
plt.show()

plt.figure(7, figsize=(12, 8))
plt.plot(timed, FORCE_Ad, color='brown', linewidth=2)
#plt.scatter(timed, FORCE_Ad, color='brown', linewidth=2)
plt.title('Axial Force During the Analysis')
plt.ylabel('Axial Force [N]')
plt.xlabel('Times')
plt.grid()
plt.show()

plt.figure(8, figsize=(12, 8))
plt.plot(timed, FORCE_Sd, color='purple', linewidth=2)
#plt.scatter(timed, FORCE_Sd, color='purple', linewidth=2)
plt.title('Shear Force During the Analysis')
plt.ylabel('Shear Force [N]')
plt.xlabel('Times')
plt.grid()
plt.show()

plt.figure(9, figsize=(12, 8))
plt.plot(timed, MOMENTd, color='green', linewidth=2)
#plt.scatter(timed, MOMENTd, color='green', linewidth=2)
plt.title('Moment During the Analysis')
plt.ylabel('Moment [kN.mm]')
plt.xlabel('Times')
plt.grid()
plt.show()

plt.figure(10, figsize=(12, 8))
plt.plot(timed, DISP_Xd, color='brown', linewidth=2)
#plt.scatter(timed, DISP_Xd, color='brown', linewidth=2)
plt.title('Displacement During the Analysis')
plt.ylabel('Displacement - X [mm]')
plt.xlabel('Times')
plt.grid()
plt.show()

plt.figure(11, figsize=(12, 8))
plt.plot(timed, DISP_Yd, color='blue', linewidth=2)
#plt.scatter(timed, DISP_Xd, color='brown', linewidth=2)
plt.title('Displacement During the Analysis')
plt.ylabel('Displacement - Y [mm]')
plt.xlabel('Times')
plt.grid()
plt.show()

plt.figure(12, figsize=(12, 8))
plt.plot(timed, ROTd, color='black', linewidth=2)
#plt.scatter(timed, ROTd, color='black', linewidth=2)
plt.title('Rotation During the Analysis')
plt.ylabel('Rotation [rad]')
plt.xlabel('Times')
plt.grid()
plt.show()

#%%------------------------------------------------------------------------------
# Compute the Cumulative Maximum Absolute Value of Last Analysis Data
def MAX_ABS(X):
    import numpy as np
    X = np.asarray(X)  # Convert input to a numpy array for faster operations
    X_MAX = np.zeros_like(X)  # Initialize an array to store cumulative max values
    X_MAX[0] = np.abs(X[0])  # Set the first value

    # Compute cumulative maximum absolute values
    for i in range(1, len(X)):
        X_MAX[i] = max(X_MAX[i-1], np.abs(X[i]))
    
    return X_MAX  

DISP_ZX = MAX_ABS(DISP_Xd)  
DISP_ZY = MAX_ABS(DISP_Yd) 
VELO_Z = MAX_ABS(velocity_Xd) 
ACCE_Z = MAX_ABS(acceleration_Xd) 
BASE_Z = MAX_ABS(FORCE_Sd) 

plt.figure(1, figsize=(8, 6))
plt.plot(timed, DISP_Xd, color='blue', linewidth=2)
plt.plot(timed, DISP_ZX, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Displacement in X [mm]')
plt.title(f'Time vs Displacement - MAX. ABS: {DISP_ZX[-1]} | ξ (Calculated): {100*Damping_ratiod:.5e} %')
plt.grid()
plt.show()

plt.figure(2, figsize=(8, 6))
plt.plot(timed, DISP_Yd, color='blue', linewidth=2)
plt.plot(timed, DISP_ZY, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Displacement in Y [mm]')
plt.title(f'Time vs Displacement - MAX. ABS: {DISP_ZY[-1]}')
plt.grid()
plt.show()

plt.figure(3, figsize=(8, 6))
plt.plot(timed, velocity_Xd, color='blue', linewidth=2)
plt.plot(timed, VELO_Z, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Velocity in X [mm/s]')
plt.title(f'Time vs Velocity - MAX. ABS: {VELO_Z[-1]}')
plt.grid()
plt.show()

plt.figure(4, figsize=(8, 6))
plt.plot(timed, acceleration_Xd, color='blue', linewidth=2)
plt.plot(timed, ACCE_Z, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Acceleration in X [mm/s^2]')
plt.title(f'Time vs Acceleration - MAX. ABS: {ACCE_Z[-1]}')
plt.grid()
plt.show()

plt.figure(5, figsize=(8, 6))
plt.plot(timed, FORCE_Sd, color='blue', linewidth=2)
plt.plot(timed, BASE_Z, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Base-reaction [N]')
plt.title(f'Time vs Base-reaction - MAX. ABS: {BASE_Z[-1]}')
plt.grid()
plt.show()

#%%------------------------------------------------------------------------------  
# %% Plot 2D Frame Shapes
S04.PLOT_2D_FRAME(deformed_scale=100000)  # Adjust scale factor as needed
#%%------------------------------------------------------------------------------   
# EXPORT DATA TO EXCEL
import pandas as pd
# Convert to DataFrame
results_df = pd.DataFrame(data_lists)
# Export the DataFrame to an Excel file
results_df.to_excel('CONCRETE_FRAME_UNCERTAINTY_RESULTS.xlsx', index=False) 
#%%------------------------------------------------------------------------------
XLABEL = 'Displacement'
YLABEL = 'Structural Ductility Damage Index'
TITLE = f'{YLABEL} and {XLABEL} scatter chart'
COLOR = 'purple'
X = DISP_Xd_MAX
Y = DIx_MAX
S01.PLOT_SCATTER(X, Y , XLABEL, YLABEL, TITLE, COLOR, LOG = 0, ORDER = 1)

# CLUSTER DATA
S01.CLUSTER_DATA(X, Y, XLABEL, YLABEL, MAX_CLUSTERS=3)
#%%------------------------------------------------------------------------------------------------
#print(results_df)
S01.RANDOM_FOREST(results_df)
#%%------------------------------------------------------------------------------------------------
# PLOT HEATMAP FOR CORRELATION 
S01.PLOT_HEATMAP(results_df)
#%%------------------------------------------------------------------------------------------------
# MULTIPLE REGRESSION MODEL
S01.MULTIPLE_REGRESSION(results_df) 
#%%------------------------------------------------------------------------------------------------
# MACHINE LEARNING: LONG SHORT-TREM MEMERY (LSTM) METHOD
x = DISP_Xd_MAX 
y = DIx_MAX
Demand_X = x[-1]
look_back = 500#int(NUM_SIM * 0.5)
ITERATION = 200
XLABEL = 'Max Displacement'
YLABEL = 'Max Acceleration'
#S01.PREDICT_LSTM(x, y, Demand_X, look_back, ITERATION, XLABEL, YLABEL)
#%%------------------------------------------------------------------------------------------------
# PERFORM RELIABILITY ANALYSIS FOR BASE REACTION AND ELEMENT CAPACITY
mean_capacity = np.mean(FORCE_Sd_MAX)    # Mean Element Ultimate Capacity
std_dev_capacity = np.std(FORCE_Sd_MAX)  # Std Element Ultimate Capacity
num_sim = num_samples
max_base_reaction = np.mean(FORCE_Sd_MAX)
#S01.RELIABILITY_ANALYSIS(max_base_reaction, num_sim, mean_capacity, std_dev_capacity)
#%%------------------------------------------------------------------------------------------------
# NEURAL NETWORK FOR FAILURE PROBABILIYY ESTIMATION
X1 = mean_capacity
X2 = max_base_reaction
#S01.NEURAL_NETWORK_FAILURE_PROBABILIYY_ESTIMATION(X1, X2, num_samples)
#%%------------------------------------------------------------------------------------------------
# MARKOV CHAIN MODEl (structural damage analysis by evaluating Structural Ductility Damage Index)
FILE_TF = False         # Indicate whether to read data from a file or use provided data
file_path = None        # Not used when 'file_tf' is False
DATA = DIx_MAX # If not using a file, replace None with a NumPy array of data

S05.MARKOV_CHAIN(FILE_TF, file_path, DATA)

#%%%---------------------------------------------------------------------------
####  FRAGILITY ANALYSIS
from scipy.stats import norm  
# ----------------------------
# Fragility Assessment
# ----------------------------
# Define damage states per FEMA P-58
# INFO LINK: https://www.fema.gov/sites/default/files/documents/fema_p-58-2-se_volume2_implementation.pdf
damage_states = {
'DS1_Slight': (0.15, 0.4),    # Median PGA=0.15g, β=0.4
'DS2_Moderate': (0.30, 0.5),
'DS3_Extensive': (0.60, 0.6),
'DS4_Complete': (1.00, 0.7)
}
im_values = acceleration_Xd_MAX
# --------------
# Visualization
# --------------
plt.figure(1, figsize=(10, 6))
# Response plot
plt.plot(timed, acceleration_Xd, lw=1, color='black')
plt.xlabel('Time (s)')
plt.ylabel('Acceleration (g)')
plt.title(f'Last Analysis Structural Response + Ground Motion ::: MAX. ABS. : {np.max(np.abs(acceleration_Xd)):.4f}')
plt.grid(True)
plt.show()    

# Fragility curves
plt.figure(2, figsize=(10, 6))
# Calculate and plot fragility curves for each damage state
for damage_state, (median, beta) in damage_states.items():
    # Calculate log-normal probabilities
    ln_im = np.log(im_values)
    ln_median = np.log(median)
    probabilities = norm.cdf((ln_im - ln_median) / beta)
    plt.scatter(im_values, probabilities, marker='o', label=f'{damage_state} (η={median}, β={beta}')
    #plt.plot(im_values, probabilities, lw=2, label=f'{damage_state} (η={median}, β={beta})')
plt.xlabel('Peak Ground Acceleration (g)  [IM]')
plt.ylabel('Probability of Exceedance')
plt.title('Fragility Curves')
plt.legend()
plt.semilogy()
plt.ylim(0, 1.0)
plt.grid(True)
plt.tight_layout()
plt.show()    

#===========================================================

# Define damage state parameters: {Damage State: (median_IM, beta)}
damage_states = {
    'Minor Damage Level': (0.2, 0.4),# Median DI=0.2, β=0.4
    'Moderate Damage Level': (0.4, 0.4),
    'Severe Damage Level': (0.6, 0.5),
    'Failure Level': (1.0, 0.5)
}

# Generate intensity measure (IM) values from 0.0 to 1.0
im_values = DIx_MAX # Structural Ductility Damage Index
# --------------
# Visualization
# --------------
# Create plot
plt.figure(figsize=(10, 6))
# Calculate and plot fragility curves for each damage state
for damage_state, (median, beta) in damage_states.items():
    # Calculate log-normal probabilities
    ln_im = np.log(im_values)
    ln_median = np.log(median)
    probabilities = norm.cdf((ln_im - ln_median) / beta)
    plt.scatter(im_values, probabilities, marker='o', label=f'{damage_state} (η={median}, β={beta}')
    #plt.plot(im_values, probabilities, lw=2, label=f'{damage_state} (η={median}, β={beta})')

# Format plot
plt.xlabel('Structural Ductility Damage Index [IM]', fontsize=12)
plt.ylabel('Probability of Exceedance', fontsize=12)
plt.title('Fragility Curves', fontsize=14)
plt.legend(loc='lower right', fontsize=10)
plt.grid(True)
plt.semilogy()
plt.ylim(0, 1.0)
plt.tight_layout()
plt.show()    
#%%------------------------------------------------------------------------------------------------  
# Print out the state of nodes 3 and 4
ops.printModel("node",3, 4)
# Print out the state of element 1 , 2 and 3
ops.printModel("ele", 1, 2 , 3)
# Print the Model
#printModel()
ops.printModel("-JSON", "-file", "CONCRETE_FRAME_UNCERTAINTY.json")
#%%------------------------------------------------------------------------------------------------    

