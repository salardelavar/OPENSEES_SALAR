###########################################################################################################
#                    >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                     #
#  PARAMETRIC NONLINEAR SEISMIC DEMAND SPECTRA OF A 2-DOF (X & Y SPRINGS) SYSTEM UNDER 2 GROUND MOTIONS   #
#                                 WITH VARIABLE DUCTILITY CAPACITY                                        #
#---------------------------------------------------------------------------------------------------------#
#                EQUIVALENT VISCOUS DAMPING RATIO: xi_eq = 100 * E_d / (4 * pi * E_s)                     #
#---------------------------------------------------------------------------------------------------------#
#          ENERGY DISSIPATION CAPACITY INDEX = 100 * E_d(earthquake) / E_d(cyclic displacement)           #
#---------------------------------------------------------------------------------------------------------#
#                        BASED ON THE SCRIPT BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                        #
#                                   EMAIL: salar.d.ghashghaei@gmail.com                                   #
###########################################################################################################
"""
 This script builds a 2-DOF OpenSeesPy structure (-ndm 2 -ndf 2) with a fixed base node 
 and a mass node, joined by one zeroLength element that holds an independent hysteretic
 spring and viscous damper in X and in Y. The ultimate displacement of each spring is
 tied to ductility through DSU = DUCT * DY (with DY = FY/Ke), and the damping ratio is
 fixed at 3%. The outer loop sweeps ductility capacity from DUCT_MIN to DUCT_MAX. For
 each level it first runs a cyclic displacement protocol on the X and Y springs, whose
 amplitude equals DSU, to get the baseline dissipated energy (loop area by the shoelace 
 formula). The inner loop then sweeps the mass, which changes the elastic period, and for
 each mass it applies two ground motions together, one along X and one along Y, in a
 nonlinear Newmark time-history analysis. For each direction it records the peak
 displacement, velocity, acceleration and base reaction, the ductility damage index,
 the stiffness and the damping ratio from logarithmic decrement. It also computes the
 Energy Dissipation Capacity Index (earthquake dissipated energy as a percentage of the
 cyclic baseline) and the equivalent viscous damping ratio (dissipated energy divided by
 4π times the strain energy). These peaks are collected into response spectra against
 period, one curve per ductility level, with mean, median, quartile and standard-deviation
 bands. The script then exports everything to Excel and produces histograms, fragility
 curves, scatter and cluster plots, random forest, a correlation heatmap and multiple
 regression, all separately for X and Y. It ends with hysteresis loops, time histories
 and an X-Y displacement orbit for the last analysis.
MODEL (OpenSeesPy):
    ops.model('basic', '-ndm', 2, '-ndf', 2)
    Node 1 : fixed base (both DOFs fixed)
    Node 2 : mass node  (DOF 1 = X, DOF 2 = Y)
    One zeroLength element carrying
        - X spring (Hysteretic or Elastic)   + X viscous damper  -> direction 1
        - Y spring (Hysteretic or Elastic)   + Y viscous damper  -> direction 2
    The X and Y springs are independent (uncoupled) and can have different FY, Ke and OSF (see SPRING).

DUCTILITY CAPACITY:
    DSU = DUCT * DY,  DY = FY / Ke   ->  the ultimate displacement of both springs follows the ductility DUCT.
    DUCT is swept from DUCT_MIN to DUCT_MAX (NUM_DUCT + 1 levels); the damping ratio is fixed (DAMPING_RATIO).
    Because DSU changes, the cyclic baseline (Ed_CP) is recomputed for every ductility level.

LOADING:
    Two UniformExcitation patterns: ground motion X on DOF 1 and ground motion Y on DOF 2
    (files Ground_Acceleration_{III_X+1}.txt and Ground_Acceleration_{III_Y+1}.txt, values in g).

PER-DIRECTION OUTPUTS (X and Y):
    displacement, velocity, acceleration, base reaction, ductility damage index (DI), damping ratio,
    stiffness, energy dissipation capacity index (EDCI) and equivalent viscous damping ratio (EVDR).

EDCI   = 100 * Ed(earthquake) / Ed(cyclic displacement protocol)   (computed separately for X and Y)
EVDR   = 100 * Ed / (4*pi*Es),  Es = 0.5*Fmax*umax
DI     = (Dd - Dy) / (Du - Dy)
"""
#%%------------------------------------------------------
import openseespy.opensees as ops
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time as ti
import SALAR_MATH as S01
import ANALYSIS_FUNCTION as S02
import DAMPING_RATIO_FUN as S04
import OPENSEEES_HYSTERETICSM_FORCE_DISP_FUN as S08
import FRAGILITY_CURVE_FUN as S077
import EQULIVALENT_VISCOUS_DAMPING_RATIO_FUN as S055
import DAMAGE_INDEX_FUN as S066
#%%------------------------------------------------------
# Define parameters (units: m, N)
NUM_T = 30                                     # Total number of structural periods (mass steps) in each simulation
NUM_DUCT = 10                                  # Total number of ductility-capacity simulations (NUM_DUCT + 1 levels)

DUCT_MIN = 2.0                                 # Minimum ductility capacity  (DUCT = DSU / DY, must be > 1)
DUCT_MAX = 20.0                                # Maximum ductility capacity
DAMPING_RATIO = 0.03                           # Fixed Structural Damping Ratio
III_X = int(25)                                # Ground motion file index for X  -> Ground_Acceleration_{III_X+1}.txt
III_Y = int(26)                                # Ground motion file index for Y  -> Ground_Acceleration_{III_Y+1}.txt
SEI = 'XY'                                     # Seismic direction: 'X', 'Y' or 'XY'

TOTAL_MASS = 100000.0                          # [kg] Total Mass of Structure (same in X and Y)
GMfact = 9.81                                  # [m/s^2] standard acceleration of gravity (ground motion files are in g)
SSF_X = 50.0                                   # Seismic Acceleration Scale Factor in X Direction
SSF_Y = 50.0                                   # Seismic Acceleration Scale Factor in Y Direction
IV0_X = 0.0005                                 # [m/s] Initial velocity in X Direction
IV0_Y = 0.0005                                 # [m/s] Initial velocity in Y Direction

SHOW_LOOP_PLOTS = False                        # False: close the hysteresis figures created inside the sweep loops
TRACK_TANGENT_PERIOD = False                   # True : eigen-analysis at every step (slow) -> tangent period history

# Spring properties in each direction  (FY [N], OSF = FU/FY, KE [N/m]);  DSU = DUCT * DY is set by the sweep
SPRING = {
    'X': dict(FY=8500000.0, OSF=1.18, KE=450000000.0),
    'Y': dict(FY=8500000.0, OSF=1.18, KE=450000000.0),
}
#%%------------------------------------------------------
# Define Analysis Properties
MAX_ITERATIONS = 5000     # Convergence iteration for test
MAX_TOLERANCE = 1.0e-6    # Convergence tolerance for test
DIR = {'X': 1, 'Y': 2}    # DOF number of each direction
MAT = {'X': 1, 'Y': 2}    # Spring material tags
CMAT = {'X': 3, 'Y': 4}   # Damper material tags
#%%------------------------------------------------------
def SPRING_PARAMS(d, DUCT):
    P = SPRING[d]
    FY, Ke = P['FY'], P['KE']
    FU = P['OSF'] * FY                                # [N] Ultimate Force
    DY = FY / Ke                                      # [m] Yield Displacement
    DSU = DUCT * DY                                   # [m] Ultimate Displacement: DSU = DUCT * DY
    DP = [DY, DSU, 1.1*DSU, 1.25*DSU]
    FP = [FY, FU, 0.20*FU, 0.10*FU]
    DN = [-x for x in DP]
    FN = [-x for x in FP]
    return dict(FY=FY, FU=FU, Ke=Ke, DY=DY, DUCT=DUCT, DSU=DSU, DP=DP, FP=FP, DN=DN, FN=FN)

PAR = {}
def SET_DUCTILITY(DUCT):
    # Update the spring properties of both directions for the ductility DUCT (DSU = DUCT * DY).
    if DUCT <= 1.0:
        raise ValueError("DUCT must be greater than 1 (DSU > DY).")
    PAR['X'] = SPRING_PARAMS('X', DUCT)
    PAR['Y'] = SPRING_PARAMS('Y', DUCT)

SET_DUCTILITY(DUCT_MIN)
#%%------------------------------------------------------
# EVALUATION OF DISSIPATED ENERGY CAPACITY INDEX
def DISSIPATED_ENERGY_FUN_WITH_PLOT(displacement, base_shear, method, title="Hysteresis Curve"):
    if method == 1:
        """
        Compute dissipated energy using convex hull and plot the hysteresis curve
        with the outer hull area shaded.
        """
        from scipy.spatial import ConvexHull
        displacement = np.asarray(displacement)
        base_shear  = np.asarray(base_shear)

        if displacement.size != base_shear.size:
            raise ValueError("Displacement and base shear arrays must have equal lengths.")

        points = np.column_stack((displacement, base_shear))
        hull = ConvexHull(points)
        area = hull.volume   # 2D hull -> area

        fig, ax = plt.subplots(figsize=(7, 6))
        ax.plot(displacement, base_shear, 'k-', linewidth=1, label="Hysteresis Curve")
        hull_pts = points[hull.vertices]
        ax.plot(hull_pts[:, 0], hull_pts[:, 1], 'r--', lw=2, label="Convex Hull")
        ax.fill(hull_pts[:, 0], hull_pts[:, 1], color='red', alpha=0.25, label="Hull Area (Energy)")
        ax.set_title(f"{title} - (Convex Hull)")
        ax.set_xlabel("Displacement (m)")
        ax.set_ylabel("Base Shear (N)")
        ax.grid(True, linestyle='--', alpha=0.5)
        ax.legend()

    if method == 2:
        disp = np.asarray(displacement, dtype=float)
        shear = np.asarray(base_shear, dtype=float)

        if disp.size != shear.size:
            raise ValueError("Displacement and base shear arrays must have the same length.")
        if disp.size < 3:
            raise ValueError("At least 3 points are required to form a closed loop.")

        # Close the loop if not already closed (important for shoelace)
        if not (disp[0] == disp[-1] and shear[0] == shear[-1]):
            disp = np.append(disp, disp[0])
            shear = np.append(shear, shear[0])

        # Dissipated energy (E_d) via Shoelace formula
        x = disp
        y = shear
        area = 0.5 * np.abs(np.dot(x[:-1], y[1:]) - np.dot(y[:-1], x[1:]))

        fig, ax = plt.subplots(figsize=(7, 6))
        idx_max = np.argmax(np.abs(disp))
        ax.plot(disp, shear, 'k-', linewidth=1.2, label="Hysteresis Loop")
        ax.scatter(disp[idx_max], shear[idx_max], color='blue', s=80,
                   zorder=5, label=r"$(u_{\rm max}, F_{\rm max})$")
        ax.fill(disp, shear, color='red', alpha=0.25, label=f"E$_d$ = {area:.3f} N·m")
        ax.set_title(title)
        ax.set_xlabel("Displacement (m)")
        ax.set_ylabel("Base Shear (N)")
        ax.grid(True, linestyle='--', alpha=0.5)
        ax.legend(loc='lower right')
        fig.tight_layout()

    return area, fig
#%%------------------------------------------------------
def DECI_ZONE(DECI):
    # Print the damage zone of the Energy Dissipation Capacity Index [%] and cap it at 100.
    zones = ["ZONE 1: VERY MINOR DAMAGE", "ZONE 2: MINOR DAMAGE", "ZONE 3: MODERATE-LOW DAMAGE",
             "ZONE 4: MODERATE DAMAGE", "ZONE 5: MODERATE-HIGH DAMAGE", "ZONE 6: SEVERE-LOW DAMAGE",
             "ZONE 7: SEVERE-MEDIUM DAMAGE", "ZONE 8: SEVERE-HIGH DAMAGE", "ZONE 9: VERY SEVERE DAMAGE",
             "ZONE 10: FAILURE DAMAGE"]
    if DECI <= 100:
        print(f'\n\tDISSIPATED ENERGY CAPACITY INDEX: {DECI:.3f} [%]\n')
        if DECI <= 0:
            print("\n\tZONE 0: NO DAMAGE\n")
        else:
            print(f"\n\t{zones[min(int(np.ceil(DECI / 10.0)), 10) - 1]}\n")
    else:
        print('\n\tFOR EVALUATION OF DISSIPATED ENERGY CAPACITY INDEX:')
        print('\n\tCHECK THE CYCLIC DISPLACEMENT ANALYSIS AND IF IT IS POSSIBLE')
        print('\t\t\tINCREASE THE DISPLACEMENT.\n')
    return min(DECI, 100.0)
#%%------------------------------------------------------
### OPENSEES FUNCTIONS
def BUILD_MODEL(DR, MAT_TYPE, M):
    # 2-DOF model: node 1 fixed, node 2 free in X (DOF 1) and Y (DOF 2).
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 2)

    ops.node(1, 0.0, 0.0)               # Fixed base
    ops.node(2, 0.0, 0.0)               # Mass node
    ops.fix(1, 1, 1)                    # Fix both DOFs of the base
    ops.mass(2, M, M)                   # Mass in X and Y

    for d in ('X', 'Y'):
        P = PAR[d]
        if MAT_TYPE == 'INELASTIC':
            S08.OPENSEEES_HYSTERETICSM_FORCE_DISP_FUN(MAT[d], P['DP'], P['FP'], P['DN'], P['FN'], PLOT=False,
                                                      X_LABEL='Displacement (m)', Y_LABEL='Force [N]',
                                                      TITLE=f'FORCE-DISPLACEMENT CURVE - {d}')
        if MAT_TYPE == 'ELASTIC':
            ops.uniaxialMaterial('Elastic', MAT[d], P['Ke'], 0.0, 0.5*P['Ke'])  # Tension and compression are not same
        omega = np.sqrt(P['Ke'] / M)
        Cd = 2 * DR * omega * M                                  # [N·s/m] Damping coefficient
        ops.uniaxialMaterial('Viscous', CMAT[d], Cd, 1.0)        # alpha = 1.0 -> linear damper

    # One element: X spring + X damper on direction 1, Y spring + Y damper on direction 2
    ops.element('zeroLength', 1, 1, 2, '-mat', MAT['X'], MAT['Y'], CMAT['X'], CMAT['Y'], '-dir', 1, 2, 1, 2)


def SETUP_ANALYSIS(MAT_TYPE):
    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.system('BandGeneral')
    ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS)
    ops.algorithm('Newton')   # Newton for both: the 'Elastic' material above is bilinear (Eneg = 0.5*E)


def TANGENT_PERIODS():
    # Tangent periods (min, max) of the 2-DOF system from eigen analysis.
    lam = ops.eigen('-fullGenLapack', 2)
    T = [2*np.pi/np.sqrt(l) if l > 1e-12 else np.nan for l in lam]
    return np.fmin(T[0], T[1]), np.fmax(T[0], T[1])


def NEW_REC():
    return {'disp': [], 'velo': [], 'acc': [], 'reaction': [], 'stiffness': [], 'DI': []}


def RECORD(REC, MAT_TYPE, dynamic):
    # Store the response of both directions at the current step.
    ops.reactions()
    for d in ('X', 'Y'):
        dof, P, r = DIR[d], PAR[d], REC[d]
        r['reaction'].append(ops.nodeReaction(1, dof))           # BASE REACTION
        r['disp'].append(ops.nodeDisp(2, dof))                   # DISPLACEMENT NODE 2
        if dynamic:
            r['velo'].append(ops.nodeVel(2, dof))                # VELOCITY NODE 2
            r['acc'].append(ops.nodeAccel(2, dof))               # ACCELERATION NODE 2
        u = abs(r['disp'][-1])
        r['stiffness'].append(abs(r['reaction'][-1]) / u if u > 1e-12 else np.nan)
        if MAT_TYPE == 'INELASTIC':
            r['DI'].append(S066.DAMAGE_INDEX_FUN(r['disp'][-1], P['DY'], P['DSU']))   # DAMAGE INDEX
        else:
            r['DI'].append(0.0)
    if TRACK_TANGENT_PERIOD:
        Tmin, Tmax = TANGENT_PERIODS()
        REC['T_MIN'].append(Tmin)
        REC['T_MAX'].append(Tmax)


def SDOF_2D(DAMPING_RATIO, MAT_TYPE, TOTAL_MASS, ANAL_TYPE, III_X=0, III_Y=0, CTRL='X'):
    """
    ANAL_TYPE: 'PERIOD' | 'STATIC' | 'PUSHOVER' | 'CYCLIC_DISPLACEMENT' | 'SEISMIC'
    CTRL     : direction ('X' or 'Y') controlled in PUSHOVER / CYCLIC_DISPLACEMENT
    """
    BUILD_MODEL(DAMPING_RATIO, MAT_TYPE, TOTAL_MASS)
    REC = {'X': NEW_REC(), 'Y': NEW_REC(), 'TIME': [], 'T_MIN': [], 'T_MAX': []}
    center_node = 2

    if ANAL_TYPE == 'PERIOD':
        Tmin, Tmax = TANGENT_PERIODS()
        ops.modalProperties("-print", "-file", f"SALAR_ModalReport_{ANAL_TYPE}.txt", "-unorm")
        return Tmin, Tmax

    # Reference load for static analyses (unit load in each direction)
    if ANAL_TYPE in ('STATIC', 'PUSHOVER', 'CYCLIC_DISPLACEMENT'):
        ops.timeSeries('Linear', 1)
        ops.pattern('Plain', 1, 1)
        if ANAL_TYPE == 'STATIC':
            ops.load(center_node, 1.0, 1.0)
        else:
            ops.load(center_node, *( (1.0, 0.0) if CTRL == 'X' else (0.0, 1.0) ))

    if ANAL_TYPE == 'STATIC':
        SETUP_ANALYSIS(MAT_TYPE)
        ops.integrator('LoadControl', 1.0)
        ops.analysis('Static')
        OK = ops.analyze(1)
        if OK != 0:
            S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS)
        RECORD(REC, MAT_TYPE, dynamic=False)
        print('\n\nSTATIC ANALYSIS DONE.\n\n')
        return REC

    if ANAL_TYPE == 'PUSHOVER':
        SETUP_ANALYSIS(MAT_TYPE)
        dof = DIR[CTRL]
        DINCR = -0.001                       # [m] Displacement increment
        DMAX = -PAR[CTRL]['DSU']             # [m] Max. displacement
        ops.integrator('DisplacementControl', center_node, dof, DINCR)
        ops.analysis('Static')
        for step in range(int(np.abs(DMAX / DINCR))):
            OK = ops.analyze(1)
            if OK != 0:
                S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS)
            RECORD(REC, MAT_TYPE, dynamic=False)
        print('\n\nPUSHOVER ANALYSIS DONE.\n\n')
        return REC

    if ANAL_TYPE == 'CYCLIC_DISPLACEMENT':
        SETUP_ANALYSIS(MAT_TYPE)
        dof = DIR[CTRL]
        DMAX = -PAR[CTRL]['DSU']             # [m] Max. displacement
        n_points = 10000
        key_disp = np.array([
             0.0,
             0.1*DMAX,   -0.1*DMAX,
             0.5*DMAX,   -0.5*DMAX,
             0.8*DMAX,   -0.8*DMAX,
             DMAX,       -DMAX,
             0.1*DMAX,   -0.1*DMAX,
             0.5*DMAX,   -0.5*DMAX,
             0.8*DMAX,   -0.8*DMAX,
             DMAX,       -DMAX,
             0.0
        ])
        disp_protocol = np.interp(np.linspace(0, len(key_disp) - 1, n_points),
                                  np.arange(len(key_disp)), key_disp)
        ops.analysis('Static')
        for target_disp in disp_protocol:
            dU = target_disp - ops.nodeDisp(center_node, dof)
            ops.integrator('DisplacementControl', center_node, dof, dU)
            OK = ops.analyze(1)
            if OK != 0:
                S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS)
            RECORD(REC, MAT_TYPE, dynamic=False)
        print(f'\n\nCYCLIC DISPLACEMENT ANALYSIS ({CTRL}) DONE.\n\n')
        return REC

    if ANAL_TYPE == 'SEISMIC':
        duration = 20.0                      # [s] Analysis duration
        dt = 0.01                            # [s] Time step

        # Define time series + UniformExcitation patterns for the two ground motions
        if SEI in ('X', 'XY'):
            gm_X = np.loadtxt(f'Ground_Acceleration_{III_X+1}.txt')            # acceleration in g
            ops.timeSeries('Path', 100, '-dt', dt, '-values', *gm_X.tolist(), '-factor', GMfact)
            ops.pattern('UniformExcitation', 100, 1, '-accel', 100, '-vel0', IV0_X, '-fact', SSF_X)   # SEISMIC-X
        if SEI in ('Y', 'XY'):
            gm_Y = np.loadtxt(f'Ground_Acceleration_{III_Y+1}.txt')            # acceleration in g
            ops.timeSeries('Path', 200, '-dt', dt, '-values', *gm_Y.tolist(), '-factor', GMfact)
            ops.pattern('UniformExcitation', 200, 2, '-accel', 200, '-vel0', IV0_Y, '-fact', SSF_Y)   # SEISMIC-Y
        print('Seismic Defined Done.')

        SETUP_ANALYSIS(MAT_TYPE)
        ops.integrator('Newmark', 0.5, 0.25)
        ops.analysis('Transient')

        stable = 0
        current_time = 0.0
        while stable == 0 and current_time < duration:
            stable = ops.analyze(1, dt)
            if stable != 0:
                S02.ANALYSIS(stable, 1, MAX_TOLERANCE, MAX_ITERATIONS)
            current_time = ops.getTime()
            REC['TIME'].append(current_time)
            RECORD(REC, MAT_TYPE, dynamic=True)
        print('\n\nSEISMIC ANALYSIS DONE.\n\n')

        # Damping ratio by logarithmic decrement (each direction)
        REC['damping_ratio'] = {d: S04.DAMPING_RATIO(REC[d]['disp']) for d in ('X', 'Y')}
        return REC
#%%------------------------------------------------------
# Analysis Durations:
starttime = ti.process_time()

# Collect Data  (each key -> {'X': [...], 'Y': [...]})
DATA = {k: {'X': [], 'Y': []} for k in range(1, 10)}
# 1: DISP, 2: VELOCITY, 3: ACCELERATION, 4: REACTION, 5: DI, 6: DAMPING RATIO, 7: STIFFNESS,
# 8: DISSIPATED ENERGY CAPACITY INDEX, 9: EQUIVALENT VISCOUS DAMPING RATIO

# Ductility levels and elastic periods of each direction for the mass sweep: T = 2*pi*sqrt(m/K)
# NOTE: mass step is written exactly as in the original script: m = TOTAL_MASS * (j + 1/NUM_T)
DUCTS = [DUCT_MIN + (DUCT_MAX - DUCT_MIN) * (i / NUM_DUCT) for i in range(NUM_DUCT + 1)]
MASSES = [TOTAL_MASS * (j + 1 / NUM_T) for j in range(NUM_T)]
T_EL = {d: np.array([2*np.pi*np.sqrt(m / PAR[d]['Ke']) for m in MASSES]) for d in ('X', 'Y')}

#%% RESPONSE SPECTRUM SWEEP OVER DUCTILITY CAPACITY
for i, DUCT in enumerate(DUCTS):

    SET_DUCTILITY(DUCT)                                  # DSU = DUCT * DY (X and Y)
    print(f"\n DUCTILITY {DUCT:.3f}:  DSU_X = {PAR['X']['DSU']:.4f} m | DSU_Y = {PAR['Y']['DSU']:.4f} m\n")

    #%% BASELINE: CYCLIC DISPLACEMENT ANALYSIS FOR EACH DIRECTION (DMAX = DSU -> depends on ductility)
    ED_CP = {}
    for d in ('X', 'Y'):
        REC_CP = SDOF_2D(DAMPING_RATIO, 'INELASTIC', TOTAL_MASS, 'CYCLIC_DISPLACEMENT', CTRL=d)
        ED_CP[d], fig_CP = DISSIPATED_ENERGY_FUN_WITH_PLOT(REC_CP[d]['disp'], REC_CP[d]['reaction'], method=2,
                                                           title=f"Cyclic Loading ({d}), DUCT = {DUCT:.2f} – Dissipated Energy")
        if SHOW_LOOP_PLOTS:
            fig_CP.show()
        else:
            plt.close(fig_CP)
        print(f"Dissipated Energy from Cyclic Displacement {d} = {ED_CP[d]:.2f} N·m")

    # Lists of maxima for this ductility
    mx = {d: {k: [] for k in ('time', 'displacement', 'velocity', 'acceleration', 'base_reaction',
                              'DI', 'K', 'dr', 'edci', 'zeta')} for d in ('X', 'Y')}

    for j, m in enumerate(MASSES):

        # SEISMIC ANALYSIS (DYNAMIC TIME-HISTORY ANALYSIS) - X AND Y GROUND MOTIONS TOGETHER
        REC_SEI = SDOF_2D(DAMPING_RATIO, 'INELASTIC', m, 'SEISMIC', III_X, III_Y)

        for d in ('X', 'Y'):
            R = REC_SEI[d]

            #%% START - Energy Dissipation Capacity Index
            Ed_SEI, fig_SEI = DISSIPATED_ENERGY_FUN_WITH_PLOT(
                R['disp'], R['reaction'], method=2,
                title=f"Earthquake Response ({d}) – Dissipated Energy")
            if SHOW_LOOP_PLOTS:
                fig_SEI.show()
            else:
                plt.close(fig_SEI)
            print(f"Dissipated Energy from Earthquake {d} = {Ed_SEI:.2f} N·m")
            DECI = DECI_ZONE(100 * Ed_SEI / ED_CP[d])
            #%% END - Energy Dissipation Capacity Index

            #%% START - Equivalent viscous damping ratio
            zeta = S055.EQULIVALENT_VISCOUS_DAMPING_RATIO_FUN(
                R['disp'], R['reaction'], 1, "Displacement [m]", "Base Reaction [N]",
                f"Equivalent viscous damping ratio ({d}) - Incremental Dynamic Analysis Hysteresis")
            print(f"Equivalent viscous damping ratio {d} = {zeta:.4f}")
            #%% END - Equivalent viscous damping ratio

            # Calculate and store the max absolute values
            mx[d]['time'].append(np.max(np.abs(REC_SEI['TIME'])))
            mx[d]['displacement'].append(np.max(np.abs(R['disp'])))
            mx[d]['velocity'].append(np.max(np.abs(R['velo'])))
            mx[d]['acceleration'].append(np.max(np.abs(R['acc'])))
            mx[d]['base_reaction'].append(np.max(np.abs(R['reaction'])))
            mx[d]['DI'].append(np.max(np.abs(R['DI'])))
            mx[d]['K'].append(np.nanmax(np.abs(R['stiffness'])))
            mx[d]['dr'].append(REC_SEI['damping_ratio'][d])
            mx[d]['edci'].append(DECI)
            mx[d]['zeta'].append(zeta)

        print(f'SEISMIC (X: {III_X+1}, Y: {III_Y+1}) - DUCT {i + 1} - STEP {j + 1} DONE')

    # Store Data
    for d in ('X', 'Y'):
        for key, name in zip(range(1, 10), ('displacement', 'velocity', 'acceleration', 'base_reaction',
                                            'DI', 'dr', 'K', 'edci', 'zeta')):
            DATA[key][d].append(mx[d][name])
    print(f'\n DUCTILITY LEVEL {i + 1} DONE \n')
else:
    print('\n\n Analysis completed successfully')

totaltime = ti.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n')

# Variables of the last analysis (last ductility / last mass) are used for the post-processing below
max_X, max_Y = mx['X'], mx['Y']
#%%------------------------------------------------------
# OUTPUT DATA FROM RESPONSE SPECTRUM ANALYSIS TO EXCEL FILE
labels = {1: "DISP", 2: "VELO", 3: "ACC", 4: "BASE_REACTION", 5: "DAMAGE_INDEX", 6: "DAMPING_RATIO",
          7: "STIFFNESS", 8: "EDCI", 9: "EVDR"}

with pd.ExcelWriter("INELASTIC_RESPONSE_SPECTRUM_SEISMIC_2DOF_DUCTILITY_RESULTS.xlsx", engine="openpyxl") as writer:
    for key, value in DATA.items():
        for d in ('X', 'Y'):
            pd.DataFrame(value[d]).to_excel(writer, sheet_name=f"{labels[key]}_{d}"[:31], index=False)
    pd.DataFrame({'T_X': T_EL['X'], 'T_Y': T_EL['Y']}).to_excel(writer, sheet_name="PERIOD", index=False)
    pd.DataFrame({'DUCT': DUCTS, 'DSU_X': [d * PAR['X']['DY'] for d in DUCTS],
                  'DSU_Y': [d * PAR['Y']['DY'] for d in DUCTS]}).to_excel(writer, sheet_name="DUCTILITY", index=False)

print("Excel file saved: INELASTIC_RESPONSE_SPECTRUM_SEISMIC_2DOF_DUCTILITY_RESULTS.xlsx")
#%%------------------------------------------------------
# PLOT THE RESPONSE SPECTRUMS (X AND Y)
def PLOT_2D(COUNT, d, XLABEL, YLABEL, TITLE):
    T = T_EL[d]
    arr = np.array(DATA[COUNT][d])                  # Shape: (NUM_DUCT+1, NUM_T)
    plt.figure(figsize=(12, 10))
    for j in range(arr.shape[0]):
        plt.plot(T, arr[j], alpha=0.4, label=f'DUCT {DUCTS[j]:.1f} | Max: {np.nanmax(np.abs(arr[j])):.4e}')

    mean_curve   = np.nanmean(arr, axis=0)
    median_curve = np.nanmedian(arr, axis=0)
    std_curve    = np.nanstd(arr, axis=0)
    q1_curve = np.nanpercentile(arr, 25, axis=0)
    q3_curve = np.nanpercentile(arr, 75, axis=0)

    plt.plot(T, mean_curve, color='navy', linestyle='-.', linewidth=3, label='Mean')
    plt.plot(T, median_curve, color='red', linestyle='--', linewidth=2, label='Median')
    plt.plot(T, q1_curve, color='green', linestyle='-.', linewidth=2, label='Q1 (25%)')
    plt.plot(T, q3_curve, color='orange', linestyle='-.', linewidth=2, label='Q3 (75%)')
    plt.fill_between(T, mean_curve - std_curve, mean_curve + std_curve,
                     color='gray', alpha=0.25, label='Mean ± Std')
    plt.xlabel(XLABEL)
    plt.ylabel(YLABEL)
    plt.title(f'{TITLE} - {d} Direction')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

SPECTRA = [
    (1, "Max Displacement (PGD) [m]",               "Displacement Response Spectrum"),
    (2, "Max Velocity (PGV) [m/s]",                 "Velocity Response Spectrum"),
    (3, "Max Acceleration (PGA) [m/s²]",            "Acceleration Response Spectrum"),
    (4, "Max Base Reaction (N)",                    "Base Reaction Response Spectrum"),
    (5, "Damage Index (DI) [%]",                    "Ductility Damage Index Spectrum"),
    (6, "Damping Ratio [%]",                        "Damping Ratio Spectrum"),
    (7, "Structural Stiffness Spectrum [N/m]",      "Structural Stiffness Spectrum"),
    (8, "Energy Dissipation Capacity Index [%]",    "Energy Dissipation Capacity Index Spectrum"),
    (9, "Equivalent viscous damping ratio [%]",     "Equivalent viscous damping ratio Spectrum"),
]
for d in ('X', 'Y'):
    for COUNT, YLABEL, TITLE in SPECTRA:
        PLOT_2D(COUNT, d, f"Structural period {d} [s]", YLABEL, TITLE)
#%%------------------------------------------------------
# Print the last results
for d in ('X', 'Y'):
    M = max_X if d == 'X' else max_Y
    print(f"\nMaximum Absolute Values - Direction {d}:")
    print("Duration:", M['time'][-1])
    print("Displacement:", M['displacement'][-1])
    print("Velocity:", M['velocity'][-1])
    print("Acceleration:", M['acceleration'][-1])
    print("Base Reaction:", M['base_reaction'][-1])
    print("Ductility Damage Index:", M['DI'][-1])
    print("Energy Dissipation Capacity Index:", M['edci'][-1])
    print("Equivalent viscous damping ratio:", M['zeta'][-1])
    print("Period:", T_EL[d][-1])
#%%------------------------------------------------------
for d in ('X', 'Y'):
    M = max_X if d == 'X' else max_Y
    S01.HISROGRAM_BOXPLOT(M['displacement'], HISTO_COLOR='blue', LABEL=f'Displacement {d}')
    S01.HISROGRAM_BOXPLOT(M['velocity'], HISTO_COLOR='purple', LABEL=f'Velocity {d}')
    S01.HISROGRAM_BOXPLOT(M['acceleration'], HISTO_COLOR='green', LABEL=f'Acceleration {d}')
    S01.HISROGRAM_BOXPLOT(M['base_reaction'], HISTO_COLOR='gold', LABEL=f'Base Reaction {d}')
    S01.HISROGRAM_BOXPLOT(M['DI'], HISTO_COLOR='pink', LABEL=f'Ductility Damage Index {d}')
    S01.HISROGRAM_BOXPLOT(T_EL[d], HISTO_COLOR='lime', LABEL=f'Structure Period {d}')
    S01.HISROGRAM_BOXPLOT(M['edci'], HISTO_COLOR='cyan', LABEL=f'Energy Dissipation Capacity Index {d}')
    S01.HISROGRAM_BOXPLOT(M['zeta'], HISTO_COLOR='brown', LABEL=f'Equivalent viscous damping ratio {d}')
#%%------------------------------------------------------
####  FRAGILITY ANALYSIS
# Define damage states per FEMA P-58
# INFO LINK: https://www.fema.gov/sites/default/files/documents/fema_p-58-2-se_volume2_implementation.pdf
damage_states_PGA = {
    'DS1_Slight': (0.15, 0.05),
    'DS2_Moderate': (0.30, 0.05),
    'DS3_Extensive': (0.60, 0.05),
    'DS4_Complete': (1.00, 0.05)
}
# Define damage state parameters: {Damage State: (median_IM, beta)}
damage_states = {
    'Minor Damage Level': (20.0, 0.4),
    'Moderate Damage Level': (40.0, 0.4),
    'Severe Damage Level': (60.0, 0.4),
    'Failure Level': (100.0, 0.4)
}

# Last-analysis acceleration time history of both directions
REC_LAST = REC_SEI
plt.figure(figsize=(10, 6))
plt.plot(REC_LAST['TIME'], REC_LAST['X']['acc'], lw=1, color='black', label='X')
plt.plot(REC_LAST['TIME'], REC_LAST['Y']['acc'], lw=1, color='red', label='Y')
plt.xlabel('Time (s)')
plt.ylabel('Acceleration (m/s²)')
plt.title(f"Last Analysis Response - MAX. ABS. X: {np.max(np.abs(REC_LAST['X']['acc'])):.4f} | "
          f"Y: {np.max(np.abs(REC_LAST['Y']['acc'])):.4f}")
plt.legend()
plt.grid(True)
plt.show()

for d in ('X', 'Y'):
    M = max_X if d == 'X' else max_Y
    S077.FRAGILITY_CURVE(M['acceleration'], damage_states_PGA,
                         f'Peak Ground Acceleration {d} (g)  [IM]', SEMILOGY=False, PLOT_KIND=True)
    S077.FRAGILITY_CURVE(M['edci'], damage_states,
                         f'Energy Dissipation Capacity Index {d} (%)  [IM]', SEMILOGY=False, PLOT_KIND=True)
    S077.FRAGILITY_CURVE(M['DI'], damage_states,
                         f'Ductility Damage Index {d} (%)  [IM]', SEMILOGY=False, PLOT_KIND=True)
    S077.FRAGILITY_CURVE(M['zeta'], damage_states,
                         f'Equivalent viscous damping ratio {d} (%)  [IM]', SEMILOGY=False, PLOT_KIND=True)
#%%------------------------------------------------------
# HYSTERESIS LOOPS AND DISPLACEMENT ORBIT OF THE LAST DYNAMIC ANALYSIS
fig, axs = plt.subplots(1, 3, figsize=(18, 6))
for ax, d in zip(axs[:2], ('X', 'Y')):
    ax.plot(REC_LAST[d]['disp'], REC_LAST[d]['reaction'], color='black')
    ax.set_xlabel(f"Displacement {d} [m]")
    ax.set_ylabel(f"Base Reaction {d} [N]")
    ax.set_title(f"Hysteresis Loop - {d}")
    ax.grid(True)
axs[2].plot(REC_LAST['X']['disp'], REC_LAST['Y']['disp'], color='navy')
axs[2].set_xlabel("Displacement X [m]")
axs[2].set_ylabel("Displacement Y [m]")
axs[2].set_title("Displacement Orbit (X-Y)")
axs[2].set_aspect('equal', adjustable='datalim')
axs[2].grid(True)
plt.tight_layout()
plt.show()
#%%------------------------------------------------------
# Cumulative Maximum Absolute Value of Last Analysis Data
def MAX_ABS(X):
    return np.maximum.accumulate(np.abs(np.asarray(X)))

TIME_SEI = REC_LAST['TIME']
for d in ('X', 'Y'):
    R = REC_LAST[d]
    for name, ylabel, key in (('Displacement', 'Displacement [m]', 'disp'),
                              ('Velocity', 'Velocity [m/s]', 'velo'),
                              ('Acceleration', 'Acceleration [m/s^2]', 'acc'),
                              ('Base-reaction', 'Base-reaction [N]', 'reaction')):
        Z = MAX_ABS(R[key])
        plt.figure(figsize=(8, 6))
        plt.plot(TIME_SEI, R[key], color='blue', linewidth=2)
        plt.plot(TIME_SEI, Z, color='red', linewidth=2)
        plt.xlabel('Time [s]')
        plt.ylabel(ylabel)
        plt.title(f'Time vs {name} ({d}) - MAX. ABS: {Z[-1]}')
        plt.grid()
        plt.show()
#%%------------------------------------------------------
# EXPORT DATA TO EXCEL
with pd.ExcelWriter('MDOF2D_DUCTILITY_RESPONSE_SPECTRUM_EDCI_EVDR_DATA_RESULTS.xlsx', engine='openpyxl') as writer:
    for d in ('X', 'Y'):
        M = max_X if d == 'X' else max_Y
        pd.DataFrame({
            'Max_displacement': M['displacement'],
            'Max_velocity': M['velocity'],
            'Max_acceleration': M['acceleration'],
            'Max_Base_Reaction': M['base_reaction'],
            'Ductility_Damage_Index': M['DI'],
            'Period': T_EL[d],
            'Energy Dissipation Capacity Index': M['edci'],
            'Equivalent viscous damping ratio': M['zeta'],
        }).to_excel(writer, sheet_name=f'RESULTS_{d}', index=False)
#%%------------------------------------------------------
# SCATTER + CLUSTER PLOTS (X AND Y)
for d in ('X', 'Y'):
    M = max_X if d == 'X' else max_Y
    for XL, YL, XD, YD, COLOR in (
        ('Displacement', 'Base Reaction', M['displacement'], M['base_reaction'], 'orange'),
        ('Velocity', 'Base Reaction', M['velocity'], M['base_reaction'], 'cyan'),
        ('Acceleration', 'Base Reaction', M['acceleration'], M['base_reaction'], 'lime'),
        ('Displacement', 'Structural Ductility Damage Index (%)', M['displacement'], M['DI'], 'purple')):
        TITLE = f'{YL} and {XL} scatter chart - {d}'
        S01.PLOT_SCATTER(XD, YD, XL, YL, TITLE, COLOR, LOG=0, ORDER=1)
        S01.CLUSTER_DATA(XD, YD, XL, YL, MAX_CLUSTERS=3)
#%%------------------------------------------------------
# PLOT THE TIME-HISTORY (X AND Y)
for d in ('X', 'Y'):
    R = REC_LAST[d]
    S01.PLOT_TIME_HISTORY(TIME_SEI, R['disp'], R['velo'], R['acc'], R['reaction'])
#%%------------------------------------------------------
# RANDOM FOREST / HEATMAP / MULTIPLE REGRESSION (X AND Y)
for d in ('X', 'Y'):
    M = max_X if d == 'X' else max_Y
    df = pd.DataFrame({
        'Max_displacement': M['displacement'],
        'Max_velocity': M['velocity'],
        'Max_acceleration': M['acceleration'],
        'Max_Base_Reaction': M['base_reaction'],
        'Ductility_Damage_Index': M['DI'],
        'Energy_Dissipation_Capacity_Index': M['edci'],
        'Equivalent_viscous_damping_ratio': M['zeta'],
    })
    print(f'\n------ DIRECTION {d} ------')
    S01.RANDOM_FOREST(df)
    S01.PLOT_HEATMAP(df)
    S01.MULTIPLE_REGRESSION(df)
