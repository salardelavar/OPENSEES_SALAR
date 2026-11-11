################################################################################################################
#                        >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                      #
#   DIRECT DISPLACEMENT-BASED DESIGN (DDBD) OF A MDOF SHEAR BUILDING WITH NONLINEAR TIME-HISTORY IN OPENSEES   #
#--------------------------------------------------------------------------------------------------------------#
#                        THIS PYTHON SCRIPT IS WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                   # 
#                                          E-mail:salar.d.ghashghaei@gmail.com                                 #
################################################################################################################
"""
This Python script implements Direct Displacement-Based Design (DDBD) for a multi-degee-of-freedom (MDOF) system and
 verifies it with nonlinear dynamic time-history analysis in OpenSeesPy. It begins by defining a target displacement,
 estimating ductility and equivalent viscous damping, then constructs a multilinear force-displacement backbone that is
 scaled so the yield strength matches the base shear required by DDBD; an effective period is found by intersecting the
 target displacement with a damped design displacement spectrum. Parallel zero-length hysteretic springs representing the
 scaled backbone, together with Rayleigh damping and a scaled ground-motion record, are assembled into an OpenSees model
 that is subjected to transient analysis; the resulting displacement history, hysteresis loops, and design spectrum are
 finally plotted to compare the achieved peak displacement against the DDBD target.
 
DDBD of an N-storey shear building (Priestley, Calvi & Kowalsky 2007) verified by nonlinear time-history.

PROCEDURE
  1. Design displacement profile  Delta_i  (inelastic shape + higher-mode drift reduction omega_theta)
  2. Substitute (equivalent SDOF) structure:   Delta_d, m_e, H_e
  3. Yield displacement  Delta_y = theta_y * H_e   ->   ductility mu   ->   equivalent damping xi_eq
  4. Effective period T_eff from the xi_eq-damped design displacement spectrum at Delta_d
  5. Effective stiffness K_e = 4 pi^2 m_e / T_eff^2   ->   design base shear V_b = K_e * Delta_d
  6. Distribute V_b to the floors:  F_i = F_t + (V_b - F_t) * m_i*Delta_i / sum(m_j*Delta_j)
     -> story shears V_i -> story strengths (STORY_COFF scale factors of the shared backbone)
  7. Verify with nonlinear time-history (OpenSeesPy, reuses MDOF_DYNAMIC_TIME_HISTORY.py), the record
     being scaled so that its xi_eq-damped linear Sd(T_eff) equals the design spectrum at T_eff.

SIMPLIFICATIONS (flagged, change them if your code requires it)
  - P-Delta is ignored (Priestley: V_b -= C*P*Delta_d/H_e)
  - Constant story height and the same yield drift theta_y in every story
  - Rayleigh damping in the THA is the ELASTIC 5 % only; the hysteretic energy is already in the spring
    (using xi_eq in the THA would count the damping twice)
"""

# PAPER: Fundamentals of Direct Displacement-Based Seismic Design and Assessment
'https://link.springer.com/chapter/10.1007/978-3-211-74214-3_8'
# PAPER: Direct Displacement-Based Seismic Design of Structures
'http://db.nzsee.org.nz/2007/Paper18.pdf' 
# PAPER: Highlighting Differences between Force-Based and Displacement-Based Design Solutions for Reinforced Concrete Frame Structures - Timothy John Sullivan
'https://www.tandfonline.com/doi/abs/10.2749/101686613X13439149156958'
# PAPER: 30 YEARS OF DIRECT DISPLACEMENT-BASED DESIGN: WHERE WE STARTED, WHAT'S BEEN DONE, WHERE WE'RE GOING
'https://ecopress.gr/wp-content/uploads/WCEE-2024-.pdf'
#YOUTUBE: Displacement-based seismic design of structures
'https://www.youtube.com/watch?v=MZUhSHmIUdI'

import numpy as np
import matplotlib.pyplot as plt
import PLOT_FORCE_DISP_FUN as S01
import MDOF_DYNAMIC_TIME_HISTORY as S03

#%% ------------------------------------------------------------------
# 1. BUILDING AND DESIGN PARAMETERS   (units: N, m, kg, s)
N_STORY = 10
H_STORY = 3.0                                  # [m] constant story height
MASS = np.full(N_STORY, 200.0e3)               # [kg] floor masses (edit for non-uniform buildings)
H = H_STORY * np.arange(1, N_STORY + 1)        # [m] floor heights
HN = H[-1]

THETA_C = 0.025                                # code/performance drift limit at the critical (first) story
THETA_Y = 0.010                                # story yield drift (RC frame: 0.5*eps_y*Lb/hb)
FT_RATIO = 0.10                                # fraction of V_b applied as extra force at the roof
XI_ELASTIC = 0.05                              # elastic viscous damping (THA Rayleigh and spectrum base)

# Design displacement spectrum at 5 % damping: linear up to the corner period TD, then constant
SD_MAX = 0.90                                  # [m] plateau displacement
TD = 5.0                                       # [s] corner period

#%% ------------------------------------------------------------------
# 2. DESIGN DISPLACEMENT PROFILE
if N_STORY <= 4:
    delta = H / HN
else:
    delta = (4.0 / 3.0) * (H / HN) * (1.0 - H / (4.0 * HN))      # inelastic shape (n > 4)

omega_theta = min(1.0, 1.15 - 0.0034 * HN)                       # higher-mode drift reduction (HN in m)
theta_d = omega_theta * THETA_C
Delta = delta * (theta_d * H[0]) / delta[0]                      # critical story = story 1 -> drift = theta_d
drift_design = np.diff(np.concatenate(([0.0], Delta)))

#%% ------------------------------------------------------------------
# 3. SUBSTITUTE (EQUIVALENT SDOF) STRUCTURE
DELTA_D = np.sum(MASS * Delta**2) / np.sum(MASS * Delta)         # design displacement
M_E = np.sum(MASS * Delta) / DELTA_D                             # effective mass
H_E = np.sum(MASS * Delta * H) / np.sum(MASS * Delta)            # effective height

DELTA_Y = THETA_Y * H_E                                          # yield displacement
MU = DELTA_D / DELTA_Y
XI_EQ = XI_ELASTIC + 0.565 * (MU - 1.0) / (MU * np.pi) if MU > 1.0 else XI_ELASTIC   # Takeda thin, RC frame

print(f"omega_theta = {omega_theta:.3f},  design drift theta_d = {theta_d*100:.2f} %")
print(f"Delta_d = {DELTA_D:.4f} m,  m_e = {M_E:.1f} kg ({M_E/np.sum(MASS)*100:.0f} % of total),  H_e = {H_E:.2f} m")
print(f"Delta_y = {DELTA_Y:.4f} m,  ductility mu = {MU:.2f},  xi_eq = {XI_EQ*100:.1f} %")

#%% ------------------------------------------------------------------
# 4. DESIGN DISPLACEMENT SPECTRUM AND EFFECTIVE PERIOD
def SD_5PCT(T):
    return SD_MAX * np.minimum(np.asarray(T, float), TD) / TD

def DAMPING_FACTOR(xi):
    # EC8-style reduction (same form as the SDOF script), bounded below by 0.55
    return 1.0 if xi <= 0.05 else max(np.sqrt(10.0 / (5.0 + 100.0 * xi)), 0.55)

def DESIGN_SD(T, xi):
    return SD_5PCT(T) * DAMPING_FACTOR(xi)

ETA = DAMPING_FACTOR(XI_EQ)
if DELTA_D > ETA * SD_MAX:
    raise ValueError(f"Delta_d = {DELTA_D:.3f} m exceeds the spectrum plateau {ETA*SD_MAX:.3f} m: "
                     "no T_eff exists. Reduce the drift limit or increase SD_MAX.")
T_EFF = DELTA_D * TD / (ETA * SD_MAX)                            # closed form on the linear branch
print(f"Effective period T_eff = {T_EFF:.3f} s")

#%% ------------------------------------------------------------------
# 5. DESIGN BASE SHEAR AND FORCE DISTRIBUTION
K_E = 4.0 * np.pi**2 * M_E / T_EFF**2
VB = K_E * DELTA_D
print(f"Effective stiffness K_e = {K_E:.4e} N/m,  design base shear V_b = {VB:.4e} N "
      f"({VB/(np.sum(MASS)*9.81)*100:.1f} % of weight)")

F_FLOOR = (1.0 - FT_RATIO) * VB * MASS * Delta / np.sum(MASS * Delta)
F_FLOOR[-1] += FT_RATIO * VB
V_DESIGN = np.cumsum(F_FLOOR[::-1])[::-1]                        # story shears, story 1 = base

#%% ------------------------------------------------------------------
# 6. STORY BACKBONES: shared shape, strength scaled by STORY_COFF
FY_REF, DY = 1000.0, THETA_Y * H_STORY                           # reference yield force and story yield drift
Force = FY_REF * np.array([1.0, 1.2, 0.25, 0.10])
Displacement = DY * np.array([1.0, 10.0, 12.0, 15.0])

DESIGN_ON_BACKBONE = True    # True : V_i is reached at the design drift of story i (secant/DDBD consistent)
                             # False: V_i is used as yield strength (same as the SDOF script)
if DESIGN_ON_BACKBONE:
    F_at_design = np.interp(drift_design, np.concatenate(([0.0], Displacement)), np.concatenate(([0.0], Force)))
    STORY_COFF = V_DESIGN / F_at_design
else:
    STORY_COFF = V_DESIGN / FY_REF

S01.PLOT_FORCE_DISP_FUN(Displacement, Force)

omega, _, _, _ = S03.MODAL_PROPERTIES(Displacement, Force, STORY_COFF, MASS)
print("Elastic periods of the designed building (first 3) [s]:", np.round(2 * np.pi / omega[:3], 3))
print("Story shears V_i [N]:", np.round(V_DESIGN, 0))

#%% ------------------------------------------------------------------
# 7. GROUND MOTION SCALED TO THE DESIGN SPECTRUM, THEN OPENSEES THA
def SD_LINEAR(ag, dt, T, xi):
    """Peak displacement of a linear SDOF (Newmark average acceleration), ag in m/s^2."""
    beta, gamma = 0.25, 0.5
    w = 2 * np.pi / T
    c, k = 2 * xi * w, w**2
    keff = k + gamma / (beta * dt) * c + 1.0 / (beta * dt**2)
    u = v = a = 0.0
    umax = 0.0
    for i in range(1, len(ag)):
        rhs = (-ag[i] + (u / (beta * dt**2) + v / (beta * dt) + (1 / (2 * beta) - 1) * a)
               + c * (gamma / (beta * dt) * u + (gamma / beta - 1) * v + dt * (gamma / (2 * beta) - 1) * a))
        un = rhs / keff
        an = (un - u) / (beta * dt**2) - v / (beta * dt) - (1 / (2 * beta) - 1) * a
        v = v + dt * ((1 - gamma) * a + gamma * an)
        u, a = un, an
        umax = max(umax, abs(u))
    return umax

GM_FILE = 'Ground_Acceleration_X.txt'
DT = 0.01
GM_UNIT_G = 9.81
AUTO_SCALE = True              # scale the record so Sd_record(T_eff, xi_eq) = design spectrum at T_eff = Delta_d
GM_SCALE = 1.0                 # used only when AUTO_SCALE = False

ag_raw = np.loadtxt(GM_FILE)
if AUTO_SCALE:
    GM_SCALE = DELTA_D / SD_LINEAR(ag_raw * GM_UNIT_G, DT, T_EFF, XI_EQ)
    print(f"Auto scale factor = {GM_SCALE:.3f}")
GM_FACTOR = GM_UNIT_G * GM_SCALE

MODES = (1, 2)
MAX_ITERATIONS, MAX_TOLERANCE = 500, 1.0e-6
OS = S03.MDOF_OPENSEES_DYNAMIC('INELASTIC', Displacement, Force, STORY_COFF, MASS, ag_raw, DT,
                               XI_ELASTIC, GM_FACTOR, MAX_ITERATIONS, MAX_TOLERANCE, MODES=MODES)

#%% ------------------------------------------------------------------
# 8. VERIFICATION
U_env = np.abs(OS['U_floor']).max(axis=0)
D_env = np.abs(OS['drift']).max(axis=0)
V_env = np.abs(OS['V_story']).max(axis=0)
DELTA_D_THA = np.sum(MASS * U_env**2) / np.sum(MASS * U_env)

print(f"\nRoof displacement  : THA = {U_env[-1]:.4f} m   design = {Delta[-1]:.4f} m   ratio = {U_env[-1]/Delta[-1]:.2f}")
print(f"Equivalent Delta_d : THA = {DELTA_D_THA:.4f} m   design = {DELTA_D:.4f} m   ratio = {DELTA_D_THA/DELTA_D:.2f}")
print(f"Max story drift    : THA = {D_env.max()/H_STORY*100:.2f} %   limit  = {THETA_C*100:.2f} %  (story {D_env.argmax()+1})")
print(f"Base shear         : THA = {V_env[0]:.4e} N   design V_b = {VB:.4e} N   ratio = {V_env[0]/VB:.2f}")

#%% ------------------------------------------------------------------
# 9. PLOTS
fig, ax = plt.subplots(2, 4, figsize=(22, 10))
story = np.arange(1, N_STORY + 1)
t = OS['t']
roof = OS['U_floor'][:, -1]

ax[0, 0].plot(np.arange(len(ag_raw)) * DT, ag_raw * GM_FACTOR, 'g')
ax[0, 0].set_title(f'Scaled ground acceleration (x{GM_SCALE:.2f})'); ax[0, 0].set_xlabel('Time [s]')
ax[0, 0].set_ylabel('ag [m/s²]')

ax[0, 1].plot(t, roof, 'b')
ax[0, 1].axhline(Delta[-1], color='r', ls='--', label=f'Design roof Δn = {Delta[-1]:.3f} m')
ax[0, 1].axhline(-Delta[-1], color='r', ls='--')
ax[0, 1].set_title('Roof displacement'); ax[0, 1].set_xlabel('Time [s]'); ax[0, 1].legend()

ax[0, 2].plot(roof, OS['base_reaction'], 'b')
ax[0, 2].set_title('Global hysteresis (DDBD designed strength)')
ax[0, 2].set_xlabel('Roof displacement [m]'); ax[0, 2].set_ylabel('Base reaction [N]')

T_plot = np.linspace(0.05, 6, 200)
ax[0, 3].plot(T_plot, SD_5PCT(T_plot), 'k--', label='5 % damping')
ax[0, 3].plot(T_plot, [DESIGN_SD(T, XI_EQ) for T in T_plot], 'k', lw=2, label=f'ξeq = {XI_EQ*100:.1f} %')
ax[0, 3].plot(T_EFF, DELTA_D, 'ro', ms=10, label=f'Teff = {T_EFF:.2f} s, Δd = {DELTA_D:.3f} m')
ax[0, 3].set_title('Design displacement spectrum'); ax[0, 3].set_xlabel('Period [s]')
ax[0, 3].set_ylabel('Sd [m]'); ax[0, 3].legend()

ax[1, 0].plot(Delta, story, 'r--o', lw=2, label='Design')
ax[1, 0].plot(U_env, story, 'k-s', lw=2, label='THA envelope')
ax[1, 0].set_title('Displacement profile'); ax[1, 0].set_xlabel('Displacement [m]'); ax[1, 0].set_ylabel('Floor')
ax[1, 0].legend()

ax[1, 1].plot(drift_design / H_STORY * 100, story, 'r--o', lw=2, label='Design')
ax[1, 1].plot(D_env / H_STORY * 100, story, 'k-s', lw=2, label='THA envelope')
ax[1, 1].axvline(THETA_C * 100, color='gray', ls=':', label='Drift limit')
ax[1, 1].set_title('Story drift profile'); ax[1, 1].set_xlabel('Drift [%]'); ax[1, 1].set_ylabel('Story')
ax[1, 1].legend()

ax[1, 2].plot(V_DESIGN, story, 'r--o', lw=2, label='Design Vi')
ax[1, 2].plot(V_env, story, 'k-s', lw=2, label='THA envelope')
ax[1, 2].set_title('Story shear profile'); ax[1, 2].set_xlabel('Story shear [N]'); ax[1, 2].set_ylabel('Story')
ax[1, 2].legend()

ax[1, 3].plot(t, OS['base_reaction'], 'b')
ax[1, 3].axhline(VB, color='r', ls='--', label=f'Design Vb = {VB:.3e} N')
ax[1, 3].axhline(-VB, color='r', ls='--')
ax[1, 3].set_title('Base shear time history'); ax[1, 3].set_xlabel('Time [s]'); ax[1, 3].legend()

for a_ in ax.ravel():
    a_.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('MDOF_DIRECT_DISPLACEMENT-BASED_DESIGN.png', dpi=150)
plt.show()
