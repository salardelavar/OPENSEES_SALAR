################################################################################################################
#                        >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                      #
#   DIRECT DISPLACEMENT-BASED DESIGN (DDBD) OF A SDOF SYSTEM WITH NONLINEAR DYNAMIC TIME-HISTORY IN OPENSSEES  #
#--------------------------------------------------------------------------------------------------------------#
#                        THIS PYTHON SCRIPT IS WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                   # 
#                                          E-mail:salar.d.ghashghaei@gmail.com                                 #
################################################################################################################
"""
This Python script implements Direct Displacement-Based Design (DDBD) for a single-degree-of-freedom (SDOF) system and
 verifies it with nonlinear dynamic time-history analysis in OpenSeesPy. It begins by defining a target displacement,
 estimating ductility and equivalent viscous damping, then constructs a multilinear force-displacement backbone that is
 scaled so the yield strength matches the base shear required by DDBD; an effective period is found by intersecting the
 target displacement with a damped design displacement spectrum. Parallel zero-length hysteretic springs representing the
 scaled backbone, together with Rayleigh damping and a scaled ground-motion record, are assembled into an OpenSees model
 that is subjected to transient analysis; the resulting displacement history, hysteresis loops, and design spectrum are
 finally plotted to compare the achieved peak displacement against the DDBD target.
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
import openseespy.opensees as ops
import PLOT_FORCE_DISP_FUN as S01
import ANALYSIS_FUNCTION as S02
#%% ------------------------------------
#  1. DDBD PARAMETERS
# Target performance
DELTA_D = 0.015          # [m] Target displacement at effective heigh --> 2.5% drift for H=1m SDOF
DELTA_Y_EST = 0.008      # [m] Estimated yield displacement
MU = DELTA_D / DELTA_Y_EST

# Equivalent viscous damping (example for RC frame - Takeda model)
XI_EQ = 0.05 + 0.565 * (MU - 1) / MU
print(f"Ductility μ = {MU:.2f}")
print(f"Equivalent damping ξeq = {XI_EQ*100:.1f}%")

# ----- Backbone (units: N, m, kg, s) -----
FY, DY = 2000.0, 0.01                 # [N], [m] Yield force & displacement
FU, DU = 1.2 * FY, DY * 10            # [N], [m] Ultimate (peak) force & displacement
 
FR1, DR1 = FY / 4,  DY * 12           # [N], [m] 1st residual / post-peak point
FR2, DR2 = FY / 10, DY * 15           # [N], [m] 2nd residual / post-peak point
 
Force = np.array([FY, FU, FR1, FR2])
Displacement = np.array([DY, DU, DR1, DR2])

"""
Coff is an array of stiffness scaling coefficients (one per parallel spring).
 Each value multiplies the shared multilinear force–displacement backbone,
 so the springs can have different initial stiffnesses and force capacities
 while keeping the same shape of the nonlinear curve. In the analytical solution
 the coefficients scale the segment stiffnesses; in OpenSeesPy they scale the
 material points of each zeroLength spring. When all entries equal 1.0 the springs
 are identical; different values produce a set of parallel springs with proportional
 strength and stiffness.
"""
Coff = np.array([1.0, 1.0, 1.0, 1.0, 1.0])   # STIFFNESS COEFFICIENTS

S01.PLOT_FORCE_DISP_FUN(Displacement, Force)

K0_initial = np.sum(Coff) * Force[0] / Displacement[0]
T0_target = 0.5
MASS = K0_initial * (T0_target / (2 * np.pi))**2
#%% ------------------------------------
# 2. Design Displacement Spectrum
def design_displacement_spectrum(T, xi, Sd_5pct_func):
    """
    Simple damping modification factor (EC8 / NZ style)
    """
    if xi <= 0.05:
        eta = 1.0
    else:
        eta = np.sqrt(10 / (5 + 100*xi))   # EC8 damping reduction
        eta = max(eta, 0.55)
    return Sd_5pct_func(T) * eta

# Simple design spectrum (you can replace with NZS 1170.5 or any code spectrum)
def Sd_5pct(T):
    # Example bilinear displacement spectrum (plateau after 2.5 s)
    # Units: meteres
    Sd_max = 0.12          # peak spectral displacement
    Tc = 2.5               # corner period
    if T <= Tc:
        return (Sd_max / Tc) * T
    else:
        return Sd_max

# Find Teff from target displacement
from scipy.optimize import fsolve
def find_Teff(delta_d, xi):
    def func(T):
        return design_displacement_spectrum(T[0], xi, Sd_5pct) - delta_d
    Teff = fsolve(func, 1.0)[0]
    return Teff

Teff = find_Teff(DELTA_D, XI_EQ)
print(f"Effective period Teff = {Teff:.3f} s")

#%% ------------------------------------
# 3. Required Strength from DDBD
Ke = 4 * np.pi**2 * MASS / Teff**2
Vb = Ke * DELTA_D          # Design base shear

print(f"Effective stiffness Ke = {Ke:.1f} N/m")
print(f"Design base shear Vb  = {Vb:.1f} N")

# Scale the backbone so that the yield force matches the required strength
# (simple approach: keep shape, scale force capacity)
scale_force = Vb / FY
Force_DDBD = Force * scale_force
print(f"Scaled yield force = {Force_DDBD[0]:.1f} N")

#%% ------------------------------------
# 4. OpenSees Model
def OPENSEES_DDBD(Force, Displacement, Coff, ag_raw, dt, mass, xi, GM_FACTOR):
    
    MAX_ITERATIONS = 500     # Convergence iteration for test
    MAX_TOLERANCE = 1.0e-6   # Convergence tolerance for test
    
    # Reset model
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)
    
    # Nodes
    ops.node(1, 0.0)
    ops.node(2, 0.0)
    ops.fix(1, 1)
    
    # Mass
    ops.mass(2, mass)

    for j in range(len(Coff)):
        tag = 1 + j
        pos, neg = [], []
        for f, d in zip(Force, Displacement):
            pos.extend([Coff[j] * f, d])
            neg.extend([-Coff[j] * f, -d])
        ops.uniaxialMaterial('HystereticSM', tag, *pos, *neg, 1.0, 1.0, 0.0, 0.0, 0.0)
        # INFO LINK: https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/HystereticSM.html
        ops.element('zeroLength', tag, 1, 2, '-mat', tag, '-dir', 1)
        # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/ZeroLength.html

    # mass-proportional damping: c = alphaM * m = 2*xi*omega0*m
    K0 = np.sum(Coff) * Force[0] / Displacement[0]
    omega0 = np.sqrt(K0 / mass)
    ops.rayleigh(2.0 * xi * omega0, 0.0, 0.0, 0.0)
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/reyleigh.html

    # ground motion
    ops.timeSeries('Path', 1, '-dt', dt, '-values', *ag_raw, '-factor', GM_FACTOR)
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/pathTs.html
    ops.pattern('UniformExcitation', 1, 1, '-accel', 1)
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/uniformExcitation.html

    ops.system('BandGeneral')
    ops.numberer('Plain')
    ops.constraints('Plain')
    ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS)
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/normDispIncr.html
    ops.algorithm('Newton')              # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/algorithm.html
    ops.integrator('Newmark', 0.5, 0.25) # INFO LINK: https://opensees.github.io/OpenSeesDocumentation/user/manual/analysis/integrator/Newmark.html
    ops.analysis('Transient')            # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/analysis.html

    t, u, f = [0.0], [0.0], [0.0]
    for i in range(1, len(ag_raw)):
        OK = ops.analyze(1, dt)
        if OK != 0:
            S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS)
        ops.reactions()
        t.append(ops.getTime())
        u.append(ops.nodeDisp(2, 1))
        f.append(ops.nodeReaction(1, 1))
    return np.array(t), np.array(u), np.array(f)

#%% ------------------------------------
# 5. Run OpenSees
GM_FILE = 'Ground_Acceleration_X.txt'
DT = 0.01
GM_UNIT_G = 9.81
GM_SCALE = 50.0                # adjust so that demand is close to target
ag_raw = np.loadtxt(GM_FILE)
GM_FACTOR = GM_UNIT_G * GM_SCALE

t_os, U_os, F_os = OPENSEES_DDBD(Force_DDBD, Displacement, Coff,
                                      ag_raw, DT, MASS, XI_EQ, GM_FACTOR)

print(f"\nMax displacement from THA = {np.max(np.abs(U_os)):.4f} m")
print(f"Target displacement Δd   = {DELTA_D:.4f} m")
print(f"Ratio (achieved / target) = {np.max(np.abs(U_os))/DELTA_D:.2f}")

#%% ------------------------------------
# 6. Plots
fig, ax = plt.subplots(2, 2, figsize=(14, 9))

ax[0,0].plot(np.arange(len(ag_raw))*DT, ag_raw*GM_FACTOR, 'g')
ax[0,0].set_title('Ground Acceleration')
ax[0,0].set_xlabel('Time [s]')
ax[0,0].grid(True)

ax[0,1].plot(t_os, U_os, 'b')
ax[0,1].axhline(DELTA_D, color='r', ls='--', label=f'Target Δd = {DELTA_D} m')
ax[0,1].axhline(-DELTA_D, color='r', ls='--')
ax[0,1].set_title('Displacement Time History')
ax[0,1].set_xlabel('Time [s]')
ax[0,1].legend()
ax[0,1].grid(True)

ax[1,0].plot(U_os, F_os, 'b')
ax[1,0].set_title('Hysteresis (DDBD designed strength)')
ax[1,0].set_xlabel('Displacement [m]')
ax[1,0].set_ylabel('Force [N]')
ax[1,0].grid(True)

# Displacement spectrum used
T_plot = np.linspace(0.1, 4, 100)
Sd_plot = [design_displacement_spectrum(T, XI_EQ, Sd_5pct) for T in T_plot]
ax[1,1].plot(T_plot, Sd_plot, 'k', lw=2)
ax[1,1].plot(Teff, DELTA_D, 'ro', ms=10, label=f'Teff = {Teff:.2f}s')
ax[1,1].set_title('Design Displacement Spectrum (damped)')
ax[1,1].set_xlabel('Period [s]')
ax[1,1].set_ylabel('Sd [m]')
ax[1,1].legend()
ax[1,1].grid(True)
plt.savefig('SDOF_DIRECT_DISPLACEMENT-BASED_DESIGN.png', dpi=150)
plt.tight_layout()
plt.show()