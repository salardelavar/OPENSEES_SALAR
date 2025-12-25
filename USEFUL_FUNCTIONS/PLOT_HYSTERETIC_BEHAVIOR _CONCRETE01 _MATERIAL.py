#**************************************************************
#                    >> IN THE NAME OF ALLAH << 
#  HYSTERETIC BEHAVIOR OF CONCRETE MATERIAL USING OPENSEES     
#  Cyclic Displacement-Controlled Analysis
#  THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
#**************************************************************

import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt

#%% ------------------------------------------------------------
# 1. INITIALIZE MODEL
ops.wipe()
ops.model('Basic', '-ndm', 1, '-ndf', 1)

#%% ------------------------------------------------------------
# 2. MATERIAL PROPERTIES (CONCRETE - CONCRETE01)
# Define parameters (units: mm, N)
# ------------------------------------------
# CONCRETE  tag   f'cec0   f'cuecu
# Parametric definitions for unconfined concrete

fcU = -35                             # [N/mm²] Unconfined concrete compressive strength
ec0U  = -0.002 * (abs(fcU)/30)**0.25  # [mm/mm] Initial strain at peak strength (semi-empirical)
fcUU  = 0.1 * fcU                     # [N/mm²] Ultimate stress (~10% of peak compressive strength)
ecuU  = ec0U * 3.5                    # [mm/mm] Ultimate strain (e.g., 3.5× ec0U)

# Parametric definitions for confined core concrete
Kc = 1.1
fcC   = fcU * Kc                      # [N/mm²] Confined strength (increased by confinement factor Kc)
ec0C  = ec0U * 1.8                    # [mm/mm] Peak strain increases with confinement
fcUC  = 0.95 * fcC                    # [N/mm²] Ultimate confined stress (~95% of peak)
ecuC  = 1.475 * ecuU * Kc             # [mm/mm] Ultimate confined strain (larger ductility due to confinement)

coreTag, coverTag = 100, 200 
ops.uniaxialMaterial('Concrete01', coreTag, fcC, ec0C, fcUC, ecuC)  # Core concrete (confined)
ops.uniaxialMaterial('Concrete01', coverTag, fcU, ec0U, fcUU, ecuU) # Cover concrete (unconfined)
# INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Concrete01_Material_--_Zero_Tensile_Strength
#%% ------------------------------------------------------------
# 3. MODEL GEOMETRY (ZERO-LENGTH ELEMENT)
ops.node(1, 0.0)
ops.node(2, 0.0)

ops.fix(1, 1)   # Fixed node
ops.fix(2, 0)   # Free node

ops.element('zeroLength', 1, 1, 2, '-mat', coreTag, '-dir', 1)

#%% ------------------------------------------------------------
# 4. CYCLIC STRAIN PROTOCOL

# Key strain points (same logic as your protocol)
key_disp = np.array([
     0.0,
     0.005*ecuC,   -0.005*ecuC,
     0.01*ecuC,   -0.01*ecuC,
     0.05*ecuC,   -0.05*ecuC,
     0.09*ecuC,   -0.09*ecuC,
     0.20*ecuC,   -0.20*ecuC,
     0.40*ecuC,   -0.40*ecuC,
     0.60*ecuC,   -0.60*ecuC,
     0.80*ecuC,   -0.80*ecuC,
     1.00*ecuC,   -1.00*ecuC,
     0.0
])

#%% ------------------------------------------------------------
# Generate 1000-point strain protocol
n_points = 1000
disp_protocol = np.interp(
    np.linspace(0, len(key_disp) - 1, n_points),
    np.arange(len(key_disp)),
    key_disp
)


#%% ------------------------------------------------------------
# 5. ANALYSIS SETUP
ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)
ops.load(2, 1.0)

ops.system('BandGeneral')
ops.numberer('Plain')
ops.constraints('Plain')
ops.algorithm('Newton')
ops.analysis('Static')

#%% ------------------------------------------------------------
# 6. RUN ANALYSIS & RECORD RESPONSE
disp = []
force = []

for target_disp in disp_protocol:
    current_disp = ops.nodeDisp(2, 1)
    dU = target_disp - current_disp

    ops.integrator('DisplacementControl', 2, 1, dU)
    ops.analyze(1)

    disp.append(ops.nodeDisp(2, 1))
    force.append(ops.eleForce(1)[0])

#%% ------------------------------------------------------------
# 7. PLOT HYSTERESIS LOOP
plt.figure(figsize=(7, 6))
plt.plot(disp, force, '-bo', linewidth=1.5)
plt.xlabel('Strain (mm/mm)')
plt.ylabel('Stress (N/mm^2)')
plt.title('Hysteresis Loop - Hysteretic Concrete Material')
plt.grid(True)
plt.tight_layout()
plt.show()
