#**************************************************************
#                    >> IN THE NAME OF ALLAH << 
#  HYSTERETIC BEHAVIOR OF STEEL MATERIAL USING OPENSEES     
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
# 2. MATERIAL PROPERTIES (STEEL - ElasticPP)
Es = 210e4        # [N/mm^2] Young's modulus
fy = 240          # [N/mm^2] Yield stress
ey = fy / Es      # [mm/mm] Yield strain
fu = 1.1818 * fy  # [N/mm^2] Ultimate stress
esu = ey * 75.2   # [mm/mm] Ultimate strain
Esh = (fu-fy)/(esu-ey)
b = Esh/Es      # strain hardening ratio


matTag = 1

# Definition of Steel01 material
ops.uniaxialMaterial('ElasticPP', matTag, Es, ey, -ey, 0.01*ey)
# INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Elastic-Perfectly_Plastic_Material
#%% ------------------------------------------------------------
# 3. MODEL GEOMETRY (ZERO-LENGTH ELEMENT)
ops.node(1, 0.0)
ops.node(2, 0.0)

ops.fix(1, 1)   # Fixed node
ops.fix(2, 0)   # Free node

ops.element('zeroLength', 1, 1, 2, '-mat', matTag, '-dir', 1)

#%% ------------------------------------------------------------
# 4. CYCLIC STRAIN PROTOCOL

# Key strain points (same logic as your protocol)
key_disp = np.array([
     0.0,
     0.1*ey,   -0.1*ey,
     0.5*ey,   -0.5*ey,
     0.8*ey,   -0.8*ey,
     ey,       -ey,
     1.5*ey,   -1.5*ey,
     2*ey,     -2*ey,
     4*ey,     -4*ey,
     6*ey,     -6*ey,
     0.5*esu,  -0.5*esu,
     0.9*esu,  -0.9*esu,
     1.0*esu,  -1.0*esu,
     1.1*esu,  -1.1*esu,
     0.0
])

#%% ------------------------------------------------------------
# Generate 1000-point strain protocol
n_points = 10000
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
plt.title('Hysteresis Loop - ElasticPP Steel Material')
plt.grid(True)
plt.tight_layout()
plt.show()
