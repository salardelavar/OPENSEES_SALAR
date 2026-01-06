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
# 2. MATERIAL PROPERTIES (STEEL - Dodd_Restrepo)
Es = 210e4        # [N/mm^2] Young's modulus
fy = 240          # [N/mm^2] Yield stress
ey = fy / Es      # [mm/mm] Yield strain
esh = ey * 5.2    # [mm/mm] Tensile strain at initiation of strain hardening
fu = 1.1818 * fy  # [N/mm^2] Ultimate stress
esu = 0.12        # [mm/mm] Ultimate strain
Esh = (fu-fy)/(esu-ey) # [N/mm^2] Tensile strain at initiation of strain hardening
eshi = (esu + 5*esh)/5   # [mm/mm] Tensile strain for a point on strain hardening curve, recommended range of values for ESHI: [ (ESU + 5*ESH)/6, (ESU + 3*ESH)/4]
fshi = fu * 0.5          # [N/mm^2] Tensile stress at point on strain hardening curve corresponding to ESHI
OmegaFac = 0.8     # hardening curve. Range: [0.75, 1.15]. Largest value tends to near a bilinear Bauschinger curve. Default = 1.0

matTag = 1

# Definition of Dodd_Restrepo material
# uniaxialMaterial Dodd_Restrepo $tag $Fy $Fsu $ESH $ESU $Youngs $ESHI $FSHI <$OmegaFac>
ops.uniaxialMaterial('Dodd_Restrepo', matTag, fy, fu, esh, esu, Es, eshi, -ey, fshi, OmegaFac)
# INFO LINK: https://opensees.berkeley.edu/wiki/index.php/DoddRestrepo
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
     #0.5*esu,  -0.5*esu,
     #0.9*esu,  -0.9*esu,
     #1.0*esu,  -1.0*esu,
     #1.1*esu,  -1.1*esu,
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
ops.test('NormDispIncr', 1e-8, 10000) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/test.html
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
plt.title('Hysteresis Loop - Dodd_Restrepo Steel Material')
plt.grid(True)
plt.tight_layout()
plt.show()
