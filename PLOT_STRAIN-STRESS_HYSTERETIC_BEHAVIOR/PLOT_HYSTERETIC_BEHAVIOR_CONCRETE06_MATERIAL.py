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
# 2. MATERIAL PROPERTIES (CONCRETE - CONCRETE06)
# Define parameters (units: mm, N)
# ------------------------------------------
# CONCRETE  tag   f'cec0   f'cuecu
# Parametric definitions for unconfined concrete

fc =  -35	               # [N/mm^2] concrete compressive strength (compression is negative)*
Ec = 4700 * np.sqrt(-fc)   # [N/mm^2] Unconfined Concrete Elastic Modulus
e0 = 2*fc/Ec               # [mm/mm] strain at compressive strength*
n = 2	                   # compressive shape factor
k = 1	                   # post-peak compressive shape factor
alpha1 = 0.32              # parameter for compressive plastic strain definition
fcr = 0.7 * np.sqrt(-fc)   # [N/mm^2] tensile strength
ecr	= -e0                  # [mm/mm] tensile strain at peak stress ($fcr)
b = 4                      # exponent of the tension stiffening curve
alpha2 = 0.8               # parameter for tensile plastic strain definition

coreTag = 1       # Confined Concrete Section Tag
# uniaxialMaterial Concrete06 $matTag $fc $e0 $n $k $alpha1 $fcr $ecr $b $alpha2
ops.uniaxialMaterial('Concrete06', coreTag, fc, e0, n, k, alpha1, fcr, ecr, b, alpha2) # build core concrete (confined)
# INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Concrete06_Material
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
     0.005*e0,   -0.005*e0,
     0.01*e0,   -0.01*e0,
     0.05*e0,   -0.05*e0,
     0.90*e0,   -0.09*e0,
     1.20*e0,   -0.20*e0,
     1.40*e0,   -0.40*e0,
     1.60*e0,   -0.60*e0,
     1.80*e0,   -0.80*e0,
     2.00*e0,   -1.00*e0,
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
plt.title('Hysteresis Loop - Hysteretic Concrete06 Material')
plt.grid(True)
plt.tight_layout()
plt.show()
