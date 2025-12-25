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
# 2. MATERIAL PROPERTIES (CONCRETE - CONCRETE02)
# Define parameters (units: mm, N)
# ------------------------------------------
# CONCRETE  tag   f'cec0   f'cuecu
# Parametric definitions for unconfined concrete

fc = -35 # [N/mm^2] Nominal concrete compressive strength
Ec = 4700 * np.sqrt(-fc) # [N/mm^2] Concrete Elastic Modulus

# confined concrete
Kfc = 1.3;			    # ratio of confined to unconfined concrete strength
fc1C = Kfc*fc;		    # CONFINED concrete (mander model), maximum stress
eps1C = 2*fc1C/Ec;	    # strain at maximum stress 
fc2C = 0.2*fc1C;		# ultimate stress
eps2C = 5*eps1C;		# strain at ultimate stress 
# unconfined concrete
fc1U = fc;			    # UNCONFINED concrete (todeschini parabolic model), maximum stress
eps1U = -0.0025;		# strain at maximum strength of unconfined concrete
fc2U = 0.2*fc1U;		# ultimate stress
eps2U = -0.012;			# strain at ultimate stress
Lambda = 0.1;			# ratio between unloading slope at $eps2 and initial slope $Ec
# tensile-strength properties
ftC = -0.55*fc1C;		# tensile strength +tension
ftU = -0.55*fc1U;		# tensile strength +tension
Ets = ftU/0.002;		# tension softening stiffness

coreTag = 1       # Confined Concrete Section Tag
coverTag = 2      # Unconfined Concrete Section Tag
ops.uniaxialMaterial('Concrete02', coreTag, fc1C, eps1C, fc2C, eps2C, Lambda, ftC, Ets) # build core concrete (confined)
ops.uniaxialMaterial('Concrete02', coverTag, fc1U, eps1U, fc2U, eps2U, Lambda, ftU, Ets) # build cover concrete (unconfined)
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
     0.005*eps2U,   -0.005*eps2U,
     0.01*eps2U,   -0.01*eps2U,
     0.05*eps2U,   -0.05*eps2U,
     0.09*eps2U,   -0.09*eps2U,
     0.20*eps2U,   -0.20*eps2U,
     0.40*eps2U,   -0.40*eps2U,
     0.60*eps2U,   -0.60*eps2U,
     0.80*eps2U,   -0.80*eps2U,
     1.00*eps2U,   -1.00*eps2U,
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
