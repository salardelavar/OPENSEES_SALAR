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
# 2. MATERIAL PROPERTIES (CONCRETE - CONCRETE04)
# Define parameters (units: mm, N)
# ------------------------------------------
# CONCRETE  tag   f'cec0   f'cuecu
# Parametric definitions for unconfined concrete

fc = -35 # [N/mm^2] Nominal concrete compressive strength

# confined concrete
Kfc = 1.3;			    # ratio of confined to unconfined concrete strength
fc1C = Kfc*fc;		    # CONFINED concrete (mander model), maximum stress
EcC = 4700 * np.sqrt(-fc1C) # [N/mm^2] Confined Concrete Elastic Modulus
eps1C = 2*fc1C/EcC;	    # strain at maximum stress 
fc2C = 0.2*fc1C;		# ultimate stress
eps2C = 5*eps1C;		# strain at ultimate stress 
# unconfined concrete
fc1U = fc;			    # UNCONFINED concrete (todeschini parabolic model), maximum stress
EcU = 4700 * np.sqrt(-fc1U) # [N/mm^2] Unconfined Concrete Elastic Modulus
eps1U = 2*fc1U/EcU;		# strain at maximum strength of unconfined concrete
fc2U = 0.2*fc1U;		# ultimate stress
eps2U = -0.012;			# strain at ultimate stress
beta = 1.01;		    # floating point value defining the exponential curve parameter to define the residual stress
# tensile-strength properties
ftC = -0.7 * np.sqrt(fc1C);		# Confined tensile strength +tension
ftU = -0.7 * np.sqrt(fc1U);		# Unconfined tensile strength +tension

EtsC = ftC/np.abs(eps1C);		# Confined softening stiffness
EtsU = ftU/np.abs(eps1U);		# Unconfined softening stiffness



coreTag = 1       # Confined Concrete Section Tag
coverTag = 2      # Unconfined Concrete Section Tag
# uniaxialMaterial Concrete04 $matTag $fc $ec $ecu $Ec <$fct $et> <$beta>
ops.uniaxialMaterial('Concrete04', coreTag, fc1C, eps1C, eps2C, EcC,ftC, EtsC, beta) # build core concrete (confined)
ops.uniaxialMaterial('Concrete04', coverTag, fc1U, eps1U, eps2U, EcU,ftU, EtsU, beta) # build cover concrete (unconfined)
# INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Concrete04_Material_--_Popovics_Concrete_Material
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
plt.title('Hysteresis Loop - Hysteretic Concrete04 Material')
plt.grid(True)
plt.tight_layout()
plt.show()
