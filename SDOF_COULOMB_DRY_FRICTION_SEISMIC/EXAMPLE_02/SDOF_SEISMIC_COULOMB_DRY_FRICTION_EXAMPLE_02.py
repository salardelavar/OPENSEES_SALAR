######################################################################################################################
#                          >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                          #
#            SIMULATE THE SEISMIC RESPONSE OF A SINGLE-DEGREE-OF-FREEDOM (SDOF) STRUCTURE WITH COULOMB               #
#                                                DRY FRICTION USING OPENSEES                                         #
#--------------------------------------------------------------------------------------------------------------------#
#                                                         EXAMPLE 02                                                 #
#--------------------------------------------------------------------------------------------------------------------#
#                              THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                            #
#                                       EMAIL: salar.d.ghashghaei@gmail.com                                          #
######################################################################################################################
"""
This Python code simulates a single-degree-of-freedom (SDOF) system with both a linear spring and a Coulomb friction damper.
The system's mass is 4000 tons, and the structural spring stiffness is 10000 kN/m. The Coulomb damper is modeled using two different approaches: a `Steel01` material for an elasto-plastic representation and a dedicated `CoulombDamper` material.
A transient analysis is performed, subjecting the system to a sinusoidal acceleration input of 0.4g at 1.5 Hz for 2000 steps with a time step of 0.01 seconds.
The code records and plots the displacement, damper forces from both models, and hysteresis loops against displacement for comparison. This allows for analysis of how each damping model behaves under seismic excitation.
"""
# WIKIPEDIA
'https://en.wikipedia.org/wiki/Friction#:~:text=Coulomb%20friction%2C%20named%20after%20Charles,to%20the%20net%20applied%20force.'
# BOOK: Differential Equations for Engineers - Wei-Chau Xie - CAMBRIDGE
'https://www.cambridge.org/core/books/differential-equations-for-engineers/1B8F1A62BF6F98EB972CFE9114FA8B84'
# portwooddigital Website by Michael H. Scott
'https://portwooddigital.com/2024/02/10/two-sprung-masses-and-some-friction-force/'
#%%------------------------------------------------------------------------------------------------
import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
#%%------------------------------------------------------------------------------------------------
# Parameters
mass = 4000.0               # [Tons]
k_structure = 10000.0       # [kN/m]
mu = 0.25                   # Coefficient of kinetic friction
g = 9.81
F_friction = mu * mass * g  # [kN] Sliding force 
weight = mass * g
#%%------------------------------------------------------------------------------------------------
ops.wipe()
ops.model('basic', '-ndm', 1, '-ndf', 1)

# Nodes
ops.node(0, 0.0)
ops.node(1, 0.0)
ops.node(2, 0.0)

ops.fix(1, 1)

ops.mass(2, mass)
# Materials
# 1. Structural Spring
ops.uniaxialMaterial('Elastic', 1, k_structure)

# 2. Friction Damper (Modeled as Elasto-Plastic)
# Steel01: matTag, yieldForce, k_initial, post_yield_stiffness_ratio
k_damper_initial = 1e6 # High stiffness to represent "stiction"
ops.uniaxialMaterial('Steel01', 2, F_friction, k_damper_initial, 0.0001)

# Elements
ops.element('zeroLength', 1, 1, 2, '-mat', 1, '-dir', 1) # Spring
ops.element('zeroLength', 2, 1, 2, '-mat', 2, '-dir', 1) # Damper

ops.uniaxialMaterial('CoulombDamper', 11, 0.0, F_friction)
 
ops.element('zeroLength', 3, 0, 2,'-mat', 11,'-dir', 1)

# Excitation
dt = 0.01
n_steps = 2000
time = np.arange(0, n_steps*dt, dt)
accel = 9.81 * 0.4 * np.sin(2*np.pi*1.5*time) # 0.4g, 1.5 Hz
ops.timeSeries('Path', 1, '-dt', dt, '-values', *accel.tolist())
ops.pattern('UniformExcitation', 1, 1, '-accel', 1)

# Analysis Setup
ops.system('BandGeneral')
ops.numberer('Plain')
ops.constraints('Plain')
ops.integrator('Newmark', 0.5, 0.25)
ops.algorithm('Newton')
ops.analysis('Transient')

# Record data
disp_list = []
force_damper_list = []
Coulomb_Damper_list = []

for i in range(n_steps):
    OK = ops.analyze(1, dt)
    disp_list.append(ops.nodeDisp(2, 1))
    # Get the basic force of element 2 (damper)
    force_damper_list.append(ops.eleResponse(2, 'force')[0])  # Get force from element 2
    Coulomb_Damper_list.append(ops.eleResponse(3, 'force')[0])# Get force from element 3
    print(i+1, disp_list[-1], force_damper_list[-1], Coulomb_Damper_list[-1])
#%%------------------------------------------------------------------------------------------------
# Plotting
plt.figure(figsize=(10, 10))
plt.subplot(4,1,1); plt.plot(time, disp_list); plt.title('Displacement (m)')
plt.subplot(4,1,2); plt.plot(time, force_damper_list, 'r'); plt.title('Damper Force (kN)')
plt.subplot(4,1,3); plt.plot(disp_list, Coulomb_Damper_list, 'black'); plt.title('Coulomb Damper Force (kN) vs Displacement Hysteresis Loop (m)')
plt.subplot(4,1,4); plt.plot(disp_list, force_damper_list, 'g'); plt.title('Damper Force (kN) vs Displacement (m) Hysteresis Loop')
plt.tight_layout()
plt.show()
#%%------------------------------------------------------------------------------------------------