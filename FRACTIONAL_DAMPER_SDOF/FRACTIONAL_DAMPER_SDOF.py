###########################################################################################################
#                                         IN THE NAME OF ALLAH                                            #
#      DYNAMIC RESPONSE ANALYSIS OF A SINGLE-DEGREE-OF-FREEDOM (SDOF) SYSTEM WITH FRACTIONAL DAMPING      #
#                               AND ELASTIC STRUCTURE UNDER SEISMIC ACCELERATION                          #
#---------------------------------------------------------------------------------------------------------#
#                          THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                     #
#                                   EMAIL: salar.d.ghashghaei@gmail.com                                   #
###########################################################################################################
import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
#------------------------------------------------------------------------------------------------
# Define parameters (units: m, kN)
# Define units
m = 1.0      # mass in kg
kN = 1.0     # force in kN
m_sec = 1.0  # length in meters, time in seconds
#------------------------------------------------------------------------------------------------
# Clear OpenSees model
ops.wipe()
#------------------------------------------------------------------------------------------------
# Model parameters
mass = 1.0 * m                  # Mass of the SDOF system
k = 10.0 * kN/m_sec             # Initial stiffness
damping_ratio = 0.05            # Damping ratio
Omega = (k/mass)**0.5
C = (damping_ratio * 2) / Omega # Viscous damping coefficient (approximating fractional damping effect)
#------------------------------------------------------------------------------------------------
# Hysteretic spring parameters (HystereticSM model)
s1p = 0.5 * kN         # Positive yield strength
e1p = 0.02 * m_sec     # Positive yield displacement
s2p = 0.8 * kN         # Second stiffness point strength
e2p = 0.05 * m_sec     # Second stiffness point displacement
s3p = 1.0 * kN         # Third stiffness point strength
e3p = 0.10 * m_sec     # Third stiffness point displacement
s1n, e1n = -s1p, -e1p  # Negative counterparts (symmetric)
s2n, e2n = -s2p, -e2p
s3n, e3n = -s3p, -e3p
pinchX = 0.8           # Pinching factor in X direction
pinchY = 0.5           # Pinching factor in Y direction
damage1 = 0.0          # Damage due to ductility
damage2 = 0.0          # Damage due to energy
beta = 0.1             # Stiffness degradation parameter

# Stiffness 
# Define hysteretic material
ops.uniaxialMaterial('Hysteretic', 1, s1p, e1p, s2p, e2p, s3p, e3p,
                     s1n, e1n, s2n, e2n, s3n, e3n, pinchX, pinchY, damage1, damage2, beta)
# INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Hysteretic_Material

# Damper
ops.uniaxialMaterial('Viscous', 2, C, 20.0)  # Material for C (alpha=1.0 for linear) 
# INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Viscous_Material                    
#------------------------------------------------------------------------------------------------
# Define the model (1D, 2 nodes, 1 DOF per node)
ops.model('basic', '-ndm', 1, '-ndf', 1)

# Nodes
ops.node(1, 0.0)  # Fixed base
ops.node(2, 0.0)  # Mass node
#------------------------------------------------------------------------------------------------
# Fixity
ops.fix(1, 1)     # Fixed at base
#------------------------------------------------------------------------------------------------
# Mass
ops.mass(2, mass)
#------------------------------------------------------------------------------------------------
# Define element (zeroLength element with hysteretic material)
ops.element('zeroLength', 1, 1, 2, '-mat', 1, 2, '-dir', 1, 1)
#------------------------------------------------------------------------------------------------
# Define damping
ops.rayleigh(C/mass, 0.0, 0.0, 0.0)  # Mass-proportional damping
#------------------------------------------------------------------------------------------------
# Ground motion (sine wave for simplicity)
time = np.linspace(0, 10.0, 1000)  # 10 seconds, 1000 points
freq = 1.0  # Hz
accel = 5.9 * np.sin(2 * np.pi * freq * time)  # Acceleration in m/s^2
#------------------------------------------------------------------------------------------------
# Time series for ground motion
ops.timeSeries('Path', 1, '-dt', time[1]-time[0], '-values', *accel)
#------------------------------------------------------------------------------------------------
# Pattern
ops.pattern('UniformExcitation', 1, 1, '-accel', 1)

# Analysis setup
ops.wipeAnalysis()
ops.constraints('Plain')
ops.numberer('Plain')
ops.system('BandGeneral')
ops.test('NormUnbalance', 1e-6, 10)
ops.algorithm('Newton')
ops.integrator('Newmark', 0.5, 0.25)  # Newmark-beta method
ops.analysis('Transient')
#------------------------------------------------------------------------------------------------
# Record results
disp = []
force = []
time_list = []

for i in range(len(time)):
    ops.analyze(1, time[1]-time[0])
    disp.append(ops.nodeDisp(2, 1))
    force.append(ops.eleForce(1, 1))
    time_list.append(ops.getTime())
#------------------------------------------------------------------------------------------------
# Plotting
plt.figure(figsize=(12, 5))

# Time history of displacement
plt.subplot(1, 2, 1)
plt.plot(time_list, disp, label='Displacement', color='purple')
plt.xlabel('Time (s)')
plt.ylabel('Displacement (m)')
plt.title('Displacement Time History')
plt.grid(True)
plt.legend()

# Hysteresis loop
plt.subplot(1, 2, 2)
plt.plot(disp, force, label='Hysteretic Response', color='black')
plt.xlabel('Displacement (m)')
plt.ylabel('Force (kN)')
plt.title('Force-Displacement Hysteresis')
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.show()
#------------------------------------------------------------------------------------------------