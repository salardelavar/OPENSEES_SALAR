#          #####################################################################################
#          #                                  IN THE NAME OF ALLAH                             #
#          #                           FREE VIBRATION ANAYSIS OF ELASTIC SPRING                #
#          #-----------------------------------------------------------------------------------#
#          #              THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)           #
#          #                       EMAIL: salar.d.ghashghaei@gmail.com                         #
#          #####################################################################################

import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt

# Set up variables
e = 210.0e3  # Modulus of elasticity (Young's modulus)
A = 0.25  # Cross-sectional area of the element
l = 1000.0  # Length of the element
k = (e*A)/l # Truss Axial Stiffness
m = 50.0  # Mass of the node
u0 = 0.001  # Initial displacement applied to the node
damping_ratio = 0.02  # Damping ratio

# Set analysis parameters
duration = 50.0 # [s] 50 Seconds
dt = 0.01 # time step

# Function to perform transient analysis
def perform_analysis(damping=False):
    # Set up the model
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)

    # Define nodes
    ops.node(1, 0, 0)
    ops.node(2, l, 0)

    # Define boundary conditions
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 0, 1, 1)

    # Define mass
    ops.mass(2, m, 0, 0)

    # Define material
    ops.uniaxialMaterial('Elastic', 1, e)

    # Define element
    ops.element('Truss', 1, 1, 2, A, 1)

    # Static analysis to apply initial displacement
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(2, 1.0, 0, 0)

    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.system('BandGeneral')
    ops.algorithm('Newton')
    ops.test('NormDispIncr', 1.0e-8, 10)
    ops.integrator('DisplacementControl', 2, 1, u0)
    ops.analysis('Static')
    ops.analyze(1)

    ops.setTime(0.0)

    # Wipe analysis and reset time
    ops.wipeAnalysis()
    ops.remove('loadPattern', 1)
    ops.system('UmfPack')

    # Dynamic analysis
    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.system('UmfPack')
    ops.test('NormDispIncr', 1.0e-8, 10)
    ops.integrator('Newmark', 0.5, 0.25)
    ops.algorithm('Newton')

    if damping:
        # Calculate Rayleigh damping factors
        omega1 = np.sqrt(k / m)
        a0 = (2 * damping_ratio * omega1) / omega1
        a1 = (2 * damping_ratio) / omega1
        # Apply Rayleigh damping
        ops.rayleigh(a0, a1, 0, 0)

    ops.analysis('Transient')

    # Perform transient analysis and store results
    time = []
    displacement = []
    velocity = []
    acceleration = []
    spring_force = []

    stable = 0
    current_time = 0.0

    while stable == 0 and current_time < duration:
        stable = ops.analyze(1, dt)
        current_time = ops.getTime()
        time.append(current_time)
        displacement.append(ops.nodeDisp(2, 1))
        velocity.append(ops.nodeVel(2, 1))
        acceleration.append(ops.nodeAccel(2, 1))
        spring_force.append(-ops.eleResponse(1, 'force')[0])

    return time, displacement, velocity, acceleration, spring_force

# Perform analyses
time_undamped, displacement_undamped, velocity_undamped, acceleration_undamped, spring_force_undamped = perform_analysis(damping=False)
time_damped, displacement_damped, velocity_damped, acceleration_damped, spring_force_damped = perform_analysis(damping=True)

### PLOT THE TIME HISTORY:
PLOT_4_CHART()
