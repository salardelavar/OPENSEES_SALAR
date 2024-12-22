###########################################################################################################
#                                                 IN THE NAME OF ALLAH                                    #
#---------------------------------------------------------------------------------------------------------#
#                                   THIS PROGRAM IS WRITTEN BY SALAR DELAVAR GHASHGHAEI                   #
#                                          EMAIL: SALAR.D.GHASHGHAEI@GMAIL.COM                            #
###########################################################################################################
#
# Title:
# Dynamic Response Analysis of a Single-Degree-of-Freedom (SDOF) System Subjected to Explosion Loading Using OpenSeesPy
#
# Target:
# This code simulates the dynamic response of a single-degree-of-freedom (SDOF) structural system subjected to explosion-induced loading.
# It achieves the following:
#
# 1. Define Explosion Loading:
#    - Implements the Friedlander equation to model pressure-time history from an explosive event.
#    - Plots the explosion pressure profile over time.
#
# 2. Structural Model:
#    - Models the SDOF system in OpenSeesPy with a linear elastic spring (stiffness `k`), mass `m`, and damping ratio.
#    - Applies the explosion-induced time-dependent loading to the mass node.
#
# 3. Dynamic Analysis:
#    - Simulates the time history response using the Newmark method for transient analysis.
#    - Tracks system responses, including displacement, velocity, acceleration, and base reactions.
#
# 4. Visualization:
#    - Generates plots for displacement, velocity, acceleration, and base reactions to evaluate the impact of the explosion loading on the structure.
#
# This simulation is useful for structural engineers studying the effects of blast loads on structures, aiding in the design and assessment of resilient systems.
import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt

# Define Friedlander equation for explosion loading
def friedlander(t, P0, t0, A):
    """
    Calculate the pressure at time t using the Friedlander equation.

    Parameters:
    t (float): Time at which to calculate the pressure.
    P0 (float): Peak pressure.
    t0 (float): Positive phase duration.
    A (float): Wave decay coefficient.

    Returns:
    float: Pressure at time t.
    """
    if t < 0:
        return 0
    return P0 * (1 - t / t0) * np.exp(-A * t / t0)

# Define structure parameters
k = 1.0e6  # Stiffness of the structure (N/m)
m = 1000.0  # Mass of the structure (kg)
damping_ratio = 0.05  # Damping ratio

# Define explosion loading parameters
P0 = 1e5  # Peak pressure (Pa)
t0 = 0.1  # Positive phase duration (s)
A = 1.3   # Wave decay coefficient

dt = 0.001  # Time step (s)
impact_duration = 2.0  # Total duration of explosion impact (s)loading 
duration = 10.0  # Total simulation duration (s)

def explosion_time_series(impact_duration, dt, P0, t0, A):
    time = np.arange(0, impact_duration, dt)
    pressure = [friedlander(t, P0, t0, A) for t in time]
    return time, pressure

# Generate explosion pressure time series
time_series, explosion_pressure = explosion_time_series(impact_duration, dt, P0, t0, A)

def plot_friedlander(P0, t0, A, dt, impact_duration):
    """
    Plot the Friedlander explosion loading formula over a given duration.

    Parameters:
    P0 (float): Peak pressure (Pa).
    t0 (float): Positive phase duration (s).
    A (float): Wave decay coefficient.
    dt (float): Time step (s).
    duration (float): Total duration for the plot (s).
    """
    time = np.arange(0, impact_duration, dt)
    pressure = [friedlander(t, P0, t0, A) for t in time]

    # Plot the pressure time history
    plt.figure(figsize=(10, 6))
    plt.plot(time, pressure, label='Explosion Pressure', color='blue')
    plt.xlabel('Time (s)')
    plt.ylabel('Pressure (Pa)')
    plt.title('Friedlander Explosion Loading')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()
    

# Plot the Friedlander explosion loading
plot_friedlander(P0, t0, A, dt, impact_duration)

# Initialize OpenSees model
ops.wipe()
ops.model('basic', '-ndm', 1, '-ndf', 1)

# Define nodes
ops.node(1, 0.0)  # Fixed base
ops.node(2, 0.0)  # Mass node

# Define boundary conditions
ops.fix(1, 1)

# Define mass
ops.mass(2, m)

# Natural frequency (rad/s)
wn = (k / m) ** 0.5
# Damping coefficient (Ns/m)
CS = 2 * wn * m * damping_ratio

# Define material properties
ops.uniaxialMaterial('Elastic', 1, k)
# Define materials for structure damper
ops.uniaxialMaterial('Elastic', 2, 0.0, CS)

# Define element
#ops.element('Truss', 1, 1, 2, 1.0, 1)
ops.element('zeroLength', 1, 1, 2, '-mat', 1, 2, '-dir', 1, 1) # DOF[1] LATERAL SPRING 

# Apply time-dependent explosion loading
time_series_tag = 1
pattern_tag = 1
ops.timeSeries('Path', time_series_tag, '-dt', dt, '-values', *explosion_pressure)
ops.pattern('Plain', pattern_tag, time_series_tag)
ops.load(2, 1.0)  # Load applied to the mass node

# Set analysis parameters
ops.constraints('Plain')
ops.numberer('Plain')
ops.system('BandGeneral')
ops.test('NormDispIncr', 1.0e-6, 10)
ops.algorithm('Newton')
ops.integrator('Newmark', 0.5, 0.25)
ops.analysis('Transient')

# Perform dynamic analysis
time = []
displacement = []
velocity = []
acceleration = []
base_reaction = []

stable = 0
current_time = 0.0
    
while stable == 0 and current_time < duration:
    stable = ops.analyze(1, dt)
    current_time = ops.getTime()
    time.append(current_time)
    displacement.append(ops.nodeDisp(2, 1))
    velocity.append(ops.nodeVel(2, 1))
    acceleration.append(ops.nodeAccel(2, 1))
    base_reaction.append(-k * displacement[-1])  # Reaction force
    #base_reaction.append(ops.nodeResponse(1, 1, 6))  # Reaction force

# Plot results
plt.figure(figsize=(18, 16))

# Displacement
plt.subplot(4, 1, 1)
plt.plot(time, displacement, label='Displacement')
plt.xlabel('Time [s]')
plt.ylabel('Displacement [m]')
plt.title('Displacement Time History')
plt.grid(True)

# Velocity
plt.subplot(4, 1, 2)
plt.plot(time, velocity, label='Velocity', color='orange')
plt.xlabel('Time [s]')
plt.ylabel('Velocity [m/s]')
plt.title('Velocity Time History')
plt.grid(True)

# Acceleration
plt.subplot(4, 1, 3)
plt.plot(time, acceleration, label='Acceleration', color='green')
plt.xlabel('Time [s]')
plt.ylabel('Acceleration [m/s^2]')
plt.title('Acceleration Time History')
plt.grid(True)

# Base Reaction Displacement
plt.subplot(4, 1, 4)
plt.plot(time, base_reaction, label='Base Reaction', color='red')
plt.xlabel('Time [s]')
plt.ylabel('Base Reaction [N]')
plt.title('Base Reaction Time History')
plt.grid(True)

plt.tight_layout()
plt.show()
