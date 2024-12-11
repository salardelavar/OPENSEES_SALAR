"""
### IN THE NAME OF ALLAH ###

Title:  
Dynamic Analysis of a 2D Structural System Under Explosion Loading Using the Friedlander Equation and OpenSeesPy

Target:  
This code aims to simulate and analyze the dynamic response of a two-dimensional structural system subjected to explosion-induced pressure loading.
 The specific objectives include:  

1. Explosion Pressure Modeling:  
   - Implements the Friedlander equation to model pressure loading from an explosion.  
   - Provides adjustable parameters for peak pressure, positive phase duration, and decay coefficient.  
   - Generates and visualizes the time-dependent pressure history.

2. Structural Model Setup:  
   - Models a 2D structural system using OpenSeesPy with stiffness in both x and y directions.  
   - Defines mass, damping, and spring properties for the system.  

3. Dynamic Loading and Analysis:  
   - Applies the generated explosion pressure as a dynamic load in both x and y directions.  
   - Conducts transient analysis using Newmark integration to compute system responses over time.

4. Response Visualization:  
   - Plots time histories of key responses, including:  
     - Displacement in x and y directions.  
     - Velocity in x and y directions.  
     - Acceleration in x and y directions.  
     - Base reactions in x and y directions.  

This simulation is useful for evaluating the behavior of structural systems under blast loads, aiding in the design and assessment of blast-resistant structures.

THIS PROGRAM IS WRITTEN BY SALAR DELAVAR GHASHGHAEI 
"""
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
k_x = 1.0e6  # Stiffness in x-direction (N/m)
k_y = 1.0e6  # Stiffness in y-direction (N/m)
m = 1000.0  # Mass of the structure (kg)
damping_ratio = 0.05  # Damping ratio

# Define explosion loading parameters
P0 = 1e5  # Peak pressure (Pa)
t0 = 0.1  # Positive phase duration (s)
A = 1.3  # Wave decay coefficient

dt = 0.001  # Time step (s)
impact_duration = 2.0  # Total duration of water impact (s)
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
ops.model('basic', '-ndm', 2, '-ndf', 2)

# Define nodes
ops.node(1, 0.0, 0.0)  # Fixed base
ops.node(2, 0.0, 0.0)  # Mass node

# Define boundary conditions
ops.fix(1, 1, 1)  # Fix both x and y directions at the base

# Define mass
ops.mass(2, m, m)  # Mass in both directions

# Define material properties
ops.uniaxialMaterial('Elastic', 1, k_x)  # Spring in x-direction
ops.uniaxialMaterial('Elastic', 2, k_y)  # Spring in y-direction

# Define elements
ops.element('zeroLength', 1, 1, 2, '-mat', 1, 2, '-dir', 1, 2)  # Springs in x and y

# Apply time-dependent explosion loading
time_series_tag = 1
pattern_tag = 1
ops.timeSeries('Path', time_series_tag, '-dt', dt, '-values', *explosion_pressure)
ops.pattern('Plain', pattern_tag, time_series_tag)
ops.load(2, 1.0, 1.0)  # Load applied in both x and y directions

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
displacement_x = []
displacement_y = []
velocity_x = []
velocity_y = []
acceleration_x = []
acceleration_y = []
base_reaction_x = []
base_reaction_y = []

stable = 0
current_time = 0.0

while stable == 0 and current_time < duration:
    stable = ops.analyze(1, dt)
    current_time = ops.getTime()
    time.append(current_time)
    displacement_x.append(ops.nodeDisp(2, 1))
    displacement_y.append(ops.nodeDisp(2, 2))
    velocity_x.append(ops.nodeVel(2, 1))
    velocity_y.append(ops.nodeVel(2, 2))
    acceleration_x.append(ops.nodeAccel(2, 1))
    acceleration_y.append(ops.nodeAccel(2, 2))
    base_reaction_x.append(ops.nodeResponse(1, 1, 6))  # Reaction force in x-direction
    base_reaction_y.append(ops.nodeResponse(1, 2, 6))  # Reaction force in y-direction

# Plot results
plt.figure(figsize=(18, 16))

# Displacement
plt.subplot(4, 2, 1)
plt.plot(time, displacement_x, label='Displacement X')
plt.xlabel('Time [s]')
plt.ylabel('Displacement [m]')
plt.title('Displacement Time History (X-direction)')
plt.grid(True)

plt.subplot(4, 2, 2)
plt.plot(time, displacement_y, label='Displacement Y')
plt.xlabel('Time [s]')
plt.ylabel('Displacement [m]')
plt.title('Displacement Time History (Y-direction)')
plt.grid(True)

# Velocity
plt.subplot(4, 2, 3)
plt.plot(time, velocity_x, label='Velocity X', color='orange')
plt.xlabel('Time [s]')
plt.ylabel('Velocity [m/s]')
plt.title('Velocity Time History (X-direction)')
plt.grid(True)

plt.subplot(4, 2, 4)
plt.plot(time, velocity_y, label='Velocity Y', color='purple')
plt.xlabel('Time [s]')
plt.ylabel('Velocity [m/s]')
plt.title('Velocity Time History (Y-direction)')
plt.grid(True)

# Acceleration
plt.subplot(4, 2, 5)
plt.plot(time, acceleration_x, label='Acceleration X', color='green')
plt.xlabel('Time [s]')
plt.ylabel('Acceleration [m/s^2]')
plt.title('Acceleration Time History (X-direction)')
plt.grid(True)

plt.subplot(4, 2, 6)
plt.plot(time, acceleration_y, label='Acceleration Y', color='red')
plt.xlabel('Time [s]')
plt.ylabel('Acceleration [m/s^2]')
plt.title('Acceleration Time History (Y-direction)')
plt.grid(True)

# Base Reaction
plt.subplot(4, 2, 7)
plt.plot(time, base_reaction_x, label='Base Reaction X', color='cyan')
plt.xlabel('Time [s]')
plt.ylabel('Base Reaction [N]')
plt.title('Base Reaction Time History (X-direction)')
plt.grid(True)

plt.subplot(4, 2, 8)
plt.plot(time, base_reaction_y, label='Base Reaction Y', color='magenta')
plt.xlabel('Time [s]')
plt.ylabel('Base Reaction [N]')
plt.title('Base Reaction Time History (Y-direction)')
plt.grid(True)

plt.tight_layout()
plt.show()
