###########################################################################################################
#                                                 IN THE NAME OF ALLAH                                    #
#---------------------------------------------------------------------------------------------------------#
#                                   THIS PROGRAM IS WRITTEN BY SALAR DELAVAR GHASHGHAEI                   #
#                                          EMAIL: SALAR.D.GHASHGHAEI@GMAIL.COM                            #
###########################################################################################################
# Title:
# Dynamic Analysis of a Structural System Under Combined Water Impact and Seismic Loading Using OpenSeesPy
# 
# Target:
# This code aims to simulate and analyze the dynamic response of a simple structural system subjected to:
# 
# 1. Water Impact Loading: Generated using a sinusoidal and cosine-based pressure time series to model water wave impact.
# 2. Seismic Loading: Applied as ground accelerations in both X and Y directions using predefined seismic time histories.
# 
# The analysis is conducted in OpenSeesPy, a finite element framework, to compute the structure's displacement, velocity, acceleration, 
# and base reactions over time. The results are visualized using Matplotlib to interpret the system's dynamic behavior under these combined loads.
# This setup is particularly relevant for designing and evaluating marine or coastal structures exposed to both hydrodynamic and seismic forces.
import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt

# Define structure parameters
k = 1.0e6  # Stiffness of the structure (N/m)
m = 1000.0  # Mass of the structure (kg)
damping_ratio = 0.05  # Damping ratio

# Define water impact pressure parameters
impact_duration = 2.0  # Total duration of water impact (s)
duration = 10.0  # Total simulation duration (s)
dt = 0.01  # Time step (s)

# Define seismic acceleration parameters
seismic_x_acceleration = np.sin(np.linspace(0, 2 * np.pi, int(duration / dt))) * 10.3  # Example X-direction acceleration
seismic_y_acceleration = np.sin(np.linspace(0, 2 * np.pi, int(duration / dt))) * 5.2  # Example Y-direction acceleration

# Water impact pressure time series generator
def water_impact_time_series(impact_duration, dt, wave_type='sin+cos', frequency_factor=1.0, amplitude=1e4):
    time = np.arange(0, impact_duration, dt)
    frequency = frequency_factor / impact_duration
    if wave_type == 'sin+cos':
        pressure = (np.sin(5 * np.pi * frequency * time) + np.cos(2 * np.pi * frequency * time)) * amplitude
    else:
        raise ValueError("Unsupported wave type. Choose 'sin+cos'.")
    return time, pressure

# Generate water impact pressure time series
time_series, water_pressure = water_impact_time_series(impact_duration, dt)

# Plot water impact loading
def plot_water_impact_time_series(time, pressure):
    plt.figure(figsize=(10, 6))
    plt.plot(time, pressure, label='Water Wave Pressure', color='blue')
    plt.xlabel('Time (s)')
    plt.ylabel('Pressure Force (N)')
    plt.title('Water Wave Pressure Loading')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

plot_water_impact_time_series(time_series, water_pressure)

# Initialize OpenSees model
ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 2)  # 2D model with 2 DOFs per node

# Define nodes
ops.node(1, 0.0, 0.0)  # Fixed base
ops.node(2, 0.0, 0.0)  # Mass node

# Define boundary conditions
ops.fix(1, 1, 1)  # Fixed in both X and Y

# Define mass
ops.mass(2, m, m)

# Define material properties
ops.uniaxialMaterial('Elastic', 1, k)

# Define element
ops.element('zeroLength', 1, 1, 2, '-mat', 1, '-dir', 1)  # Lateral spring in X direction
ops.element('zeroLength', 2, 1, 2, '-mat', 1, '-dir', 2)  # Lateral spring in Y direction

# Apply water impact loading
water_time_series_tag = 1
water_pattern_tag = 1
ops.timeSeries('Path', water_time_series_tag, '-dt', dt, '-values', *water_pressure)
ops.pattern('Plain', water_pattern_tag, water_time_series_tag)
ops.load(2, 1.0, 0.0)  # Load applied to the mass node in the X direction

# Apply seismic accelerations
seismic_x_time_series_tag = 2
seismic_y_time_series_tag = 3
ops.timeSeries('Path', seismic_x_time_series_tag, '-dt', dt, '-values', *seismic_x_acceleration)
ops.timeSeries('Path', seismic_y_time_series_tag, '-dt', dt, '-values', *seismic_y_acceleration)

ops.pattern('UniformExcitation', 2, 1, '-accel', seismic_x_time_series_tag)  # Seismic in X
ops.pattern('UniformExcitation', 3, 2, '-accel', seismic_y_time_series_tag)  # Seismic in Y

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
    base_reaction_x.append(-k * displacement_x[-1])  # Reaction force in x-direction
    base_reaction_y.append(-k * displacement_y[-1])  # Reaction force in y-direction
    #base_reaction_x.append(ops.nodeResponse(1, 1, 6))  # Reaction force in x-direction
    #base_reaction_y.append(ops.nodeResponse(1, 2, 6))  # Reaction force in y-direction

# Plot results
plt.figure(figsize=(18, 28))

# Displacement
plt.subplot(8, 1, 1)
plt.plot(time, displacement_x, label='Displacement X', color='blue')
plt.xlabel('Time [s]')
plt.ylabel('Displacement X [m]')
plt.title('Displacement Time History')
plt.legend()
plt.grid(True)

plt.subplot(8, 1, 2)
plt.plot(time, displacement_y, label='Displacement Y', color='green')
plt.xlabel('Time [s]')
plt.ylabel('Displacement Y [m]')
plt.title('Displacement Time History')
plt.legend()
plt.grid(True)

# Velocity
plt.subplot(8, 1, 3)
plt.plot(time, velocity_x, label='Velocity X', color='orange')
plt.xlabel('Time [s]')
plt.ylabel('Velocity X [m/s]')
plt.title('Velocity Time History')
plt.legend()
plt.grid(True)

plt.subplot(8, 1, 4)
plt.plot(time, velocity_y, label='Velocity Y', color='purple')
plt.xlabel('Time [s]')
plt.ylabel('Velocity Y [m/s]')
plt.title('Velocity Time History')
plt.legend()
plt.grid(True)

# Acceleration
plt.subplot(8, 1, 5)
plt.plot(time, acceleration_x, label='Acceleration X', color='red')
plt.xlabel('Time [s]')
plt.ylabel('Acceleration X [m/s^2]')
plt.title('Acceleration Time History')
plt.legend()
plt.grid(True)

plt.subplot(8, 1, 6)
plt.plot(time, acceleration_y, label='Acceleration Y', color='brown')
plt.xlabel('Time [s]')
plt.ylabel('Acceleration Y [m/s^2]')
plt.title('Acceleration Time History')
plt.legend()
plt.grid(True)

# Base Reaction Displacement
plt.subplot(8, 1, 7)
plt.plot(time, base_reaction_x, label='Base Reaction X', color='red')
plt.xlabel('Time [s]')
plt.ylabel('Base Reaction X [N]')
plt.title('Base Reaction Time History')
plt.legend()
plt.grid(True)

plt.subplot(8, 1, 8)
plt.plot(time, base_reaction_x, label='Base Reaction Y', color='red')
plt.xlabel('Time [s]')
plt.ylabel('Base Reaction Y [N]')
plt.title('Base Reaction Time History')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()
