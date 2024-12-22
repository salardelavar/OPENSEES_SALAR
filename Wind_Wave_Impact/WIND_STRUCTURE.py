###########################################################################################################
#                                                 IN THE NAME OF ALLAH                                    #
#---------------------------------------------------------------------------------------------------------# 
#                                   THIS PROGRAM IS WRITTEN BY SALAR DELAVAR GHASHGHAEI                   #
#                                          EMAIL: SALAR.D.GHASHGHAEI@GMAIL.COM                            #
###########################################################################################################
#
# Title:  
# Dynamic Response Analysis of a Single-Degree-of-Freedom (SDOF) System Subjected to Wind Impact Loading Using OpenSeesPy
# Generate a combined time series of turbulent wind pressure, gust loads, and random fluctuations.
#
# Target:  
# This code models and analyzes the dynamic response of an SDOF structural system under time-varying wind impact pressure.
# The objectives include:  
#
# 1. Wind Pressure Time Series Generation:  
#    - Implements a flexible waveform generator for wind pressure loading, supporting sine, cosine, combined sine+cosine, and turbulent fluctuations.  
#    - Visualizes the generated wind impact pressure time series.
#
# 2. Structural Model Development:  
#    - Constructs an SDOF system in OpenSeesPy with specified stiffness (`k`), mass (`m`), and damping ratio.  
#    - Applies the wind pressure as a dynamic load to the system.  
#
# 3. Dynamic Analysis:  
#    - Conducts transient analysis using the Newmark integration method.  
#    - Monitors key responses, including displacement, velocity, acceleration, and base reaction forces over time.  
#
# 4. Visualization:  
#    - Provides time-history plots for displacement, velocity, acceleration, and base reaction forces to understand the structural response to wind impact loading.
#
# This simulation aids engineers and researchers in evaluating the structural resilience of systems subjected to varying wind impact pressures,
# applicable in contexts such as wind-induced vibrations, flutter analysis, and high-rise building performance.


import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt

# Define structure parameters
k = 2.0e6  # Stiffness of the structure (N/m)
m = 1500.0  # Mass of the structure (kg)
damping_ratio = 0.02  # Damping ratio

# Define wind impact pressure parameters
impact_duration = 5.0  # Total duration of wind impact (s)
duration = 15.0  # Total simulation duration (s)
dt = 0.01  # Time step (s)

def generate_wind_pressure_time_series(
    impact_duration, dt, wave_type='sin+cos', frequency_factor=1.0, amplitude=500.0, turbulence_intensity=100.0):
    """
    Generate a combined time series of wind pressure including gusts and turbulence.

    Parameters:
    - impact_duration: Duration of the impact event (seconds).
    - dt: Time step for the time series (seconds).
    - wave_type: Type of waveform ('sine', 'cos', 'sin+cos'). Default is 'sin+cos'.
    - frequency_factor: Factor to modify the frequency of the waveform. Default is 1.0.
    - amplitude: Amplitude of the waveform. Default is 500.0.
    - turbulence_intensity: Intensity of the random turbulence. Default is 100.0.

    Returns:
    - time: Array of time values.
    - pressure: Array of composite pressure values.
    """
    time = np.arange(0, impact_duration, dt)
    frequency = frequency_factor / impact_duration  # Adjust frequency based on the impact duration

    # Base waveform
    if wave_type == 'sine':
        base_pressure = np.sin(2 * np.pi * frequency * time) * amplitude
    elif wave_type == 'cos':
        base_pressure = np.cos(2 * np.pi * frequency * time) * amplitude
    elif wave_type == 'sin+cos':
        base_pressure = (np.sin(2 * np.pi * frequency * time) + np.cos(1.5 * np.pi * frequency * time)) * amplitude
    else:
        raise ValueError("Unsupported wave type. Choose from 'sine', 'cos', or 'sin+cos'.")

    # Gust-induced load
    gust_pressure = 0.8 * amplitude * np.sin(0.3 * np.pi * time / impact_duration)

    # Random turbulent fluctuations
    turbulent_pressure = np.random.normal(0, turbulence_intensity, size=time.shape)

    # Better simulation of wind pressure with exponential decay factor for gusts
    decay_factor = np.exp(-0.5 * time / impact_duration)
    gust_pressure_with_decay = gust_pressure * decay_factor

    # Combine all loads
    total_pressure = base_pressure + gust_pressure_with_decay + turbulent_pressure

    return time, total_pressure

# Generate composite loading
time_series, wind_pressure = generate_wind_pressure_time_series(
    impact_duration, dt, wave_type='sin+cos', frequency_factor=0.5, amplitude=500.0, turbulence_intensity=200.0)

# Plot the composite loading
def plot_wind_pressure_time_series(time, pressure):
    # Plot the pressure time history
    plt.figure(figsize=(10, 6))
    plt.plot(time, pressure, label='Wind Pressure', color='purple')
    plt.xlabel('Time (s)')
    plt.ylabel('Pressure Force (N)')
    plt.title('Wind Pressure Loading')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

plot_wind_pressure_time_series(time_series, wind_pressure)

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
# Define materials for structural damper
ops.uniaxialMaterial('Elastic', 2, 0.0, CS)

# Define element
ops.element('zeroLength', 1, 1, 2, '-mat', 1, 2, '-dir', 1, 1)  # DOF[1] LATERAL SPRING

# Apply time-dependent wind loading
time_series_tag = 1
pattern_tag = 1
# Replace the wind pressure with composite loading in OpenSees
ops.timeSeries('Path', time_series_tag, '-dt', dt, '-values', *wind_pressure)
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

# Re-run dynamic analysis with new loading
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
plt.ylabel('Acceleration [m/s²]')
plt.title('Acceleration Time History')
plt.grid(True)

# Base Reaction Force
plt.subplot(4, 1, 4)
plt.plot(time, base_reaction, label='Base Reaction', color='red')
plt.xlabel('Time [s]')
plt.ylabel('Base Reaction [N]')
plt.title('Base Reaction Time History')
plt.grid(True)

plt.tight_layout()
plt.show()
