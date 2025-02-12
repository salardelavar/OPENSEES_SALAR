"""
   #######################################################################
   #                            IN THE NAME OF ALLAH                     #
   #     Using Combined Kanai-Tajimi model and Increasing Sinusoidal     #         
   #        Ground Acceleration involves simulating a non-stationary     #
   # stochastic process to represent earthquake ground motion            #
   # With Uncertainty Conditions Using Probability Distribution Function #
   # Monte-Carlo Method                                                  #
   #---------------------------------------------------------------------#
   #            THIS PROGRAM WRITTEN BY SALAR DELAVAR QASHQAI            #
   #                 EMAIL: salar.d.ghashghaei@gmail.com                 #
   #######################################################################

Generate ground acceleration using the Kanai-Tajimi model involves simulating a non-stationary stochastic process:

The Kanai-Tajimi model is a widely used approach to simulate the random ground acceleration experienced during an earthquake.
It models the ground acceleration as a filtered Gaussian process, where a white noise process (representing bedrock acceleration)
is filtered through a viscoelastic system that represents the soil layer. This model is particularly useful for analyzing structures
subjected to seismic events.

The classical Kanai-Tajimi model consists of a linear oscillator attached to the bedrock, which moves with an acceleration modeled
as a Gaussian white noise process. The model parameters, such as the damping ratio (ζ_g) and the natural frequency (ω_g) of the soil
deposit, are calibrated based on historical earthquake data, including zero crossings and other statistical properties.

In advanced versions, such as the fractional Kanai-Tajimi model, the purely viscous element in the Kelvin-Voigt element (a spring in
parallel with a dashpot) is replaced with a springpot. This modification allows for a more accurate representation of the viscoelastic
behavior of the ground, ensuring that the number of zero crossings of the absolute acceleration at the free field remains finite.

The non-stationary nature of earthquake ground motion is typically accounted for by applying an envelope or shape function to the
stationary process. This adjusts the intensity of the ground motion over time, reflecting the way an actual earthquake's intensity
changes during the event.

This program implements the Kanai-Tajimi model to generate artificial earthquake records. It combines the stochastic ground motion
from the Kanai-Tajimi model with a deterministic increasing sinusoidal function to simulate more complex ground motion scenarios.
The resulting ground acceleration is saved to a text file and visualized using matplotlib.

The program also includes uncertainty modeling using probability distribution functions (PDFs) and the Monte-Carlo method to account
for variability in soil properties and ground motion intensity.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import time as ti
from SALAR_MATH import BETA_PDF, HISROGRAM_BOXPLOT

# Kanai-Tajimi model differential equations
def kanai_tajimi_model(Y, t, omega_g, zeta_g, white_noise, dt):
    x, v = Y
    index = min(int(t / dt), len(white_noise) - 1)  # Ensure index is within bounds
    dxdt = v
    dvdt = -2 * zeta_g * omega_g * v - omega_g**2 * x + white_noise[index]
    return [dxdt, dvdt]

# Kanai-Tajimi Power Spectral Density function
def kanai_tajimi_psd(omega, omega_g, zeta_g, S_0):
    return S_0 * (1 + (2 * zeta_g * omega / omega_g)**2) / ((1 - (omega / omega_g)**2)**2 + (2 * zeta_g * omega / omega_g)**2)

# Time settings
t_max = 10                  # Total time for simulation (s)
NUM_SIM = 6000              # Total number for simulation (set to 1 for simplicity)
dt = 0.01                   # Time step (s)
# Define the frequency range for the simulation
freq = np.linspace(dt, t_max, 1000)
ZETA_G = BETA_PDF(0.5, 0.7, 1, 1, NUM_SIM) # Damping ratio of the ground
s_0 = BETA_PDF(0.9, 1.1, 1, 1, NUM_SIM)    # Intensity of the ground motion
# Define the Kanai-Tajimi model parameters ground motion  
#  Z and W are variables that play a role in defining the natural frequency of the ground (𝜔𝑔) and, by extension, influence the earthquake simulation         
ZZ = BETA_PDF(1.8, 2.2,1, 1, NUM_SIM)      #Constant Value
WW = BETA_PDF(2.3, 2.7,1, 1, NUM_SIM)      #Variable Value 

AF = 0.005 # Amplitude Factor 

HISROGRAM_BOXPLOT(ZETA_G, HISTO_COLOR='blue', LABEL='Damping ratio of the ground')
HISROGRAM_BOXPLOT(s_0, HISTO_COLOR='purple', LABEL='Velocity')
HISROGRAM_BOXPLOT(ZZ, HISTO_COLOR='green', LABEL='Constant Value')
HISROGRAM_BOXPLOT(WW, HISTO_COLOR='gold', LABEL='Variable Value')
                  
TIME = np.arange(dt, t_max, dt)

# Analysis Durations:
starttime = ti.process_time()
    
for I in range(NUM_SIM):
    # Define the Kanai-Tajimi model parameters
    Z = ZZ[I]
    W = WW[I]
    omega_g = Z * np.pi * W  # Natural frequency of the ground (rad/s)
    omega = Z * np.pi * freq
    zeta_g = ZETA_G[I]
    S_0 = s_0[I]
    # Calculate the PSD
    psd = kanai_tajimi_psd(omega, omega_g, zeta_g, S_0)
    # Find the psd
    max_psd = np.max(np.abs(psd))
    
    # White noise generation
    white_noise = BETA_PDF(-np.sqrt(s_0[I]), np.sqrt(s_0[I]), 0.5, 0.5, len(TIME))

    # Initial conditions
    Y0 = [0, 0]

    # Solve the differential equations
    sol = odeint(kanai_tajimi_model, Y0, TIME, args=(omega_g, zeta_g, white_noise, dt))

    # Extract the ground acceleration time history
    ground_acceleration = sol[:, 1]

    # Add an increasing sinusoidal function to the ground acceleration
    sinusoidal_amplitude = AF * TIME  # Amplitude increases over time
    sinusoidal_frequency = 2 * np.pi  # Frequency of the sinusoidal function
    sinusoidal_component = sinusoidal_amplitude * np.sin(sinusoidal_frequency * TIME)

    # Combine Kanai-Tajimi and sinusoidal components
    combined_acceleration = ground_acceleration + sinusoidal_component

    # Find the maximum absolute combined acceleration
    max_combined_acceleration = np.max(np.abs(combined_acceleration))

    # Print combined acceleration to a text file
    with open(f'Combined_Acceleration_{I+1}.txt', 'w') as file:
        for acc in combined_acceleration:
            file.write(f'{acc:.6f}\n')

    print(f'{I+1} Maximum absolute combined acceleration: {max_combined_acceleration:.6f} m/s²')
    print(f'      Combined acceleration data has been printed to combined_acceleration_{I+1}.txt')

    # Plot the combined acceleration
    plt.figure(figsize=(10, 6))
    plt.plot(TIME, combined_acceleration, label=f'Combined Acceleration: {np.max(np.abs(combined_acceleration)): .3f}', color='black')
    plt.xlabel('Time (s)')
    plt.ylabel('Acceleration (m/s²)')
    plt.title(f'Combined Kanai-Tajimi and Increasing Sinusoidal Ground Acceleration -  Amplitude Factor: {AF:.3f}')
    plt.legend()
    plt.grid(True)
    plt.show()

totaltime = ti.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n')