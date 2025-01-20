"""
   #######################################################################
   #                            IN THE NAME OF ALLAH                     #
   # Using the Kanai-Tajimi model involves simulating a non-stationary   #
   # stochastic process to represent earthquake ground motion            #
   # Nonlinear Dynamic Analysis on Concrete Confined Section Column      #
   # With Uncertainty Conditions Using Probability Distribution Function #
   # Monte-Carlo Method                                                  #
   #---------------------------------------------------------------------#
   #     THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)      #
   #                 EMAIL: salar.d.ghashghaei@gmail.com                 #
   #######################################################################

Generate ground acceleration using the Kanai-Tajimi model involves simulating a non-stationary stochastic process:

Using the Kanai-Tajimi model involves simulating a non-stationary stochastic process to represent
earthquake ground motion. The Kanai-Tajimi model is a widely used approach to simulate the random ground acceleration
experienced during an earthquake. It is based on the concept of a filtered Gaussian process, where the ground acceleration
is modeled as a white noise process (acceleration at the bedrock) that is filtered through a viscoelastic system representing
the soil layer.

The classical Tajimi-Kanai model consists of a linear oscillator attached to the bedrock, which moves with an acceleration
modeled as a Gaussian white noise process. This model is often used for the analysis of structures subjected to earthquakes.
The parameters of the model, such as the damping ratio ( \zeta_g ) and the frequency ( \omega_g ) of the soil deposit,
are calibrated by observing zero crossings and other statistics of historical earthquakes1.

In more advanced versions, such as the fractional Tajimi-Kanai model, the purely viscous element in the Kelvin-Voigt element
(a spring in parallel with a dashpot) is replaced with a springpot, which exhibits behavior between purely elastic and purely
viscous. This modification allows for a more accurate representation of the viscoelastic behavior of the ground and ensures
that the number of zero crossings of the absolute acceleration at the free field remains finite1.

The non-stationary aspect of the process is typically accounted for by applying an envelope or shape function to the stationary
process, which adjusts the intensity of the ground motion over time, reflecting the way an actual earthquake’s intensity
changes2.
Overall, the Kanai-Tajimi model and its variations provide a mathematical framework to generate artificial earthquake records
that can be used for the design and analysis of structures to withstand seismic events.
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
NUM_SIM = 6000             # Total number for simulation
dt = 0.01                   # Time step (s)
# Define the frequency range for the simulation
freq = np.linspace(dt, t_max, 1000)
# The selection of alpha and beta coefficients in the beta probability distribution is crucial.
#  In uncertainty analysis, careful consideration must also be given to the numerical interval (maximum and minimum) and the alpha and beta coefficients.
ZETA_G = BETA_PDF(0.5, 0.7, 1, 1, NUM_SIM) # Damping ratio of the ground
s_0 = BETA_PDF(0.9, 1.1, 1, 1, NUM_SIM)    # Intensity of the ground motion
# Define the Kanai-Tajimi model parameters ground motion  
#  Z and W are variables that play a role in defining the natural frequency of the ground (𝜔𝑔) and, by extension, influence the earthquake simulation         
ZZ = BETA_PDF(1.8, 2.2,1, 1, NUM_SIM)      #Constant Value
WW = BETA_PDF(2.3, 2.7,1, 1, NUM_SIM)      #Variable Value 

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

    # Find the maximum absolute ground acceleration
    max_ground_acceleration = np.max(np.abs(ground_acceleration))

    # Print ground acceleration to a text file
    with open(f'Ground_Acceleration_{I+1}.txt', 'w') as file:
        new_ground_acceleration = ground_acceleration[1:]
        for acc in new_ground_acceleration:
            file.write(f'{acc:.6f}\n')

    print(f'{I+1} Maximum absolute ground acceleration: {max_ground_acceleration:.6f} m/s²')
    print(f'      Ground acceleration data has been printed to ground_acceleration_{I+1} .txt')
    
totaltime = ti.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n')   
