###########################################################################################################
#                                          IN THE NAME OF ALLAH                                           #
#                               FRAGILITY ANALYSIS BASED ON ACCELERATION                                  #
#---------------------------------------------------------------------------------------------------------#
#                          THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                     #
#                                   EMAIL: salar.d.ghashghaei@gmail.com                                   #
########################################################################################################### 
"""
PAPER:
https://www.jackwbaker.com/Publications/Baker_(2015)_fragility_fitting_EQS_final.pdf

1. IM_levels: Peak Ground Acceleration (PGA) as the Intensity Measure
The array IM_levels lists the Intensity Measure (IM) levels used in the analysis, with units in g (gravitational acceleration).
 In this case, the IM is Peak Ground Acceleration (PGA), which measures the maximum acceleration the ground experiences during an earthquake.
 Think of PGA as the “shake factor”—it’s a go-to metric in earthquake engineering to gauge how intense the ground motion is.

Range and Increments: The values span from 0.1g to 0.8g, stepping up by 0.1g each time. This range covers a spectrum of seismic intensities:
0.1g: Minor shaking, typical of small earthquakes or distant events—structures should handle this without breaking a sweat.
0.8g: Severe shaking, characteristic of strong, near-fault earthquakes that can push structures to their limits.
Why This Matters: These levels let us test the structure across a gradient of earthquake intensities, from “barely noticeable” to “potential disaster.
” It’s a systematic way to see how performance degrades as the ground motion amps up
---------------------------------------------
2. n_analyses: Number of Analyses per PGA Level
The array n_analyses tells us how many structural analyses were run at each PGA level. For every value in IM_levels, we’ve got 20 simulations.

What’s Happening Here: Each analysis likely involves running a structural model (say, a finite element model of a building) through a unique ground motion
 record scaled to the corresponding PGA. Imagine 20 different earthquake “movies,” each with the same peak acceleration but starring different frequency content, durations,
 and wiggles.
Why 20?: Earthquakes are messy—two ground motions with the same PGA can stress a structure differently due to their spectral shapes or energy distribution.
 By using 20 analyses per level, we’re capturing that variability, ensuring our results aren’t skewed by a fluke record. This is standard practice in multiple stripe analysis (MSA), a technique we’ll revisit later.
---------------------------------------------
3. n_failures: Number of Failures per PGA Level
The array n_failures shows how many of those 20 analyses resulted in the structure exceeding a damage state at each PGA level.
 “Failure” here could mean collapse, yielding, excessive drift, or another performance limit defined by the engineer—whatever threshold says, “This structure’s in trouble.”

Breaking It Down:
PGA = 0.1g: 0 failures out of 20. The structure sails through low-intensity shaking unscathed.
PGA = 0.2g: 1 failure. Damage starts creeping in, but it’s rare.
PGA = 0.4g: 6 failures. We’re hitting a tipping point—30% of the simulations show problems.
PGA = 0.8g: 19 failures. At this brutal intensity, 95% of the analyses flag the structure as failing.
The Trend: As PGA climbs, failures skyrocket. This makes sense—stronger shaking pushes the structure past its capacity more often. It’s a classic dose-response curve
 for seismic vulnerability.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.optimize import minimize
#------------------------------------------------------------------------
# Simulated structural analysis data
IM_levels = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])  # PGA in g
n_analyses = np.array([20, 20, 20, 20, 20, 20, 20, 20])  # Number of analyses per level
n_failures = np.array([0, 1, 3, 6, 10, 14, 17, 19])  # Number of failures per level
#------------------------------------------------------------------------
# Define the lognormal fragility function
def fragility_function(IM, theta, beta):
    """Compute the probability of exceeding a damage state using lognormal CDF."""
    return norm.cdf((np.log(IM) - theta) / beta)
#------------------------------------------------------------------------
# Define the negative log-likelihood function for MLE
def neg_log_likelihood(params, IM, n_fail, n_tot):
    """Calculate the negative log-likelihood to be minimized."""
    theta, beta = params
    if beta <= 0:  # Ensure beta is positive
        return np.inf
    p = fragility_function(IM, theta, beta)
    p = np.clip(p, 1e-10, 1 - 1e-10)  # Prevent log(0) issues
    likelihood = n_fail * np.log(p) + (n_tot - n_fail) * np.log(1 - p)
    return -np.sum(likelihood)
#------------------------------------------------------------------------
# Initial parameter guesses: theta = ln(median), beta = log standard deviation
initial_params = [np.log(0.5), 0.4]
#------------------------------------------------------------------------
# Fit the fragility curve by minimizing the negative log-likelihood
result = minimize(
    neg_log_likelihood,
    initial_params,
    args=(IM_levels, n_failures, n_analyses),
    method='Nelder-Mead',
    bounds=[(-5, 5), (0.01, 2)]  # Bounds for theta and beta
)
#------------------------------------------------------------------------
# Extract fitted parameters
theta_fit, beta_fit = result.x
print(f"Fitted parameters: theta = {theta_fit:.3f}, beta = {beta_fit:.3f}")
print(f"Median PGA = {np.exp(theta_fit):.3f} g")
#------------------------------------------------------------------------
# Generate data for the fragility curve plot
IM_range = np.linspace(0.01, 1.0, 100)
fragility_curve = fragility_function(IM_range, theta_fit, beta_fit)
#------------------------------------------------------------------------
# Plot the results
plt.figure(figsize=(8, 6))
plt.plot(IM_range, fragility_curve, label='Fitted Fragility Curve', color='blue')
plt.scatter(IM_levels, n_failures / n_analyses, color='red', label='Observed Data', zorder=5)
plt.xlabel('Peak Ground Acceleration (PGA, g)')
plt.ylabel('Probability of Exceeding Damage State')
plt.title('Fragility Curve Fitting Based on Baker (2015)')
plt.legend()
plt.grid(True)
plt.show()
#------------------------------------------------------------------------