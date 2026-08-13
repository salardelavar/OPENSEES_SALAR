def FRAGILITY_CURVE_FUN(IMs, EDPs):
    """
    Compute and plot fragility curves for four damage states.
    IMs : array-like, intensity measures
    EDPs : array-like, engineering demand parameters (same length as IMs)
    
    THIS PYTHON SCRIPT WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
    """
    # LINK INFO: https://www.sciencedirect.com/topics/engineering/fragility-curves#definition
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from scipy.stats import lognorm
    from scipy.optimize import curve_fit
    
    # Damage state thresholds (example values, adjust as needed)
    damage_states = {
        'Slight': 10.0,
        'Moderate': 30.0,
        'Extensive': 70.0,
        'Complete': 100.0
    }
    fragility_results = {}
    unique_IMs = np.linspace(min(IMs), max(IMs), 20)
    tolerance = 0.05 * (max(IMs) - min(IMs))  # adaptive tolerance

    for ds, threshold in damage_states.items():
        prob_exceed = []
        for im in unique_IMs:
            # Select EDPs whose IM is within tolerance
            edp_subset = EDPs[np.abs(IMs - im) < tolerance]
            if len(edp_subset) > 0:
                prob = np.mean(edp_subset >= threshold)
            else:
                prob = 0.0
            prob_exceed.append(prob)

        prob_array = np.array(prob_exceed)

        # Skip if probability is flat (all 0 or all 1) – no meaningful fit
        if np.all(prob_array == 0) or np.all(prob_array == 1):
            print(f" Skipping '{ds}' — flat probability distribution (all 0s or all 1s)")
            fragility_results[ds] = (np.nan, np.nan)
            continue

        # Define lognormal CDF for fitting
        def lognorm_cdf(x, mu, sigma):
            return lognorm.cdf(x, s=sigma, scale=np.exp(mu))

        # Initial guess: mu from median IM, sigma=0.5
        try:
            popt, _ = curve_fit(lognorm_cdf, unique_IMs, prob_exceed,
                                p0=[np.log(np.mean(unique_IMs)), 0.5],
                                bounds=([-10, 0.01], [10, 5]))  # reasonable bounds
            mu, sigma = popt
            fragility_results[ds] = (mu, sigma)
            print(f" Fitted '{ds}': mu = {mu:.3f}, sigma = {sigma:.3f}")
        except Exception as e:
            print(f" Failed to fit for '{ds}': {e}")
            fragility_results[ds] = (np.nan, np.nan)

    # Plotting
    plt.figure(figsize=(10, 6))
    x_vals = np.linspace(min(IMs), max(IMs), 300)
    for ds, (mu, sigma) in fragility_results.items():
        if np.isnan(mu):
            continue
        plt.plot(x_vals, lognorm.cdf(x_vals, s=sigma, scale=np.exp(mu)),
                 label=ds, linewidth=2)

    plt.ylim(0, 1)
    plt.xlabel("Spectral Acceleration Sa(T1) [g]", fontsize=12)
    plt.ylabel("Probability of Exceedance", fontsize=12)
    plt.title("Seismic Fragility Curves by Damage State (Synthetic Data)", fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend(title="Damage State")
    plt.tight_layout()
    plt.savefig("FRAGILITY_CURVE_FUN.png", dpi=300)
    plt.show()
    print("Fragility plot saved as 'FRAGILITY_CURVE_FUN.png'")

#%% ------------------------------------------------------------
# 1. Generate synthetic random data
import numpy as np

n_samples = 200
# IMs: lognormally distributed (mean ~0.5g, std ~0.3g)
IMs = np.random.lognormal(mean=-0.7, sigma=0.6, size=n_samples)
IMs = np.clip(IMs, 0.05, 2.0)  # keep within reasonable range

# EDPs: power-law with scatter
a = 50.0      # scale factor
b = 1.2       # exponent
beta = 0.4    # lognormal dispersion
edp_median = a * (IMs ** b)
epsilon = np.random.lognormal(mean=0, sigma=beta, size=n_samples)
EDPs = edp_median * epsilon

# Ensure thresholds are in the same units as EDPs; if not, scale thresholds.
# Here we use the thresholds as given (10,30,70,100) – they should be
# comparable to the EDP range.
# We'll scale the thresholds to match the typical EDP values.
# For demonstration, we keep them as is, but you can adjust.
# In practice, you should set thresholds based on physical limits.
# Let's just use them directly; they will likely be exceeded for high IM.

#%% ------------------------------------------------------------
# 2. Call the fragility function

FRAGILITY_CURVE_FUN(IMs, EDPs)
