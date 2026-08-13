def FRAGILITY_CURVE_FUN_02(IMs, DAMAGE_STATES, XLABEL, SCATTER, SEMILOGY, A, B, BETA):
    """
    Fragility Curve Generator for Seismic Risk Analysis
    This script calculates and visualizes fragility curves for structural
    components or systems based on nonlinear response data from a suite of
    ground motions.
    #
    Fragility curves represent the conditional probability that a structure
    will exceed a predefined damage state for a given ground motion intensity
    measure (IM), typically spectral acceleration (Sa). These curves are
    essential tools for probabilistic seismic risk assessment and
    performance-based earthquake engineering.
    #
    The script uses a lognormal cumulative distribution function (CDF) to fit
    observed probabilities of exceedance derived from drift thresholds and
    structural response data. Drift ratios are compared against damage state
    thresholds to calculate these probabilities at various IM bins.

    Parameters
    ----------
    IMs : TYPE
        DESCRIPTION.
    DAMAGE_STATES : TYPE
        DESCRIPTION.
    XLABEL : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.
        
    -------
    A is the median EDP (e.g., drift ratio, acceleration) when the intensity measure 
    IM=1 (in the units you choose for IM, e.g., g for spectral acceleration).
    It scales the entire demand curve vertically.
    -------
    B = 1 → linear relationship (elastic behaviour).
    B > 1 → demand grows faster than linearly (common for inelastic structures).
    B < 1 → demand grows slower (e.g., base‑shear‑controlled responses).
    A and B are often derived from nonlinear dynamic analyses (e.g., regression of drifts vs. Sa).
    BETA typically ranges from 0.3 to 0.6 for structural components, reflecting the inherent variability in response.
    -------
    THIS PYTHON SCRIPT WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI) 
    """
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from scipy.stats import lognorm
    from scipy.optimize import curve_fit

    IMs = np.abs(IMs)
    
    # EDPs: power-law with scatter
    edp_median = A * (IMs ** B)
    n_samples = len(IMs)
    epsilon = np.random.lognormal(mean=0, sigma=BETA, size=n_samples)
    EDPs = edp_median * epsilon

    fragility_results = {}
    unique_IMs = np.linspace(min(IMs), max(IMs), 20)
    tolerance = 0.05

    for ds, threshold in DAMAGE_STATES.items():
        prob_exceed = []
        for im in unique_IMs:
            edp_subset = EDPs[np.abs(IMs - im) < tolerance]
            if len(edp_subset) > 0:
                prob = np.mean(edp_subset >= threshold)
            else:
                prob = 0
            prob_exceed.append(prob)

        def lognorm_cdf(x, mu, sigma):
            return lognorm.cdf(x, s=sigma, scale=np.exp(mu))

        prob_array = np.array(prob_exceed)
        if np.all(prob_array == 0) or np.all(prob_array == 1):
            print(f" Skipping '{ds}' — flat probability distribution (all 0s or all 1s)")
            fragility_results[ds] = (np.nan, np.nan)
            continue

        try:
            popt, _ = curve_fit(lognorm_cdf, unique_IMs, prob_exceed, p0=[np.log(np.mean(unique_IMs)), 0.5])
            mu, sigma = popt
            fragility_results[ds] = (mu, sigma)
            print(f" Fitted '{ds}': mu = {mu:.3f}, sigma = {sigma:.3f}")
        except Exception as e:
            print(f" Failed to fit for '{ds}': {e}")
            fragility_results[ds] = (np.nan, np.nan)

    plt.figure(figsize=(10, 6))
    x_vals = np.linspace(min(IMs), max(IMs), 300)

    for ds, (mu, sigma) in fragility_results.items():
        if np.isnan(mu):
            continue
        y_vals = lognorm.cdf(x_vals, s=sigma, scale=np.exp(mu))
        if SCATTER == True:
            # Plot every 5th point to keep scatter readable
            plt.scatter(x_vals[::5], y_vals[::5], label=ds, marker='o', s=40)
        else:
            plt.plot(x_vals, y_vals, label=ds, linewidth=2)

    plt.ylim(0, 1)
    plt.xlabel(XLABEL, fontsize=12)
    plt.ylabel("Probability of Exceedance", fontsize=12)
    plt.title(f"Seismic Fragility Curves by Damage State from {XLABEL}", fontsize=14)
    #plt.title("Seismic Fragility Curves by Damage State (Synthetic Data)", fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend(title="Damage State")
    plt.tight_layout()
    if SEMILOGY == True:
        plt.semilogy()
    plt.savefig(f"FRAGILITY_CURVE_FUN_{XLABEL}.png", dpi=300)
    plt.show()
    print(f" Fragility plot saved as 'FRAGILITY_CURVE_FUN_{XLABEL}.png'")
#%% ------------------------------------------------------------

# Ensure thresholds are in the same units as EDPs; if not, scale thresholds.
# Here we use the thresholds as given (10,30,70,100) – they should be
# comparable to the EDP range.
# We'll scale the thresholds to match the typical EDP values.
# For demonstration, we keep them as is, but you can adjust.
# In practice, you should set thresholds based on physical limits.
# Let's just use them directly; they will likely be exceeded for high IM.
#%% ------------------------------------------------------------
# 1. Generate synthetic random data
import numpy as np

n_samples = 200
# IMs: lognormally distributed (mean ~0.5g, std ~0.3g)
IMs = np.random.lognormal(mean=-0.7, sigma=0.6, size=n_samples)
IMs = np.clip(IMs, 0.05, 2.0)  # keep within reasonable range
#%% ------------------------------------------------------------

# 2. Call the fragility function

DAMAGE_STATES = {
    'Minor Damage Level': 0.10, # Median Value = 20%
    'Moderate Damage Level': 0.30,
    'Severe Damage Level': 0.70,
    'Failure Level': 1.00
}


A = 0.5       # Scale factor (median EDP at IM = 1)
B = 1.2       # Power‑law exponent (elastic or inelastic response trend)
BETA = 0.4    # Dispersion (logarithmic standard deviation of EDP)

IM, XLABEL, SCATTER, SEMILOGY = IMs, 'Spectral Acceleration Sa(T1) [g]', False, False
FRAGILITY_CURVE_FUN_02(IM, DAMAGE_STATES, XLABEL, SCATTER, SEMILOGY, A, B, BETA)
