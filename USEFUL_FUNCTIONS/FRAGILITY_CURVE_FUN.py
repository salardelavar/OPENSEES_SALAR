def FRAGILITY_CURVE(im_values, damage_states, X_LABEL, SEMILOGY=True, PLOT_KIND=True):
    # THIS FUNCTION WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
    from scipy.stats import norm
    import numpy as np
    import matplotlib.pyplot as plt

    # Fragility curves
    plt.figure(1, figsize=(12, 8))
    # Calculate and plot fragility curves for each damage state
    for damage_state, (median, beta) in damage_states.items():
        # Calculate log-normal probabilities
        ln_im = np.log(im_values)
        ln_median = np.log(median)
        probabilities = norm.cdf((ln_im - ln_median) / beta)
        if PLOT_KIND == False:
            plt.scatter(im_values, probabilities, marker='o', label=f'{damage_state} (η={median}, β={beta}')
        if PLOT_KIND == True:
            plt.plot(im_values, probabilities, lw=2, label=f'{damage_state} (η={median}, β={beta})')
    plt.xlabel(X_LABEL)
    plt.ylabel('Probability of Exceedance')
    plt.title('Fragility Curves')
    plt.legend()
    if PLOT_KIND == True:
        plt.semilogy()
    plt.ylim(0, 1.0)
    plt.grid(True)
    plt.tight_layout()
    plt.show()