####  FRAGILITY ANALYSIS
def FRAGILITY_ANALYSIS(damage_states, im_values, TITLE, SCATTER='False', SEMI_LOG='False'):
    # THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI) 
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.stats import norm
    # --------------
    # Visualization
    # --------------
    # Create plot
    plt.figure(figsize=(10, 6))
    # Calculate and plot fragility curves for each damage state
    for damage_state, (median, beta) in damage_states.items():
        # Calculate log-normal probabilities
        ln_im = np.log(im_values)
        ln_median = np.log(median)
        probabilities = norm.cdf((ln_im - ln_median) / beta)
        if SCATTER == 'True':
            plt.scatter(im_values, probabilities, marker='o', label=f'{damage_state} (η={median}, β={beta})')
        else:    
            plt.plot(im_values, probabilities, lw=2, label=f'{damage_state} (η={median}, β={beta})')

    # Format plot
    plt.xlabel(TITLE, fontsize=12)
    plt.ylabel('Probability of Exceedance', fontsize=12)
    plt.title(f'Fragility Curves - {TITLE}', fontsize=14)
    plt.legend(loc='lower right', fontsize=10)
    plt.grid(True)
    if SEMI_LOG == 'True':
        plt.semilogy()
    plt.ylim(0, 1.0)
    plt.tight_layout()
    plt.show()    