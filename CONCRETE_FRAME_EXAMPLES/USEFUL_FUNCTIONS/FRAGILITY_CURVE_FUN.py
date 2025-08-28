def FRAGILITY_CURVE(im_values, damage_states, X_LABEL_TEXT, PLOT):
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.stats import norm
        
    # Generate intensity measure (IM) values from 0.0 to 1.0
    #im_values = max_DI # Structural Ductility Damage Index

    # Create plot
    plt.figure(figsize=(10, 6))
    # Calculate and plot fragility curves for each damage state
    for damage_state, (median, beta) in damage_states.items():
        # Calculate log-normal probabilities
        ln_im = np.log(im_values)
        ln_median = np.log(median)
        probabilities = norm.cdf((ln_im - ln_median) / beta)
        if PLOT == False:
            plt.scatter(im_values, probabilities, marker='o', label=f'{damage_state} (η={median}, β={beta}')
        if PLOT == True:
            plt.plot(im_values, probabilities, lw=2, label=f'{damage_state} (η={median}, β={beta})')
    
    # --------------
    # Visualization
    # --------------
    # Format plot
    plt.xlabel(f'{X_LABEL_TEXT} [IM]', fontsize=12)
    plt.ylabel('Probability of Exceedance', fontsize=12)
    plt.title('Fragility Curves', fontsize=14)
    plt.legend(loc='lower right', fontsize=10)
    plt.grid(True)
    plt.semilogy()
    plt.ylim(0, 1.0)
    plt.tight_layout()
    plt.show()

"""
####  FRAGILITY ANALYSIS BASED ON ACCELERATION : 
damage_states = {
'DS1_Slight': (0.15, 0.4),    # Median PGA=0.15g, β=0.4
'DS2_Moderate': (0.30, 0.5),
'DS3_Extensive': (0.60, 0.6),
'DS4_Complete': (1.00, 0.7)
}

####  FRAGILITY ANALYSIS BASED ON STRUCTURAL DUCTILITY DAMAGE INDEX:

# Define damage state parameters: {Damage State: (median_IM, beta)}
damage_states = {
    'Minor Damage Level': (0.2, 0.4),# Median DI=0.2, β=0.4
    'Moderate Damage Level': (0.4, 0.4),
    'Severe Damage Level': (0.6, 0.5),
    'Failure Level': (1.0, 0.5)
}
FRAGILITY_CURVE(im_values, damage_states, 'Ductility Damage Index', PLOT=false)
"""
