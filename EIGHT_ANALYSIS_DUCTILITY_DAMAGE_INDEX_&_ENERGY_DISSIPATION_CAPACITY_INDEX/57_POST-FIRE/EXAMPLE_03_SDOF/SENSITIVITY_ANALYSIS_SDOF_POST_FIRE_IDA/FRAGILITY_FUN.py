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


def INTERPRET_FRAGILITY(damage_states, im_values):
    # THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
    import numpy as np
    from scipy.stats import norm
    print("="*60)
    print("FRAGILITY ANALYSIS INTERPRETATION")
    print("="*60)
    
    # 1. Basic Statistics
    print(f"\n[1] INPUT SUMMARY:")
    print(f"    - Number of IM samples: {len(im_values)}")
    print(f"    - IM range: {min(im_values):.2f}% to {max(im_values):.2f}%")
    print(f"    - Average IM: {np.mean(im_values):.2f}%")
    
    # 2. Damage State Interpretation
    print(f"\n[2] DAMAGE STATE PARAMETERS:")
    for state, (median, beta) in damage_states.items():
        print(f"    - {state}:")
        print(f"        Median (50% exceedance) at DI = {median}%")
        print(f"        Dispersion (uncertainty) = {beta}%")
        if beta < 30:
            print(f"        -> Low uncertainty (confident prediction)")
        elif beta < 50:
            print(f"        -> Moderate uncertainty")
        else:
            print(f"        -> High uncertainty (significant variability)")
    
    # 3. Probability Analysis at Mean IM
    mean_im = np.mean(im_values)
    print(f"\n[3] PROBABILITIES AT MEAN IM ({mean_im:.2f}%):")
    for state, (median, beta) in damage_states.items():
        ln_im = np.log(mean_im) if mean_im > 0 else 0
        ln_median = np.log(median)
        prob = norm.cdf((ln_im - ln_median) / beta) if mean_im > 0 else 0
        print(f"    - {state}: {prob*100:.1f}% probability of exceedance")
    
    # 4. Critical Damage Assessment
    print(f"\n[4] CRITICAL FINDINGS:")
    # Find probability of failure at max IM
    max_im = max(im_values)
    if max_im > 0:
        failure_median = damage_states['Failure Level'][0]
        failure_beta = damage_states['Failure Level'][1]
        p_failure = norm.cdf((np.log(max_im) - np.log(failure_median)) / failure_beta)
        print(f"    - At peak DI ({max_im:.2f}%):")
        print(f"        -> Failure probability = {p_failure*100:.1f}%")
        if p_failure > 0.5:
            print(f"        [!]  HIGH RISK: Structure likely to FAIL")
        elif p_failure > 0.2:
            print(f"        [!]  MODERATE RISK: Failure possible")
        else:
            print(f"        [V] LOW RISK: Structure likely safe")
    
    # 5. Damage Hierarchy
    print(f"\n[5] DAMAGE HIERARCHY (Easiest -> Hardest to occur):")
    sorted_states = sorted(damage_states.items(), key=lambda x: x[1][0])
    for i, (state, _) in enumerate(sorted_states, 1):
        print(f"    {i}. {state}")
    
    print("\n" + "="*60)
    print("INTERPRETATION COMPLETE")
    print("="*60)
    
    return {
        'mean_im': np.mean(im_values),
        'max_im': max(im_values),
        'probabilities': {state: norm.cdf((np.log(np.mean(im_values)) - np.log(median))/beta) 
                         if np.mean(im_values) > 0 else 0 
                         for state, (median, beta) in damage_states.items()}
    }    