def FRAGILITY_ANALYSIS_TWO_PARAMETERS_FUN(damage_states, im1_grid, im2_grid, TITLE,
                          SCATTER='False', CONTOUR='True', SURFACE='True'):
    """
    1. Bivariate fragility modeling: The function evaluates the probability of exceeding 
    predefined damage states given two correlated intensity measures (e.g., spectral acceleration 
    and duration) using a bivariate lognormal cumulative distribution.
    2. Input parameters: Each damage state is defined by medians (η1, η2), logarithmic standard 
    deviations (β1, β2), and the correlation coefficient ρ between the logarithms of the two IMs.
    3. Probability computation: It transforms the IM grids to standardized normal space, then 
    integrates the bivariate standard normal density over the exceedance region 
    z1 ≤ (ln(im1)-ln(η1))/β1 and z2 ≤ (ln(im2)-ln(η2))/β2 using scipy.stats.multivariate_normal.cdf.
    4. Output visualization: The resulting joint probability surface is displayed as both contour 
    lines (showing iso-probability levels 0.1, 0.5, 0.9) and a 3D surface plot, directly mapping 
    IM pairs to damage exceedance likelihood.
    5. Engineering value: This approach captures the interactive effect of two IMs on structural 
    performance, moving beyond univariate fragility curves to a more realistic risk assessment 
    when damage is driven by multiple loading characteristics.
    
    # THIS PYTHON SCRIPT WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.stats import multivariate_normal

    # Create a meshgrid for the two IM dimensions
    IM1, IM2 = np.meshgrid(im1_grid, im2_grid)
    ln_im1 = np.log(IM1)
    ln_im2 = np.log(IM2)

    # Prepare plot
    fig = plt.figure(figsize=(12, 5))

    for state_name, params in damage_states.items():
        # Unpack parameters for this damage state
        med1, beta1 = params['median1'], params['beta1']
        med2, beta2 = params['median2'], params['beta2']
        rho = params['rho']

        # Standardized log values
        z1 = (ln_im1 - np.log(med1)) / beta1
        z2 = (ln_im2 - np.log(med2)) / beta2

        # Compute bivariate CDF at each grid point
        # The covariance matrix is [[1, rho], [rho, 1]]
        cov = np.array([[1, rho], [rho, 1]])
        rv = multivariate_normal(mean=[0, 0], cov=cov)
        # Flatten, compute CDF, reshape back
        Z = rv.cdf(np.dstack((z1, z2)))   # shape (len(im2_grid), len(im1_grid))

        # Plot options
        if CONTOUR == 'True':
            ax = fig.add_subplot(1, 2, 1)
            CS = ax.contour(IM1, IM2, Z, levels=[0.1, 0.2, 0.5, 0.8, 0.9])
            ax.clabel(CS, inline=True, fontsize=8)
            ax.set_xlabel('IM 1')
            ax.set_ylabel('IM 2')
            ax.set_title('Fragility Contours')
            ax.grid(True)

        if SURFACE == 'True':
            ax = fig.add_subplot(1, 2, 2, projection='3d')
            ax.plot_surface(IM1, IM2, Z, cmap='viridis', alpha=0.8)
            ax.set_xlabel('IM 1')
            ax.set_ylabel('IM 2')
            ax.set_zlabel('P(Exceedance)')
            ax.set_title('Fragility 3D')
            #ax.set_title(f'{state_name}')

    #plt.suptitle(f'Bivariate Fragility - {TITLE}', fontsize=14)
    plt.tight_layout()
    plt.show()


#%%----------------------------------------------------
import numpy as np
damage = {
    'Minor Damage Level': {'median1': 0.2, 'beta1': 0.4, 'median2': 5.0, 'beta2': 0.5, 'rho': 0.3},
    'Moderate Damage Level': {'median1': 0.5, 'beta1': 0.4, 'median2': 10.0, 'beta2': 0.5, 'rho': 0.3},
    'Severe Damage Level': {'median1': 0.7, 'beta1': 0.4, 'median2': 5.0, 'beta2': 0.5, 'rho': 0.3},
    'Failure Level': {'median1': 1.0, 'beta1': 0.4, 'median2': 5.0, 'beta2': 0.5, 'rho': 0.3},
}

im1_vals = np.linspace(0.05, 1.5, 50)
im2_vals = np.linspace(1, 25, 50)

FRAGILITY_ANALYSIS_TWO_PARAMETERS_FUN(damage, im1_vals, im2_vals, 'Sa vs Duration', CONTOUR='True', SURFACE='True') 
#%%----------------------------------------------------   
