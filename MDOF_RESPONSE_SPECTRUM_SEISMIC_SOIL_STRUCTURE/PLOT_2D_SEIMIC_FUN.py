def PLOT_2D_SEIMIC(NUM_SEISMIC, COUNT, XLABEL, YLABEL, TITLE, SEMILOGY=False):
    import matplotlib.pyplot as plt
    plt.figure(COUNT, figsize=(12, 10))

    # Plot all simulations
    for j in range(NUM_SEISMIC):
        #plt.plot(T, DATA[COUNT][j], alpha=0.4, label=f'Sim {j+1} | Max: {np.max(np.abs(DATA[COUNT][j])):.4e}')
        plt.plot(T, DATA[COUNT][j], alpha=0.4)
    
    # Convert to NumPy array for statistical calculations
    arr = np.array(DATA[COUNT])    # Shape: (NUM_SIM, len(T))

    # Compute statistical metrics
    mean_curve   = np.mean(arr, axis=0)
    median_curve = np.median(arr, axis=0)
    std_curve    = np.std(arr, axis=0)
    std_curveDOWN = mean_curve - np.std(arr, axis=0)
    std_curveUP = mean_curve + np.std(arr, axis=0)
    q1_curve = np.percentile(arr, 25, axis=0)
    q3_curve = np.percentile(arr, 75, axis=0)

    # Plot statistical curves
    plt.plot(T, mean_curve, color='navy', linestyle='-.', linewidth=3, label='Mean')
    
    #plt.plot(T, std_curveDOWN, color='steelblue', linestyle='-', linewidth=2, label='Mean - Std')
    #plt.plot(T, std_curveUP,   color='deepskyblue', linestyle='-', linewidth=2, label='Mean + Std')
    
    plt.plot(T, median_curve, color='red', linestyle='--', linewidth=2, label='Median')
    
    plt.plot(T, q1_curve, color='green', linestyle='-.', linewidth=2, label='Q1 (25%)')
    plt.plot(T, q3_curve, color='orange', linestyle='-.', linewidth=2, label='Q3 (75%)')

    # Standard deviation band around the mean
    plt.fill_between(T, mean_curve - std_curve, mean_curve + std_curve,
                     color='gray', alpha=0.25, label='Mean ± Std')

    # Labels and title
    plt.xlabel(XLABEL)
    plt.ylabel(YLABEL)
    plt.title(TITLE)
    
    if SEMILOGY == True:
        plt.semilogy()
        
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()