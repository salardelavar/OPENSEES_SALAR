def PLOT_FORCE_DISP(positive_cycle, negative_cycle):
    import matplotlib.pyplot as plt
    forces_pos = positive_cycle[::2]  # Extract force values
    displacements_pos = positive_cycle[1::2]  # Extract displacement values
    forces_neg = negative_cycle[::2]
    displacements_neg = negative_cycle[1::2]
    # Create the plot
    plt.figure(figsize=(10, 6))
    plt.plot(displacements_pos, forces_pos, 'b-o', label='Positive Loading', linewidth=2)
    plt.plot(displacements_neg, forces_neg, 'r-s', label='Negative Loading', linewidth=2)
    
    # Add labels and title
    plt.xlabel('Displacement (mm)', fontsize=12)
    plt.ylabel('Force (kN)', fontsize=12)
    plt.title('Force-Displacement Hysteresis Loop', fontsize=14)
    
    # Add grid and legend
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend(fontsize=12)
    
    # Show the plot
    plt.tight_layout()
    plt.show() 