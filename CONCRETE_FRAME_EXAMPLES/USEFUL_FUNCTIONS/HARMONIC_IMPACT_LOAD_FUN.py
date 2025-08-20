def HARMONIC_IMPACT_LOAD(force_amplitude, impact_time, DT, omega_n, PLOT):
    import numpy as np
    import matplotlib.pyplot as plt
    # Harmonic Load Durations
    num_steps = int(impact_time / DT)
    load_time = np.linspace(0, 0.1 * impact_time, num_steps) # Force time 10 is percent of total time
    
    target_frequency = 0.65 * omega_n  # Target excitation frequency
    
    harmonic_load = force_amplitude * np.sin(target_frequency * load_time)
    if PLOT == True:
        # Plot harmonic loading
        plt.figure(figsize=(10, 6))
        plt.plot(load_time, harmonic_load, label='Harmonic Loading', color='purple',linewidth=5)
        plt.title('Harmonic Loading Over Time')
        plt.xlabel('Time (s)')
        plt.ylabel('Force (N)')
        plt.grid(True)
        plt.legend()
        plt.show()
    return  load_time, harmonic_load  
