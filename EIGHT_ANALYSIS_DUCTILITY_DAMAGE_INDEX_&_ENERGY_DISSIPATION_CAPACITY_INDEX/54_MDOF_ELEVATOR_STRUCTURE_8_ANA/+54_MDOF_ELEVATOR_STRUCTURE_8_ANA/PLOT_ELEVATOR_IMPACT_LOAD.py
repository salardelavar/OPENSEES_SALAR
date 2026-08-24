def PLOT_ELEVATOR_IMPACT_LOAD(
    W_kg=2300,          # Total weight (cabin + capacity + counterweight) in kg
    v_mps=2.0,          # [m/s] Impact velocity
    S_m=0.45,           # [m] Buffer stroke
    duration=0.5,       # [s] Impact duration
    dt= 0.01,
    g=9.81
):
    import numpy as np
    import matplotlib.pyplot as plt
    
    W = W_kg * g                     # Weight in Newtons
    # Maximum impact force using buffer formula
    R_max = W * (1 + (v_mps**2) / (2 * g * S_m))
    R_max_N = R_max
    N_STEPS  = int(duration / dt) 
    # Time array for simulation (Half-Sine pulse - common for impact loads)
    t = np.linspace(0, duration, N_STEPS)
    force = R_max * np.sin(np.pi * t / duration) 
    #force = R_max * np.sin(np.pi * t / duration) + np.cos(np.pi * t / duration)
    #force = R_max * np.sin(np.pi * t / duration) * np.exp(0.05 * t/ duration)
    force[t > duration] = 0
    
    # Display key information
    print(f"Total Weight: {W_kg} kg")
    print(f"Maximum Impact Force: {R_max_N:.1f} N")
    print(f"Impact Velocity: {v_mps} m/s | Buffer Stroke: {S_m} m")
    
    # Plot the impact load
    plt.figure(figsize=(10, 6))
    plt.plot(t, force, 'b-', linewidth=2.5, label='Elevator Impact Load')
    plt.fill_between(t, force, alpha=0.3, color='blue')
    
    plt.title('Dynamic Elevator Impact Load (Buffer Impact)', fontsize=14, fontweight='bold')
    plt.xlabel('Time (seconds)', fontsize=12)
    plt.ylabel('Force (N)', fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.axhline(y=R_max_N, color='r', linestyle='--', alpha=0.7, 
                label=f'Max Force = {R_max_N:.1f} N')
    plt.legend(fontsize=11)
    
    # Add text box with parameters
    plt.text(0.02, 0.88, f'Total Weight: {W_kg} kg\n'
                        f'Approximate Impact Factor: {R_max/W:.2f}', 
             transform=plt.gca().transAxes, 
             bbox=dict(facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    plt.show()
    
    return R_max_N


# Run the visualization
if __name__ == "__main__":
    print("=== Elevator Impact Load Visualization ===\n")
    max_force = PLOT_ELEVATOR_IMPACT_LOAD(
        W_kg=2300,      # Change as needed
        v_mps=1.8,      # Impact velocity
        S_m=0.4,        # Buffer stroke
        duration=0.5,   # Impact duration
        dt=0.01
    )