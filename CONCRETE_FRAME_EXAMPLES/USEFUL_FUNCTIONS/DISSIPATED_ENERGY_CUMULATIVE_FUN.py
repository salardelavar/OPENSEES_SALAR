def CUMULATIVE_DISSIPATED_ENERGY_FUN(disp, force):
    import numpy as np
    import matplotlib.pyplot as plt
    cum_energy = np.zeros(len(disp))
    for i in range(1, len(disp)):
        d_disp = disp[i] - disp[i-1]
        avg_force = 0.5 * (force[i] + force[i-1])
        cum_energy[i] = cum_energy[i-1] + abs(avg_force * d_disp)
    
    plt.figure()
    plt.plot(disp, cum_energy, 'r-')
    plt.xlabel('Displacement')
    plt.ylabel('Cumulative Dissipated Energy')
    plt.grid(True)
    plt.show()