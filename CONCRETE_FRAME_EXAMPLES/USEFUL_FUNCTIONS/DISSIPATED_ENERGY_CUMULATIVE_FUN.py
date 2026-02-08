def CUMULATIVE_DISSIPATED_ENERGY_FUN(disp, force, TITLE, COLOR):
    import numpy as np
    import matplotlib.pyplot as plt
    disp = np.array(disp)
    force = np.array(force)
    cum_energy = np.zeros(len(disp))
    for i in range(1, len(disp)):
        d_disp = disp[i] - disp[i-1]
        avg_force = 0.5 * (force[i] + force[i-1])
        cum_energy[i] = cum_energy[i-1] + abs(avg_force * d_disp)
    
    plt.figure()
    plt.plot(disp, cum_energy, color=COLOR)
    plt.xlabel('Displacement')
    plt.ylabel(TITLE)
    plt.grid(True)
    plt.show()
    return cum_energy

