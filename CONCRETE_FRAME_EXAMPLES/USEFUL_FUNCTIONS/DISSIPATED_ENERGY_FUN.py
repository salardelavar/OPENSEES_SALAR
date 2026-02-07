def DISSIPATED_ENERGY_FUN(force, disp):
    import numpy as np
    disp = np.array(disp)
    force = np.array(force)
    energy_cumulative = np.cumsum(np.abs(np.diff(disp) * (force[:-1] + force[1:]) / 2))
    print("Dissipated Energy =", energy_cumulative)
    return energy_cumulative