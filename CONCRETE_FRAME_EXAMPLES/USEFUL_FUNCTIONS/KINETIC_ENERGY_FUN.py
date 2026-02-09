def KINETIC_ENERGY_FUN(mass, velocity):
    import numpy as np
    # kinetic energy
    kinetic_energy = 0.5 * mass * velocity**2
    Ek_min = np.min(kinetic_energy)
    Ek_mean = np.mean(kinetic_energy)
    Ek_median = np.median(kinetic_energy)
    Ek_max = np.max(kinetic_energy)
    return Ek_min, Ek_mean, Ek_median, Ek_max