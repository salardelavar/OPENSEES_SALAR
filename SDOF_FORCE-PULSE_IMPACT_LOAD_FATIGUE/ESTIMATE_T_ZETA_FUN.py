# Function to estimate T and zeta
def ESTIMATE_T_ZETA(disp, dt):
    from scipy.optimize import fsolve
    import numpy as np
    n = len(disp)
    freqs = np.fft.rfftfreq(n, d=dt)
    spectrum = np.abs(np.fft.rfft(disp))
    peak_idx = np.argmax(spectrum[1:]) + 1
    fn = freqs[peak_idx]
    T = 1.0 / fn

    peaks, _ = find_peaks(np.abs(disp))
    if len(peaks) > 1:
        x1 = np.abs(disp[peaks[0]])
        x2 = np.abs(disp[peaks[1]])
        delta = np.log(x1 / x2)
    else:
        delta = np.nan
        
    # Define equation for root finding
    def EQUATION(x, delta):
       if np.any(x == 0):  # avoid divide by zero
           return np.inf
       return x**2 - 1 + ((2 * np.pi * x) / np.mean(delta)) ** 2  
     
    # Solve for damping ratio
    x0 = 0.0
    solution = fsolve(EQUATION, x0, args=(delta))   
    zeta = solution[0]
    
    return T, zeta