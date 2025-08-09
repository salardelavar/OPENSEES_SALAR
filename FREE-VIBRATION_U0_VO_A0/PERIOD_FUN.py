def PERIOD_FUN(disp, dt):
    from scipy.fft import fft, fftfreq
    import numpy as np
    # --- FFT to find dominant frequency ---
    N = len(disp)
    yf = fft(disp)
    xf = fftfreq(N, dt)[:N//2]
    
    # Find peak frequency in positive frequencies
    idx_peak = np.argmax(np.abs(yf[:N//2]))
    freq_peak = xf[idx_peak]
    period_nonlinear = 1.0 / freq_peak
    
    print(f"\nEstimated from FFT: {period_nonlinear:.4f} [s]\n\n")
    return period_nonlinear