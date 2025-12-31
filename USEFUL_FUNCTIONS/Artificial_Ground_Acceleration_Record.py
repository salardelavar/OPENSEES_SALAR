###########################################################################################################
#                     >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                    #
#                                 GENERATE ARTIFICIAL GROUND ACCELERATION                                 #
#---------------------------------------------------------------------------------------------------------#
#                          THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                     #
#                                   EMAIL: salar.d.ghashghaei@gmail.com                                   #
###########################################################################################################
import numpy as np
import matplotlib.pyplot as plt

# Generate Artificial Ground Acceleration Record
def GENERATE_ARTIFICIAL_ACCEL(duration, dt, max_accel=0.3*9.81):
    """
    Generate artificial acceleration record
    """
    import numpy as np
    npts = int(duration/dt)
    t = np.linspace(0, duration, npts)
    
    # Combine multiple sine waves with different frequencies
    acc = np.zeros(npts)
    frequencies = [0.5, 1.0, 2.0, 3.0, 5.0]  # Hz
    
    for freq in frequencies:
        amplitude = max_accel / len(frequencies)
        phase = np.random.random() * 2 * np.pi
        acc += amplitude * np.sin(2 * np.pi * freq * t + phase)
    
    # Add random noise
    noise = 0.1 * max_accel * np.random.randn(npts)
    acc += noise
    
    # Windowing
    window = np.ones(npts)
    ramp_up = int(0.1 * npts)
    ramp_down = int(0.9 * npts)
    
    window[:ramp_up] = np.linspace(0, 1, ramp_up)
    window[ramp_down:] = np.linspace(1, 0, npts - ramp_down)
    
    acc *= window
    
    return t, acc
print("Generating artificial acceleration record...")
t_acc, acc_record = GENERATE_ARTIFICIAL_ACCEL(duration=20.0, dt=0.01, max_accel=0.3*9.81)

plt.figure(figsize=(10, 4))
plt.plot(t_acc, acc_record, color='black', linewidth=1)
plt.xlabel('Time [s]')
plt.ylabel('Earthquake Ground Acceleration [m/s²]')
plt.title('Artificial Acceleration Record for Response Spectrum Analysis')
plt.grid(True)
plt.tight_layout()
plt.show()
