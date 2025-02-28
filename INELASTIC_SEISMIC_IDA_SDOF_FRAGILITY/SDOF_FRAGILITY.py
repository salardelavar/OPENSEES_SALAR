"""
Nonlinear Dynamic Analysis and Fragility Assessment of a Single-Degree-of-Freedom System Using OpenSees:

Performs nonlinear dynamic analysis of a Single Degree of Freedom (SDOF) system subjected to ground motion, 
followed by 'fragility analysis' to assess the probability of exceeding predefined damage states.
 Key steps include:  
1. Model Setup: A 2D SDOF system with a bilinear material model is created using OpenSeesPy.  
2. Dynamic Analysis: Ground motion is applied, and the system's absolute acceleration response is computed.  
3. Fragility Analysis: Damage states (e.g., Slight, Moderate) are defined using lognormal distributions,
 and probabilities of exceedance are calculated for varying intensity measures (PGA).  
4. Visualization: The structural response and fragility curves are plotted for interpretation.  
This framework is robust, modular, and suitable for performance-based earthquake engineering applications.

Written By Salar Delavar Ghashghaei (Qashqai)
"""
import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import lognorm

# ================================
# 1. Structural Response Analysis
# ================================

def DYNAMIC_ANALYSIS(gm_accels, dt, k=200e3, m=100):
    
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    
    # Nodes - fixed base at (0,0), mass at (0,3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 0.0, 3.0)
    ops.fix(1, 1, 1, 1)
    ops.mass(2, m, 1e-9, 1e-9)  # Negligible rotational mass

    # Nonlinear material definition (Bilinear with 3% strain hardening)
    fy = 15000  # Yield force (N)
    E0 = k      # Initial stiffness
    b = 0.03    # Strain hardening ratio
    ops.uniaxialMaterial('Steel01', 1, fy, E0, b)

    # Zero-length element for nonlinear spring
    ops.element('zeroLength', 1, 1, 2, '-mat', 1, '-dir', 1)

    # Analysis configuration
    ops.timeSeries('Path', 1, '-dt', dt, '-values', *gm_accels.tolist(), '-factor', 9.81)
    ops.pattern('UniformExcitation', 1, 1, '-accel', 1)
    
    ops.system('BandSPD')
    ops.numberer('Plain')
    ops.constraints('Plain')
    ops.integrator('Newmark', 0.5, 0.25)
    ops.analysis('Transient')
    
    # Run analysis and collect results
    time = np.arange(0, len(gm_accels)*dt, dt)
    absolute_accels = []
    for step in range(len(gm_accels)):  # Use step counter instead of time value
        ops.analyze(1, dt)
        # Get total acceleration = structural + ground motion
        absolute_accels.append(ops.nodeAccel(2, 1) + gm_accels[step])  
        
    return time, np.array(absolute_accels), np.max(np.abs(absolute_accels))

# =================================
# 2. Fragility Analysis Framework
# =================================

class FragilityModel:
    def __init__(self, damage_states):
        """
        Initialize fragility model
        Args:
            damage_states (dict): {ds_name: (median_edp, beta)}
        """
        self.damage_states = damage_states
        
    def calculate_probability(self, im_values):
        """
        Calculate damage probabilities for IM values
        Args:
            im_values (np.array): Intensity measure values (PGA in g)
        Returns:
            dict: {ds_name: array of probabilities}
        """
        prob_dict = {}
        for ds, (median, beta) in self.damage_states.items():
            prob = lognorm(s=beta, scale=median).cdf(im_values)
            prob_dict[ds] = prob
        return prob_dict

# ========================
# 3. Application Example
# ========================

if __name__ == "__main__":
    # ----------------------------
    # Ground Motion Selection
    # ----------------------------
    gm_data = np.loadtxt('Ground_Acceleration_1.txt')  # Assumes acceleration in m/s²
    gm_data /= 9.81  # Convert to g units
    dt = 0.02  # Time step from record
    
    # ----------------------------
    # Structural Analysis
    # ----------------------------
    analysis_time, response_accels, pga = DYNAMIC_ANALYSIS(gm_data, dt)
    print(f"Peak Structural Acceleration: {pga:.2f}g")
    
    # ----------------------------
    # Fragility Assessment
    # ----------------------------
    # Define damage states per FEMA P-58
    damage_params = {
        'DS1_Slight': (0.15, 0.4),    # Median PGA=0.15g, β=0.4
        'DS2_Moderate': (0.30, 0.5),
        'DS3_Extensive': (0.60, 0.6),
        'DS4_Complete': (1.00, 0.7)
    }
    
    fragility_model = FragilityModel(damage_params)
    im_values = np.linspace(0.05, 2.0, 100)
    probabilities = fragility_model.calculate_probability(im_values)
    
    # ----------------------------
    # Visualization
    # ----------------------------
    plt.figure(figsize=(12, 6))
    
    # Response plot
    plt.subplot(1, 2, 1)
    plt.plot(analysis_time, response_accels, lw=1, color='black')
    plt.xlabel('Time (s)')
    plt.ylabel('Absolute Acceleration (g)')
    plt.title('Structural Response\nGround Motion')
    plt.grid(True)
    
    # Fragility curves
    plt.subplot(1, 2, 2)
    for ds, prob in probabilities.items():
        plt.plot(im_values, prob, lw=2, label=ds)
    plt.xlabel('Peak Ground Acceleration (g)')
    plt.ylabel('Probability of Exceedance')
    plt.title('Fragility Curves')
    plt.legend()
    plt.grid(True)
    
    plt.tight_layout()
    plt.show()