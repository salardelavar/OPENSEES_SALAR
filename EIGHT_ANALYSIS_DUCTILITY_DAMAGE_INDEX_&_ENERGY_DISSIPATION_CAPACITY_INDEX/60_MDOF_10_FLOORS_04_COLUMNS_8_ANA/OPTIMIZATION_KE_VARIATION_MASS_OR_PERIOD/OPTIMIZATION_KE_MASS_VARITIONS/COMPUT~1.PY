#%% EQUIVALENT SDOF SYSTEM DERIVATION VIA DISPLACEMENT-BASED SEISMIC DESIGN PROCEDURE WITH PUSHOVER ANALYSIS
# Change MDOF to SDOF System with Displacement Based Design Concept
# THIS PYTHON SCRIPT IS WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
"""
This script implements a displacement-based pushover transformation,
 converting a multi-degree-of-freedom (MDOF) system into an equivalent
 single-degree-of-freedom (SDOF) system for seismic assessment.
 It calculates effective modal properties—displacement, mass, and
 stiffness—by weighting element forces and nodal displacements according
 to a presumed deformed shape.
 The derived effective period provides a simplified dynamic characteristic
 for performance-based engineering. The visualizations effectively track the
 evolution of these equivalent parameters throughout the nonlinear analysis steps.
"""
#---------------------------------------------------------------------------
import numpy as np 

displacement_X_1 = np.array(list(node_displacements.values())[1])   # DOF 02    
displacement_X_2 = np.array(list(node_displacements.values())[2])   # DOF 03 
displacement_X_3 = np.array(list(node_displacements.values())[3])   # DOF 04 
displacement_X_4 = np.array(list(node_displacements.values())[4])   # DOF 05 
displacement_X_5 = np.array(list(node_displacements.values())[5])   # DOF 06    
displacement_X_6 = np.array(list(node_displacements.values())[6])   # DOF 07 
displacement_X_7 = np.array(list(node_displacements.values())[7])   # DOF 08 
displacement_X_8 = np.array(list(node_displacements.values())[8])   # DOF 09 
displacement_X_9 = np.array(list(node_displacements.values())[9])   # DOF 10 
displacement_X_10 = np.array(list(node_displacements.values())[10])   # DOF 11 

ele_force_01 = np.array(list(ele_force.values())[0])   # ELEMENT 01 
ele_force_02 = np.array(list(ele_force.values())[1])   # ELEMENT 02 
ele_force_03 = np.array(list(ele_force.values())[2])   # ELEMENT 03 
ele_force_04 = np.array(list(ele_force.values())[3])   # ELEMENT 04 
ele_force_05 = np.array(list(ele_force.values())[4])   # ELEMENT 05 
ele_force_06 = np.array(list(ele_force.values())[5])   # ELEMENT 06 
ele_force_07 = np.array(list(ele_force.values())[6])   # ELEMENT 07 
ele_force_08 = np.array(list(ele_force.values())[7])   # ELEMENT 08 
ele_force_09 = np.array(list(ele_force.values())[8])   # ELEMENT 09 
ele_force_10 = np.array(list(ele_force.values())[9])   # ELEMENT 10 

STIFF_X_1 = np.abs(ele_force_01 / displacement_X_1)
STIFF_X_2 = np.abs(ele_force_02 / displacement_X_2)
STIFF_X_3 = np.abs(ele_force_03 / displacement_X_3)
STIFF_X_4 = np.abs(ele_force_04 / displacement_X_4)
STIFF_X_5 = np.abs(ele_force_05 / displacement_X_5)
STIFF_X_6 = np.abs(ele_force_06 / displacement_X_6)
STIFF_X_7 = np.abs(ele_force_07 / displacement_X_7)
STIFF_X_8 = np.abs(ele_force_08 / displacement_X_8)
STIFF_X_9 = np.abs(ele_force_09 / displacement_X_9)
STIFF_X_10 = np.abs(ele_force_10 / displacement_X_10)

MX2 = (MASS[0] * np.square(displacement_X_1) + 
       MASS[1] * np.square(displacement_X_2) + 
       MASS[2] * np.square(displacement_X_3) +
       MASS[3] * np.square(displacement_X_4) +
       MASS[4] * np.square(displacement_X_5) + 
       MASS[5] * np.square(displacement_X_6) +
       MASS[6] * np.square(displacement_X_7) + 
       MASS[7] * np.square(displacement_X_8) + 
       MASS[8] * np.square(displacement_X_9) +
       MASS[9] * np.square(displacement_X_10))

MX = (MASS[0] * np.array(displacement_X_1) + 
      MASS[1] * np.array(displacement_X_2) + 
      MASS[2] * np.array(displacement_X_3) + 
      MASS[3] * np.array(displacement_X_4) +
      MASS[4] * np.array(displacement_X_5) + 
      MASS[5] * np.array(displacement_X_6) + 
      MASS[6] * np.array(displacement_X_7) +
      MASS[7] * np.array(displacement_X_8) + 
      MASS[8] * np.array(displacement_X_9) + 
      MASS[9] * np.array(displacement_X_10))
EFFECTIVE_DISP_X = MX2 / np.abs(MX)

EFFECTIVE_MASS_X = np.abs(MX) / EFFECTIVE_DISP_X


# Effective Stiffness
KX = (np.array(STIFF_X_1) * np.array(displacement_X_1) + 
      np.array(STIFF_X_2) * np.array(displacement_X_2) + 
      np.array(STIFF_X_3) * np.array(displacement_X_3) + 
      np.array(STIFF_X_4) * np.array(displacement_X_4) + 
      np.array(STIFF_X_5) * np.array(displacement_X_5) + 
      np.array(STIFF_X_6) * np.array(displacement_X_6) + 
      np.array(STIFF_X_7) * np.array(displacement_X_7) + 
      np.array(STIFF_X_8) * np.array(displacement_X_8) + 
      np.array(STIFF_X_9) * np.array(displacement_X_9) + 
      np.array(STIFF_X_10) * np.array(displacement_X_10))

EFFECTIVE_STIFF_X = np.abs(KX) / EFFECTIVE_DISP_X


# Effective Period
EFFECTIVE_PERIOD_X = 2 * np.pi / np.sqrt(EFFECTIVE_STIFF_X/EFFECTIVE_MASS_X)


# Effective Damping Ratio
DRX = (damping_ratio_02 * np.array(displacement_X_1) + 
      damping_ratio_03 * np.array(displacement_X_2) + 
      damping_ratio_04 * np.array(displacement_X_3) + 
      damping_ratio_05 * np.array(displacement_X_4) +
      damping_ratio_06 * np.array(displacement_X_5) + 
      damping_ratio_07 * np.array(displacement_X_6) + 
      damping_ratio_08 * np.array(displacement_X_7) +
      damping_ratio_09 * np.array(displacement_X_8) + 
      damping_ratio_10 * np.array(displacement_X_9) + 
      damping_ratio_11 * np.array(displacement_X_10))

EFFECTIVE_DR_X = np.abs(DRX / EFFECTIVE_DISP_X)

MED_DISP = np.median(EFFECTIVE_DISP_X)
MED_MASS = np.median(EFFECTIVE_MASS_X)
MED_STIFF = np.median(EFFECTIVE_STIFF_X)
MED_PERIOD = np.median(EFFECTIVE_PERIOD_X)
MED_DR = np.median(EFFECTIVE_DR_X)

print('Median Effective Displacement:           ', MED_DISP)
print('Median Effective Mass:                   ', MED_MASS)
print('Median Effective Stiffness:              ', MED_STIFF)
print('Median Effective Period:                 ', MED_PERIOD)
print('Median Effective Damping Ratio:          ', MED_DR)

# Create a figure with two subplots (Effective Displacement and Effective Mass)
fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, 1, figsize=(18, 12))

ax1.plot(time, EFFECTIVE_DISP_X, color='green', linewidth=3)
ax1.set_title(f'Effective Displacement - Median: {np.median(EFFECTIVE_DISP_X): .5f}')
ax1.set_xlabel('Time [s]')
ax1.set_ylabel('Effective Displacement')
#ax1.legend(loc='upper right')
ax1.grid(True)

ax2.plot(time, EFFECTIVE_MASS_X, color='magenta', linewidth=3)
ax2.set_title(f'Effective Mass - Median: {np.median(EFFECTIVE_MASS_X): .5f}')
ax2.set_xlabel('Time [s]')
ax2.set_ylabel('Effective Mass')
#ax2.legend(loc='upper right')
ax2.grid(True)

ax3.plot(time, EFFECTIVE_STIFF_X, color='cyan', linewidth=3)
ax3.set_title(f'Effective Stiffness - Median: {np.median(EFFECTIVE_STIFF_X): .5f}')
ax3.set_xlabel('Time [s]')
ax3.set_ylabel('Effective Stiffness')
#ax3.legend(loc='upper right')
ax3.grid(True)

ax4.plot(time, EFFECTIVE_PERIOD_X, color='black', linewidth=3)
ax4.set_title(f'Effective Period - Median: {np.median(EFFECTIVE_PERIOD_X): .5f}')
ax4.set_xlabel('Time [s]')
ax4.set_ylabel('Effective Period')
#ax4.legend(loc='upper right')
#ax4.semilogy()
ax4.grid(True)

ax5.plot(time, EFFECTIVE_DR_X, color='red', linewidth=3)
ax5.set_title(f'Effective Damping Ratio - Median: {np.median(EFFECTIVE_DR_X): .2f} [%]')
ax5.set_xlabel('Time [s]')
ax5.set_ylabel('Effective Damping Ratio')
#ax5.legend(loc='upper right')
#ax5.semilogy()
ax5.grid(True)

plt.tight_layout()
plt.show()