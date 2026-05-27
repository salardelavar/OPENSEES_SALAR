#%% EQUIVALENT SDOF SYSTEM DERIVATION VIA DISPLACEMENT-BASED SEISMIC DESIGN PROCEDURE WITH PUSHOVER ANALYSIS
# Change MDOF to SDOF System with Displacement Based Design Concept
# THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
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
displacement_X_1 = np.array(list(node_displacements.values())[0])   # DOF 02    
displacement_X_2 = np.array(list(node_displacements.values())[1])   # DOF 03 
displacement_X_3 = np.array(list(node_displacements.values())[2])   # DOF 04 
displacement_X_4 = np.array(list(node_displacements.values())[3])   # DOF 05 

ele_force_01 = np.array(list(node_reactions.values())[0])   # ELEMENT 01 
ele_force_02 = np.array(list(node_reactions.values())[1])   # ELEMENT 02 
ele_force_03 = np.array(list(node_reactions.values())[2])   # ELEMENT 03 
ele_force_04 = np.array(list(node_reactions.values())[3])   # ELEMENT 04 

STIFF_X_1 = np.abs(ele_force_01 / displacement_X_1)
STIFF_X_2 = np.abs(ele_force_02 / displacement_X_2)
STIFF_X_3 = np.abs(ele_force_03 / displacement_X_3)
STIFF_X_4 = np.abs(ele_force_04 / displacement_X_4)

MX2 = (MASS[0] * np.square(displacement_X_1) + 
       MASS[1] * np.square(displacement_X_2) + 
       MASS[2] * np.square(displacement_X_3) +
       MASS[3] * np.square(displacement_X_4))

MX = (MASS[0] * np.array(displacement_X_1) + 
      MASS[1] * np.array(displacement_X_2) + 
      MASS[2] * np.array(displacement_X_3) + 
      MASS[3] * np.array(displacement_X_4))
EFFECTIVE_DISP_X = MX2 / np.abs(MX)

EFFECTIVE_MASS_X = np.abs(MX) / EFFECTIVE_DISP_X


# Effective Stiffness
KX = (np.array(STIFF_X_1) * np.array(displacement_X_1) + 
      np.array(STIFF_X_2) * np.array(displacement_X_2) + 
      np.array(STIFF_X_3) * np.array(displacement_X_3) + 
      np.array(STIFF_X_4) * np.array(displacement_X_4))

EFFECTIVE_STIFF_X = np.abs(KX) / EFFECTIVE_DISP_X


# Effective Period
EFFECTIVE_PERIOD_X = 2 * np.pi / np.sqrt(EFFECTIVE_STIFF_X/EFFECTIVE_MASS_X)

print('Median Effective Displacement:           ', np.median(EFFECTIVE_DISP_X))
print('Median Effective Mass:                   ', np.median(EFFECTIVE_MASS_X))
print('Median Effective Stiffness:              ', np.median(EFFECTIVE_STIFF_X))
print('Median Effective Period:                 \n', np.median(EFFECTIVE_PERIOD_X))

# Create a figure with two subplots (Effective Displacement and Effective Mass)
fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(12, 10))

ax1.plot(EFFECTIVE_DISP_X, color='green', linewidth=3)
ax1.set_title(f'Effective Displacement - Median: {np.median(EFFECTIVE_DISP_X): .5f}')
ax1.set_xlabel('Step')
ax1.set_ylabel('Effective Displacement')
#ax1.legend(loc='upper right')
ax1.grid(True)

ax2.plot(EFFECTIVE_MASS_X, color='magenta', linewidth=3)
ax2.set_title(f'Effective Mass - Median: {np.median(EFFECTIVE_MASS_X): .5f}')
ax2.set_xlabel('Step')
ax2.set_ylabel('Effective Mass')
#ax2.legend(loc='upper right')
ax2.grid(True)

ax3.plot(EFFECTIVE_STIFF_X, color='cyan', linewidth=3)
ax3.set_title(f'Effective Stiffness - Median: {np.median(EFFECTIVE_STIFF_X): .5f}')
ax3.set_xlabel('Step')
ax3.set_ylabel('Effective Stiffness')
#ax3.legend(loc='upper right')
ax3.grid(True)

ax4.plot(EFFECTIVE_PERIOD_X, color='black', linewidth=3)
ax4.set_title(f'Effective Period - Median: {np.median(EFFECTIVE_PERIOD_X): .5f}')
ax4.set_xlabel('Step')
ax4.set_ylabel('Effective Period')
#ax4.legend(loc='upper right')
ax4.semilogy()
ax4.grid(True)

plt.tight_layout()
plt.show()