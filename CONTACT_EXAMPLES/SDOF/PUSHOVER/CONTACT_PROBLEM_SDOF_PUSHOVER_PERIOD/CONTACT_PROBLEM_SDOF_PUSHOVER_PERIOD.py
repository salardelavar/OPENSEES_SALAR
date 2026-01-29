######################################################################################################################
#                          >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                          #
#         CONTACT-DRIVEN PUSHOVER ANALYSIS OF INELASTIC SDOF SYSTEMS: MONITORING PERIOD SHIFTS DURING                #
#                                       SECONDARY SPRING ACTIVATION IN OPENSEES                                      #
#--------------------------------------------------------------------------------------------------------------------#
#                              THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                            #
#                                       EMAIL: salar.d.ghashghaei@gmail.com                                          #
######################################################################################################################
"""
This script simulates the nonlinear pushover response of a single-degree-of-freedom system with a
 contact/gap mechanism. The structure has a primary spring (elastic or hysteretic) that activates
 immediately, while a secondary parallel spring engages only when displacement exceeds a specified 
 gap distance. This models structural components that come into contact only after certain
 deformation thresholds, such as gap-opening in masonry infills, pounding between adjacent structures,
 or secondary bracing systems activating during strong seismic events.

The analysis tracks force-displacement response, stiffness degradation, and period elongation
 as damage accumulates. The eigenvalue analysis at each step captures how the natural period
 increases with structural softening, a critical indicator of seismic vulnerability during
 progressive damage. Contact activation causes a sudden stiffness increase when the gap closes,
 followed by further period evolution as the system yields.
"""
import openseespy.opensees as ops
import ANALYSIS_FUNCTION as S01
import EIGENVALUE_ANALYSIS_FUN as S02
import numpy as np
import matplotlib.pyplot as plt
import time as TI

#%%------------------------------------------------------------------------------------------------
# Define  Structural Properties
FY = 85000.0                                     # [N] Yield Force of Structure
FU = 1.5 * FY                                    # [N] Ultimate Force of Structure
Ke = 4500000.0                                   # [N/m] Spring Elastic Stiffness
DY = FY / Ke                                     # [m] Yield Displacement
DSU = 0.36                                       # [m] Ultimate Displacement
Ksh = (FU - FY) / (DSU - DY)                     # [N/m] Displacement Hardening Modulus
Kp = FU / DSU                                    # [N/m] Spring Plastic Stiffness
b = Ksh / Ke                                     # Displacement Hardening Ratio

M = 150000.0                                     # [kg] Mass
zi = 0.05                                        # Damping ratio
#%%------------------------------------------------------------------------------------------------
# Positive branch points
pos_disp = [0, DY, DSU, 1.1*DSU, 1.25*DSU]
pos_force = [0, FY, FU, 0.2*FU, 0.1*FU]
KP = np.array([FY, DY, FU, DSU, 0.2*FU, 1.1*DSU, 0.1*FU, 1.25*DSU])

# Negative branch points
neg_disp = [0, -DY, -DSU, -1.1*DSU, -1.25*DSU]
neg_force = [0, -FY, -FU, -0.2*FU, -0.1*FU]
KN = np.array([-FY, -DY, -FU, -DSU, -0.2*FU, -1.1*DSU, -0.1*FU, -1.25*DSU])

# Plot
plt.plot(pos_disp, pos_force, marker='o', color='red')
plt.plot(neg_disp, neg_force, marker='o', color='black')

plt.xlabel("Displacement [m]")
plt.ylabel("Force [N]")
plt.title("Force–Displacement Curve")
plt.grid(True)
plt.axhline(0, linewidth=0.5)
plt.axvline(0, linewidth=0.5)
plt.show()
#%% DEFINE ANALYSIS PROPERTIES
MAX_ITERATIONS = 20000    # Convergence iteration for test
MAX_TOLERANCE = 1.0e-10   # Convergence tolerance for test
#%%------------------------------------------------------------------------------------------------
def PUSHOVER_CONTACT(SPRING, DMAX, DINCR, GAP_DIST):
    #%% Wipe and create model
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)
    #%% Create nodes (both at x=0 for zeroLength)
    ops.node(1, 0.0)
    ops.node(2, 0.0)
    
    ops.node(3, 0.0) # Contact Node
    
    #%% Define Constrain
    ops.fix(1, 1)
    
    ops.fix(3, 1) # Contact Node
    
    #%% Define Mass
    ops.mass(2, M) 
    
    Ke = FY/DY # [N/mm] Elastic Stiffness
    #%% Calculate natural frequency and damping coefficient
    omega = np.sqrt(Ke / M)
    C = 2 * zi * omega * M  # Damping Coefficient
    MatTag01, MatTag02 = 1, 2
    #%% Define uniaxial materials for each spring
    if SPRING == 'ELASTIC':
        ops.uniaxialMaterial('Elastic', MatTag01, Ke, C)    # Spring 1
        # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/ElasticUni.html
        # Define Spring 1 (always active)
        ops.element('zeroLength', 1, 1, 2, '-mat', MatTag01, '-dir', 1)
        # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/ZeroLength.html
    if SPRING == 'INELASTIC':
        ops.uniaxialMaterial('HystereticSM', MatTag01, '-posEnv', *KP.flatten(), '-negEnv', *KN.flatten(), '-pinch', 1, 1)
        # INFO LINK: https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/HystereticSM.html
        ops.uniaxialMaterial('Viscous', MatTag02, C, 1.0)  # Material for C (alpha=1.0 for linear)
        # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/Viscous.html
        # Define Spring 1 (always active)
        ops.element('zeroLength', 1, 1, 2, '-mat', MatTag01, MatTag02, '-dir', 1, 1)
            
    #%% Create a dummy load pattern for displacement control
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(2, 1.0)
        
    #%% Displacement control integrator (node 2, dof 1, increment = DINCR)
    ops.integrator('DisplacementControl', 2, 1, DINCR) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/integrator.html
    ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/test.html
    ops.algorithm('Newton')   # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/algorithm.html
    ops.system('BandGeneral') # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/system.html
    ops.numberer('Plain')     # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/numberer.html
    ops.constraints('Plain')  # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/constraints.html
    ops.analysis('Static')    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/analysis.html
    
    
    spring2_added = False
    disp = []; force1 = []; force2 = []; total_reaction = [];
    PERIOD_MIN, PERIOD_MAX = [], []
    stiffness = []
    
    # Initialize lists for each node's displacement
    node_displacements = {
        2: [],  # DISP02
    }
    
    node_reactions = {
        2: [],  # REACTION 02
    }
    Nsteps =  int(np.abs(DMAX/ DINCR))
    #print(Nsteps)
    
    for step in range(Nsteps):
        #%% Perform one load increment
        OK = ops.analyze(1)
        S01.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS
        u = ops.nodeDisp(2,1); disp.append(u)
        # Get force in Spring 1 (element 1)
        f1 = ops.eleResponse(1, 'force')[1]    # queries 'force' of zeroLength
        force1.append(f1) # Internal Force for element 1
        
        #%% Check activation condition
        if (u >= GAP_DIST) and not spring2_added:
            # Create Spring 2 in parallel (element tag 2)
            if SPRING == 'ELASTIC':
                ops.element('zeroLength', 2, 3, 2, '-mat', MatTag01, '-dir', 1)
                ops.domainChange()    # update the solver with new element
            if SPRING == 'INELASTIC':    
                ops.element('zeroLength', 2, 3, 2, '-mat', MatTag01, MatTag02, '-dir', 1, 1)
                ops.domainChange()    # update the solver with new element
            # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/domainChange.html
            spring2_added = True
            
        #%% If Spring 2 is active, record its force
        if spring2_added:
            f2 = ops.eleResponse(2, 'force')[0]
        else:
            f2 = 0.0
        force2.append(f2) # Internal Force for element 2
        
        total_reaction.append(np.abs((f1+f2))) # Total Reaction Force
        stiffness.append(np.abs(total_reaction[-1]) / np.abs(disp[-1])) # Structural Stiffness During Analysis
        #%% IN EACH STEP, STRUCTURE PERIOD GOING TO BE CALCULATED
        PERIODmin, PERIODmax = S02.EIGENVALUE_ANALYSIS(1, PLOT=True)
        PERIOD_MIN.append(PERIODmin)
        PERIOD_MAX.append(PERIODmax)
        
        # Store displacements and reacions
        for node_id in node_displacements.keys():
            node_displacements[node_id].append(ops.nodeDisp(node_id, 1))
            node_reactions[node_id].append(total_reaction[-1])  # Reaction force
            
        print(step+1, disp[-1], total_reaction[-1])
    else:
        print('Analysis is Finish Successfully \n\n')
        
    # Compute modal properties
    ops.modalProperties("-print", "-file", "SALAR_ModalReport.txt", "-unorm") 
    
    DATA = (disp, force1, force2, total_reaction, stiffness,
            np.array(PERIOD_MIN), np.array(PERIOD_MAX),
            node_displacements, node_reactions)
    return DATA    

#%%------------------------------------------------------------------------------------------------
# Analysis Durations for Static Analysis:
starttime = TI.process_time()

SPRING = 'INELASTIC'
DMAX = 0.35            # [m] Max. Pushover Incremental Displacement
DINCR = 0.001          # [m] Pushover Increment
GAP_DIST = 0.25        # [m] Gap Distance (Conact Distance)

DATA = PUSHOVER_CONTACT(SPRING, DMAX, DINCR, GAP_DIST)
(disp, force1, force2, total_reaction, stiffness,
 period_min, period_max,
 node_displacements, node_reactions) = DATA

totaltime = TI.process_time() - starttime
print(f'\nTotal Analysis Durations (s): {totaltime:.4f} \n\n')
#%%------------------------------------------------------------------------------------------------
plt.figure(1, figsize=(12, 8))
plt.plot(disp, force1, label='Spring 1', color='black', linewidth=4)
plt.plot(disp, force2, label='Spring 2', color='red', linewidth=4)
plt.axvline(GAP_DIST, color='purple', label='Gap Distance', linestyle='--', linewidth=4)
plt.xlabel('Displacement (Node 2) [m]')
plt.ylabel('Spring Force [N]')
plt.title('Elements Force vs Displacement During Pushover Analysis')
plt.legend(); plt.grid(True); plt.show()
#%%------------------------------------------------------------------------------------------------
plt.figure(2, figsize=(12, 8))
plt.plot(disp, total_reaction, color='purple', linewidth=4)
plt.xlabel('Displacement (Node 2) [m]')
plt.ylabel('Total Reaction Spring [N]')
plt.legend(); plt.grid(True); plt.show()
#%%------------------------------------------------------------------------------------------------
# PLOT STRUCTURAL PERIOD DURING THE ANALYSIS
plt.figure(0, figsize=(12, 8))
plt.plot(disp, period_min, color='black', linewidth=4)
plt.plot(disp, period_max, color='red', linewidth=4)
plt.title('Period of Structure vs Displacement During Pushover Analysis')
plt.ylabel('Structural Period [s]')
plt.xlabel('Displacement [m]')
#plt.semilogy()
plt.grid()
plt.legend([f'PERIOD - MIN VALUES: Min: {np.min(period_min):.3f} (s) - Mean: {np.mean(period_min):.3f} (s) - Max: {np.max(period_min):.3f} (s)', 
            f'PERIOD - MAX VALUES:  Min: {np.min(period_max):.3f} (s) - Mean: {np.mean(period_max):.3f} (s) - Max: {np.max(period_max):.3f} (s)',
            ])
plt.show()

plt.figure(1, figsize=(12, 8))
plt.plot(100*np.array(disp)/DSU, period_min, color='black', linewidth=4)
plt.plot(100*np.array(disp)/DSU, period_max, color='red', linewidth=4)
plt.title('Period of Structure vs Damage Level During Pushover Analysis')
plt.ylabel('Structural Period [s]')
plt.xlabel('Damage Level [%]')
#plt.xlabel('Displacement [m]')
#plt.semilogy()
plt.grid()
plt.legend([f'PERIOD - MIN VALUES: Min: {np.min(period_min):.3f} (s) - Mean: {np.mean(period_min):.3f} (s) - Max: {np.max(period_min):.3f} (s)', 
            f'PERIOD - MAX VALUES:  Min: {np.min(period_max):.3f} (s) - Mean: {np.mean(period_max):.3f} (s) - Max: {np.max(period_max):.3f} (s)',
            ])
plt.show()

plt.figure(2, figsize=(12, 8))
plt.plot(stiffness, period_min, color='black', linewidth=4)
plt.plot(stiffness, period_max, color='red', linewidth=4)
plt.title('Period of Structure vs Structural Stiffness During Pushover Analysis')
plt.ylabel('Structural Period [s]')
plt.xlabel('Structural Stiffness [N/m]')
#plt.semilogy()
plt.grid()
plt.legend([f'PERIOD - MIN VALUES: Min: {np.min(period_min):.3f} (s) - Mean: {np.mean(period_min):.3f} (s) - Max: {np.max(period_min):.3f} (s)', 
            f'PERIOD - MAX VALUES:  Min: {np.min(period_max):.3f} (s) - Mean: {np.mean(period_max):.3f} (s) - Max: {np.max(period_max):.3f} (s)',
            ])
plt.show()
#%%------------------------------------------------------------------------------
# Plotting Nodal Base-Reaction and Displacements
def PLOT_BASEREACTION_DISPLACEMENT(displacements_dict, reactions_dict, TITLE):
    plt.figure(figsize=(10, 6))
    
    for node_id, disp_values in displacements_dict.items():
        # Check if node exists in reactions dictionary and both arrays have same length
        if node_id in reactions_dict and len(disp_values) == len(reactions_dict[node_id]):
            plt.plot(disp_values, reactions_dict[node_id], 
                    label=f'Node {node_id} - MAX. ABS. : {np.max(np.abs(disp_values)): 0.4e}', 
                    linewidth=2)
        else:
            print(f"Warning: Node {node_id} not found in reactions dict or array length mismatch")
    
    plt.xlabel('Displacement [m]')
    plt.ylabel('Base Reaction [N]')
    plt.title(TITLE)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()


#PLOT_BASEREACTION_DISPLACEMENT(node_displacementsE, node_reactionsE, TITLE="Node Displacements vs Base Reactions for Elastic Structure") # ELASTIC STRUCTURE
PLOT_BASEREACTION_DISPLACEMENT(node_displacements, node_reactions, TITLE="Node Displacements vs Base Reactions for Inelastic Structure") # INELASTIC STRUCTURE
#%%------------------------------------------------------------------------------
# Print out the state of all nodes
ops.printModel("node",1, 2)
# Print out the state of all elements
ops.printModel("ele", 1)
# Print the Model
#printModel()
ops.printModel("-JSON", "-file", "CONTACT_PROBLEM_SDOF_PUSHOVER_PERIOD.json")
#%%-------------------------------------------------------------------------------