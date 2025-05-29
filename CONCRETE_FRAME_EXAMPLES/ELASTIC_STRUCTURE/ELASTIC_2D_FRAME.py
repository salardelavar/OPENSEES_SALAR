###########################################################################################################
#                                         IN THE NAME OF ALLAH                                            #
#                GRAVITY AND CYCLIC STATIC ANALYSIS OF ELASTIC CONCRETE STRUCTURE                         # 
#---------------------------------------------------------------------------------------------------------#
#                         THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                      #
#                                      EMAIL: salar.d.ghashghaei@gmail.com                                #
###########################################################################################################
"""
 Model and analyze a 2D elastic concrete frame under gravity and cyclic lateral loads.
 It defines a simple frame with beam-column elements, material and geometric properties,
 and load patterns. Two separate analyses are performed: gravity loading and cyclic
 lateral loading, recording displacements and internal forces at each step.
 Results are plotted for displacements, rotations, and base reactions over time.
 Finally, the deformed shape of the frame is visualized using a custom plotting function.
"""
import openseespy.opensees as ops
import ANALYSIS_FUNCTION as S02
import PLOT_2D as S04
# Clear any existing model
ops.wipe()

# Create a 2D model with 3 DOF per node
ops.model('Basic', '-ndm', 2, '-ndf', 3)
Bc = 500     # [mm] Column Section Width
Hc = 500     # [mm] Column Section Height
Bb = 500     # [mm] Beam Section Width
Hb = 300     # [mm] Beam Section Height
Ea = 2.1e5      # [N/mm^1] Young's Modulus
LENGTH_COL = 3000 # [mm] Column Length 
LENGTH_BM = 7000  # [mm] Beam Length 
PY = 200.0      # [N] Applied Load In Y-Direction
MZ = 160.0      # [N.mm] Applied Moment In Z-Direction
MASS = 57500    # [kg] Nodal Mass
DIS_LOAD = -0.3 #[N/mm] Distributed Load on Beam
numSteps = 1000
# Define Analysis Properties
MAX_ITERATIONS = 1000000   # Convergence iteration limit
MAX_TOLERANCE = 1.0e-10    # Convergence tolerance

GLA = 'False'
CLLA = 'True'

AREAc = Bc*Hc
IZc = Bc*(Hc**3) / 12
AREAb = Bb*Hb
IZb = Bb*(Hb**3) / 12

# Create nodes
ops.node(1, 0.0, 0.0)
ops.node(2, LENGTH_BM, 0.0)
ops.node(3, 0.0, LENGTH_COL)
ops.node(4, LENGTH_BM, LENGTH_COL)

# Fix base nodes (1 and 2) - fixed in all DOFs
ops.fix(1, 1, 1, 1)
ops.fix(2, 1, 1, 1)

# Add masses to nodes 3 and 4
ops.mass(3, MASS, MASS, 0.0)
ops.mass(4, MASS, MASS, 0.0)

# Geometric transformation
#ops.geomTransf('Linear', 1)
ops.geomTransf('Corotational', 1)


# Create elements
# Note: Corrected duplicate element tag (1 was used twice in original)
# element elasticBeamColumn $eleTag $iNode $jNode $A $E $Iz $transfTag <-mass $massDens> <-cMass>
ops.element('elasticBeamColumn', 1, 1, 3, AREAc, Ea, IZc, 1) # Column 01
ops.element('elasticBeamColumn', 2, 2, 4, AREAc, Ea, IZc, 1) # Column 02
ops.element('elasticBeamColumn', 3, 3, 4, AREAb, Ea, IZb, 1) # Beam 01

# Create a load pattern
ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)

# Apply loads
ops.load(3, 0.0, -PY, -MZ)
ops.load(4, 0.0, -PY, MZ)
# Apply distributed load to all beams in the structure
# eleLoad -ele $eleTag1 <$eleTag2 ....> -type -beamUniform $Wz <$Wx>
ops.eleLoad('-ele', 3, '-type', '-beamUniform', DIS_LOAD)

print("\nModel Done.")

#%% ------------------------------------------------------------------------------------
# Gravity Load Analysis
# Create analysis components
import matplotlib.pyplot as plt
import numpy as np
def GRAVITY_ANALYSIS(GLA=True, numSteps=10):
    if GLA:
        # Initialize data storage
        time_history = []
        disp_X, disp_Y, rot_Z = [], [], []
        base_shear, base_axial, base_moment = [], [], []
        
        # Create analysis components
        ops.constraints('Transformation')
        ops.numberer('Plain')
        ops.system('BandGeneral')
        ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS, 2)
        ops.algorithm('Newton')
        
        # Set up integrator and analysis type
        ops.integrator('LoadControl', 1/numSteps)
        ops.analysis('Static')
        
        # Perform the analysis in numSteps
        for i in range(numSteps):
            stable = ops.analyze(1)
            S02.ANALYSIS(stable, 1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS
            
            # Record responses at each step
            time_history.append(ops.getTime())
            disp_X.append(ops.nodeDisp(3, 1))  # X-displacement of node 3
            disp_Y.append(ops.nodeDisp(3, 2))  # Y-displacement of node 3
            rot_Z.append(ops.nodeDisp(3, 3))   # Rotation of node 3
            
            # Get element forces (columns 1 and 2)
            f1 = ops.eleResponse(1, 'force')
            f2 = ops.eleResponse(2, 'force')
            base_shear.append(-f1[0] - f2[0])    # Sum of X forces (shear)
            base_axial.append(-f1[1] - f2[1])    # Sum of Y forces (axial)
            base_moment.append(-f1[2] - f2[2])   # Sum of moments
        
        print('\nGravity Analysis Done. \n')    
        
        # Convert to numpy arrays
        time = np.array(time_history)
        dX = np.array(disp_X)
        dY = np.array(disp_Y)
        rZ = np.array(rot_Z)
        b_sh = np.array(base_shear)
        b_ax = np.array(base_axial)
        b_mom = np.array(base_moment)
                
        # Create plots
        plt.figure(figsize=(18, 12))
        
        # 1. Displacement Time History
        plt.subplot(3, 2, 1)
        plt.plot(time, dX, label='X-disp', color='black')
        plt.xlabel('Load Factor')
        plt.ylabel('Displacement (mm)')
        plt.title('Displacement-X Time History')
        plt.legend()
        plt.grid(True)
        
        # 2. Displacement Time History
        plt.subplot(3, 2, 3)
        plt.plot(time, dY, label='Y-disp', color='purple')
        plt.xlabel('Load Factor')
        plt.ylabel('Displacement (mm)')
        plt.title('Displacement-Y Time History')
        plt.grid(True)
        
        # 3. Rotation Time History
        plt.subplot(3, 2, 5)
        plt.plot(time, rZ)
        plt.xlabel('Load Factor')
        plt.ylabel('Rotation (rad)')
        plt.title('Rotation Time History')
        plt.grid(True)
        
        # 4. Base Shear Force History
        plt.subplot(3, 2, 2)
        plt.plot(time, b_sh)
        plt.xlabel('Load Factor')
        plt.ylabel('Base Shear Force (N)')
        plt.title('Base Shear Force History')
        plt.grid(True)
        
        # 5. Base Axial Force History
        plt.subplot(3, 2, 4)
        plt.plot(time, b_ax)
        plt.xlabel('Load Factor')
        plt.ylabel('Base Axial Force (N)')
        plt.title('Base Axial Force History')
        plt.grid(True)
        
        # 6. Base Moment Force History
        plt.subplot(3, 2, 6)
        plt.plot(time, b_mom)
        plt.xlabel('Load Factor')
        plt.ylabel('Base Moment Force (N)')
        plt.title('Base Moment Force History')
        plt.grid(True)
        
        plt.tight_layout()
        plt.show()
        
        return {
            'time': time,
            'disp_X': dX,
            'disp_Y': dY,
            'rot_Z': rZ,
            'shear': b_sh,
            'axial': b_ax,
            'moment': b_mom
        }

# Run the analysis
results = GRAVITY_ANALYSIS(GLA=True, numSteps=100)
#%% ------------------------------------------------------------------------------------
# Cyclic Lateral Load Analysis
# Set gravity loads constant and reset time
import matplotlib.pyplot as plt
import numpy as np

def CYCLIC_ANALYSIS(CLLA=True, px=1, numSteps=10):
    if CLLA:
        # Reset and set analysis parameters
        ops.loadConst('-time', 0.0)
        ops.constraints('Transformation')
        ops.numberer('Plain')
        ops.system('BandGeneral')
        ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS)
        ops.algorithm('Newton')
        ops.analysis('Static')
        
        # Lateral load pattern
        ops.pattern('Plain', 2, 1)
        ops.load(3, px, 0.0, 0.0)
        ops.load(4, px, 0.0, 0.0)
        
        # Initialize data storage
        time_history = []
        disp_X, disp_Y, rot_Z = [], [], []
        base_shear, base_axial, base_moment = [], [], []
        
        # Cyclic loading parameters
        cycles = [
            (numSteps, 1/numSteps),    # Push right
            (numSteps, -1/numSteps),   # Pull left
            (numSteps, -1/numSteps),   # Pull left
            (numSteps, 1/numSteps),    # Push right
            (numSteps, 1/numSteps)     # Push right
        ]
        
        # Perform cyclic analysis
        for n_steps, step_size in cycles:
            ops.integrator('LoadControl', step_size)
            for _ in range(n_steps):
                stable = ops.analyze(1)
                S02.ANALYSIS(stable, 1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS
                
                # Record responses
                time_history.append(ops.getTime())
                disp_X.append(ops.nodeDisp(3, 1))
                disp_Y.append(ops.nodeDisp(3, 2))
                rot_Z.append(ops.nodeDisp(3, 3))
                
                # Get element forces (columns 1 and 2)
                f1 = ops.eleResponse(1, 'force')
                f2 = ops.eleResponse(2, 'force')
                base_shear.append(f1[0] + f2[0])
                base_axial.append(f1[1] + f2[1])
                base_moment.append(f1[2] + f2[2])
                
        print('\nCyclic Analysis Done. \n')        
        
        # Convert to numpy arrays
        time = np.array(time_history)
        dX = np.array(disp_X)
        dY = np.array(disp_Y)
        rZ = np.array(rot_Z)
        b_sh = np.array(base_shear)
        b_ax = np.array(base_axial)
        b_mom = np.array(base_moment)
        
        # Create plots
        plt.figure(figsize=(18, 12))
        
        # 1. Displacement Time History
        plt.subplot(3, 2, 1)
        plt.plot(time, dX, label='X-disp', color='black')
        plt.xlabel('Load Factor')
        plt.ylabel('Displacement (mm)')
        plt.title('Displacement-X Time History')
        plt.legend()
        plt.grid(True)
        
        # 2. Displacement Time History
        plt.subplot(3, 2, 3)
        plt.plot(time, dY, label='Y-disp', color='purple')
        plt.xlabel('Load Factor')
        plt.ylabel('Displacement (mm)')
        plt.title('Displacement-Y Time History')
        plt.grid(True)
        
        # 3. Rotation Time History
        plt.subplot(3, 2, 5)
        plt.plot(time, rZ)
        plt.xlabel('Load Factor')
        plt.ylabel('Rotation (rad)')
        plt.title('Rotation Time History')
        plt.grid(True)
        
        # 4. Base Shear Force History
        plt.subplot(3, 2, 2)
        plt.plot(time, b_sh)
        plt.xlabel('Load Factor')
        plt.ylabel('Base Shear Force (N)')
        plt.title('Base Shear Force History')
        plt.grid(True)
        
        # 5. Base Axial Force History
        plt.subplot(3, 2, 4)
        plt.plot(time, b_ax)
        plt.xlabel('Load Factor')
        plt.ylabel('Base Axial Force (N)')
        plt.title('Base Axial Force History')
        plt.grid(True)
        
        # 6. Base Moment Force History
        plt.subplot(3, 2, 6)
        plt.plot(time, b_mom)
        plt.xlabel('Load Factor')
        plt.ylabel('Base Moment Force (N)')
        plt.title('Base Moment Force History')
        plt.grid(True)
        
        plt.tight_layout()
        plt.show()
        
        return {
            'time': time,
            'disp_X': dX,
            'disp_Y': dY,
            'rot_Z': rZ,
            'shear': b_sh,
            'axial': b_ax,
            'moment': b_mom
        }

# Run the analysis
results = CYCLIC_ANALYSIS(CLLA=True, px=10, numSteps=100)
#%% ------------------------------------------------------------------------------------ 
# %% Plot 2D Frame Shapes
S04.PLOT_2D_FRAME(deformed_scale=100)  # Adjust scale factor as needed
#%%------------------------------------------------------------------------------    
    
    
