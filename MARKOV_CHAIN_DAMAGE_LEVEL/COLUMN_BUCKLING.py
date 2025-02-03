import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt

# -------------------------------
# Model: Post-Buckling of a Column
# -------------------------------

# Clear any previous model
ops.wipe()

# Define model dimensions: 2D frame (ndm=2, ndf=3)
ops.model('basic', '-ndm', 2, '-ndf', 3)

# Geometry and material properties
L  = 120.0      # column height in inches
E  = 29000.0    # modulus of elasticity (ksi)
A  = 10.0       # cross-sectional area (in^2)
I  = 500.0      # moment of inertia (in^4)
rho = 0.1       # mass per unit length (lb/in)

# Create nodes: bottom fixed, top free
ops.node(1, 0.0, 0.0)
ops.node(2, 0.0, L)

# Apply boundary conditions: fix node 1 in all DOFs
ops.fix(1, 1, 1, 1)

# Define geometric transformation
ops.geomTransf('Linear', 1)

# Define a uniaxial material for flexural behavior.
# (Here we use an elastic-perfectly plastic material to allow for plastic hinges.)
Fy = 5.0       # yield moment (kip-in)
ops.uniaxialMaterial('Steel01', 1, Fy, E, 0.02)

# Define a fiber section using a simplified approach:
# For demonstration, we use a layered section approximated with two fibers.
# (A more refined fiber discretization is recommended for accurate post-buckling response.)
ops.section('Fiber', 1)
# Add fibers: one for tension and one for compression, each representing half the area
ops.patch('rect', 1, 1, 1, -A/2.0, -L/20.0, A/2.0, L/20.0)
# Note: In practice, you may define a more detailed fiber section.

# Define the beam-column element using force-based formulation (which is good for inelastic buckling)
# Use forceBeamColumn element with 5 integration points
ops.geomTransf('Corotational', 1)
numIntPts = 5
ops.element('forceBeamColumn', 1, 1, 2, 1, numIntPts, 1)

# Define self-weight (lumped at the nodes) if desired (optional)
ops.mass(2, rho*L, 0.0, 0.0)

# -------------------------------
# Analysis: Nonlinear static (load-control) analysis
# -------------------------------

# Apply a vertical compressive load at the top node (node 2) in the y-direction.
P = -20.0   # applied load in kips (compressive)

# Define a plain load pattern (for static analysis)
ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)
ops.load(2, 0.0, P, 0.0)

# Set up recorders to capture displacement at node 2
disp_rec_file = 'node2_disp.out'
ops.recorder('Node', '-file', disp_rec_file, '-time', '-node', 2, '-dof', 2, 'disp')

# Define analysis parameters
ops.system('BandGeneral')
ops.numberer('RCM')
ops.constraints('Plain')
ops.integrator('DisplacementControl', 2, 2, -0.05)  # control displacement in y DOF of node 2, step size -0.05 in
ops.algorithm('Newton')
ops.test('NormDispIncr', 1e-6, 10)
ops.analysis('Static')

# Perform analysis: apply load steps to push beyond buckling
numSteps = 80
timeStep = 1.0  # pseudo-time step

for i in range(numSteps):
    ops.analyze(1, timeStep)
    currentDisp = ops.nodeDisp(2, 2)
    print(f"Step {i+1}, node 2 vertical disp = {currentDisp:.4f}")
    # Check if large post-buckling displacement reached
    if abs(currentDisp) > 10.0:
        print("Post-buckling regime reached.")
        break

# -------------------------------
# Post-Processing: Plotting deformed shape
# -------------------------------

# Get original coordinates and deformed coordinates of nodes
# For simplicity, we assume only node 2 moves (column pivot)
x0 = 0.0; y0 = 0.0   # Node 1 (fixed)
xL = 0.0; yL = L     # Node 2 original

# Get displacement at node 2
disp2 = ops.nodeDisp(2, 1), ops.nodeDisp(2, 2)
scaleFactor = 10.0  # exaggeration factor for plotting

# Deformed coordinates (only add horizontal displacement to illustrate buckling)
xL_def = xL + scaleFactor * disp2[0]
yL_def = yL + scaleFactor * disp2[1]

plt.figure(figsize=(6,8))
# Plot original column (line from node1 to node2)
plt.plot([x0, xL], [y0, yL], 'k--', label='Original')

# Plot deformed column shape
plt.plot([x0, xL_def], [y0, yL_def], 'r-', linewidth=2, label='Deformed (post-buckled)')

plt.xlabel('Horizontal displacement (inches)')
plt.ylabel('Vertical position (inches)')
plt.title('Post-Buckling Analysis of a Column')
plt.legend()
plt.grid(True)
plt.axis('equal')
plt.show()

# Clean up
ops.wipe()
