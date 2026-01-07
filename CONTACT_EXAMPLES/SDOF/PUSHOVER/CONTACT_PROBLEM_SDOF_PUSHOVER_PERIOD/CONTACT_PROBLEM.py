'https://portwooddigital.com/2025/10/19/positive-opensees-contact/'
'https://openseespydoc.readthedocs.io/en/stable/src/zeroLengthContact2D.html'
import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt

# Define units
kip = 1
ft = 1
inch = ft/12.0

# Parameters
L = 8*ft
EI = 5120*kip*ft**2
P = 9*kip
gap = 0.6*inch

# Wipe and create model
ops.wipe()
ops.model('basic','-ndm',2,'-ndf',3)

# Create nodes
ops.node(1,0,gap); ops.fix(1,1,1,1)
ops.node(2,L,gap)
ops.node(3,0,0); ops.fix(3,1,1,0)
ops.node(4,L,0)
ops.node(5,2*L,0); ops.fix(5,0,1,0)

# Transformation
ops.geomTransf('Linear',1)

# Elements
ops.element('elasticBeamColumn',1,1,2,1,1,EI,1)
ops.element('elasticBeamColumn',2,3,4,1,1,EI,1)
ops.element('elasticBeamColumn',3,4,5,1,1,EI,1)

# Contact element
kn = 1e6*EI/L**3
ops.element('zeroLengthContact2D',4,4,2,kn,0,0,'-normal',0,-1)

# Loading
ops.timeSeries('Constant',1)
ops.pattern('Plain',1,1)
ops.load(2,0,-P,0)

# Analysis
ops.analysis('Static','-noWarnings')
ops.analyze(1)
ops.reactions()

# Get results
R3 = ops.nodeReaction(3,2)
R5 = ops.nodeReaction(5,2)
u2 = ops.nodeDisp(2,2)
u4 = ops.nodeDisp(4,2)

# Print results
print("Results:")
print(f"Vertical displacement at node 2 (u2): {u2:.6f} ft")
print(f"Vertical displacement at node 4 (u4): {u4:.6f} ft")
print(f"Reaction at node 3 (R3): {R3:.4f} kip")
print(f"Reaction at node 5 (R5): {R5:.4f} kip")

# Plot the deformed shape
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Plot 1: Structure geometry and deformation
ax1.set_title('Structure Geometry and Deformation')
ax1.set_xlabel('X (ft)')
ax1.set_ylabel('Y (ft)')
ax1.grid(True, alpha=0.3)
ax1.set_aspect('equal')

# Get node coordinates
nodes = [1, 2, 3, 4, 5]
x_coords = []
y_coords = []
for node in nodes:
    x_coords.append(ops.nodeCoord(node, 1))
    y_coords.append(ops.nodeCoord(node, 2))

# Get deformed coordinates
x_def = []
y_def = []
for node in nodes:
    x_def.append(ops.nodeCoord(node, 1) + ops.nodeDisp(node, 1))
    y_def.append(ops.nodeCoord(node, 2) + ops.nodeDisp(node, 2))

# Plot undeformed shape
ax1.plot([x_coords[0], x_coords[1]], [y_coords[0], y_coords[1]], 'k-', linewidth=2, label='Undeformed')
ax1.plot([x_coords[2], x_coords[3]], [y_coords[2], y_coords[3]], 'k-', linewidth=2)
ax1.plot([x_coords[3], x_coords[4]], [y_coords[3], y_coords[4]], 'k-', linewidth=2)

# Plot deformed shape (scaled for visibility)
scale_factor = 5  # Scale displacement for visualization
x_def_scaled = []
y_def_scaled = []
for node in nodes:
    x_def_scaled.append(ops.nodeCoord(node, 1) + scale_factor * ops.nodeDisp(node, 1))
    y_def_scaled.append(ops.nodeCoord(node, 2) + scale_factor * ops.nodeDisp(node, 2))

ax1.plot([x_def_scaled[0], x_def_scaled[1]], [y_def_scaled[0], y_def_scaled[1]], 'r--', linewidth=2, label=f'Deformed (x{scale_factor})')
ax1.plot([x_def_scaled[2], x_def_scaled[3]], [y_def_scaled[2], y_def_scaled[3]], 'r--', linewidth=2)
ax1.plot([x_def_scaled[3], x_def_scaled[4]], [y_def_scaled[3], y_def_scaled[4]], 'r--', linewidth=2)

# Plot nodes
ax1.scatter(x_coords, y_coords, color='blue', s=100, zorder=5, label='Nodes')
ax1.scatter(x_def_scaled, y_def_scaled, color='red', s=50, zorder=5, alpha=0.5, label='Deformed nodes')

# Add annotations
for i, node in enumerate(nodes):
    ax1.annotate(f'N{node}', (x_coords[i], y_coords[i]), xytext=(5, 5), textcoords='offset points')
    ax1.annotate(f'u={ops.nodeDisp(node,2):.3f}', (x_coords[i], y_coords[i]), 
                 xytext=(5, -15), textcoords='offset points', fontsize=8)

ax1.legend()
ax1.set_xlim(-1, 17)
ax1.set_ylim(-0.5, 1)

# Plot 2: Bar chart of reactions and displacements
ax2.set_title('Analysis Results')
categories = ['R3', 'R5', 'u2', 'u4']
values = [R3, R5, u2, u4]
units = ['kip', 'kip', 'ft', 'ft']

bars = ax2.bar(categories, values, color=['blue', 'green', 'orange', 'red'])
ax2.set_ylabel('Value')
ax2.grid(True, alpha=0.3, axis='y')

# Add value labels on bars
for bar, value, unit in zip(bars, values, units):
    height = bar.get_height()
    ax2.text(bar.get_x() + bar.get_width()/2., height,
             f'{value:.4f}\n{unit}',
             ha='center', va='bottom')

ax2.set_ylim(min(values)*1.2 if min(values) < 0 else 0, max(values)*1.2)

plt.tight_layout()
plt.show()

# Print summary
print("\nSummary:")
print("-" * 40)
print(f"Applied load P = {P} kip")
print(f"Initial gap = {gap:.4f} ft ({0.6} in)")
print(f"Total vertical reaction = {R3+R5:.4f} kip")
print(f"Load equilibrium check: Applied P = {P}, Reactions sum = {R3+R5:.4f}")
print(f"Difference: {abs(P - (R3+R5)):.6f} kip")