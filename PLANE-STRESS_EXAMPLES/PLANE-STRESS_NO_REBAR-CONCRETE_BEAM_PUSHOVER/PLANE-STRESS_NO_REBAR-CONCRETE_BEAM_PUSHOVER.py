######################################################################################################################
#                                               IN THE NAME OF ALLAH                                                 #
#        NONLINEAR PLANE-STRESS PUSHOVER ANALYSIS OF A REINFORCED CONCRETE BEAM WITHOUT REBARS USING OPENSEES        #
#--------------------------------------------------------------------------------------------------------------------#
#                              THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                            #
#                                       EMAIL: salar.d.ghashghaei@gmail.com                                          #
######################################################################################################################
"""
Nonlinear Plane-stress Pushover Analysis of a Reinforced Concrete Beam without rebars Using OpenSeesPy

This script performs an advanced nonlinear pushover analysis of a simply supported reinforced concrete (RC)
 beam using OpenSeesPy, incorporating a 2D plane stress model with damage mechanics.
 The beam has dimensions of 4800 mm (length) × 600 mm (height) × 400 mm (thickness)
 and is discretized into a 20×15 triangular mesh using `tri31` elements.

The concrete material is modeled using the ASDConcrete3D nonlinear constitutive model,
 incorporating tensile and compressive damage mechanics. 
 Material properties are defined as follows:

 -Compressive strength: 30 MPa
 -Elastic modulus: E = 4700  * sqrt(fc')
 -Tensile strength: ft = fc' / 10

The analysis applies displacement-controlled pushover loading at the mid-span of the beam using
 a master node linked rigidly to key nodes at the top and bottom fiber. The loading proceeds up
 to a target displacement of 20 mm, with adaptive load stepping that automatically adjusts increments
 (minimum 0.0001 mm) to ensure convergence. Newton–Raphson iteration is employed for solving the nonlinear system.

Boundary Conditions:

 -Left bottom node: pinned support (fixed in both X and Y directions)
 -Right bottom and both top corner nodes: roller supports (fixed in Y only)

Output and Visualization:

The model tracks the full structural response, including:

 -Displacement history at the control node
 -Internal forces: axial force, shear, and bending moment at a critical mid-span element
 -Stress components from the damage model

All results are exported to CSV files for post-processing. A separate visualization routine is included to
 display the mesh, boundary conditions, and loading setup, providing clear graphical insight into the model behavior.

The script is modular and assumes the existence of an external `ANALYSIS_FUNCTION` responsible for convergence checking
 and step size control. This design allows for robust and customizable analysis loops, making the script well-suited for
 advanced structural engineers investigating the nonlinear behavior of RC beams under monotonic loading.
"""
import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
import ANALYSIS_FUNCTION as S02

# Define parameters (units: mm, N)
# ------------------------------------------
# Beam dimensions and mesh parameters
L = 4800.0    # [mm] Beam length
H = 600.0     # [mm] Beam height
B = 400       # [mm] Element thickness
Nx = 20       # Elements along length
Ny = 15       # Elements along height
#%%------------------------------------------------------
# Material properties (concrete with nonlinear behavior)
fc = 30.0             # [MPa] Compressive strength
v = 0.2               # Poisson's ratio
Ec = 4700*np.sqrt(fc) # [MPa] Elastic modulus
ft = fc/10.0          # [MPa] Tensile strength

# Pushover parameters
MAX_DISP = 20.0       # [mm] Maximum displacement
dU = 0.1              # [mm] Initial displacement increment
min_du = 1e-4         # [mm] Minimum displacement increment

# Define Analysis Properties
MAX_ITERATIONS = 20000     # Convergence iteration for test
MAX_TOLERANCE = 1.0e-4     # Convergence tolerance for test
#%%------------------------------------------------------
def build_model():
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 2)
        
    # Material definition with damage
    ec = fc/Ec  # Compressive yield strain
    et = ft/Ec  # Tensile yield strain
    MatTag = 1
    ops.nDMaterial('ASDConcrete3D', MatTag, Ec, v,
        '-Ce', 0.0, ec, ec+0.01,
        '-Cs', 0.0, fc, 0.1*fc,
        '-Cd', 0.0, 0.0, 1.0,
        '-Te', 0.0, et, et+0.01,
        '-Ts', 0.0, ft, 0.1*ft,
        '-Td', 0.0, 0.0, 1.0)
    # INFO LINK: https://opensees.github.io/OpenSeesDocumentation/user/manual/material/ndMaterials/ASDConcrete3D.html
    ops.nDMaterial('PlaneStress', 2, MatTag)

    # Create nodes
    node_id = 1
    nodes = {}
    for i in range(Nx + 1):
        for j in range(Ny + 1):
            x = i * L / Nx
            y = j * H / Ny
            ops.node(node_id, x, y)
            nodes[(i, j)] = node_id
            node_id += 1

    # Create elements
    ele_id = 1
    critical_element = None
    for i in range(Nx):
        for j in range(Ny):
            # Get node IDs for current quad
            n1 = nodes[(i, j)]
            n2 = nodes[(i, j+1)]
            n3 = nodes[(i+1, j+1)]
            n4 = nodes[(i+1, j)]
            
            # Split quad into two triangles
            # element tri31 $eleTag $iNode $jNode $kNode $thick $type $matTag <$pressure $rho $b1 $b2>
            # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/tri31.html
            ops.element('tri31', ele_id, n1, n2, n4, B, 'PlaneStress', 2)
            if i == Nx//2 - 1 and j == 0:  # Mid-span bottom element
                critical_element = ele_id
            ele_id += 1
            
            ops.element('tri31', ele_id, n2, n3, n4, B, 'PlaneStress', 2)
            ele_id += 1

    # Boundary conditions (simply supported)
    # Left support (pin)
    left_bottom = nodes[(0, 0)]
    left_top = nodes[(0, Ny)]
    ops.fix(left_bottom, 1, 1)  # Fixed x and y
    ops.fix(left_top, 0, 1)     # Fixed x only
    
    # Right support (roller)
    right_bottom = nodes[(Nx, 0)]
    right_top = nodes[(Nx, Ny)]
    ops.fix(right_bottom, 0, 1)  # Fixed x only
    ops.fix(right_top, 0, 1)     # Fixed x only
    
    # Create master node for pushover control
    master_node = node_id
    ops.node(master_node, L/2, H/2)
    
    # Identify mid-span nodes
    mid_bot = nodes[(Nx//2, 0)]
    mid_top = nodes[(Nx//2, Ny)]
    
    # Create rigid links to master node
    ops.rigidLink('bar', master_node, mid_bot)
    ops.rigidLink('bar', master_node, mid_top)
    
    return nodes, master_node, left_bottom, left_top, right_bottom, right_top, critical_element
#%%------------------------------------------------------
def RUN_PUSHOVER_ANALYSIS(MAX_DISP, dU):
    angle_deg = 90
    angle = angle_deg * np.pi/180.0  # Convert to radians
    
    # Build the model
    nodes, master_node, left_bot, left_top, right_bot, right_top, crit_elem = build_model()
    
    # Setup displacement control
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.sp(master_node, 1, np.cos(angle))
    ops.sp(master_node, 2, np.sin(angle))

    # Analysis settings
    ops.constraints('Transformation')
    ops.numberer('Plain')
    ops.system('FullGeneral')
    ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS, 0)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 0.0)  # Will be set during analysis
    ops.analysis('Static')
    
    current_disp = 0.0
    max_steps = int(np.abs(MAX_DISP)/np.abs(dU))
    # Results storage
    results = {
        'dispX': [],
        'dispY': [],
        'force_x': [],
        'force_y': [],
        'axial': [],
        'shear': [],
        'moment': [],
        'stress_x': [],
        'stress_y': []
    }
    
    print(f"Starting pushover analysis for {angle_deg}° direction")
    step = 0

    while current_disp < MAX_DISP and step < max_steps:
        # Set displacement increment
        ops.integrator('LoadControl', dU)
        OK = ops.analyze(1)
        S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS
        
        if OK == 0:
            # Successful step - record results
            step += 1
            current_disp = ops.getTime()
            Dx = ops.nodeDisp(master_node, 1)
            Dy = ops.nodeDisp(master_node, 2)
            
            # Get reaction forces at master node
            ops.reactions()
            Fx = ops.nodeReaction(master_node, 1)
            Fy = ops.nodeReaction(master_node, 2)
            
            # Get reactions at supports
            R_left_bottom = ops.nodeReaction(left_bot)
            R_left_top = ops.nodeReaction(left_top)
            
            # Calculate base reactions
            axial = R_left_bottom[0]
            shear = R_left_bottom[1] + R_left_top[1]
            moment = R_left_top[1] * (0.5*L)
            
            # Get stresses at critical element
            stress = ops.eleResponse(crit_elem, 'stress')

            sx = stress[0] if stress else 0.0
            sy = stress[1] if len(stress) > 1 else 0.0
            
            # Store results
            results['dispX'].append(Dx)
            results['dispY'].append(Dy)
            results['force_x'].append(Fx)
            results['force_y'].append(Fy)
            results['axial'].append(axial)
            results['shear'].append(shear)
            results['moment'].append(moment)
            results['stress_x'].append(sx)
            results['stress_y'].append(sy)
            
            # Print progress
            if step % 10 == 0:
                print(f"  Step {step}: DISP = {Dy:.3f} (mm), BASE-SHEAR = {shear:.2f} (N), BASE-AXIAL = {axial:.2f} (N)")
        else:
            # Reduce step size if analysis fails
            dU /= 2.0
            print(f"  Reduced step size to {dU:.4f} mm")
            if dU < min_du:
                print("  Minimum step size reached - stopping analysis")
                break
                
    print(f"Pushover analysis completed at {Dy:.2f} mm displacement")
    return results
#%%------------------------------------------------------
def PLOT_PUSHOVER_RESULTS(results):
    plt.figure(figsize=(15, 10))
    angle_deg = 90
    # Calculate force in push direction
    angle_rad = angle_deg * np.pi/180.0
    force_push = [
        fx*np.cos(angle_rad) + fy*np.sin(angle_rad) 
        for fx, fy in zip(results['force_x'], results['force_y'])
    ]
    
    
    # Base reactions
    plt.subplot(2, 2, 1)
    plt.plot(results['dispX'], results['axial'], 'r', label='Axial', markersize=10)
    plt.xlabel('Displacement in X Dir. (mm)')
    plt.ylabel('Base Axiaal Reaction')
    plt.title('Base Axial Reaction vs Displacement')
    plt.legend()
    plt.grid(True)
    
    plt.subplot(2, 2, 2)
    plt.plot(results['dispY'], results['shear'], 'g', label='Shear', markersize=10)
    plt.xlabel('Displacement in Y Dir. (mm)')
    plt.ylabel('Base Shear Reaction')
    plt.title('Base Shear Reaction vs Displacement')
    plt.legend()
    plt.grid(True)
    
    plt.subplot(2, 2, 3)
    plt.plot(results['dispY'], results['moment'], 'purple', label='Moment', markersize=10)
    plt.xlabel('Displacement in Y Dir. (mm)')
    plt.ylabel('Moment in the Middle of Beam')
    plt.title('Moment vs Displacement')
    plt.legend()
    plt.grid(True)
    
    # Stresses
    plt.subplot(2, 2, 4)
    plt.plot(results['dispY'], results['stress_x'], 'r-', label='Stress X', markersize=10)
    plt.plot(results['dispY'], results['stress_y'], 'b-', label='Stress Y', markersize=10)
    plt.xlabel('Displacement in Y Dir. (mm)')
    plt.ylabel('Stress (MPa)')
    plt.title('Critical Element Stresses')
    plt.legend()
    plt.grid(True)

    plt.tight_layout()
    plt.savefig(f'PLANE-STRESS_NO_REBAR-CONCRETE_BEAM_PUSHOVER.png', dpi=300)
    plt.show()
    
    # Save results to CSV
    import csv
    with open(f'PLANE-STRESS_NO_REBAR-CONCRETE_BEAM_PUSHOVER.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Disp_X (mm)', 'Disp_Y (mm)', 'Axial Base-reaction (N)', 'Shear Base-reaction (N)', 
                         'Moment in midle of beam (N.mm)', 'Stress_X (MPa)', 'Stress_Y (MPa)'])
        
        for i in range(len(results['dispX'])):
            writer.writerow([
                results['dispX'][i],
                results['dispY'][i],
                results['axial'][i],
                results['shear'][i],
                results['moment'][i],
                results['stress_x'][i],
                results['stress_y'][i]
            ])

#%%------------------------------------------------------
results = RUN_PUSHOVER_ANALYSIS(MAX_DISP, dU)
PLOT_PUSHOVER_RESULTS(results)

print("All pushover analyses completed!")

#%%------------------------------------------------------
def PLOT_BEAM():
    # Create nodes
    #angle = 90.0
    node_coords = []
    nodes = {}
    node_id = 1
    for i in range(Nx + 1):
        for j in range(Ny + 1):
            x = i * L / Nx
            y = j * H / Ny
            node_coords.append([x, y])
            nodes[(i, j)] = node_id
            node_id += 1
    
    # Create element connectivity for triangular elements
    elements = []
    ele_id = 1
    for i in range(Nx):
        for j in range(Ny):
            n1 = nodes[(i, j)]
            n2 = nodes[(i, j+1)]
            n3 = nodes[(i+1, j+1)]
            n4 = nodes[(i+1, j)]
            elements.append([n1-1, n2-1, n4-1])  # -1 for zero-based indexing
            ele_id += 1
            elements.append([n2-1, n3-1, n4-1])
            ele_id += 1
    
    # Convert to numpy arrays
    node_coords = np.array(node_coords)
    elements = np.array(elements)
    
    # Identify constrained nodes
    left_bottom = nodes[(0, 0)] - 1  # Zero-based indexing
    left_top = nodes[(0, Ny)] - 1
    right_bottom = nodes[(Nx, 0)] - 1
    right_top = nodes[(Nx, Ny)] - 1
    
    # Identify loaded nodes
    mid_bot = nodes[(Nx//2, 0)] - 1
    mid_top = nodes[(Nx//2, Ny)] - 1
    
    # Plot the mesh
    plt.figure(figsize=(12, 3))
    #plt.title(f'Beam Mesh with Constraints and Loads (Angle = {angle*180/np.pi:.1f}°)')
    plt.xlabel('X (mm)')
    plt.ylabel('Y (mm)')
    
    # Plot elements (triangles)
    for element in elements:
        x = node_coords[element, 0]
        y = node_coords[element, 1]
        x = np.append(x, x[0])
        y = np.append(y, y[0])
        plt.plot(x, y, 'b-', linewidth=0.5)
    
    # Plot nodes
    plt.plot(node_coords[:, 0], node_coords[:, 1], 'ro', markersize=2)
    
    # Plot constraints
    # Pinned support at left bottom (triangle)
    x, y = node_coords[left_bottom]
    plt.plot(x, y, 'g^', markersize=15, label='Pinned Support (Fixed X, Y)')
    
    # Roller supports (vertical arrows)
    for node in [left_top, right_bottom, right_top]:
        x, y = node_coords[node]
        plt.arrow(x, y-50, 0, 20, head_width=30, head_length=20, fc='g', ec='g') #
    plt.plot(node_coords[left_top, 0], node_coords[left_top, 1], 'gs', markersize=10, label='Roller Support (Fixed Y)')
    
    # Plot loads at mid-span
    arrow_scale = 200  # Scale factor for arrow visibility
    for node in [mid_bot, mid_top]:
        x, y = node_coords[node]
        #plt.arrow(x, y, 1*arrow_scale, 1*arrow_scale, head_width=50, head_length=50, fc='m', ec='m')
    plt.plot(node_coords[mid_bot, 0], node_coords[mid_bot, 1], 'm*', markersize=15, label='Applied Load')
    
    # Set aspect ratio, legend, and grid
    plt.axis('equal')
    plt.legend(loc='upper right', fontsize=8)
    plt.grid(True)
    plt.tight_layout()
    plt.show()

PLOT_BEAM()    
#%%------------------------------------------------------
