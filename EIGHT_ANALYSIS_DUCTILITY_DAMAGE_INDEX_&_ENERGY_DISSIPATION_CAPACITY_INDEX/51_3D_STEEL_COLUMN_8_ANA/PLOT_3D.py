def PLOT_3D_FRAME(deformed_scale=1.0):
    """
    Plot 3D frame geometry with undeformed (black) and deformed (red dashed) shapes from OpenSees model.
    
    Parameters
    ----------
    deformed_scale : float, optional
        Scale factor for visualizing deformations (defaults to 1.0)
        
    THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)    
    """
    import openseespy.opensees as ops
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    # Initialize figure and 3D axes
    fig = plt.figure(figsize=(14, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Gather node coordinates
    nodes = ops.getNodeTags()
    node_coords = {node: ops.nodeCoord(node) for node in nodes}

    # Plot undeformed shape (black solid lines)
    for i, ele in enumerate(ops.getEleTags()):
        node1, node2 = ops.eleNodes(ele)
        x1, y1, z1 = node_coords[node1]
        x2, y2, z2 = node_coords[node2]
        ax.plot([x1, x2], [y1, y2], [z1, z2],
                'k-', linewidth=2,
                label='Undeformed' if i == 0 else "")

    # Plot deformed shape (scaled, red dashed lines)
    for i, ele in enumerate(ops.getEleTags()):
        node1, node2 = ops.eleNodes(ele)
        x1, y1, z1 = node_coords[node1]
        x2, y2, z2 = node_coords[node2]
        ux1, uy1, uz1, _, _, _ = ops.nodeDisp(node1)
        ux2, uy2, uz2, _, _, _ = ops.nodeDisp(node2)
        ax.plot(
            [x1 + deformed_scale * ux1, x2 + deformed_scale * ux2],
            [y1 + deformed_scale * uy1, y2 + deformed_scale * uy2],
            [z1 + deformed_scale * uz1, z2 + deformed_scale * uz2],
            'r--', linewidth=2,
            label='Deformed' if i == 0 else ""
        )

    # Annotate node tags (undeformed & deformed)
    for node, (x, y, z) in node_coords.items():
        ux, uy, uz, _, _, _ = ops.nodeDisp(node)
        # Original position
        ax.text(x, y, z, f"{node}", color='blue', fontsize=10, ha='center', va='center')
        # Deformed position
        ax.text(x + deformed_scale * ux,
                y + deformed_scale * uy,
                z + deformed_scale * uz,
                f"{node}", color='purple', fontsize=10, ha='center', va='center')

    # Configure axes and plot appearance
    ax.set_xlabel('X [mm]', fontsize=12)
    ax.set_ylabel('Y [mm]', fontsize=12)
    ax.set_zlabel('Z [mm]', fontsize=12)
    ax.set_title(f"3D Frame: Undeformed vs Deformed (Scale = {deformed_scale:.2f})", fontsize=14)
    ax.legend()
    ax.grid(True)

    # Make aspect ratio close to equal (important for geometry visualization)
    all_coords = np.array(list(node_coords.values()))
    ranges = [all_coords[:, i].max() - all_coords[:, i].min() for i in range(3)]
    max_range = max(ranges)
    mid = [np.mean(all_coords[:, i]) for i in range(3)]
    ax.set_xlim(mid[0] - max_range / 2, mid[0] + max_range / 2)
    ax.set_ylim(mid[1] - max_range / 2, mid[1] + max_range / 2)
    ax.set_zlim(mid[2] - max_range / 2, mid[2] + max_range / 2)

    plt.tight_layout()
    plt.show()
