def PLOT_2D_TRI31(deformed_scale=1.0):
    """
    Plot undeformed and deformed shapes of a mesh with tri31 elements in 2D.
    Arguments:
        deformed_scale (float): Scaling factor for visualizing deformation shape
    """
    import openseespy.opensees as ops
    import numpy as np
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(12, 10))

    # Node coordinates dictionary
    node_tags = ops.getNodeTags()
    node_coords = {node: ops.nodeCoord(node) for node in node_tags}

    # Extract undeformed and deformed coordinates for tri31 elements
    for ele in ops.getEleTags():
        node1, node2, node3 = ops.eleNodes(ele)
        x = np.array([node_coords[node][0] for node in [node1, node2, node3, node1]])
        y = np.array([node_coords[node][1] for node in [node1, node2, node3, node1]])

        # Undeformed shape (solid black lines)
        ax.plot(x, y, 'k-', lw=1.2, label='Undeformed' if ele == ops.getEleTags()[0] else "")

        # Deformed shape (scaled displacements)
        ux = np.array([ops.nodeDisp(node)[0] for node in [node1, node2, node3, node1]])
        uy = np.array([ops.nodeDisp(node)[1] for node in [node1, node2, node3, node1]])
        ax.plot(x + deformed_scale * ux, y + deformed_scale * uy, 'r--', lw=1.2, label='Deformed' if ele == ops.getEleTags()[0] else "")

    # Annotate node tags
    for node, (x, y) in node_coords.items():
        ux, uy = ops.nodeDisp(node)
        ax.text(x, y, f"{node}", color='blue', fontsize=9, ha='center', va='center')
        ax.text(x + deformed_scale * ux, y + deformed_scale * uy, f"{node}", color='purple', fontsize=9, ha='center', va='center')

    # Plot settings
    ax.set_xlabel('X [mm]')
    ax.set_ylabel('Y [mm]')
    ax.set_title(f'Undeformed vs Deformed Shape (tri31 Elements)\nSCALE = {deformed_scale:.2f}')
    ax.grid(True)
    ax.set_aspect('equal', 'box')
    ax.legend()
    plt.tight_layout()
    plt.show()
