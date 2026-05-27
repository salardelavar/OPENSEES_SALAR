def PLOT_2D_QUAD(deformed_scale=1.0):
    """
    Plot undeformed and deformed shapes of a mesh with quad elements in 2D.
    Arguments:
        deformed_scale (float): Scaling factor for visualizing deformation shape
    """
    import openseespy.opensees as ops
    import numpy as np
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(12, 10))

    node_tags = ops.getNodeTags()
    node_coords = {node: ops.nodeCoord(node) for node in node_tags}

    ele_tags = ops.getEleTags()
    if not ele_tags:
        print("No elements to plot.")
        return
    first_ele = ele_tags[0]

    for ele in ele_tags:
        nodes = ops.eleNodes(ele)                 # list of node tags
        nodes_closed = nodes + [nodes[0]]         # close the polygon (list + list)

        # Undeformed shape
        x_undef = np.array([node_coords[node][0] for node in nodes_closed])
        y_undef = np.array([node_coords[node][1] for node in nodes_closed])
        ax.plot(x_undef, y_undef, 'k-', lw=1.2,
                label='Undeformed' if ele == first_ele else "")

        # Deformed shape
        ux = np.array([ops.nodeDisp(node)[0] for node in nodes_closed])
        uy = np.array([ops.nodeDisp(node)[1] for node in nodes_closed])
        ax.plot(x_undef + deformed_scale * ux,
                y_undef + deformed_scale * uy,
                'r--', lw=1.2,
                label='Deformed' if ele == first_ele else "")

    # Annotate node tags
    for node, (x, y) in node_coords.items():
        ux, uy = ops.nodeDisp(node)
        ax.text(x, y, f"{node}", color='blue', fontsize=9, ha='center', va='center')
        ax.text(x + deformed_scale * ux, y + deformed_scale * uy,
                f"{node}", color='purple', fontsize=9, ha='center', va='center')

    ax.set_xlabel('X [mm]')
    ax.set_ylabel('Y [mm]')
    ax.set_title(f'Undeformed vs Deformed Shape (Quad Elements)\nSCALE = {deformed_scale:.2f}')
    ax.grid(True)
    ax.set_aspect('equal', 'box')
    ax.legend()
    plt.tight_layout()
    plt.show()