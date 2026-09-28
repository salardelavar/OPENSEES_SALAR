def PLOT_1D_SPRING(deformed_scale=1.0,
                        virtual_spring_length=100.0,
                        virtual_spring_angle=45.0):
    """
    Plot the undeformed and deformed shapes of a 1D,
    including support for zero-length springs.

    Parameters
    ----------
    deformed_scale : float
        Scale factor for the deformed shape.
    virtual_spring_length : float
        Virtual length used to visualize zero-length springs.
    virtual_spring_angle : float
        Angle in degrees of the virtual spring length with respect to the X axis.
    """
    import openseespy.opensees as ops
    import numpy as np
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, figsize=(20, 16))

    # --- Node coordinates (first two dimensions only) ---
    nodes = ops.getNodeTags()
    node_coords = {n: np.array(ops.nodeCoord(n)[:2], dtype=float) for n in nodes}

    # --- Element data and zero-length detection ---
    ele_info = {}                 # ele -> (n1, n2, is_zero_length)
    zero_length_eles = []

    for ele in ops.getEleTags():
        n1, n2 = ops.eleNodes(ele)
        c1, c2 = node_coords[n1], node_coords[n2]
        is_zero = np.linalg.norm(c1 - c2) < 1e-8
        ele_info[ele] = (n1, n2, is_zero)
        if is_zero:
            zero_length_eles.append(ele)

    # --- Compute visualization offsets for zero-length spring nodes ---
    theta = np.deg2rad(virtual_spring_angle)
    direction = np.array([np.cos(theta), np.sin(theta)])

    node_offsets = {n: np.zeros(2) for n in nodes}

    # For each zero-length spring, move its two nodes apart by half the virtual length
    for ele in zero_length_eles:
        n1, n2, _ = ele_info[ele]
        node_offsets[n1] -= 0.5 * virtual_spring_length * direction
        node_offsets[n2] += 0.5 * virtual_spring_length * direction

    # Visualization coordinates (not the real coordinates)
    plot_coords = {n: node_coords[n] + node_offsets[n] for n in nodes}

    first_ele = ops.getEleTags()[0]

    # --- Plot undeformed shape ---
    for ele, (n1, n2, is_zero) in ele_info.items():
        c1, c2 = plot_coords[n1], plot_coords[n2]
        style = 'k-' if not is_zero else 'g-'   # zero-length springs in green
        lbl = 'Undeformed'
        if is_zero:
            lbl = 'Undeformed (zero-length spring)'

        ax.plot([c1[0], c2[0]], [c1[1], c2[1]], style, lw=1.2,
                label=lbl if ele == first_ele or (is_zero and ele == zero_length_eles[0]) else "")

    # --- Plot deformed shape ---
    for ele, (n1, n2, is_zero) in ele_info.items():
        c1, c2 = plot_coords[n1], plot_coords[n2]

        d1 = np.array(ops.nodeDisp(n1)[:2], dtype=float)
        d2 = np.array(ops.nodeDisp(n2)[:2], dtype=float)

        # Visualization offset + real displacement
        p1 = c1 + deformed_scale * d1
        p2 = c2 + deformed_scale * d2

        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], 'r--', lw=1.2,
                label='Deformed' if ele == first_ele else "")

    # --- Node labels ---
    for n in nodes:
        c = plot_coords[n]
        u = np.array(ops.nodeDisp(n)[:2], dtype=float)

        ax.text(c[0], c[1], f"{n}", color='blue', fontsize=10, ha='center',
                label='Node Tags' if n == nodes[0] else "")
        ax.text(c[0] + deformed_scale * u[0],
                c[1] + deformed_scale * u[1],
                f"{n}", color='purple', fontsize=10, ha='center')

    # --- Final plot settings ---
    ax.set_xlabel('X [mm]')
    ax.set_ylabel('Y [mm]')
    ax.set_title(
        f'Undeformed & Deformed Shapes — scale = {deformed_scale:.2f}  |  '
        f'virtual spring length = {virtual_spring_length:.1f} @ {virtual_spring_angle}°'
    )
    ax.legend(loc='best')
    ax.grid(True)
    ax.set_aspect('equal', adjustable='datalim')
    plt.tight_layout()
    plt.show()