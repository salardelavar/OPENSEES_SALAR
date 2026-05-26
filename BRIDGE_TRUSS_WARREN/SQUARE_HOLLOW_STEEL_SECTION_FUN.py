def SQUARE_HOLLOW_STEEL_SECTION_FUN(
    secTag,
    matTag,
    B,      # [mm] outer width of the square section
    t,      # [mm] wall thickness
    nFib=20,
    plot=True
):
    import openseespy.opensees as ops
    import matplotlib.pyplot as plt
    """
    Fiber section definition for a Square Hollow Steel Section (SHS)
    and optional visualization of the section geometry.
    """

    # -------------------------
    # OpenSees fiber section
    # -------------------------
    ops.section('Fiber', secTag)

    # Outer coordinates
    y1 = -B / 2.0
    y2 =  B / 2.0
    z1 = -B / 2.0
    z2 =  B / 2.0

    # Inner (hollow) coordinates
    yi1 = y1 + t
    yi2 = y2 - t
    zi1 = z1 + t
    zi2 = z2 - t

    # Top plate
    ops.patch('quad', matTag, nFib, 1,
              y1,  z2,
              y2,  z2,
              yi2, zi2,
              yi1, zi2)

    # Bottom plate
    ops.patch('quad', matTag, nFib, 1,
              yi1, zi1,
              yi2, zi1,
              y2,  z1,
              y1,  z1)

    # Right plate
    ops.patch('quad', matTag, 1, nFib,
              yi2, zi1,
              y2,  z1,
              y2,  z2,
              yi2, zi2)

    # Left plate
    ops.patch('quad', matTag, 1, nFib,
              y1,  z1,
              yi1, zi1,
              yi1, zi2,
              y1,  z2)

    # -------------------------
    # Plot section geometry
    # -------------------------
    if plot:
        fig, ax = plt.subplots()

        # Outer square
        outer_y = [y1, y2, y2, y1, y1]
        outer_z = [z1, z1, z2, z2, z1]

        # Inner square (hole)
        inner_y = [yi1, yi2, yi2, yi1, yi1]
        inner_z = [zi1, zi1, zi2, zi2, zi1]

        ax.plot(outer_y, outer_z, 'k-', linewidth=2, label='Outer boundary')
        ax.plot(inner_y, inner_z, 'r--', linewidth=2, label='Inner boundary')

        ax.set_aspect('equal')
        ax.set_xlabel('y (mm)')
        ax.set_ylabel('z (mm)')
        ax.set_title('Square Hollow Steel Section (SHS)')
        ax.grid(True)
        ax.legend()

        plt.show()