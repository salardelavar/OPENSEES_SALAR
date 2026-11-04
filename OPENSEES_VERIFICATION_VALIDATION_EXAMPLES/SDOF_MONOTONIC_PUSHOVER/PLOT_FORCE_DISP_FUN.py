def PLOT_FORCE_DISP_FUN(Displacement, Force):
    import numpy as np
    import matplotlib.pyplot as plt
    # Sort by displacement just in case
    idx = np.argsort(Displacement)
    D = Displacement[idx]
    F = Force[idx]
    # ---------- 1. Elastic stiffness ----------
    # Slope from origin to the first (yield) point
    K_elastic = F[0] / D[0]

    # ---------- 2. Plastic stiffness ----------
    # Slope of the plastic branch: from yield point to the final point
    K_plastic = (F[-1] - 0.0) / (D[-1] - 0.0)

    # ---------- 3. Tangent stiffness at each point ----------
    # Central differences for interior points, one-sided at the ends
    K_tangent = np.zeros_like(D, dtype=float)

    # Forward difference at first point
    K_tangent[0] = (F[1] - F[0]) / (D[1] - D[0])
    # Backward difference at last point
    K_tangent[-1] = (F[-1] - F[-2]) / (D[-1] - D[-2])
    # Central difference for interior points
    for i in range(1, len(D) - 1):
        K_tangent[i] = (F[i + 1] - F[i - 1]) / (D[i + 1] - D[i - 1])

    # ---------- Report ----------
    print(f"Elastic  stiffness  Ke = {K_elastic:.4f} ")
    print(f"Plastic  stiffness  Kp = {K_plastic:.4f} ")
    print("\nTangent stiffness at each point:")
    for d, f, kt in zip(D, F, K_tangent):
        print(f"  D = {d:6.2f} ,  F = {f:6.2f}  ->  Kt = {kt:8.4f}")

    # ---------- Plot ----------
    fig, ax1 = plt.subplots(figsize=(9, 6))

    # Force–Displacement curve
    ax1.plot(D, F, 'o-', color='steelblue', lw=2, ms=8, label='F–D curve')

    # Elastic line (origin -> first point)
    d_line = np.linspace(0, D[0], 20)
    ax1.plot(d_line, K_elastic * d_line, '--', color='green', label=f'Elastic stiffness = {K_elastic:.2f}')

    # Plastic line (first point -> last point)
    ax1.plot([0, D[-1]], [0, F[-1]], '--', color='red', label=f'Plastic stiffness = {K_plastic:.2f}')

    ax1.set_xlabel('Displacement')
    ax1.set_ylabel('Force')
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='upper right')
    plt.title('Force–Displacement: Elastic, Plastic and Tangent Stiffness')
    plt.tight_layout()
    plt.show()