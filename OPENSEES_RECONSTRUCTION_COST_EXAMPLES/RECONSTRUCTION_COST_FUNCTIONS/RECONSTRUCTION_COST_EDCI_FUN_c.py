
# ================================================================
# HYSTERETIC ENERGY AND ENERGY DISSIPATION CAPACITY INDEX (DECI)
# AUTHOR: SALAR DELAVAR GHASHGHAEI (QASHQAI)
# ================================================================

import numpy as np
import matplotlib.pyplot as plt


def DISSIPATED_ENERGY_FUN_WITH_PLOT(
    displacement,
    base_shear,
    method=2,
    XLABEL="Displacement [m]",
    YLABEL="Base Reaction [N]",
    title="HYSTERETIC RESPONSE",
    PLOT=True
):
    """
    HYSTERETIC WORK / ENERGY EVALUATION

    Parameters
    ----------
    displacement : array-like
        Chronological displacement response [m].

    base_shear : array-like
        Corresponding restoring-force or base-reaction response [N].

    method : int
        1 : Convex-hull area (geometric diagnostic only).
        2 : Force-displacement path integral (recommended).

    XLABEL, YLABEL : str
        Plot axis labels.

    title : str
        Figure title.

    PLOT : bool
        If True, generate the response plot.

    Returns
    -------
    energy : float
        Method 1: convex-hull area [J], not dissipated energy.
        Method 2: absolute net force-displacement work [J].

    fig : matplotlib.figure.Figure or None
        Figure object, or None if PLOT=False.

    Engineering note
    ----------------
    For a complete, closed hysteresis cycle, the absolute line
    integral of force with respect to displacement equals the
    enclosed loop area.

    For a non-closed earthquake response, the integral includes
    the net change in recoverable strain energy. It must not
    automatically be interpreted as pure hysteretic dissipation.

    A convex hull is not a valid measure of cumulative hysteretic
    energy.
    """

    # ------------------------------------------------------------
    # 1. INPUT VALIDATION
    # ------------------------------------------------------------
    disp = np.asarray(displacement, dtype=float).reshape(-1)
    force = np.asarray(base_shear, dtype=float).reshape(-1)

    if disp.size != force.size:
        raise ValueError(
            "Displacement and force arrays must have equal lengths."
        )

    if disp.size < 2:
        raise ValueError("At least two response points are required.")

    if not np.all(np.isfinite(disp)):
        raise ValueError("Displacement contains NaN or infinite values.")

    if not np.all(np.isfinite(force)):
        raise ValueError("Force contains NaN or infinite values.")

    if method not in (1, 2):
        raise ValueError("method must be 1 or 2.")

    # ------------------------------------------------------------
    # 2. ENERGY / WORK CALCULATION
    # ------------------------------------------------------------
    points = np.column_stack((disp, force))

    if method == 1:
        # Geometric envelope; NOT physical hysteretic dissipation.
        from scipy.spatial import ConvexHull, QhullError

        if np.unique(points, axis=0).shape[0] < 3:
            raise ValueError(
                "At least three distinct points are required "
                "for convex-hull evaluation."
            )

        try:
            hull = ConvexHull(points)
        except QhullError as exc:
            raise ValueError(
                "The response points do not define a valid "
                "two-dimensional convex hull."
            ) from exc

        energy = float(hull.volume)
        hull_pts = points[hull.vertices]

    else:
        # Chronological trapezoidal line integral: integral(F du).
        # No artificial closing segment is introduced.
        delta_disp = np.diff(disp)
        mean_force = 0.5 * (force[:-1] + force[1:])

        signed_work = float(np.sum(mean_force * delta_disp))
        energy = abs(signed_work)

    # ------------------------------------------------------------
    # 3. RESPONSE PLOT
    # ------------------------------------------------------------
    fig = None

    if PLOT:
        fig, ax = plt.subplots(figsize=(8, 6))

        ax.plot(
            disp,
            force,
            color="black",
            linewidth=1.2,
            label="Force-Displacement Response"
        )

        if method == 1:
            # Draw the hull without implying it is dissipated energy.
            closed_hull = np.vstack((hull_pts, hull_pts[0]))

            ax.plot(
                closed_hull[:, 0],
                closed_hull[:, 1],
                "r--",
                linewidth=1.8,
                label="Convex Hull (Geometric Envelope)"
            )

            ax.fill(
                closed_hull[:, 0],
                closed_hull[:, 1],
                color="red",
                alpha=0.15
            )

            ax.set_title(f"{title} — CONVEX HULL DIAGNOSTIC")

        else:
            ax.set_title(f"{title} — FORCE-DISPLACEMENT WORK")

            ax.scatter(
                [disp[0], disp[-1]],
                [force[0], force[-1]],
                s=35,
                zorder=3,
                label="Initial / Final Response"
            )

        ax.set_xlabel(XLABEL)
        ax.set_ylabel(YLABEL)
        ax.grid(True, linestyle="--", alpha=0.4)
        ax.legend()
        fig.tight_layout()

    return energy, fig


def RECONSTRUCTION_COST_EDCI_FUN(
    disp_SEI,
    reaction_SEI,
    disp_CP,
    reaction_CP,
    RECONSTRUCTION_COST=None,
    PLOT=True
):
    """
    COMPARE EARTHQUAKE-RESPONSE WORK WITH CYCLIC-TEST CAPACITY.

    SEI : Seismic / earthquake response history.
    CP  : Cyclic loading reference history.

    DECI = 100 * E_SEI / E_CP

    DECI is a relative energy-utilization index, not automatically
    a validated damage index or reconstruction-cost ratio.

    Returns
    -------
    DECI : float
        Energy utilization index [%].

    """

    # ------------------------------------------------------------
    # 1. EARTHQUAKE RESPONSE WORK
    # ------------------------------------------------------------
    E_SEI, fig_SEI = DISSIPATED_ENERGY_FUN_WITH_PLOT(
        displacement=disp_SEI,
        base_shear=reaction_SEI,
        method=2,
        XLABEL="Displacement [m]",
        YLABEL="Base Reaction [N]",
        title="EARTHQUAKE RESPONSE",
        PLOT=PLOT
    )

    if fig_SEI is not None:
        plt.show()

    print(f"\nEARTHQUAKE RESPONSE WORK = {E_SEI:,.3f} J")

    # ------------------------------------------------------------
    # 2. CYCLIC REFERENCE WORK
    # ------------------------------------------------------------
    E_CP, fig_CP = DISSIPATED_ENERGY_FUN_WITH_PLOT(
        displacement=disp_CP,
        base_shear=reaction_CP,
        method=2,
        XLABEL="Displacement [m]",
        YLABEL="Base Reaction [N]",
        title="CYCLIC REFERENCE RESPONSE",
        PLOT=PLOT
    )

    if fig_CP is not None:
        plt.show()

    print(f"CYCLIC REFERENCE WORK = {E_CP:,.3f} J")

    # ------------------------------------------------------------
    # 3. REFERENCE CAPACITY VALIDATION
    # ------------------------------------------------------------
    if E_CP <= np.finfo(float).eps:
        raise ValueError(
            "The cyclic reference work is zero or numerically "
            "negligible. DECI cannot be evaluated."
        )

    # ------------------------------------------------------------
    # 4. ENERGY DISSIPATION CAPACITY INDEX
    # ------------------------------------------------------------
    DECI = 100.0 * E_SEI / E_CP

    print("\n" + "=" * 55)
    print(" ENERGY DISSIPATION CAPACITY ASSESSMENT")
    print("=" * 55)
    print(f"Earthquake Response Work : {E_SEI:,.3f} J")
    print(f"Cyclic Reference Work    : {E_CP:,.3f} J")
    print(f"Energy Utilization Index : {DECI:.3f} %")

    if DECI < 100.0:
        print("Status: Below the reference work value.")
    elif np.isclose(DECI, 100.0):
        print("Status: Equal to the reference work value.")
    else:
        print("Status: Exceeds the reference work value.")
        print(
            "Check whether the cyclic reference history adequately "
            "represents the required deformation and loading protocol."
        )

    # ------------------------------------------------------------
    # 5. OPTIONAL ECONOMIC SCENARIO
    # ------------------------------------------------------------
    if RECONSTRUCTION_COST is not None:
        if (
            not np.isscalar(RECONSTRUCTION_COST)
            or not np.isreal(RECONSTRUCTION_COST)
            or not np.isfinite(RECONSTRUCTION_COST)
            or RECONSTRUCTION_COST < 0
        ):
            raise ValueError(
                "RECONSTRUCTION_COST must be a finite, "
                "non-negative real scalar."
            )

        # This is a scenario estimate, not a validated loss model.
        cost_scenario = (
            float(RECONSTRUCTION_COST) * DECI / 100.0
        )

        print(
            f"Linear Cost Scenario     : ${cost_scenario:,.2f}"
        )
        print(
            "Warning: The energy ratio alone does not establish "
            "actual reconstruction cost or structural damage."
        )

    print("=" * 55 + "\n")

    return DECI


"""
DECI = RECONSTRUCTION_COST_EDCI_FUN(
    disp_SEI=disp_SEI,
    reaction_SEI=reaction_SEI,
    disp_CP=disp_CP,
    reaction_CP=reaction_CP,
    RECONSTRUCTION_COST=1_000_000,
    PLOT=True
)

print(f"DECI = {DECI:.3f} %")
"""