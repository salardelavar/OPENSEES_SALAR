
def RECONSTRUCTION_COST_DI_FUN(
    DISP,
    DY,
    DSU,
    RECONSTRUCTION_COST,
    PRINT_REPORT=True
):
    """
    DISPLACEMENT-BASED STRUCTURAL DAMAGE AND RECONSTRUCTION COST ASSESSMENT

    Author:
        Salar Delavar Ghashghaei (Qashqai)

    Description:
        Evaluates a displacement-based damage index, classifies the
        corresponding structural damage state, and estimates the
        reconstruction cost using a linear damage-cost relationship.

    Parameters
    ----------
    DISP : float
        Maximum absolute structural displacement [m].

    DY : float
        Yield displacement capacity [m].

    DSU : float
        Ultimate displacement capacity [m].
        Must be greater than DY.

    RECONSTRUCTION_COST : float
        Reference reconstruction cost [$].
        Must be non-negative.

    PRINT_REPORT : bool, optional
        If True, prints the damage assessment report.

    Returns
    -------
    rc : int
        Estimated reconstruction cost [$].

    Notes
    -----
    Damage Index (DI):
        DI = 100 * (abs(DISP) - DY) / (DSU - DY)

    The damage index is bounded between 0% and 100%.

    The estimated reconstruction cost is:
        RC = RECONSTRUCTION_COST * DI / 100

    This linear cost model is a simplified engineering assumption,
    not a universally validated relationship between displacement
    damage and actual economic loss.
    """

    import numpy as np

    # ---------------------------------------------------------
    # 1. INPUT VALIDATION
    # ---------------------------------------------------------
    values = {
        "DISP": DISP,
        "DY": DY,
        "DSU": DSU,
        "RECONSTRUCTION_COST": RECONSTRUCTION_COST
    }

    for name, value in values.items():
        if not np.isscalar(value) or not np.isreal(value):
            raise TypeError(f"{name} must be a real scalar.")

        if not np.isfinite(value):
            raise ValueError(f"{name} must be finite.")

    DISP = float(DISP)
    DY = float(DY)
    DSU = float(DSU)
    RECONSTRUCTION_COST = float(RECONSTRUCTION_COST)

    if DY < 0.0:
        raise ValueError("DY must be non-negative.")

    if DSU <= DY:
        raise ValueError("DSU must be greater than DY.")

    if RECONSTRUCTION_COST < 0.0:
        raise ValueError(
            "RECONSTRUCTION_COST must be non-negative."
        )

    # ---------------------------------------------------------
    # 2. DISPLACEMENT-BASED DAMAGE INDEX
    # ---------------------------------------------------------
    DISP_ABS = abs(DISP)

    DI_RAW = 100.0 * (DISP_ABS - DY) / (DSU - DY)

    # Bound the damage index to [0, 100] percent.
    DI = float(np.clip(DI_RAW, 0.0, 100.0))

    # ---------------------------------------------------------
    # 3. STRUCTURAL DAMAGE-STATE CLASSIFICATION
    # ---------------------------------------------------------
    DAMAGE_STATES = [
        (0.0,  "NO DAMAGE"),
        (10.0, "VERY MINOR DAMAGE"),
        (20.0, "MINOR DAMAGE"),
        (30.0, "MODERATE-LOW DAMAGE"),
        (40.0, "MODERATE DAMAGE"),
        (50.0, "MODERATE-HIGH DAMAGE"),
        (60.0, "SEVERE-LOW DAMAGE"),
        (70.0, "SEVERE-MEDIUM DAMAGE"),
        (80.0, "SEVERE-HIGH DAMAGE"),
        (90.0, "VERY SEVERE DAMAGE"),
        (100.0, "FAILURE DAMAGE")
    ]

    if DI == 0.0:
        ZONE = 0
        DAMAGE_STATE = DAMAGE_STATES[0][1]
    else:
        ZONE = min(int(np.ceil(DI / 10.0)), 10)
        DAMAGE_STATE = DAMAGE_STATES[ZONE][1]

    # ---------------------------------------------------------
    # 4. RECONSTRUCTION COST ESTIMATION
    # ---------------------------------------------------------
    RC_RAW = RECONSTRUCTION_COST * DI / 100.0

    # Preserve full precision internally; return integer dollars.
    rc = int(round(RC_RAW))

    # ---------------------------------------------------------
    # 5. ENGINEERING ASSESSMENT REPORT
    # ---------------------------------------------------------
    if PRINT_REPORT:
        print("\n" + "=" * 58)
        print(" STRUCTURAL DAMAGE AND RECONSTRUCTION COST ASSESSMENT")
        print("=" * 58)

        print(f"Maximum Displacement  : {DISP_ABS:.6f} m")
        print(f"Yield Displacement    : {DY:.6f} m")
        print(f"Ultimate Displacement : {DSU:.6f} m")
        print(f"Damage Index          : {DI:.2f} %")
        print(f"Damage Zone           : {ZONE}")
        print(f"Damage State          : {DAMAGE_STATE}")
        print(f"Reference Cost        : ${RECONSTRUCTION_COST:,.2f}")
        print(f"Estimated Cost        : ${rc:,.0f}")

        print("=" * 58 + "\n")

    return rc


"""
DISP = 0.18           # Maximum displacement [m]
DY = 0.04             # Yield displacement [m]
DSU = 0.36            # Ultimate displacement [m]
RECONSTRUCTION_COST = 1_000_000  # [$]

RC = RECONSTRUCTION_COST_DI_FUN(
    DISP=DISP,
    DY=DY,
    DSU=DSU,
    RECONSTRUCTION_COST=RECONSTRUCTION_COST
)

print(f"Reconstruction Cost = ${RC:,.0f}")
"""
