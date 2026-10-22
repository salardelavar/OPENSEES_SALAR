# PUSHOVER ANALYSIS OF A MULTI-DEGREE-OF-FREEDOM STRUCTURE VIA EVALUATION OF A MULTILINEAR FITTING CURVE
This Python-based OpenSees framework performs nonlinear static pushover analysis of a ten-story, multi-degree-of-freedom (MDOF) structural system comprising four columns per floor, with the objective of investigating global lateral-force resistance and deformation capacity under progressively increasing lateral loading. By implementing a displacement-controlled analysis procedure, the model can be used to establish the structural capacity curve, characterize the transition from initial elastic response to nonlinear behavior, and examine the implications of stiffness degradation, yielding, and inelastic deformation for overall structural performance, subject to the constitutive models and response quantities implemented in the script. From an advanced structural engineering perspective, the framework provides a computational basis for evaluating global load–displacement behavior, identifying potential structural performance limitations, and supporting subsequent seismic assessment using appropriate ductility-based damage indices and energy-dissipation measures. The reliability of these evaluations depends on the adopted material models, element formulations, boundary conditions, gravity-load state, and convergence criteria, while the resulting pushover response should be interpreted as a nonlinear static approximation rather than a complete representation of earthquake-induced dynamic behavior.

# تحلیل پوش‌آور سازهٔ چنددرجه‌آزادی با ارزیابی منحنی برازش چندخطی
![alt text](https://github.com/salardelavar/OPENSEES_SALAR/blob/main/EIGHT_ANALYSIS_DUCTILITY_DAMAGE_INDEX_%26_ENERGY_DISSIPATION_CAPACITY_INDEX/60_MDOF_10_FLOORS_04_COLUMNS_8_ANA/PUSHOVER_MDOF_10_FLOORS_04_COLUMNS/COVER.png)

# EQUIVALENT SDOF SYSTEM DERIVATION VIA DISPLACEMENT-BASED SEISMIC DESIGN PROCEDURE WITH PUSHOVER ANALYSIS

![alt text](https://github.com/salardelavar/OPENSEES_SALAR/blob/main/EIGHT_ANALYSIS_DUCTILITY_DAMAGE_INDEX_%26_ENERGY_DISSIPATION_CAPACITY_INDEX/60_MDOF_10_FLOORS_04_COLUMNS_8_ANA/PUSHOVER_MDOF_10_FLOORS_04_COLUMNS/COVER_DISPLACEMENT_BASED_PUSHOVER.png)

This Python script performs nonlinear pushover analysis of a 10-story MDOF
 building frame using OpenSeesPy to derive an equivalent SDOF system for
 displacement-based seismic design.
1. It models a 10-DOF lumped-mass shear building with 4 parallel zeroLength 
springs per floor, each combining a hysteretic uniaxial material (elastic or inelastic)
 and a viscous damper.
 
2. Spring properties (yield force, ultimate force, elastic stiffness, ultimate displacement) 
are defined per column, from which yield displacement, hardening ratio, and damping coefficient are computed.

3. Multiple analysis types are supported: STATIC, PUSHOVER, CYCLIC_DISPLACEMENT, FREE-VIBRATION,
 SEISMIC, and external time-dependent loading (static or dynamic).

4. In the PUSHOVER case, displacement control is applied at the top node with small increments
 up to twice the ultimate displacement, recording base reaction, nodal displacement, element
 forces, and damage index at each step.

5. At every step, eigenvalue analysis (via `EIGENVALUE_ANALYSIS_FUN`) recomputes the structure's
 period (min/max) to track stiffness degradation.

6. The damage index is computed using a ductility-based function (`DAMAGE_INDEX_FUN`)
 comparing current displacement to yield and ultimate displacements.

7. Results are visualized: base shear vs. displacement, period evolution, element forces,
 nodal displacements, and damage index.

8. A `PLOT_1D_SPRING` routine optionally renders a virtual 1D spring representation of the
 deformed structure.

9. The pushover curve is then fitted with multilinear idealizations bilinear, trilinear,
 quadrilinear, and pentalinear using an area-preserving piecewise-linear fitting algorithm
 (`MULTILINEAR_CURVE`).

10. Fitting outputs include elastic stiffness, plastic stiffness, tangent stiffness, ductility ratio,
 and over-strength factor, plotted alongside the original pushover curve for comparison.
