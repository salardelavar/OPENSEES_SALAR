# COMPREHENSIVE NONLINEAR SEISMIC ASSESSMENT OF A COMPOSITE SECTION FRAME : AN OPENSEES FRAMEWORK FOR STATIC PUSHOVER, CYCLIC DEGRADATION, STATIC TIME-HISTORY AND DYNAMIC TIME-HISTORY ANALYSIS

ارزیابی جامع غیرخطی سازه با المان مرکب فلزی بتنی : چارچوبی مبتنی بر اوپنسیس در پایتون برای تحلیل پوش‌آور استاتیکی، تحلیل چرخه‌ای همراه با کاهش سختی، تحلیل تاریخچه زمانی استاتیکی و دینامیکی، با در نظر گرفتن شاخص‌های آسیب شکل‌پذیری در المان و سازه و ارزیابی شاخص ظرفیت انرژی اتلاف‌شده


![alt text](https://github.com/salardelavar/OPENSEES_SALAR/blob/main/EIGHT_ANALYSIS_DUCTILITY_DAMAGE_INDEX_%26_ENERGY_DISSIPATION_CAPACITY_INDEX/10_COMPOSITE_FRAME_8_ANA/+10_COMPOSITE_FRAME_8_ANA/COVER.png) 

# ADVANCED NONLINEAR STRUCTURAL ANALYSIS OF COMPOSITE FRAMES USING OPENSEESPY

## 1. Overview

This directory is part of the **EIGHT ANALYSIS DUCTILITY DAMAGE INDEX & ENERGY DISSIPATION CAPACITY INDEX** project in the [OPENSEES_SALAR repository](https://github.com/salardelavar/OPENSEES_SALAR).

The objective of this example is to provide a computational framework for investigating the structural response of composite frame systems using Python and OpenSeesPy. Composite structural systems combine different construction materials, such as structural steel and reinforced concrete, to achieve a desired balance between strength, stiffness, ductility, and deformation capacity. Their nonlinear response depends on material constitutive behavior, member interaction, connection characteristics, boundary conditions, gravity loading, and the applied loading history.

The framework is intended to support the numerical investigation of structural response, nonlinear behavior, deformation capacity, and energy-related performance indicators within a unified structural-analysis workflow.

The implementation should be interpreted according to the constitutive models, element formulations, loading protocols, and response quantities actually defined in each Python script. The presence of an analysis category in this documentation does not imply that every possible analysis procedure is implemented in every model.

## 2. Engineering Objectives

The principal objectives of the composite-frame example are to:

- Develop a computational model of a composite structural system using OpenSeesPy.
- Investigate the relationship between structural stiffness, strength, and deformation demand.
- Examine the transition from elastic behavior to nonlinear response when supported by the selected material and element formulations.
- Evaluate the structural capacity under appropriate static or dynamic loading protocols.
- Investigate displacement demand, force response, and deformation-related performance indicators.
- Quantify hysteretic energy dissipation when the implemented loading protocol and recorded response permit a reliable calculation.
- Establish a foundation for further parametric studies, uncertainty analysis, and nonlinear seismic assessment.

## 3. Computational Framework

The analysis is based on the finite-element modeling capabilities of OpenSees and its Python interface, OpenSeesPy.

A nonlinear structural model generally consists of the following components:

| Component | Engineering purpose |
|---|---|
| Model geometry | Defines the structural configuration, nodes, connectivity, and degrees of freedom. |
| Material constitutive laws | Describe the elastic and inelastic behavior of the modeled materials. |
| Element formulations | Represent the stiffness and force-deformation behavior of structural members. |
| Boundary conditions | Define supports, restraints, and other kinematic constraints. |
| Mass representation | Establishes the inertial properties required for dynamic analyses. |
| Gravity loading | Establishes the initial force state when gravity analysis is included. |
| Analysis configuration | Specifies the solution algorithm, convergence tests, constraints, and integration procedure. |
| Response recorders | Store the quantities needed for post-processing and performance assessment. |

The accuracy of the resulting response depends on the consistency of these modeling assumptions and the suitability of the selected numerical formulations.

## 4. Nonlinear Structural Response

Composite frames can exhibit nonlinear behavior through material yielding, cracking, stiffness degradation, local deformation mechanisms, and connection response, depending on the physical behavior represented by the model.

A nonlinear analysis can be used to examine:

- The relationship between applied load and structural displacement.
- The development of inelastic deformation.
- Changes in effective structural stiffness.
- The distribution of internal forces and deformation demands.
- The influence of loading history on residual deformation and hysteretic response.
- The evolution of energy dissipation during cyclic or dynamic loading.

These mechanisms should only be attributed to the numerical model when they are explicitly represented by the selected material models, elements, and connection idealizations.

## 5. Analysis Procedures

The broader project is organized around multiple structural-analysis procedures. The following descriptions identify the principal analytical categories relevant to a composite-frame model; the actual availability of each procedure must be confirmed against the individual scripts.

### 5.1. Structural Period and Modal Properties

Eigenvalue analysis can be used to obtain natural frequencies, vibration periods, and mode shapes for a model with appropriate mass and stiffness definitions.

For a linearized system, the eigenvalue problem is

\[
\left[\mathbf K-\omega_i^2\mathbf M\right]\boldsymbol\phi_i=\mathbf 0
\]

where:

- \(\mathbf K\) is the relevant stiffness matrix.
- \(\mathbf M\) is the mass matrix.
- \(\omega_i\) is the natural circular frequency of mode \(i\).
- \(\boldsymbol\phi_i\) is the corresponding mode shape.

The natural period is calculated as

\[
T_i=\frac{2\pi}{\omega_i}
\]

These properties depend on the modeled mass distribution, boundary conditions, and stiffness assumptions. Periods from an eigenvalue analysis should not automatically be interpreted as the effective periods of a structure undergoing substantial nonlinear deformation.

### 5.2. Gravity Load Analysis

Where implemented, gravity analysis establishes the initial equilibrium state before subsequent loading.

The resulting internal forces and deformations provide the initial conditions for later analyses. Gravity loads should be applied consistently with the intended physical model, and the converged gravity state should be retained when subsequent analysis procedures require it.

### 5.3. Nonlinear Static Pushover Analysis

Pushover analysis evaluates the relationship between an imposed lateral deformation pattern and the corresponding structural restoring force.

A typical capacity curve is expressed as

\[
V_b=f(\Delta_r)
\]

where \(V_b\) denotes base shear and \(\Delta_r\) denotes the selected control displacement.

The resulting curve can be used to investigate initial lateral stiffness, strength development, nonlinear deformation, post-yield response, and the evolution of the modeled structural mechanism.

Interpretation requires attention to the imposed load pattern, control node, displacement increment, convergence criteria, and termination condition. A pushover curve alone does not establish the probability of damage or collapse under earthquake loading.

### 5.4. Cyclic Loading and Hysteretic Response

When cyclic loading is implemented, the force-displacement history can be used to investigate hysteresis, unloading and reloading behavior, residual deformation, and changes in response across repeated cycles.

The enclosed area of a complete force-displacement loop represents the energy dissipated during that cycle, subject to the assumptions of the modeled force and deformation quantities.

For a closed cycle,

\[
E_{\mathrm{diss}}=\oint F\,du
\]

where \(F\) is the relevant resisting force and \(u\) is the corresponding deformation.

The integration procedure must preserve the chronological order of the recorded response. The accumulated signed work, loop-by-loop dissipation, and geometrical area of a plotted curve are not interchangeable quantities.

### 5.5. Time-History Analysis

Where dynamic time-history analysis is provided, the structural equations of motion can be expressed as

\[
\mathbf M\ddot{\mathbf u}
+\mathbf C\dot{\mathbf u}
+\mathbf f_{\mathrm{int}}(\mathbf u,\mathbf z)
=-\mathbf M\mathbf r\,a_g(t)
\]

where:

- \(\mathbf M\) is the mass matrix.
- \(\mathbf C\) is the damping matrix.
- \(\mathbf u\), \(\dot{\mathbf u}\), and \(\ddot{\mathbf u}\) are the displacement, velocity, and acceleration vectors.
- \(\mathbf f_{\mathrm{int}}\) represents the internal restoring forces, potentially depending on deformation history through internal variables \(\mathbf z\).
- \(\mathbf r\) is the influence vector for ground excitation.
- \(a_g(t)\) is the ground acceleration time history.

A dynamic analysis can provide time-dependent response histories, including displacement, velocity, acceleration, and restoring-force quantities, when these are recorded by the model.

The reliability of seismic response estimates depends on ground-motion scaling, time-step selection, numerical convergence, damping assumptions, mass representation, and the treatment of nonlinear material behavior.

## 6. Ductility and Damage Assessment

Ductility describes the capacity of a structural component or system to undergo deformation beyond its elastic range. A common displacement-based definition is

\[
\mu=\frac{\Delta_{\max}}{\Delta_y}
\]

where \(\Delta_{\max}\) is the maximum displacement demand and \(\Delta_y\) is a consistently defined yield displacement.

The definition must be appropriate to the response quantity under consideration. Global displacement ductility, member deformation ductility, and material strain ductility are different measures and should not be used interchangeably.

A ductility-based damage index can provide a normalized measure of deformation demand relative to a selected reference capacity. However, its physical interpretation depends on the exact mathematical formulation, yield definition, deformation limits, and calibration procedure.

If the project computes a Ductility Damage Index (DDI), the README and accompanying scripts should explicitly identify:

1. The displacement or deformation quantity used.
2. The method used to determine yield deformation.
3. The normalization or limiting value.
4. Whether the index is calculated at the member or structural level.
5. Whether cumulative cyclic effects are included.

A normalized index is not automatically a damage probability or a code-defined performance level. Such interpretations require a suitable physical or probabilistic basis.

## 7. Energy Dissipation Capacity

Energy-based assessment provides a complementary perspective to displacement-based ductility measures.

The energy dissipated during a loading history can be obtained from the work associated with the relevant force-deformation response. For a discrete response history, incremental work can be approximated by

\[
\Delta W_i\approx
\frac{F_i+F_{i-1}}{2}
(u_i-u_{i-1})
\]

and accumulated work can be evaluated as

\[
W=\sum_{i=1}^{n}\Delta W_i
\]

This trapezoidal integration represents work along the recorded path. Its interpretation as dissipated energy requires the appropriate treatment of conservative elastic contributions, unloading, sign conventions, and the physical meaning of the recorded force.

If an Energy Dissipation Capacity Index (EDCI) is implemented, its numerator, denominator, reference capacity, and integration limits must be documented explicitly. In particular, the energy demand from a time-history analysis and the energy capacity from a cyclic or pushover-based procedure must be defined consistently before their ratio is interpreted.

An energy ratio is a model-dependent performance indicator rather than a universal measure of structural safety.

## 8. Equivalent Viscous Damping

Equivalent viscous damping can be used to describe the energy dissipation of a nonlinear system in terms of an equivalent linearized representation, provided that the selected formulation is appropriate to the response being analyzed.

For a stable, fully developed hysteresis loop, a commonly used expression is

\[
\xi_{\mathrm{eq}}=
\frac{E_D}{4\pi E_S}
\]

where:

- \(E_D\) is the energy dissipated during the cycle.
- \(E_S\) is the maximum recoverable strain energy associated with the adopted convention.

The exact definition of \(E_S\), the treatment of residual deformation, and the loop-selection procedure must follow the method used by the corresponding script.

Equivalent viscous damping is an interpretation of energy dissipation, not necessarily the same as the viscous damping ratio assigned to the OpenSees model.

## 9. Parametric Studies and Uncertainty Analysis

The composite-frame model can serve as a basis for systematic investigations of parameters that affect stiffness, strength, deformation demand, and energy dissipation.

Potential research extensions include:

- Material strength and elastic modulus.
- Member stiffness and section properties.
- Connection stiffness and rotational restraint.
- Mass distribution and structural period.
- Damping assumptions.
- Ground-motion characteristics and intensity.
- Yield deformation and ductility capacity.
- Model parameters associated with nonlinear degradation.

For a defensible parametric study, vary selected parameters while maintaining a clearly defined baseline model. Record the actual values used in each analysis and distinguish input uncertainty from numerical variability.

Where multiple parameters are varied simultaneously, their effects may be confounded. Sensitivity analysis, factorial designs, regression, or variance-based methods can help separate individual and interaction effects when the available data and sampling strategy support those methods.

## 10. Numerical Reliability and Verification

Nonlinear finite-element results require numerical checks before engineering interpretation.

Recommended checks include:

1. **Equilibrium:** Verify applied loads, support reactions, and residual equilibrium.
2. **Convergence:** Examine force and displacement convergence criteria and document failed steps.
3. **Time-step sensitivity:** For dynamic analyses, compare results using appropriate integration time steps.
4. **Energy consistency:** Check external work, kinetic energy, recoverable strain energy, and dissipative contributions where these quantities are available.
5. **Model sensitivity:** Assess the influence of constitutive parameters, boundary conditions, and element formulations.
6. **Independent comparison:** Compare selected response quantities with analytical solutions, benchmark models, experimental evidence, or independently implemented calculations where suitable reference data exist.
7. **Reproducibility:** Record software versions, input parameters, loading definitions, and analysis settings.

Successful numerical convergence does not by itself demonstrate physical accuracy. Verification of the numerical implementation and validation against independent evidence are separate tasks.

## 11. Software Requirements

The examples are intended for a Python environment with OpenSeesPy and the numerical or visualization libraries required by the individual scripts.

Typical dependencies for structural-analysis workflows include:

- Python
- OpenSeesPy
- NumPy
- Matplotlib
- Pandas
- SciPy, when required by the implementation
- OpenPyXL, when Excel output is generated

The actual requirements should be confirmed from the imports used by each Python file. Install compatible package versions in a dedicated virtual environment to improve reproducibility.

## 12. Running the Examples

1. Download or clone the [OPENSEES_SALAR repository](https://github.com/salardelavar/OPENSEES_SALAR).
2. Navigate to `EIGHT_ANALYSIS_DUCTILITY_DAMAGE_INDEX_&_ENERGY_DISSIPATION_CAPACITY_INDEX/10_COMPOSITE_FRAME_8_ANA/+10_COMPOSITE_FRAME_8_ANA`.
3. Review the Python script and identify its input parameters, material definitions, element formulations, and loading protocol.
4. Install the dependencies required by that script.
5. Execute the selected script from the configured Python environment.
6. Inspect the generated output files and plots, if provided.
7. Check convergence messages and verify that the recorded response quantities correspond to the intended engineering definitions.

Input parameters should be reviewed before changing the structural geometry, material properties, loading history, or numerical integration settings.

## 13. Scope and Limitations

This repository provides educational and computational examples for structural engineering research and development. Results depend on the assumptions and limitations of the implemented model.

The presence of nonlinear constitutive laws does not guarantee that all relevant composite-member mechanisms are represented. Interface slip, local buckling, composite action, connection fracture, residual stresses, geometric imperfections, and other phenomena require explicit modeling when they are relevant to the problem.

The scripts should not be treated as independently validated design software or as a substitute for applicable design standards, qualified engineering judgment, and project-specific verification.

## 14. Further Research Directions

Potential extensions to this composite-frame framework include:

- Incremental Dynamic Analysis (IDA) using a suite of ground motions.
- Fragility analysis based on explicitly defined engineering demand parameters and damage-state thresholds.
- Probabilistic modeling of material and connection uncertainties.
- Sensitivity analysis of structural strength, stiffness, ductility, and energy-related indicators.
- Comparison of global response with local member deformation demands.
- Verification against experimental composite-frame results or established numerical benchmarks.
- Performance-based assessment using clearly defined acceptance criteria.

These extensions should be implemented and verified individually before their results are presented as findings of the repository.

## 15. Disclaimer

The Python scripts are provided free of charge for educational, research, and computational development purposes. They may contain errors, simplifications, or limitations in numerical implementation and physical modeling.

Users are responsible for checking the model assumptions, validating the implementation, reviewing the output, and determining whether the methods are suitable for their intended application. No structural design decision should rely solely on these scripts without independent engineering verification.

## 16. Repository

**Author:** Salar Delavar Qashqai

**Main repository:** [OPENSEES_SALAR](https://github.com/salardelavar/OPENSEES_SALAR)

**Parent project:** [EIGHT ANALYSIS DUCTILITY DAMAGE INDEX & ENERGY DISSIPATION CAPACITY INDEX](https://github.com/salardelavar/OPENSEES_SALAR/tree/main/EIGHT_ANALYSIS_DUCTILITY_DAMAGE_INDEX_%26_ENERGY_DISSIPATION_CAPACITY_INDEX)

The project is intended to facilitate transparent, reproducible, and code-based exploration of structural and earthquake engineering concepts through OpenSeesPy.
