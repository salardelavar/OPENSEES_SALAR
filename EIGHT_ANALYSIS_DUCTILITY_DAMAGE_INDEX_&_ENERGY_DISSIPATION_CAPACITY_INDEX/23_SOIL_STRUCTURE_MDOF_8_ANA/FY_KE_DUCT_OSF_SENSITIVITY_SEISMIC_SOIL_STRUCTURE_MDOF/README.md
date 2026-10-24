# SENSITIVITY ANALYSIS OF SOIL'S DAMPING RATIO, YIELD STRENGTH, ELASTIC STIFFNESS, DUCTILITY RATIO AND OVER-STRENGTH FACTOR ON THE SEISMIC PERFORMANCE OF SOIL–STRUCTURE MDOF SYSTEMS USING OPENSEES: CYCLIC DISPLACEMENT AND SEISMIC TIME-HISTORY ANALYSES WITH EVALUATION OF PERIOD ELONGATION,ENERGY DISSIPATION CAPACITY INDEX, EQUIVALENT VISCOUS DAMPING RATIO, DUCTILITY DAMAGE INDEX AND RECONSTRUCTION COST 

# تحلیل حساسیت نسبت میرایی خاک، مقاومت تسلیم، سختی الاستیک، نسبت شکل‌پذیری و ضریب مقاومت افزون بر عملکرد لرزه‌ای سیستم‌های چنددرجه‌آزادی خاک–سازه با استفاده از OpenSees: تحلیل‌های جابجایی چرخه‌ای و تاریخچه زمانی لرزه‌ای همراه با ارزیابی طولانی‌شدن دوره تناوب، شاخص ظرفیت اتلاف انرژی، نسبت میرایی ویسکوز معادل، شاخص خسارت شکل‌پذیری و هزینه بازسازی

This OpenSeesPy-based Python framework performs a comprehensive multi-parameter
 sensitivity analysis of a two-degree-of-freedom nonlinear soil–structure system subjected
 to both cyclic displacement protocols and selected ground-motion time histories.
 By systematically varying soil yield strength, elastic stiffness, ductility ratio,
 over-strength factor, and viscous damping ratios of both soil and structure across a
 discrete design space, the code assembles a high-fidelity hysteretic model using multi-linear
 OpenSees hysteretic materials and parallel viscous dampers, extracts instantaneous eigenvalue
 periods at every analysis step to quantify period elongation, and computes a suite of performance
 metrics—including peak displacement, velocity and acceleration demands, energy-dissipation capacity 
 index, equivalent viscous damping ratio from hysteresis loops, ductility-based damage index, and 
 reconstruction cost. Extensive post-processing then generates correlation scatter plots,
 three-dimensional response surfaces, Pearson sensitivity heatmaps, Random-Forest safety
 classification, and both main-effect and interaction ANOVA tables, thereby ranking the relative
 influence of each soil and structural parameter on seismic performance and providing a rigorous
 quantitative basis for performance-based design and risk assessment of soil–structure systems.
 
![alt text](https://github.com/salardelavar/OPENSEES_SALAR/blob/main/EIGHT_ANALYSIS_DUCTILITY_DAMAGE_INDEX_%26_ENERGY_DISSIPATION_CAPACITY_INDEX/23_SOIL_STRUCTURE_MDOF_8_ANA/FY_KE_DUCT_OSF_SENSITIVITY_SEISMIC_SOIL_STRUCTURE_MDOF/COVER.png)

# 200 SEISMIC GROUND MOTIONS ANALYSIS
![alt text](SEISMIC_GROUND_MOTION_ANALYSIS.png)

THIS PYTHON SCRIPT IS WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)

Please note that the content may not be entirely free of errors or inaccuracies.





# FY–KE–DUCT–OSF SENSITIVITY ANALYSIS OF A SEISMIC SOIL–STRUCTURE INTERACTION MDOF SYSTEM

**Nonlinear Seismic Analysis | Soil–Structure Interaction (SSI) | Parameter Sensitivity | OpenSeesPy | Python**

## 1. Overview

This project focuses on the sensitivity analysis of a multi-degree-of-freedom (MDOF) soil–structure interaction system subjected to seismic excitation. The analysis investigates how variations in selected mechanical and structural parameters influence the nonlinear dynamic response of the coupled soil–structure system.

The investigated parameters are:

- **FY — Yield Strength:** Represents the yield-strength parameter assigned to the relevant nonlinear component or structural model.
- **KE — Elastic Stiffness:** Represents the initial elastic stiffness of the modeled component or subsystem.
- **DUCT — Ductility Capacity:** Represents the ductility parameter used to characterize deformation capacity.
- **OSF — Overstrength Factor:** Represents the overstrength parameter used in the structural model.

By systematically varying these parameters and evaluating the resulting seismic responses, the framework can be used to identify influential modeling variables, investigate nonlinear response mechanisms, and examine the sensitivity of structural performance to parameter uncertainty.

The project is part of the `OPENSEES_SALAR` repository, which contains Python-based computational tools and structural engineering analysis examples.

## 2. Project Objectives

The primary objectives of this project are to:

1. Investigate the influence of yield strength on the seismic response of a soil–structure interaction system.
2. Evaluate the sensitivity of nonlinear response to changes in elastic stiffness.
3. Examine how ductility capacity affects deformation demand and structural response.
4. Investigate the role of the overstrength factor in the modeled system.
5. Compare the response obtained from different combinations of the investigated parameters.
6. Identify parameter combinations associated with comparatively high or low response demands.
7. Establish a computational framework for examining parameter sensitivity in nonlinear seismic analysis.
8. Support the interpretation of coupled soil–structure behavior under earthquake excitation.

## 3. Theoretical Background

### 3.1. Soil–Structure Interaction

Soil–structure interaction describes the mutual influence between a structure and its supporting soil. Unlike a fixed-base structural model, an SSI model can account for the influence of soil flexibility and the mechanical characteristics of the foundation–soil system on structural response.

Depending on the adopted modeling assumptions, SSI can influence:

- Natural periods and vibration modes.
- Effective system stiffness.
- Seismic displacement demands.
- Acceleration and force distributions.
- Nonlinear deformation demands.
- Energy transfer and dissipation within the coupled system.

The magnitude and direction of these effects depend on the soil properties, structural characteristics, foundation representation, excitation, and nonlinear constitutive assumptions.

### 3.2. Multi-Degree-of-Freedom Formulation

The dynamic equilibrium equation of a general MDOF system can be expressed as

\[
\mathbf{M}\ddot{\mathbf{u}}(t)
+\mathbf{C}\dot{\mathbf{u}}(t)
+\mathbf{f}_{\mathrm{int}}(\mathbf{u},\mathbf{z})
=
-\mathbf{M}\mathbf{r}a_g(t)
\]

where:

- \(\mathbf{M}\) is the mass matrix.
- \(\mathbf{C}\) is the damping matrix.
- \(\mathbf{u}(t)\) is the displacement vector relative to the ground.
- \(\dot{\mathbf{u}}(t)\) is the relative velocity vector.
- \(\ddot{\mathbf{u}}(t)\) is the relative acceleration vector.
- \(\mathbf{f}_{\mathrm{int}}\) is the internal restoring-force vector.
- \(\mathbf{z}\) represents internal variables associated with nonlinear material behavior.
- \(\mathbf{r}\) is the influence vector for ground excitation.
- \(a_g(t)\) is the imposed ground acceleration.

For nonlinear systems, the restoring-force vector depends on the current deformation and the relevant constitutive history variables.

The exact formulation depends on the implemented structural model, soil representation, boundary conditions, and excitation method.

### 3.3. Yield Strength (FY)

Yield strength influences the onset of inelastic behavior in a nonlinear component.

Changes in yield strength can affect:

- The transition from elastic to inelastic response.
- The magnitude and distribution of plastic deformation.
- Residual deformation.
- Force redistribution.
- Hysteretic energy dissipation.

The effect of FY cannot be interpreted independently of the adopted stiffness, deformation capacity, and constitutive law. In particular, changing FY while keeping other parameters fixed may also change the yield displacement or yield acceleration.

### 3.4. Elastic Stiffness (KE)

Elastic stiffness controls the initial force–deformation relationship of the modeled component.

For a linear elastic range,

\[
F=K_Eu
\]

where \(F\) is the restoring force and \(u\) is the corresponding deformation.

Changes in KE can influence the system's natural periods, force distribution, and dynamic amplification. In a coupled soil–structure system, changes in component stiffness can also modify the interaction between the structural and soil subsystems.

The precise influence depends on whether KE represents structural stiffness, foundation stiffness, soil stiffness, or a specific nonlinear component's initial tangent stiffness.

### 3.5. Ductility Capacity (DUCT)

Ductility describes the capacity of a component or structural system to sustain deformation beyond the onset of yielding.

A conventional displacement ductility ratio is

\[
\mu=\frac{|u_{\max}|}{u_y}
\]

where \(u_{\max}\) is the maximum deformation and \(u_y\) is the yield deformation.

This response measure should be distinguished from a prescribed ductility-capacity parameter. In this project, DUCT denotes the parameter used by the implemented model; its exact definition must be taken from the corresponding Python script.

Increasing the prescribed ductility capacity does not necessarily reduce displacement demand. Its influence depends on the constitutive model, post-yield behavior, strength, stiffness, and earthquake record.

### 3.6. Overstrength Factor (OSF)

The overstrength factor characterizes the relationship between a reference strength level and a higher strength level defined by the adopted modeling convention.

A general representation is

\[
\Omega=\frac{F_{\mathrm{reference}}^{*}}
{F_{\mathrm{design}}}
\]

where the numerator and denominator must be defined consistently with the implementation.

In this project, OSF refers to the model's overstrength parameter. Its exact mathematical relationship with FY and other strength parameters should be verified from the code before interpreting its effects.

## 4. Sensitivity Analysis Framework

The analysis considers four input parameters:

| Parameter | Description | Primary area of investigation |
|---|---|---|
| FY | Yield strength | Strength-dependent nonlinear response |
| KE | Elastic stiffness | Initial stiffness and dynamic characteristics |
| DUCT | Ductility capacity | Deformation capacity and inelastic behavior |
| OSF | Overstrength factor | Strength relationships and nonlinear response |

A parameter sensitivity study typically involves the following procedure:

1. Define the reference model and baseline parameter values.
2. Establish the parameter ranges or discrete values to be investigated.
3. Generate parameter combinations according to the adopted sampling strategy.
4. Construct or update the soil–structure interaction model for each combination.
5. Apply the selected earthquake ground-motion record or records.
6. Perform nonlinear time-history analysis.
7. Extract the response quantities required for comparison.
8. Compare the response across the parameter combinations.
9. Quantify sensitivity using an appropriate statistical or numerical method.
10. Identify important trends, interactions, and potential limitations.

The actual sampling strategy, number of simulations, and response quantities should be documented from the implemented script.

## 5. Seismic Input and Nonlinear Dynamic Analysis

Earthquake ground-motion analysis provides a means of investigating the dynamic behavior of the modeled system under time-varying excitation.

For a reproducible study, the following information should be documented where applicable:

- Ground-motion record identification and source.
- Sampling interval and record duration.
- Acceleration units and preprocessing.
- Ground-motion scaling procedure.
- Damping model and damping parameters.
- Numerical integration algorithm.
- Nonlinear solution algorithm and convergence tolerances.
- Boundary conditions and support constraints.
- Material models and component force–deformation relationships.
- Treatment of soil and foundation flexibility.

If multiple ground motions are used, the record-to-record variability should be distinguished from sensitivity caused by the four model parameters.

## 6. Response Quantities and Performance Evaluation

The sensitivity analysis can be evaluated using response quantities that are actually calculated and stored by the implementation.

Potential quantities include:

| Response quantity | Engineering interpretation |
|---|---|
| Maximum displacement | Peak deformation demand |
| Interstory drift ratio | Relative deformation demand between adjacent floors |
| Peak acceleration | Maximum absolute or relative acceleration, according to the output definition |
| Base shear | Global lateral-force demand |
| Residual displacement | Permanent deformation after excitation |
| Ductility demand | Deformation relative to yield deformation |
| Hysteretic energy | Energy dissipated through cyclic inelastic behavior |
| Natural period | Dynamic characteristic of the modeled system |
| Damage index | A model-specific measure of accumulated or deformation-related damage |

These quantities are candidate evaluation metrics, not a claim that every metric is calculated by the current script.

A meaningful comparison requires consistent units, response definitions, parameter bounds, and analysis settings.

### 6.1. Damage and Energy-Based Measures

If the implementation evaluates ductility-based damage or energy-dissipation indices, the corresponding definitions and normalization procedures should be documented explicitly.

For example, a displacement-based ductility measure can be defined as

\[
\mu_d=\frac{|u_{\max}|}{u_y}
\]

where \(u_y\) must be defined consistently for the nonlinear component or system under investigation.

Energy-based measures require a clear definition of the energy quantity being evaluated. The area enclosed by a force–displacement hysteresis loop, the cumulative work of the restoring force, and the total energy dissipated by the complete damped system are different quantities and should not be treated as interchangeable.

Any damage index should be interpreted according to its mathematical definition and calibration. A numerical damage index is not automatically equivalent to a probability of failure or a code-defined performance state.

## 7. Interpretation of Parameter Effects

The four parameters may interact rather than produce independent effects.

Examples of questions that can be investigated include:

- Does increasing FY delay yielding under a given earthquake record?
- How does changing KE alter the dynamic characteristics and displacement response?
- How does the prescribed ductility capacity affect the modeled nonlinear force–deformation relationship?
- Does changing OSF modify the strength available to resist seismic demand?
- Are parameter effects consistent across different ground motions?
- Which parameters have the greatest influence on the selected response quantity?
- Do combinations of parameters produce behavior that cannot be explained by one-factor-at-a-time comparisons?

These are research questions rather than predetermined conclusions. The direction and magnitude of each effect must be established from the computed results.

## 8. Recommended Statistical Analysis

Where the simulation dataset supports it, several approaches can be used to quantify sensitivity.

### 8.1. One-Factor-at-a-Time Analysis

Vary one parameter while holding the others at their reference values.

This approach is useful for preliminary investigation but does not capture interactions between parameters.

### 8.2. Factorial Analysis

Use combinations of parameter levels to investigate main effects and interaction effects.

For example, an FY–KE interaction can indicate whether the influence of yield strength depends on the elastic stiffness level.

### 8.3. Analysis of Variance (ANOVA)

ANOVA can be applied to investigate whether differences in a selected response are associated with categorical parameter levels or appropriately binned continuous parameters.

The model specification should account for relevant interactions and, where appropriate, earthquake-record variability. Statistical significance alone does not establish engineering significance.

### 8.4. Regression-Based Sensitivity

A regression model can approximate the relationship between the input parameters and a response quantity.

Depending on the dataset, candidate models include linear regression, polynomial regression, and other response-surface methods.

Model quality should be evaluated using suitable diagnostics and validation procedures.

### 8.5. Global Sensitivity Analysis

Where parameter distributions and sampling are adequately defined, global methods such as variance-based sensitivity analysis can evaluate the contribution of input uncertainty to response variability.

The selected method should be consistent with the sample design, computational budget, and treatment of correlated inputs.

## 9. Computational Implementation

The project is organized within the following repository directory:

```text
OPENSEES_SALAR/
└── EIGHT_ANALYSIS_DUCTILITY_DAMAGE_INDEX_&_ENERGY_DISSIPATION_CAPACITY_INDEX/
    └── 23_SOIL_STRUCTURE_MDOF_8_ANA/
        └── FY_KE_DUCT_OSF_SENSITIVITY_SEISMIC_SOIL_STRUCTURE_MDOF/
```

The directory name identifies the project as a parameter-sensitivity study involving a seismic soil–structure interaction MDOF model.

The exact Python entry point, supporting modules, output filenames, and dependency versions should be listed here after checking the repository contents.

### Recommended Software

- **Python** for model orchestration, numerical processing, and data analysis.
- **OpenSeesPy**, if used by the implementation, for finite-element-based structural analysis.
- **NumPy** for numerical arrays and mathematical operations.
- **Pandas** for tabular data management.
- **Matplotlib** for response and sensitivity plots.
- **SciPy** for additional numerical and statistical methods, where required.

Only dependencies imported by the actual scripts need to be installed.

## 10. Reproducibility and Numerical Verification

For reproducible computational experiments, the following practices are recommended:

1. Record the baseline values and units of FY, KE, DUCT, and OSF.
2. Document the parameter ranges and sampling method.
3. Keep earthquake records and preprocessing procedures identifiable.
4. Store the analysis configuration alongside the response data.
5. Check nonlinear convergence and retain information about failed analyses.
6. Verify that all parameter combinations satisfy the intended physical and constitutive constraints.
7. Confirm that comparisons use consistent response definitions.
8. Record software versions and random seeds when stochastic sampling is used.
9. Distinguish numerical failures from genuine large-response or collapse outcomes.
10. Validate selected cases against independent calculations or benchmark models where possible.

A sensitivity study is only as reliable as its underlying model, numerical settings, and interpretation of the resulting data.

## 11. Limitations and Scope

The conclusions of this analysis are conditional on the implemented model and its assumptions.

Important considerations include:

- The adopted representation of soil–structure interaction.
- The calibration and limitations of the nonlinear constitutive models.
- The selected ground-motion records and their scaling.
- The ranges and physical relationships of the investigated parameters.
- The numerical integration and convergence settings.
- The definition of each response and damage metric.
- The statistical treatment of record-to-record variability.

Parameter sensitivity should not be interpreted as proof of causality outside the adopted computational model. Likewise, conclusions obtained for one structural configuration or ground-motion dataset should not automatically be generalized to all soil–structure systems.

## 12. Intended Applications

This computational framework can support:

- Educational demonstrations of nonlinear seismic analysis.
- Investigation of soil–structure interaction effects.
- Parametric studies of strength and stiffness.
- Exploration of deformation capacity and overstrength assumptions.
- Comparison of nonlinear response across model configurations.
- Development and testing of Python-based sensitivity-analysis tools.
- Preliminary assessment of influential parameters for further investigation.

It is a computational analysis example, not a substitute for project-specific structural assessment, geotechnical evaluation, or code-compliant engineering design.

## 13. Author and Repository

**Author:** Salar Delavar Ghashghaei (Qashqai)

**GitHub:** [salardelavar](https://github.com/salardelavar)

**Repository:** [OPENSEES_SALAR](https://github.com/salardelavar/OPENSEES_SALAR)

This project is intended to support the development, study, and educational use of computational methods in structural engineering, nonlinear analysis, and earthquake engineering.

## 14. Disclaimer

This repository provides computational examples for educational and research purposes. The scripts may contain errors, numerical limitations, or modeling assumptions that require independent verification.

Users should review the implementation, validate the numerical results, and confirm the suitability of the adopted models before applying the methods to engineering decisions. The author does not guarantee that the scripts are suitable for a specific design, assessment, or safety-critical application.

---

**Keywords:** OpenSeesPy, Python, soil–structure interaction, SSI, MDOF system, nonlinear time-history analysis, seismic response, sensitivity analysis, yield strength, elastic stiffness, ductility capacity, overstrength factor, earthquake engineering, structural dynamics.
