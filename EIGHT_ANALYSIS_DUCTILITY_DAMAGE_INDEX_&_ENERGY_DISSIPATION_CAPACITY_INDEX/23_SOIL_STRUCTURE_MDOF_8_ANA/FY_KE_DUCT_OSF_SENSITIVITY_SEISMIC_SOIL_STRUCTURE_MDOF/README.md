![alt text](https://github.com/salardelavar/OPENSEES_SALAR/blob/main/EIGHT_ANALYSIS_DUCTILITY_DAMAGE_INDEX_%26_ENERGY_DISSIPATION_CAPACITY_INDEX/23_SOIL_STRUCTURE_MDOF_8_ANA/FY_KE_DUCT_OSF_SENSITIVITY_SEISMIC_SOIL_STRUCTURE_MDOF/COVER.png)

# SENSITIVITY ANALYSIS OF SOIL'S DAMPING RATIO, YIELD STRENGTH, ELASTIC STIFFNESS, DUCTILITY RATIO AND OVER-STRENGTH FACTOR ON THE SEISMIC PERFORMANCE OF SOIL–STRUCTURE MDOF SYSTEMS USING OPENSEES: CYCLIC DISPLACEMENT AND SEISMIC TIME-HISTORY ANALYSES WITH EVALUATION OF PERIOD ELONGATION,ENERGY DISSIPATION CAPACITY INDEX, EQUIVALENT VISCOUS DAMPING RATIO, DUCTILITY DAMAGE INDEX AND RECONSTRUCTION COST 

# تحلیل حساسیت نسبت میرایی خاک، مقاومت تسلیم، سختی الاستیک، نسبت شکل‌پذیری و ضریب مقاومت افزون بر عملکرد لرزه‌ای سیستم‌های چنددرجه‌آزادی خاک–سازه با استفاده از OpenSees: تحلیل‌های جابجایی چرخه‌ای و تاریخچه زمانی لرزه‌ای همراه با ارزیابی طولانی‌شدن دوره تناوب، شاخص ظرفیت اتلاف انرژی، نسبت میرایی ویسکوز معادل، شاخص خسارت شکل‌پذیری و هزینه بازسازی


**A Computational Framework for Evaluating Soil Damping, Yield Strength, Elastic Stiffness, Ductility Capacity, Overstrength, Period Elongation, Energy Dissipation, Structural Damage, and Reconstruction Cost**

[![Python](https://img.shields.io/badge/Python-3.x-blue.svg)](https://www.python.org/)
[![OpenSeesPy](https://img.shields.io/badge/Analysis-OpenSeesPy-orange.svg)](https://openseespydoc.readthedocs.io/)
[![Engineering](https://img.shields.io/badge/Domain-Earthquake%20Engineering-darkred.svg)](https://github.com/salardelavar/OPENSEES_SALAR)

## 1. Project Overview

This project presents a Python-based computational framework for investigating the nonlinear seismic response of a multi-degree-of-freedom (MDOF) soil–structure interaction system using OpenSeesPy.

The framework examines how variations in soil and structural parameters influence the dynamic characteristics, nonlinear deformation demands, energy-dissipation behavior, and damage-related response of the coupled system. The analysis incorporates cyclic displacement loading and earthquake ground-motion time histories, with the seismic analysis component identified as involving **200 ground motions**.

The investigated parameters include soil damping ratio, yield strength, elastic stiffness, ductility capacity, and overstrength factor. Their effects are evaluated through numerical simulations and post-processing procedures designed to characterize changes in structural response and identify influential parameters.

The computational workflow combines nonlinear structural analysis with response extraction, damage and energy-related calculations, visualization, and statistical sensitivity assessment. Depending on the implemented analysis configuration, the resulting dataset can be used to investigate individual parameter effects, interactions between parameters, and variations in response across earthquake records.

The project is part of the [`OPENSEES_SALAR`](https://github.com/salardelavar/OPENSEES_SALAR) repository, which provides Python-based computational examples and tools for structural engineering, nonlinear analysis, and earthquake engineering.

![Soil–structure interaction MDOF model](https://github.com/salardelavar/OPENSEES_SALAR/blob/main/EIGHT_ANALYSIS_DUCTILITY_DAMAGE_INDEX_%26_ENERGY_DISSIPATION_CAPACITY_INDEX/23_SOIL_STRUCTURE_MDOF_8_ANA/FY_KE_DUCT_OSF_SENSITIVITY_SEISMIC_SOIL_STRUCTURE_MDOF/COVER.png)

## 2. Analysis Objectives

The primary objectives are to:

1. Investigate the influence of soil damping on the dynamic response of the coupled soil–structure system.
2. Evaluate the sensitivity of nonlinear response to variations in yield strength and elastic stiffness.
3. Examine the effects of ductility capacity and overstrength on deformation demands and nonlinear behavior.
4. Investigate changes in the system's dynamic characteristics, including period elongation.
5. Evaluate hysteretic energy dissipation and equivalent viscous damping using the implemented calculation procedures.
6. Quantify ductility-related damage indicators and reconstruction-cost measures.
7. Compare responses obtained under different parameter combinations and earthquake ground motions.
8. Apply statistical and machine-learning-based techniques, where implemented, to investigate parameter importance and response variability.
9. Develop a reproducible computational workflow for nonlinear seismic sensitivity analysis.

These objectives describe the scope of the computational investigation; the actual influence of each parameter must be established from the numerical results.

## 3. Investigated Parameters

| Parameter | Description | Engineering significance |
|---|---|---|
| Soil damping ratio | Viscous damping assigned to the soil subsystem | Influences energy dissipation and dynamic response |
| Structural damping ratio | Viscous damping assigned to the structural subsystem | Influences the dynamic response and energy dissipation |
| **FY** | Yield strength | Controls the onset of inelastic behavior and influences strength-dependent response |
| **KE** | Elastic stiffness | Influences initial force–deformation behavior and dynamic characteristics |
| **DUCT** | Ductility-capacity parameter | Defines a deformation-capacity parameter in the adopted nonlinear model |
| **OSF** | Overstrength factor | Represents a strength-related parameter according to the model's definition |

The precise mathematical relationships, units, parameter ranges, and dependencies must be taken from the corresponding Python implementation. In particular, changing yield strength or overstrength may also change other properties of the force–deformation relationship, depending on how the model is formulated.

## 4. Theoretical Background

### 4.1. Soil–Structure Interaction

Soil–structure interaction (SSI) describes the coupled behavior of a structure and its supporting soil. Soil flexibility and damping can modify the effective dynamic properties of the system compared with an idealized fixed-base structural model.

Potentially affected response quantities include:

- Natural periods and vibration characteristics.
- Displacement and acceleration demands.
- Internal forces and force redistribution.
- Nonlinear deformation and yielding.
- Energy transfer and dissipation.
- Residual deformation and damage-related indicators.

The magnitude and direction of these effects depend on the soil and structural properties, foundation representation, boundary conditions, constitutive models, and earthquake excitation.

### 4.2. Governing Dynamic Equilibrium

A general equation of motion for a nonlinear MDOF system subjected to ground acceleration is

\[
\mathbf{M}\ddot{\mathbf{u}}(t)
+\mathbf{C}\dot{\mathbf{u}}(t)
+\mathbf{f}_{\mathrm{int}}(\mathbf{u},\mathbf{z})
=-\mathbf{M}\mathbf{r}a_g(t)
\]

where:

- \(\mathbf{M}\) is the mass matrix.
- \(\mathbf{C}\) is the damping matrix.
- \(\mathbf{u}(t)\) is the displacement vector relative to the ground.
- \(\dot{\mathbf{u}}(t)\) and \(\ddot{\mathbf{u}}(t)\) are the relative velocity and acceleration vectors.
- \(\mathbf{f}_{\mathrm{int}}\) is the internal restoring-force vector.
- \(\mathbf{z}\) represents the internal variables associated with nonlinear material behavior.
- \(\mathbf{r}\) is the ground-motion influence vector.
- \(a_g(t)\) is the ground acceleration.

This equation provides the general theoretical framework for nonlinear seismic analysis. The specific formulation is determined by the implemented model, including its degrees of freedom, soil representation, boundary conditions, and damping formulation.

### 4.3. Yield Strength and Elastic Stiffness

For a component with an initial linear elastic range, the force–deformation relationship can be represented as

\[
F=K_Eu
\]

where \(K_E\) is the initial elastic stiffness and \(u\) is the corresponding deformation.

Yield strength determines the force level associated with the onset of yielding under the adopted constitutive model. Together, stiffness and strength influence the yield deformation, post-yield response, and distribution of inelastic deformation.

Consequently, parameter sensitivity should be interpreted in terms of the complete force–deformation relationship rather than treating strength and stiffness as entirely independent physical effects.

### 4.4. Ductility Capacity and Overstrength

A conventional displacement ductility demand is defined as

\[
\mu_d=\frac{|u_{\max}|}{u_y}
\]

where \(u_{\max}\) is the maximum displacement and \(u_y\) is the reference yield displacement.

Ductility demand is distinct from a prescribed ductility-capacity parameter. The relationship between the two depends on the nonlinear model and the definitions adopted in the code.

Similarly, the overstrength factor represents a ratio between strength levels defined by the adopted modeling convention. Its exact relationship with FY and other strength parameters must be verified from the implementation before interpreting its influence on seismic performance.

## 5. Analysis Methodology

The overall computational workflow consists of model definition, parameter variation, nonlinear analysis, response extraction, and sensitivity assessment.

1. **Model definition:** Establish the soil–structure system, material properties, masses, boundary conditions, and damping formulation.
2. **Parameter specification:** Define the investigated values or ranges of soil damping, FY, KE, DUCT, and OSF.
3. **Parameter combinations:** Generate the combinations required for the sensitivity study.
4. **Cyclic displacement analysis:** Evaluate the force–deformation response under the cyclic displacement protocol used in the model.
5. **Seismic time-history analysis:** Apply the selected earthquake ground motions to investigate dynamic response.
6. **Response extraction:** Obtain the required displacement, velocity, acceleration, force, and other response histories.
7. **Dynamic-characteristic evaluation:** Calculate period-related quantities using the implemented eigenvalue-analysis procedure.
8. **Damage and energy assessment:** Evaluate the implemented energy-dissipation, equivalent-damping, ductility-damage, and reconstruction-cost measures.
9. **Statistical post-processing:** Analyze relationships between model parameters and selected response quantities using the implemented statistical methods.
10. **Visualization and comparison:** Generate plots and summary results for interpreting parameter sensitivity.

The number of parameter combinations, analysis steps, ground-motion records successfully completed, and treatment of failed simulations should be verified from the output data when reporting final results.

## 6. Seismic Ground-Motion Analysis

The seismic component of this project is identified as a **200-ground-motion analysis**.

Using multiple ground motions allows the response to be examined under different excitation histories rather than relying on a single earthquake record. It can also help distinguish parameter-related response changes from record-to-record variability, provided that the statistical analysis treats the records appropriately.

For reproducibility, the following information should be documented alongside the analysis results:

- Ground-motion record identifiers and sources.
- Sampling intervals, record durations, and acceleration units.
- Record preprocessing and scaling procedures.
- Soil and structural damping definitions.
- Nonlinear material models and parameter values.
- Numerical integration and solution algorithms.
- Convergence criteria and handling of unsuccessful analyses.
- The number of records successfully analyzed for each parameter combination.

The use of 200 ground motions does not, by itself, establish that the records are statistically representative of every seismic hazard or site condition. Their suitability depends on record selection, scaling, and the objectives of the study.

![Seismic ground-motion analysis](SEISMIC_GROUND_MOTION_ANALYSIS.png)

## 7. Response Evaluation and Performance Indicators

The framework investigates several response and performance indicators. Their exact definitions should follow the relevant functions in the Python implementation.

| Indicator | Purpose |
|---|---|
| Displacement response | Characterizes deformation demand |
| Velocity response | Describes the system's motion during excitation |
| Acceleration response | Characterizes acceleration demand |
| Base reaction or restoring force | Supports evaluation of global force response, subject to the recorded quantity's definition |
| Period elongation | Examines changes in the calculated dynamic period |
| Energy-dissipation capacity index (EDCI) | Represents an implemented energy-related performance measure |
| Equivalent viscous damping ratio (EVDR) | Expresses hysteretic behavior through an equivalent damping measure |
| Ductility damage index (DDI) | Represents a deformation-related damage indicator |
| Reconstruction cost | Estimates a cost-related measure using the implemented cost function |

These indicators are not interchangeable. Their engineering interpretation depends on their mathematical definitions, normalization, and calibration.

### 7.1. Period Elongation

Period elongation refers to an increase in the calculated vibration period relative to a defined reference condition. In nonlinear analysis, changes in effective stiffness can influence the system's dynamic characteristics.

The interpretation of period elongation depends on how the period is calculated, whether the analysis uses instantaneous or reference eigenvalues, and how changes in the modeled properties are handled. Comparisons should therefore use a consistent period definition and clearly stated reference state.

### 7.2. Energy Dissipation and Equivalent Viscous Damping

Energy-based response measures require a clear distinction between different physical quantities.

The area enclosed by a force–displacement hysteresis loop represents energy dissipated over that cycle under the relevant assumptions. Cumulative restoring-force work, hysteretic energy, and total energy dissipated by the complete damped system are not necessarily equivalent.

Likewise, equivalent viscous damping is a model-dependent representation of energy dissipation. The calculation should identify the relevant hysteresis loop, energy definition, and normalization procedure.

### 7.3. Ductility Damage Index

A ductility-based indicator can be used to characterize deformation relative to a reference yield displacement or another specified deformation limit.

The interpretation of DDI depends on the exact implemented equation, its normalization, and whether it represents peak deformation, cumulative damage, or another response measure. It should not automatically be interpreted as a probability of failure or as a code-defined damage state.

### 7.4. Reconstruction Cost

The reconstruction-cost function provides a cost-related indicator based on the assumptions encoded in the implementation.

Its engineering use requires documentation of the cost equation, reference deformation levels, normalization, and any limits imposed on the calculated value. Unless calibrated against appropriate economic or empirical data, the result should be interpreted as a model-based estimate rather than a validated prediction of actual repair expenditure.

## 8. Statistical and Machine-Learning-Based Sensitivity Analysis

Statistical post-processing can help identify relationships between the investigated parameters and the selected response quantities.

The project description identifies the following analytical approaches.

### 8.1. Correlation Analysis

Correlation analysis can be used to investigate the direction and strength of statistical relationships between input parameters and response variables.

Correlation alone does not establish causation and may fail to capture nonlinear relationships or parameter interactions.

### 8.2. Three-Dimensional Response Surfaces

Response surfaces can illustrate how a selected response varies with combinations of two input parameters while other parameters or conditions are held fixed or otherwise accounted for.

Their interpretation depends on the available parameter combinations, the range of the plotted variables, and whether the displayed surface is based on actual simulations or an interpolated approximation.

### 8.3. Pearson Correlation Heatmaps

Pearson correlation heatmaps provide a visual summary of linear associations between variables. They can help identify preliminary trends but should be supplemented by other methods when relationships are nonlinear or influenced by interactions.

### 8.4. Random-Forest Safety Classification

A Random Forest classifier can be used to explore how model parameters relate to a defined safety or performance class.

Such a classification requires an explicit target-label definition and a documented training and validation procedure. Model predictions should not be interpreted as validated safety probabilities unless the classification framework has been calibrated and independently evaluated for that purpose.

### 8.5. Analysis of Variance (ANOVA)

ANOVA can be used to examine the statistical significance of parameter effects and interactions when the dataset and model specification support the method.

The analysis should document the response variable, factor definitions, binning of continuous predictors where applicable, interaction terms, and treatment of repeated earthquake records. Statistical significance should be distinguished from engineering significance.

These techniques provide complementary perspectives. No single correlation, classifier, or statistical test is sufficient to establish a universal ranking of parameter importance across all soil–structure systems.

## 9. Computational Environment

The project is located in the following directory of the repository:

```text
OPENSEES_SALAR/
└── EIGHT_ANALYSIS_DUCTILITY_DAMAGE_INDEX_&_ENERGY_DISSIPATION_CAPACITY_INDEX/
    └── 23_SOIL_STRUCTURE_MDOF_8_ANA/
        └── FY_KE_DUCT_OSF_SENSITIVITY_SEISMIC_SOIL_STRUCTURE_MDOF/
```

The analysis is developed around Python-based computational workflows. The exact dependencies should be taken from the imports used in the scripts.

Potential dependencies include:

- **Python** — scripting and workflow orchestration.
- **OpenSeesPy** — nonlinear finite-element and structural analysis.
- **NumPy** — numerical arrays and mathematical operations.
- **Pandas** — tabular data processing.
- **Matplotlib** — visualization of response histories and sensitivity results.
- **SciPy** — numerical and statistical calculations, where required.
- **Scikit-learn** — machine-learning methods, if used by the implementation.

Install only the packages required by the actual scripts. The repository should also identify the main execution script, required input files, supporting modules, and expected output files to facilitate reproducibility.

## 10. Verification and Reproducibility

Reliable sensitivity analysis requires consistent model definitions, numerical settings, and data processing.

Recommended verification practices include:

1. Document the units, baseline values, and ranges of all investigated parameters.
2. Verify the relationships between FY, KE, DUCT, and OSF in the nonlinear material model.
3. Confirm that soil and structural damping are assigned as intended.
4. Record ground-motion identifiers, preprocessing, and scaling factors.
5. Check numerical convergence and retain information about unsuccessful analyses.
6. Verify the definitions and units of all response quantities.
7. Ensure that comparisons use consistent analysis settings and response definitions.
8. Record software versions and random seeds where stochastic sampling is used.
9. Distinguish numerical failures from physically meaningful high-response outcomes.
10. Validate selected simulations against independent calculations, benchmark models, or established reference solutions where appropriate.

A sensitivity ranking is meaningful only when the model, input parameters, response measures, and statistical procedures are sufficiently well defined.

## 11. Scope and Limitations

The results of this project are conditional on the implemented soil–structure model, parameter ranges, constitutive assumptions, ground-motion dataset, and numerical procedures.

Important limitations to consider include:

- The adopted representation of soil and foundation flexibility.
- The calibration and applicability of the nonlinear material models.
- The physical consistency of the parameter combinations.
- The selection and scaling of the earthquake records.
- The treatment of damping and period calculations.
- The definitions and calibration of damage and reconstruction-cost indicators.
- The handling of numerical convergence failures.
- The statistical treatment of record-to-record variability and parameter interactions.

The computed sensitivities should not be generalized automatically to different structural configurations, soil conditions, earthquake datasets, or modeling assumptions. Engineering conclusions should be supported by the underlying simulation results and independent verification.

## 12. Intended Applications

This computational project may support:

- Investigation of nonlinear soil–structure interaction.
- Parametric studies of strength, stiffness, ductility, and overstrength.
- Evaluation of seismic response under multiple ground motions.
- Exploration of period-related changes and hysteretic behavior.
- Comparative analysis of energy- and ductility-related indicators.
- Development of statistical and machine-learning-based sensitivity tools.
- Educational demonstrations of Python and OpenSeesPy in earthquake engineering.
- Further development of performance-based structural analysis workflows.

The framework is intended for computational study and educational or research use. It does not replace project-specific structural assessment, geotechnical investigation, independent model validation, or code-compliant engineering design.

## 13. Author and Repository

**Author:** Salar Delavar Ghashghaei (Qashqai)

**GitHub Profile:** [github.com/salardelavar](https://github.com/salardelavar)

**Main Repository:** [OPENSEES_SALAR](https://github.com/salardelavar/OPENSEES_SALAR)

The repository shares computational examples and Python-based tools for structural engineering, nonlinear analysis, and earthquake engineering.

## 14. Disclaimer

This project is provided for educational and research purposes. Although care has been taken in developing the computational procedures, the scripts may contain errors, numerical limitations, or modeling assumptions that require independent verification.

Users are responsible for reviewing the implementation, validating the numerical results, and determining whether the adopted methods are appropriate for their intended applications. The author makes no guarantee that the scripts are suitable for a specific design, assessment, or safety-critical engineering decision.

## Keywords

OpenSeesPy, Python, soil–structure interaction, SSI, MDOF system, nonlinear structural analysis, earthquake engineering, seismic time-history analysis, ground-motion analysis, sensitivity analysis, yield strength, elastic stiffness, ductility capacity, overstrength factor, period elongation, energy dissipation, equivalent viscous damping, ductility damage index, reconstruction cost, ANOVA, Random Forest, structural performance.

 

# 200 SEISMIC GROUND MOTIONS ANALYSIS
![alt text](SEISMIC_GROUND_MOTION_ANALYSIS.png)


