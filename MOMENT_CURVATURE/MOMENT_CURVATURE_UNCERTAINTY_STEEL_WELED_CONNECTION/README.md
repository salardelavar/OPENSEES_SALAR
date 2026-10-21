# PROBABILISTIC MOMENT–CURVATURE ANALYSIS OF STEEL WELDED CONNECTIONS UNDER UNCERTAINTY

![alt text](https://github.com/salardelavar/OPENSEES_SALAR/blob/main/MOMENT_CURVATURE/MOMENT_CURVATURE_UNCERTAINTY_STEEL_WELED_CONNECTION/MOMENT_CURVATURE_UNCERTAINTY_STEEL_WELDED_CONNECTION.png)

### OpenSees | Python | Monte Carlo Simulation | Structural Reliability

**Author:** Salar Delavar Ghashghaei (Qashqai)  
**Repository:** [OPENSEES_SALAR](https://github.com/salardelavar/OPENSEES_SALAR)  
**Project Directory:** `MOMENT_CURVATURE/MOMENT_CURVATURE_UNCERTAINTY_STEEL_WELED_CONNECTION`

---

## 1. Overview

This project focuses on the numerical investigation of the moment–curvature response of steel welded connections under uncertainty using Python and OpenSees.

Moment–curvature analysis provides a fundamental framework for examining the nonlinear flexural response of structural sections and connection models. When uncertainties in material properties, geometric parameters, connection characteristics, or other model inputs are considered, the resulting response may no longer be adequately represented by a single deterministic curve.

This project explores a probabilistic approach to structural response analysis by combining nonlinear numerical modeling with Monte Carlo simulation. The objective is to investigate how uncertainty in model parameters can influence the predicted moment–curvature relationship and the associated structural response characteristics.

The computational framework is intended for educational and research applications in structural mechanics, nonlinear finite-element analysis, uncertainty quantification, and steel connection behavior.

## 2. Engineering Motivation

Steel welded connections are important components of structural load paths. Their mechanical response can influence the stiffness, strength, deformation capacity, and overall behavior of steel structures.

A conventional deterministic analysis evaluates the response using one specified set of input parameters. However, real structural systems exhibit variability arising from material properties, fabrication tolerances, geometric imperfections, welding-related characteristics, and uncertainty in modeling assumptions.

A probabilistic analysis can help investigate the implications of these variations by evaluating multiple realizations of the model rather than relying exclusively on one nominal configuration.

The moment–curvature relationship is particularly useful because it describes the evolution of bending resistance with curvature and can reveal important nonlinear response characteristics.

## 3. Main Objectives

The project is intended to support the following investigations:

- Develop a computational framework for moment–curvature analysis of steel welded connection models.
- Investigate the influence of uncertain input parameters on nonlinear structural response.
- Apply Monte Carlo simulation to evaluate multiple realizations of uncertain model inputs.
- Examine the variability of moment–curvature curves across simulated realizations.
- Investigate statistical characteristics of response quantities derived from the analysis.
- Explore the potential application of probabilistic response analysis to structural performance assessment and uncertainty quantification.
- Provide a Python-based computational example that can be extended for further structural-engineering research.

The parameters investigated and the specific outputs available depend on the implementation in the accompanying source files.

## 4. Theoretical Background

### 4.1 Moment–Curvature Relationship

The moment–curvature relationship describes the relationship between the bending moment \(M\) and the curvature \(\phi\) of a structural section or an appropriately defined connection model.

\[
M = f(\phi)
\]

For small deformations in the linear-elastic range, the relationship can be approximated by

\[
M \approx EI\phi
\]

where:

- \(M\) is the bending moment.
- \(E\) is Young's modulus.
- \(I\) is the relevant second moment of area.
- \(\phi\) is the curvature.

The expression above represents the classical linear-elastic section relationship under its usual assumptions. The response of a welded connection model may require additional constitutive relationships and connection-specific modeling assumptions.

Beyond the elastic range, nonlinear behavior may arise from material yielding, changes in tangent stiffness, local deformation mechanisms, and other effects represented in the numerical model.

The resulting moment–curvature response can be used to investigate initial stiffness, nonlinear response, peak resistance, and post-peak behavior when these characteristics are captured by the adopted model and analysis procedure.

### 4.2 Uncertainty Quantification

Let the uncertain input parameters be represented by the random vector

\[
\mathbf{X} =
\left[X_1,X_2,\ldots,X_n\right]
\]

The calculated response can then be expressed as

\[
M(\phi,\mathbf{X})
\]

where the response depends on curvature and the realization of the uncertain input parameters.

Different realizations of \(\mathbf{X}\) may produce different moment–curvature curves. Their variation provides information about the sensitivity of the predicted response to the specified uncertainties.

The statistical interpretation depends on the assumed probability distributions, parameter dependencies, numerical model, and sampling procedure.

### 4.3 Monte Carlo Simulation

Monte Carlo simulation is a numerical technique for propagating uncertainty through a computational model.

For a simulation involving \(N\) realizations, the general procedure is:

1. Define the uncertain input parameters and their probability distributions.
2. Generate a realization of the input parameters.
3. Construct or update the numerical model using that realization.
4. Execute the moment–curvature analysis.
5. Extract the response quantities of interest.
6. Repeat the procedure for the prescribed number of realizations.
7. Evaluate the statistical distribution and variability of the calculated responses.

For a response quantity \(Y\), the sample mean can be estimated as

\[
\bar{Y}=\frac{1}{N}\sum_{i=1}^{N}Y_i
\]

and the sample standard deviation as

\[
s_Y=
\sqrt{
\frac{1}{N-1}
\sum_{i=1}^{N}(Y_i-\bar{Y})^2
}
\]

These estimators are applicable when the corresponding response values have been successfully obtained from the simulation.

The adequacy of the simulation should be evaluated using convergence checks, sufficient sampling, and appropriate treatment of failed or incomplete analyses.

## 5. Computational Framework

The project combines numerical structural analysis and probabilistic computation.

### 5.1 Model Definition

The structural model defines the geometry, material behavior, connection representation, boundary conditions, and loading conditions required for the analysis.

The actual model assumptions must be established from the source code. In particular, the representation of the welded connection should be examined to determine which deformation mechanisms and failure modes are explicitly included.

### 5.2 Nonlinear Moment–Curvature Analysis

The numerical analysis evaluates the relationship between bending moment and curvature according to the implemented loading protocol and constitutive models.

Important implementation details include:

- Material constitutive relationships.
- Connection and section idealization.
- Loading and boundary conditions.
- Curvature increments or displacement-control parameters.
- Nonlinear solution algorithms and convergence tolerances.
- Treatment of numerical nonconvergence.

The reliability of the resulting response depends on these modeling and numerical choices.

### 5.3 Probabilistic Simulation

The uncertainty-analysis procedure generates multiple model realizations and evaluates their responses.

For a reproducible study, the following information should be documented explicitly:

- The uncertain parameters.
- The probability distribution assigned to each parameter.
- The distribution parameters and units.
- Any assumed correlations between uncertain variables.
- The number of Monte Carlo realizations.
- The random-number seed, if controlled.
- The criteria used to identify successful and failed analyses.

These details are essential for interpreting the statistical results and reproducing the analysis.

## 6. Potential Response Quantities

Depending on the available numerical outputs, the results can be examined in terms of the following quantities.

| Response quantity | Engineering interpretation |
|---|---|
| Initial tangent stiffness | Initial resistance to curvature |
| Moment–curvature response | Evolution of bending resistance with curvature |
| Yield-related response | Transition from predominantly elastic to nonlinear behavior, subject to the adopted definition |
| Maximum bending moment | Peak resistance attained during the simulated loading protocol |
| Curvature at a specified response point | Deformation measure associated with a defined response criterion |
| Post-peak response | Reduction or evolution of resistance after the peak, if captured |
| Mean response | Average across the successful simulation realizations |
| Standard deviation | Measure of dispersion in the simulated response |
| Response percentiles | Distribution-based characterization of response variability |

These quantities are potential outputs, not a claim that every item is calculated by the current implementation.

## 7. Interpretation of Probabilistic Results

A probabilistic moment–curvature analysis should be interpreted in terms of both the response of individual realizations and the statistical behavior of the complete simulation set.

**Response variability:** A spread of moment–curvature curves indicates that the specified input uncertainty influences the predicted response.

**Stiffness variability:** Differences in the initial slopes may indicate sensitivity to uncertain elastic properties or other model parameters affecting initial stiffness.

**Strength variability:** Differences in peak moments may reveal sensitivity of the predicted resistance to the uncertain parameters represented in the model.

**Deformation variability:** Differences in curvature at a defined response point can help characterize uncertainty in deformation-related performance measures.

**Post-peak variability:** If the model captures post-peak behavior, the simulations may help examine the variability of the descending response and the numerical stability of that regime.

These interpretations must remain consistent with the adopted connection model. A numerical moment–curvature response alone does not establish that all physical failure mechanisms of a real welded connection have been represented.

## 8. Suggested Statistical Visualization

For a comprehensive probabilistic investigation, the following figures may be useful where the required data are available:

1. Individual moment–curvature curves for the Monte Carlo realizations.
2. The nominal or reference response compared with the ensemble of simulated responses.
3. The estimated mean moment–curvature curve.
4. Response dispersion represented by standard-deviation bands or selected percentiles.
5. Histograms of peak moment or other extracted response parameters.
6. Empirical cumulative distribution functions for selected engineering demand or resistance measures.
7. Convergence of estimated statistical quantities as the number of simulations increases.

A mean curve should not automatically be interpreted as the response of a physically realizable connection. Likewise, percentile bands should be calculated at consistent curvature values and from clearly defined populations of successful analyses.

## 9. Applications in Structural Engineering

The computational approach may support research in:

- Probabilistic assessment of nonlinear structural response.
- Sensitivity studies involving uncertain material or connection parameters.
- Statistical characterization of moment–curvature behavior.
- Numerical investigations of steel connection mechanics.
- Monte Carlo methods for structural engineering.
- Development and verification of uncertainty-aware analysis procedures.
- Educational studies of deterministic and probabilistic modeling.

Further applications to reliability assessment require an explicitly defined limit-state function, a suitable failure criterion, and justified probabilistic input models.

## 10. Software and Dependencies

The project is associated with the Python and OpenSees ecosystem.

The required software and packages should be verified against the source files before execution.

Potential dependencies include:

- Python.
- OpenSeesPy, if the model is implemented through its Python interface.
- NumPy for numerical calculations.
- Matplotlib for plotting.
- Pandas for tabular data processing, if used by the scripts.
- SciPy for statistical or numerical operations, if used by the scripts.

Not every listed package is necessarily required by the current implementation.

### General execution procedure

1. Download or clone the repository.
2. Navigate to the project directory.
3. Install the dependencies imported by the relevant source file.
4. Inspect the model parameters, units, uncertainty definitions, and simulation settings.
5. Execute the analysis script or notebook.
6. Review the generated outputs and confirm that the numerical analyses completed successfully.
7. Check convergence and statistical stability before interpreting the results.

The precise execution command depends on the filenames and organization of the actual project files.

## 11. Verification and Limitations

Monte Carlo simulation does not, by itself, guarantee accurate structural predictions. Its results are conditional on the quality of the numerical model and the assumptions used to define uncertainty.

Before using the results in research or engineering assessments, consider the following checks:

- Verify material properties, geometric dimensions, and unit consistency.
- Confirm that the connection idealization is appropriate for the intended investigation.
- Check the sensitivity of the results to numerical discretization and solution settings.
- Investigate nonconvergence and distinguish numerical failure from a physically meaningful limit state.
- Examine the convergence of statistical estimates with increasing sample size.
- Confirm that the assumed probability distributions are supported by experimental data, standards, or defensible engineering assumptions.
- Avoid extrapolating beyond the range of behavior represented by the numerical model.

**Important:** This project is provided for educational and research purposes. Its results should not be treated as validated predictions of real connection strength or failure probability without appropriate verification, validation, and supporting evidence.

## 12. Reproducibility and Future Development

The framework can be extended by documenting the uncertain variables, their distributions, and their physical justification; conducting sensitivity and convergence studies; and comparing numerical predictions against analytical solutions, experimental measurements, or independently verified finite-element models where appropriate.

Potential future developments include:

- Sensitivity analysis to identify influential uncertain parameters.
- Confidence intervals for selected statistical estimates.
- Reliability calculations based on clearly defined limit states.
- Comparison of alternative uncertainty-propagation techniques.
- Validation against relevant experimental moment–curvature or connection-response data.

Such extensions should be considered research directions unless they are already implemented and verified in the project.

## 13. Author and Open-Source Educational Work

**Salar Delavar Ghashghaei (Qashqai)** develops and shares computational examples involving structural engineering, earthquake engineering, nonlinear analysis, numerical methods, and uncertainty quantification.

The broader repository contains Python, OpenSees, and related computational examples intended to support learning, experimentation, and knowledge-sharing within the engineering community.

- **GitHub:** https://github.com/salardelavar/OPENSEES_SALAR
- **Technical weblog:** https://salardqashqai.blogfa.com/

## 14. Disclaimer

This code is distributed for educational and research purposes. The implementation may contain errors, simplifying assumptions, or numerical limitations. Users are responsible for independently verifying the source code, model assumptions, numerical procedures, and results before relying on them.

The results should not be used as the sole basis for structural design, safety assessments, or engineering decisions without appropriate independent validation.

---

**Keywords:** OpenSeesPy, Python, moment–curvature analysis, steel welded connections, uncertainty quantification, Monte Carlo simulation, nonlinear structural analysis, probabilistic mechanics, structural reliability, finite-element modeling.



