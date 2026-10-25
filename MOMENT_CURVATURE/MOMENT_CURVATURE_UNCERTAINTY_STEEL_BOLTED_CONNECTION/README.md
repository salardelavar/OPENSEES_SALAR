# PROBABILISTIC MOMENT–CURVATURE ANALYSIS OF STEEL BOLTED CONNECTIONS UNDER UNCERTAINTY

**Author:** Salar Delavar Ghashghaei (Qashqai)  
**Repository:** [OPENSEES_SALAR](https://github.com/salardelavar/OPENSEES_SALAR)  
**Project Category:** Structural Engineering | Nonlinear Analysis | Uncertainty Quantification | Python | OpenSees

---

## 1. Project Overview

This project focuses on the numerical investigation of moment–curvature behavior in steel bolted connections considering uncertainty in the parameters that govern their mechanical response. By combining nonlinear structural analysis with probabilistic simulation, the framework is intended to explore how variations in material properties, connection characteristics, and modeling parameters may affect the predicted flexural response.

Moment–curvature analysis provides a useful basis for investigating the relationship between bending moment and curvature, identifying changes in stiffness, evaluating nonlinear response, and examining the sensitivity of structural behavior to uncertain input parameters.

The project is part of the broader **OPENSEES_SALAR** collection, which provides Python and OpenSees-based educational implementations for structural analysis, earthquake engineering, numerical simulation, and engineering computation.

## 2. Engineering Objectives

The principal objectives of this project are:

- Investigate the nonlinear moment–curvature response of steel bolted connection models.
- Examine the influence of uncertain input parameters on the calculated structural response.
- Develop a computational framework for repeated analyses using Python and OpenSees, where implemented by the associated scripts.
- Explore the variability of moment–curvature curves under different realizations of uncertain parameters.
- Identify parameters that may substantially influence initial stiffness, yielding behavior, peak moment resistance, and post-yield response.
- Support comparative and sensitivity studies for improving the understanding of connection behavior.
- Provide an educational foundation for further development of probabilistic connection analysis.

The objectives represent the intended scope of the study. Specific capabilities and outputs depend on the implementation contained in the project files.

## 3. Mechanical Background

### 3.1 Moment–Curvature Relationship

The moment–curvature relationship describes the flexural response of a structural section or an appropriately formulated connection model.

\[
M = f(\phi)
\]

where:

- \(M\) = bending moment.
- \(\phi\) = curvature.
- \(f(\phi)\) = the modeled moment–curvature relationship.

The response may include an initial approximately linear range, nonlinear behavior, yielding, stiffness degradation, and post-peak response, depending on the constitutive laws, connection idealization, and loading protocol.

The initial tangent stiffness of a moment–curvature curve can be expressed as

\[
K_{\phi,0}=\left.\frac{dM}{d\phi}\right|_{\phi\rightarrow 0}
\]

This quantity characterizes the initial flexural response of the modeled system. It should not automatically be interpreted as the rotational stiffness of an entire bolted joint unless the model's kinematics and deformation measures establish that equivalence.

### 3.2 Modeling Steel Bolted Connections

The response of a bolted connection may depend on several interacting mechanisms, including:

- Elastic and inelastic behavior of the connected steel components.
- Bolt strength and deformation characteristics.
- Plate geometry, thickness, and material properties.
- Contact, bearing, slip, and separation between connected components.
- Connection geometry and boundary conditions.
- The adopted representation of bolts, plates, interfaces, and load transfer.

Not every model explicitly captures all these mechanisms. The interpretation of the results must therefore remain consistent with the actual numerical formulation.

A simplified section-level model may characterize flexural response efficiently, whereas a detailed component or contact model may be required to investigate local bolt forces, bearing deformation, slip, or other connection-specific limit states.

## 4. Uncertainty Quantification

Structural response is influenced by uncertainty in material properties, geometric parameters, connection behavior, and numerical idealization. Probabilistic analysis provides a framework for investigating the consequences of these variations.

Let the uncertain parameter vector be

\[
\mathbf{X}=[X_1,X_2,\ldots,X_n]
\]

where each \(X_i\) represents an uncertain model input.

The resulting response can be represented as

\[
M(\phi,\mathbf{X})
\]

For each realization of \(\mathbf{X}\), the numerical model may produce a different moment–curvature curve. The collection of these curves can be used to investigate response variability over the sampled parameter space.

### 4.1 Potential Sources of Uncertainty

Depending on the implemented model, relevant uncertain parameters may include:

| Parameter group | Representative variables |
|---|---|
| Material properties | Steel yield strength, elastic modulus, ultimate strength |
| Geometry | Plate thickness, connection dimensions, component eccentricities |
| Bolt properties | Bolt strength, diameter, or other modeled mechanical parameters |
| Connection behavior | Rotational stiffness, slip-related parameters, or contact properties |
| Constitutive response | Post-yield stiffness and other material-model parameters |
| Modeling assumptions | Selected idealizations and boundary-condition parameters |

These are candidate uncertainty sources, not a claim that every variable is randomized in the current implementation.

### 4.2 Monte Carlo Simulation

If Monte Carlo simulation is used, each analysis realization is generated from a specified joint probability model for the uncertain inputs.

For \(N\) realizations,

\[
\mathbf{X}^{(i)}\sim p_{\mathbf{X}}(\mathbf{x}),
\qquad i=1,\ldots,N
\]

The numerical model then evaluates

\[
M^{(i)}(\phi)=M(\phi,\mathbf{X}^{(i)})
\]

The resulting ensemble can be used to estimate response statistics, including the sample mean,

\[
\overline{M}(\phi)=
\frac{1}{N}\sum_{i=1}^{N}M^{(i)}(\phi)
\]

and sample standard deviation,

\[
s_M(\phi)=
\sqrt{\frac{1}{N-1}
\sum_{i=1}^{N}
\left[M^{(i)}(\phi)-\overline{M}(\phi)\right]^2}
\]

These expressions assume that the responses are evaluated at a consistent curvature value and that the relevant realizations have valid results at that value.

Response percentiles can also be calculated to describe the spread of the simulated results. Such percentiles are statistical summaries, not automatically design resistance values or code-prescribed characteristic strengths.

## 5. Potential Response Quantities

The following quantities can be investigated if the model and post-processing routines provide the required data.

### 5.1 Initial Flexural Stiffness

The initial slope of the moment–curvature response can be used to compare the initial stiffness of different realizations.

### 5.2 Yielding Response

An appropriate yield criterion can be used to estimate the onset of yielding. The definition should be stated explicitly, particularly when the response does not exhibit a sharp yield point.

### 5.3 Peak Moment Resistance

The maximum moment reached within the analyzed loading range is

\[
M_{\max}^{(i)}=\max_{\phi\in\mathcal{D}_i}M^{(i)}(\phi)
\]

where \(\mathcal{D}_i\) is the curvature domain successfully analyzed for realization \(i\).

This quantity is the maximum computed moment over the analyzed domain; it is not necessarily the physical ultimate resistance if the simulation terminates before the governing limit state.

### 5.4 Post-Yield Behavior

The post-yield response may be examined to investigate changes in tangent stiffness, deformation capacity, and possible loss of resistance.

The reliability of these interpretations depends on the constitutive model, convergence behavior, loading range, and representation of connection failure mechanisms.

### 5.5 Response Variability

Statistical summaries of the simulated curves can be used to quantify uncertainty in the predicted response. A narrow response distribution indicates less variability for the selected inputs and modeling assumptions, whereas a wider distribution indicates greater variability.

## 6. Recommended Computational Workflow

The following workflow provides a reproducible structure for probabilistic moment–curvature analysis.

1. **Define the structural model.** Specify the geometry, material constitutive laws, connection idealization, boundary conditions, and loading procedure.

2. **Identify uncertain parameters.** Select variables that are physically justified and relevant to the intended investigation.

3. **Define probability distributions.** Specify the distribution type, statistical parameters, units, and dependencies for each uncertain variable.

4. **Generate input realizations.** Sample the uncertain parameters using the selected simulation method.

5. **Execute nonlinear analyses.** Run the structural model for each realization while recording convergence status and analysis limits.

6. **Extract response quantities.** Store the moment–curvature history and any additional quantities supported by the implementation.

7. **Evaluate statistical variability.** Calculate relevant statistics and compare response distributions.

8. **Perform sensitivity analysis.** Investigate which uncertain parameters have the greatest influence on selected response measures.

9. **Validate the numerical procedure.** Check representative results against analytical solutions, benchmark examples, independent calculations, or experimental evidence when available.

10. **Document limitations.** Report assumptions, failed analyses, parameter ranges, and restrictions on interpretation.

## 7. Numerical Reliability and Validation

Probabilistic nonlinear analysis requires more than repeated execution of a numerical model. The reliability of the resulting statistics depends on the validity of individual analyses.

Recommended checks include:

- Confirm that the baseline model produces a mechanically reasonable response.
- Verify the sign conventions and units used for moment and curvature.
- Confirm that sampled parameters remain within physically meaningful ranges.
- Distinguish numerical nonconvergence from a physically meaningful connection failure.
- Record failed realizations rather than silently discarding them.
- Use a consistent curvature range when comparing curves.
- Check whether the number of simulations is sufficient for stable estimates of the quantities of interest.
- Examine sensitivity to sampling assumptions and numerical solution settings.
- Compare representative results with independent calculations or suitable benchmarks.

If failed realizations are excluded, the resulting statistics may be biased, particularly when failure to converge is associated with severe nonlinear response. The treatment of unsuccessful analyses should therefore be reported explicitly.

## 8. Interpretation and Scope

This project is intended to support computational investigation and engineering education. Its results should be interpreted in the context of the model assumptions and the uncertain variables actually considered.

A probabilistic distribution of moment–curvature response does not, by itself, establish the probability of connection failure. A failure-probability calculation requires a clearly defined limit-state function and a suitable treatment of uncertainty.

For example, a limit state may be formulated as

\[
g(\mathbf{X})=R(\mathbf{X})-S(\mathbf{X})
\]

where \(R\) denotes resistance and \(S\) denotes demand, provided these quantities and their relationship are defined appropriately for the connection being assessed.

The failure probability is then

\[
P_f=P[g(\mathbf{X})\leq 0]
\]

This formulation is only applicable when the selected resistance, demand, and failure criterion are physically meaningful and consistently evaluated. It should not be presented as an output of this repository unless the corresponding limit-state analysis has actually been implemented.

## 9. Software and Reproducibility

The project belongs to the OPENSEES_SALAR repository, which contains educational structural analysis scripts written primarily in Python, with OpenSees-based implementations across the collection.

For reproducible studies, users should document:

- Python and OpenSeesPy versions.
- Required third-party packages.
- Units and sign conventions.
- Random seed and sampling procedure, where applicable.
- Input distributions and parameter correlations.
- Number of simulation realizations.
- Nonlinear solution algorithms and convergence criteria.
- Output file formats and post-processing procedures.

The exact installation commands and file-specific execution instructions should be added after confirming the dependencies and entry-point scripts in this directory.

## 10. Suggested Future Development

Potential extensions include:

- Automated statistical summaries of moment–curvature curves.
- Confidence intervals and percentile-based response envelopes.
- Sensitivity analysis for material, geometric, and connection parameters.
- Convergence diagnostics for Monte Carlo estimates.
- Comparison of alternative probability distributions.
- Correlation modeling between uncertain material and geometric properties.
- Benchmarking against experimental moment–rotation or moment–curvature data.
- Reliability analysis using explicitly defined connection limit states.
- Improved reporting of failed simulations and numerical uncertainty.

These extensions are recommendations for future development and should not be interpreted as existing features unless supported by the project code.

## 11. Author and Repository

**Author:** Salar Delavar Ghashghaei (Qashqai)

**GitHub repository:** [OPENSEES_SALAR](https://github.com/salardelavar/OPENSEES_SALAR)

**Project directory:** [MOMENT_CURVATURE_UNCERTAINTY_STEEL_BOLTED_CONNECTION](https://github.com/salardelavar/OPENSEES_SALAR/tree/main/MOMENT_CURVATURE/MOMENT_CURVATURE_UNCERTAINTY_STEEL_BOLTED_CONNECTION)

The repository provides educational examples for exploring numerical methods and structural engineering applications through Python and OpenSees.

## 12. Disclaimer

This project is provided for educational, research, and numerical-exploration purposes. The scripts and their results may contain errors, simplifications, or modeling limitations. Users should independently verify the numerical implementation, constitutive assumptions, boundary conditions, and interpretation of results before applying the methods to engineering decisions.

The examples should not be treated as a substitute for validated engineering software, applicable design standards, experimental verification, or professional engineering judgment.

---

**Keywords:** Moment–Curvature Analysis, Steel Bolted Connections, Structural Mechanics, Nonlinear Analysis, Uncertainty Quantification, Monte Carlo Simulation, Probabilistic Structural Analysis, Sensitivity Analysis, Python, OpenSees, OpenSeesPy, Structural Engineering.

