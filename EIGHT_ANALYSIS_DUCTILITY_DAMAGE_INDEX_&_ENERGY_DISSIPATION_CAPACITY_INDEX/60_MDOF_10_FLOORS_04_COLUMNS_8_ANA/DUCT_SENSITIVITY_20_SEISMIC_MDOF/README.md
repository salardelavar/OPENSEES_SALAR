# SENSITIVITY ANALYSIS OF COLUMN DUCTILITY RATIO  AND OPENSEES VIA 20 SEISMIC GROUND MOTIONS
# تحلیل حساسیت نسبت شکل‌پذیری ستون و اوپنسیس از طریق ۲۰ حرکت زمین‌لرزه‌



![alt text](https://github.com/salardelavar/OPENSEES_SALAR/blob/main/EIGHT_ANALYSIS_DUCTILITY_DAMAGE_INDEX_%26_ENERGY_DISSIPATION_CAPACITY_INDEX/60_MDOF_10_FLOORS_04_COLUMNS_8_ANA/DUCT_SENSITIVITY_20_SEISMIC_MDOF/COVER.png) 

# EQUIVALENT SDOF SYSTEM DERIVATION VIA DISPLACEMENT-BASED SEISMIC DESIGN PROCEDURE WITH PUSHOVER ANALYSIS

![alt text](https://github.com/salardelavar/OPENSEES_SALAR/blob/main/EIGHT_ANALYSIS_DUCTILITY_DAMAGE_INDEX_%26_ENERGY_DISSIPATION_CAPACITY_INDEX/60_MDOF_10_FLOORS_04_COLUMNS_8_ANA/DUCT_SENSITIVITY_20_SEISMIC_MDOF/COVER_DISPLACEMENT_BASED_PUSHOVER.png) 

Assume that a multi-degree-of-freedom structure is subjected to nonlinear dynamic analysis under twenty different ground motion records. To investigate the influence of column ductility on seismic performance, a sensitivity analysis is performed by systematically varying each column’s ductility ratio from 20 to 50 in 21 steps. For each ductility level, a pushover analysis provides the capacity curve and equivalent SDOF parameters, a cyclic displacement analysis establishes a reference hysteretic energy, and twenty nonlinear time-history analyses are conducted using different ground motions. From these analyses, the Energy Dissipation Capacity Index (EDCI) and other response metrics such as maximum displacement, velocity, acceleration, damage index, over-strength factor, ductility ratio, and equivalent viscous damping ratio are computed. The median response over the twenty ground motions is then evaluated for each ductility level. Finally, the results are assessed through trend plots, 3D contour surfaces, correlation heatmaps, Random Forest regression, and ANOVA sensitivity analysis to identify the optimal column ductility ratio that maximizes energy dissipation capacity while maintaining structural safety.

# SENSITIVITY ANALYSIS BY CHANGING COLUMN DUCTILITY RATIO
1. Sets MAT_TYPE = 'INELASTIC' and sweeps the COLUMN DUCTILITY RATIO from 20.0 to 50.0
   in 21 linear steps (DUCT_MIN -> DUCT_MAX).  COL_DUCT stores each swept value.
2. Declares ~20 accumulator lists: raw per-record responses, plus *_MED lists that
   will hold one MEDIAN value per ductility level.
3. Starts a CPU timer.  Outer loop `for JJ in range(0, 21)` computes DUCT and resets
   the per-ductility raw accumulators (DISP, VELO, ACC, EDCI, DII, OMEGA, MU, RR, EDVR).
4. STEP 1 -- PUSHOVER: MDOF(DUCT, MAT_TYPE, TOTAL_MASS, 'PUSHOVER', II) returns the
   capacity curve (disp_PUSH, reaction_PUSH) and the equivalent-SDOF system
   (mass, stiffness, period) via displacement-based pushover reduction.
5. Stores SDOF_ef_* arrays and PERIOD_PUSH = max(PERIOD_MAX_PUSH)  (secant period).
6. STEP 2 -- CYCLIC_DISPLACEMENT: benchmark constant-amplitude hysteresis loop.
   S055.EQULIVALENT_VISCOUS_DAMPING_RATIO_FUN returns zeta_CP (Jacobsen area method).
7. STEP 3 -- SEISMIC LOOP `for II in range(40, 60)`: 20 nonlinear time-history runs.
8. Each run unpacks time/reaction/disp/velo/acc/DI/stiffness/period arrays from MDOF().
9. S12.DUCTILITY_DAMAGE_INDEX_FUN fits the pushover curve and returns Park-Ang style
   DIx, over-strength Omega_0, ductility mu, and behaviour coefficient R.
10. S10.ENERGY_DISSIPATION_CAPACITY_INDEX compares seismic vs cyclic dissipated energy
   -> EDCI (%).  Peak DISP/VELO/ACC also appended.
11. Ratio metrics appended: PERIOD_RATIO, DI_RATIO, REACTION_RATIO, DISP_RATIO
   (each = SEISMIC / PUSHOVER).
12. S055 again computes zeta_SEI for the seismic loop -> EDVR and EDVR_RATIO = zeta_SEI/zeta_CP.
13. STEP 4 -- after 20 records, np.median() reduces each raw list to ONE robust value
   (DISP_MED, DI_RATIO_MED, EDCI_MED, EDVR_MED, ...), the design-point response.
14. Outer loop ends after 21 design points.  Total CPU time is printed.
15. PLOTTING: ~17 calls to S01.PLOT_SCATTER(COL_DUCT, Y, ...) draw scatter + polynomial
   fit of order 1/3/7 for every response vs. the ductility ratio.
16. 5 calls to S11.PLOT_CONTOUR_3D_2D_FUN build 2D contours + 3D surfaces pairing two
   responses (e.g. DISP_RATIO vs DI_RATIO) against the ductility ratio axis.
17. MACHINE LEARNING: builds a pandas DataFrame with features DISP, VEL, ACC,
   REACTION_RATIO, EDVR_RATIO, DI_RATIO, then calls S01.RANDOM_FOREST(df) to train a
   RandomForestClassifier (safe/unsafe) + Regressor (safety likelihood) and print
   accuracy, MSE, R2, and feature importances.
18. S01.PLOT_HEATMAP(df) draws the Pearson correlation matrix of the six features.
19. S099.SENSITIVITY_HEATMAP_FUN(X=[DISP,VEL,ACC], Y=[EDVR_RATIO,DUCT,DI_RATIO],
   method='pearson') returns a 3x3 sensitivity-coefficient heatmap.
20. S100.ANOVA_SENSITIVITY_FUN(df_sens, output_col="DISP", param_cols=[VEL,ACC,
   EDVR_RATIO,DI_RATIO], n_bins=4) bins continuous predictors and runs classical ANOVA.
21. Two ANOVA runs: (a) MAIN EFFECTS only, (b) MAIN + TWO-WAY INTERACTIONS; each prints
   a SumSq / df / F / p-value table and shows the diagnostic plot.

22. Engineering reading: plots reveal whether higher column ductility reduces damage
   (DII), damping (EDVR), energy index (EDCI) and peak demand -- the core question of
   performance-based seismic design.

23. ML + correlation + ANOVA together rank which EDP (disp, vel or acc) drives damage.

# 200 SEISMIC GROUND MOTIONS ANALYSIS
![alt text](https://github.com/salardelavar/OPENSEES_SALAR/blob/main/EIGHT_ANALYSIS_DUCTILITY_DAMAGE_INDEX_%26_ENERGY_DISSIPATION_CAPACITY_INDEX/60_MDOF_10_FLOORS_04_COLUMNS_8_ANA/DUCT_SENSITIVITY_20_SEISMIC_MDOF/SEISMIC_GROUND_MOTION_ANALYSIS.png)  


# COLUMN DUCTILITY SENSITIVITY ANALYSIS UNDER 20 EARTHQUAKE GROUND MOTIONS

## 1. Overview

This project investigates the influence of column ductility on the nonlinear seismic response of a ten-story, multi-degree-of-freedom (MDOF) structural system using **OpenSeesPy and Python**.

The framework combines nonlinear structural analysis, earthquake time-history simulations, displacement-based pushover analysis, energy dissipation evaluation, damage assessment, and statistical sensitivity analysis. Its primary objective is to examine how variations in column ductility influence global structural response and seismic performance under multiple ground-motion records.

The repository also includes supporting numerical functions, visualization utilities, ground-motion preprocessing scripts, and supplementary technical documentation.

## 2. Research Objectives

The principal objectives are:

- Investigate the sensitivity of nonlinear structural response to column ductility.
- Evaluate the seismic behavior of a ten-story MDOF structural model.
- Perform nonlinear time-history analyses under multiple earthquake ground motions.
- Quantify structural displacement, base shear, ductility demand, and damage-related response parameters.
- Investigate hysteretic energy dissipation and equivalent viscous damping.
- Evaluate effective structural properties using pushover and free-vibration-based procedures.
- Apply statistical techniques to investigate relationships between ductility and structural response.
- Produce numerical results and visualizations suitable for engineering interpretation and further research.

## 3. Structural Model

The principal model represents a ten-story, multi-degree-of-freedom structure with four columns, as defined in the main OpenSeesPy analysis script.

The structural system is analyzed using nonlinear numerical procedures to investigate its response under prescribed loading and earthquake excitation.

The principal modeling components include:

- Structural mass and stiffness representation.
- Nonlinear material and hysteretic behavior.
- Column ductility parameters.
- Modal and eigenvalue analysis.
- Nonlinear static pushover analysis.
- Nonlinear dynamic time-history analysis.
- Response-history extraction and post-processing.

**Important:** The precise interpretation of the four-column configuration, element connectivity, boundary conditions, floor constraints, and material properties should be established from the model definition in the Python source code. The number of columns alone does not uniquely define the structural topology.

## 4. Analysis Framework

The workflow is organized into the following stages.

### 4.1. Model Definition and Modal Analysis

The structural model is established in OpenSeesPy. Its dynamic properties are evaluated using the relevant eigenvalue-analysis procedures.

Depending on the implemented model, these calculations can provide natural frequencies, vibration periods, and mode shapes.

### 4.2. Column Ductility Parameter Study

The main script investigates a range of column ductility parameters.

The ductility parameter is varied systematically to investigate its relationship with structural response. The implemented parameter values and number of cases should be taken directly from the current version of the main Python script.

For a controlled sensitivity study, other model parameters should remain consistent unless their variation is explicitly part of the research design.

### 4.3. Nonlinear Pushover Analysis

A nonlinear static analysis is used to obtain the structural capacity relationship, typically expressed as base shear versus control-node displacement.

The capacity curve can be used to investigate:

- Maximum base shear.
- Displacement capacity.
- Yield and post-yield response.
- Effective stiffness.
- Displacement ductility.
- Equivalent structural properties, where the adopted conversion method supports them.

The accuracy of these quantities depends on the structural model, loading pattern, control node, convergence criteria, and post-processing definitions.

### 4.4. Earthquake Time-History Analysis

The framework includes nonlinear dynamic analysis under multiple earthquake ground motions.

For each ground-motion case, the structural response is evaluated using the acceleration history and the model's dynamic characteristics.

Depending on the record and implemented analysis settings, the response histories can be used to obtain peak displacement, velocity, acceleration, base reaction, and other response quantities.

The main sensitivity analysis uses its configured set of earthquake cases. The separate ground-motion preprocessing utility supports a broader set of input files; its total number of records should not be confused with the number actually used in the sensitivity study.

### 4.5. Damage and Energy Dissipation Evaluation

Supporting functions are provided for damage-related and hysteretic response calculations.

These include utilities for:

- Ductility-based damage indices.
- General damage-index evaluation.
- Hysteretic force–displacement relationships.
- Dissipated energy estimation.
- Equivalent viscous damping ratio calculation.
- Effective structural property estimation.

The interpretation of each index depends on its mathematical definition, normalization, and the response quantities supplied to it. These functions should be reviewed before their outputs are compared across different models or studies.

### 4.6. Statistical Sensitivity Analysis

The project includes supporting statistical utilities for investigating the influence of model parameters on calculated responses.

The available tools include ANOVA-related sensitivity analysis, heatmaps, and plotting functions.

These methods can help identify response parameters associated with changes in column ductility. However, statistical association should not automatically be interpreted as physical causation.

## 5. Ground-Motion Data

The repository includes a separate utility:

`SEISMIC_GROUND_MOTION_ANALYSIS.py`

This script processes ground-motion files named according to the following pattern:

`Ground_Acceleration_1.txt`  
`Ground_Acceleration_2.txt`  
`...`  
`Ground_Acceleration_200.txt`

The actual number of files processed depends on the script configuration and the files available in the expected input directory.

The utility calculates peak absolute acceleration and descriptive statistics for the records it processes.

Before running the analysis, verify the following:

1. The ground-motion files are available in the expected directory.
2. Each file contains the acceleration data in the format expected by the script.
3. The acceleration units are identified correctly.
4. The sampling interval is consistent with the record metadata.
5. Any baseline correction, filtering, or resampling is documented.
6. The records are scaled consistently with the intended engineering objective.

**Important:** Peak ground acceleration alone does not fully characterize earthquake severity or spectral compatibility. For rigorous seismic comparisons, consider record duration, frequency content, response spectra, and the adopted record-selection and scaling criteria.

## 6. Repository Contents

| File | Purpose |
|---|---|
| `DUCT_SENSITIVITY_20_SEISMIC_MDOF.py` | Main ductility sensitivity analysis. |
| `SEISMIC_GROUND_MOTION_ANALYSIS.py` | Ground-motion data processing and descriptive statistics. |
| `ANALYSIS_FUNCTION.py` | Supporting analysis procedures. |
| `ANOVA_SENSITIVITY_FUN.py` | ANOVA-based sensitivity analysis utilities. |
| `BILINEAR_CURVE.py` | Bilinear curve processing utility. |
| `COMPUTE_EFFECTIVE_PROPERTIES_FUN_FREE_VIBRATION.py` | Effective-property calculations based on free-vibration response. |
| `COMPUTE_EFFECTIVE_PROPERTIES_FUN_PUSHOVER.py` | Effective-property calculations associated with pushover response. |
| `DAMAGE_INDEX_FUN.py` | Damage-index calculations. |
| `DAMPING_RATIO_FUN.py` | Damping-related calculations. |
| `DISSIPATED_ENERGY_WITH_PLOT_FUN.py` | Dissipated-energy calculation and visualization. |
| `DUCTILITY_DAMAGE_INDEX_FUN.py` | Ductility-based damage-index calculation. |
| `EIGENVALUE_ANALYSIS_FUN.py` | Eigenvalue and modal-analysis utilities. |
| `EQULIVALENT_VISCOUS_DAMPING_RATIO_FUN.py` | Equivalent viscous damping ratio calculation. |
| `FRAGILITY_CURVE_FUN.py` | Fragility-curve analysis utilities. |
| `OPENSEEES_HYSTERETICSM_FORCE_DISP_FUN.py` | Hysteretic force–displacement analysis. |
| `PERIOD_FUN.py` | Structural period-related calculations. |
| `PLOT_1D_SPRING.py` | One-dimensional spring visualization. |
| `PLOT_CONTOUR_3D_2D_FUN.py` | Two- and three-dimensional contour visualization. |
| `RAYLEIGH_DAMPING_FUN.py` | Rayleigh damping utilities. |
| `SALAR_MATH.py` | Supporting mathematical functions. |
| `SENSITIVITY_HEATMAP_FUN.py` | Sensitivity heatmap generation. |
| `SENSITIVITY_PLOTS_ANOVA.py` | Visualization of ANOVA sensitivity results. |
| `200_SEISMIC_TIME_HISTORY.rar` | Compressed ground-motion data archive. |
| `PDF_DUCT_SENSITIVITY_20_SEISMIC_MDOF.pdf` | Supplementary technical documentation. |
| `PPT_DUCT_SENSITIVITY_20_SEISMIC_MDOF.pptx` | Presentation material. |
| `COVER.png` | Project cover image. |
| `COVER_DISPLACEMENT_BASED_PUSHOVER.png` | Pushover-related cover image. |
| `SEISMIC_GROUND_MOTION_ANALYSIS.png` | Ground-motion analysis illustration. |
| `EQULIVALENT_VISCOUS_DAMPING_RATIO_FUN.png` | Equivalent viscous damping illustration. |

File availability and implementation details may change as the repository develops.

## 7. Requirements

The main analysis uses Python and OpenSeesPy, together with numerical, statistical, and visualization libraries.

Install the commonly required packages using:

```bash
pip install openseespy numpy pandas matplotlib scipy
```

Some functions may require additional packages depending on the current implementation. Review the imports in the relevant scripts if Python reports a missing dependency.

## 8. How to Run

Clone the repository:

```bash
git clone https://github.com/salardelavar/OPENSEES_SALAR.git
```

Navigate to the project directory:

```bash
cd OPENSEES_SALAR/EIGHT_ANALYSIS_DUCTILITY_DAMAGE_INDEX_\&_ENERGY_DISSIPATION_CAPACITY_INDEX/60_MDOF_10_FLOORS_04_COLUMNS_8_ANA/DUCT_SENSITIVITY_20_SEISMIC_MDOF
```

Run the main sensitivity analysis:

```bash
python DUCT_SENSITIVITY_20_SEISMIC_MDOF.py
```

Run the ground-motion preprocessing utility separately:

```bash
python SEISMIC_GROUND_MOTION_ANALYSIS.py
```

The scripts may depend on input files, relative paths, and supporting functions. Confirm that the required data files and modules are accessible from the execution environment before running them.

## 9. Engineering Validation and Reproducibility

For meaningful engineering interpretation, the following checks are recommended.

### Structural model verification

- Confirm element connectivity, boundary conditions, floor constraints, and mass distribution.
- Verify material models, section properties, and nonlinear constitutive parameters.
- Check the initial stiffness and natural periods against independent calculations where possible.
- Confirm that the applied loading and control-node displacement represent the intended analysis.

### Dynamic analysis verification

- Check the ground-motion time step, units, scaling factors, and input direction.
- Verify the damping formulation and the modes used to determine damping coefficients.
- Review convergence criteria, integration method, and analysis time step.
- Confirm that the applied record duration and analysis termination are appropriate.

### Response and damage calculations

- Check that displacement and force units are consistent.
- Confirm the definitions of yield displacement, ultimate displacement, and ductility.
- Verify the energy-dissipation calculation for the relevant loading protocol.
- Ensure that equivalent viscous damping calculations use a clearly defined energy input and reference elastic energy.
- Distinguish peak response quantities from cumulative response measures.

### Statistical analysis

- Document the number of earthquake records and ductility cases included in each calculation.
- Identify how missing, failed, or nonconvergent analyses are treated.
- Consider record-to-record variability when summarizing seismic response.
- Check statistical assumptions before interpreting ANOVA results.
- Avoid data leakage if machine-learning models are used for prediction or sensitivity assessment.

## 10. Interpretation and Limitations

The results describe the behavior of the implemented numerical model under the selected parameter values and ground-motion inputs.

They should not be interpreted as universal relationships between column ductility and seismic performance.

In particular:

- Numerical results depend on the structural configuration, material properties, damping model, and analysis settings.
- Ground-motion selection and scaling can substantially affect calculated response.
- Different damage indices may represent different physical mechanisms and are not necessarily interchangeable.
- Equivalent SDOF properties depend on the adopted transformation and normalization.
- Statistical sensitivity results depend on the sampled parameter range and the variability of the earthquake records.
- Numerical convergence does not, by itself, establish that a structural model is physically accurate.

Independent verification and validation are recommended before using the results for design decisions, structural assessment, or safety-critical applications.

## 11. Suggested Research Extensions

Potential extensions include:

1. Compare median, mean, and upper-percentile seismic responses across ground motions.
2. Investigate the sensitivity of results to ground-motion scaling and record selection.
3. Compare ductility-based damage indices with energy-based and deformation-based measures.
4. Examine the influence of damping assumptions on nonlinear response.
5. Evaluate alternative definitions of effective period and equivalent SDOF properties.
6. Quantify uncertainty in the estimated sensitivity relationships.
7. Compare ANOVA-based results with regression or other validated sensitivity-analysis methods.
8. Develop a reproducible benchmark dataset for comparing nonlinear structural response models.

## 12. References and Software

- **OpenSees:** Open System for Earthquake Engineering Simulation  
  https://opensees.berkeley.edu/

- **OpenSeesPy:** Python interface to OpenSees  
  https://openseespydoc.readthedocs.io/

- **NumPy:** Numerical computing in Python  
  https://numpy.org/

- **SciPy:** Scientific computing and numerical methods  
  https://scipy.org/

- **Pandas:** Data analysis and data structures  
  https://pandas.pydata.org/

- **Matplotlib:** Scientific plotting and visualization  
  https://matplotlib.org/

## 13. Author

**Salar Delavar Ghashghaei (Qashqai)**  
Civil and Structural Engineering | OpenSeesPy | Nonlinear Structural Analysis | Earthquake Engineering

GitHub: https://github.com/salardelavar

## 14. Disclaimer

This repository is intended for educational, research, and engineering-development purposes. The scripts and numerical results are provided as-is and may contain errors, simplifications, or implementation limitations.

Users are responsible for independently verifying the theoretical formulations, source code, input data, numerical results, and applicability of the methods before using them in research publications, engineering assessments, or professional decisions.

**The availability of a numerical result does not establish the validity of the underlying structural model or its suitability for engineering practice.**

