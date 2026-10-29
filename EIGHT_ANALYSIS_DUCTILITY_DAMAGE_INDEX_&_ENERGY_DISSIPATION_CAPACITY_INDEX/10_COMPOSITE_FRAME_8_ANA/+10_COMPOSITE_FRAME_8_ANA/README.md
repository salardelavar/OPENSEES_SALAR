# COMPREHENSIVE NONLINEAR SEISMIC ASSESSMENT OF A COMPOSITE SECTION FRAME : AN OPENSEES FRAMEWORK FOR STATIC PUSHOVER, CYCLIC DEGRADATION, STATIC TIME-HISTORY AND DYNAMIC TIME-HISTORY ANALYSIS

# ارزیابی جامع غیرخطی سازه با المان مرکب فلزی بتنی : چارچوبی مبتنی بر اوپنسیس در پایتون برای تحلیل پوش‌آور استاتیکی، تحلیل چرخه‌ای همراه با کاهش سختی، تحلیل تاریخچه زمانی استاتیکی و دینامیکی، با در نظر گرفتن شاخص‌های آسیب شکل‌پذیری در المان و سازه و ارزیابی شاخص ظرفیت انرژی اتلاف‌شده


![alt text](https://github.com/salardelavar/OPENSEES_SALAR/blob/main/EIGHT_ANALYSIS_DUCTILITY_DAMAGE_INDEX_%26_ENERGY_DISSIPATION_CAPACITY_INDEX/10_COMPOSITE_FRAME_8_ANA/+10_COMPOSITE_FRAME_8_ANA/COVER.png) 



## 1. Overview

This directory belongs to the **EIGHT_ANALYSIS_DUCTILITY_DAMAGE_INDEX_&_ENERGY_DISSIPATION_CAPACITY_INDEX** collection within the [OPENSEES_SALAR](https://github.com/salardelavar/OPENSEES_SALAR) repository.

It provides a complete, research-grade computational framework for the nonlinear seismic assessment of **composite steel-concrete frames** using **OpenSeesPy**. The models combine structural steel and reinforced (or FRP-confined) concrete to capture the interaction of strength, stiffness, ductility, and energy dissipation under monotonic, cyclic, and dynamic loading.

The implementation supports **eight complementary analysis protocols** that enable:

- Evaluation of fundamental periods and modal characteristics
- Establishment of the gravity-load state
- Generation of capacity curves and plastic mechanism identification (pushover)
- Quantification of hysteretic degradation, pinching, and cyclic energy dissipation
- Static and dynamic response under time-dependent external forces
- Free-vibration response after initial conditions
- Computation of **Ductility Damage Index (DDI)** at member and system levels
- Computation of **Energy Dissipation Capacity Index (EDCI)**

All results are obtained from rigorously defined fiber-section composite members, consistent constitutive models, and carefully controlled numerical procedures.

> **Important**: The presence of an analysis category in this documentation does not guarantee that every possible option is activated in every script. Always inspect the specific Python file for the exact constitutive models, element formulations, recorders, and post-processing routines that are implemented.

---

## 2. Engineering Objectives

1. Construct high-fidelity nonlinear models of composite frames using fiber-section formulations.
2. Quantify the transition from elastic to inelastic behavior under realistic loading histories.
3. Generate complete capacity curves and identify plastic mechanisms via displacement-controlled pushover analysis.
4. Capture cyclic degradation, pinching, and residual deformation under symmetric cyclic displacement protocols.
5. Evaluate structural response under harmonic or exponentially decaying time-dependent loading (static and dynamic).
6. Compute system- and member-level ductility demand and damage indices.
7. Quantify hysteretic energy dissipation and derive an Energy Dissipation Capacity Index.
8. Provide a verified, reproducible baseline for parametric studies, sensitivity analysis, and fragility assessment.

---

## 3. The Eight Analysis Protocols

| # | Protocol | Description |
|---|----------|-------------|
| 1 | **PERIOD** | Eigenvalue analysis – fundamental periods, mode shapes, and modal participation |
| 2 | **STATIC (Gravity)** | Nonlinear static analysis under self-weight and permanent loads to establish the initial stress state |
| 3 | **PUSHOVER** | Displacement-controlled monotonic pushover – capacity curve, yield displacement, ultimate displacement, plastic mechanism |
| 4 | **CYCLIC_DISPLACEMENT** | Symmetric cyclic displacement protocol – hysteresis loops, stiffness/strength degradation, pinching, cumulative energy dissipation |
| 5 | **STATIC_EXTERNAL_TIME-DEPENDENT_LOADING** | Quasi-static analysis under \( P(t) = P_0 \sin(\omega t) \) or \( P(t) = P_0 e^{-0.05\omega t}\sin(\omega t) \) |
| 6 | **DYNAMIC_EXTERNAL_TIME-DEPENDENT_LOADING** | Nonlinear dynamic analysis under the same time-dependent external forces |
| 7 | **FREE_VIBRATION** | Free vibration after initial displacement/velocity (or after a loading phase) |
| 8 | **DAMAGE & ENERGY INDICES** | Post-processing for Ductility Damage Index and Energy Dissipation Capacity Index |

---

## 4. Composite Section Library

The framework includes a rich set of composite fiber-section generators (both 2-D and 3-D visualizations):

- Rectangular and circular concrete cores with steel plates / I-sections
- Confined concrete (stirrup or FRP confinement)
- Composite I-sections with different flange widths/heights
- Double-I and specialized steel sections
- Composite shear-wall sections
- FRP-confined circular sections

Key scripts:

- `COMPOSITE_REC_I_SECTION_FUN.py`
- `COMPOSITE_CIRCULAR_CONFINED_CONCRETE_SECTION_FUN_EXTRA.py` (+ 3-D versions)
- `COMPOSITE_FRP_*` series
- `COMPOSITE_SECTION_DIFF_WIDTH_HEIGHT_FUN_EXTRA_TWO.py`
- `DOUBLE_I_SECTION_FIBER_FUN.py` / `DOUBLE_I_SECTION_QUAD_FUN.py`
- `CONFINED_CONCRETE_WITH_PLATE_SECTION_FUN.py`
- `FLANGE_DIFFERENT_I_STEEL_SECTION_FUN.py`

All sections are defined with fiber discretization suitable for accurate moment-curvature and axial-flexure interaction response.

---

## 5. Core Computational Components

| File | Purpose |
|------|---------|
| `P(t)_COMPOSITE_SECTION_FRAME_8_ANA.py` | Main analysis driver implementing the eight protocols |
| `ANALYSIS_FUNCTION.py` | Reusable analysis control, solvers, and convergence utilities |
| `BILINEAR_CURVE.py` | Idealized bilinear approximation of capacity curves |
| `PERIOD_FUN.py` / `EIGENVALUE_ANALYSIS_FUN.py` | Modal analysis utilities |
| `DAMPING_RATIO_FUN.py` / `RAYLEIGH_DAMPING_FUN.py` | Rayleigh and equivalent viscous damping |
| `EQULIVALENT_VISCOUS_DAMPING_RATIO_FUN.py` | Calculation of \(\xi_{eq}\) from hysteretic loops |
| `FRAGILITY_FUN.py` | Fragility-related post-processing helpers |
| `PLOT_2D.py` | Standardized 2-D response plotting |
| `Ground_Acceleration_X.txt` / `Y.txt` | Sample ground-motion records (when used) |

Supporting documentation:

- `PDF_P(t)_COMPOSITE_SECTION_FRAME_8_ANA.pdf`
- `PPT_P(t)_COMPOSITE_SECTION_FRAME_8_ANA.pptx`
- `COVER.png` and numerous section visualization figures

---

## 6. Ductility Damage Index (DDI)

The Ductility Damage Index is evaluated at both **member** and **system** levels.  

A generic normalized form is:

\[
DDI = \frac{\mu - 1}{\mu_u - 1}
\]

or more advanced cumulative formulations that incorporate cyclic plastic deformation demand.

**Users must document**:

- Definition of yield displacement / curvature (\(\Delta_y\) or \(\phi_y\))
- Ultimate capacity used for normalization
- Whether the index is peak-based or cumulative
- Member-level vs. global (roof drift) evaluation

The index is a **model-dependent performance indicator**, not a direct probability of collapse.

---

## 7. Energy Dissipation Capacity Index (EDCI)

Hysteretic energy is computed via trapezoidal integration of the force-deformation response:

\[
\Delta W_i \approx \frac{F_i + F_{i-1}}{2}(u_i - u_{i-1}), \qquad
W = \sum_{i=1}^{n} \Delta W_i
\]

The Energy Dissipation Capacity Index is typically expressed as a ratio of dissipated energy demand to a reference energy capacity (from cyclic or pushover analysis).  

**Critical documentation requirements**:

- Exact integration limits and sign conventions
- Treatment of elastic recoverable energy
- Definition of the reference capacity
- Consistency between the loading protocol used for demand and for capacity

Equivalent viscous damping is obtained from:

\[
\xi_{eq} = \frac{E_D}{4\pi E_S}
\]

where \(E_D\) is the area of the hysteresis loop and \(E_S\) is the maximum recoverable strain energy (definition must match the script implementation).

---

## 8. Software Requirements

```bash
Python ≥ 3.8
openseespy
numpy
scipy
matplotlib
pandas          # recommended for post-processing
