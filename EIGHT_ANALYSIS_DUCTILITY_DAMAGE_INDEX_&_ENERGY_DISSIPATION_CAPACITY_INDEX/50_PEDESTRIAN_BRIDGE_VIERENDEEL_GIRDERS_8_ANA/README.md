# COMPREHENSIVE NONLINEAR SEISMIC ASSESSMENT OF A PEDESTRIAN BRIDGE STEEL VIERENDEEL GIRDERS AN OPENSEES FRAMEWORK FOR STATIC PUSHOVER, CYCLIC DEGRADATION, STATIC TIME-HISTORY AND DYNAMIC TIME-HISTORY ANALYSIS

![alt text](https://github.com/salardelavar/OPENSEES_SALAR/blob/main/EIGHT_ANALYSIS_DUCTILITY_DAMAGE_INDEX_%26_ENERGY_DISSIPATION_CAPACITY_INDEX/50_PEDESTRIAN_BRIDGE_VIERENDEEL_GIRDERS_8_ANA/COVER.png) 

# Pedestrian Bridge – Steel Vierendeel Girders  
### Comprehensive Nonlinear Seismic Assessment Framework in OpenSeesPy

**An advanced OpenSeesPy framework for static pushover, cyclic degradation, static & dynamic time-history analysis, free vibration, and seismic evaluation of a pedestrian bridge with steel Vierendeel girders.**

**Author:** Salar Delavar Ghashghaei (Qashqai)  
**Email:** salar.d.ghashghaei@gmail.com  
**Repository Path:**  
`EIGHT_ANALYSIS_DUCTILITY_DAMAGE_INDEX_&_ENERGY_DISSIPATION_CAPACITY_INDEX / 50_PEDESTRIAN_BRIDGE_VIERENDEEL_GIRDERS_8_ANA`

---

## Overview

This project implements a high-fidelity **2D nonlinear finite-element model** of a pedestrian bridge composed of **steel Vierendeel girders**. The framework systematically executes **eight complementary analysis protocols** and quantifies structural performance through advanced indices:

- **Ductility Damage Index (DDI)**
- **Energy Dissipation Capacity Index (EDCI)**
- **Equivalent Viscous Damping Ratio (EVDR)**

The code supports both material nonlinearity (elastic-perfectly plastic or hysteretic steel models) and geometric nonlinearity (large displacements, P-Δ effects, and intentional deck camber).

It is designed for **performance-based earthquake engineering (PBEE)** applications, educational use, research on Vierendeel systems, and comparative studies of energy dissipation mechanisms in bridge structures.

---

## Structural Model

| Feature                        | Description                                      |
|--------------------------------|--------------------------------------------------|
| Structure type                 | Pedestrian bridge with steel Vierendeel girders |
| Dimensionality                 | 2D (`-ndm 2 -ndf 3`)                             |
| Geometry                       | Top & bottom chords + vertical web members       |
| Camber                         | Sinusoidal deck camber (`ε = 0.0033 × Length`)   |
| Nodes                          | 60 nodes (30 bottom + 30 top)                    |
| Elements                       | `forceBeamColumn` or `dispBeamColumn` with fiber sections |
| Sections                       | Custom I-shaped steel plate sections (variable width/height) |
| Material models                | Elastic-Perfectly Plastic / Hysteretic / Steel02 |
| Boundary conditions            | Simply supported (pinned–roller)                 |
| Mass distribution              | Consistent or lumped mass                        |

---

## Eight Analysis Protocols

The script executes the following standardized analysis suite:

| #  | Protocol                              | Description                                                                 | Key Outputs |
|----|---------------------------------------|-----------------------------------------------------------------------------|-------------|
| 1  | **PERIOD**                            | Eigenvalue analysis & modal properties                                      | Periods, mode shapes, modal participation |
| 2  | **STATIC**                            | Gravity / dead-load analysis                                                | Base reactions, vertical displacements, internal forces |
| 3  | **PUSHOVER**                          | Displacement-controlled nonlinear static analysis                           | Capacity curve, yield & ultimate points, plastic mechanism |
| 4  | **CYCLIC_DISPLACEMENT**               | Symmetric cyclic displacement protocol                                      | Full hysteresis loops, pinching, strength/stiffness degradation, cyclic energy |
| 5  | **STATIC_EXTERNAL_TIME-DEPENDENT**    | Static analysis under time-varying load \( P(t) \)                          | Response under harmonic or damped-harmonic loading |
| 6  | **DYNAMIC_EXTERNAL_TIME-DEPENDENT**   | Dynamic analysis under the same external loads                              | Transient response, resonance effects |
| 7  | **FREE-VIBRATION**                    | Free vibration with initial conditions                                      | Logarithmic decrement, damping ratio estimation |
| 8  | **SEISMIC**                           | Nonlinear time-history analysis under multi-directional ground motions      | Displacement/acceleration histories, energy demand, damage indices |

**Loading functions available:**
\[
P(t) = P_0 \sin(\omega t)
\]
\[
P(t) = P_0 e^{-0.05\omega t} \sin(\omega t)
\]

Rayleigh damping (typically 3 % critical) is applied for dynamic analyses.

---

## Performance Metrics

### 1. Ductility Damage Index (DDI)
\[
DDI = \frac{D_d - D_y}{D_u - D_y}
\]
- \( D_d \): Peak displacement demand  
- \( D_y \): Yield displacement (from bilinear idealization of pushover curve)  
- \( D_u \): Ultimate displacement  

### 2. Energy Dissipation Capacity Index (EDCI)
\[
EDCI = 100 \times \frac{E_d^{\text{earthquake}}}{E_d^{\text{cyclic}}}
\]
Compares the energy dissipated under seismic excitation with the maximum energy dissipation capacity obtained from cyclic testing.

### 3. Equivalent Viscous Damping Ratio (EVDR)
\[
\xi_{eq} = 100 \times \frac{E_d}{4\pi E_s}
\]
- \( E_d \): Dissipated energy (area of hysteresis loop)  
- \( E_s \): Strain energy (area under the backbone curve)

---

## Key Features

- Automatic bilinear idealization of capacity curves
- Fiber-section modeling of I-shaped steel plates with different flange/web dimensions
- Section curvature and stress–strain monitoring
- Comprehensive post-processing and publication-quality plotting
- Modal properties report generation
- Fragility-related helper functions
- Support for both force-based and displacement-based beam-column elements
- Robust convergence handling with adaptive algorithms

---

## File Structure



THIS PYTHON SCRIPT IS WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
