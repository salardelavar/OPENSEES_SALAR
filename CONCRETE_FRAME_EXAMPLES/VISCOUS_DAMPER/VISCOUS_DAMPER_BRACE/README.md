# Viscous Damper Brace – Reinforced Concrete Frame Examples

**Nonlinear Static & Dynamic Analysis of a 2D RC Frame with Viscous Damper Brace using OpenSeesPy**

**Author:** Salar Delavar Ghashghaei (Qashqai)  
**Email:** salar.d.ghashghaei@gmail.com  

---

## Overview

This directory contains a set of educational and research-oriented OpenSeesPy examples that demonstrate the seismic and dynamic performance of a single-bay reinforced concrete frame equipped with a **diagonal viscous damper brace**.

The examples cover three complementary analysis cases:

| Folder | Description |
|--------|-------------|
| `VISCOUS_DAMPER_FREE_VIBRATION` | Free-vibration analysis (elastic vs inelastic) with viscous damper |
| `VISCOUS_DAMPER_IMPACT_LOAD` | Harmonic / impact loading analysis with viscous damper |
| `VISCOUS_DAMPER_SEISMIC` | Nonlinear time-history analysis under ground motion records |

All scripts use the same base RC frame geometry and material models so that the effect of the viscous damper can be clearly compared across different loading scenarios.

---

## Structural Model

- **Geometry**: Single-bay, single-story 2D RC frame  
  - Column length = 3000 mm  
  - Beam length = 7000 mm  
- **Sections**:
  - Columns: 500 × 500 mm  
  - Beam: 500 × 300 mm  
- **Materials**:
  - Concrete: `Concrete01` (confined core + unconfined cover)  
  - Steel: bilinear or hysteretic models (with/without hardening)  
- **Elements**: `nonlinearBeamColumn` (distributed plasticity) or elastic beam-column  
- **Geometric nonlinearity**: Corotational transformation  
- **Damper**: Modeled as a truss element with viscous damping properties connecting the base of one column to the top of the opposite column (diagonal brace configuration)

### Viscous Damper Modeling Notes

The damper is implemented via a custom truss element that can represent both linear and nonlinear viscous behavior:

\[
F = C \cdot |v|^{\alpha-1} \cdot v
\]

Typical parameters used in the examples:
- Damping ratio target for damper ≈ 15 %
- Section area of damper element = 35² mm²

---

## Analysis Types

### 1. Free Vibration (`VISCOUS_DAMPER_FREE_VIBRATION`)
- Initial displacement / velocity / acceleration applied
- Comparison of elastic vs inelastic frame response
- Estimation of damping ratio via logarithmic decrement
- Period elongation and stiffness degradation tracking

### 2. Harmonic / Impact Loading (`VISCOUS_DAMPER_IMPACT_LOAD`)
- Sinusoidal force applied at the top of the frame
- Study of resonance and energy dissipation under sustained dynamic loading
- Force–displacement hysteresis of the damper and frame

### 3. Seismic Analysis (`VISCOUS_DAMPER_SEISMIC`)
- Uniform excitation using ground acceleration records (X and/or Y)
- Pushover analysis for capacity curve and ductility assessment
- Ductility Damage Index (DDI)
- Overstrength factor (Ω₀) and response modification factor (R)
- Rayleigh damping calibrated from eigenvalue analysis

---

## Key Features

- Comparative elastic vs inelastic response
- Variable confinement enhancement ratios in columns and beams
- Rayleigh damping + viscous damper damping
- Eigenvalue analysis and period elongation monitoring
- Bilinear idealization of pushover curves
- Damage index calculation
- Comprehensive plotting (time histories, hysteresis loops, deformed shapes, stiffness degradation)
- Results exported to Excel

---

## Requirements

```bash
pip install openseespy numpy matplotlib pandas
