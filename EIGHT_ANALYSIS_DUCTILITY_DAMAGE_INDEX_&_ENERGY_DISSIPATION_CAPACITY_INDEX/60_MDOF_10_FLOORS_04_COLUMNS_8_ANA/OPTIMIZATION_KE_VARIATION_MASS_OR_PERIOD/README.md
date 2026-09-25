# OPTIMIZATION_KE_VARIATION_MASS_OR_PERIOD

**Optimum Elastic Stiffness (Ke) of Columns in a 10-Story MDOF Frame**  
Using Finite-Difference Newton–Raphson Iteration and OpenSees

Author: **Salar Delavar Ghashghaei (Qashqai)**  
Email: salar.d.ghashghaei@gmail.com

---

## Overview

This folder contains two complementary optimization studies for a **10-story, 4-column-per-floor** multi-degree-of-freedom (MDOF) frame.

In both cases the goal is the same:  
**find the optimum initial elastic stiffness (Ke) of the columns** so that the structural period matches a prescribed target.

| Subfolder | Mass | Target Period | What is optimized |
|-----------|------|---------------|-------------------|
| `OPTIMIZATION_KE_MASS_VARITIONS` | **Varied** | **Constant** | Ke |
| `OPTIMIZATION_KE_PERIOD_VARITIONS` | **Constant** | **Varied** | Ke |

The optimization uses a **finite-difference Newton–Raphson** solver.

---

## Structure Model

- 10 floors  
- 4 columns per floor  
- Beams with infinite flexural rigidity (`EI = ∞`)  
- Elastic or inelastic material models  
- Equivalent SDOF system obtained via displacement-based design (pushover)

---

## Key Features

- Newton–Raphson optimization of column elastic stiffness  
- Structural period analysis (`ANAL_TYPE = 'PERIOD'`)  
- Ductility Damage Index  
- Energy Dissipation Capacity Index (EDCI)  
- Equivalent Viscous Damping Ratio  
- Multiple analysis protocols (period, pushover, cyclic, free-vibration, seismic, …)

---

## How to Run

### 1. Constant Period – Vary Mass
```bash
cd OPTIMIZATION_KE_MASS_VARITIONS
python OPTIMIZATION_KE_MASS_VARITIONS.py

### 2. Constant Mass – Vary Period
```bash
cd OPTIMIZATION_KE_PERIOD_VARITIONS
python OPTIMIZATION_KE_PERIOD_VARITIONS.py
