# OPENSEES_CABLE_SUSPENSION_BRIDGE_04

**Simplified Numerical Modeling and Analysis of the Presidente Ibáñez Cable Suspension Bridge using Python & OpenSeesPy**

![Cover](COVER.png)

---

## Overview

This repository provides a complete, research-oriented Python implementation for the **2D nonlinear finite-element analysis** of a cable suspension bridge inspired by the historic **Puente Presidente Ibáñez** (Aysén Region, Chile).

The bridge is Chile’s longest suspension bridge (main span ≈ 210 m). It features two steel towers (~25 m high) supporting eight main steel cables on each side, vertical hangers, and a concrete deck stiffened by longitudinal and transverse beams.

The model deliberately simplifies the originally truss-shaped deck into equivalent truss elements while retaining the essential cable–hanger–deck interaction, geometric nonlinearity (corotational formulation), and optional material nonlinearity. It is designed for educational, research, and preliminary design studies in structural and earthquake engineering.

**Author:** Salar Delavar Ghashghaei (Qashqai)  
**Email:** salar.d.ghashghaei@gmail.com  
**Repository path:** [`CABLE_SUSPENSION_BRIDGE/OPENSEES_CABLE_SUSPENSION_BRIDGE_04`](https://github.com/salardelavar/OPENSEES_SALAR/tree/main/CABLE_SUSPENSION_BRIDGE/OPENSEES_CABLE_SUSPENSION_BRIDGE_04)

---

## Bridge Description (Reference Structure)

| Parameter                    | Value                          |
|-----------------------------|--------------------------------|
| Main span                   | 210 m                          |
| Tower height                | ≈ 25 m                         |
| Main cables                 | 8 cables per side              |
| Hangers                     | 22 hangers per side            |
| Deck width                  | ≈ 7 m                          |
| Year of construction        | 1961–1966 (Krupp Rheinhausen)  |
| Status                      | National Monument of Chile     |

**Wikipedia (Spanish):** [Puente Presidente Ibáñez](https://es.wikipedia.org/wiki/Puente_Presidente_Ib%C3%A1%C3%B1ez)

---

## Key Features of the Model

- **2D plane-frame idealization** (`-ndm 2 -ndf 2`)
- **Corotational truss elements** (`corotTruss`) for large-displacement geometric nonlinearity
- Three distinct cable groups:
  - Top main cables
  - Bottom/deck longitudinal elements
  - Vertical hangers
- Optional **linear elastic** or **nonlinear hysteretic** material models
- Mass lumped at deck nodes
- Support conditions: fixed at both ends of cables and deck (simply-supported idealization)
- Comprehensive post-processing and visualization

---

## Analysis Capabilities

| Analysis Type              | Description                                                                 | Key Outputs |
|---------------------------|-----------------------------------------------------------------------------|-------------|
| **Pushover**              | Displacement-controlled vertical loading at mid-span                       | Capacity curve, deformed shapes, base reactions |
| **Free Vibration**        | Initial displacement / velocity followed by free oscillation               | Time histories (disp, vel, accel), damping estimation |
| **Eigenvalue Analysis**   | Modal analysis (natural frequencies & mode shapes)                         | Periods, frequencies, mode shapes |
| **Dynamic / Seismic**     | Time-history analysis under ground accelerations (X & Y components)        | Response histories, base shears |

Supporting modules handle Rayleigh damping, bilinear curve fitting, damping-ratio estimation, and robust convergence algorithms.

---

## Repository Structure
OPENSEES_CABLE_SUSPENSION_BRIDGE_04/
├── OPENSEES_CABLE_SUSPENSION_BRIDGE_04.py   # Main driver script
├── ANALYSIS_FUNCTION.py                     # Robust convergence helper
├── BILINEAR_CURVE.py                        # Bilinear approximation utilities
├── DAMPING_RATIO_FUN.py                     # Damping ratio calculation
├── EIGENVALUE_ANALYSIS_FUN.py               # Modal analysis
├── RAYLEIGH_DAMPING_FUN.py                  # Rayleigh damping coefficients
├── PLOT_2D.py                               # 2-D plotting helpers
├── SALAR_MATH.py                            # Mathematical utilities
├── Ground_Acceleration_X.txt                # Ground motion (X)
├── Ground_Acceleration_Y.txt                # Ground motion (Y)
├── COVER.png
├── OPENSEES_CABLE_SUSPENSION_BRIDGE_04.png
├── PDF_OPENSEES_CABLE_SUSPENSION_BRIDGE_04.pdf
├── PPT_OPENSEES_CABLE_SUSPENSION_BRIDGE_04.pptx
└── README.md


---

## Model Parameters (Default Values)

```python
L            = 210000.0      # [mm]  Span length
H1           = 25000.0       # [mm]  Height of top cable
arc_depth    = 20000.0       # [mm]  Sag of main cable
E_cable      = 210e5         # [N/mm²] Modulus of elasticity
Cable_Dia_01 = 145600        # [mm]  Equivalent diameter – bottom/deck
Cable_Dia_02 = 300           # [mm]  Equivalent diameter – top cables
Cable_Dia_03 = 50            # [mm]  Equivalent diameter – hangers
num_nodes    = 22            # Number of nodes along the span


![alt text](https://github.com/salardelavar/OPENSEES_SALAR/blob/main/CABLE_SUSPENSION_BRIDGE/OPENSEES_CABLE_SUSPENSION_BRIDGE_04/COVER.png) 

THIS PYTHON SCRIPT WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
