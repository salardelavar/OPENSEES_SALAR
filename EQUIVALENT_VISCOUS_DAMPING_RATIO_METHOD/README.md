The equivalent viscous damping ratio is a linearization concept that quantifies the energy dissipated
 by a nonlinear damping or hysteretic mechanism as if it were produced by an ideal viscous damper,
 enabling simplified linear dynamic analysis. It is derived by equating the actual energy lost per
 cycle of steady-state harmonic motion—arising from sources like friction, material yielding, or
 aerodynamic drag to the energy that a linear viscous damper would dissipate over the same cycle
 at a given amplitude and frequency. This yields an amplitude- and frequency-dependent damping ratio (ζ_eq),
 often expressed as ζ_eq = (E_dissipated) / (4π E_stored) for a single-degree-of-freedom system.
 The resulting value allows engineers to approximate a nonlinear system's damping behavior within a
 linear framework, such as modal analysis, response spectrum methods, or equivalent linearization of
 inelastic structures, where it is iteratively updated to match the expected response level.
 

In Each OpenSees Project Script File, Five Analysis Protocols Are Implemented:    
(1) [CYCLIC_DISPLACEMENT] : Symmetric cyclic displacement protocol capturing hysteresis,
 pinching behavior, and energy dissipation degradation
(2) [STATIC_EXTERNAL_TIME-DEPENDENT_LOADING] : Static Analysis of External time-dependent loading P(t) = P0 sin(wt) or P(t) = P0 exp(-0.05wt) sin(wt)  
(3) [DYNAMIC_EXTERNAL_TIME-DEPENDENT_LOADING] : Dynamic Analysis of External time-dependent loading P(t) = P0 sin(wt) or P(t) = P0 exp(-0.05wt) sin(wt)  
(4) [FREE-VIBRATION] : Free-vibration with initial conditions extracting damping ratios
 via logarithmic decrement
(5) [SEISMIC] : Multi-directional seismic excitation with Rayleigh damping (3% ratio)

On This Page, There Are Python and OpenSees codes, which are written by Salar Delavar Ghashghaei (Qashqai).
Please note that the content may not be entirely free of errors or inaccuracies.
