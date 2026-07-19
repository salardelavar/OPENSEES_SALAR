COMPREHENSIVE NONLINEAR SEISMIC ASSESSMENT OF MULTI-MODE POST-BUCKLING PHENOMENA IN A POST-FIRE STEEL COLUMN WITH SEMI-RIGID CONNECTIONS: AN OPENSEES FRAMEWORK FOR STATIC PUSHOVER, CYCLIC DEGRADATION, STATIC TIME-HISTORY, AND DYNAMIC TIME-HISTORY ANALYSES

Assume we have a structure where, after a fire, certain regions of the structural members have lost their initial material properties (yield strength and modulus of elasticity), and the entire column member has undergone lateral buckling. Now, we assess the strength of this structure by performing eight different analyses.

In Each OpenSees Project Script File, Eight Analysis Protocols Are Implemented:    
(1) [PERIOD] : Structural Period
(2) [STATIC] : Gravity load analysis establishing dead load state
(3) [PUSHOVER] : Displacement-controlled pushover generating full capacity curves
 and plastic mechanism identification
(4) [CYCLIC_DISPLACEMENT] : Symmetric cyclic displacement protocol capturing hysteresis,
 pinching behavior, and energy dissipation degradation
(5) [STATIC_EXTERNAL_TIME-DEPENDENT_LOADING] : Static Analysis of External time-dependent loading P(t) = P0 sin(wt) or P(t) = P0 exp(-0.05wt) sin(wt)  
(6) [DYNAMIC_EXTERNAL_TIME-DEPENDENT_LOADING] : Dynamic Analysis of External time-dependent loading P(t) = P0 sin(wt) or P(t) = P0 exp(-0.05wt) sin(wt)  
(7) [FREE-VIBRATION] : Free-vibration with initial conditions extracting damping ratios
 via logarithmic decrement
(8) [SEISMIC] : Multi-directional seismic excitation with Rayleigh damping (3% ratio)

On This Page, There Are Python and OpenSees Python Scripts, which are written by Salar Delavar Ghashghaei (Qashqai).
Please note that the content may not be entirely free of errors or inaccuracies.