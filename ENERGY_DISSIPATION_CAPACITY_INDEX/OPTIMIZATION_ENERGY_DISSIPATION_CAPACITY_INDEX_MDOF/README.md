OPTIMIZATION PROBLEM FOR ENERGY DISSIPATION CAPACITY INDEX (EDCI)
-----------------------------------------------------------------
1. Problem Formulation – We invert the conventional design process:
   instead of checking a given ultimate displacement DSU, we find the 'minimum'
   DSU that makes the structure's hysteretic energy dissipation under
   a specific earthquake exactly equal to a target fraction (20%) of
   its maximum dissipation capacity measured from a full cyclic test.

2. Energy‑Based Performance Metric – The Energy Dissipation Capacity
   Index (EDCI) = E_seismic / E_cyclic is a scalar that quantifies
   how close the seismic demand comes to saturating the structure's
   plastic energy absorption capability – a direct proxy for cumulative
   damage and collapse margin.

3. Nonlinear Forward Model – Each evaluation of DSU requires solving
   the equation of motion: M·ü + C·u̇ + R(u, DSU) = f_ext(t) with a
   rate‑independent hysteretic material (elastic‑perfectly‑plastic or
   with hardening), where the restoring force R depends sensitively on
   DSU through the yield surface.

4. Root‑Finding via Newton–Raphson – The residual g(DSU) = EDCI(DSU)-0.20
   is driven to zero using a secant‑like Newton iteration, where the
   Jacobian g'(DSU) is approximated by central finite differences
   [g(DSU+ε) - g(DSU-ε)]/(2ε), with ε tuned to balance truncation and
   conditioning errors.

5. Numerical Challenges – The function g is non‑smooth due to abrupt
   yielding, pinching, and stiffness degradation; hence the solver is
   augmented with a line‑search and a fallback bisection step when the
   Newton update overshoots, ensuring robust convergence despite
   discontinuities.

6. Cyclic Capacity Envelope – The denominator E_cyclic is obtained from
   a quasi‑static displacement‑controlled cyclic protocol that
   progressively increases amplitude to capture the full hysteresis
   loop, including the post‑peak softening branch – defining the
   structure's ultimate energy dissipation reservoir.

7. Seismic Demand Evaluation – The numerator E_seismic is computed from
   a full dynamic time‑history analysis with Rayleigh damping (3%) and
   a scaled ground motion, integrating the instantaneous power
   ∫ R·du over the entire shaking duration.

8. Physical Insight – The converged DSU corresponds to the strength
   threshold that allows the structure to exploit its full ductility
   without exceeding the damage index associated with the target EDCI;
   below this value, the system would enter the "failure" zone
   (EDCI > 100%) under the same excitation.

9. Sensitivity and Stability – The finite‑difference step ε is set to
   1e-3·DSU to avoid numerical ill‑conditioning, and the solver terminates
   when the relative change in DSU falls below 1e-6, ensuring that the
   design strength is determined with sub‑percent accuracy.

10. Engineering Decision – This optimization yields a rational,
    risk‑consistent design point that directly links seismic hazard
    (energy input) to structural capacity (energy absorption), enabling
    performance‑based earthquake engineering without empirical
    overstrength factors – a move towards collapse‑prevention
    verification based on explicit energy balance.

On This Page, There Are Python and OpenSees codes, which are written by Salar Delavar Ghashghaei (Qashqai).
Please note that the content may not be entirely free of errors or inaccuracies.
