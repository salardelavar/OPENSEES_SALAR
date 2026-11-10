# NUMERICAL NEWTON-RAPHSON VERIFICATION OF THE OPENSEES CONTACT-DRIVEN SDOF MONOTINIC PUSHOVER MODEL

This Python script implements a Numerical Newton-Raphson solver to independently verify an OpenSees contact-driven single-degree-of-freedom (SDOF) monotonic pushover analysis. It defines a symmetric piecewise-linear force-displacement backbone, assembles internal forces and tangent stiffness for two parallel springs (one always active, the second activated after a gap closes), and performs displacement-controlled pushover analysis under two interpretations of the contact spring: an “OPENSEES” mode that matches the original model’s strain definition and a more physical “gap” mode. Force-controlled Newton-Raphson is used as an inverse check, periods are computed from the tangent stiffness, and results (pushover curves and period evolution) are compared with an optional OpenSeesPy run that dynamically adds the second zeroLength element once the gap is closed.

 ![alt text](COVER_SDOF_CONTACT_NUMERICAL_OPENSEESPY.jpg)

