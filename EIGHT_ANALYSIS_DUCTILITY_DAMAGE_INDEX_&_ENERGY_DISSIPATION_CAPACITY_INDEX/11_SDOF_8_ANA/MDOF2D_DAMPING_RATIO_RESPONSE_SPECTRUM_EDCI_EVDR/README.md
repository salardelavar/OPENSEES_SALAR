# PARAMETRIC NONLINEAR SEISMIC DEMAND SPECTRA OF A 2-DOF (X & Y SPRINGS) SYSTEM UNDER 2 GROUND MOTIONS WITH VARIABLE DAMPING RATIO CAPACITY

This script builds a 2-DOF OpenSeesPy structure (`-ndm 2 -ndf 2`) made of a fixed base node and a mass node, linked by one `zeroLength` element that holds an independent hysteretic spring and viscous damper in X and in Y, each with its own yield force, stiffness and ultimate displacement. It first runs a cyclic displacement protocol on each spring to get a baseline dissipated energy (the area of the force-displacement loop, computed with the shoelace formula). It then sweeps damping ratio (1% to 10%) and mass, which changes the period, and for each combination it subjects the structure to two ground motions applied at the same time, one along X and one along Y, using a nonlinear Newmark time-history analysis. For each direction it extracts the peak displacement, velocity, acceleration and base reaction, a ductility damage index (how far the displacement has gone from yield toward ultimate), the stiffness and the damping ratio from logarithmic decrement. It also computes the Energy Dissipation Capacity Index (earthquake dissipated energy as a percentage of the cyclic baseline) and the equivalent viscous damping ratio (dissipated energy divided by 4π times the strain energy). These peaks are collected into response spectra against period, with mean, median, quartile and standard-deviation bands, then exported to Excel and used for fragility curves, histograms, scatter and cluster plots, random forest, correlation heatmap and multiple regression. Everything is done separately for X and Y, and the script finishes with hysteresis loops, time histories and an X-Y displacement orbit for the last analysis.

MODEL (OpenSeesPy):
    ops.model('basic', '-ndm', 2, '-ndf', 2)
    Node 1 : fixed base (both DOFs fixed)
    Node 2 : mass node  (DOF 1 = X, DOF 2 = Y)
    One zeroLength element carrying
        - X spring (Hysteretic or Elastic)   + X viscous damper  -> direction 1
        - Y spring (Hysteretic or Elastic)   + Y viscous damper  -> direction 2
    The X and Y springs are independent (uncoupled) and can have different FY, Ke, OSF and DSU (see SPRING).

LOADING:
    Two UniformExcitation patterns: ground motion X on DOF 1 and ground motion Y on DOF 2
    (files Ground_Acceleration_{III_X+1}.txt and Ground_Acceleration_{III_Y+1}.txt, values in g).

PER-DIRECTION OUTPUTS (X and Y):
    displacement, velocity, acceleration, base reaction, ductility damage index (DI), damping ratio,
    stiffness, energy dissipation capacity index (EDCI) and equivalent viscous damping ratio (EVDR).

EDCI   = 100 * Ed(earthquake) / Ed(cyclic displacement protocol)   (computed separately for X and Y)
EVDR   = 100 * Ed / (4*pi*Es),  Es = 0.5*Fmax*umax
DI     = (Dd - Dy) / (Du - Dy)
