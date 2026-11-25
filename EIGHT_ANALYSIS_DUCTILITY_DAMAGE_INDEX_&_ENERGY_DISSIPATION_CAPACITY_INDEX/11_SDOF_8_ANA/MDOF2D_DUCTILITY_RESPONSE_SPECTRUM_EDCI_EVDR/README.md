# PARAMETRIC NONLINEAR SEISMIC DEMAND SPECTRA OF A 2-DOF (X & Y SPRINGS) SYSTEM UNDER 2 GROUND MOTIONS WITH VARIABLE DUCTILITY CAPACITY 

 طیف‌های تقاضای لرزه‌ای غیرخطی پارامتریک یک سیستم دو درجه آزادی (با فنرهای X و Y) تحت دو شتاب‌نگاشت زلزله با ظرفیت نسبت میرایی متغیر

![alt text](COVER.jpg)

This script builds a 2-DOF OpenSeesPy structure (-ndm 2 -ndf 2) with a fixed base node and a mass node, joined by one zeroLength element that holds an independent hysteretic spring and viscous damper in X and in Y. The ultimate displacement of each spring is tied to ductility through DSU = DUCT * DY (with DY = FY/Ke), and the damping ratio is fixed at 3%. The outer loop sweeps ductility capacity from DUCT_MIN to DUCT_MAX. For each level it first runs a cyclic displacement protocol on the X and Y springs, whose amplitude equals DSU, to get the baseline dissipated energy (loop area by the shoelace formula). The inner loop then sweeps the mass, which changes the elastic period, and for each mass it applies two ground motions together, one along X and one along Y, in a nonlinear Newmark time-history analysis. For each direction it records the peak displacement, velocity, acceleration and base reaction, the ductility damage index, the stiffness and the damping ratio from logarithmic decrement. It also computes the Energy Dissipation Capacity Index (earthquake dissipated energy as a percentage of the cyclic baseline) and the equivalent viscous damping ratio (dissipated energy divided by 4π times the strain energy). These peaks are collected into response spectra against period, one curve per ductility level, with mean, median, quartile and standard-deviation bands. The script then exports everything to Excel and produces histograms, fragility curves, scatter and cluster plots, random forest, a correlation heatmap and multiple regression, all separately for X and Y. It ends with hysteresis loops, time histories and an X-Y displacement orbit for the last analysis.

MODEL (OpenSeesPy):
    ops.model('basic', '-ndm', 2, '-ndf', 2)
    Node 1 : fixed base (both DOFs fixed)
    Node 2 : mass node  (DOF 1 = X, DOF 2 = Y)
    One zeroLength element carrying
        - X spring (Hysteretic or Elastic)   + X viscous damper  -> direction 1
        - Y spring (Hysteretic or Elastic)   + Y viscous damper  -> direction 2
    The X and Y springs are independent (uncoupled) and can have different FY, Ke and OSF (see SPRING).

DUCTILITY CAPACITY:
    DSU = DUCT * DY,  DY = FY / Ke   ->  the ultimate displacement of both springs follows the ductility DUCT.
    DUCT is swept from DUCT_MIN to DUCT_MAX (NUM_DUCT + 1 levels); the damping ratio is fixed (DAMPING_RATIO).
    Because DSU changes, the cyclic baseline (Ed_CP) is recomputed for every ductility level.

LOADING:
    Two UniformExcitation patterns: ground motion X on DOF 1 and ground motion Y on DOF 2
    (files Ground_Acceleration_{III_X+1}.txt and Ground_Acceleration_{III_Y+1}.txt, values in g).

PER-DIRECTION OUTPUTS (X and Y):
    displacement, velocity, acceleration, base reaction, ductility damage index (DI), damping ratio,
    stiffness, energy dissipation capacity index (EDCI) and equivalent viscous damping ratio (EVDR).

EDCI   = 100 * Ed(earthquake) / Ed(cyclic displacement protocol)   (computed separately for X and Y)

EVDR   = 100 * Ed / (4*pi*Es),  Es = 0.5*Fmax*umax

DI     = (Dd - Dy) / (Du - Dy)
