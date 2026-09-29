# SENSITIVITY ANALYSIS OF STRUCTURE DUCTILITY RATIO  AND OPENSEES VIA 20 SEISMIC GROUND MOTIONS

![alt text](https://github.com/salardelavar/OPENSEES_SALAR/blob/main/EIGHT_ANALYSIS_DUCTILITY_DAMAGE_INDEX_%26_ENERGY_DISSIPATION_CAPACITY_INDEX/11_SDOF_8_ANA/DUCT_SENSITIVITY_20_SEISMIC_SDOF/COVER.png) 

Assume that a single-degree-of-freedom structure is subjected to nonlinear dynamic analysis under
 twenty different ground motion records. To investigate the influence of structure ductility on
 seismic performance, a sensitivity analysis is performed by systematically varying each column’s
 ductility ratio from 20 to 50 in 21 steps. For each ductility level, a pushover analysis provides a cyclic displacement analysis establishes a reference hysteretic energy, and twenty nonlinear time-history analyses are conducted using different ground motions. From these analyses, the Energy Dissipation Capacity Index (EDCI) and other response metrics such as maximum displacement, velocity, acceleration, damage index, over-strength factor, ductility ratio, and equivalent viscous damping ratio are computed.
 The median response over the twenty ground motions is then evaluated for each ductility level.
 Finally, the results are assessed through trend plots, 3D contour surfaces, correlation heatmaps,
 Random Forest regression, and ANOVA sensitivity analysis to identify the optimal column ductility
 ratio that maximizes energy dissipation capacity while maintaining structural safety.

# SENSITIVITY ANALYSIS BY CHANGING STRUCTURE DUCTILITY RATIO
1. Sets MAT_TYPE = 'INELASTIC' and sweeps the COLUMN DUCTILITY RATIO from 20.0 to 50.0
   in 21 linear steps (DUCT_MIN -> DUCT_MAX).  COL_DUCT stores each swept value.
2. Declares ~20 accumulator lists: raw per-record responses, plus *_MED lists that
   will hold one MEDIAN value per ductility level.
3. Starts a CPU timer.  Outer loop `for JJ in range(0, 21)` computes DUCT and resets
   the per-ductility raw accumulators (DISP, VELO, ACC, EDCI, DII, OMEGA, MU, RR, EDVR).
4. STEP 1 -- PUSHOVER: SDOF(DUCT, MAT_TYPE, TOTAL_MASS, 'PUSHOVER', II) returns the
   capacity curve (disp_PUSH, reaction_PUSH).
5. Stores SDOF_ef_* arrays and PERIOD_PUSH = max(PERIOD_MAX_PUSH)  (secant period).
6. STEP 2 -- CYCLIC_DISPLACEMENT: benchmark constant-amplitude hysteresis loop.
   S055.EQULIVALENT_VISCOUS_DAMPING_RATIO_FUN returns zeta_CP (Jacobsen area method).
7. STEP 3 -- SEISMIC LOOP `for II in range(40, 60)`: 20 nonlinear time-history runs.
8. Each run unpacks time/reaction/disp/velo/acc/DI/stiffness/period arrays from SDOF().
9. S12.DUCTILITY_DAMAGE_INDEX_FUN fits the pushover curve and returns Park-Ang style
   DIx, over-strength Omega_0, ductility mu, and behaviour coefficient R.
10. S10.ENERGY_DISSIPATION_CAPACITY_INDEX compares seismic vs cyclic dissipated energy
   -> EDCI (%).  Peak DISP/VELO/ACC also appended.
11. Ratio metrics appended: PERIOD_RATIO, DI_RATIO, REACTION_RATIO, DISP_RATIO
   (each = SEISMIC / PUSHOVER).
12. S055 again computes zeta_SEI for the seismic loop -> EDVR and EDVR_RATIO = zeta_SEI/zeta_CP.
13. STEP 4 -- after 20 records, np.median() reduces each raw list to ONE robust value
   (DISP_MED, DI_RATIO_MED, EDCI_MED, EDVR_MED, ...), the design-point response.
14. Outer loop ends after 21 design points.  Total CPU time is printed.
15. PLOTTING: ~17 calls to S01.PLOT_SCATTER(COL_DUCT, Y, ...) draw scatter + polynomial
   fit of order 1/3/7 for every response vs. the ductility ratio.
16. 5 calls to S11.PLOT_CONTOUR_3D_2D_FUN build 2D contours + 3D surfaces pairing two
   responses (e.g. DISP_RATIO vs DI_RATIO) against the ductility ratio axis.
17. MACHINE LEARNING: builds a pandas DataFrame with features DISP, VEL, ACC,
   REACTION_RATIO, EDVR_RATIO, DI_RATIO, then calls S01.RANDOM_FOREST(df) to train a
   RandomForestClassifier (safe/unsafe) + Regressor (safety likelihood) and print
   accuracy, MSE, R2, and feature importances.
18. S01.PLOT_HEATMAP(df) draws the Pearson correlation matrix of the six features.
19. S099.SENSITIVITY_HEATMAP_FUN(X=[DISP,VEL,ACC], Y=[EDVR_RATIO,DUCT,DI_RATIO],
   method='pearson') returns a 3x3 sensitivity-coefficient heatmap.
20. S100.ANOVA_SENSITIVITY_FUN(df_sens, output_col="DISP", param_cols=[VEL,ACC,
   EDVR_RATIO,DI_RATIO], n_bins=4) bins continuous predictors and runs classical ANOVA.
21. Two ANOVA runs: (a) MAIN EFFECTS only, (b) MAIN + TWO-WAY INTERACTIONS; each prints
   a SumSq / df / F / p-value table and shows the diagnostic plot.

22. Engineering reading: plots reveal whether higher column ductility reduces damage
   (DII), damping (EDVR), energy index (EDCI) and peak demand -- the core question of
   performance-based seismic design.

23. ML + correlation + ANOVA together rank which EDP (disp, vel or acc) drives damage.
