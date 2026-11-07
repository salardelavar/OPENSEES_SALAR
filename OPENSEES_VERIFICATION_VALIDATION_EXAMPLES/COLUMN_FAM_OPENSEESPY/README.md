# MONOTONIC PUSHOVER ANALYSIS OF STRUCTURE WITH FORCE ANALOGY METHOD (FAM) AND OPENSEESPY VERIFICATION
# تحلیل پوش‌اور یکنواخت سازه با روش آنالوژی نیرو و صحت‌سنجی با اوپنسیس

This python script performs a displacement-controlled pushover analysis of a fixed-guided column 
using the Force Analogy Method (FAM). 
1. The global stiffness of the structure remains constant at its initial elastic state (K_elastic).
2. Material nonlinearity (hinge formation) is captured by solving for plastic rotations (theta_p)
   at the base and top, which act as analog strain-like parameters reducing internal forces.
3. The moment and shear equations are decomposed into elastic and plastic components.
4. A return-mapping algorithm tracks the onset of yielding and computes the plastic rotations.
5. Three subplots are generated: Pushover Curve, Moment Evolution, and Plastic Rotation Evolution.

![alt text](COVER.png)

# COMPOSITE_ELASTIC_SECTION_ANALYSIS_FUN
![alt text](COMPOSITE_ELASTIC_SECTION_ANALYSIS_FUN.png)
