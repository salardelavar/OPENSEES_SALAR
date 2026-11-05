# MONOTONIC PUSHOVER ANALYSIS OF A ONE-BAY, TWO-STORY, SHEAR-CRITICAL CONCRETE FRAME IN OPENSEES

This OpenSeesPy program implements a nonlinear monotonic displacement-controlled pushover analysis of a one-bay, two-story shear-critical reinforced concrete plane frame, based on the experimental and analytical research presented in Serhan Güner’s Ph.D. thesis, *Performance Assessment of Shear-Critical Reinforced Concrete Plane Frames* (University of Toronto, 2008), and related publications by Serhan Güner and Frank J. Vecchio. The numerical model consists of a two-dimensional frame with fixed base supports, nonlinear beam-column elements with five integration points, P-Delta geometric transformation to account for second-order effects, and fiber-based reinforced concrete sections incorporating material nonlinearity and confinement effects. The model includes 300 × 400 mm column and beam sections with 20 mm reinforcing bars, confined concrete with a compressive strength of 30 MPa, and reinforcing steel models with or without strain hardening and ultimate strain. Following the application of gravity loads, lateral loads are applied simultaneously at the two upper-story nodes using displacement control up to approximately 83 mm. Throughout the analysis, base shear reactions, axial forces, bending moments, lateral displacements, rotations, and corresponding stiffness parameters are recorded. The resulting base-shear–displacement response is subsequently idealized using area-preserving bilinear, trilinear, quadrilinear, and pentalinear curves to evaluate the structural elastic and post-yield stiffness, ductility ratio, and overstrength factor. Furthermore, the numerical pushover response is compared with experimental data, and the calculated results are exported to Excel for further processing and evaluation. This provides a numerical framework for investigating the strength, ductility, stiffness degradation, and overall lateral response of shear-critical reinforced concrete frames. However, an important modeling distinction should be considered: although the fiber-based OpenSees formulation captures nonlinear material behavior and second-order effects, its explicit ability to reproduce shear-related mechanisms and shear failure depends on the implementation of the section formulation in `CONCRETE_SECTION_FUN_TWO.py`. This differs from Güner’s VecTor5/Modified Compression Field Theory (MCFT)-based approach, which was specifically developed to represent both flexural and shear-related response mechanisms.
The main reason VecTor5 is generally more suitable for shear-critical reinforced concrete frame analysis is that its formulation explicitly couples shear, flexural, and axial responses at the sectional level. A conventional fiber-based beam-column model, such as the one in your current OpenSeesPy implementation, primarily represents axial and flexural nonlinearities through uniaxial material fibers. Unless additional shear constitutive models are introduced, it cannot fully reproduce the interaction between diagonal cracking, aggregate interlock, compression softening, reinforcement response, and shear-induced deformation.

VecTor5 was specifically developed to address this limitation using the Modified Compression Field Theory (MCFT) and subsequent refinements associated with the Disturbed Stress Field Model (DSFM). These approaches account for the behavior of cracked concrete under combined stress conditions, making it possible to simulate shear degradation, shear deformation, and interactions between different failure mechanisms. Güner and Vecchio's 2010 verification study applied their procedure to 33 previously tested specimens, approximately two-thirds of which were shear-critical, and reported good agreement with experimental strengths, stiffnesses, ductilities, and failure modes.


# All information has been extracted from the PhD Thesis:
"PERFORMANCE ASSESSMENT OF SHEAR-CRITICAL REINFORCED CONCRETE PLANE FRAMES - Serhan Güner - University of Toronto - 2008"

 PhD Thesis, Department of Civil Engineering, University of Toronto, 2008. Current analysis methods for new
 reinforced concrete structures rely primarily on linear-elastic principles. In many practical situations,
 however, a more accurate prediction of structural behaviour is required. Nonlinear analysis procedures can
 provide this accuracy, but the specialized software needed for such analyses is limited, and most existing
 programs fail to capture shear-related mechanisms adequately. As a result, they often significantly overestimate
 both strength and ductility in shear-critical structures.This study develops and verifies an analytical procedure
 for the nonlinear analysis of frame structures that accounts for shear effects together with flexural and axial behaviour.
 The work builds on the existing program VecTor5, originally developed at the University of Toronto in the early 1980s.
 VecTor5 is based on the Modified Compression Field Theory (MCFT) and was previously capable of nonlinear frame analysis
 under temperature and monotonic loading. Although it produced generally satisfactory results, several deficiencies remained
 in its computational algorithms.
 The research comprises three main parts:
 Improvement of the original analysis procedure for monotonic loading,
 Extension of the procedure to general loading conditions, including cyclic and reversed-cyclic loading, and
 Further development for dynamic loading, covering time-varying base accelerations, impulse, impact and blast
 forces, initial mass velocities, and constant mass accelerations.
 Each part is supported by extensive verification studies using a large and diverse set of previously tested structures
 reported in the literature. Practical considerations for nonlinear modelling are also discussed, providing guidance for
 general applications.Analyses of 63 previously tested structures (approximately half of which were shear-critical)
 demonstrate that the developed procedure successfully reproduces experimental behaviour. Key response parameters—including
 load–deflection curves, reinforcement strains, crack widths, failure mode, failure displacement, total energy dissipation,
 displacement ductility, and post-peak vibrational characteristics—are captured with a high degree of accuracy.
 
'https://collectionscanada.gc.ca/obj/thesescanada/vol2/002/NR58000.PDF'
