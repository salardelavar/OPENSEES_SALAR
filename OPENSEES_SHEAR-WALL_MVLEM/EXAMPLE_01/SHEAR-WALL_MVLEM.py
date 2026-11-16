############################################################################################################
#                    >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<                      #
#                       CYCLIC ANALYSIS OF RC WALL RW2 USING MVLEM IN OPENSEESPY                           #
#----------------------------------------------------------------------------------------------------------#
#   ORIGINAL TCL EXAMPLE: MVLEM - Multiple-Vertical-Line-Element-Model for RC Walls                        #
#   https://opensees.berkeley.edu/wiki/index.php/MVLEM_-_Multiple-Vertical-Line-Element-Model_for_RC_Walls #
#   Example 1: Simulation of wall flexural response using MVLEM                                            #
#   Specimen: RW2 (Thomsen and Wallace, 1995)                                                              #
#   Created by: Kristijan Kolozvari (kkolozvari@fullerton.edu)  dated: July/2015                           #
#   Rewritted by Salar Delavar Ghashghaei (salar.d.ghashghaei@gmail.com) dated: Dec. 2026                  #     
#   UNITS: N, mm                                                                                           #                                        
############################################################################################################
"""
 This script simulates a real RC wall test (specimen RW2, Thomsen & Wallace, 1995)
 in OpenSeesPy. The wall is 3657.6 mm tall and 101.6 mm thick, and it is modeled
 with 16 stacked MVLEM elements, each made of 8 vertical fibers with steel, confined
 concrete and unconfined concrete materials. First, it applies a constant axial load
 of about 378 kN at the top. Then it pushes and pulls the top of the wall sideways,
 in small displacement steps, through eight drift levels (up to about 2%), with two full
 cycles at each level. If a step fails to converge, the solver switches to backup methods.
 At every step the code records the top displacement, base shear, base curvature, and
 the stress and strain of a boundary fiber. At the end it plots the hysteresis loops,
 the loading protocol, the peak envelope, the curvature, and the steel and concrete
 stress-strain curves.
"""
# PEER REPORT: Analytical Modeling of Reinforced Concrete Walls for Predicting Flexural and Coupled–Shear-Flexural Responses
'https://peer.berkeley.edu/sites/default/files/web_peer607_k._orakcal_l._massone_j._wallace_.pdf'
# PAPER: 
'https://www.researchgate.net/figure/Specimen-RW2-tested-by-Thomsen-and-Wallace-at-Clackson-University-3_fig10_350836242'    

#%% Imports
# ---------------------------------------
import argparse
import os
import time
import numpy as np
import matplotlib.pyplot as plt
import openseespy.opensees as ops
#%% Unit conversion factors (original model: kip, in, ksi)
# ---------------------------------------
IN = 25.4                      # 1 in   = 25.4 mm
KIP = 4448.2216152605          # 1 kip  = 4448.2216 N
KSI = KIP / IN**2              # 1 ksi  = 6.894757 MPa  (N/mm2)

#%% Input data
# ---------------------------------------
MODEL_NAME = "RW2"
DATA_DIR = f"MVLEM_{MODEL_NAME}_N_mm"

H = 144.0 * IN          # [mm] wall height            (144 in)
T = 4.0 * IN            # [mm] wall thickness         (4 in)
N_AXIAL = 85.0 * KIP    # [N] axial load               (85 kip)

CTRL_NODE = 17          # CTRL_NODE is the control node. Node 17 is the top node of the wall, at y = H = 3657.6 mm. The cyclic displacement is imposed on it, and the lateral reference load is applied there.
CTRL_DOF = 1            # CTRL_DOF is the degree of freedom at that node whose displacement is controlled. With ndf = 3 in this 2D model, DOF 1 is horizontal translation (x), DOF 2 is vertical translation (y) and DOF 3 is in-plane rotation.

# analysis parameters (end of MVLEM_RW2.tcl)
DRIFT_PEAKS = [0.000330792, 0.001104233, 0.002925758, 0.004558709,
               0.006625238, 0.010816268, 0.014985823, 0.019655056]   # drift ratios (dimensionless)

DINCR = 0.001 * IN      # displacement increment [mm] (0.001 in = 0.0254 mm)
CYCLE_TYPE = "Full"     # Full / Push / Half
NCYCLES = 2             # NCYCLES is the number of times the full cycle is repeated at each drift level. At each of the 8 peaks in DRIFT_PEAKS, the wall goes through 0 → +peak → 0 → −peak → 0 twice before moving on to the next, larger peak.
FACT = H                # drift ratio -> displacement [mm]

NODE_Y = [y * IN for y in (0, 1.25, 1.75, 2.0, 5.25, 9.0, 18.0, 30.0, 36.0,
                           37.25, 37.75, 54.0, 72.0, 90.0, 108.0, 126.0)] + [H]   # [mm]
# NODE_Y is the list of heights (y-coordinates) of the 17 nodes along the wall, measured from the base. It is in mm; the original Tcl values were in inches, which is why each is multiplied by IN.

#%% Model generation
# ---------------------------------------
def build_model():
    """
    build_model() creates the 17 nodes, fixes the base, defines the steel,
    confined and unconfined concrete and shear materials, and adds the 16
    MVLEM wall elements.
    """
    # Reset model
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)

    # Nodes
    for i, y in enumerate(NODE_Y, start=1):
        ops.node(i, 0.0, float(y))
        
    # Contrain    
    ops.fix(1, 1, 1, 1)

    # steel: SteelMPF tag fyp fyn E0 bp bn R0 a1 a2      [MPa]
    fyYbp, bybp, fyYbn, bybn = 57.3 * KSI, 0.0185, 63.0 * KSI, 0.02     # boundary
    fyYwp, bywp, fyYwn, bywn = 48.8 * KSI, 0.035, 65.0 * KSI, 0.02      # web
    Es, R0, a1, a2 = 29000.0 * KSI, 20.0, 0.925, 0.0015
    ops.uniaxialMaterial("SteelMPF", 1, fyYbp, fyYbn, Es, bybp, bybn, R0, a1, a2)
    ops.uniaxialMaterial("SteelMPF", 2, fyYwp, fyYwn, Es, bywp, bywn, R0, a1, a2)
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/SteelMPF.html

    rouYb, rouYw = 0.029333, 0.003333        # reinforcing ratios (boundary / web)

    # concrete: ConcreteCM tag fpcc epcc Ec rc xcrn ft et rt xcrp -GapClose gap   [MPa]
    fpc, ec0, ft, et, Ec = 6.2 * KSI, -0.0021, 0.295 * KSI, 0.00008, 4500.0 * KSI
    xcrnu, xcrp, ru, rt = 1.039, 10000.0, 7.0, 1.2
    fpcc, ec0c, Ecc, xcrnc, rc = 6.9036 * KSI, -0.0033, 5091.3 * KSI, 1.0125, 7.3049
    ops.uniaxialMaterial("ConcreteCM", 3, -fpcc, ec0c, Ecc, rc, xcrnc, ft, et, rt, xcrp, "-GapClose", 1)  # confined
    ops.uniaxialMaterial("ConcreteCM", 4, -fpc, ec0, Ec, ru, xcrnu, ft, et, rt, xcrp, "-GapClose", 1)     # unconfined
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/ConcreteCM.html

    # shear: large elastic stiffness (flexure-only response)   GAs [N]
    Ashweb, G = 192.0 * IN**2, 1875000.0 * KSI      # [mm2], [MPa]
    ops.uniaxialMaterial("Elastic", 5, G * Ashweb)

    # MVLEM elements: element MVLEM eleTag Dens iNode jNode m c -thick .. -width .. -rho .. ...
    widths = [w * IN for w in (7.5, 1.5, 7.5, 7.5, 7.5, 7.5, 1.5, 7.5)]    # [mm]
    rho = [rouYb, 0.0, rouYw, rouYw, rouYw, rouYw, 0.0, rouYb]
    mat_c = [3, 4, 4, 4, 4, 4, 4, 3]
    mat_s = [1, 2, 2, 2, 2, 2, 2, 1]
    for e in range(1, 17):
        ops.element("MVLEM", e, 0.0, e, e + 1, 8, 0.4,
                    "-thick", *([T] * 8),
                    "-width", *widths,
                    "-rho", *rho,
                    "-matConcrete", *mat_c,
                    "-matSteel", *mat_s,
                    "-matShear", 5)
        # INFO LINK: https://openseespydoc.readthedocs.io/en/stable/src/MVLEM.html
    ops.initialize()  # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/InitialStateAnalysis.html


#%% Gravity
# ---------------------------------------
def gravity(MAX_ITERATIONS, MAX_TOLERANCE):
    """
    gravity() applies the 378,100 N axial load at the top node in 10
    load-control steps, then locks it in with loadConst.
    """
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(CTRL_NODE, 0.0, -N_AXIAL, 0.0)

    ops.integrator("LoadControl", 0.1) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/loadControl.html
    ops.system("BandGeneral")          # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/BandGen.html
    ops.test("NormDispIncr", MAX_TOLERANCE, 100, 0) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/normDispIncr.html
    ops.numberer("RCM")                                        # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/RCM.html
    ops.constraints("Transformation")                          # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/constraints.html
    ops.algorithm("Newton")
    ops.analysis("Static")
    ok = ops.analyze(10)
    if ok != 0:
        raise RuntimeError("Gravity analysis failed")
    ops.loadConst("-time", 0.0)


#%% GeneratePeaks
# ---------------------------------------
def generate_peaks(dmax, dincr=0.01, cycle_type="Full", fact=1.0):
    # List of absolute displacement targets [mm] for one cycle (starts with two zeros)
    """
    generate_peaks(dmax, dincr, cycle_type, fact) returns the list of target displacements
    for one cycle: 0 → +peak → 0 → −peak → 0, in steps of dincr.
    """
    disp = 0.0
    steps = [disp, disp]
    dmax = dmax * fact
    dx = -dincr if dmax < 0 else dincr
    n = int(abs(dmax) / dincr)

    for _ in range(n):                       # 0 -> peak
        disp += dx
        steps.append(disp)
    if cycle_type != "Push":
        for _ in range(n):                   # peak -> 0
            disp -= dx
            steps.append(disp)
        if cycle_type != "Half":
            for _ in range(n):               # 0 -> -peak
                disp -= dx
                steps.append(disp)
            for _ in range(n):               # -peak -> 0
                disp += dx
                steps.append(disp)
    return steps


#%% Static analysis parameters
# ---------------------------------------
def set_static_parameters(dincr, MAX_ITERATIONS, MAX_TOLERANCE):
    """
    set_static_parameters() sets the solver for the cyclic stage: constraints,
    numberer, system, convergence test, Newton algorithm, displacement-control integrator
    and the static analysis.
    """
    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system("BandGeneral")
    ops.test("NormDispIncr", MAX_TOLERANCE, 100, 0)
    ops.algorithm("Newton")
    ops.integrator("DisplacementControl", CTRL_NODE, CTRL_DOF, dincr)
    ops.analysis("Static")


def restore_defaults(MAX_ITERATIONS, MAX_TOLERANCE):
    """
    restore_defaults() puts the convergence test and algorithm back to the normal
    settings after a fallback strategy has changed them.
    """
    ops.test("NormDispIncr", MAX_TOLERANCE, 100, 0)
    ops.algorithm("Newton")


def analyze_one_step(MAX_ITERATIONS, MAX_TOLERANCE):
    """
    analyze_one_step() runs one analysis step and, if it fails, retries with Newton using
    the initial tangent and then with Broyden.
    """
    ok = ops.analyze(1)
    if ok != 0:
        print("Trying Newton with Initial Tangent ..")
        ops.test("NormDispIncr", MAX_TOLERANCE, 2000, 0)
        ops.algorithm("Newton", "-initial")
        ok = ops.analyze(1)
        restore_defaults(MAX_ITERATIONS, MAX_TOLERANCE)
    if ok != 0:
        print("Trying Broyden ..")
        ops.algorithm("Broyden", 500)
        ok = ops.analyze(1)
        restore_defaults(MAX_ITERATIONS, MAX_TOLERANCE)
    return ok


#%% Displacement-controlled cyclic analysis
# ---------------------------------------
def cyclic_analysis(MAX_ITERATIONS, MAX_TOLERANCE):
    """
    cyclic_analysis() adds the lateral reference load, loops over the drift peaks and cycles
    applying each displacement increment, and records displacement, base shear, curvature
    and fiber stress and strain after every step.
    """
    ops.timeSeries("Linear", 200)
    ops.pattern("Plain", 200, 200)
    ops.load(CTRL_NODE, 1.0 * KIP, 0.0, 0.0)

    set_static_parameters(DINCR, MAX_ITERATIONS, MAX_TOLERANCE)

    rec = {k: [] for k in ("step", "level", "disp", "shear", "curv", "shear_def",
                           "eps_b", "sig_s_b", "sig_c_b")}
    step_no = 0
    ok = 0
    for level, dmax in enumerate(DRIFT_PEAKS):
        targets = generate_peaks(dmax, DINCR, CYCLE_TYPE, FACT)
        for cyc in range(NCYCLES):
            d0 = 0.0
            for d1 in targets:
                ops.integrator("DisplacementControl", CTRL_NODE, CTRL_DOF, d1 - d0)
                ops.analysis("Static")
                ok = analyze_one_step(MAX_ITERATIONS, MAX_TOLERANCE)
                if ok != 0:
                    print(f"PROBLEM Cyclic analysis: CtrlNode {CTRL_NODE:03d}, dof {CTRL_DOF}, "
                          f"Disp={ops.nodeDisp(CTRL_NODE, CTRL_DOF):.4f} mm")
                    return rec, False
                d0 = d1
                step_no += 1

                # ---- recorders (replaced by in-memory lists) ----
                ops.reactions()
                rec["step"].append(step_no)
                rec["level"].append(level)
                rec["disp"].append(ops.nodeDisp(CTRL_NODE, CTRL_DOF))      # [mm]
                rec["shear"].append(-ops.nodeReaction(1, 1))               # base shear [N]
                rec["curv"].append(ops.eleResponse(1, "Curvature")[0])     # [1/mm]
                rec["shear_def"].append(ops.eleResponse(1, "shear_force_deformation"))
                eps = ops.eleResponse(1, "fiber_strain")
                sc = ops.eleResponse(1, "fiber_stress_concrete")           # [MPa]
                ss = ops.eleResponse(1, "fiber_stress_steel")              # [MPa]
                rec["eps_b"].append(eps[0])             # fiber 1 = left boundary element
                rec["sig_c_b"].append(sc[0])
                rec["sig_s_b"].append(ss[0])
                if step_no % 5000 == 0:
                    print(f"Load Step: {step_no}   disp = {rec['disp'][-1]:.2f} mm")
    print(f"DONE Cyclic analysis: CtrlNode {CTRL_NODE:03d}, dof {CTRL_DOF}, "
          f"Disp={ops.nodeDisp(CTRL_NODE, CTRL_DOF):.4f} mm")
    return rec, True


#%% Plots
# ---------------------------------------
def plot_results(rec, outfile="RW2_results_N_mm.png"):
    """
    plot_results(rec) draws the six plots: hysteresis loops, displacement protocol,
    force envelope, base curvature, boundary-steel stress–strain and confined-concrete
    stress–strain.
    """
    d = np.array(rec["disp"])                 # [mm]
    V = np.array(rec["shear"]) / 1000.0       # N -> kN for readable axes
    drift = d / H * 100.0
    lev = np.array(rec["level"])
    stp = np.array(rec["step"])

    fig, ax = plt.subplots(2, 3, figsize=(18, 10))

    # 1) hysteresis
    ax[0, 0].plot(d, V, lw=0.9, color="tab:blue")
    ax[0, 0].axhline(0, color="k", lw=0.5); ax[0, 0].axvline(0, color="k", lw=0.5)
    ax[0, 0].set_xlabel("Top displacement [mm]"); ax[0, 0].set_ylabel("Base shear [kN]")
    ax[0, 0].set_title("Base shear - top displacement")
    sec = ax[0, 0].secondary_xaxis("top", functions=(lambda x: x / H * 100, lambda x: x * H / 100))
    sec.set_xlabel("Drift ratio [%]")
    ax[0, 0].grid(alpha=0.3)

    # 2) protocol
    ax[0, 1].plot(stp, drift, lw=0.8, color="tab:red")
    ax[0, 1].axhline(0, color="k", lw=0.5)
    ax[0, 1].set_xlabel("Analysis step"); ax[0, 1].set_ylabel("Drift ratio [%]")
    ax[0, 1].set_title("Displacement protocol (2 cycles per level)")
    ax[0, 1].grid(alpha=0.3)

    # 3) envelope: max / min shear in each drift level
    xs_p, ys_p, xs_n, ys_n = [], [], [], []
    for L in np.unique(lev):
        m = lev == L
        ip = np.argmax(d[m]); im = np.argmin(d[m])
        xs_p.append(drift[m][ip]); ys_p.append(V[m][ip])
        xs_n.append(drift[m][im]); ys_n.append(V[m][im])
    ax[0, 2].plot(xs_p, ys_p, "-o", color="tab:green", label="push (+)")
    ax[0, 2].plot(np.abs(xs_n), np.abs(ys_n), "-s", color="tab:orange", label="pull (-)")
    ax[0, 2].set_xlabel("|Peak drift| [%]"); ax[0, 2].set_ylabel("|Base shear| [kN]")
    ax[0, 2].set_title("Envelope (peak of each drift level)")
    ax[0, 2].legend(); ax[0, 2].grid(alpha=0.3)

    # 4) base curvature vs drift
    ax[1, 0].plot(drift, np.array(rec["curv"]) * 1e6, lw=0.9, color="tab:purple")
    ax[1, 0].axhline(0, color="k", lw=0.5); ax[1, 0].axvline(0, color="k", lw=0.5)
    ax[1, 0].set_xlabel("Drift ratio [%]"); ax[1, 0].set_ylabel("Curvature, element 1 [1/mm x 1e-6] (OpenSees sign)")
    ax[1, 0].set_title("Base-element curvature")
    ax[1, 0].grid(alpha=0.3)

    # 5) steel stress-strain, boundary fiber (element 1, fiber 1)
    ax[1, 1].plot(np.array(rec["eps_b"]) * 100, rec["sig_s_b"], lw=0.9, color="tab:brown")
    ax[1, 1].axhline(0, color="k", lw=0.5); ax[1, 1].axvline(0, color="k", lw=0.5)
    ax[1, 1].set_xlabel("Strain [%]"); ax[1, 1].set_ylabel("Stress [MPa]")
    ax[1, 1].set_title("Boundary steel (SteelMPF), element 1, fiber 1")
    ax[1, 1].grid(alpha=0.3)

    # 6) confined concrete stress-strain, boundary fiber
    ax[1, 2].plot(np.array(rec["eps_b"]) * 100, rec["sig_c_b"], lw=0.9, color="tab:gray")
    ax[1, 2].axhline(0, color="k", lw=0.5); ax[1, 2].axvline(0, color="k", lw=0.5)
    ax[1, 2].set_xlim(min(rec["eps_b"]) * 100 * 1.3, 0.5)     # zoom on the compression range
    ax[1, 2].set_xlabel("Strain [%]"); ax[1, 2].set_ylabel("Stress [MPa]")
    ax[1, 2].set_title("Confined concrete (ConcreteCM), element 1, fiber 1")
    ax[1, 2].grid(alpha=0.3)

    fig.suptitle(f"MVLEM wall specimen {MODEL_NAME} (Thomsen & Wallace 1995) - cyclic response [N, mm, MPa]",
                 fontsize=14)
    fig.tight_layout()
    fig.savefig(outfile, dpi=200)
    print(f"figure saved -> {outfile}")


def save_out(rec):
    """
    save_out(rec) writes displacement, base shear and curvature to text files,
    like the Tcl recorders did.
    """
    os.makedirs(DATA_DIR, exist_ok=True)
    np.savetxt(f"{DATA_DIR}/MVLEM_Dtop_mm.txt", np.column_stack([rec["step"], rec["disp"]]))
    np.savetxt(f"{DATA_DIR}/MVLEM_BaseShear_N_Disp_mm.txt", np.column_stack([rec["disp"], rec["shear"]]))
    np.savetxt(f"{DATA_DIR}/MVLEM_Curvature_1_per_mm.txt", np.column_stack([rec["step"], rec["curv"]]))


#%% MAIN
# ---------------------------------------
"""
The #%% MAIN block at the bottom is not a function. It parses --save-out and --no-show,
 builds the model, applies gravity, runs the cycles and plots the results.
"""
if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--save-out", action="store_true", help="write text files like the Tcl recorders")
    ap.add_argument("--no-show", action="store_true", help="do not open the plot window")
    args = ap.parse_args()

    MAX_ITERATIONS = 2000   # Maximum number of iterations
    MAX_TOLERANCE = 1.0e-4  # Specified tolerance for convergence
    
    t0 = time.time()
    build_model()
    gravity(MAX_ITERATIONS, MAX_TOLERANCE)
    print("Model generated and gravity load applied successfully")
    rec, ok = cyclic_analysis(MAX_ITERATIONS, MAX_TOLERANCE)
    print(f"elapsed: {time.time() - t0:.1f} s, steps recorded: {len(rec['step'])}")

    if len(rec["step"]):
        d = np.array(rec["disp"]); V = np.array(rec["shear"])
        print(f"max displacement {d.max():.2f} mm ({d.max()/H*100:.2f} % drift), "
              f"max base shear {V.max()/1000:.1f} kN ({V.max():.0f} N), min {V.min()/1000:.1f} kN")
        if args.save_out:
            save_out(rec)
        plot_results(rec)
        if not args.no_show:
            plt.show()