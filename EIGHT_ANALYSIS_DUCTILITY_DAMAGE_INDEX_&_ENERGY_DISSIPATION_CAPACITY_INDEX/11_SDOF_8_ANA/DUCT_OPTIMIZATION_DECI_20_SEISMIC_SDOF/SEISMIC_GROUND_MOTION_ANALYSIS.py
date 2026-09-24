# =============================================================================
#     >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<
#                         SEISMIC GROUND MOTION ANALYSIS
#------------------------------------------------------------------------------
#  Purpose : Compute the peak absolute acceleration (PGA proxy) for 200
#            seismic records and visualize their distribution.
#
#  Input   : Ground_Acceleration_1.txt ... Ground_Acceleration_200.txt
#  Output  : - Console: statistical summary (mean, std, min, Q1-Q4, IQR)
#            - File   : max_acceleration_histogram.png
#------------------------------------------------------------------------------
# THIS PYTHON SCRIPT IS WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
# =============================================================================
"""
The code reads 200 seismic acceleration text files named `Ground_Acceleration_1.txt` through 
`Ground_Acceleration_200.txt`, loads each file with `np.loadtxt` (falling back to comma-delimited 
parsing if needed), selects the acceleration column if the data has multiple columns, computes the
 maximum absolute acceleration for each record, and stores the results in a NumPy array while skipping
 any missing files or NaN values. It then calculates key statistics—sample mean, standard deviation,
 minimum, quartiles Q1, Q2 (median), Q3, Q4 (maximum), and the interquartile range—and prints them to
 the console. Finally, it creates a histogram of all 200 peak values, overlays vertical reference lines
 for the mean, Q1, Q2, Q3, minimum, and Q4/maximum, adds a shaded mean ± 1σ band, inserts a text box
 summarizing the statistics, labels the axes and legend, and saves the figure as `max_acceleration_histogram.png`
 while displaying it on screen.
"""
# -----------------------------------------------------------------------------
# 1. IMPORTS & GLOBAL CONFIGURATION
# -----------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

FOLDER  = "."                       # folder containing the .txt files
N_FILES = 200                       # number of records to process
PREFIX  = "Ground_Acceleration_"    # filename prefix
N_BINS  = 20                        # histogram bin count


# -----------------------------------------------------------------------------
# 2. DATA INGESTION — READ FILES & EXTRACT PEAK |ACCELERATION|
# -----------------------------------------------------------------------------
max_vals = []
file_ids = []

for i in range(1, N_FILES + 1):
    filepath = Path(FOLDER) / f"{PREFIX}{i}.txt"

    if not filepath.exists():
        print(f"⚠ Missing file: {filepath}")
        continue

    try:
        data = np.loadtxt(filepath)
    except ValueError:
        data = np.loadtxt(filepath, delimiter=",")

    acc = data[:, 1] if data.ndim > 1 else data      # acceleration column

    max_vals.append(np.max(np.abs(acc)))
    file_ids.append(i)

max_vals = np.array(max_vals)
file_ids = np.array(file_ids)


# -----------------------------------------------------------------------------
# 3. STATISTICAL SUMMARY — MEAN, STD, MIN, Q1, Q2, Q3, Q4, IQR
# -----------------------------------------------------------------------------
mean   = np.mean(max_vals)
std    = np.std(max_vals, ddof=1)                    # sample standard deviation
vmin   = np.min(max_vals)
vmax   = np.max(max_vals)
Q1, Q2, Q3, Q4 = np.percentile(max_vals, [25, 50, 75, 100])
IQR    = Q3 - Q1

print("=" * 50)
print(f"  N      = {len(max_vals)}")
print(f"  Mean   = {mean:.4f}")
print(f"  Std    = {std:.4f}")
print(f"  Min    = {vmin:.4f}")
print(f"  Q1     = {Q1:.4f}")
print(f"  Q2     = {Q2:.4f}   (= median)")
print(f"  Q3     = {Q3:.4f}")
print(f"  Q4     = {Q4:.4f}   (= max)")
print(f"  IQR    = {IQR:.4f}")
print("=" * 50)


# -----------------------------------------------------------------------------
# 4. BAR CHART — PEAK |a| PER RECORD
# -----------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(14, 6))

ax.bar(file_ids, max_vals, color="steelblue",
       edgecolor="black", width=0.8)
ax.plot(file_ids, max_vals, color="darkred", linewidth=1, marker="o",
        markersize=3, alpha=0.7, label="Max |a|")

# Highlight the global maximum
imax = int(np.argmax(max_vals))
ax.annotate(f"Max = {max_vals[imax]:.4f}\n(file #{file_ids[imax]})",
            xy=(file_ids[imax], max_vals[imax]),
            xytext=(file_ids[imax], max_vals[imax] * 1.1),
            ha="center", fontsize=10, color="darkred",
            arrowprops=dict(arrowstyle="->", color="darkred"))

ax.set_xlabel("Ground Motion Record #", fontsize=12)
ax.set_ylabel("Max |Acceleration|", fontsize=12)
ax.set_title("Peak Absolute Acceleration for 200 Seismic Records",
             fontsize=13, fontweight="bold")
ax.grid(True, linestyle="--", alpha=0.5)
ax.set_xlim(0, N_FILES + 1)              # ← was: n_files (lowercase, undefined)
ax.legend()
plt.tight_layout()
plt.savefig("max_acceleration_chart.png", dpi=200)
plt.show()


# -----------------------------------------------------------------------------
# 5. FIGURE SETUP — HISTOGRAM
# -----------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(12, 7))


# -----------------------------------------------------------------------------
# 6. HISTOGRAM OF PEAK ABSOLUTE ACCELERATION
# -----------------------------------------------------------------------------
counts, bins, patches = ax.hist(
    max_vals,
    bins=N_BINS,
    color="steelblue",
    edgecolor="black",
    alpha=0.85,
    label=f"Max |a|  (N = {len(max_vals)})",
)

ymax_plot = counts.max() * 1.18                      # headroom for annotations


# -----------------------------------------------------------------------------
# 7. VERTICAL REFERENCE LINES — MEAN, QUARTILES, MIN/MAX
# -----------------------------------------------------------------------------
reference_lines = [
    (mean, "red",        "-",  f"Mean = {mean:.3f}",        1.00),
    (Q2,   "green",      "--", f"Q2 (median) = {Q2:.3f}",   0.88),
    (Q1,   "darkorange", ":",  f"Q1 = {Q1:.3f}",            0.76),
    (Q3,   "darkorange", ":",  f"Q3 = {Q3:.3f}",            0.76),
    (vmin, "purple",     "-.", f"Min = {vmin:.3f}",         0.64),
    (Q4,   "black",      "-.", f"Q4 (max) = {Q4:.3f}",      0.64),
]

for x, color, ls, label, frac in reference_lines:
    ax.axvline(x, color=color, linestyle=ls, linewidth=2,
               label=label, ymax=frac)


# -----------------------------------------------------------------------------
# 8. MEAN ± 1σ SHADED BAND
# -----------------------------------------------------------------------------
ax.axvspan(mean - std, mean + std, color="red", alpha=0.08,
           label=f"Mean ± 1σ ({std:.3f})")
ax.axvline(mean - std, color="red", linestyle="--", linewidth=1, alpha=0.5)
ax.axvline(mean + std, color="red", linestyle="--", linewidth=1, alpha=0.5)


# -----------------------------------------------------------------------------
# 9. STATISTICS TEXT BOX
# -----------------------------------------------------------------------------
stats_txt = (
    f"Statistics\n"
    f"──────────────\n"
    f"N     = {len(max_vals)}\n"
    f"Mean  = {mean:.4f}\n"
    f"Std   = {std:.4f}\n"
    f"Min   = {vmin:.4f}\n"
    f"Q1    = {Q1:.4f}\n"
    f"Q2    = {Q2:.4f}\n"
    f"Q3    = {Q3:.4f}\n"
    f"Q4    = {Q4:.4f}\n"
    f"IQR   = {IQR:.4f}"
)

ax.text(0.985, 0.97, stats_txt,
        transform=ax.transAxes,
        ha="right", va="top",
        fontsize=10, family="monospace",
        bbox=dict(boxstyle="round,pad=0.5",
                  facecolor="lightyellow",
                  edgecolor="gray", alpha=0.95))


# -----------------------------------------------------------------------------
# 10. LABELS, TITLES & COSMETICS
# -----------------------------------------------------------------------------
ax.set_xlabel("Max |Acceleration|", fontsize=12)
ax.set_ylabel("Frequency (count)", fontsize=12)
ax.set_title("Histogram of Peak Absolute Acceleration — 200 Seismic Records",
             fontsize=13, fontweight="bold")
ax.set_ylim(0, ymax_plot)
ax.grid(True, linestyle="--", alpha=0.4, axis="y")
ax.legend(loc="upper left", fontsize=9, framealpha=0.9)


# -----------------------------------------------------------------------------
# 11. EXPORT & DISPLAY
# -----------------------------------------------------------------------------
plt.tight_layout()
plt.savefig("max_acceleration_histogram.png", dpi=200)
plt.show()