#%% SENSITIVITY PLOTS + ANOVA
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pandas.plotting import parallel_coordinates
from scipy.stats import pearsonr, spearmanr
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm

#%% 1. Tornado Diagram
def plot_tornado(df, output_col, param_cols=None,
                 quantile_low=0.1, quantile_high=0.9,
                 figsize=(10, 6), title=None, xlabel=None, color="steelblue"):
    """Horizontal bar chart showing the effect of each parameter
    when moving from its low to high quantile on the output.
    """
    if param_cols is None:
        param_cols = [c for c in df.columns if c != output_col]

    results = []
    for p in param_cols:
        low_thr  = df[p].quantile(quantile_low)
        high_thr = df[p].quantile(quantile_high)

        low_mask  = df[p] <= low_thr
        high_mask = df[p] >= high_thr

        # Guard against empty groups
        low_mean  = df.loc[low_mask,  output_col].mean()  if low_mask.any()  else np.nan
        high_mean = df.loc[high_mask, output_col].mean()  if high_mask.any() else np.nan
        effect    = high_mean - low_mean

        results.append({"Parameter": p, "Low": low_mean, "High": high_mean, "Range": effect})

    tornado_df = (pd.DataFrame(results)
                    .sort_values("Range", ascending=True)
                    .reset_index(drop=True))

    fig, ax = plt.subplots(figsize=figsize)
    ax.barh(tornado_df["Parameter"], tornado_df["Range"], color=color)
    ax.set_xlabel(xlabel or f"Change in {output_col}")
    ax.set_title(title or f"Tornado Diagram – Sensitivity of {output_col}")
    ax.axvline(0, color="black", linewidth=0.8)
    ax.grid(axis="x", linestyle="--", alpha=0.7)
    fig.tight_layout()
    return fig, ax


#%% 2. Parallel Coordinates
def plot_parallel_coordinates(df, target_col, cols_to_plot=None,
                              class_col=None, n_classes=4,
                              figsize=(12, 6), colormap="viridis"):
    """Parallel-coordinates plot coloured by the target (or by an explicit class column)."""
    if cols_to_plot is None:
        cols_to_plot = [c for c in df.columns if c != target_col]

    plot_df = df[cols_to_plot + [target_col]].copy()

    if class_col is None:
        plot_df["class"] = pd.qcut(plot_df[target_col], q=n_classes, labels=False, duplicates="drop")
        class_col = "class"
    else:
        plot_df["class"] = plot_df[class_col]

    fig, ax = plt.subplots(figsize=figsize)
    parallel_coordinates(plot_df, class_col, cols=cols_to_plot, colormap=colormap, ax=ax)
    ax.set_title(f"Parallel Coordinates – coloured by {target_col}")
    ax.grid(True, linestyle="--", alpha=0.5)
    fig.tight_layout()
    return fig, ax


#%% 3. Pairplot
def plot_pairplot(df, vars_to_plot=None, hue_col=None,
                  diag_kind="hist", plot_kws=None, figsize=(12, 12)):
    """Seaborn pairplot wrapper."""
    if vars_to_plot is None:
        vars_to_plot = df.columns.tolist()

    g = sns.pairplot(
        df,
        vars=vars_to_plot,
        hue=hue_col,
        diag_kind=diag_kind,
        plot_kws=plot_kws or {},
        height=2.5,
        aspect=1,
    )
    g.fig.set_size_inches(*figsize)
    g.fig.tight_layout()
    return g


#%% 4. Faceted Grid
def plot_facet_grid(df, x_var, y_var, row_var, col_var,
                    kind="scatter", hue_var=None,
                    row_bins=3, col_bins=3, figsize=(12, 10)):
    """FacetGrid with quantile-binned row/column variables."""
    df_plot = df.copy()
    df_plot["row_bin"] = pd.qcut(df_plot[row_var], q=row_bins, labels=False, duplicates="drop")
    df_plot["col_bin"] = pd.qcut(df_plot[col_var], q=col_bins, labels=False, duplicates="drop")

    g = sns.FacetGrid(
        df_plot,
        row="row_bin",
        col="col_bin",
        height=3.5,
        aspect=1.2,
        margin_titles=True,
    )

    if kind == "scatter":
        g.map_dataframe(sns.scatterplot, x=x_var, y=y_var, hue=hue_var, alpha=0.6)
    elif kind == "line":
        g.map_dataframe(sns.lineplot, x=x_var, y=y_var, hue=hue_var, marker="o")

    g.set_axis_labels(x_var, y_var)
    g.add_legend()
    g.fig.set_size_inches(*figsize)
    g.fig.tight_layout()
    return g


#%% 5. Correlation / Sensitivity Bars
def plot_sensitivity_bars(df, output_col, param_cols=None,
                          method="pearson", figsize=(8, 5), title=None):
    """Bar chart of Pearson or Spearman correlations with the output."""
    if param_cols is None:
        param_cols = [c for c in df.columns if c != output_col]

    corr_values = []
    for p in param_cols:
        if method == "pearson":
            r, _ = pearsonr(df[p], df[output_col])
        else:
            r, _ = spearmanr(df[p], df[output_col])
        corr_values.append(r)

    corr_df = (pd.DataFrame({"Parameter": param_cols, "Correlation": corr_values})
                 .sort_values("Correlation", ascending=False)
                 .reset_index(drop=True))

    fig, ax = plt.subplots(figsize=figsize)
    ax.bar(corr_df["Parameter"], corr_df["Correlation"], color="darkorange")
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylabel(f"{method.capitalize()} Correlation with {output_col}")
    ax.set_title(title or f"Sensitivity ranking ({method})")
    ax.grid(axis="y", linestyle="--", alpha=0.5)
    fig.tight_layout()
    return fig, ax


# -----------------------------------------------------------------------------
# 6. ANOVA Sensitivity (robust)
# -----------------------------------------------------------------------------
def anova_sensitivity(df, output_col, param_cols=None,
                      n_bins=4, include_interactions=False, plot=True):
    """
    ANOVA with automatic binning of continuous predictors.
    Returns the ANOVA table (with % contribution) and optionally a bar plot.
    """
    if param_cols is None:
        param_cols = [c for c in df.columns if c != output_col]

    df_anova = df.copy()

    # Bin continuous variables into a small number of ordered categories
    for p in param_cols:
        df_anova[p] = pd.qcut(df_anova[p], q=n_bins, labels=False, duplicates="drop")
        df_anova[p] = df_anova[p].astype("category")

    # Build formula
    terms = list(param_cols)
    if include_interactions:
        for i, p1 in enumerate(param_cols):
            for p2 in param_cols[i + 1 :]:
                terms.append(f"{p1}:{p2}")

    formula = f"{output_col} ~ " + " + ".join(terms)

    model = ols(formula, data=df_anova).fit()
    anova_table = anova_lm(model, typ=2)

    # Percentage contribution
    ss = anova_table["sum_sq"]
    total_ss = ss.sum()
    anova_table["percent"] = (ss / total_ss) * 100

    if plot:
        plot_data = (anova_table[anova_table.index != "Residual"]
                       .sort_values("percent", ascending=False))

        fig, ax = plt.subplots(figsize=(10, 6))
        bars = ax.bar(range(len(plot_data)), plot_data["percent"],
                      color="skyblue", edgecolor="black")
        ax.set_xticks(range(len(plot_data)))
        ax.set_xticklabels(plot_data.index, rotation=45, ha="right")
        ax.set_ylabel("Percentage of Sum of Squares (%)")
        ax.set_title(f"ANOVA Sensitivity – {output_col}")
        ax.grid(axis="y", linestyle="--", alpha=0.7)
        fig.tight_layout()
        return anova_table, (fig, ax)

    return anova_table, (None, None)


#%% DATA GENERATION
#np.random.seed(42)
N = 1000                                   # realistic sample size

UI_MAX   = np.random.uniform(0.003, 0.015, N)
KI_MAX   = np.random.uniform(0.8, 1.2, N) * 4_500_000.0
DUCT_MAX = np.random.uniform(6.0, 8.0, N)
OSF_MAX  = np.random.uniform(1.05, 1.25, N)
DPR_MAX  = np.random.uniform(0.01, 0.05, N)

# Synthetic outputs (same physics-inspired formulas)
noise = np.random.normal(0, 0.15, N)

zeta_base = (0.02 * DUCT_MAX) + (0.5 * DPR_MAX) - (1e-7 * KI_MAX) + 0.03
zeta_MAX  = np.clip(zeta_base * (1 + 0.2 * noise), 0.01, 0.40)

disp_base   = (0.8 * UI_MAX) * (DUCT_MAX / 5.0) / (KI_MAX / 4_500_000.0)
disp_FV_MAX = np.clip(disp_base * (1 + 0.2 * np.random.normal(0, 1, N)), 0.001, 0.08)

vel_base    = disp_FV_MAX * np.sqrt(KI_MAX / 5000.0)
velo_FV_MAX = np.clip(vel_base * (1 + 0.25 * np.random.normal(0, 1, N)), 0.001, 1.5)

acc_base    = disp_FV_MAX * (KI_MAX / 5000.0)
acc_FV_MAX  = np.clip(acc_base * (1 + 0.25 * np.random.normal(0, 1, N)), 0.1, 15.0)

period_base        = 2 * np.pi / np.sqrt(KI_MAX / 5000.0) * (1 + 0.1 * (DUCT_MAX - 5.0))
PERIOD_MAX_FV_MAX  = np.clip(period_base * (1 + 0.1 * np.random.normal(0, 1, N)), 0.2, 3.0)

df_sens = pd.DataFrame({
    "UI":     UI_MAX,
    "KI":     KI_MAX,
    "DUCT":   DUCT_MAX,
    "OSF":    OSF_MAX,
    "DPR":    DPR_MAX,
    "ZETA":   zeta_MAX,
    "DISP":   disp_FV_MAX,
    "VEL":    velo_FV_MAX,
    "ACC":    acc_FV_MAX,
    "PERIOD": PERIOD_MAX_FV_MAX,
})

print(f" DataFrame shape: {df_sens.shape}")

#%%------------------------------------------------
# CALL ALL PLOTTING FUNCTIONS
param_list = ["KI", "DUCT", "OSF", "DPR", "UI"]

# 1. Tornado for ZETA
plot_tornado(df_sens, output_col="ZETA", param_cols=param_list)
plt.show()

# 2. Parallel coordinates (optional – can be slow / dense with N=1000)
# plot_parallel_coordinates(df_sens, target_col="ZETA", cols_to_plot=param_list)
# plt.show()

# 3. Pairplot (optional)
# df_sens["ZETA_class"] = pd.qcut(df_sens["ZETA"], q=4, labels=["Q1", "Q2", "Q3", "Q4"])
# plot_pairplot(df_sens, vars_to_plot=["KI", "DUCT", "OSF", "DPR", "ZETA"],
#               hue_col="ZETA_class")
# plt.show()

# 4. Facet grid
plot_facet_grid(df_sens,
                x_var="DUCT", y_var="ZETA",
                row_var="OSF", col_var="DPR",
                hue_var="KI", kind="scatter")
plt.show()

# 5. Sensitivity bar plot (Pearson)
plot_sensitivity_bars(df_sens, output_col="ZETA",
                      param_cols=param_list, method="pearson")
plt.show()

# 6. ANOVA – main effects only
anova_table, _ = anova_sensitivity(
    df_sens, output_col="ZETA",
    param_cols=param_list,
    n_bins=4,
    include_interactions=False,
    plot=True,
)
plt.show()
print("\nANOVA table (main effects):")
print(anova_table.round(4))

# 7. ANOVA – with two-way interactions
anova_table_inter, _ = anova_sensitivity(
    df_sens, output_col="ZETA",
    param_cols=param_list,
    n_bins=4,
    include_interactions=True,
    plot=True,
)
plt.show()
print("\nANOVA table (with interactions):")
print(anova_table_inter.round(4))
#%%------------------------------------------------