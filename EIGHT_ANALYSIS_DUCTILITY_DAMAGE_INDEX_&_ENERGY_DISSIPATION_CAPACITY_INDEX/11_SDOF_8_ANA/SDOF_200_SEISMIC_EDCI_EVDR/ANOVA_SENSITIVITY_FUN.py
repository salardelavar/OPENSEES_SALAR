def ANOVA_SENSITIVITY_FUN(df, output_col, param_cols=None,
                      n_bins=4, include_interactions=False, plot=True):
    """
    ANOVA with automatic binning of continuous predictors.
    Returns the ANOVA table (with % contribution) and optionally a bar plot.
    
    THIS PYTHON SCRIPT WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI) 
    """
    import matplotlib.pyplot as plt 
    import pandas as pd
    from statsmodels.formula.api import ols
    from statsmodels.stats.anova import anova_lm
    
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

#%%------------------------------------------
"""
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

param_list = ["KI", "DUCT", "OSF", "DPR", "UI"]
# ANOVA – main effects only
anova_table, _ = ANOVA_SENSITIVITY(
    df_sens, output_col="ZETA",
    param_cols=param_list,
    n_bins=4,
    include_interactions=False,
    plot=True,
)
plt.show()
print("\nANOVA table (main effects):")
print(anova_table.round(4))

# ANOVA – with two-way interactions
anova_table_inter, _ = ANOVA_SENSITIVITY(
    df_sens, output_col="ZETA",
    param_cols=param_list,
    n_bins=4,
    include_interactions=True,
    plot=True,
)
plt.show()
print("\nANOVA table (with interactions):")
print(anova_table_inter.round(4))
"""