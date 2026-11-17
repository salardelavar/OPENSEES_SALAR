# INCREMENTAL DYNAMIC ANALYSIS (IDA) OF A MULTI NONLINEAR SHEAR BUILDING





PLOT_FRAGILITY_FUN.py turns your IDA collapse results into a fragility curve, which gives the probability of collapse as a function of ground-motion intensity. FIT_LOGNORMAL_FRAGILITY counts, at each intensity level, how many of the 200 records have collapsed, then finds the two parameters of a lognormal curve, P = Φ(ln(x/θ)/β), that best explain those counts. It does this by maximum likelihood: it minimizes the negative log-likelihood, treating the number of collapses at each level as a binomial outcome. The parameters are optimized as logarithms so that the median θ and the dispersion β always stay positive, and the starting guess comes from each record’s first collapse level. Because it works with counts per level, records that never collapse up to the highest intensity are handled automatically. PLOT_FRAGILITY_FUN calls the fit and draws the observed collapse fractions as dots, the fitted curve as a smooth line, a dashed vertical line at the median θ (the intensity where 50% of records collapse), and a legend showing θ and β. It returns theta and beta. In the __main__ block it loads IDA_DATA.npz from your earlier run, so you can produce the fragility curve without repeating the analysis.

