# Characterise clusters with external (biological / clinical) variables

For a clustering, summarise how each external variable differs across
the clusters. Continuous variables get \\n\\, mean, SD, median, IQR and
range *per cluster*, an overall Kruskal-Wallis test, a rank \\\eta^2\\
effect size, and **gradient detection** (are the cluster means
monotonically ordered?). Categorical variables get per-cluster counts
and column percentages, a chi-squared / Fisher test and Cramer's V. The
emphasis is on effect sizes and descriptive summaries rather than
p-values alone.

## Usage

``` r
characterize_clusters(labels, metadata, continuous = NULL, categorical = NULL)
```

## Arguments

- labels:

  Cluster label of each object (vector; `NA` allowed).

- metadata:

  A data frame of external variables, one row per object, in the same
  order as `labels`.

- continuous, categorical:

  Optional character vectors naming the columns to treat as continuous /
  categorical. If `NULL`, numeric columns are taken as continuous and
  the rest as categorical.

## Value

A list with `continuous` (per-cluster summary, long format),
`continuous_tests` (one row per variable: `KW_p`, `eta2`, `gradient` –
the clusters ordered by mean – and a monotone-trend
`trend_rho`/`trend_p`), `categorical` (counts and column percentages)
and `categorical_tests` (`p`, `cramers_v`).

## See also

[`plot_cluster_profile`](https://lucp12891.github.io/MosaiClusteR/reference/plot_cluster_profile.md),
[`cluster_biology_panel`](https://lucp12891.github.io/MosaiClusteR/reference/cluster_biology_panel.md)

## Examples

``` r
data(mosaic_toy)
fit <- M_ABC(mosaic_toy$List, numsim = 40, NC = 3)
#> M_ABC: source 1 of 2 (weighting = var) ...
#> M_ABC: source 2 of 2 (weighting = var) ...
#> M_ABC: building consensus dissimilarity (C++ kernel) ...
cl  <- mosaic_labels(fit, 3)
meta <- data.frame(score = rnorm(length(cl)), grp = sample(c("A", "B"), length(cl), TRUE))
ch  <- characterize_clusters(cl, meta)
ch$continuous_tests
#>       variable      KW_p eta2                           gradient  trend_rho
#> score    score 0.7583901    0 3(-0.0548) < 1(-0.0238) < 2(0.178) 0.09311533
#>         trend_p
#> score 0.4791718
```
