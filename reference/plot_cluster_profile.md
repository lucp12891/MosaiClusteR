# Plot one external variable by cluster

Boxplot (continuous) or grouped bar plot of column proportions
(categorical) of an external variable across clusters – the per-cluster
distribution and gradient the manuscript should show.

## Usage

``` r
plot_cluster_profile(labels, x, name = deparse(substitute(x)), col = NULL)
```

## Arguments

- labels:

  Cluster labels.

- x:

  The variable (numeric -\> boxplot; factor/character -\> barplot).

- name:

  Axis / title label for the variable.

- col:

  Optional vector of cluster colours.

## Value

Invisibly `NULL`; called for the plot.

## See also

[`characterize_clusters`](https://lucp12891.github.io/MosaiClusteR/reference/characterize_clusters.md),
[`cluster_biology_panel`](https://lucp12891.github.io/MosaiClusteR/reference/cluster_biology_panel.md)
