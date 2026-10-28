# Panel of cluster characterisation plots

Lay out
[`plot_cluster_profile`](https://lucp12891.github.io/MosaiClusteR/reference/plot_cluster_profile.md)
for every external variable in a single figure – boxplots for continuous
variables, bar plots for categorical ones – ready for the manuscript.

## Usage

``` r
cluster_biology_panel(labels, metadata, vars = NULL, ncol = 3, base_size = 13)
```

## Arguments

- labels:

  Cluster labels.

- metadata:

  Data frame of external variables (rows = objects).

- vars:

  Optional subset of columns to show (default: all).

- ncol:

  Number of columns in the panel grid.

- base_size:

  Base font size passed to the plot theme.

## Value

Invisibly `NULL`.

## See also

[`characterize_clusters`](https://lucp12891.github.io/MosaiClusteR/reference/characterize_clusters.md)
