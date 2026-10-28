# Justify the number of clusters (silhouette + elbow / scree)

For a fused dissimilarity (or a fitted result, or a `List` of sources)
compute, over a range of \\k\\, the average silhouette width and the
total within-cluster dispersion (an elbow / scree curve), and plot them
so the chosen \\k\\ can be justified in the manuscript.

## Usage

``` r
cluster_selection_plot(
  x,
  k_range = 2:8,
  linkage = "ward.D2",
  select = c("both", "elbow", "silhouette"),
  plot = TRUE,
  point_size = 2.6,
  base_size = 13
)
```

## Arguments

- x:

  A fitted result (with `$DistM`), a dissimilarity matrix, or a `List`
  of data sources (objects in rows).

- k_range:

  Integer vector of candidate cluster counts.

- linkage:

  Agglomeration method used to cut the dissimilarity.

- select:

  Which criterion to show *and* use for `best_k`: `"both"` (default)
  draws both panels and picks \\k\\ from the elbow; `"elbow"` draws only
  the within-cluster dispersion (scree) curve and uses its knee;
  `"silhouette"` draws only the silhouette curve and uses its maximum.
  Passing a single criterion therefore shows a single panel.

- plot:

  Logical; draw the figure.

- point_size:

  Marker size for the curve points.

- base_size:

  Base font size for the ggplot2 theme.

## Value

Invisibly a data frame with `k`, `silhouette` and `within_dispersion`,
plus attributes `"best_k"`, `"elbow_k"`, `"silhouette_k"` and `"plot"`
(the ggplot object).

## See also

[`embedding_plot`](https://lucp12891.github.io/MosaiClusteR/reference/embedding_plot.md),
[`SelectnrClusters`](https://lucp12891.github.io/MosaiClusteR/reference/SelectnrClusters.md)

## Examples

``` r
data(mosaic_toy)
fit <- M_ABC(mosaic_toy$List, numsim = 40, NC = 3)
#> M_ABC: source 1 of 2 (weighting = var) ...
#> M_ABC: source 2 of 2 (weighting = var) ...
#> M_ABC: building consensus dissimilarity (C++ kernel) ...
cluster_selection_plot(fit, k_range = 2:6, select = "elbow")   # one panel
```
