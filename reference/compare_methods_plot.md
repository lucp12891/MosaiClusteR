# Visually compare clustering solutions with ComparePlot (no ground truth)

Aligns several clustering solutions of the same objects and draws them
with
[`ComparePlot`](https://lucp12891.github.io/MosaiClusteR/reference/ComparePlot.md)
so their agreement can be read off directly – which objects the methods
place together and where they disagree – without any ground truth. Each
solution may be a fitted result (with `$Clust` / `$DistM`) or a plain
label vector; label vectors are turned into a co-membership
dissimilarity so every solution is comparable.

## Usage

``` r
compare_methods_plot(
  solutions,
  nrclusters,
  cols = NULL,
  cex.labels = 0.5,
  margins = c(5, 11, 3, 2),
  ...
)
```

## Arguments

- solutions:

  A named list of fitted results and/or label vectors, all over the same
  objects in the same order.

- nrclusters:

  Number of clusters to cut each solution into.

- cols:

  Optional cluster colours; defaults to a distinct palette.

- cex.labels, margins:

  Passed to
  [`ComparePlot`](https://lucp12891.github.io/MosaiClusteR/reference/ComparePlot.md).

- ...:

  Further arguments to
  [`ComparePlot`](https://lucp12891.github.io/MosaiClusteR/reference/ComparePlot.md).

## Value

Invisibly `NULL`; called for the plot. Needs circlize and plotrix.

## See also

[`clustering_concordance`](https://lucp12891.github.io/MosaiClusteR/reference/clustering_concordance.md),
[`ComparePlot`](https://lucp12891.github.io/MosaiClusteR/reference/ComparePlot.md)

## Examples

``` r
data(mosaic_toy)
sols <- list(
  var    = M_ABC(mosaic_toy$List, numsim = 40, NC = 3, weighting = c("var", "var")),
  nugget = M_ABC(mosaic_toy$List, numsim = 40, NC = 3, weighting = c("nugget", "nugget")),
  NEMO   = NEMO(mosaic_toy$List, k = 3)$cluster)
#> M_ABC: source 1 of 2 (weighting = var) ...
#> M_ABC: source 2 of 2 (weighting = var) ...
#> M_ABC: building consensus dissimilarity (C++ kernel) ...
#> M_ABC: source 1 of 2 (weighting = nugget) ...
#> M_ABC: source 2 of 2 (weighting = nugget) ...
#> M_ABC: building consensus dissimilarity (C++ kernel) ...
if (FALSE) compare_methods_plot(sols, nrclusters = 3) # \dontrun{}
```
