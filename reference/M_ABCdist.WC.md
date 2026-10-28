# Multi-source M-ABC fused with WeightedClust

Runs
[`ABCdist.SingleInMultiple`](https://lucp12891.github.io/MosaiClusteR/reference/ABCdist.SingleInMultiple.md)
on each source, then feeds the per-source dissimilarities to
[`WeightedClust`](https://lucp12891.github.io/MosaiClusteR/reference/WeightedClust.md)
(`type = "dist"`) to explore convex combinations over a weight grid and
return the focus-weight clustering.

## Usage

``` r
M_ABCdist.WC(
  List,
  distmeasure = rep("euclidean", length(List)),
  weighting = rep("var", length(List)),
  normalize = rep(list(NULL), length(List)),
  gr = c(),
  bag = rep(TRUE, length(List)),
  numsim = 1000L,
  numvar = rep(100L, length(List)),
  linkage = rep("ward.D", length(List)),
  NC = NULL,
  NC2 = NULL,
  StopRange = FALSE,
  wc_weight = seq(1, 0, -0.1),
  wc_weightclust = 0.5,
  wc_linkage = "ward",
  wc_alpha = 0.625,
  mds = FALSE,
  nugget_type = "between",
  nugget_args = list()
)
```

## Arguments

- List:

  A list of \\K\\ data matrices over the same objects, **objects
  (samples) in rows** and features in columns.

- distmeasure, weighting, normalize, bag, numvar, linkage:

  Per-source vectors (recycled to \\K\\) forwarded to
  [`ABCpp.SingleInMultiple`](https://lucp12891.github.io/MosaiClusteR/reference/ABCpp.SingleInMultiple.md).

- gr:

  Unused; kept for compatibility.

- numsim:

  Iterations per source.

- NC:

  Base-cluster count per source (`NULL` = data-driven).

- NC2:

  Optional number of clusters to cut the final tree into.

- StopRange:

  Passed to
  [`WeightedClust`](https://lucp12891.github.io/MosaiClusteR/reference/WeightedClust.md).

- wc_weight:

  Weight grid for
  [`WeightedClust`](https://lucp12891.github.io/MosaiClusteR/reference/WeightedClust.md).

- wc_weightclust:

  Focus weight whose result is returned in `Clust`.

- wc_linkage, wc_alpha:

  Final clustering controls for
  [`WeightedClust`](https://lucp12891.github.io/MosaiClusteR/reference/WeightedClust.md).

- mds:

  Logical; MDS of the consensus dissimilarities.

- nugget_type, nugget_args:

  Forwarded to the nugget weighting.

## Value

A list of class `"M_ABCdist_WC"` with `source_D`, the full `WC` object,
the focus `DistM` and `Clust`, and the `weight_grid`.

## References

(Amaratunga et al. 2008)

## See also

[`M_ABCdist`](https://lucp12891.github.io/MosaiClusteR/reference/M_ABCdist.md),
[`WeightedClust`](https://lucp12891.github.io/MosaiClusteR/reference/WeightedClust.md)

## Examples

``` r
data(mosaic_toy)
fit <- M_ABCdist.WC(mosaic_toy$List, numsim = 40,
                    wc_weight = seq(1, 0, -0.25), wc_weightclust = 0.5)
#> M_ABCdist.WC: source 1 of 2 (weighting=var) ...
#> M_ABCdist.WC: source 2 of 2 (weighting=var) ...
#> Warning: longer object length is not a multiple of shorter object length
dim(fit$DistM)
#> [1] 60 60
```
