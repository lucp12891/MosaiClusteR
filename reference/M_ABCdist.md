# Multi-source M-ABC with distance accumulation

Runs
[`ABCdist.SingleInMultiple`](https://lucp12891.github.io/MosaiClusteR/reference/ABCdist.SingleInMultiple.md)
on each source to obtain a `[0, 1]` dissimilarity, combines them as a
weighted mean \\D = \sum_k w_k D_k\\ and applies a final agglomerative
clustering.

## Usage

``` r
M_ABCdist(
  List,
  source_weights = NULL,
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
  final_linkage = "ward",
  alpha = 0.625,
  mds = FALSE,
  nfeat = NULL,
  nugget_type = "between",
  nugget_args = list()
)
```

## Arguments

- List:

  A list of \\K\\ data matrices over the same objects, **objects
  (samples) in rows** and features in columns.

- source_weights:

  Optional non-negative source weights (normalised to sum to one);
  `NULL` = equal weights.

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

- final_linkage, alpha:

  Final clustering controls.

- mds:

  Logical; MDS of the consensus dissimilarities.

- nfeat:

  Feature-subsample size forwarded to
  [`ABCpp.SingleInMultiple`](https://lucp12891.github.io/MosaiClusteR/reference/ABCpp.SingleInMultiple.md)
  (`NULL` = \\\lfloor\sqrt{G}\rfloor\\).

- nugget_type, nugget_args:

  Forwarded to the nugget weighting.

## Value

A list of class `"M_ABCdist"` with `source_D`, `DistM`, `weights` and
`Clust`.

## References

(Amaratunga et al. 2008)

## See also

[`M_ABC`](https://lucp12891.github.io/MosaiClusteR/reference/M_ABC.md),
[`ABCdist.SingleInMultiple`](https://lucp12891.github.io/MosaiClusteR/reference/ABCdist.SingleInMultiple.md),
[`M_ABCdist.WC`](https://lucp12891.github.io/MosaiClusteR/reference/M_ABCdist.WC.md)

## Examples

``` r
data(mosaic_toy)
fit <- M_ABCdist(mosaic_toy$List, numsim = 40, weighting = c("var", "nugget"))
#> M_ABCdist: source 1 of 2 (weight=0.5, weighting=var) ...
#> M_ABCdist: source 2 of 2 (weight=0.5, weighting=nugget) ...
dim(fit$DistM)
#> [1] 60 60
```
