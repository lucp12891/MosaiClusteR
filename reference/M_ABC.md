# Multi-source Aggregating Bundles of Clusters (M-ABC), C++-accelerated

Extends the single-source ABC algorithm of (Amaratunga et al. 2008) to
\\K \ge 2\\ sources: runs
[`ABCpp.SingleInMultiple`](https://lucp12891.github.io/MosaiClusteR/reference/ABCpp.SingleInMultiple.md)
per source, row-stacks the per-iteration label matrices, and builds a
single consensus dissimilarity via
[`f.clustABC.MultiSource`](https://lucp12891.github.io/MosaiClusteR/reference/f.clustABC.MultiSource.md)
(C++ kernel when available).

## Usage

``` r
M_ABC(
  List,
  distmeasure = rep("euclidean", length(List)),
  weighting = rep("var", length(List)),
  normalize = rep(FALSE, length(List)),
  gr = c(),
  bag = rep(TRUE, length(List)),
  numsim = 1000,
  numvar = rep(100, length(List)),
  linkage = rep("ward.D2", length(List)),
  NC = NULL,
  dissimilarity = "mabc",
  final_linkage = "ward",
  alpha = 0.625,
  mds = FALSE,
  nfeat = NULL,
  nugget_type = "between",
  nugget_args = list(),
  seed = NULL
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

- dissimilarity:

  Consensus dissimilarity, `"mabc"` or `"abc"`.

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

- seed:

  Optional integer RNG seed. `M_ABC` is stochastic (it bootstraps
  objects and subsamples features), so set `seed` (or call
  [`set.seed`](https://rdrr.io/r/base/Random.html) beforehand) for
  reproducible results.

## Value

A list of class `"M_ABC"` with `DistM` and `Clust`.

## Data-nugget weighting

Set `weighting = "nugget"` (per source, recycled) to weight features by
the data-nugget between-group variance instead of the classic variance;
see
[`ABCpp.SingleInMultiple`](https://lucp12891.github.io/MosaiClusteR/reference/ABCpp.SingleInMultiple.md)
and
[`create_data_nuggets`](https://lucp12891.github.io/MosaiClusteR/reference/create_data_nuggets.md).

## References

Amaratunga D, Cabrera J, Kovtun V (2008). “Microarray Learning with
ABC.” *Biostatistics*, **9**(1), 128–136.
[doi:10.1093/biostatistics/kxm017](https://doi.org/10.1093/biostatistics/kxm017)
.
<https://academic.oup.com/biostatistics/article-abstract/9/1/128/254198>.

## See also

[`ABCpp.SingleInMultiple`](https://lucp12891.github.io/MosaiClusteR/reference/ABCpp.SingleInMultiple.md),
[`f.clustABC.MultiSource`](https://lucp12891.github.io/MosaiClusteR/reference/f.clustABC.MultiSource.md),
[`create_data_nuggets`](https://lucp12891.github.io/MosaiClusteR/reference/create_data_nuggets.md)

## Examples

``` r
data(mosaic_toy)
fit <- M_ABC(mosaic_toy$List, numsim = 50, NC = 3,
               weighting = c("var", "nugget"))
#> M_ABC: source 1 of 2 (weighting = var) ...
#> M_ABC: source 2 of 2 (weighting = nugget) ...
#> M_ABC: building consensus dissimilarity (C++ kernel) ...
table(stats::cutree(stats::as.hclust(fit$Clust), 3))
#> 
#>  1  2  3 
#> 20 20 20 
```
