# Single-source Aggregating Bundles of Clusters (ABC), C++-accelerated engine

The canonical single-source ABC engine of (Amaratunga et al. 2008) .
Over `numsim` iterations it (optionally) bootstraps the objects, draws a
weighted subsample of features, clusters the sub-matrix with
[`fastcluster::hclust`](https://rdrr.io/pkg/fastcluster/man/hclust.html),
and records each selected object's base-cluster label.
[`M_ABC`](https://lucp12891.github.io/MosaiClusteR/reference/M_ABC.md)
turns the stacked label matrix into a consensus dissimilarity using the
C++ co-clustering kernel.

## Usage

``` r
ABCpp.SingleInMultiple(
  data,
  distmeasure = "euclidean",
  weighting = "var",
  normalize = NULL,
  gr = c(),
  bag = TRUE,
  numsim = 1000,
  numvar = 100,
  linkage = "ward.D2",
  NC = NULL,
  mds = FALSE,
  nfeat = NULL,
  nugget_type = "between",
  nugget_args = list()
)
```

## Arguments

- data:

  Numeric matrix with **objects (samples) in rows** and features in
  columns (the MosaiClusteR convention; transposed internally to the
  features-by-objects layout of the original engine).

- distmeasure:

  Distance measure passed to
  [`Distance`](https://lucp12891.github.io/MosaiClusteR/reference/Distance.md).

- weighting:

  `"var"`, `"cv"`, `"nugget"` or equal.

- normalize:

  Optional normalisation (see
  [`Normalization`](https://lucp12891.github.io/MosaiClusteR/reference/Normalization.md)).

- gr:

  Unused; kept for interface compatibility.

- bag:

  Logical; bootstrap the objects each iteration (Step 2).

- numsim:

  Number of iterations \\R\\.

- numvar:

  Kept for backward compatibility; superseded by `nfeat`.

- linkage:

  Agglomeration method for
  [`fastcluster::hclust`](https://rdrr.io/pkg/fastcluster/man/hclust.html).

- NC:

  Base-cluster count (`NULL` = \\\lfloor\sqrt{N}\rfloor\\).

- mds:

  Unused placeholder.

- nfeat:

  Optional feature-subsample size. `NULL` (default) uses the classic
  \\g=\lfloor\sqrt{G}\rfloor\\; a value \\\ge 1\\ sets an absolute count
  (raise it for low-dimensional sources where most features matter).

- nugget_type, nugget_args:

  Passed to the data-nugget weighting when `weighting = "nugget"`.

## Value

A data frame, `numsim` rows by one column per object; entries are
base-cluster labels (`0` = object not selected that iteration).

## Feature weighting (incl. data nuggets)

Features are subsampled with probability proportional to a
rank-transformed importance score:

- `"var"` – feature variance (the classic ABC choice);

- `"cv"` – coefficient of variation;

- `"nugget"` – *NEW*: the data-nugget between-group weighted variance
  ([`create_data_nuggets`](https://lucp12891.github.io/MosaiClusteR/reference/create_data_nuggets.md)),
  a robust, Big-Data-friendly alternative to the raw variance;

- anything else – equal weights.

## References

Amaratunga D, Cabrera J, Kovtun V (2008). “Microarray Learning with
ABC.” *Biostatistics*, **9**(1), 128–136.
[doi:10.1093/biostatistics/kxm017](https://doi.org/10.1093/biostatistics/kxm017)
.
<https://academic.oup.com/biostatistics/article-abstract/9/1/128/254198>.

## See also

[`M_ABC`](https://lucp12891.github.io/MosaiClusteR/reference/M_ABC.md),
[`create_data_nuggets`](https://lucp12891.github.io/MosaiClusteR/reference/create_data_nuggets.md)

## Examples

``` r
data(mosaic_toy)
lab <- ABCpp.SingleInMultiple(mosaic_toy$List[[1]], numsim = 40, NC = 3,
                              weighting = "nugget")
dim(lab)
#> [1] 40 60
```
