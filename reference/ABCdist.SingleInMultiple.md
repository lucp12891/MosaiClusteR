# Single-source ABC with distance accumulation

A variant of
[`ABCpp.SingleInMultiple`](https://lucp12891.github.io/MosaiClusteR/reference/ABCpp.SingleInMultiple.md)
that accumulates the mean scaled pairwise distance between co-selected
objects across iterations, yielding a `[0, 1]` dissimilarity matrix that
[`M_ABCdist`](https://lucp12891.github.io/MosaiClusteR/reference/M_ABCdist.md)
fuses across sources.

## Usage

``` r
ABCdist.SingleInMultiple(
  data,
  distmeasure = "euclidean",
  weighting = "var",
  normalize = NULL,
  gr = c(),
  bag = TRUE,
  numsim = 1000,
  numvar = 100,
  linkage = "ward.D",
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

A list with `D` (\\N \times N\\ dissimilarity), `NC` and `ng`
(per-iteration feature count); if `mds = TRUE` also `mds_coords`.

## References

(Amaratunga et al. 2008)

## See also

[`M_ABCdist`](https://lucp12891.github.io/MosaiClusteR/reference/M_ABCdist.md),
[`ABCpp.SingleInMultiple`](https://lucp12891.github.io/MosaiClusteR/reference/ABCpp.SingleInMultiple.md)

## Examples

``` r
data(mosaic_toy)
out <- ABCdist.SingleInMultiple(mosaic_toy$List[[1]], numsim = 40)
dim(out$D)
#> [1] 60 60
```
