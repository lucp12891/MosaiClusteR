# Compute a distance matrix

Dispatcher computing an object-by-object distance matrix under a range
of measures suitable for continuous, binary and mixed-type data.

## Usage

``` r
Distance(
  Data,
  distmeasure = c("tanimoto", "jaccard", "euclidean", "hamming", "cont tanimoto",
    "MCA_coord", "gower", "chi.squared", "cosine"),
  normalize = NULL,
  method = NULL
)

Distance_v2(Data, distmeasure = "tanimoto", normalize = NULL)
```

## Arguments

- Data:

  A data matrix; rows are the objects.

- distmeasure:

  One of `"tanimoto"`, `"jaccard"`, `"euclidean"` (true Euclidean, via
  [`stats::dist`](https://rdrr.io/r/stats/dist.html)), `"gower"`
  (mixed-type, via
  [`FD::gowdis`](https://rdrr.io/pkg/FD/man/gowdis.html)), `"hamming"`,
  `"cont tanimoto"`, `"MCA_coord"`, `"chi.squared"`, `"cosine"`.

- normalize:

  Optional normalisation method (see
  [`Normalization`](https://lucp12891.github.io/MosaiClusteR/reference/Normalization.md)),
  applied for the `"euclidean"` / `"cont tanimoto"` measures.

- method:

  Ignored; retained for backward compatibility.

## Value

A symmetric numeric distance matrix.

## Details

`"cosine"`, `"chi.squared"` and `"MCA_coord"` require the suggested
packages lsa, analogue and FactoMineR respectively.

## Examples

``` r
m <- matrix(rbinom(40, 1, 0.5), nrow = 8)
D <- Distance(m, distmeasure = "tanimoto")
dim(D)
#> [1] 8 8
```
