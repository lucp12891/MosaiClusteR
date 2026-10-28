# Single-source base clustering

Cluster one data (or distance) matrix with agglomerative hierarchical
clustering. This is the Layer-2 baseline primitive of the MoSaIC
framework: run it on each tile on its own to see what each modality
recovers before integrating. It also underpins
[`HierarchicalEnsembleClustering`](https://lucp12891.github.io/MosaiClusteR/reference/HierarchicalEnsembleClustering.md).

## Usage

``` r
Cluster(
  Data,
  type = c("data", "dist"),
  distmeasure = "euclidean",
  normalize = FALSE,
  method = NULL,
  clust = "agnes",
  linkage = "ward",
  alpha = 0.625,
  StopRange = TRUE,
  ...
)
```

## Arguments

- Data:

  A data matrix (objects in rows) or a precomputed distance matrix.

- type:

  `"data"` to compute a distance first, or `"dist"` if `Data` is already
  a dissimilarity matrix.

- distmeasure:

  Distance measure (see
  [`Distance`](https://lucp12891.github.io/MosaiClusteR/reference/Distance.md)),
  used when `type = "data"`.

- normalize:

  Feature normalisation applied before the distance is computed. Either
  a logical (`FALSE`, the default, for none) or a
  [`Normalization`](https://lucp12891.github.io/MosaiClusteR/reference/Normalization.md)
  method name, one of `"Quantile"`, `"Fisher-Yates"`, `"Standardize"` or
  `"Range"`. Note that
  [`Distance`](https://lucp12891.github.io/MosaiClusteR/reference/Distance.md)
  applies feature normalisation only for `distmeasure = "euclidean"`.

- method:

  Normalisation method used when `normalize = TRUE` (default
  `"Standardize"`). Ignored when `normalize` is `FALSE` or already names
  a method.

- clust:

  Clustering backend; currently `"agnes"`.

- linkage:

  Agglomeration method. Default `"ward"`.

- alpha:

  Flexible-linkage parameter for
  [`cluster::agnes`](https://rdrr.io/pkg/cluster/man/agnes.html).

- StopRange:

  Logical; if `FALSE` and the dissimilarities fall outside `[0, 1]` they
  are range-normalised for comparability.

- ...:

  Ignored; accepted for backward compatibility (e.g. `gap`, `maxK`).

## Value

A list with `DistM` (the dissimilarity matrix) and `Clust` (an `agnes`
object).

## See also

[`Distance`](https://lucp12891.github.io/MosaiClusteR/reference/Distance.md),
[`M_ABC`](https://lucp12891.github.io/MosaiClusteR/reference/M_ABC.md)

## Examples

``` r
data(mosaic_toy)
# objects are in rows, so no transpose is needed
base1 <- Cluster(mosaic_toy$List[[1]], type = "data",
                 distmeasure = "euclidean")
mosaic_labels(base1, k = 3)[1:6]
#> obj01 obj02 obj03 obj04 obj05 obj06 
#>     1     1     1     1     1     1 
```
