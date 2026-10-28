# Data-nugget clustering

Cluster a (potentially large) data set by first compressing its
observations into weighted data nuggets
([`create_data_nuggets`](https://lucp12891.github.io/MosaiClusteR/reference/create_data_nuggets.md)),
clustering the nugget centres with a weighted clusterer, and mapping the
nugget labels back to the original observations. This is the clustering
analogue of the data-nugget weighting used by
[`M_ABC`](https://lucp12891.github.io/MosaiClusteR/reference/M_ABC.md)
and follows `DN.Wkmeans` / `DN.Whclust` from the WCluster package.

## Usage

``` r
nugget_cluster(
  data,
  k,
  method = c("kmeans", "ward"),
  max_nuggets = NULL,
  nugget_args = list()
)
```

## Arguments

- data:

  A numeric matrix (objects in rows) *or* a precomputed
  [`create_data_nuggets`](https://lucp12891.github.io/MosaiClusteR/reference/create_data_nuggets.md)`{data_nugget}`
  object. For multi-source data, column-bind the per-source matrices
  first (the direct/`ADC` strategy) or pass a single fused matrix.

- k:

  Number of clusters.

- method:

  Weighted clusterer for the nugget centres: `"kmeans"` (WWCSS) or
  `"ward"` (weighted Ward).

- max_nuggets, nugget_args:

  Passed to
  [`create_data_nuggets`](https://lucp12891.github.io/MosaiClusteR/reference/create_data_nuggets.md)
  when `data` is a raw matrix.

## Value

A list of class `"nugget_cluster"` with `cluster` (labels for the
original objects), `nugget_cluster` (labels for the nuggets) and the
`nuggets` object.

## References

Cherasia et al. (2023); Dey et al. (2025), WCluster.

## See also

[`create_data_nuggets`](https://lucp12891.github.io/MosaiClusteR/reference/create_data_nuggets.md),
[`Wkmeans`](https://lucp12891.github.io/MosaiClusteR/reference/Wkmeans.md)

## Examples

``` r
data(mosaic_toy)
X <- do.call(cbind, mosaic_toy$List)        # fuse sources (objects in rows)
fit <- nugget_cluster(X, k = 3, max_nuggets = 30)
table(fit$cluster)
#> 
#>  1  2  3 
#> 20 20 20 
```
