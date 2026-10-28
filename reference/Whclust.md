# Weighted hierarchical clustering (weighted Ward)

Agglomerative clustering of weighted observations using the weighted
Ward merge cost. Intended for the small set of data-nugget centres.

## Usage

``` r
Whclust(x, weights = rep(1, nrow(x)))
```

## Arguments

- x:

  Numeric matrix, observations (e.g. nugget centres) in rows.

- weights:

  Non-negative weight per observation.

## Value

An [`hclust`](https://rdrr.io/r/stats/hclust.html) object.

## References

Dey et al. (2025). WCluster.

## See also

[`Wkmeans`](https://lucp12891.github.io/MosaiClusteR/reference/Wkmeans.md),
[`nugget_cluster`](https://lucp12891.github.io/MosaiClusteR/reference/nugget_cluster.md)

## Examples

``` r
set.seed(1)
x <- rbind(matrix(rnorm(8 * 3, 0), 8), matrix(rnorm(8 * 3, 5), 8))
stats::cutree(Whclust(x, weights = rep(2, 16)), k = 2)[1:4]
#> 1 2 3 4 
#> 1 1 1 1 
```
