# Weighted k-means (minimises WWCSS)

k-means in which each observation carries a non-negative weight. Cluster
centres are weighted means and the objective is the Weighted
Within-Cluster Sum of Squares \\\sum_c \sum\_{i \in c} w_i \lVert x_i -
\mu_c \rVert^2\\.

## Usage

``` r
Wkmeans(x, k, weights = rep(1, nrow(x)), nstart = 10, max_iter = 100)
```

## Arguments

- x:

  Numeric matrix, observations in rows.

- k:

  Number of clusters.

- weights:

  Non-negative weight per observation (default all 1).

- nstart:

  Number of random restarts (best WWCSS kept).

- max_iter:

  Maximum Lloyd iterations per start.

## Value

A list with `cluster`, `centers` and `wwcss`.

## References

Dey, T., Duan, Y., Cabrera, J. & Cheng, X. (2025). WCluster: Clustering
and PCA with Weights, and Data Nuggets Clustering.

## See also

[`nugget_cluster`](https://lucp12891.github.io/MosaiClusteR/reference/nugget_cluster.md),
[`Whclust`](https://lucp12891.github.io/MosaiClusteR/reference/Whclust.md)

## Examples

``` r
set.seed(1)
x <- rbind(matrix(rnorm(20 * 3, 0), 20), matrix(rnorm(20 * 3, 5), 20))
Wkmeans(x, k = 2, weights = runif(40, 1, 5))$cluster[1:5]
#> [1] 1 1 1 1 1
```
