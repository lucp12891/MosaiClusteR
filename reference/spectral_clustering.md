# Spectral clustering of an affinity matrix

Normalised (Ng-Jordan-Weiss) spectral clustering. Given a symmetric,
non-negative affinity matrix it forms the normalised Laplacian, embeds
the objects in the space of its leading eigenvectors and applies
k-means.

## Usage

``` r
spectral_clustering(affinity, k = NULL, max_k = 10, nstart = 10)
```

## Arguments

- affinity:

  A symmetric \\n \times n\\ non-negative similarity/affinity matrix
  over the objects.

- k:

  Number of clusters. If `NULL`, chosen by the largest eigengap of the
  normalised Laplacian (searched over `2:max_k`).

- max_k:

  Upper bound for the eigengap search when `k = NULL`.

- nstart:

  Number of k-means restarts.

## Value

A list with `cluster` (named integer labels), `k` and the spectral
`embedding`.

## References

Ng, A. Y., Jordan, M. I. & Weiss, Y. (2002). On spectral clustering:
analysis and an algorithm. NIPS 14.

## See also

[`NEMO`](https://lucp12891.github.io/MosaiClusteR/reference/NEMO.md),
[`SNF`](https://lucp12891.github.io/MosaiClusteR/reference/SNF.md)

## Examples

``` r
set.seed(1)
X <- rbind(matrix(rnorm(20 * 5, 0), 20), matrix(rnorm(20 * 5, 6), 20))
A <- exp(-as.matrix(dist(X))^2 / 2)
spectral_clustering(A, k = 2)$cluster[1:5]
#> 1 2 3 4 5 
#> 1 1 1 1 1 
```
