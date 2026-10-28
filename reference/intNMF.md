# Integrative NMF clustering (intNMF)

Joint non-negative matrix factorisation across data sources: a single
shared basis matrix \\W\\ (objects \\\times\\ rank) and source-specific
coefficient matrices \\H_k\\ minimise \\\sum_k \alpha_k \lVert X_k - W
H_k \rVert_F^2\\ with all factors non-negative. Sources are weighted by
\\\alpha_k = \max_j \overline{SS}\_j / \overline{SS}\_k\\. Object
clusters are obtained by consensus of the dominant-factor assignment
across restarts.

## Usage

``` r
intNMF(List, k, nstart = 10, max_iter = 200, seed = NULL)
```

## Arguments

- List:

  A list of data matrices over the same objects, **objects in rows**.
  Matrices are shifted to be non-negative internally.

- k:

  Number of clusters (also the factorisation rank).

- nstart:

  Number of random restarts forming the consensus.

- max_iter:

  Maximum multiplicative-update iterations per restart.

- seed:

  Optional RNG seed.

## Value

A list of class `"intNMF"` with `cluster`, the consensus `DistM`, an
`agnes` `Clust`, the best-fit basis `W` and source weights `alpha`.

## References

Chalise, P. & Fridley, B. L. (2017). Integrative clustering of
multi-level omic data based on non-negative matrix factorization. PLoS
ONE 12:e0176278. Reviewed in Zhang et al. (2022), WIREs Comput Stat
14:e1553.

## See also

[`M_ABC`](https://lucp12891.github.io/MosaiClusteR/reference/M_ABC.md),
[`SNF`](https://lucp12891.github.io/MosaiClusteR/reference/SNF.md)

## Examples

``` r
data(mosaic_toy)
fit <- intNMF(mosaic_toy$List, k = 3, nstart = 5, seed = 1)
table(fit$cluster)
#> 
#>  1  2  3 
#> 20 20 20 
```
