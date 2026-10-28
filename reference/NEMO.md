# NEMO: neighborhood-based multi-omics clustering

NEMO (Wang et al. 2014) integrates several omics by building a
locally-scaled, k-nearest-neighbour affinity network for each one and
averaging them, for every pair of objects, over only the omics in which
both objects are measured. Partial / missing modalities are therefore
handled without imputation. The integrated affinity is clustered
spectrally.

## Usage

``` r
NEMO(List, k = NULL, NN = 20, normalize = TRUE)
```

## Arguments

- List:

  A list of data matrices over the same objects, **objects in rows**.
  Sources may cover different subsets of objects (matched by row name)
  for the partial-data setting.

- k:

  Number of clusters (`NULL` = eigengap estimate).

- NN:

  Number of neighbours for the affinity graphs and local scaling.

- normalize:

  Logical; row-normalise each affinity to a transition matrix before
  integration (NEMO default `TRUE`).

## Value

A list of class `"NEMO"` with `FusedM` (integrated affinity), `DistM`
(`1 - FusedM`), `cluster` (labels) and `Clust` (an `agnes` object on
`DistM` for interoperability).

## References

(Wang et al. 2014) ; Rappoport, N. & Shamir, R. (2019). NEMO: cancer
subtyping by integration of partial multi-omic data. Bioinformatics
35:3348-3356.

## See also

[`SNF`](https://lucp12891.github.io/MosaiClusteR/reference/SNF.md),
[`spectral_clustering`](https://lucp12891.github.io/MosaiClusteR/reference/spectral_clustering.md)

## Examples

``` r
data(mosaic_toy)
fit <- NEMO(mosaic_toy$List, k = 3, NN = 15)
table(fit$cluster)
#> 
#>  1  2  3 
#> 20 20 20 
```
