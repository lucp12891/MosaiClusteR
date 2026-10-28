# Consensus clustering from stacked ABC label matrices

Turns a row-stacked matrix of per-iteration cluster labels (\\K \cdot
R\\ rows by \\N\\ columns) into a co-clustering dissimilarity and
applies [`cluster::agnes`](https://rdrr.io/pkg/cluster/man/agnes.html).
\\S\_{ij}\\ is the co-clustering probability over iterations selecting
both objects (B6). The co-clustering counts use the compiled C++ kernel
`count_coclustering_cpp` (from `src/mcpp_sampling.cpp`) when available,
otherwise a pure-R loop.

## Usage

``` r
f.clustABC.MultiSource(
  res,
  numclust = NULL,
  dissimilarity = "mabc",
  linkage = "ward",
  alpha = 0.625,
  mds = FALSE
)
```

## Arguments

- res:

  Integer label matrix (`0` = object not selected).

- numclust:

  Optional number of clusters to cut the tree into.

- dissimilarity:

  `"mabc"` (arc-cosine, default) or `"abc"` (linear).

- linkage:

  AGNES agglomeration method.

- alpha:

  Flexible-linkage parameter.

- mds:

  Logical; draw an MDS of the dissimilarities.

## Value

A list with `DistM`, `Clust` and (if `numclust` given) `cut`.

## References

(Amaratunga et al. 2008)

## See also

[`M_ABC`](https://lucp12891.github.io/MosaiClusteR/reference/M_ABC.md)
