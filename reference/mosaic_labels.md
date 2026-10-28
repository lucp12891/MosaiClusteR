# Extract a flat partition from a MosaiClusteR result

Cuts the hierarchical clustering carried in a MosaiClusteR result object
into `k` groups. Works with any method returning a `Clust` element (an
`agnes`/`hclust` object) and/or a `DistM`.

## Usage

``` r
mosaic_labels(object, k, linkage = "ward")
```

## Arguments

- object:

  A MosaiClusteR result (e.g. from
  [`M_ABC`](https://lucp12891.github.io/MosaiClusteR/reference/M_ABC.md)),
  an `agnes`/`hclust` object, or a distance matrix.

- k:

  Number of clusters: a single positive integer, at most the number of
  objects.

- linkage:

  Linkage used when `object` is a bare distance matrix.

## Value

An integer vector of cluster labels, named by object.

## Examples

``` r
set.seed(1)
s <- paste0("o", 1:20)
L <- list(matrix(rnorm(20 * 40), 20, dimnames = list(s, NULL)),
          matrix(rnorm(20 * 50), 20, dimnames = list(s, NULL)))
fit <- M_ABC(L, numsim = 30)
#> M_ABC: source 1 of 2 (weighting = var) ...
#> M_ABC: source 2 of 2 (weighting = var) ...
#> M_ABC: building consensus dissimilarity (C++ kernel) ...
mosaic_labels(fit, k = 3)
#>  o1  o2  o3  o4  o5  o6  o7  o8  o9 o10 o11 o12 o13 o14 o15 o16 o17 o18 o19 o20 
#>   1   2   1   3   1   1   3   1   2   2   2   2   2   3   1   1   2   3   2   2 
```
