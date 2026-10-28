# Two-dimensional cluster visualisation (PCA / UMAP / t-SNE)

Project the objects to two dimensions and colour them by cluster.
`"pca"` uses classical MDS of the fused dissimilarity (always
available); `"umap"` needs uwot and `"tsne"` needs Rtsne.

## Usage

``` r
embedding_plot(
  x,
  labels,
  method = c("pca", "umap", "tsne"),
  col = NULL,
  main = NULL
)
```

## Arguments

- x:

  A fitted result (`$DistM`), a dissimilarity, or a `List`.

- labels:

  Cluster label of each object.

- method:

  `"pca"` (default), `"umap"` or `"tsne"`.

- col:

  Optional cluster colours.

- main:

  Plot title.

## Value

Invisibly the 2-column coordinate matrix.

## See also

[`cluster_selection_plot`](https://lucp12891.github.io/MosaiClusteR/reference/cluster_selection_plot.md)

## Examples

``` r
data(mosaic_toy)
fit <- M_ABC(mosaic_toy$List, numsim = 40, NC = 3)
#> M_ABC: source 1 of 2 (weighting = var) ...
#> M_ABC: source 2 of 2 (weighting = var) ...
#> M_ABC: building consensus dissimilarity (C++ kernel) ...
embedding_plot(fit, mosaic_labels(fit, 3), method = "pca")
```
