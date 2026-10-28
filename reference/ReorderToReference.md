# Relabel and reorder clusterings to a reference

`ReorderToReference` relabels and reorders the clusters of a series of
clustering results so that the cluster numbers of every method are
matched, via a Gale-Shapley stable-matching scheme, to those of the
first (reference) method. The result is a matrix in which the columns
represent the objects (ordered as for the reference method) and the rows
represent the methods. Each cell contains the number of the cluster the
object belongs to for that method, relabelled to the reference. It is a
possibility that two or more clusters are fused together compared to the
reference; in that case the function alerts the user and asks to set
`fusionsLog` to `TRUE`.

## Usage

``` r
ReorderToReference(
  List,
  nrclusters = NULL,
  fusionsLog = FALSE,
  weightclust = FALSE,
  names = NULL
)
```

## Arguments

- List:

  A list of clustering outputs to be compared. The first element of the
  list will be used as the reference.

- nrclusters:

  The number of clusters to cut the dendrogram in. Default is NULL.

- fusionsLog:

  Logical. Indicator for the fusion of clusters. Default is FALSE.

- weightclust:

  Logical. To be used for the outputs of CEC, WeightedClust or
  WeightedSimClust. If TRUE, only the result of the Clust element is
  considered. Default is FALSE.

- names:

  Optional. A character vector with the names of the methods.

## Value

A matrix of which the cells indicate to what cluster the objects belong
to according to the rearranged methods.

## Note

The `ReorderToReference` function was optimized for the situations
presented by the data sets at hand. It is noted that the function might
fail in a particular situation which results in an infinite loop.

## References

Van Moerbeke, M. (2017). MoSaiClusteR: integrative clustering of
multi-source data. Hasselt University.

## Examples

``` r
if (FALSE) { # \dontrun{
data(mosaic_toy)
L <- mosaic_toy$List
res <- list(Cluster(L[[1]], type = "data", distmeasure = "euclidean"),
            Cluster(L[[2]], type = "data", distmeasure = "euclidean"))
M <- ReorderToReference(res, nrclusters = 5, fusionsLog = FALSE,
                        weightclust = FALSE, names = c("S1", "S2"))
} # }
```
