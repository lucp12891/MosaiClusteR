# Find a selection of objects in the output of `ReorderToReference`

`FindCluster` selects the objects belonging to a cluster after the
results of the methods have been rearranged by `ReorderToReference`.

## Usage

``` r
FindCluster(
  List,
  nrclusters = NULL,
  select = c(1, 1),
  fusionsLog = TRUE,
  weightclust = TRUE,
  names = NULL
)
```

## Arguments

- List:

  A list of the clustering outputs to be compared. The first element of
  the list will be used as the reference in `ReorderToReference`.

- nrclusters:

  The number of clusters to cut the dendrogram in. Default is NULL.

- select:

  The row (the method) and the number of the cluster to select. Default
  is c(1,1).

- fusionsLog:

  Logical. To be handed to `ReorderToReference`: indicator for the
  fusion of clusters. Default is TRUE.

- weightclust:

  Logical. To be handed to `ReorderToReference`: to be used for the
  outputs of CEC, WeightedClust or WeightedSimClust. If TRUE, only the
  result of the Clust element is considered. Default is TRUE.

- names:

  Optional. Names of the methods. Default is NULL.

## Value

A character vector containing the names of the objects in the selected
cluster.

## Examples

``` r
if (FALSE) { # \dontrun{
data(mosaic_toy)
L <- mosaic_toy$List
res <- list(Cluster(L[[1]], type = "data", distmeasure = "euclidean"),
            Cluster(L[[2]], type = "data", distmeasure = "euclidean"))
Comps <- FindCluster(List = res, nrclusters = 5, select = c(1, 2),
                     weightclust = FALSE)
} # }
```
