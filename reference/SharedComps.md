# Objects shared across clusterings

`SharedComps` determines, per cluster, which objects are shared across
all listed clustering results. The clusterings are first rearranged to a
common reference with `ReorderToReference` and the intersection of the
lead objects of every method is taken per cluster.

## Usage

``` r
SharedComps(
  List,
  nrclusters = NULL,
  fusionsLog = FALSE,
  weightclust = FALSE,
  names = NULL
)
```

## Arguments

- List:

  A list of the clustering outputs to be compared. The first element of
  the list will be used as the reference in `ReorderToReference`.

- nrclusters:

  If List is the output of several clustering methods, it has to be
  provided in how many clusters to cut the dendrograms in. Default is
  NULL.

- fusionsLog:

  Logical. To be handed to `ReorderToReference`: indicator for the
  fusion of clusters. Default is FALSE.

- weightclust:

  Logical. To be handed to `ReorderToReference`: to be used for the
  outputs of CEC, WeightedClust or WeightedSimClust. If TRUE, only the
  result of the Clust element is considered. Default is FALSE.

- names:

  Names of the methods or clusters. Default is NULL.

## Value

A list containing the shared objects of all listed elements per cluster.

## Examples

``` r
if (FALSE) { # \dontrun{
data(mosaic_toy)
L <- mosaic_toy$List
res <- list(Cluster(L[[1]], type = "data", distmeasure = "euclidean"),
            Cluster(L[[2]], type = "data", distmeasure = "euclidean"))
Comps <- SharedComps(List = res, nrclusters = 5, fusionsLog = FALSE,
                     weightclust = FALSE, names = c("S1", "S2"))
} # }
```
