# Track a cluster or a selection of objects across multiple methods

`TrackCluster` follows a selection of objects or a specific cluster
across a sequence of clustering results, recording how the objects are
distributed over the clusters of each method and optionally producing
tracking plots.

## Usage

``` r
TrackCluster(
  List,
  Selection,
  nrclusters = NULL,
  followMaxComps = FALSE,
  followClust = TRUE,
  fusionsLog = TRUE,
  weightclust = TRUE,
  names = NULL,
  selectionPlot = FALSE,
  table = FALSE,
  completeSelectionPlot = FALSE,
  ClusterPlot = FALSE,
  cols = NULL,
  legendposx = 0.5,
  legendposy = 2.4,
  plottype = "sweave",
  location = NULL
)
```

## Arguments

- List:

  A list of clustering outputs to be compared. The first element of the
  list is used as the reference in `ReorderToReference`.

- Selection:

  Either a character vector of objects, or a numeric cluster number.

- nrclusters:

  The number of clusters to cut the dendrogram in. Default is NULL.

- followMaxComps:

  Logical. Whether to follow the cluster with the maximum number of
  objects. Default is FALSE.

- followClust:

  Logical. Whether to follow the cluster of interest. Default is TRUE.

- fusionsLog:

  Logical. Indicator for the fusion of clusters. Default is TRUE.

- weightclust:

  Logical. For the outputs of CEC, WeightedClust or WeightedSimClust, if
  TRUE only the result of the Clust element is considered. Default is
  TRUE.

- names:

  Optional. Names of the methods. Default is NULL.

- selectionPlot:

  Logical. Whether to produce the selection plot. Default is FALSE.

- table:

  Logical. Whether to produce a table of the tracking. Default is FALSE.

- completeSelectionPlot:

  Logical. Whether to produce the complete selection plot. Default is
  FALSE.

- ClusterPlot:

  Logical. Whether to produce the cluster plot. Default is FALSE.

- cols:

  The colors to use in the plots. Default is NULL.

- legendposx:

  The x-coordinate of the legend. Default is 0.5.

- legendposy:

  The y-coordinate of the legend. Default is 2.4.

- plottype:

  Should be one of "pdf","new" or "sweave". Default is "sweave".

- location:

  If plottype is "pdf", a location should be provided. Default is NULL.

## Value

A list with one element per method describing how the selection is
distributed across the clusters, and optionally a tracking table.

## Examples

``` r
if (FALSE) { # \dontrun{
data(fingerprintMat)
data(targetMat)

MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)

L=list(MCF7_F,MCF7_T)
MCF7_Track=TrackCluster(List=L,Selection=4,nrclusters=7,followMaxComps=FALSE,
followClust=TRUE,selectionPlot=TRUE,table=TRUE,cols=NULL)
} # }
```
