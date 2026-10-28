# Colouring clusters in a dendrogram

Plot a dendrogram with leaves colored by a result of choice.

## Usage

``` r
ClusterPlot(
  Data1,
  Data2 = NULL,
  nrclusters = NULL,
  cols = NULL,
  colorComps = NULL,
  hangdend = 0.02,
  plottype = "new",
  location = NULL,
  ...
)
```

## Arguments

- Data1:

  The resulting clustering of a method which contains the dendrogram to
  be colored.

- Data2:

  Optional. The resulting clustering of another method , i.e. the
  resulting clustering on which the colors should be based. Default is
  NULL.

- nrclusters:

  Optional. The number of clusters to cut the dendrogram in. If not
  specified the dendrogram will be drawn without colours to discern the
  different clusters. Default is NULL.

- cols:

  The colours for the clusters if nrclusters is specified. Default is
  NULL.

- colorComps:

  If only a specific set of objects needs to be highlighted, this can be
  specified here. The objects should be given in a character vector. If
  specified, all other compound labels will be colored black. Default is
  NULL.

- hangdend:

  A specification for the length of the brances of the dendrogram.
  Default is 0.02.

- plottype:

  Should be one of "pdf","new" or "sweave". If "pdf", a location should
  be provided in "location" and the figure is saved there. If "new" a
  new graphic device is opened and if "sweave", the figure is made
  compatible to appear in a sweave or knitr document. Default is "new".

- location:

  If plottype is "pdf", a location should be provided in "location" and
  the figure is saved there. Default is NULL.

- ...:

  Other options which can be given to the plot function.

## Value

A plot of the dendrogram of the first clustering result with colored
leaves. If a second clustering result is given in Data2, the colors are
based on this clustering result.

## Examples

``` r
if (FALSE) { # \dontrun{
data(fingerprintMat)
data(targetMat)
data(Colors1)

MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)

ClusterPlot(MCF7_T ,nrclusters=7,cols=Colors1,plottype="new",location=NULL,
main="Clustering on Target Predictions: Dendrogram",ylim=c(-0.1,1.8))
} # }
```
