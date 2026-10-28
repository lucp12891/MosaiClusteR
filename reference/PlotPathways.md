# A GO plot of a pathway analysis output.

The `PlotPathways` function takes an output of the `PathwayAnalysis`
function and plots a GO graph with the help of the `plotGOgraph`
function of the MLP package.

## Usage

``` r
PlotPathways(
  Pathways,
  nRow = 5,
  main = NULL,
  plottype = "new",
  location = NULL
)
```

## Arguments

- Pathways:

  One element of the output list returned by `PathwayAnalysis` or
  `Geneset.intersect`.

- nRow:

  Number of GO IDs for which to produce the plot. Default is 5.

- main:

  Title of the plot. Default is NULL.

- plottype:

  Should be one of "pdf","new" or "sweave". If "pdf", a location should
  be provided in "location" and the figure is saved there. If "new" a
  new graphic device is opened and if "sweave", the figure is made
  compatible to appear in a sweave or knitr document. Default is "new".

- location:

  If plottype is "pdf", a location should be provided in "location" and
  the figure is saved there. Default is NULL.

## Value

The output is a GO graph.

## Examples

``` r
if (FALSE) { # \dontrun{
data(fingerprintMat)
data(targetMat)
data(geneMat)
data(GeneInfo)

MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)

L=list(MCF7_F,MCF7_T)
names=c('FP','TP')

MCF7_PathsFandT=PathwayAnalysis(List=L, geneExpr = geneMat, nrclusters = 7, method = c("limma",
"MLP"), geneInfo = GeneInfo, geneSetSource = "GOBP", topP = NULL,
topG = NULL, GENESET = NULL, sign = 0.05,niter=2,fusionsLog = TRUE, weightclust = TRUE,
 names =names,seperatetables=FALSE,separatepvals=FALSE)

PlotPathways(MCF7_PathsFandT$FP$"Cluster 1"$Pathways,nRow=5,main=NULL)
} # }
```
