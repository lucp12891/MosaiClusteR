# Plotting gene profiles

In `ProfilePlot`, the gene profiles of the significant genes for a
specific cluster are shown on 1 plot. Therefore, each gene is normalized
by subtracting its the mean.

## Usage

``` r
ProfilePlot(
  Genes,
  Comps,
  geneExpr = NULL,
  raw = FALSE,
  orderLab = NULL,
  colorLab = NULL,
  nrclusters = NULL,
  cols = NULL,
  addLegend = TRUE,
  margins = c(8.1, 4.1, 1.1, 6.5),
  extra = 5,
  plottype = "new",
  location = NULL
)
```

## Arguments

- Genes:

  The genes to be plotted.

- Comps:

  The objects to be plotted or to be separated from the other objects.

- geneExpr:

  The gene expression matrix or ExpressionSet of the objects.

- raw:

  Logical. Should raw p-values be plotted? Default is FALSE.

- orderLab:

  Optional. If the objects are to set in a specific order of a specific
  method. Default is NULL.

- colorLab:

  The clustering result that determines the color of the labels of the
  objects in the plot. Default is NULL.

- nrclusters:

  Optional. The number of clusters to cut the dendrogram in.

- cols:

  Optional. The color to use for the objects in Clusters for each
  method.

- addLegend:

  Optional. Whether a legend of the colors should be added to the plot.

- margins:

  Optional. Margins to be used for the plot. Default is
  margins=c(8.1,4.1,1.1,6.5).

- extra:

  The space between the plot and the legend. Default is 5.

- plottype:

  Should be one of "pdf","new" or "sweave". If "pdf", a location should
  be provided in "location" and the figure is saved there. If "new" a
  new graphic device is opened and if "sweave", the figure is made
  compatible to appear in a sweave or knitr document, i.e. no new device
  is opened and the plot appears in the current device or document.
  Default is "new".

- location:

  If plottype is "pdf", a location should be provided in "location" and
  the figure is saved there. Default is NULL.

## Value

A plot which contains multiple gene profiles. A distinction is made
between the values for the objects in Comps and the others.

## Examples

``` r
if (FALSE) { # \dontrun{
data(fingerprintMat)
data(targetMat)
data(geneMat)

MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
  method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
  method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)

L=list(MCF7_F,MCF7_T)
names=c('FP','TP')

MCF7_FT_DE = DiffGenes(List=L,geneExpr=geneMat,nrclusters=7,method="limma",sign=0.05,topG=10,
fusionsLog=TRUE,weightclust=TRUE)

Comps=SharedComps(list(MCF7_FT_DE$`Method 1`$"Cluster 1",MCF7_FT_DE$`Method 2`$"Cluster 1"))[[1]]

MCF7_SharedGenes=FindGenes(dataLimma=MCF7_FT_DE,names=c("FP","TP"))

Genes=names(MCF7_SharedGenes[[1]])[-c(2,4,5)]

colscl=ColorPalette(colors=c("red","green","purple","brown","blue","orange"),ncols=9)

ProfilePlot(Genes=Genes,Comps=Comps,geneExpr=geneMat,raw=FALSE,orderLab=MCF7_F,
colorLab=NULL,nrclusters=7,cols=colscl,addLegend=TRUE,margins=c(16.1,6.1,1.1,13.5),
extra=4,plottype="sweave",location=NULL)
} # }
```
