# Comparison of clustering results for the single and multiple source clustering.

The function `CompareSvsM` plots the comparison of the single source
clustering results on the left and that of the multiple source
clustering results on the right such that a visual comparison is
possible.

## Usage

``` r
CompareSvsM(
  ListS,
  ListM,
  nrclusters = NULL,
  cols = NULL,
  fusionsLogS = FALSE,
  fusionsLogM = FALSE,
  weightclustS = FALSE,
  weightclustM = FALSE,
  namesS = NULL,
  namesM = NULL,
  margins = c(8.1, 3.1, 3.1, 4.1),
  plottype = "new",
  location = NULL
)
```

## Arguments

- ListS:

  A list of the outputs from the single source clusterings to be
  compared.

- ListM:

  A list of the outputs from the multiple source clusterings to be
  compared.

- nrclusters:

  The number of clusters to cut the dendrogram in. Default is NULL.

- cols:

  A character vector with the names of the colours. Default is NULL.

- fusionsLogS:

  The fusionslog parameter for the elements in ListS. Default is FALSE.

- fusionsLogM:

  The fusionsLog parameter for the elements in ListM. Default is FALSE.

- weightclustS:

  The weightclust parameter for the elements in ListS. Default is FALSE.

- weightclustM:

  The weightclust parameter for the elements in ListM. Default is FALSE.

- namesS:

  Optional. Names of the single source clusterings. Default is NULL.

- namesM:

  Optional. Names of the multiple source clusterings. Default is NULL.

- margins:

  Optional. Margins to be used for the plot. Default is
  c(8.1,3.1,3.1,4.1).

- plottype:

  Should be one of "pdf","new" or "sweave". Default is "new".

- location:

  If plottype is "pdf", a location should be provided here. Default is
  NULL.

## Value

A plot with on the left the comparison over the objects in ListS and on
the right a comparison over the objects in ListM.

## Examples

``` r
if (FALSE) { # \dontrun{
data(fingerprintMat)
data(targetMat)
data(Colors1)

MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)

L=list(fingerprintMat,targetMat)

MCF7_W=WeightedClust(List=L,type="data", distmeasure=c("tanimoto","tanimoto"),
normalize=c(FALSE,FALSE),method=c(NULL,NULL),weight=seq(1,0,-0.1),weightclust=0.5
,clust="agnes",linkage="ward",StopRange=FALSE)

ListM=list(MCF7_W)
namesM=seq(1.0,0.0,-0.1)

ListS=list(MCF7_F,MCF7_T)
namesS=c("FP","TP")

CompareSvsM(ListS,ListM,nrclusters=7,cols=Colors1,fusionsLogS=FALSE,
fusionsLogM=FALSE,weightclustS=FALSE,weightclustM=FALSE,namesS,
namesM,plottype="new",location=NULL)
} # }
```
