# Interactive comparison of single and multiple source clustering results

The function `CompareInteractive` produces an interactive plot of the
multiple source clustering results in `ListM`. By identifying a cluster
or a method, a comparison with the single source clustering results in
`ListS` is shown in a new plotting window.

## Usage

``` r
CompareInteractive(
  ListM,
  ListS,
  nrclusters = NULL,
  cols = NULL,
  fusionsLogM = FALSE,
  fusionsLogS = FALSE,
  weightclustM = FALSE,
  weightclustS = FALSE,
  namesM = NULL,
  namesS = NULL,
  marginsM = c(2, 2.5, 2, 2.5),
  marginsS = c(8, 2.5, 2, 2.5),
  Interactive = TRUE,
  n = 1
)
```

## Arguments

- ListM:

  A list of the outputs from the multiple source clusterings to be
  compared.

- ListS:

  A list of the outputs from the single source clusterings to be
  compared.

- nrclusters:

  The number of clusters to cut the dendrogram in. Default is NULL.

- cols:

  A character vector with the names of the colours. Default is NULL.

- fusionsLogM:

  The fusionsLog parameter for the elements in ListM. Default is FALSE.

- fusionsLogS:

  The fusionsLog parameter for the elements in ListS. Default is FALSE.

- weightclustM:

  The weightclust parameter for the elements in ListM. Default is FALSE.

- weightclustS:

  The weightclust parameter for the elements in ListS. Default is FALSE.

- namesM:

  Optional. Names of the multiple source clusterings. Default is NULL.

- namesS:

  Optional. Names of the single source clusterings. Default is NULL.

- marginsM:

  Optional. Margins to be used for the ListM plot. Default is
  c(2,2.5,2,2.5).

- marginsS:

  Optional. Margins to be used for the ListS plot. Default is
  c(8,2.5,2,2.5).

- Interactive:

  Optional. Do you want an interactive plot? Default is TRUE.

- n:

  The number of methods/clusters you want to identify. Default is 1.

## Value

A plot of the comparison of the elements of ListM, on which multiple
clusters and/or methods can be identified to compare them to the
elements in ListS.

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

MCF7_W=WeightedClust(List=L,type="data",distmeasure=c("tanimoto","tanimoto"),
normalize=c(FALSE,FALSE),method=c(NULL,NULL),weight=seq(1,0,-0.1),weightclust=0.5,
clust="agnes",linkage="ward",StopRange=FALSE)

ListM=list(MCF7_W)
namesM=c(seq(1.0,0.0,-0.1))

ListS=list(MCF7_F,MCF7_T)
namesS=c("FP","TP")

CompareInteractive(ListM,ListS,nrclusters=7,cols=Colors1,fusionsLogM=FALSE,
fusionsLogS=FALSE,weightclustM=FALSE,weightclustS=TRUE,namesM,namesS,
marginsM=c(2,2.5,2,2.5),marginsS=c(8,2.5,2,2.5),Interactive=TRUE,n=1)
} # }
```
