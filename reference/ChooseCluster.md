# Interactive plot to determine DE Genes and DE features for a specific cluster

If desired, the function produces a dendrogram of a clustering result.
One or multiple clusters can be indicated by a mouse click. From these
clusters DE genes and characteristic features are determined. It is also
possible to provide the objects of interest without producing the plot.

## Usage

``` r
ChooseCluster(
  Interactive = TRUE,
  leadCpds = NULL,
  clusterResult = NULL,
  colorLab = NULL,
  binData = NULL,
  contData = NULL,
  datanames = c("FP"),
  geneExpr = NULL,
  topChar = 20,
  topG = 20,
  sign = 0.05,
  nrclusters = NULL,
  cols = NULL,
  n = 1
)
```

## Arguments

- Interactive:

  Logical. Whether an interactive plot should be made. Defaults to TRUE.

- leadCpds:

  A list of the objects of the clusters of interest. Defaults to NULL.

- clusterResult:

  The output of one of the aggregated cluster functions. Default is
  NULL.

- colorLab:

  The clustering result the dendrogram should be colored after. Default
  is NULL.

- binData:

  A list of the binary feature data matrices. Default is NULL.

- contData:

  A list of continuous data sets of the objects. Default is NULL.

- datanames:

  A vector with the names of the data matrices. Default is c("FP").

- geneExpr:

  A gene expression matrix, may also be an ExpressionSet. Default is
  NULL.

- topChar:

  The number of top characteristics to return. Default is 20.

- topG:

  The number of top genes to return. Default is 20.

- sign:

  The significance level. Default is 0.05.

- nrclusters:

  Optional. The number of clusters to cut the dendrogram in. Default is
  NULL.

- cols:

  The colors to use in the dendrogram. Default is NULL.

- n:

  The number of clusters one wants to identify by a mouse click. Default
  is 1.

## Value

A list with one element per cluster of interest with elements objects,
Characteristics and (optionally) Genes.

## Details

This function may require the suggested packages 'a4Base', 'a4Core' and
'limma' when a gene expression matrix is supplied.

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

MCF7_Interactive=ChooseCluster(Interactive=TRUE,leadCpds=NULL,clusterResult=MCF7_T,
colorLab=MCF7_F,binData=list(fingerprintMat),datanames=c("FP"),geneExpr=geneMat,
topChar = 20, topG = 20,nrclusters=7,n=1)
} # }
```
