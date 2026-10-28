# Determine the characteristic features of a cluster or selection

`FeatSelection` determines the characteristic features of a selection of
objects or of a specific cluster across the methods. Binary features are
evaluated with Fisher's exact test, continuous features with the t-test.

## Usage

``` r
FeatSelection(
  List,
  Selection = NULL,
  binData = NULL,
  contData = NULL,
  datanames = NULL,
  nrclusters = NULL,
  topChar = NULL,
  sign = 0.05,
  fusionsLog = TRUE,
  weightclust = TRUE
)
```

## Arguments

- List:

  A list of clustering outputs. Default is NULL.

- Selection:

  Either a character vector of objects, or a numeric cluster number.
  Default is NULL.

- binData:

  A list of the binary feature data matrices. Default is NULL.

- contData:

  A list of continuous feature data matrices. Default is NULL.

- datanames:

  A vector with the names of the data matrices. Default is NULL.

- nrclusters:

  The number of clusters to cut the dendrogram in. Default is NULL.

- topChar:

  The number of top characteristics to return. Default is NULL.

- sign:

  The significance level. Default is 0.05.

- fusionsLog:

  Logical. Indicator for the fusion of clusters. Default is TRUE.

- weightclust:

  Logical. For the outputs of CEC, WeightedClust or WeightedSimClust, if
  TRUE only the result of the Clust element is considered. Default is
  TRUE.

## Value

A list with the found (top) characteristics of the feature data.

## Examples

``` r
if (FALSE) { # \dontrun{
data(fingerprintMat)
data(targetMat)

MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)

L=list(MCF7_F)
MCF7_FeatSel=FeatSelection(List=L,Selection=4,binData=list(fingerprintMat),
datanames=c("FP"),nrclusters=7,topChar=20)
} # }
```
