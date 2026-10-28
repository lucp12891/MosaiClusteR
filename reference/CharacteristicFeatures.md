# Determine the characteristic features of clusters

`CharacteristicFeatures` identifies the features that are characteristic
for the objects of each cluster of one or more clustering results.
Binary features are evaluated with Fisher's exact test, continuous
features with the t-test. Multiplicity correction is included.

## Usage

``` r
CharacteristicFeatures(
  List,
  Selection = NULL,
  binData = NULL,
  contData = NULL,
  datanames = NULL,
  nrclusters = NULL,
  sign = 0.05,
  topChar = NULL,
  fusionsLog = TRUE,
  weightclust = TRUE,
  names = NULL
)
```

## Arguments

- List:

  A list of clustering outputs to be compared. The first element of the
  list is used as the reference in `ReorderToReference`.

- Selection:

  If a selection of objects is provided, the function only investigates
  the features of this selection. Default is NULL.

- binData:

  A list of the binary feature data matrices. Default is NULL.

- contData:

  A list of continuous feature data matrices. Default is NULL.

- datanames:

  A vector with the names of the data matrices. Default is NULL.

- nrclusters:

  The number of clusters to cut the dendrogram in. Default is NULL.

- sign:

  The significance level. Default is 0.05.

- topChar:

  The number of top characteristics to return. If NULL, only the
  significant characteristics are saved. Default is NULL.

- fusionsLog:

  Logical. Indicator for the fusion of clusters. Default is TRUE.

- weightclust:

  Logical. For the outputs of CEC, WeightedClust or WeightedSimClust, if
  TRUE only the result of the Clust element is considered. Default is
  TRUE.

- names:

  Optional. Names of the methods. Default is NULL.

## Value

A list with for each method the found (top) characteristics of the
feature data per cluster.

## References

Perualila-Tan N. et al. (2016). Weighted similarity-based clustering of
chemical structures and bioactivity data in early drug discovery.

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
MCF7_Char=CharacteristicFeatures(List=L,Selection=NULL,binData=list(fingerprintMat),
contData=NULL,datanames=c("FP"),nrclusters=7,topChar=20)
} # }
```
