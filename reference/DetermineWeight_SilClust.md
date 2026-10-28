# Determines an optimal weight for weighted clustering by silhouettes widths.

The function `DetermineWeight_SilClust` determines an optimal weight for
weighted similarity clustering by calculating silhouettes widths. For
each given weight, a linear combination of the distance matrices of the
single data sources is obtained, medoid clustering is performed and the
silhouette widths are regressed against the cluster memberships. A
statistic is derived and a p-value obtained via bootstrapping.

## Usage

``` r
DetermineWeight_SilClust(
  List,
  type = c("data", "dist", "clusters"),
  distmeasure = c("tanimoto", "tanimoto"),
  normalize = c(FALSE, FALSE),
  method = c(NULL, NULL),
  weight = seq(0, 1, by = 0.01),
  nrclusters = NULL,
  names = NULL,
  nboot = 10,
  StopRange = FALSE,
  plottype = "new",
  location = NULL
)
```

## Arguments

- List:

  A list of matrices of the same type. It is assumed the rows are
  corresponding with the objects.

- type:

  indicates whether the provided matrices in "List" are either data
  matrices, distance matrices or clustering results obtained from the
  data. Type should be one of "data", "dist" or "clusters".

- distmeasure:

  A vector of the distance measures to be used on each data matrix.
  Defaults to c("tanimoto","tanimoto").

- normalize:

  Logical. Indicates whether to normalize the distance matrices. Default
  is c(FALSE, FALSE).

- method:

  A method of normalization. Default is c(NULL,NULL).

- weight:

  Optional. A list of different weight combinations for the data sets in
  List. Defaults to seq(0,1,by=0.01).

- nrclusters:

  The number of clusters to cut the dendrogram in. Default is NULL.

- names:

  The labels to give to the elements in List. Default is NULL.

- nboot:

  Number of bootstraps to be run. Default is 10.

- StopRange:

  Logical. Indicates whether the distance matrices with values not
  between zero and one should be standardized to have so. Default is
  FALSE.

- plottype:

  Should be one of "pdf","new" or "sweave". Default is "new".

- location:

  If plottype is "pdf", a location should be provided here. Default is
  NULL.

## Value

Two plots and a list with two elements:

- Result:

  A data frame with the statistic for each weight combination

- Weight:

  The optimal weight

.

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

MC7_Weight=DetermineWeight_SilClust(List=L,type="clusters",distmeasure=
c("tanimoto","tanimoto"),normalize=c(FALSE,FALSE),method=c(NULL,NULL),
weight=seq(0,1,by=0.01),nrclusters=c(7,7),names=c("FP","TP"),nboot=10,
StopRange=FALSE,plottype="new",location=NULL)
} # }
```
