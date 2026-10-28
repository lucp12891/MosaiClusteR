# Determines an optimal weight for weighted clustering by similarity weighted clustering.

The function `DetermineWeight_SimClust` determines an optimal weight for
performing weighted similarity clustering. For each given weight, each
separate clustering is compared to the clustering on a weighted
dissimilarity matrix and a Jaccard coefficient is calculated. The ratio
of the Jaccard coefficients closest to one indicates an optimal weight.

## Usage

``` r
DetermineWeight_SimClust(
  List,
  type = c("data", "dist", "clusters"),
  distmeasure = c("tanimoto", "tanimoto"),
  normalize = c(FALSE, FALSE),
  method = c(NULL, NULL),
  weight = seq(0, 1, by = 0.01),
  nrclusters = NULL,
  clust = "agnes",
  linkage = c("flexible", "flexible"),
  linkageF = "ward",
  alpha = 0.625,
  gap = FALSE,
  maxK = 15,
  names = NULL,
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

- clust:

  Choice of clustering function (character). Defaults to "agnes".

- linkage:

  Choice of inter group dissimilarity for the individual clusterings.
  Defaults to c("flexible","flexible").

- linkageF:

  Choice of inter group dissimilarity for the final clustering. Defaults
  to "ward".

- alpha:

  The parameter alpha to be used in the "flexible" linkage of the agnes
  function. Defaults to 0.625.

- gap:

  Logical. Whether or not to calculate the gap statistic. Only if
  type="data". Default is FALSE.

- maxK:

  The maximal number of clusters to consider in calculating the gap
  statistic. Default is 15.

- names:

  The labels to give to the elements in List. Default is NULL.

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

The returned value is a list with three elements:

- ClustSep:

  The result of `Cluster` for each single element of List

- Result:

  A data frame with the Jaccard coefficients and their ratios for each
  weight

- Weight:

  The optimal weight

.

## References

Perualila-Tan, N. et al. (2016). Weighted similarity-based clustering of
chemical structures and bioactivity data in early drug discovery.
Journal of Bioinformatics and Computational Biology.

## Examples

``` r
if (FALSE) { # \dontrun{
data(fingerprintMat)
data(targetMat)

MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="flexible",alpha=0.625,gap=FALSE,maxK=55,StopRange=FALSE)
MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="flexible",alpha=0.625,gap=FALSE,maxK=55,StopRange=FALSE)

L=list(MCF7_F,MCF7_T)

MCF7_Weight=DetermineWeight_SimClust(List=L,type="clusters",weight=seq(0,1,by=0.01),
nrclusters=c(7,7),distmeasure=c("tanimoto","tanimoto"),normalize=c(FALSE,FALSE),
method=c(NULL,NULL),clust="agnes",linkage=c("flexible","flexible"),linkageF="ward",
alpha=0.625,gap=FALSE,maxK=50,names=c("FP","TP"),StopRange=FALSE,plottype="new",location=NULL)
} # }
```
