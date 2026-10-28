# Compares medoid clustering results based on silhouette widths

The function `CompareSilCluster` compares the results of two medoid
clusterings. The null hypothesis is that the clustering is identical. A
test statistic is calculated and a p-value obtained with bootstrapping.

## Usage

``` r
CompareSilCluster(
  List,
  type = c("data", "dist"),
  distmeasure = c("tanimoto", "tanimoto"),
  normalize = c(FALSE, FALSE),
  method = c(NULL, NULL),
  nrclusters = NULL,
  names = NULL,
  nboot = 100,
  plottype = "new",
  location = NULL
)
```

## Arguments

- List:

  A list of data matrices. It is assumed the rows are corresponding with
  the objects.

- type:

  indicates whether the provided matrices in "List" are either data
  matrices or distance matrices. Type should be one of "data" or "dist".

- distmeasure:

  A vector of the distance measures to be used on each data matrix.
  Defaults to c("tanimoto","tanimoto").

- normalize:

  Logical. Indicates whether to normalize the distance matrices. Default
  is c(FALSE, FALSE).

- method:

  A method of normalization. Default is c(NULL,NULL).

- nrclusters:

  The number of clusters to cut the dendrogram in. Default is NULL.

- names:

  The labels to give to the elements in List. Default is NULL.

- nboot:

  Number of bootstraps to be run. Default is 100.

- plottype:

  Should be one of "pdf","new" or "sweave". Default is "new".

- location:

  If plottype is "pdf", a location should be provided here. Default is
  NULL.

## Value

A plot of the density of the statistic under the null hypothesis.
Further a list with two elements:

- Observed Statistic:

  The observed statistical value

- P-Value:

  The P-value of the obtained statistic retrieved after bootstrapping

.

## Examples

``` r
if (FALSE) { # \dontrun{
data(fingerprintMat)
data(targetMat)

List=list(fingerprintMat,targetMat)

Comparison=CompareSilCluster(List=List,type="data",
distmeasure=c("tanimoto","tanimoto"),normalize=c(FALSE,FALSE),method=c(NULL,NULL),
nrclusters=7,names=NULL,nboot=100,plottype="new",location=NULL)

Comparison
} # }
```
