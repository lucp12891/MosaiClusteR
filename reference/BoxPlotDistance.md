# Box plots of one distance matrix categorized against another distance matrix.

Given two distance matrices, the function categorizes one distance
matrix and produces a box plot from the other distance matrix against
the created categories. The option is available to choose one of the
plots or to have both plots. The function also works on outputs from
ADEC and CEC functions which do not have distance matrices but incidence
matrices.

## Usage

``` r
BoxPlotDistance(
  Data1,
  Data2,
  type = c("data", "dist", "clusters"),
  distmeasure = c("tanimoto", "tanimoto"),
  normalize = c(FALSE, FALSE),
  method = c(NULL, NULL),
  lab1,
  lab2,
  limits1 = NULL,
  limits2 = NULL,
  plot = 1,
  StopRange = FALSE,
  plottype = "new",
  location = NULL
)
```

## Arguments

- Data1:

  The first data matrix, cluster outcome or distance matrix to be
  plotted.

- Data2:

  The second data matrix, cluster outcome or distance matrix to be
  plotted.

- type:

  indicates whether the provided matrices in "List" are either data
  matrices, distance matrices or clustering results obtained from the
  data. If type="dist" the calculation of the distance matrices is
  skipped and if type="clusters" the single source clustering is
  skipped. Type should be one of "data", "dist" or "clusters".

- distmeasure:

  A vector of the distance measures to be used on each data matrix.
  Should be one of "tanimoto", "euclidean", "jaccard", "hamming".
  Defaults to c("tanimoto","tanimoto").

- normalize:

  Logical. Indicates whether to normalize the distance matrices or not,
  defaults to c(FALSE, FALSE) for two data sets. This is recommended if
  different distance types are used. More details on normalization in
  `Normalization`.

- method:

  A method of normalization. Should be one of "Quantile","Fisher-Yates",
  "standardize","Range" or any of the first letters of these names.
  Default is c(NULL,NULL) for two data sets.

- lab1:

  The label to plot for Data1.

- lab2:

  The label to plot for Data2.

- limits1:

  The limits for the categories of Data1.

- limits2:

  The limits for the categories of Data2.

- plot:

  The type of plots: 1 - Plot the values of Data1 versus the categories
  of Data2. 2 - Plot the values of Data2 versus the categories of Data1.
  3 - Plot both types 1 and 2.

- StopRange:

  Logical. Indicates whether the distance matrices with values not
  between zero and one should be standardized to have so. If FALSE the
  range normalization is performed. See `Normalization`. If TRUE, the
  distance matrices are not changed. This is recommended if different
  types of data are used such that these are comparable. Default is
  FALSE.

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

One/multiple box plots.

## Examples

``` r
if (FALSE) { # \dontrun{
data(fingerprintMat)
data(targetMat)

MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)


BoxPlotDistance(MCF7_F,MCF7_T,type="cluster",lab1="FP",lab2="TP",limits1=c(0.3,0.7),
limits2=c(0.3,0.7),plot=1,StopRange=FALSE,plottype="new", location=NULL)
} # }
```
