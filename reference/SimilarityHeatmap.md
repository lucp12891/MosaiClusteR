# A heatmap of similarity values between objects

The function `SimilarityHeatmap` plots the similarity values between
objects. The darker the shade, the more similar objects are. The option
is available to set a cutoff value to highlight the most similar
objects.

## Usage

``` r
SimilarityHeatmap(
  Data,
  type = c("data", "clust", "sim", "dist"),
  distmeasure = "tanimoto",
  normalize = FALSE,
  method = NULL,
  linkage = "flexible",
  cutoff = NULL,
  percentile = FALSE,
  plottype = "new",
  location = NULL
)
```

## Arguments

- Data:

  The data of which a heatmap should be drawn.

- type:

  indicates whether the provided matrices in "List" are either data
  matrices, distance matrices or clustering results obtained from the
  data. If type="dist" the calculation of the distance matrices is
  skipped and if type="clusters" the single source clustering is
  skipped. Type should be one of "data", "dist" ,"sim" or "clusters".

- distmeasure:

  The distance measure. Should be one of "tanimoto", "euclidean",
  "jaccard", "hamming". Defaults to "tanimoto".

- normalize:

  Logical. Indicates whether to normalize the distance matrices or not,
  defaults to c(FALSE, FALSE) for two data sets. This is recommended if
  different distance types are used. More details on normalization in
  `Normalization`.

- method:

  A method of normalization. Should be one of "Quantile","Fisher-Yates",
  "standardize","Range" or any of the first letters of these names.
  Default is NULL.

- linkage:

  Choice of inter group dissimilarity (character). Defaults to
  "flexible".

- cutoff:

  Optional. If a cutoff value is specified, all values lower are put to
  zero while all other values are kept. This helps to highlight the most
  similar objects. Default is NULL.

- percentile:

  Logical. The cutoff value can be a percentile. If one want the cutoff
  value to be the 90th percentile of the data, one should specify cutoff
  = 0.90 and percentile = TRUE. Default is FALSE.

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

A heatmap with the names of the objects on the right and bottom and a
dendrogram of the clustering at the left and top.

## Details

If data is of type "clust", the distance matrix is extracted from the
result and transformed to a similarity matrix. Possibly a range
normalization is performed. If data is of type "dist", it is also
transformed to a similarity matrix and cluster is performed on the
distances. If data is of type "sim", the data is tranformed to a
distance matrix on which clustering is performed. Once the similarity
mattrix is obtained, the cutoff value is applied and a heatmap is drawn.
If no cutoff value is desired, one can leave the default NULL
specification.

## Examples

``` r

if (FALSE) { # \dontrun{
data(fingerprintMat)

MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55)

SimilarityHeatmap(Data=MCF7_F,type="clust",cutoff=0.90,percentile=TRUE)
SimilarityHeatmap(Data=MCF7_F,type="clust",cutoff=0.75,percentile=FALSE)

} # }
```
