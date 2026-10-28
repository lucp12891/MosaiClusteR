# List all features present in a selected cluster of objects

The function `FeaturesOfCluster` lists the number of features objects of
the cluster have in common. A threshold can be set selecting among how
many objects of the cluster the features should be shared. An optional
plot of the features is available.

## Usage

``` r
FeaturesOfCluster(
  leadCpds,
  data,
  threshold = 1,
  plot = TRUE,
  cexB = 1.1,
  cexR = 1.1,
  cexN = 1.1,
  plottype = "new",
  location = NULL
)
```

## Arguments

- leadCpds:

  A character vector containing the objects one wants to investigate in
  terms of features.

- data:

  The data matrix.

- threshold:

  The number of objects the features at least should be shared amongst.
  Default is set to 1.

- plot:

  Logical. Indicates whether or not a BinFeaturesPlot should be set up.
  Default is TRUE.

- cexB:

  Font size of the bottom (object) labels in the plot. Default is 1.1.

- cexR:

  Font size of the right (feature) labels in the plot. Default is 1.1.

- cexN:

  Font size of the data-name label in the plot. Default is 1.1.

- plottype:

  Should be one of "pdf","new" or "sweave". Default is "new".

- location:

  If plottype is "pdf", a location should be provided. Default is NULL.

## Value

A list with 2 items: the number of shared features among the objects,
and a character vector of the plotted features.

## Examples

``` r
if (FALSE) { # \dontrun{
data(fingerprintMat)

Lead=rownames(fingerprintMat)[1:5]

FeaturesOfCluster(leadCpds=Lead,data=fingerprintMat,
threshold=1,plot=TRUE,plottype="new",location=NULL)
} # }
```
