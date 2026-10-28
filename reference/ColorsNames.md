# Function that annotates colors to their names

The `ColorsNames` function is used on the output of the
`ReorderToReference` and matches the cluster numbers indicated by the
cell with the names of the colors. This is necessary to produce the plot
of the `ComparePlot` function and is therefore an internal function of
this function but can also be applied separately.

## Usage

``` r
ColorsNames(matrixColors, cols = NULL)
```

## Arguments

- matrixColors:

  The output of the ReorderToReference function.

- cols:

  A character vector with the names of the colours to be used. Default
  is NULL.

## Value

A matrix containing the hex code of the color that corresponds to each
cell of the matrix to be colored. This function is called upon by the
`ComparePlot` function.

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

L=list(MCF7_F,MCF7_T)
names=c("FP","TP")

MatrixColors=ReorderToReference(List=L,nrclusters=7,fusionsLog=TRUE,weightclust=TRUE,
names=names)

Names=ColorsNames(matrixColors=MatrixColors,cols=Colors1)
} # }
```
