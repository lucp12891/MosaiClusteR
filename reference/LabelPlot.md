# Coloring specific leaves of a dendrogram

The function plots a dendrogrmam of which specific leaves are coloured.

## Usage

``` r
LabelPlot(Data, sel1, sel2 = NULL, col1 = NULL, col2 = NULL)
```

## Arguments

- Data:

  The result of a method which contains the dendrogram to be colored.

- sel1:

  The selection of objects to be colored. Default is NULL.

- sel2:

  An optional second selection to be colored. Default is NULL.

- col1:

  The color for the first selection. Default is NULL.

- col2:

  The color for the optional second selection. Default is NULL.

## Value

A plot of the dendrogram of which the leaves of the selection(s) are
colored.

## Examples

``` r
if (FALSE) { # \dontrun{
data(fingerprintMat)
MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)

ClustF_6=cutree(MCF7_F$Clust,6)

SelF=rownames(fingerprintMat)[ClustF_6==6]
SelF

LabelPlot(Data=MCF7_F,sel1=SelF,sel2=NULL,col1='darkorchid')
} # }
```
