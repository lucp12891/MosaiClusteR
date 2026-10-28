# Plot of continuous features

The function `ContFeaturesPlot` plots the values of continuous features.
It is possible to separate between objects of interest and the other
objects.

## Usage

``` r
ContFeaturesPlot(
  leadCpds,
  data,
  nrclusters = NULL,
  orderLab = NULL,
  colorLab = NULL,
  cols = NULL,
  ylab = "features",
  addLegend = TRUE,
  margins = c(5.5, 3.5, 0.5, 8.7),
  plottype = "new",
  location = NULL
)
```

## Arguments

- leadCpds:

  A character vector containing the objects one wants to separate from
  the others.

- data:

  The data matrix.

- nrclusters:

  Optional. The number of clusters to consider if colorLab is specified.
  Default is NULL.

- orderLab:

  Optional. If the objects are to set in a specific order of a specific
  method. Default is NULL.

- colorLab:

  The clustering result that determines the color of the labels of the
  objects in the plot. If NULL, the labels are black. Default is NULL.

- cols:

  The colors for the labels of the objects. Default is NULL.

- ylab:

  The lable of the y-axis. Default is "features".

- addLegend:

  Logical. Indicates whether a legend should be added to the plot.
  Default is TRUE.

- margins:

  Optional. Margins to be used for the plot. Default is
  c(5.5,3.5,0.5,8.7).

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

A plot in which the values of the features of the leadCpds are separeted
from the others.

## Examples

``` r
if (FALSE) { # \dontrun{
data(Colors1)
Comps=c("Cpd1", "Cpd2", "Cpd3", "Cpd4", "Cpd5")

Data=matrix(sample(15, size = 50*5, replace = TRUE), nrow = 50, ncol = 5)
colnames(Data)=colnames(Data, do.NULL = FALSE, prefix = "col")
rownames(Data)=rownames(Data, do.NULL = FALSE, prefix = "row")
for(i in 1:50){
  rownames(Data)[i]=paste("Cpd",i,sep="")
}

ContFeaturesPlot(leadCpds=Comps,orderLab=rownames(Data),colorLab=NULL,data=Data,
nrclusters=7,cols=Colors1,ylab="features",addLegend=TRUE,margins=c(5.5,3.5,0.5,8.7),
plottype="new",location=NULL)
} # }
```
