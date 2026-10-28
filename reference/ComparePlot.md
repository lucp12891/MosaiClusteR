# Comparison of clustering results over multiple methods

A visual comparison of the clustering results of several methods. The
function relies on `ReorderToReference` and `ColorsNames` and renders
the resulting matrix either as a rectangular heatmap (using the
`plotrix` package) or in a circular format (using the `circlize`
package).

## Usage

``` r
ComparePlot(
  List,
  nrclusters = NULL,
  cols = NULL,
  fusionsLog = FALSE,
  weightclust = FALSE,
  names = NULL,
  margins = c(8.1, 3.1, 3.1, 4.1),
  circle = FALSE,
  canvaslims = c(-1, 1, -1, 1),
  Highlight = NULL,
  substr = NULL,
  cex.highlight = 1,
  cex.labels = 0.7,
  trackheightdend = 0.5,
  plottype = "new",
  location = NULL
)
```

## Arguments

- List:

  A list of the outputs from the methods to be compared. The first
  element of the list is used as the reference in `ReorderToReference`.

- nrclusters:

  The number of clusters to cut the dendrogram in. Default is NULL.

- cols:

  A character vector with the colours to be used. Default is NULL.

- fusionsLog:

  Logical. To be handed to `ReorderToReference`. Default is FALSE.

- weightclust:

  Logical. To be handed to `ReorderToReference`. Default is FALSE.

- names:

  Optional. Names of the methods to be used as labels. Default is NULL.

- margins:

  Optional. Margins for the plot. Default is c(8.1,3.1,3.1,4.1).

- circle:

  Logical. Whether the figure should be circular (TRUE) or a rectangle
  (FALSE). Default is FALSE.

- canvaslims:

  The limits for the circular dendrogram. Default is
  c(-1.0,1.0,-1.0,1.0).

- Highlight:

  Optional. A list of character vectors of objects to be highlighted.
  Default is NULL.

- substr:

  Optional. A vector of length two giving start and stop positions to
  shorten labels. Default is NULL.

- cex.highlight:

  Magnification for the highlight text. Default is 1.

- cex.labels:

  Magnification for the labels. Default is 0.7.

- trackheightdend:

  The height of the dendrogram track in the circular plot. Default is
  0.5.

- plottype:

  Should be one of "pdf","new" or "sweave". Default is "new".

- location:

  Optional. If plottype is "pdf", a location should be provided. Default
  is NULL.

## Value

A plot which translates the matrix output of `ReorderToReference`.
