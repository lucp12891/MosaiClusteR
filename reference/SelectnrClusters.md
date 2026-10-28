# Select the number of clusters via silhouette widths

For each element of `List`, partitions around medoids
([`cluster::pam`](https://rdrr.io/pkg/cluster/man/pam.html)) are
computed over a range of cluster numbers and the average silhouette
width is recorded. The number of clusters maximising the average
silhouette width is reported per source and for the across-source
average. The analytical result (the silhouette table and the chosen
number of clusters) is returned and is testable without producing any
plot.

## Usage

``` r
SelectnrClusters(
  List,
  type = c("data", "dist", "pam"),
  distmeasure = c("tanimoto", "tanimoto"),
  normalize = c(FALSE, FALSE),
  method = c(NULL, NULL),
  nrclusters = seq(5, 25, 1),
  names = NULL,
  StopRange = FALSE,
  plottype = "none",
  location = NULL
)
```

## Arguments

- List:

  A list of data matrices (objects in rows), distance matrices, or
  [`cluster::pam`](https://rdrr.io/pkg/cluster/man/pam.html) outputs,
  depending on `type`.

- type:

  One of `"data"` (compute a distance first via
  [`Distance`](https://lucp12891.github.io/MosaiClusteR/reference/Distance.md)),
  `"dist"` (elements are dissimilarity matrices) or `"pam"` (elements
  are `pam` objects).

- distmeasure:

  Distance measure(s) used when `type = "data"`.

- normalize, method:

  Normalisation controls passed to
  [`Distance`](https://lucp12891.github.io/MosaiClusteR/reference/Distance.md).

- nrclusters:

  Sequence of cluster numbers to evaluate.

- names:

  Optional names for the sources; used in output labels.

- StopRange:

  Logical; if `FALSE` and dissimilarities fall outside `[0, 1]` they are
  range-normalised for comparability.

- plottype:

  One of `"none"` (default; no plot, safe for headless runs), `"new"`,
  `"pdf"` or `"sweave"`. Any value other than `"none"` attempts to draw
  the silhouette curve; drawing is guarded so that a missing graphics
  device does not abort the computation.

- location:

  File location (without extension) used when `plottype = "pdf"`.
  Default `NULL`.

## Value

A list with two elements: `Silhoutte_Widths` (a matrix of average
silhouette widths per source and their mean) and
`Optimal_Nr_of_CLusters` (a one-row data frame with the silhouette
optimal number of clusters per source and overall).

## References

Kaufman, L. and Rousseeuw, P. J. (1990). Finding Groups in Data: An
Introduction to Cluster Analysis. Wiley, New York.

## Examples

``` r
if (FALSE) { # \dontrun{
data(mosaic_toy)
sel <- SelectnrClusters(mosaic_toy$List, type = "data",
                        distmeasure = c("euclidean", "euclidean"),
                        nrclusters = seq(2, 5), plottype = "none")
sel$Optimal_Nr_of_CLusters
} # }
```
