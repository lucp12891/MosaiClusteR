# Evidence accumulation clustering (co-association)

`EvidenceAccumulation` builds a co-association (similarity) matrix from
the individual source partitions and extracts a consensus partition by
graph partitioning. Three graph-partitioning options are available: a
minimum-spanning-tree procedure (`"MST"`, requires the igraph package),
a thresholded single-link procedure (`"SL"`), and an agnes single-link
procedure (`"SL_agnes"`).

## Usage

``` r
EvidenceAccumulation(
  List,
  type = c("data", "dist", "clust"),
  distmeasure = c("tanimoto", "tanimoto"),
  normalize = c(FALSE, FALSE),
  method = c(NULL, NULL),
  clust = "agnes",
  linkage = c("flexible", "flexible"),
  alpha = 0.625,
  nrclusters = c(7, 7),
  gap = FALSE,
  maxK = 15,
  graphPartitioning = c("MST", "SL", "SL_agnes"),
  t = NULL
)
```

## Arguments

- List:

  A list of data matrices. It is assumed the rows correspond with the
  objects.

- type:

  Indicates whether the provided matrices in "List" are either data
  matrices ("data"), distance matrices ("dist") or clustering results
  ("clust").

- distmeasure:

  A vector of the distance measures to be used on each data matrix.
  Should be one of "tanimoto", "euclidean", "jaccard", "hamming".
  Defaults to c("tanimoto","tanimoto").

- normalize:

  Logical. Indicates whether to normalize the distance matrices or not,
  defaults to c(FALSE, FALSE).

- method:

  A method of normalization. Should be one of "Quantile","Fisher-Yates",
  "standardize","Range" or any of the first letters of these names.
  Default is c(NULL,NULL).

- clust:

  Choice of clustering function (character). Defaults to "agnes".

- linkage:

  Choice of inter group dissimilarity (character) for each data set.
  Defaults to c("flexible","flexible").

- alpha:

  The parameter alpha to be used in the "flexible" linkage of the agnes
  function. Defaults to 0.625.

- nrclusters:

  The number of clusters to divide each individual dendrogram in.
  Default is c(7,7).

- gap:

  Logical. Whether the optimal number of clusters should be determined
  with the gap statistic. Defaults to FALSE.

- maxK:

  The maximal number of clusters to investigate in the gap statistic.
  Default is 15.

- graphPartitioning:

  The graph partitioning method to be performed: "MST", "SL" or
  "SL_agnes".

- t:

  A threshold used by the "MST" and "SL" partitioning procedures.
  Defaults to NULL.

## Value

A list of two elements:

- DistM:

  A NULL object

- Clust:

  The resulting clustering

The value has class 'Ensemble' (or 'Single Clustering' for "SL_agnes").

## Details

The co-association matrix `S` records the fraction of source partitions
in which two objects are co-clustered. `"MST"` prunes a minimum spanning
tree of the similarity graph (optionally at threshold `t`); `"SL"`
performs a thresholded single-link grouping; and `"SL_agnes"` runs agnes
single linkage on the dissimilarity `1-S`. The `"MST"` path requires the
suggested package igraph.

## References

Fred A.L.N. and Jain A.K. (2005). Combining multiple clusterings using
evidence accumulation. IEEE Transactions on Pattern Analysis and Machine
Intelligence, 27(6), 835-850.

## Examples

``` r
if (FALSE) { # \dontrun{
data(mosaic_toy)
L <- mosaic_toy$List
res <- EvidenceAccumulation(List = L, type = "data",
  distmeasure = c("euclidean","euclidean"), nrclusters = c(5,5),
  graphPartitioning = "SL_agnes")
} # }
```
