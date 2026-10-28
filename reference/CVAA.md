# Cumulative Voting Aggregation (CVAA / W-CVAA)

`CVAA` performs cumulative voting consensus clustering. A reference
partition is used to align the individual source partitions, which are
then aggregated into a consensus stochastic membership matrix. Two
variants are available: cumulative voting (`"CVAA"`) and the
information-weighted variant (`"W-CVAA"`).

## Usage

``` r
CVAA(
  Reference = NULL,
  nrclustersR = 7,
  List,
  typeL = c("data", "dist", "clust"),
  distmeasure = c("tanimoto", "tanimoto"),
  normalize = c(FALSE, FALSE),
  method = c(NULL, NULL),
  clust = "agnes",
  linkage = c("flexible", "flexible"),
  alpha = 0.625,
  nrclusters = c(7, 7),
  gap = FALSE,
  maxK = 15,
  votingMethod = c("CVAA", "W-CVAA"),
  optimalk = nrclustersR
)
```

## Arguments

- Reference:

  A reference partition (a clustering result, optionally an object
  produced by a weighted/CEC method). Required.

- nrclustersR:

  The number of clusters of the reference partition. Default is 7.

- List:

  A list of data matrices. It is assumed the rows correspond with the
  objects.

- typeL:

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

- votingMethod:

  The method to be performed: "CVAA" or "W-CVAA".

- optimalk:

  An estimate of the final optimal number of clusters. Defaults to
  nrclustersR.

## Value

A list of two elements:

- DistM:

  A NULL object

- Clust:

  The resulting clustering (order, order.lab and Clusters)

The value has class 'Ensemble'.

## Details

Cumulative voting aggregation aligns each source membership matrix to a
running consensus by least-squares voting, then averages cumulatively.
The weighted variant orders the partitions by average information
content and weights their contribution accordingly. When `optimalk` is
smaller than `nrclustersR` the JS-ALink agglomeration is used to merge
reference clusters.

## References

Ayad H.G. and Kamel M.S. (2008). Cumulative voting consensus method for
partitions with variable number of clusters. IEEE Transactions on
Pattern Analysis and Machine Intelligence, 30(1), 160-173.

## Examples

``` r
if (FALSE) { # \dontrun{
data(mosaic_toy)
L <- mosaic_toy$List
Ref <- Cluster(L[[1]], type = "data", distmeasure = "euclidean")
attr(Ref, "method") <- "Single"
res <- CVAA(Reference = Ref, nrclustersR = 5, List = L, typeL = "data",
  distmeasure = c("euclidean","euclidean"), votingMethod = "CVAA", optimalk = 5)
} # }
```
