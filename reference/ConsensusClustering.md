# Voting-based consensus clustering

Consensus Clustering performs voting-based consensus on a set of
single-source partitions. Each data matrix is clustered separately and
the resulting label matrix is reconciled into a single consensus
partition by an iterative voting scheme. Three voting methods are
available: Iterative Voting Consensus (IVC), Iterative Probabilistic
Voting Consensus (IPVC) and Iterative Pairwise Consensus (IPC).

## Usage

``` r
ConsensusClustering(
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
  votingMethod = c("IVC", "IPVC", "IPC"),
  optimalk = 7
)
```

## Arguments

- List:

  A list of data matrices. It is assumed the rows are corresponding with
  the objects.

- type:

  Indicates whether the provided matrices in "List" are either data
  matrices ("data"), distance matrices ("dist") or clustering results
  obtained from the data ("clust").

- distmeasure:

  A vector of the distance measures to be used on each data matrix.
  Should be one of "tanimoto", "euclidean", "jaccard", "hamming".
  Defaults to c("tanimoto","tanimoto").

- normalize:

  Logical. Indicates whether to normalize the distance matrices or not,
  defaults to c(FALSE, FALSE) for two data sets. More details on
  normalization in `Normalization`.

- method:

  A method of normalization. Should be one of "Quantile","Fisher-Yates",
  "standardize","Range" or any of the first letters of these names.
  Default is c(NULL,NULL) for two data sets.

- clust:

  Choice of clustering function (character). Defaults to "agnes".

- linkage:

  Choice of inter group dissimilarity (character) for each data set.
  Defaults to c("flexible", "flexible") for two data sets.

- alpha:

  The parameter alpha to be used in the "flexible" linkage of the agnes
  function. Defaults to 0.625 and is only used if the linkage is set to
  "flexible".

- nrclusters:

  The number of clusters to divide each individual dendrogram in.
  Default is c(7,7) for two data sets.

- gap:

  Logical. Whether the optimal number of clusters should be determined
  with the gap statistic. Defaults to FALSE.

- maxK:

  The maximal number of clusters to investigate in the gap statistic.
  Default is 15.

- votingMethod:

  The voting method to be performed: one of "IVC", "IPVC" or "IPC".

- optimalk:

  The number of clusters for the final consensus partition. Defaults to
  7.

## Value

The returned value is a list of two elements:

- DistM:

  A NULL object

- Clust:

  The resulting clustering

The value has class 'Ensemble'.

## References

Nguyen N. and Caruana R. (2007). Consensus Clusterings. Seventh IEEE
International Conference on Data Mining (ICDM 2007), 607-612.

## Examples

``` r
if (FALSE) { # \dontrun{
data(mosaic_toy)
fit <- ConsensusClustering(List = mosaic_toy$List, type = "data",
  distmeasure = c("euclidean", "euclidean"), normalize = c(FALSE, FALSE),
  method = c(NULL, NULL), clust = "agnes", linkage = c("ward", "ward"),
  nrclusters = c(3, 3), gap = FALSE, votingMethod = "IVC", optimalk = 3)
} # }
```
