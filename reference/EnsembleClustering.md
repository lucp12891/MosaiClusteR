# Ensemble clustering via graph partitioning (CSPA / HGPA / MCLA)

Cluster-based ensemble clustering that combines several single-source
partitions into a single consensus partition through the
graph-partitioning algorithms of Strehl and Ghosh: the Cluster-based
Similarity Partitioning Algorithm (CSPA), the HyperGraph Partitioning
Algorithm (HGPA) and the Meta-CLustering Algorithm (MCLA), as well as a
"Best" selection strategy.

**External executable required.** This method does not run in pure R. It
writes the stacked label matrix to disk and delegates the actual graph
partitioning to an **external graph-partitioning executable** (e.g. a
compiled CSPA/HGPA/MCLA binary backed by METIS / `gpmetis` / `shmetis`)
or to a **MATLAB** session. These binaries are **not bundled** with the
package; the user must supply and configure them. The path to the
executable is selected through the `executable` argument. When
`executable=FALSE` (the default), or when the configured executable
cannot be located, the function stops with an informative error rather
than attempting to launch an unavailable program.

## Usage

``` r
EnsembleClustering(
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
  ensembleMethod = c("CSPA", "HGPA", "MCLA", "Best"),
  waitingtime = 300,
  file_number = 0,
  executable = FALSE
)
```

## Arguments

- List:

  A list of data matrices. It is assumed the rows are corresponding with
  the objects.

- type:

  indicates whether the provided matrices in "List" are either data
  matrices, distance matrices or clustering results obtained from the
  data. If type="dist" the calculation of the distance matrices is
  skipped and if type="clust" the single source clustering is skipped.
  Type should be one of "data", "dist" or "clust".

- distmeasure:

  A vector of the distance measures to be used on each data matrix.
  Should be one of "tanimoto", "euclidean", "jaccard", "hamming".
  Defaults to c("tanimoto","tanimoto").

- normalize:

  Logical. Indicates whether to normalize the distance matrices or not,
  defaults to c(FALSE, FALSE) for two data sets. This is recommended if
  different distance types are used. More details on normalization in
  `Normalization`.

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
  with the gap statistic. Default is FALSE.

- maxK:

  The maximal number of clusters to investigate in the gap statistic.
  Default is 15.

- ensembleMethod:

  A character string indicating the consensus function to use. One of
  "CSPA", "HGPA", "MCLA" or "Best". Defaults to "CSPA".

- waitingtime:

  The maximum number of seconds to wait for the external program to
  write its result before stopping. Defaults to 300.

- file_number:

  An identifier appended to the temporary files exchanged with the
  external program. Defaults to 00.

- executable:

  Either FALSE (the default) or the path/flag identifying the external
  graph-partitioning executable (or MATLAB) to use. When FALSE the
  function stops, since no consensus can be computed in pure R.
  Precompiled binaries are not bundled with the package.

## Value

The returned value is a list of two elements:

- DistM:

  A NULL object

- Clust:

  The resulting consensus clustering

The value has class 'Ensemble'.

## References

Strehl A. and Ghosh J. (2002). Cluster ensembles - a knowledge reuse
framework for combining multiple partitions. Journal of Machine Learning
Research, 3, 583-617.

## Examples

``` r
if (FALSE) { # \dontrun{
data(mosaic_toy)
L <- mosaic_toy$List

# Requires an external graph-partitioning executable / MATLAB to be
# configured through the 'executable' argument.
res <- EnsembleClustering(List = L, type = "data",
  distmeasure = c("euclidean", "euclidean"), normalize = c(FALSE, FALSE),
  method = c(NULL, NULL), clust = "agnes",
  linkage = c("flexible", "flexible"), nrclusters = c(7, 7),
  gap = FALSE, maxK = 15, ensembleMethod = "CSPA",
  executable = "./EnsembleClusteringC")
} # }
```
