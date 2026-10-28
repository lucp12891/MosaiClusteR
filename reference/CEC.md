# Complementary ensemble clustering

Complementary Ensemble Clustering (CEC) Complementary Ensemble
Clustering (CEC, Fodeh2013) shows similarities with ADEC. However,
instead of merging the data matrices, ensemble clustering is performedon
each data matrix separately. The resulting incidence matrices for each
data sets are combined in a weighted linear equation. The weighted
incidence matrix is the input for the final clustering algorithm.
Similarly as ADEC, there are versions depending of the specification of
the number of features to sample and the number of clusters.

## Usage

``` r
CEC(
  List,
  distmeasure = c("tanimoto", "tanimoto"),
  normalize = c(FALSE, FALSE),
  method = c(NULL, NULL),
  t = 10,
  r = NULL,
  nrclusters = NULL,
  weight = NULL,
  clust = "agnes",
  linkage = c("flexible", "flexible"),
  alpha = 0.625,
  weightclust = 0.5
)
```

## Arguments

- List:

  A list of data matrices. It is assumed the rows are corresponding with
  the objects.

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

- t:

  The number of iterations. Defaults to 10.

- r:

  A vector with the number of features to take for the random sample for
  each element in List. If NULL (default), all features are considered.

- nrclusters:

  A list with a sequence of numbers of clusters to cut the dendrogram in
  for each element in List. If NULL (default), the function stops.

- weight:

  The weights for the weighted linear combination.

- clust:

  Choice of clustering function (character). Defaults to "agnes".

- linkage:

  Choice of inter group dissimilarity (character) for each data set.
  Defaults to c("flexible", "flexible") for two data sets.

- alpha:

  The parameter alpha to be used in the "flexible" linkage of the agnes
  function. Defaults to 0.625 and is only used if the linkage is set to
  "flexible"

- weightclust:

  A weight for which the result will be put aside of the other results.
  This was done for comparative reason and easy access.

## Value

The returned value is a list of four elements:

- DistM:

  The resulting incidence matrix

- Results:

  The hierarchical clustering result for each element in WeightedDist

- Clust:

  The result for the weight specified in Clustweight

The value has class 'CEC'.

## Details

If r is specified and nrclusters is a fixed number, only a random
sampling of the features will be performed for the t iterations (CECa).
If r is NULL and the nrclusters is a sequence, the clustering is
performedon all features and the dendrogam is divided into clusters for
the values of nrclusters (CECb). If both r is specified and nrclusters
is a sequence, the combination is performed (CECc). After every
iteration, either be random sampling, multiple divisions of the
dendrogram or both, an incidence matrix is set up. All incidence
matrices are summed and represent the distance matrix on which a final
clustering is performed.

## References

Fodeh SJ, Punch W, Tan P (2013). “On Unifying Multi-source Information
to Improve Subcellular Protein Localization Prediction.” *Proteins:
Structure, Function, and Bioinformatics*, **81**, 1298–1310.
