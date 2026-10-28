# Aggregated data ensemble clustering

Aggregated Data Ensemble Clustering (ADEC) is a direct clustering
multi-source technique. ADEC is an iterative procedure which starts with
the merging of the data sets. In each iteration, a random sample of the
features is selected and/or a resulting dendrogram is divided into k
clusters for a range of values of k.

## Usage

``` r
ADEC(
  List,
  distmeasure = "tanimoto",
  normalize = FALSE,
  method = NULL,
  t = 10,
  r = NULL,
  nrclusters = NULL,
  clust = "agnes",
  linkage = "flexible",
  alpha = 0.625,
  seed = NULL
)
```

## Arguments

- List:

  A list of data matrices of the same type. It is assumed the rows are
  corresponding with the objects.

- distmeasure:

  Choice of metric for the dissimilarity matrix (character). Should be
  one of "tanimoto", "euclidean", "jaccard", "hamming". Defaults to
  "tanimoto".

- normalize:

  Logical. Indicates whether to normalize the distance matrices or not,
  defaults to FALSE. This is recommended if different distance types are
  used. More details on normalization in `Normalization`.

- method:

  A method of normalization. Should be one of "Quantile","Fisher-Yates",
  "standardize","Range" or any of the first letters of these names.
  Default is NULL.

- t:

  The number of iterations. Defaults to 10.

- r:

  The number of features to take for the random sample. If NULL
  (default), all features are considered.

- nrclusters:

  A sequence of numbers of clusters to cut the dendrogram in. If NULL
  (default), the function stops.

- clust:

  Choice of clustering function (character). Defaults to "agnes".

- linkage:

  Choice of inter group dissimilarity (character). Defaults to
  "flexible".

- alpha:

  The parameter alpha to be used in the "flexible" linkage of the agnes
  function. Defaults to 0.625 and is only used if the linkage is set to
  "flexible".

- seed:

  Optional integer RNG seed. ADEC is stochastic when `r` is set (it
  takes a random sample of the features each iteration), so set `seed`
  for reproducible results.

## Value

The returned value is a list with the following three elements.

- AllData:

  Fused data matrix of the data matrices

- DistM:

  The resulting co-association matrix

- Clust:

  The resulting clustering

The value has class 'ADEC'. The Clust element will be of interest for
further applications.

## Details

If r is specified and nrclusters is a fixed number, only a random
sampling of the features will be performed for the t iterations (ADECa).
If r is NULL and the nrclusters is a sequence, the clustering is
performed on all features and the dendrogam is divided into clusters for
the values of nrclusters (ADECb). If both r is specified and nrclusters
is a sequence, the combination is performed (ADECc). After every
iteration, either be random sampling, multiple divisions of the
dendrogram or both, an incidence matrix is set up. All incidence
matrices are summed and represent the distance matrix on which a final
clustering is performed.

## References

Fodeh, S. J., Punch, W. & Tan, P.-N. (2013). On unifying multi-source
information from probabilistically generated functional data. Proteins
81:1298-1310.

## Examples

``` r
if (FALSE) { # \dontrun{
data(mosaic_toy)
L <- mosaic_toy$List
res_ADEC <- ADEC(List=L, distmeasure="euclidean", normalize=FALSE, method=NULL,
                 t=5, r=50, nrclusters=5, clust="agnes", linkage="flexible", alpha=0.625)
} # }
```
