# Aggregated data clustering

Aggregated Data Clustering (ADC) is a direct clustering multi-source
technique. The data matrices are merged into a single fused data matrix
on which a dissimilarity matrix is computed. Hierarchical clustering is
then performed on this dissimilarity matrix.

## Usage

``` r
ADC(
  List,
  distmeasure = "tanimoto",
  normalize = FALSE,
  method = NULL,
  clust = "agnes",
  linkage = "flexible",
  alpha = 0.625
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

- clust:

  Choice of clustering function (character). Defaults to "agnes".

- linkage:

  Choice of inter group dissimilarity (character). Defaults to
  "flexible".

- alpha:

  The parameter alpha to be used in the "flexible" linkage of the agnes
  function. Defaults to 0.625 and is only used if the linkage is set to
  "flexible".

## Value

The returned value is a list with the following three elements.

- AllData:

  Fused data matrix of the data matrices

- DistM:

  The resulting distance matrix

- Clust:

  The resulting clustering

The value has class 'ADC'. The Clust element will be of interest for
further applications.

## References

Fodeh, S. J., Punch, W. & Tan, P.-N. (2013). On unifying multi-source
information from probabilistically generated functional data. Proteins
81:1298-1310.

## Examples

``` r
if (FALSE) { # \dontrun{
data(mosaic_toy)
L <- mosaic_toy$List
res_ADC <- ADC(List=L, distmeasure="euclidean", normalize=FALSE, method=NULL,
               clust="agnes", linkage="flexible", alpha=0.625)
} # }
```
