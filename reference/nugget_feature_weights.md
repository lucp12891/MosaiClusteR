# Feature weights derived from data nuggets

Computes a per-feature importance score from a data-nugget
([`create_data_nuggets`](https://lucp12891.github.io/MosaiClusteR/reference/create_data_nuggets.md))
object built over the *observations* (samples). This score is the
data-nugget analogue of the feature variance used by the ABC / M-ABC
family: it rewards features that separate the representative groups
while discounting within-nugget noise, giving a robust,
Big-Data-friendly weighting.

## Usage

``` r
nugget_feature_weights(nuggets, type = c("between", "ratio", "inv_scale"))
```

## Arguments

- nuggets:

  A `"data_nugget"` object from
  [`create_data_nuggets`](https://lucp12891.github.io/MosaiClusteR/reference/create_data_nuggets.md)
  (built with samples in the rows).

- type:

  One of:

  `"between"`

  : (default) size-weighted between-nugget variance of each feature:
    \\\sum_i w_i (c\_{ij} - \bar c_j)^2 / \sum_i w_i\\.

  `"ratio"`

  : a separability / pseudo-\\F\\ score, \\\mathrm{between}\_j /
    (\mathrm{between}\_j + \overline{\mathrm{within}}\_j)\\.

  `"inv_scale"`

  : inverse mean within-nugget variance, which down-weights noisy
    features.

## Value

A named numeric vector of length `ncol(centers)` of non-negative feature
scores.

## References

Cherasia KE, Cabrera J, Fernholz LT, Fernholz R (2023). “Data Nuggets in
Supervised Learning.” In *Robust and Multivariate Statistical Methods:
Festschrift in Honor of David E. Tyler*, 429–449. Springer.

## See also

[`create_data_nuggets`](https://lucp12891.github.io/MosaiClusteR/reference/create_data_nuggets.md),
[`M_ABC`](https://lucp12891.github.io/MosaiClusteR/reference/M_ABC.md)

## Examples

``` r
set.seed(1)
x <- cbind(signal = c(rnorm(50, 0), rnorm(50, 5)),  # separates groups
           noise  = rnorm(100))                      # pure noise
dn <- create_data_nuggets(x, max_nuggets = 20)
nugget_feature_weights(dn, type = "between")
#>    signal     noise 
#> 6.9017623 0.8824925 
```
