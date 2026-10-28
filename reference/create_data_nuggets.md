# Create data nuggets from a data matrix

`create_data_nuggets` compresses `n` observations into `M <= n`
representative "nuggets". Each nugget stores a centre (location), a
weight (the number of member observations) and a scale (within-nugget
variability). Data nuggets preserve the covariance structure of the
data - including the periphery - far better than a random sub-sample,
which makes them a robust, memory-light surrogate for the full data set
when clustering or computing feature weights.

## Usage

``` r
create_data_nuggets(
  x,
  max_nuggets = NULL,
  min_size = 1L,
  engine = c("native", "datanugget"),
  scale_features = TRUE,
  seed = NULL,
  ...
)
```

## Arguments

- x:

  A numeric matrix or data frame with observations in the *rows* and
  variables in the *columns*.

- max_nuggets:

  Target (maximum) number of nuggets \\M\\. Defaults to a data-driven
  value (see Details). Ignored when `engine = "datanugget"`.

- min_size:

  Minimum number of observations per nugget (native engine only).
  Defaults to `1`.

- engine:

  Either `"native"` (a self-contained k-means based partition, no extra
  dependencies) or `"datanugget"` (delegates to the datanugget package
  via
  [`datanugget::create.DN`](https://rdrr.io/pkg/datanugget/man/create.DN.html)
  /
  [`datanugget::refine.DN`](https://rdrr.io/pkg/datanugget/man/refine.DN.html)
  when it is installed).

- scale_features:

  Logical; if `TRUE` (default) the partition is computed on standardised
  features so that no single high-variance variable dominates the nugget
  assignment. Centres and scales are always reported on the original
  scale.

- seed:

  Optional integer seed for reproducible partitioning.

- ...:

  Additional arguments forwarded to
  [`datanugget::create.DN`](https://rdrr.io/pkg/datanugget/man/create.DN.html)
  when `engine = "datanugget"`.

## Value

An object of class `"data_nugget"`: a list with elements

- centers:

  \\M \times p\\ matrix of nugget centres.

- weights:

  Length-\\M\\ integer vector of nugget sizes \\n_i\\.

- scales:

  \\M \times p\\ matrix of within-nugget variances.

- membership:

  Length-\\n\\ integer vector mapping each observation to its nugget.

- x_scale:

  Length-\\M\\ vector of total within-nugget variability
  \\\mathrm{tr}(\mathrm{Cov})\\.

- engine:

  The engine used.

## Details

With the native engine the number of nuggets defaults to
`min(nrow(x), max(2, floor(nrow(x) / 10)))` when `n > 50`, and to `n`
otherwise (in which case every observation is its own nugget and the
nugget-weighted variance reduces to the ordinary variance - a safe
fallback for small samples). The partition is obtained with
[`stats::kmeans`](https://rdrr.io/r/stats/kmeans.html); missing values
are mean-imputed for the partition step only.

## References

Cherasia KE, Cabrera J, Fernholz LT, Fernholz R (2023). “Data Nuggets in
Supervised Learning.” In *Robust and Multivariate Statistical Methods:
Festschrift in Honor of David E. Tyler*, 429–449. Springer.

Dey T, Duan Y, Cabrera J, Cheng X (2025). *WCluster: Clustering and PCA
with Weights, and Data Nuggets Clustering*. R package version 1.3.0.

## See also

[`nugget_feature_weights`](https://lucp12891.github.io/MosaiClusteR/reference/nugget_feature_weights.md),
[`M_ABC`](https://lucp12891.github.io/MosaiClusteR/reference/M_ABC.md)

## Examples

``` r
set.seed(1)
x <- rbind(matrix(rnorm(200 * 5,  0), ncol = 5),
           matrix(rnorm(200 * 5,  4), ncol = 5))
dn <- create_data_nuggets(x, max_nuggets = 40)
dn
#> <data_nugget>
#>   engine     : native 
#>   observations: 400 compressed into 40 nuggets
#>   reduction  : 90.0% 
#>   weight range: 3 - 18 (obs per nugget)
```
