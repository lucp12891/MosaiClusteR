# LUCID: latent unknown clustering integrating omics with an outcome

Thin wrapper around LUCIDus' `estimate_lucid`. LUCID is a *model-based,
supervised* integrative clustering method: it defines a categorical
latent variable (the clusters) jointly from upstream exposures `G`,
multi-omics data `Z` and an outcome `Y` via an EM algorithm (a
quasi-mediation model). Unlike the unsupervised methods in MosaiClusteR
it requires an outcome and returns posterior cluster membership.

## Usage

``` r
LUCID(
  G,
  Z,
  Y,
  K,
  lucid_model = c("early", "parallel", "serial"),
  family = c("normal", "binary"),
  ...
)
```

## Arguments

- G:

  Exposure/predictor matrix (objects in rows).

- Z:

  Omics matrix, or a list of matrices for the parallel/serial models
  (objects in rows).

- Y:

  Outcome vector (length = number of objects).

- K:

  Number of latent clusters.

- lucid_model:

  One of `"early"`, `"parallel"`, `"serial"`.

- family:

  Outcome family, `"normal"` or `"binary"`.

- ...:

  Further arguments forwarded to
  [`LUCIDus::estimate_lucid`](https://rdrr.io/pkg/LUCIDus/man/estimate_lucid.html).

## Value

A list of class `"LUCID"` with the fitted `model` and the posterior
`cluster` assignment per object.

## References

Zhao, Y., Jia, K., Goodrich, J. & Conti, D. V. (2024). LUCIDus: an R
package for implementing LUCID with phenotypic traits. The R Journal
16(2).

## See also

[`M_ABC`](https://lucp12891.github.io/MosaiClusteR/reference/M_ABC.md),
[`intNMF`](https://lucp12891.github.io/MosaiClusteR/reference/intNMF.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# requires the 'LUCIDus' package
data(mosaic_toy)
Z <- mosaic_toy$List[[1]]
G <- mosaic_toy$List[[2]][, 1:3]
Y <- as.numeric(mosaic_toy$truth)
fit <- LUCID(G = G, Z = Z, Y = Y, K = 3)
table(fit$cluster)
} # }
```
