# Simulate data for evaluating feature-weighting schemes (variance vs CV)

Generates a multi-source data set with a known cluster structure in
which the discriminative signal is carried either by
high-absolute-variance features (`regime = "variance"`) or by
high-*coefficient-of-variation* features buried under high-abundance,
high-variance decoys (`regime = "cv"`). It is intended for benchmarking
the `"var"`, `"cv"` and `"nugget"` weightings of the ABC family (and any
other feature-weighting method) on the same footing.

## Usage

``` r
simulate_weighting_regime(
  regime = c("cv", "variance"),
  n_samples = 120,
  n_clusters = 3,
  n_sources = 2,
  n_signal = 40,
  n_decoy = 500,
  n_noise = 60,
  effect = NULL,
  seed = NULL
)
```

## Arguments

- regime:

  `"cv"` (default) or `"variance"`; see Details.

- n_samples:

  Number of objects (rows), split evenly across clusters.

- n_clusters:

  Number of clusters.

- n_sources:

  Number of data sources (independent draws of the same design).

- n_signal, n_decoy, n_noise:

  Numbers of discriminative, decoy (non-discriminative but high-variance
  in the CV regime) and pure-noise features per source.

- effect:

  Between-cluster mean shift on the signal features (defaults: 1.4 for
  the variance regime, 0.9 for the cv regime).

- seed:

  Optional RNG seed.

## Value

A list with `List` (the sources, objects in rows), `truth` (cluster
labels), `signal_features` (indices of the discriminative features),
`regime` and `params`.

## See also

[`M_ABC`](https://lucp12891.github.io/MosaiClusteR/reference/M_ABC.md)
(the `weighting` argument),
[`mosaic_sim`](https://lucp12891.github.io/MosaiClusteR/reference/mosaic_sim.md)

## Examples

``` r
cv  <- simulate_weighting_regime("cv", seed = 1)
# variance weighting is misled by high-abundance decoys, CV recovers the truth:
fit_var <- M_ABC(cv$List, numsim = 80, NC = 3, weighting = c("var", "var"))
#> M_ABC: source 1 of 2 (weighting = var) ...
#> M_ABC: source 2 of 2 (weighting = var) ...
#> M_ABC: building consensus dissimilarity (C++ kernel) ...
fit_cv  <- M_ABC(cv$List, numsim = 80, NC = 3, weighting = c("cv",  "cv"))
#> M_ABC: source 1 of 2 (weighting = cv) ...
#> M_ABC: source 2 of 2 (weighting = cv) ...
#> M_ABC: building consensus dissimilarity (C++ kernel) ...
c(var = cluster_agreement(mosaic_labels(fit_var, 3), cv$truth)[["ARI"]],
  cv  = cluster_agreement(mosaic_labels(fit_cv,  3), cv$truth)[["ARI"]])
#>         var          cv 
#> 0.002749162 0.002861820 
```
