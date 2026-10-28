# Simulate a multi-source data set with planted clusters

Generates `K` feature-by-sample matrices that share a common set of `n`
samples drawn from `g` ground-truth groups. Each source has its own
informative features (shifted between groups) plus noise features, so
the clusters are only partially visible in any single source - the
canonical setting in which multi-source integration helps.

## Usage

``` r
mosaic_sim(
  n = 60,
  g = 3,
  p = c(120, 150),
  informative = c(20, 25),
  effect = 2.5,
  K = 2,
  seed = NULL
)
```

## Arguments

- n:

  Number of samples (objects).

- g:

  Number of ground-truth groups.

- p:

  Per-source feature counts (recycled to length `K`).

- informative:

  Per-source number of group-informative features.

- effect:

  Between-group mean shift for informative features.

- K:

  Number of sources.

- seed:

  Optional RNG seed.

## Value

A list with `List` (the `K` matrices, **objects in rows**, features in
columns) and `truth` (the integer group label per object).

## Examples

``` r
sim <- mosaic_sim(n = 30, g = 3, K = 2, seed = 1)
lengths(lapply(sim$List, dim))
#> source1 source2 
#>       2       2 
```
