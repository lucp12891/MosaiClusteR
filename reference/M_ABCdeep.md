# Multi-source M-ABC with a deep-learning (autoencoder) Step 4 (M-ABCdeep)

The deep-learning counterpart of
[`M_ABC`](https://lucp12891.github.io/MosaiClusteR/reference/M_ABC.md):
runs
[`ABCdeep.SingleInMultiple`](https://lucp12891.github.io/MosaiClusteR/reference/ABCdeep.SingleInMultiple.md)
per source (autoencoder embedding + latent clustering for the base
partitions), row-stacks the per-iteration label matrices, and builds one
consensus dissimilarity via
[`f.clustABC.MultiSource`](https://lucp12891.github.io/MosaiClusteR/reference/f.clustABC.MultiSource.md).
The ABC algorithm of (Amaratunga et al. 2008) is preserved end to end;
only Step 4 is deep. There is no Python/torch dependency and no PCA
fallback (see
[`ABCdeep.SingleInMultiple`](https://lucp12891.github.io/MosaiClusteR/reference/ABCdeep.SingleInMultiple.md)).

## Usage

``` r
M_ABCdeep(
  List,
  weighting = rep("var", length(List)),
  normalize = rep(FALSE, length(List)),
  gr = c(),
  bag = rep(TRUE, length(List)),
  numsim = 500,
  numvar = rep(100, length(List)),
  NC = NULL,
  dissimilarity = "mabc",
  final_linkage = "ward",
  alpha = 0.625,
  mds = FALSE,
  topk_frac = 0.2,
  latent_dim = "auto",
  hidden_units = 32L,
  ae_epochs = 80L,
  ae_lr = 0.001,
  cluster_in_latent = "ward",
  linkage = rep("ward.D2", length(List)),
  nugget_type = "between",
  nugget_args = list(),
  min_samples = 40L,
  seed = NULL
)
```

## Arguments

- List:

  A list of \\K\\ data matrices over the same objects, **objects
  (samples) in rows** and features in columns.

- weighting, normalize, bag, numvar:

  Per-source vectors (recycled to \\K\\) forwarded to
  [`ABCdeep.SingleInMultiple`](https://lucp12891.github.io/MosaiClusteR/reference/ABCdeep.SingleInMultiple.md).

- gr:

  Unused; kept for compatibility.

- numsim:

  Iterations per source.

- NC:

  Base-cluster count per source (`NULL` = data-driven).

- dissimilarity:

  Consensus dissimilarity, `"mabc"` or `"abc"`.

- final_linkage, alpha:

  Final clustering controls.

- mds:

  Logical; MDS of the consensus dissimilarities.

- topk_frac:

  Fraction of top-weighted features the autoencoder focuses on (default
  0.2); forwarded to
  [`ABCdeep.SingleInMultiple`](https://lucp12891.github.io/MosaiClusteR/reference/ABCdeep.SingleInMultiple.md).

- latent_dim, hidden_units, ae_epochs, ae_lr, cluster_in_latent,
  linkage:

  Autoencoder / latent-clustering controls (scalar or length-\\K\\),
  forwarded to
  [`ABCdeep.SingleInMultiple`](https://lucp12891.github.io/MosaiClusteR/reference/ABCdeep.SingleInMultiple.md).

- nugget_type, nugget_args:

  Forwarded to the nugget weighting.

- min_samples:

  Warn (per source) when \\N\\ is below this; the autoencoder is
  underdetermined at very small \\N\\ (default 40).

- seed:

  Optional RNG seed.

## Value

A list of class `"M_ABC"` with `DistM` and `Clust`.

## References

Amaratunga D, Cabrera J, Kovtun V (2008). “Microarray Learning with
ABC.” *Biostatistics*, **9**(1), 128–136.
[doi:10.1093/biostatistics/kxm017](https://doi.org/10.1093/biostatistics/kxm017)
.
<https://academic.oup.com/biostatistics/article-abstract/9/1/128/254198>.

## See also

[`ABCdeep.SingleInMultiple`](https://lucp12891.github.io/MosaiClusteR/reference/ABCdeep.SingleInMultiple.md),
[`M_ABC`](https://lucp12891.github.io/MosaiClusteR/reference/M_ABC.md),
[`f.clustABC.MultiSource`](https://lucp12891.github.io/MosaiClusteR/reference/f.clustABC.MultiSource.md)

## Examples

``` r
data(mosaic_toy)
fit <- M_ABCdeep(mosaic_toy$List, numsim = 40, NC = 3, ae_epochs = 40, seed = 1)
#> M_ABCdeep: source 1 of 2 (weighting = var, latent = auto) ...
#> ABCdeep: training autoencoder (G=120, top-K=24, N=60, hidden=32, latent=8, epochs=40) ...
#> M_ABCdeep: source 2 of 2 (weighting = var, latent = auto) ...
#> ABCdeep: training autoencoder (G=150, top-K=30, N=60, hidden=32, latent=8, epochs=40) ...
#> M_ABCdeep: building consensus dissimilarity (C++ kernel) ...
table(stats::cutree(stats::as.hclust(fit$Clust), 3))
#> 
#>  1  2  3 
#> 19 25 16 
```
