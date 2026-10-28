# Single-source ABC with a deep-learning (autoencoder) Step 4 (ABCdeep)

A deep-learning variant of
[`ABCpp.SingleInMultiple`](https://lucp12891.github.io/MosaiClusteR/reference/ABCpp.SingleInMultiple.md)
that keeps the ABC algorithm of (Amaratunga et al. 2008) intact and
replaces **only Step 4**: instead of clustering the raw feature
sub-matrix, a neural autoencoder embeds the objects into a latent space
and the latent vectors are clustered into `NC` base clusters. The
stacked label matrix is fused across sources by
[`f.clustABC.MultiSource`](https://lucp12891.github.io/MosaiClusteR/reference/f.clustABC.MultiSource.md),
exactly as in
[`M_ABC`](https://lucp12891.github.io/MosaiClusteR/reference/M_ABC.md).

## Usage

``` r
ABCdeep.SingleInMultiple(
  data,
  weighting = "var",
  normalize = NULL,
  gr = c(),
  bag = TRUE,
  numsim = 1000,
  numvar = 100,
  NC = NULL,
  topk_frac = 0.2,
  latent_dim = "auto",
  hidden_units = 32L,
  ae_epochs = 80L,
  ae_lr = 0.001,
  cluster_in_latent = "ward",
  linkage = "ward.D2",
  nugget_type = "between",
  nugget_args = list(),
  min_samples = 40L,
  seed = NULL
)
```

## Arguments

- data:

  Numeric matrix, **objects (samples) in rows**, features in columns
  (the MosaiClusteR convention).

- weighting:

  Feature weighting used to rank features for the top-K focus: `"var"`,
  `"cv"`, `"nugget"` or equal (as in
  [`ABCpp.SingleInMultiple`](https://lucp12891.github.io/MosaiClusteR/reference/ABCpp.SingleInMultiple.md)).

- normalize:

  Optional normalisation applied before the internal min-max scaling
  (see
  [`Normalization`](https://lucp12891.github.io/MosaiClusteR/reference/Normalization.md));
  `NULL`/`FALSE` = none.

- gr:

  Unused; kept for interface compatibility.

- bag:

  Logical; bootstrap the objects each iteration (Step 2).

- numsim:

  Number of iterations \\R\\.

- numvar:

  Kept for backward compatibility; unused.

- NC:

  Base-cluster count (`NULL` = \\\lfloor\sqrt{N}\rfloor\\).

- topk_frac:

  Fraction of the highest-weighted features the autoencoder is
  trained/embedded on (default 0.2). Focusing on the informative
  features drops the noise that would dilute a full-feature embedding on
  sparse signals, while keeping enough features for dense signals; the
  count is floored at `min(G, 20)` and capped at `min(G, 4000)`.

- latent_dim:

  Autoencoder bottleneck size: an integer, or `"auto"` (default) =
  `max(NC + 2, 8)`, capped at the hidden width and \\K-1\\.

- hidden_units:

  Hidden-layer width of the encoder/decoder.

- ae_epochs:

  Training epochs for the one-time autoencoder fit.

- ae_lr:

  Adam learning rate.

- cluster_in_latent:

  Latent clustering (Step 4b): `"ward"` (default, hierarchical) or
  `"kmeans"`.

- linkage:

  Agglomeration method when `cluster_in_latent = "ward"`.

- nugget_type, nugget_args:

  Passed to the nugget weighting when `weighting = "nugget"`.

- min_samples:

  Emit a warning when `nrow(data)` is below this (the autoencoder is
  underdetermined at very small \\N\\); default 40.

- seed:

  Optional RNG seed (controls the autoencoder init and the bootstrap).

## Value

A data frame, `numsim` rows by one column per object; entries are
base-cluster labels (`0` = object not selected that iteration) – the
same contract as
[`ABCpp.SingleInMultiple`](https://lucp12891.github.io/MosaiClusteR/reference/ABCpp.SingleInMultiple.md).

## Speed (trained once, forward pass per iteration)

The autoencoder is trained a single time on the source's top-`topk_frac`
most-informative features; each of the `numsim` iterations then performs
only a forward pass over the bootstrapped objects. The network is a
self-contained, vectorised R autoencoder (He initialisation, ReLU, Adam,
MSE) – there is no Python, TensorFlow or torch dependency, and **no PCA
or other fallback**: a single code path.

## References

Amaratunga D, Cabrera J, Kovtun V (2008). “Microarray Learning with
ABC.” *Biostatistics*, **9**(1), 128–136.
[doi:10.1093/biostatistics/kxm017](https://doi.org/10.1093/biostatistics/kxm017)
.
<https://academic.oup.com/biostatistics/article-abstract/9/1/128/254198>.

## See also

[`M_ABCdeep`](https://lucp12891.github.io/MosaiClusteR/reference/M_ABCdeep.md),
[`ABCpp.SingleInMultiple`](https://lucp12891.github.io/MosaiClusteR/reference/ABCpp.SingleInMultiple.md),
[`f.clustABC.MultiSource`](https://lucp12891.github.io/MosaiClusteR/reference/f.clustABC.MultiSource.md)

## Examples

``` r
data(mosaic_toy)
lab <- ABCdeep.SingleInMultiple(mosaic_toy$List[[1]], numsim = 40, NC = 3,
                                ae_epochs = 40, seed = 1)
#> ABCdeep: training autoencoder (G=120, top-K=24, N=60, hidden=32, latent=8, epochs=40) ...
dim(lab)
#> [1] 40 60
```
