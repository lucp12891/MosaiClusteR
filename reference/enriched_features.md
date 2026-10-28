# Tidy CharacteristicFeatures() output for plotting

Turns the output of
[`CharacteristicFeatures()`](https://lucp12891.github.io/MosaiClusteR/reference/CharacteristicFeatures.md)
into a per-cluster list of compact tables (`term` / `pvalue` / `Count`)
ready for
[`plot_pathways()`](https://lucp12891.github.io/MosaiClusteR/reference/plot_pathways.md).
For each cluster it keeps the characteristic features that pass a
significance cut and, by default, only those that are *enriched* in the
cluster (more frequent among its objects than in the remaining objects)
rather than depleted. Empty clusters are dropped.

## Usage

``` r
enriched_features(
  char,
  data,
  source = 1,
  method = 1,
  sign = 0.05,
  enriched_only = TRUE,
  clean_names = TRUE
)
```

## Arguments

- char:

  Output of
  [`CharacteristicFeatures()`](https://lucp12891.github.io/MosaiClusteR/reference/CharacteristicFeatures.md).

- data:

  The binary data matrix whose features were tested (objects in rows,
  features in columns) - e.g. the target matrix. Used to set the
  enrichment direction and the per-cluster feature counts.

- source:

  Which tested data source to use within each cluster, by name (e.g.
  `"Target"`) or index. Default: the first source.

- method:

  Which clustering result in `char` to use, by name or index. Default
  `1`.

- sign:

  Keep features with adjusted p-value below this. Default `0.05`.

- enriched_only:

  If `TRUE` (default) keep only features more frequent in the cluster
  than in the rest of the objects; if `FALSE` keep every significant
  characteristic feature (enriched or depleted).

- clean_names:

  If `TRUE` (default) tidy the feature names for display - reversing the
  [`make.names()`](https://rdrr.io/r/base/make.names.html) mangling
  common in binary descriptor matrices (drop a leading `X` before a
  digit, turn `.`/`_` runs into spaces), e.g. `"X3.oxo.5.alpha_steroid"`
  becomes `"3 oxo 5 alpha steroid"`. Set `FALSE` to keep the raw column
  names.

## Value

A named list - one element per non-empty cluster - of data frames with
columns `term`, `pvalue` (adjusted p-value) and `Count` (number of the
cluster's objects carrying the feature), suitable for
[`plot_pathways()`](https://lucp12891.github.io/MosaiClusteR/reference/plot_pathways.md).

## See also

[`CharacteristicFeatures()`](https://lucp12891.github.io/MosaiClusteR/reference/CharacteristicFeatures.md),
[`plot_pathways()`](https://lucp12891.github.io/MosaiClusteR/reference/plot_pathways.md)

## Examples

``` r
if (FALSE) { # \dontrun{
Char <- CharacteristicFeatures(List = list(Weighted = sol),
                               binData = list(Target = targetMat),
                               datanames = "Target", nrclusters = 4)
tabs <- enriched_features(Char, targetMat)
plot_pathways(tabs, style = "circular")
} # }
```
