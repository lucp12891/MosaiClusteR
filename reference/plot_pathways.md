# Modern Gene Ontology / pathway-enrichment plots

Publication-quality replacements for the classic green GO DAG and the
horizontal bar plot. Given the output of
[`PathwayAnalysis()`](https://lucp12891.github.io/MosaiClusteR/reference/PathwayAnalysis.md)
(or a single MLP result table) it draws either a *dot plot* that
compares the enriched terms across clusters (the enrichment-map /
`clusterProfiler::dotplot` idiom) or a *lollipop* plot for a single
cluster. Dot area encodes the number of genes in the term and colour
encodes \\-\log\_{10}\\ p on a perceptually-uniform viridis scale.

## Usage

``` r
plot_pathways(
  pathways,
  style = c("dotplot", "lollipop", "circular"),
  top = 8,
  cluster = NULL,
  source = 1,
  title = NULL,
  palette = "D",
  max_chars = 42,
  wrap = 20,
  cluster_cols = NULL
)
```

## Arguments

- pathways:

  Output of
  [`PathwayAnalysis()`](https://lucp12891.github.io/MosaiClusteR/reference/PathwayAnalysis.md)
  (a list of sources, each a list of per-cluster results carrying a
  `$Pathways` table), a single source's per-cluster list, or a single
  MLP/`data.frame` result table.

- style:

  `"dotplot"` to compare clusters side by side, `"lollipop"` for one
  cluster, or `"circular"` for a radial (sunburst-style) bar plot that
  groups the terms by cluster around a central hub; bar length encodes
  \\-\log\_{10}\\ p, colour encodes the cluster, and each term is a
  single-line radial label that stays legible as the number of terms
  grows.

- top:

  Number of top terms (by p-value) to show per cluster.

- cluster:

  For `"lollipop"` (or to subset a dot plot): which cluster to show, by
  name (e.g. `"Cluster 1"`) or index. `NULL` shows all.

- source:

  Which source (index or name) of a full
  [`PathwayAnalysis()`](https://lucp12891.github.io/MosaiClusteR/reference/PathwayAnalysis.md)
  object to plot. Default `1`.

- title:

  Optional plot title.

- palette:

  viridis option passed to
  [`ggplot2::scale_colour_viridis_c()`](https://ggplot2.tidyverse.org/reference/scale_viridis.html)
  (`"D"` viridis, `"C"` plasma, `"B"` inferno, `"G"` mako, ...).

- max_chars:

  Term descriptions longer than this are wrapped onto two lines for
  legibility (dot plot / lollipop).

- wrap:

  For `style = "circular"`, the character width at which long term
  labels are wrapped onto successive lines, so long names stay compact
  and legible around the ring. Smaller values wrap more aggressively.

- cluster_cols:

  For `style = "circular"`, an optional vector of fill colours for the
  clusters (cluster bars and legend). Supply the same palette used
  elsewhere - e.g. the `cols` passed to
  [`Cyclogram()`](https://lucp12891.github.io/MosaiClusteR/reference/Cyclogram.md) -
  so the clusters share one colour scheme across figures. Positional by
  cluster order, or a vector named by cluster. Default `NULL` uses the
  package palette.

## Value

A ggplot object (invisibly); it is also printed so the function works
directly on an open graphics device.

## See also

[`PathwayAnalysis()`](https://lucp12891.github.io/MosaiClusteR/reference/PathwayAnalysis.md),
[`PlotPathways()`](https://lucp12891.github.io/MosaiClusteR/reference/PlotPathways.md)

## Examples

``` r
if (FALSE) { # \dontrun{
P <- PathwayAnalysis(List = L, geneExpr = ge, nrclusters = 3,
                     method = c("limma", "MLP"), geneInfo = geneInfo,
                     geneSetSource = "GOBP", sign = 0.05, niter = 2)
plot_pathways(P, style = "dotplot", top = 8)          # compare clusters
plot_pathways(P, style = "lollipop", cluster = 1)     # one cluster
} # }
```
