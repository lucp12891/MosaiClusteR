# Biological interpretation: genes and pathways

For omics data, the question after clustering is biological: which genes
separate the clusters, and which pathways do those genes implicate.
MosaiClusteR connects clusters to genes and gene sets through a set of
functions that build on Bioconductor.

## Dependencies

This layer uses Bioconductor packages (`limma`, `MLP`, `biomaRt`,
`org.Hs.eg.db`), declared under *Suggests*. Install the ones a function
names:

``` r

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("limma", "MLP", "biomaRt", "org.Hs.eg.db"))
```

Each function raises a clear error naming the exact install command if a
package is missing, so you only install what you use.

## Differential features per cluster

[`DiffGenes()`](https://lucp12891.github.io/MosaiClusteR/reference/DiffGenes.md)
finds the genes that distinguish the clusters, using `limma` by default.
[`DiffGenesSelection()`](https://lucp12891.github.io/MosaiClusteR/reference/DiffGenesSelection.md)
restricts the test to a chosen set of objects, and
[`FindGenes()`](https://lucp12891.github.io/MosaiClusteR/reference/FindGenes.md)
looks up genes of interest.

``` r

library(MosaiClusteR)
# geneExpr is a genes-by-samples expression matrix aligned to the objects in L
dg <- DiffGenes(List = L, geneExpr = geneExpr, nrclusters = 3, method = "limma")
```

## Pathway enrichment

From the differential genes, test which gene sets are enriched per
cluster:

| Function | Purpose |
|----|----|
| `PathwayAnalysis`, `Pathways` | enrichment per cluster |
| `PathwaysIter`, `PathwaysSelection` | iterate over clusters, or a chosen selection |
| `PreparePathway` | assemble inputs for the enrichment step |
| `plot_pathways`, `PlotPathways` | visualise enriched terms |
| `enriched_features` | turn characterisation output into per-cluster term tables |

``` r

pa <- PathwayAnalysis(List = L, geneExpr = geneExpr, nrclusters = 3,
                      method = c("limma", "MLP"), geneSetSource = "GOBP")
plot_pathways(pa, style = "circular")     # also "dotplot" and "lollipop"
```

[`plot_pathways()`](https://lucp12891.github.io/MosaiClusteR/reference/plot_pathways.md)
offers three layouts. The `"dotplot"` and `"lollipop"` styles rank terms
in a column; the `"circular"` style groups terms by cluster around a
central hub, with bar length encoding significance and colour encoding
the cluster. Pass `cluster_cols` to match the colours used elsewhere in
your figures.

## Shared genes and pathways across solutions

When you have several clusterings, ask what they agree on biologically:

| Function | Purpose |
|----|----|
| `Geneset.intersect`, `Geneset.intersectSelection` | gene sets shared across solutions |
| `SharedGenesPathsFeat` | genes, pathways and features shared across solutions |
| `SharedSelection`, `SharedSelectionLimma`, `SharedSelectionMLP` | shared results for a chosen selection |

``` r

shared <- Geneset.intersect(PathwaysOutput = pa)
```

## A short interpretation workflow

``` r

lab <- mosaic_labels(
  SNF(L, type = "data", distmeasure = rep("euclidean", length(L)),
      normalize = rep(FALSE, length(L)), method = rep(list(NULL), length(L))), 3)

dg <- DiffGenes(List = L, geneExpr = geneExpr, nrclusters = 3, method = "limma")
pa <- PathwayAnalysis(List = L, geneExpr = geneExpr, nrclusters = 3,
                      method = c("limma", "MLP"), geneSetSource = "GOBP")
plot_pathways(pa, style = "dotplot")
```

Cluster first, confirm the clusters are stable and separated, then read
their biology from the differential genes and the enriched pathways.
