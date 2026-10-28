# Visualisation and cluster characterisation

A clustering is only useful once you can see it and say what each
cluster is. This page covers the plotting suite and the characterisation
tools. The heavier plots declare their dependencies under *Suggests*
(`circlize`, `plotrix`, `gplots`, `ggplot2`, `uwot`, `Rtsne`) and raise
an informative error if one is missing.

``` r

library(MosaiClusteR)
data(mosaic_toy)
L   <- mosaic_toy$List
fit <- SNF(L, type = "data", distmeasure = rep("euclidean", length(L)),
           normalize = rep(FALSE, length(L)), method = rep(list(NULL), length(L)))
lab <- mosaic_labels(fit, 3)
```

## The fitted hierarchy

Every fitted result carries a hierarchy in `$Clust` that you can draw
with base R, no extra package needed:

``` r

hc <- stats::as.hclust(fit$Clust)
plot(hc, labels = FALSE, main = "SNF consensus", xlab = "", sub = "")
stats::rect.hclust(hc, k = 3, border = 2:4)
```

## Comparing and drawing solutions

| Function | Shows | Needs |
|----|----|----|
| `ComparePlot` | several solutions aligned side by side | `circlize`, `plotrix` |
| `CompareInteractive` | the same, with interactive selection | `circlize`, `plotrix` |
| `ClusterPlot`, `Cyclogram` | one solution as a coloured tree or cyclogram | `circlize` |
| `SimilarityHeatmap`, `HeatmapPlot`, `HeatmapSelection` | object similarity and feature heatmaps | `gplots` |
| `BoxPlotDistance` | within- versus between-cluster distances | base |
| `embedding_plot` | objects in a 2-D PCA, UMAP or t-SNE embedding | `uwot`, `Rtsne` |

``` r

fit_a <- ADC(L, distmeasure = "euclidean", normalize = FALSE)
ComparePlot(list(ADC = fit_a, SNF = fit), nrclusters = 3)
SimilarityHeatmap(fit$DistM)
embedding_plot(L[[1]], lab, method = "pca")     # or "umap", "tsne"
```

Consistent colours across plots come from the colour helpers:
[`ColorPalette()`](https://lucp12891.github.io/MosaiClusteR/reference/ColorPalette.md),
[`ColorsNames()`](https://lucp12891.github.io/MosaiClusteR/reference/ColorsNames.md),
[`ClusterCols()`](https://lucp12891.github.io/MosaiClusteR/reference/ClusterCols.md)
and
[`LabelCols()`](https://lucp12891.github.io/MosaiClusteR/reference/LabelCols.md).
Reuse the same palette so a cluster keeps its colour from one figure to
the next.

## Per-cluster feature profiles

To see what varies across clusters, plot the features:

| Function | Shows |
|----|----|
| `ProfilePlot`, `plot_cluster_profile` | mean feature profile per cluster |
| `ContFeaturesPlot` | continuous features across clusters |
| `BinFeaturesPlot_SingleData`, `BinFeaturesPlot_MultipleData` | binary features (one or several sources) |

``` r

ProfilePlot(L[[1]], lab)
ContFeaturesPlot(leadCpds = rownames(L[[1]]), data = L[[1]], nrclusters = 3)
```

## Characterising clusters

Describe why clusters differ, from external variables and from the
features themselves.

- **External variables.**
  [`characterize_clusters()`](https://lucp12891.github.io/MosaiClusteR/reference/characterize_clusters.md)
  relates a label vector to a metadata frame of continuous and
  categorical variables, returning per-variable tests and summaries.
  [`cluster_biology_panel()`](https://lucp12891.github.io/MosaiClusteR/reference/cluster_biology_panel.md)
  draws the panel.
- **Driving features.**
  [`CharacteristicFeatures()`](https://lucp12891.github.io/MosaiClusteR/reference/CharacteristicFeatures.md)
  finds the features that define each cluster;
  [`FeatSelection()`](https://lucp12891.github.io/MosaiClusteR/reference/FeatSelection.md)
  and
  [`FeaturesOfCluster()`](https://lucp12891.github.io/MosaiClusteR/reference/FeaturesOfCluster.md)
  focus on one source or one cluster.

``` r

meta <- data.frame(score = rnorm(length(lab)),
                   grp   = sample(c("A", "B"), length(lab), TRUE))
characterize_clusters(lab, meta, continuous = "score", categorical = "grp")
CharacteristicFeatures(L, lab)
FeaturesOfCluster(L[[1]], lab, cluster = 1)
```

## Navigating and tracking

| Function | Purpose |
|----|----|
| `FindCluster`, `FindElement` | locate a cluster or an object in a result |
| `TrackCluster` | follow a cluster across several solutions |
| `SharedComps`, `SimilarityMeasure` | shared members and similarity between solutions |
| `ReorderToReference` | align one solution’s labels to a reference |

``` r

TrackCluster(list(ADC = fit_a, SNF = fit), nrclusters = 3)
SharedComps(mosaic_labels(fit_a, 3), lab)
```
