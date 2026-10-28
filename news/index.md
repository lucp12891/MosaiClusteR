# Changelog

## MosaiClusteR 1.0.0

First stable release. This version focuses on correctness: several
functions that silently ignored an argument, mishandled an edge case, or
failed to align multi-source inputs are fixed, and the documentation and
vignette examples are corrected to match. Changes since 0.1.4:

### Correctness fixes

- [`Cluster()`](https://lucp12891.github.io/MosaiClusteR/reference/Cluster.md)
  now forwards the normalisation method to
  [`Distance()`](https://lucp12891.github.io/MosaiClusteR/reference/Distance.md).
  Previously the `method` argument was documented and accepted but
  silently ignored, and `normalize = TRUE` did not select a method.
  `normalize` now accepts either a logical (`FALSE` for none) or a
  normalisation method name, and `method` sets the method when
  `normalize = TRUE`.
- Documentation examples now use the objects-in-rows orientation
  consistently (the
  [`Cluster()`](https://lucp12891.github.io/MosaiClusteR/reference/Cluster.md)
  example no longer transposes the bundled `mosaic_toy` data).
- [`cluster_agreement()`](https://lucp12891.github.io/MosaiClusteR/reference/cluster_agreement.md)
  edge cases. Identical partitions (up to relabelling) now return 1 for
  ARI, NMI and Jaccard, including the degenerate single-cluster and
  all-singleton cases that previously returned 0. Missing (`NA`) labels
  now raise a clear error instead of a cryptic “missing value where
  TRUE/FALSE needed”.
- [`ADC()`](https://lucp12891.github.io/MosaiClusteR/reference/ADC.md)
  and
  [`ADEC()`](https://lucp12891.github.io/MosaiClusteR/reference/ADEC.md)
  validate and align the sources by object identifier before combining
  them. Reordered rows are aligned to the first source; sources without
  row names are combined positionally (previously this failed with “need
  at least two non-NA values to interpolate”); and mismatched object
  sets, differing row counts or duplicated identifiers now raise clear
  errors instead of a cryptic “subscript out of bounds”.
- [`mosaic_labels()`](https://lucp12891.github.io/MosaiClusteR/reference/mosaic_labels.md)
  validates `k` (a single positive integer) with a clear message,
  instead of silently truncating a non-integer `k` or failing
  cryptically on `NA`.
- [`M_ABC()`](https://lucp12891.github.io/MosaiClusteR/reference/M_ABC.md)
  and
  [`ADEC()`](https://lucp12891.github.io/MosaiClusteR/reference/ADEC.md)
  gained a `seed` argument for reproducible results, consistent with
  [`M_ABCdeep()`](https://lucp12891.github.io/MosaiClusteR/reference/M_ABCdeep.md)
  and
  [`intNMF()`](https://lucp12891.github.io/MosaiClusteR/reference/intNMF.md).
  Their documentation states that they are stochastic and that `seed`
  (or [`set.seed()`](https://rdrr.io/r/base/Random.html)) fixes the
  result.
- [`M_ABC()`](https://lucp12891.github.io/MosaiClusteR/reference/M_ABC.md)
  now aligns its sources by object identifier before building the
  consensus, so a source with reordered rows no longer changes the
  result; it also evaluates its `List` argument before seeding so an
  argument that uses the RNG cannot offset the internal random stream.
  [`CEC()`](https://lucp12891.github.io/MosaiClusteR/reference/CEC.md)
  uses the same alignment and now accepts sources without row names
  (combined positionally) instead of failing with “need at least two
  non-NA values to interpolate”.
- Error messages for missing optional Bioconductor packages (used by the
  pathway and gene-set functions) now give the exact installation
  command, for example `BiocManager::install("MLP")`, instead of only
  “needed for this function to work”. The CRAN-package messages already
  named
  [`install.packages()`](https://rdrr.io/r/utils/install.packages.html).

### Documentation and examples

- The vignette preprocessing example now computes distances from the
  normalised data, rather than normalising and then discarding the
  result.
- The vignette method-comparison example reports each method’s name,
  success or failure status, runtime and a useful error message, instead
  of returning a bare `NA` on failure.
- [`SelectnrClusters()`](https://lucp12891.github.io/MosaiClusteR/reference/SelectnrClusters.md)
  returns one recommendation per criterion; the vignette now reduces
  these to a single validated integer before using it as a cluster
  count, rather than passing the whole table.
- The five integration paradigms are stated consistently.
  Factor/low-rank and spectral methods (`intNMF`, `spectral_clustering`,
  `LUCID`) are presented as additional methods rather than a sixth
  paradigm.

## MosaiClusteR 0.1.4

CRAN release: 2026-10-21

Third public release. It consolidates the development carried out after
0.1.1, which had not previously been published.

### Breaking change

- **`M_ABCpp()` renamed to
  [`M_ABC()`](https://lucp12891.github.io/MosaiClusteR/reference/M_ABC.md).**
  The multi-source Aggregating Bundles of Clusters function is now
  [`M_ABC()`](https://lucp12891.github.io/MosaiClusteR/reference/M_ABC.md);
  its arguments, behaviour and return value are unchanged. Update any
  calls from `M_ABCpp(...)` to `M_ABC(...)`.

### New features

- The example dataset `mosaic_toy` is bundled in `data/` and used by the
  examples and tests, so `R CMD check` (and the r-universe builds) no
  longer fail with “object ‘mosaic_toy’ not found”.
- **`plot_pathways(style = "circular")`.** A radial (sunburst-style) bar
  plot that groups the enriched terms by cluster around a central hub.
  Bar length encodes p and colour encodes the cluster; each term is a
  single-line radial label that stays legible as the number of terms
  grows. Complements the existing `"dotplot"` and `"lollipop"` styles.
  It also gained `cluster_cols` (an explicit per-cluster palette,
  e.g. to match
  [`Cyclogram()`](https://lucp12891.github.io/MosaiClusteR/reference/Cyclogram.md))
  and `wrap` (radial label wrapping).
- **[`enriched_features()`](https://lucp12891.github.io/MosaiClusteR/reference/enriched_features.md).**
  A helper that turns
  [`CharacteristicFeatures()`](https://lucp12891.github.io/MosaiClusteR/reference/CharacteristicFeatures.md)
  output into a per-cluster list of compact `term`/`pvalue`/`Count`
  tables, ready to hand straight to
  [`plot_pathways()`](https://lucp12891.github.io/MosaiClusteR/reference/plot_pathways.md).
- [`FeaturesOfCluster()`](https://lucp12891.github.io/MosaiClusteR/reference/FeaturesOfCluster.md)
  and
  [`BinFeaturesPlot_SingleData()`](https://lucp12891.github.io/MosaiClusteR/reference/BinFeaturesPlot_SingleData.md)
  gained `cexB`, `cexR` and `cexN` arguments for the object, feature and
  name label sizes;
  [`Cyclogram()`](https://lucp12891.github.io/MosaiClusteR/reference/Cyclogram.md)
  now defaults to a larger `cex.labels`.

### Bug fixes

- **[`CharacteristicFeatures()`](https://lucp12891.github.io/MosaiClusteR/reference/CharacteristicFeatures.md).**
  When called with `binData` (the common case, as in the documented
  example) an inner loop index overwrote the cluster index, so only the
  first cluster’s characteristic features were returned, and a
  `datanames` mis-index left the per-source result list unnamed. Both
  are fixed: features are now computed and named for every cluster.
- **`Distance(distmeasure = "euclidean")`** now computes a true
  Euclidean distance (via
  [`stats::dist`](https://rdrr.io/r/stats/dist.html)). Previously it
  routed to [`FD::gowdis`](https://rdrr.io/pkg/FD/man/gowdis.html)
  (Gower), which returned `NaN` when a feature was constant across the
  selected objects, producing NA distance matrices. Gower is still
  available via `distmeasure = "gower"`.
- **[`ComparePlot()`](https://lucp12891.github.io/MosaiClusteR/reference/ComparePlot.md)**
  (and
  [`Cyclogram()`](https://lucp12891.github.io/MosaiClusteR/reference/Cyclogram.md))
  now build a default colour palette when `cols` is not supplied, and
  check up front that all results describe the same objects, raising a
  clear message otherwise.

### Metadata

- The licence is **GPL-3**, consistent with 0.1.0 and 0.1.1. An interim
  local change to CC BY-SA 4.0 is reverted, since CC licences are
  intended for content rather than software.
- “MoSaIC” expands to *Multi-Omics Similarity Aggregation and
  Integrative Clustering in R*, and the README framework overview is a
  scalable SVG figure (`man/figures/framework.svg`) that reads on both
  light and dark backgrounds.

## MosaiClusteR 0.1.1

CRAN release: 2026-08-21

- Documentation and installation updates to the README, and a package
  version bump from 0.1.0.

## MosaiClusteR 0.1.0

CRAN release: 2026-07-29

First public release. MosaiClusteR is an umbrella framework for
multi-source and multi-omics clustering: it gathers a broad catalogue of
integration methods behind one consistent interface, together with the
tooling needed to preprocess, select, evaluate, compare, characterise
and biologically interpret the results.

### Framework

- **Unified interface.** Every method takes the same input — a `List` of
  numeric matrices over the same objects (objects in rows) — and returns
  the same `{DistM, Clust}` contract, so methods are interchangeable and
  directly comparable. Precomputed dissimilarities are accepted via
  `type = "dist"`.
- **Four-layer workflow.** Preprocess → baseline → integrate → evaluate.

### Integration methods (five paradigms)

- **Direct:**
  [`ADC()`](https://lucp12891.github.io/MosaiClusteR/reference/ADC.md),
  [`ADEC()`](https://lucp12891.github.io/MosaiClusteR/reference/ADEC.md).
- **Similarity-based:**
  [`WeightedClust()`](https://lucp12891.github.io/MosaiClusteR/reference/WeightedClust.md),
  [`WonM()`](https://lucp12891.github.io/MosaiClusteR/reference/WonM.md),
  [`SNF()`](https://lucp12891.github.io/MosaiClusteR/reference/SNF.md),
  [`NEMO()`](https://lucp12891.github.io/MosaiClusteR/reference/NEMO.md).
- **Graph-based:**
  [`EnsembleClustering()`](https://lucp12891.github.io/MosaiClusteR/reference/EnsembleClustering.md)
  (CSPA/HGPA/MCLA),
  [`HBGF()`](https://lucp12891.github.io/MosaiClusteR/reference/HBGF.md),
  [`ClusteringAggregation()`](https://lucp12891.github.io/MosaiClusteR/reference/ClusteringAggregation.md).
- **Voting / consensus:**
  [`CEC()`](https://lucp12891.github.io/MosaiClusteR/reference/CEC.md),
  [`CVAA()`](https://lucp12891.github.io/MosaiClusteR/reference/CVAA.md),
  [`ConsensusClustering()`](https://lucp12891.github.io/MosaiClusteR/reference/ConsensusClustering.md),
  [`EvidenceAccumulation()`](https://lucp12891.github.io/MosaiClusteR/reference/EvidenceAccumulation.md),
  [`LinkBasedClustering()`](https://lucp12891.github.io/MosaiClusteR/reference/LinkBasedClustering.md)
  (CTS/SRS/ASRS), and the Aggregating Bundles of Clusters family
  [`M_ABC()`](https://lucp12891.github.io/MosaiClusteR/reference/M_ABC.md)
  /
  [`M_ABCdist()`](https://lucp12891.github.io/MosaiClusteR/reference/M_ABCdist.md)
  /
  [`M_ABCdist.WC()`](https://lucp12891.github.io/MosaiClusteR/reference/M_ABCdist.WC.md)
  /
  [`M_ABCdeep()`](https://lucp12891.github.io/MosaiClusteR/reference/M_ABCdeep.md)
  (a deep-learning, autoencoder-based variant).
- **Hierarchy-based:**
  [`HierarchicalEnsembleClustering()`](https://lucp12891.github.io/MosaiClusteR/reference/HierarchicalEnsembleClustering.md),
  [`EHC()`](https://lucp12891.github.io/MosaiClusteR/reference/EHC.md).
- **Factor / low-rank / spectral:**
  [`intNMF()`](https://lucp12891.github.io/MosaiClusteR/reference/intNMF.md),
  [`spectral_clustering()`](https://lucp12891.github.io/MosaiClusteR/reference/spectral_clustering.md),
  [`LUCID()`](https://lucp12891.github.io/MosaiClusteR/reference/LUCID.md)
  (supervised, outcome-guided).
- [`Cluster()`](https://lucp12891.github.io/MosaiClusteR/reference/Cluster.md)
  provides the single-source baseline primitive.

### Preprocessing, selection and evaluation

- **Preprocessing:**
  [`Normalization()`](https://lucp12891.github.io/MosaiClusteR/reference/Normalization.md)
  and
  [`Distance()`](https://lucp12891.github.io/MosaiClusteR/reference/Distance.md)
  (Euclidean, Manhattan, Tanimoto/Jaccard for binary data,
  correlation-based, and more).
- **Cluster-number selection:**
  [`SelectnrClusters()`](https://lucp12891.github.io/MosaiClusteR/reference/SelectnrClusters.md),
  [`ChooseCluster()`](https://lucp12891.github.io/MosaiClusteR/reference/ChooseCluster.md).
- **Evaluation:**
  [`mosaic_labels()`](https://lucp12891.github.io/MosaiClusteR/reference/mosaic_labels.md),
  [`cluster_agreement()`](https://lucp12891.github.io/MosaiClusteR/reference/cluster_agreement.md)
  (ARI / NMI / Jaccard / purity),
  [`compare_clusterings()`](https://lucp12891.github.io/MosaiClusteR/reference/compare_clusterings.md),
  [`CompareSilCluster()`](https://lucp12891.github.io/MosaiClusteR/reference/CompareSilCluster.md),
  [`CompareSvsM()`](https://lucp12891.github.io/MosaiClusteR/reference/CompareSvsM.md),
  and data-driven source weighting
  ([`DetermineWeight_SilClust()`](https://lucp12891.github.io/MosaiClusteR/reference/DetermineWeight_SilClust.md),
  [`DetermineWeight_SimClust()`](https://lucp12891.github.io/MosaiClusteR/reference/DetermineWeight_SimClust.md)).

### Feature weighting and data nuggets

- A data-nugget compression provides a robust, big-data-friendly feature
  weighting as an alternative to variance:
  [`create_data_nuggets()`](https://lucp12891.github.io/MosaiClusteR/reference/create_data_nuggets.md),
  [`nugget_feature_weights()`](https://lucp12891.github.io/MosaiClusteR/reference/nugget_feature_weights.md),
  [`nugget_cluster()`](https://lucp12891.github.io/MosaiClusteR/reference/nugget_cluster.md),
  [`Wkmeans()`](https://lucp12891.github.io/MosaiClusteR/reference/Wkmeans.md),
  [`Whclust()`](https://lucp12891.github.io/MosaiClusteR/reference/Whclust.md).

### Characterisation, visualisation and interpretation

- **Characterisation / navigation:**
  [`CharacteristicFeatures()`](https://lucp12891.github.io/MosaiClusteR/reference/CharacteristicFeatures.md),
  [`FeatSelection()`](https://lucp12891.github.io/MosaiClusteR/reference/FeatSelection.md),
  [`FeaturesOfCluster()`](https://lucp12891.github.io/MosaiClusteR/reference/FeaturesOfCluster.md),
  [`FindCluster()`](https://lucp12891.github.io/MosaiClusteR/reference/FindCluster.md),
  [`FindElement()`](https://lucp12891.github.io/MosaiClusteR/reference/FindElement.md),
  [`TrackCluster()`](https://lucp12891.github.io/MosaiClusteR/reference/TrackCluster.md),
  [`ReorderToReference()`](https://lucp12891.github.io/MosaiClusteR/reference/ReorderToReference.md),
  [`SharedComps()`](https://lucp12891.github.io/MosaiClusteR/reference/SharedComps.md),
  [`SimilarityMeasure()`](https://lucp12891.github.io/MosaiClusteR/reference/SimilarityMeasure.md).
- **Visualisation:**
  [`ComparePlot()`](https://lucp12891.github.io/MosaiClusteR/reference/ComparePlot.md),
  [`ClusterPlot()`](https://lucp12891.github.io/MosaiClusteR/reference/ClusterPlot.md),
  [`Cyclogram()`](https://lucp12891.github.io/MosaiClusteR/reference/Cyclogram.md),
  [`SimilarityHeatmap()`](https://lucp12891.github.io/MosaiClusteR/reference/SimilarityHeatmap.md),
  [`HeatmapPlot()`](https://lucp12891.github.io/MosaiClusteR/reference/HeatmapPlot.md),
  [`ProfilePlot()`](https://lucp12891.github.io/MosaiClusteR/reference/ProfilePlot.md),
  [`ContFeaturesPlot()`](https://lucp12891.github.io/MosaiClusteR/reference/ContFeaturesPlot.md),
  [`BinFeaturesPlot_SingleData()`](https://lucp12891.github.io/MosaiClusteR/reference/BinFeaturesPlot_SingleData.md)/`_MultipleData()`,
  [`BoxPlotDistance()`](https://lucp12891.github.io/MosaiClusteR/reference/BoxPlotDistance.md),
  and colour helpers
  ([`ColorPalette()`](https://lucp12891.github.io/MosaiClusteR/reference/ColorPalette.md),
  [`ColorsNames()`](https://lucp12891.github.io/MosaiClusteR/reference/ColorsNames.md),
  [`ClusterCols()`](https://lucp12891.github.io/MosaiClusteR/reference/ClusterCols.md),
  [`LabelCols()`](https://lucp12891.github.io/MosaiClusteR/reference/LabelCols.md)).
- **Biological interpretation (omics):**
  [`DiffGenes()`](https://lucp12891.github.io/MosaiClusteR/reference/DiffGenes.md),
  [`PathwayAnalysis()`](https://lucp12891.github.io/MosaiClusteR/reference/PathwayAnalysis.md),
  [`Pathways()`](https://lucp12891.github.io/MosaiClusteR/reference/Pathways.md),
  [`Geneset.intersect()`](https://lucp12891.github.io/MosaiClusteR/reference/Geneset.intersect.md),
  [`SharedGenesPathsFeat()`](https://lucp12891.github.io/MosaiClusteR/reference/SharedGenesPathsFeat.md)
  and related functions connect clusters to genes and pathways.

### Data and tests

- Bundled `mosaic_toy` dataset and the
  [`mosaic_sim()`](https://lucp12891.github.io/MosaiClusteR/reference/mosaic_sim.md)
  structured-data generator.
- A `testthat` suite covering distances, weighting schemes, the
  integration methods, evaluation metrics and data nuggets.

### Notes

- A compiled kernel (`src/`) accelerates the consensus step; building
  from source needs a C/C++ toolchain (Rtools on Windows). Heavier and
  specialised dependencies — `SNFtool`, `circlize`, `plotrix`,
  `ggplot2`, `igraph`, `limma`, `MLP`, `biomaRt`, `org.Hs.eg.db` and
  others — are declared under *Suggests* and loaded only when the
  corresponding method is called.
- A few methods rely on an external graph-partitioning executable
  (e.g. METIS for
  [`EnsembleClustering()`](https://lucp12891.github.io/MosaiClusteR/reference/EnsembleClustering.md),
  the METIS path of
  [`EHC()`](https://lucp12891.github.io/MosaiClusteR/reference/EHC.md),
  and `LinkBasedClustering(linkBasedMethod = "asrs")`) and raise an
  informative error when it is absent; pure-R alternatives
  ([`HBGF()`](https://lucp12891.github.io/MosaiClusteR/reference/HBGF.md),
  [`EvidenceAccumulation()`](https://lucp12891.github.io/MosaiClusteR/reference/EvidenceAccumulation.md),
  CTS/SRS) cover the same need.
