# MosaiClusteR 1.0.0 (development)

Work toward the first stable release, on branch `dev-1.0.0`. This section
collects changes since 0.1.4.

## Correctness fixes

* `Cluster()` now forwards the normalisation method to `Distance()`. Previously
  the `method` argument was documented and accepted but silently ignored, and
  `normalize = TRUE` did not select a method. `normalize` now accepts either a
  logical (`FALSE` for none) or a normalisation method name, and `method` sets
  the method when `normalize = TRUE`.
* Documentation examples now use the objects-in-rows orientation consistently
  (the `Cluster()` example no longer transposes the bundled `mosaic_toy` data).

# MosaiClusteR 0.1.4

Third public release. It consolidates the development carried out after 0.1.1,
which had not previously been published.

## Breaking change

* **`M_ABCpp()` renamed to `M_ABC()`.** The multi-source Aggregating Bundles of
  Clusters function is now `M_ABC()`; its arguments, behaviour and return value
  are unchanged. Update any calls from `M_ABCpp(...)` to `M_ABC(...)`.

## New features

* The example dataset `mosaic_toy` is bundled in `data/` and used by the
  examples and tests, so `R CMD check` (and the r-universe builds) no longer
  fail with "object 'mosaic_toy' not found".
* **`plot_pathways(style = "circular")`.** A radial (sunburst-style) bar plot
  that groups the enriched terms by cluster around a central hub. Bar length
  encodes \eqn{-\log_{10}} p and colour encodes the cluster; each term is a
  single-line radial label that stays legible as the number of terms grows.
  Complements the existing `"dotplot"` and `"lollipop"` styles. It also gained
  `cluster_cols` (an explicit per-cluster palette, e.g. to match `Cyclogram()`)
  and `wrap` (radial label wrapping).
* **`enriched_features()`.** A helper that turns `CharacteristicFeatures()`
  output into a per-cluster list of compact `term`/`pvalue`/`Count` tables,
  ready to hand straight to `plot_pathways()`.
* `FeaturesOfCluster()` and `BinFeaturesPlot_SingleData()` gained `cexB`, `cexR`
  and `cexN` arguments for the object, feature and name label sizes; `Cyclogram()`
  now defaults to a larger `cex.labels`.

## Bug fixes

* **`CharacteristicFeatures()`.** When called with `binData` (the common case,
  as in the documented example) an inner loop index overwrote the cluster index,
  so only the first cluster's characteristic features were returned, and a
  `datanames` mis-index left the per-source result list unnamed. Both are fixed:
  features are now computed and named for every cluster.
* **`Distance(distmeasure = "euclidean")`** now computes a true Euclidean
  distance (via `stats::dist`). Previously it routed to `FD::gowdis` (Gower),
  which returned `NaN` when a feature was constant across the selected objects,
  producing NA distance matrices. Gower is still available via
  `distmeasure = "gower"`.
* **`ComparePlot()`** (and `Cyclogram()`) now build a default colour palette when
  `cols` is not supplied, and check up front that all results describe the same
  objects, raising a clear message otherwise.

## Metadata

* The licence is **GPL-3**, consistent with 0.1.0 and 0.1.1. An interim local
  change to CC BY-SA 4.0 is reverted, since CC licences are intended for content
  rather than software.
* "MoSaIC" expands to *Multi-Omics Similarity Aggregation and Integrative
  Clustering in R*, and the README framework overview is a scalable SVG figure
  (`man/figures/framework.svg`) that reads on both light and dark backgrounds.

# MosaiClusteR 0.1.1

* Documentation and installation updates to the README, and a package version
  bump from 0.1.0.

# MosaiClusteR 0.1.0

First public release. MosaiClusteR is an umbrella framework for multi-source and
multi-omics clustering: it gathers a broad catalogue of integration methods
behind one consistent interface, together with the tooling needed to preprocess,
select, evaluate, compare, characterise and biologically interpret the results.

## Framework

* **Unified interface.** Every method takes the same input — a `List` of numeric
  matrices over the same objects (objects in rows) — and returns the same
  `{DistM, Clust}` contract, so methods are interchangeable and directly
  comparable. Precomputed dissimilarities are accepted via `type = "dist"`.
* **Four-layer workflow.** Preprocess → baseline → integrate → evaluate.

## Integration methods (five paradigms)

* **Direct:** `ADC()`, `ADEC()`.
* **Similarity-based:** `WeightedClust()`, `WonM()`, `SNF()`, `NEMO()`.
* **Graph-based:** `EnsembleClustering()` (CSPA/HGPA/MCLA), `HBGF()`,
  `ClusteringAggregation()`.
* **Voting / consensus:** `CEC()`, `CVAA()`, `ConsensusClustering()`,
  `EvidenceAccumulation()`, `LinkBasedClustering()` (CTS/SRS/ASRS), and the
  Aggregating Bundles of Clusters family `M_ABC()` / `M_ABCdist()` /
  `M_ABCdist.WC()` / `M_ABCdeep()` (a deep-learning, autoencoder-based variant).
* **Hierarchy-based:** `HierarchicalEnsembleClustering()`, `EHC()`.
* **Factor / low-rank / spectral:** `intNMF()`, `spectral_clustering()`,
  `LUCID()` (supervised, outcome-guided).
* `Cluster()` provides the single-source baseline primitive.

## Preprocessing, selection and evaluation

* **Preprocessing:** `Normalization()` and `Distance()` (Euclidean, Manhattan,
  Tanimoto/Jaccard for binary data, correlation-based, and more).
* **Cluster-number selection:** `SelectnrClusters()`, `ChooseCluster()`.
* **Evaluation:** `mosaic_labels()`, `cluster_agreement()` (ARI / NMI / Jaccard /
  purity), `compare_clusterings()`, `CompareSilCluster()`, `CompareSvsM()`,
  and data-driven source weighting (`DetermineWeight_SilClust()`,
  `DetermineWeight_SimClust()`).

## Feature weighting and data nuggets

* A data-nugget compression provides a robust, big-data-friendly feature
  weighting as an alternative to variance: `create_data_nuggets()`,
  `nugget_feature_weights()`, `nugget_cluster()`, `Wkmeans()`, `Whclust()`.

## Characterisation, visualisation and interpretation

* **Characterisation / navigation:** `CharacteristicFeatures()`,
  `FeatSelection()`, `FeaturesOfCluster()`, `FindCluster()`, `FindElement()`,
  `TrackCluster()`, `ReorderToReference()`, `SharedComps()`, `SimilarityMeasure()`.
* **Visualisation:** `ComparePlot()`, `ClusterPlot()`, `Cyclogram()`,
  `SimilarityHeatmap()`, `HeatmapPlot()`, `ProfilePlot()`, `ContFeaturesPlot()`,
  `BinFeaturesPlot_SingleData()`/`_MultipleData()`, `BoxPlotDistance()`, and
  colour helpers (`ColorPalette()`, `ColorsNames()`, `ClusterCols()`, `LabelCols()`).
* **Biological interpretation (omics):** `DiffGenes()`, `PathwayAnalysis()`,
  `Pathways()`, `Geneset.intersect()`, `SharedGenesPathsFeat()` and related
  functions connect clusters to genes and pathways.

## Data and tests

* Bundled `mosaic_toy` dataset and the `mosaic_sim()` structured-data generator.
* A `testthat` suite covering distances, weighting schemes, the integration
  methods, evaluation metrics and data nuggets.

## Notes

* A compiled kernel (`src/`) accelerates the consensus step; building from source
  needs a C/C++ toolchain (Rtools on Windows). Heavier and specialised
  dependencies — `SNFtool`, `circlize`, `plotrix`, `ggplot2`, `igraph`, `limma`,
  `MLP`, `biomaRt`, `org.Hs.eg.db` and others — are declared under *Suggests* and
  loaded only when the corresponding method is called.
* A few methods rely on an external graph-partitioning executable (e.g. METIS for
  `EnsembleClustering()`, the METIS path of `EHC()`, and
  `LinkBasedClustering(linkBasedMethod = "asrs")`) and raise an informative error
  when it is absent; pure-R alternatives (`HBGF()`, `EvidenceAccumulation()`,
  CTS/SRS) cover the same need.
