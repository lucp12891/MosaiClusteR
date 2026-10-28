# MosaiClusteR 0.1.5

* The example dataset `mosaic_toy` is bundled in `data/` and used by the
  examples and tests, so `R CMD check` (and the r-universe builds) no longer
  fail with "object 'mosaic_toy' not found".

* **`CharacteristicFeatures()` bug fix.** When called with `binData` (the common
  case, as in the documented example) an inner loop index overwrote the cluster
  index, so only the first cluster's characteristic features were returned, and
  a `datanames` mis-index left the per-source result list unnamed. Both are
  fixed: features are now computed and named for every cluster.
* **`plot_pathways(style = "circular")`.** New radial (sunburst-style) bar plot
  that groups the enriched terms by cluster around a central hub. Bar length
  encodes \eqn{-\log_{10}} p and colour encodes the cluster; each term is a
  single-line radial label that stays legible as the number of terms grows.
  Complements the existing `"dotplot"` and `"lollipop"` styles.
* **`enriched_features()`** - new helper that turns [CharacteristicFeatures()]
  output into a per-cluster list of compact `term`/`pvalue`/`Count` tables
  (keeping the significant, enriched features and tidying their names), ready to
  hand straight to `plot_pathways()`.
* `plot_pathways(style = "circular")` gained `cluster_cols` to supply an explicit
  per-cluster palette (e.g. to match the colours used by `Cyclogram()`), and
  `wrap` to control radial label wrapping.
* `FeaturesOfCluster()` and `BinFeaturesPlot_SingleData()` gained `cexB`, `cexR`
  and `cexN` arguments to control the object, feature and name label sizes;
  `Cyclogram()` now defaults to a larger `cex.labels`.

# MosaiClusteR 0.1.2

* `Distance(distmeasure = "euclidean")` now computes a **true Euclidean**
  distance (via `stats::dist`). Previously it routed to `FD::gowdis` (Gower),
  which returned `NaN` when a feature was constant across the selected objects,
  producing NA distance matrices. Gower is still available via
  `distmeasure = "gower"`.
* `ComparePlot()` (and `Cyclogram()`) now build a default colour palette when
  `cols` is not supplied. Previously the default `cols = NULL` caused an error
  ("number of items to replace is not a multiple of replacement length") from
  `ColorsNames()`, so a palette had to be passed explicitly.
* `ComparePlot()` now checks up front that all results describe the **same
  objects** and raises a clear message otherwise, instead of a cryptic
  "missing value where TRUE/FALSE needed" error. (Results from different data
  sets / sample sizes cannot be compared.)

# MosaiClusteR 0.1.1

* **`M_ABCpp()` renamed to `M_ABC()`.** The multi-source Aggregating Bundles of
  Clusters function is now `M_ABC()`; its arguments, behaviour and return value
  are unchanged. Update any calls from `M_ABCpp(...)` to `M_ABC(...)`.
* Metadata and documentation: the "MoSaIC" name now expands to *Multi-Omics
  Similarity Aggregation and Integrative Clustering in R*, the package licence is
  **CC BY-SA 4.0**, and the four-layer framework overview in the README is now a
  scalable SVG figure (`man/figures/framework.svg`) that reads on both light and
  dark backgrounds.

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
