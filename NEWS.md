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
  Aggregating Bundles of Clusters family `M_ABCpp()` / `M_ABCdist()` /
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
