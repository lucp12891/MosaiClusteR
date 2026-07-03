# MosaiClusteR

<!-- badges: start -->
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R >= 4.0](https://img.shields.io/badge/R-%3E%3D4.0-blue)](https://www.r-project.org/)
[![Lifecycle: maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/)
![Status](https://img.shields.io/badge/tests-passing-brightgreen.svg)
<!-- badges: end -->

> **An umbrella framework for multi-source & multi-omics clustering in R.**
> *One interface, many methods, a complete workflow — from raw modalities to interpretable subgroups.*

**MosaiClusteR** — *“**MoSaIC**” = **M**ulti-**O**mics **S**ource-**A**gnostic **I**ntegration **C**lustering in R* — treats every data modality (transcriptome, methylome, proteome, fingerprints, targets, imaging, clinical…) as one **tile** of a larger picture, and assembles those tiles into a single coherent *mosaic* of sample structure. It unifies a large family of integrative clustering algorithms behind **one consistent list‑of‑matrices interface**, and wraps them in a complete analysis framework: preprocessing → single‑source baselines → integration across five paradigms → evaluation & interpretation.

No single algorithm is best across all data regimes, so MosaiClusteR gathers many methods on the same footing and gives you the tooling to *choose between them*: cluster‑number selection, agreement metrics, visual comparison, cluster characterisation, and — for omics — pathway and gene‑set interpretation. It also introduces a **data‑nugget feature‑weighting** scheme that offers a robust, Big‑Data‑friendly alternative to variance weighting.

---

## Table of contents

- [Why MosaiClusteR](#why-mosaicluster)
- [Installation](#installation)
- [60‑second quick start](#60-second-quick-start)
- [The MoSaIC framework of use](#the-mosaic-framework-of-use)
- [Method catalogue](#method-catalogue) — *every method, and when to use it*
- [Spotlight: M‑ABC and data nuggets](#spotlight-m-abc-and-data-nuggets)
- [Evaluation & comparison toolkit](#evaluation--comparison-toolkit)
- [Datasets](#datasets)
- [Design principles](#design-principles)
- [Citation](#citation)
- [References](#references)

---

## Why MosaiClusteR

Multi‑omics and other multi‑source studies measure the *same objects* (patients, tumours, compounds, cells) through several heterogeneous lenses. No single modality tells the whole story, and the many integration algorithms in the literature live in scattered packages with incompatible inputs and outputs. MosaiClusteR fixes this:

- **One input format.** Every method takes a `List` of matrices over the *same shared objects*, with **objects (samples) in rows and features in columns** (`n × m`) — the same convention across the entire package, including the M‑ABC family. Distance‑matrix inputs (`type = "dist"`) are orientation‑free.
- **One output contract.** Every method returns at least `DistM` (a dissimilarity matrix) and `Clust` (a hierarchical‑clustering object), so results are directly comparable and composable.
- **Five paradigms, one roof.** Direct, similarity‑based, graph‑based, voting‑based consensus, and hierarchy‑based integration — pick the right tool without leaving the package.
- **A full workflow, not just algorithms.** Normalisation, distances, partition extraction, agreement metrics (ARI/NMI/Jaccard), and method comparison are built in.

---

## Installation

```r
# install.packages("remotes")
remotes::install_github("bosangir/MosaiClusteR")
```

MosaiClusteR ships C++ kernels for the M‑ABC family, so a C/C++ toolchain is
needed to build from source (on Windows install
[Rtools](https://cran.r-project.org/bin/windows/Rtools/); on macOS the Xcode
command‑line tools). Every method has a pure‑R fallback, so the package still
runs if a kernel is unavailable. A few methods rely on *suggested* packages that
are only loaded when you call them:

```r
# optional, enables extra methods / distances
install.packages(c("SNFtool", "datanugget", "WCluster", "lsa", "analogue", "FactoMineR"))
```

| If you want…                              | Install (Suggests)          |
|-------------------------------------------|-----------------------------|
| Similarity Network Fusion (`SNF`)         | `SNFtool`                   |
| The reference data‑nugget engine          | `datanugget`, `WCluster`    |
| Cosine / χ² / MCA distances               | `lsa`, `analogue`, `FactoMineR` |

---

## 60‑second quick start

```r
library(MosaiClusteR)
data(mosaic_toy)                       # two sources, 60 shared objects, 3 planted groups
L     <- mosaic_toy$List               # a list of matrices, objects in rows
truth <- mosaic_toy$truth

## ---- Integrate: every method shares one interface — swap to change paradigm
fit <- SNF(L, type = "data", distmeasure = c("euclidean", "euclidean"),
           normalize = c(FALSE, FALSE), method = list(NULL, NULL))

labels <- mosaic_labels(fit, k = 3)    # pull a flat partition from any fitted result
cluster_agreement(labels, truth)       # ARI / NMI / Jaccard / purity vs ground truth

## ---- Same score, a different paradigm — e.g. a neighbourhood method
cluster_agreement(NEMO(L, k = 3)$cluster, truth)

## ---- Compare two integrations directly
compare_clusterings(fit$DistM,
                    ADC(L, distmeasure = "euclidean", normalize = FALSE)$DistM, NC = 3)
```

That is the whole loop: **integrate → extract → evaluate → compare** — and the
same three lines work for any of the methods below.

---

## The MoSaIC framework of use

MosaiClusteR is organised as a four‑layer pipeline. You can enter at any layer,
and every layer speaks the same `List → {DistM, Clust}` language.

```
            ┌─────────────────────────────────────────────────────────────┐
  LAYER 1   │  PREPROCESSING                                               │
  prepare   │  Normalization()  ·  Distance()                             │
            │  heterogeneous, differently-scaled modalities → comparable  │
            └─────────────────────────────────────────────────────────────┘
                              │  List of matrices (features × objects)
                              ▼
            ┌─────────────────────────────────────────────────────────────┐
  LAYER 2   │  SINGLE-SOURCE BASELINES                                     │
  explore   │  cluster each tile alone → see what each modality recovers  │
            │  motivates integration: agreements vs. modality-specific    │
            └─────────────────────────────────────────────────────────────┘
                              │
                              ▼
            ┌─────────────────────────────────────────────────────────────┐
  LAYER 3   │  INTEGRATION ENGINE  (five paradigms)                        │
  integrate │  Direct · Similarity · Graph · Voting-consensus · Hierarchy │
            │  M_ABCpp · M_ABCdist · WeightedClust · SNF · CEC · HEC · …   │
            └─────────────────────────────────────────────────────────────┘
                              │  {DistM, Clust}
                              ▼
            ┌─────────────────────────────────────────────────────────────┐
  LAYER 4   │  EVALUATION & INTERPRETATION                                 │
  conclude  │  mosaic_labels() · cluster_agreement() · compare_clusterings│
            │  ARI / NMI / Jaccard · cross-method consensus               │
            └─────────────────────────────────────────────────────────────┘
```

**Choosing a layer‑3 method** (decision shortcut):

- *Want a fast, transparent baseline?* → **Direct** (`ADC`) or **Weighted** (`WeightedClust`).
- *Sources are noisy and complementary, genomic scale?* → **Similarity** (`SNF`).
- *Want stability under resampling / a novel robust ensemble?* → **Voting consensus** (`M_ABCpp`, `CEC`, `EvidenceAccumulation`).
- *Care about the hierarchy / dendrogram structure itself?* → **Hierarchy** (`HEC`).
- *Have very large N or heavy‑tailed data?* → use **`weighting = "nugget"`** inside the M‑ABC family.

---

## Method catalogue

> **Legend** —  ✅ implemented & tested in this package · † requires an external
> graph‑partitioning executable (METIS) or MATLAB, not bundled · 🔭 interoperable
> with / on the roadmap via a reference package. All implemented methods share the
> `List → {DistM, Clust}` contract and take the same objects‑in‑rows input.

### 1 · Direct clustering — *combine first, then cluster*

| Method | Function | What it does | When it’s useful | Ref |
|---|---|---|---|---|
| **ADC** — Aggregated Data Clustering | `ADC` ✅ | Concatenate all tiles `[D₁│…│Dₗ]`, one distance, one hierarchy | Simplest baseline; when modalities are commensurate and you want a transparent reference | Fodeh 2013 |
| **ADEC** — Aggregated Data *Ensemble* Clustering | `ADEC` ✅ | ADC + resampling of features / cut‑points → co‑clustering matrix | When ADC is unstable and you want a robustness layer | Fodeh 2013 |

### 2 · Similarity‑based — *fuse object‑by‑object similarity*

| Method | Function | What it does | When it’s useful | Ref |
|---|---|---|---|---|
| **Weighted Clustering** | `WeightedClust` ✅ | Convex combination `Σ wₖ Dₖ` over a weight grid | When you want to *tune and inspect* the trade‑off between sources | Perualila‑Tan 2016 |
| **SNF** — Similarity Network Fusion | `SNF` ✅ | Cross‑diffuse per‑source similarity networks to one fused network | Noisy, complementary, genomic‑scale sources; the multi‑omics workhorse | Wang 2014 |
| **WonM** — Weighting on Membership | `WonM` ✅ | Consensus co‑membership across many cut‑points *k* | When you don’t want to commit to a single *k* | — |
| **NEMO** — Neighborhood‑based Multi‑Omics | `NEMO` ✅ | Locally‑scaled kNN affinities averaged over the omics measuring each pair → spectral clustering | When samples are **missing some modalities** (no imputation needed) | Rappoport 2019 |
| **Spectral clustering** | `spectral_clustering` ✅ | Normalised‑Laplacian embedding + k‑means on any affinity | Final clusterer for fused/similarity matrices (SNF, NEMO, Spectrum) | Ng 2002 |
| **ab‑SNF** | 🔭 | SNF with association/variance feature weights | When a few features carry the signal and should drive the network | Zhang 2022 |
| **CIMLR / rMKL‑LPP** | 🔭 | Multi‑kernel similarity learning | Many heterogeneous kernels; learn their weights | Zhang 2022 |

### 3 · Graph‑based — *partition an ensemble (hyper)graph*

| Method | Function | What it does | When it’s useful | Ref |
|---|---|---|---|---|
| **CSPA / HGPA / MCLA** | `EnsembleClustering` ✅† | Cluster‑ensemble via (hyper)graph partitioning | Combining many base partitions into one robust consensus | Strehl 2002 |
| **HBGF** — Hybrid Bipartite Graph Formulation | `HBGF` ✅ | Partition an object↔cluster bipartite graph (SVD + k‑means) | Ensembles where objects and clusters co‑embed naturally | Fern 2004 |
| **Clustering Aggregation** (Balls/Agglo/Furthest) | `ClusteringAggregation` ✅ | Minimise pairwise disagreement cost | When you have several partitions and want the “median” one | Gionis 2007 |

### 4 · Voting‑based consensus — *let partitions vote*

| Method | Function | What it does | When it’s useful | Ref |
|---|---|---|---|---|
| **M‑ABC** — Multi‑source Aggregating Bundles of Clusters | `M_ABCpp` ✅ | Per‑source bootstrap of *feature bundles* → consensus co‑clustering dissimilarity | High‑dimensional sources, where signal hides in feature subsets; robust by design; supports **data‑nugget weighting** | Amaratunga 2008 |
| **M‑ABCdist** | `M_ABCdist` ✅ | M‑ABC that accumulates *distances* (keeps geometry) then fuses | When you want continuous dissimilarities rather than co‑clustering counts | Amaratunga 2008 |
| **M‑ABCdeep** — deep‑learning ABC | `M_ABCdeep` / `ABCdeep.SingleInMultiple` ✅ | M‑ABC whose base clustering is a neural autoencoder embedding + latent clustering; self‑contained R network, no Python/torch | Nonlinear/complex per‑source structure; when a learned embedding beats a raw distance | Amaratunga 2008 |
| **CEC** — (Consensus) Ensemble Clustering | `CEC` ✅ | Incidence accumulation across sources & cut‑points | A solid, weight‑aware consensus baseline | Fodeh 2013 |
| **CVAA / W‑CVAA** | `CVAA` ✅ | Cumulative voting aggregation (and weighted) | Aligning and merging many partitions by voting | Saeed 2012/2014 |
| **IVC / IPVC / IPC** | `ConsensusClustering` ✅ | Iterative (probabilistic) voting consensus | Iteratively refine a consensus labelling | Nguyen 2007 |
| **EA** — Evidence Accumulation | `EvidenceAccumulation` ✅ | Co‑association matrix across partitions (SL / SL‑agnes / MST) | Classic, parameter‑light consensus | Fred 2005 |
| **CTS / SRS / ASRS** | `LinkBasedClustering` ✅ | Link‑based cluster ensembles (CTS/SRS pure‑R; ASRS needs MATLAB) | When pairwise *link* similarity refines the consensus | Iam‑On 2010 |

### 5 · Hierarchy‑based — *consense the trees themselves*

| Method | Function | What it does | When it’s useful | Ref |
|---|---|---|---|---|
| **HEC** — Hierarchical Ensemble Clustering | `HierarchicalEnsembleClustering` ✅ | Aggregate cophenetic distances → closest ultrametric (Floyd–Warshall) | When the *dendrogram structure* across sources is the object of interest | Zheng 2014 |
| **EHC** — Ensemble for Hierarchical Clustering | `EHC` ✅† | Graph aggregation of dendrogram association strengths | Graph‑flavoured hierarchical consensus | Hossain 2012 |

### 6 · Model‑based & factor integration

| Method | Function | What it does | When it’s useful | Ref |
|---|---|---|---|---|
| **intNMF** — integrative NMF | `intNMF` ✅ | Shared non‑negative basis `W` + per‑omic `Hₖ`; consensus clustering of `W` | Low‑rank, parts‑based integration; non‑negative data | Chalise 2017 |
| **LUCID** | `LUCID` ✅‡ | Latent clusters jointly from exposures, omics **and** an outcome (quasi‑mediation, EM) | Supervised subtyping where an outcome should shape the clusters | Zhao 2024 |
| **iCluster / iClusterPlus / iClusterBayes** | 🔭 | Joint latent‑variable integration, mixed data types, feature selection | Zhang 2022 |
| **moCluster / MOFA / JIVE** | 🔭 | Low‑rank / factor decomposition of shared + individual variation | Zhang 2022 |

‡ `LUCID` wraps the **LUCIDus** package (a *Suggests* dependency).

### Data reduction & weighting

| Tool | Function | Role |
|---|---|---|
| **Data nuggets** | `create_data_nuggets`, `nugget_feature_weights` ✅ | Compress N observations into M weighted representatives; derive robust **feature weights** for M‑ABC |
| **Data‑nugget clustering** | `nugget_cluster`, `Wkmeans`, `Whclust` ✅ | Cluster the compressed representatives directly (weighted k‑means / weighted Ward), then map back to objects — Big‑Data clustering | 

---

## Spotlight: M‑ABC and data nuggets

**M‑ABC** generalises the single‑source *Aggregating Bundles of Clusters* algorithm (Amaratunga et al. 2008) to many heterogeneous sources. For each source it repeats, `numsim` times: bootstrap the samples, draw a **weighted subsample of features**, cluster, and record co‑clustering. Per‑source evidence is accumulated into one consensus dissimilarity, then clustered.

The feature subsample is the heart of the method, and *how features are weighted* matters:

| `weighting` | Score behind the selection probability | Best when |
|---|---|---|
| `"var"` | feature **variance** (classic) | signal lives in high‑variance features |
| `"cv"`  | coefficient of variation | scale‑free relative dispersion |
| `"nugget"` | **data‑nugget between‑group weighted variance** | large N, heavy tails, outliers — robust signal |
| `"equal"` | flat | ablation / no prior |

**Why data nuggets?** A data nugget summarises a group of observations by a
*centre*, a *weight* (how many observations it represents) and a *scale* (its
internal spread). Building nuggets over the samples and asking *“which features
separate the nugget centres?”* yields a feature‑importance score that rewards
real between‑group signal while discounting within‑nugget noise — and it scales
to very large N because it works on `M ≪ N` representatives.

```r
# Inspect the nugget compression and the weights it induces
dn <- create_data_nuggets(t(L[[1]]), max_nuggets = 25)
dn                                   # 60 obs -> 25 nuggets (reduction reported)
w  <- nugget_feature_weights(dn, type = "between")   # robust feature weights

# Drop straight into M-ABC
M_ABCpp(L, weighting = c("nugget", "nugget"),
      nugget_args = list(max_nuggets = 25, seed = 1))
```

If the reference `datanugget` package is installed you can swap in its engine
with `create_data_nuggets(x, engine = "datanugget")`; otherwise the built‑in
native engine is used (no extra dependency).

---

## Evaluation & comparison toolkit

```r
mosaic_labels(fit, k = 3)                 # extract a partition from any result
cluster_agreement(a, b)                   # ARI, NMI, pairwise Jaccard
compare_clusterings(D1, D2, NC = 3)       # cluster two dissimilarities & compare
```

`cluster_agreement()` rule of thumb (ARI): `>0.90` excellent · `0.75–0.90` good
· `0.50–0.75` moderate · `<0.50` weak.

### Downstream analysis & visualisation

The package also ships a full downstream analysis suite (Layer 4). These
guard their heavier dependencies via *Suggests* and raise an informative error
if one is missing.

| Purpose | Functions | Needs |
|---|---|---|
| Cross‑method comparison | `ReorderToReference`, `SimilarityMeasure`, `CompareSilCluster`, `CompareSvsM`, `CompareInteractive`, `SelectnrClusters` | base / `plotrix` |
| Source weighting | `DetermineWeight_SilClust`, `DetermineWeight_SimClust` | base |
| Cluster navigation | `FindCluster`, `FindElement`, `SharedComps`, `TrackCluster`, `ChooseCluster` | base |
| Characterisation | `CharacteristicFeatures`, `FeatSelection`, `FeaturesOfCluster` | base |
| Dendrograms / colours | `ComparePlot`, `Cyclogram`, `ClusterPlot`, `ColorPalette`, `LabelPlot` | `circlize`, `plotrix` |
| Heatmaps / boxplots | `HeatmapPlot`, `HeatmapSelection`, `SimilarityHeatmap`, `BoxPlotDistance` | `gplots`, `ggplot2`, `gridExtra` |
| Feature / profile plots | `BinFeaturesPlot_SingleData`, `BinFeaturesPlot_MultipleData`, `ContFeaturesPlot`, `ProfilePlot` | `plotrix`, `Biobase` |
| Differential expression | `DiffGenes`, `DiffGenesSelection`, `FindGenes` | `limma` |
| Pathway analysis | `PathwayAnalysis`, `Pathways`, `PathwaysIter`, `PathwaysSelection`, `PreparePathway`, `PlotPathways`, `Geneset.intersect` | `MLP`, `biomaRt`, `org.Hs.eg.db` |
| Shared genes/paths/features | `SharedGenesPathsFeat`, `SharedSelection`, `SharedSelectionLimma`, `SharedSelectionMLP` | base / `limma` / `MLP` |

---

## Datasets

- **`mosaic_toy`** — a bundled two‑source example (60 shared samples, 3 planted
  groups) for instant experimentation and the test suite.
- **`mosaic_sim()`** — generate your own multi‑source data with controllable
  group structure, informative‑feature counts and effect size.

The methods are designed for, and have been evaluated on, multi‑omics
benchmarks in the spirit of the **LUCIDus** simulated‑HELIX data, **TCGA**
subtype cohorts, and the **MCF7** compound–fingerprint–transcriptome panel used
in the companion manuscript.

---

## Design principles

1. **Source‑agnostic.** Methods do not care whether a tile is RNA, methylation, fingerprints or clinical scores.
2. **Composable.** A uniform `{DistM, Clust}` return makes any two methods comparable and any result re‑usable.
3. **Honest dependencies.** The core is pure R; heavyweight methods degrade gracefully when an optional package is absent.
4. **Tested.** A `testthat` suite covers the data‑nugget engine, every weighting scheme, the M‑ABC family, distances and metrics.

---

## Citation

```bibtex
@Manual{MosaiClusteR2026,
  title  = {MosaiClusteR: An Umbrella Framework for Multi-Source and
            Multi-Omics Clustering},
  author = {Osang'ir, Bernard Isekah and Cabrera, Javier and
            Amaratunga, Dhammika and Shkedy, Ziv},
  year   = {2026},
  note   = {R package version 0.1.0},
  url    = {https://github.com/bosangir/MosaiClusteR}
}
```

## References

- Amaratunga, Cabrera & Shkedy (2008). *Exploration and Analysis of DNA Microarray and Other High‑Dimensional Data.* Wiley.
- Cherasia, Cabrera, Fernholz & Fernholz (2023). *Data Nuggets in Supervised Learning.* In *Robust and Multivariate Statistical Methods*, Springer, 429–449.
- Wang et al. (2014). *Similarity Network Fusion…* **Nature Methods** 11:333–337.
- Perualila‑Tan et al. (2016). *Weighted similarity‑based clustering…*
- Fodeh, Punch & Tan (2013). *On unifying multi‑source information…* **Proteins**.
- Zheng, Chen & Zheng (2014). *A novel approach for hierarchical ensemble clustering.* **IEEE TPAMI**.
- Strehl & Ghosh (2002); Fern & Brodley (2004); Gionis et al. (2007); Nguyen & Caruana (2007); Fred & Jain (2005); Iam‑On et al. (2010); Saeed et al. (2012/2014); Hossain et al. (2012).
- Zhang, Zhou, Xu & Liu (2022). *Integrative clustering methods for multi‑omics data.* **WIREs Comput Stat** 14:e1553. *(SNF, NEMO, iCluster family, ab‑SNF, CIMLR review.)*
- Zhao, Jia, Goodrich & Conti (2024). *LUCIDus…* **The R Journal** 16(2).

---

*MosaiClusteR is an open, extensible framework — contributions of new integration methods, distances, and evaluation tools are very welcome.*
