# Choosing a method: a side-by-side comparison

No single integration method is best across all data regimes.
MosaiClusteR gathers many methods behind one interface so you can try
several on the same footing and let the evidence decide. This page is
the shortcut: it lays out what each method assumes, what it accepts, and
when to reach for it.

Every method takes the same input, a `List` of matrices over the **same
shared objects** with **objects in rows and features in columns**, and
returns the same `{DistM, Clust}` contract. So the choice below is only
about behaviour, never about reformatting your data.

## The five paradigms at a glance

| Paradigm | Idea in one line | Representative functions | Reach for it when |
|----|----|----|----|
| **Direct** | Combine the sources, then cluster once | `ADC`, `ADEC` | You want a transparent baseline and the modalities are commensurate |
| **Similarity-based** | Fuse per-source object similarities | `SNF`, `WeightedClust`, `WonM`, `NEMO` | Sources are noisy and complementary, at genomic scale |
| **Graph-based** | Partition a graph built from many clusterings | `EnsembleClustering`, `HBGF`, `ClusteringAggregation` | You already have several base partitions to consense |
| **Voting / consensus** | Let many partitions vote on co-membership | `M_ABC`, `CEC`, `CVAA`, `ConsensusClustering`, `EvidenceAccumulation`, `LinkBasedClustering` | You want stability under resampling, or signal hidden in feature subsets |
| **Hierarchy-based** | Consense the dendrograms themselves | `HierarchicalEnsembleClustering`, `EHC` | The tree structure across sources is the object of interest |
| **Factor / low-rank** | Cluster a shared latent representation | `intNMF`, `spectral_clustering`, `LUCID` | You want a low-rank model, or an outcome should guide the clusters |

## Property matrix

How to read the columns:

- **Input**: `data` accepts raw matrices (a distance is computed
  internally); `dist` accepts a precomputed dissimilarity
  (`type = "dist"`).
- **Missing a modality**: whether a sample measured in only some sources
  can be included without imputation.
- **Reproducibility**: `deterministic` gives the same result every run;
  [`set.seed()`](https://rdrr.io/r/base/Random.html) or a `seed =`
  argument fixes the result of a stochastic method.
- **Extra dependency**: beyond base R. A dagger (†) marks an external
  executable (METIS) or MATLAB that is not bundled.

| Method | Paradigm | Input | Missing a modality | Reproducibility | Extra dependency |
|----|----|----|----|----|----|
| `ADC` | Direct | data, dist | No (needs the same objects) | Deterministic | none |
| `ADEC` | Direct | data, dist | No | `seed =` | none |
| `WeightedClust` | Similarity | data, dist | No | Deterministic | none |
| `SNF` | Similarity | data | No | Deterministic | `SNFtool` |
| `WonM` | Similarity | data, dist | No | Deterministic | none |
| `NEMO` | Similarity | data | **Yes** (averages over the omics that measured each pair) | [`set.seed()`](https://rdrr.io/r/base/Random.html) (k-means step) | none |
| `spectral_clustering` | Factor | affinity matrix | n/a | [`set.seed()`](https://rdrr.io/r/base/Random.html) (k-means step) | none |
| `EnsembleClustering` | Graph | base partitions | No | Deterministic | METIS † |
| `HBGF` | Graph | base partitions | No | [`set.seed()`](https://rdrr.io/r/base/Random.html) (k-means step) | none |
| `ClusteringAggregation` | Graph | base partitions | No | Deterministic | none |
| `CEC` | Voting | data, dist | No | Deterministic | none |
| `CVAA` | Voting | base partitions | No | Deterministic | none |
| `ConsensusClustering` | Voting | base partitions | No | [`set.seed()`](https://rdrr.io/r/base/Random.html) | none |
| `EvidenceAccumulation` | Voting | base partitions | No | Deterministic | none |
| `LinkBasedClustering` | Voting | base partitions | No | Deterministic (CTS, SRS); ASRS needs MATLAB † | none / MATLAB † |
| `M_ABC` | Voting | data | No | `seed =` | none (C++ kernel built in) |
| `M_ABCdist` | Voting | data | No | [`set.seed()`](https://rdrr.io/r/base/Random.html) | none |
| `M_ABCdeep` | Voting | data | No | `seed =` | none (self-contained R network) |
| `HierarchicalEnsembleClustering` | Hierarchy | data, dist | No | Deterministic | none |
| `EHC` | Hierarchy | data, dist | No | Deterministic | METIS † |
| `intNMF` | Factor | data (non-negative) | No | `seed =` | none |
| `LUCID` | Factor | data plus an outcome | No | [`set.seed()`](https://rdrr.io/r/base/Random.html) | `LUCIDus` |

Two things hold for every method except NEMO: the sources must describe
the same objects (MosaiClusteR aligns them by row name and tells you
clearly if they do not match), and a stochastic method gives a
reproducible result once you set a seed.

## Data types and distances

Most methods work on continuous matrices, but the input does not have to
be continuous.
[`Distance()`](https://lucp12891.github.io/MosaiClusteR/reference/Distance.md)
turns a source into a dissimilarity with a measure that suits the data:

- continuous: `"euclidean"`, correlation-based measures;
- binary (fingerprints, presence or absence): `"tanimoto"`, `"jaccard"`,
  `"hamming"`;
- mixed-type: `"gower"`.

Because the integration methods accept a precomputed dissimilarity
(`type = "dist"`), you can mix data types across sources: compute a
Tanimoto distance on a fingerprint source and a Euclidean distance on an
expression source, then integrate the two distance matrices with the
same function.

## A short decision path

1.  **Start with a baseline.** Run `ADC` (or `WeightedClust` if you want
    to inspect the source trade-off). It is fast and transparent, and it
    gives you a reference ARI to beat.
2.  **Some samples are missing a modality?** Use `NEMO`. It is the one
    method here that does not need every sample in every source.
3.  **Sources are noisy and complementary?** Use `SNF`. It is the
    multi-omics workhorse for genomic-scale data.
4.  **Want stability under resampling, or signal in feature subsets?**
    Use the voting family, `M_ABC` above all, and set
    `weighting = "nugget"` for large N or heavy-tailed data.
5.  **Already have several partitions to combine?** Use a graph or
    voting consensus (`EnsembleClustering`, `EvidenceAccumulation`,
    `CEC`).
6.  **The dendrogram structure is the point?** Use
    `HierarchicalEnsembleClustering`.

Then compare. Extract a partition from each result with
[`mosaic_labels()`](https://lucp12891.github.io/MosaiClusteR/reference/mosaic_labels.md),
score them with
[`cluster_agreement()`](https://lucp12891.github.io/MosaiClusteR/reference/cluster_agreement.md),
and line them up with
[`compare_clusterings()`](https://lucp12891.github.io/MosaiClusteR/reference/compare_clusterings.md)
or
[`ComparePlot()`](https://lucp12891.github.io/MosaiClusteR/reference/ComparePlot.md).
The workflow is the same whichever methods you picked.

``` r

library(MosaiClusteR)
data(mosaic_toy)
L <- mosaic_toy$List

fits <- list(
  ADC  = ADC(L, distmeasure = "euclidean", normalize = FALSE),
  SNF  = SNF(L, type = "data", distmeasure = c("euclidean", "euclidean"),
             normalize = c(FALSE, FALSE), method = list(NULL, NULL)),
  MABC = M_ABC(L, NC = 3, seed = 1)
)

labs  <- lapply(fits, mosaic_labels, k = 3)
sapply(labs, cluster_agreement, b = mosaic_toy$truth)
```
