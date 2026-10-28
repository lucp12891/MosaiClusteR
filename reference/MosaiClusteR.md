# MosaiClusteR: An Umbrella Framework for Multi-Source and Multi-Omics Clustering

**MosaiClusteR** ("MoSaIC" = Multi-Omics Source-Agnostic Integration
Clustering in R) unifies a large family of multi-source clustering
methodologies behind a single, consistent list-of-matrices interface,
and wraps them in a complete analysis framework: preprocessing,
single-source baselines, integrative clustering across five paradigms,
and downstream comparison/evaluation.

## The MoSaIC framework (how to use the package)

1.  **Preprocess** each source with
    [`Normalization`](https://lucp12891.github.io/MosaiClusteR/reference/Normalization.md)
    and
    [`Distance`](https://lucp12891.github.io/MosaiClusteR/reference/Distance.md).

2.  **Baseline** each source on its own (any hierarchical clusterer) to
    motivate integration.

3.  **Integrate** with one or more methods, e.g.
    [`M_ABC`](https://lucp12891.github.io/MosaiClusteR/reference/M_ABC.md),
    [`M_ABCdist`](https://lucp12891.github.io/MosaiClusteR/reference/M_ABCdist.md),
    [`WeightedClust`](https://lucp12891.github.io/MosaiClusteR/reference/WeightedClust.md),
    [`SNF`](https://lucp12891.github.io/MosaiClusteR/reference/SNF.md),
    [`CEC`](https://lucp12891.github.io/MosaiClusteR/reference/CEC.md),
    [`HierarchicalEnsembleClustering`](https://lucp12891.github.io/MosaiClusteR/reference/HierarchicalEnsembleClustering.md).

4.  **Evaluate** and compare solutions with
    [`compare_clusterings`](https://lucp12891.github.io/MosaiClusteR/reference/compare_clusterings.md)
    and
    [`cluster_agreement`](https://lucp12891.github.io/MosaiClusteR/reference/cluster_agreement.md).

Every clustering method returns a list containing at least `DistM` (a
dissimilarity matrix) and `Clust` (a hierarchical-clustering object), so
results are directly comparable and composable.

## Feature weighting with data nuggets

The M-ABC family weights features for subsampling. Beyond the classic
variance (`"var"`) and coefficient-of-variation (`"cv"`) schemes,
MosaiClusteR adds `weighting = "nugget"`, a robust, Big-Data-friendly
alternative based on
[`create_data_nuggets`](https://lucp12891.github.io/MosaiClusteR/reference/create_data_nuggets.md).
See
[`nugget_feature_weights`](https://lucp12891.github.io/MosaiClusteR/reference/nugget_feature_weights.md)
and
[`ABCpp.SingleInMultiple`](https://lucp12891.github.io/MosaiClusteR/reference/ABCpp.SingleInMultiple.md).

## See also

Useful links:

- <https://github.com/lucp12891/MosaiClusteR>

- Report bugs at <https://github.com/lucp12891/MosaiClusteR/issues>

## Author

**Maintainer**: Bernard Isekah Osang'ir <b.osangir@gmail.com>
([ORCID](https://orcid.org/0000-0002-5557-3602))

Authors:

- Bernard Isekah Osang'ir <b.osangir@gmail.com>
  ([ORCID](https://orcid.org/0000-0002-5557-3602))

- Marijke Van Moerbeke

Other contributors:

- Ziv Shkedy \[contributor\]

- Surya Gupta \[contributor\]

- Jürgen Claesen \[contributor\]
