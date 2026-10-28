# Evaluation, with and without ground truth

Once you have a clustering, you need to judge it. The tools differ
depending on whether you have a known grouping to compare against. This
page covers both cases, and how to compare methods against each other.

``` r

library(MosaiClusteR)
data(mosaic_toy)
L     <- mosaic_toy$List
truth <- mosaic_toy$truth
fit   <- SNF(L, type = "data", distmeasure = rep("euclidean", length(L)),
             normalize = rep(FALSE, length(L)), method = rep(list(NULL), length(L)))
lab   <- mosaic_labels(fit, 3)
```

## With ground truth

When you know the true grouping (a benchmark, a simulation, a validated
subtype),
[`cluster_agreement()`](https://lucp12891.github.io/MosaiClusteR/reference/cluster_agreement.md)
scores a clustering against it with four complementary measures:

``` r

cluster_agreement(lab, truth)   # ARI, NMI, Jaccard, purity
```

- **ARI** (Adjusted Rand Index): agreement on which pairs share a
  cluster, corrected for chance. The headline number.
- **NMI** (Normalised Mutual Information): shared information between
  the two labellings.
- **Jaccard**: overlap of co-membership, uncorrected for chance.
- **Purity**: how cleanly each cluster maps to a single true group.

Rule of thumb for ARI: above 0.90 excellent, 0.75 to 0.90 good, 0.50 to
0.75 moderate, below 0.50 weak. All of ARI, NMI and Jaccard return 1 for
identical partitions up to relabelling, including the single-cluster and
all-singleton cases, and missing labels raise a clear error rather than
a cryptic one.

Labels are arbitrary names, so before comparing two solutions by eye,
align one to the other with
[`ReorderToReference()`](https://lucp12891.github.io/MosaiClusteR/reference/ReorderToReference.md):

``` r

ReorderToReference(reference = truth, tocompare = lab)
```

## Without ground truth

On real data you usually have no truth to compare against. Judge a
clustering by how well separated it is and by how reliably it
reproduces.

**Internal separation.** The silhouette measures how much better each
object fits its own cluster than the nearest other cluster.
[`CompareSilCluster()`](https://lucp12891.github.io/MosaiClusteR/reference/CompareSilCluster.md)
compares candidate solutions by silhouette, and
[`cluster_k_sweep()`](https://lucp12891.github.io/MosaiClusteR/reference/cluster_k_sweep.md)
reports silhouette and within-dispersion across `k`.

``` r

CompareSilCluster(fit$DistM, nrclusters = 2:6)
```

**Single source versus multi-source.**
[`CompareSvsM()`](https://lucp12891.github.io/MosaiClusteR/reference/CompareSvsM.md)
contrasts a single-source clustering with the integrated one, which
shows whether integration actually added structure over the best single
source.

**Stability.** A solution you cannot reproduce is not worth
interpreting. Re-run a stochastic method with different seeds, or
subsample the objects, and check that the labels agree with
[`cluster_agreement()`](https://lucp12891.github.io/MosaiClusteR/reference/cluster_agreement.md).
The stability checks in the choosing-k guide apply directly here.

## Comparing two integrations

[`compare_clusterings()`](https://lucp12891.github.io/MosaiClusteR/reference/compare_clusterings.md)
clusters two dissimilarities at the same `k` and reports their
agreement, so you can ask how much two methods actually differ:

``` r

fit_a <- ADC(L, distmeasure = "euclidean", normalize = FALSE)
compare_clusterings(fit_a$DistM, fit$DistM, NC = 3)   # ARI, NMI, Jaccard + crosstab
```

[`clustering_concordance()`](https://lucp12891.github.io/MosaiClusteR/reference/clustering_concordance.md)
summarises agreement across several solutions at once, and
[`compare_methods_plot()`](https://lucp12891.github.io/MosaiClusteR/reference/compare_methods_plot.md)
and
[`ComparePlot()`](https://lucp12891.github.io/MosaiClusteR/reference/ComparePlot.md)
show the solutions side by side so disagreements are visible, not just
numeric.

## A comparison-first workflow

Because every method shares one interface, the recommended practice is
to run several, score them the same way, and prefer structure that
reproduces:

``` r

methods <- list(
  ADC = function() mosaic_labels(ADC(L, distmeasure = "euclidean", normalize = FALSE), 3),
  SNF = function() mosaic_labels(fit, 3),
  ABC = function() mosaic_labels(M_ABC(L, NC = 3, seed = 1), 3)
)
sapply(methods, function(f) cluster_agreement(f(), truth)[["ARI"]])
```

With a known truth, read the ARI column. Without one, lean on silhouette
and on agreement across methods and seeds.
