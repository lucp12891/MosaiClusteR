# Choosing the number of clusters, and stability

How many clusters are there? MosaiClusteR separates this question from
the integration step: a method builds a dissimilarity and a hierarchy,
and you choose `k` when you extract labels. This page shows how to pick
`k` and how to check that the answer is stable.

## k is chosen at label extraction

`mosaic_labels(fit, k)` cuts any fitted result into `k` clusters. You
can try several values of `k` on one fitted result without re-running
the integration:

``` r

library(MosaiClusteR)
data(mosaic_toy)
L <- mosaic_toy$List
fit <- SNF(L, type = "data", distmeasure = rep("euclidean", length(L)),
           normalize = rep(FALSE, length(L)), method = rep(list(NULL), length(L)))

lab3 <- mosaic_labels(fit, k = 3)
lab4 <- mosaic_labels(fit, k = 4)
```

## Scoring a range of k

[`SelectnrClusters()`](https://lucp12891.github.io/MosaiClusteR/reference/SelectnrClusters.md)
scores a range of `k` with several criteria and returns one
recommendation per criterion, as a one-row data frame. Reduce it to a
single integer, for example by majority vote:

``` r

S  <- length(L)
sel <- SelectnrClusters(L, type = "data", distmeasure = rep("euclidean", S),
                        normalize = rep(FALSE, S), method = rep(list(NULL), S),
                        nrclusters = 2:6)
sel$Optimal_Nr_of_CLusters                 # one column per criterion

recs  <- as.integer(unlist(sel$Optimal_Nr_of_CLusters))
k_hat <- as.integer(names(sort(table(recs), decreasing = TRUE))[1])
k_hat
```

Treat the recommendation as a starting point, not a verdict. When
criteria disagree, inspect the candidates directly rather than trusting
the vote.

## Sweeping and plotting

[`cluster_k_sweep()`](https://lucp12891.github.io/MosaiClusteR/reference/cluster_k_sweep.md)
evaluates a fitted result across a range of `k` and returns a tidy data
frame of silhouette and within-dispersion values, with the best `k`
recorded as an attribute.
[`cluster_selection_plot()`](https://lucp12891.github.io/MosaiClusteR/reference/cluster_selection_plot.md)
draws the criterion against `k` so an elbow or a silhouette peak is easy
to see.

``` r

sweep <- cluster_k_sweep(fit, k_range = 2:8)
sweep                                  # k, silhouette, within_dispersion
cluster_selection_plot(fit, k_range = 2:8, select = "elbow")
```

How to read them: the **silhouette** rises to a peak at the `k` that
separates the objects best; the **within-dispersion** curve bends at an
elbow where adding clusters stops paying off. A `k` that both criteria
favour is a confident choice.

## Cutting a fitted tree interactively

[`ChooseCluster()`](https://lucp12891.github.io/MosaiClusteR/reference/ChooseCluster.md)
helps pick a cut height on a fitted hierarchy, which is useful when you
want to see the tree before committing to a `k`.

## Checking stability

A cluster count is more trustworthy when it survives perturbation. Two
simple checks:

- **Across methods.** Extract labels from several integration methods at
  the same `k` and compare them with
  [`cluster_agreement()`](https://lucp12891.github.io/MosaiClusteR/reference/cluster_agreement.md).
  Structure that several paradigms agree on is more likely to be real.
- **Across resamples.** Re-run a stochastic method (the M-ABC family,
  `intNMF`) with different seeds, or subsample the objects, and compare
  the resulting labels. High agreement means the solution is stable at
  that `k`.

``` r

labs <- list(
  ADC = mosaic_labels(ADC(L, distmeasure = "euclidean", normalize = FALSE), 3),
  SNF = mosaic_labels(fit, 3),
  ABC = mosaic_labels(M_ABC(L, NC = 3, seed = 1), 3)
)
# pairwise agreement between solutions at k = 3
outer(seq_along(labs), seq_along(labs),
      Vectorize(function(i, j) cluster_agreement(labs[[i]], labs[[j]])[["ARI"]]))
```

If the agreement is high across methods and across seeds, the `k` you
picked is stable. If it drops, try neighbouring values of `k` and prefer
the one whose structure reproduces.
