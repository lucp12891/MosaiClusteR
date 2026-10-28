# MosaiClusteR: a complete end-to-end analysis
# https://lucp12891.github.io/MosaiClusteR/
#
# Four layers: preprocess -> single-source baselines -> integrate -> evaluate,
# then characterise the clusters. Uses the bundled toy data; swap in your own
# list of matrices (objects in rows) to run it on real data.

library(MosaiClusteR)

data(mosaic_toy)
L     <- mosaic_toy$List
truth <- mosaic_toy$truth
S     <- length(L)
DM    <- rep("euclidean", S)
NM    <- rep(FALSE, S)
ME    <- rep(list(NULL), S)

## Layer 1: preprocess -------------------------------------------------------
Xn <- Normalization(L[[1]], method = "Quantile")   # make a source comparable
D1 <- Distance(Xn, distmeasure = "euclidean")       # distances from normalised data

## Layer 2: single-source baselines -----------------------------------------
b1 <- Cluster(L[[1]], distmeasure = "euclidean", normalize = FALSE)
b2 <- Cluster(L[[2]], distmeasure = "euclidean", normalize = FALSE)
c(source1 = cluster_agreement(mosaic_labels(b1, 3), truth)[["ARI"]],
  source2 = cluster_agreement(mosaic_labels(b2, 3), truth)[["ARI"]])

## Choose the number of clusters --------------------------------------------
sel   <- SelectnrClusters(L, type = "data", distmeasure = DM, normalize = NM,
                          method = ME, nrclusters = 2:6)
recs  <- as.integer(unlist(sel$Optimal_Nr_of_CLusters))
k_hat <- as.integer(names(sort(table(recs), decreasing = TRUE))[1])

## Layer 3: integrate --------------------------------------------------------
fit <- SNF(L, type = "data", distmeasure = DM, normalize = NM, method = ME)
lab <- mosaic_labels(fit, k_hat)

## Layer 4: evaluate ---------------------------------------------------------
cluster_agreement(lab, truth)                       # with ground truth
fit_a <- ADC(L, distmeasure = "euclidean", normalize = FALSE)
compare_clusterings(fit_a$DistM, fit$DistM, NC = k_hat)   # two integrations

## Characterise the clusters -------------------------------------------------
CharacteristicFeatures(L, lab)                      # what defines each cluster
