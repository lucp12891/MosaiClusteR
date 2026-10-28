## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", warning = FALSE,
                      message = FALSE, fig.width = 7, fig.height = 4.5)
set.seed(1)


## ----setup--------------------------------------------------------------------
library(MosaiClusteR)


## ----data---------------------------------------------------------------------
data(mosaic_toy)
L     <- mosaic_toy$List      # a list of two matrices, objects in rows
truth <- mosaic_toy$truth     # the known group of each object (for illustration)
S     <- length(L)            # number of sources
K     <- 3                    # number of groups
vapply(L, dim, integer(2))    # 60 objects x (120, 150) features


## ----sim, eval = FALSE--------------------------------------------------------
# sim <- mosaic_sim(n = 90, g = 3, seed = 1)   # n objects, g groups
# str(sim$List, max.level = 1)


## ----preprocess---------------------------------------------------------------
Xn <- Normalization(L[[1]], method = "Quantile")   # optional per-source scaling
D1 <- Distance(L[[1]], distmeasure = "euclidean")   # object-by-object distance
round(as.matrix(D1)[1:4, 1:4], 2)


## ----baseline-----------------------------------------------------------------
b1 <- Cluster(L[[1]], distmeasure = "euclidean", normalize = FALSE, nrclusters = 3)
b2 <- Cluster(L[[2]], distmeasure = "euclidean", normalize = FALSE, nrclusters = 3)
c(source1 = cluster_agreement(mosaic_labels(b1, 3), truth)[["ARI"]],
  source2 = cluster_agreement(mosaic_labels(b2, 3), truth)[["ARI"]])


## ----integrate----------------------------------------------------------------
DM <- rep("euclidean", S); NM <- rep(FALSE, S); ME <- rep(list(NULL), S)

methods <- list(
  "ADC (direct)"           = function() mosaic_labels(
      ADC(L, distmeasure = "euclidean", normalize = FALSE), K),
  "WeightedClust (simil.)" = function() mosaic_labels(
      WeightedClust(L, type = "data", distmeasure = DM, normalize = NM,
                    method = ME, weight = seq(1, 0, -0.25)), K),
  "WonM (similarity)"      = function() mosaic_labels(
      WonM(L, type = "data", distmeasure = DM, normalize = NM, method = ME,
           nrclusters = rep(list(2:6), S), linkage = rep("ward", S)), K),
  "SNF (graph fusion)"     = function() mosaic_labels(
      SNF(L, type = "data", distmeasure = DM, normalize = NM, method = ME,
          NN = 20, T = 20), K),
  "NEMO (neighbourhood)"   = function() NEMO(L, k = K, NN = 20)$cluster,
  "CEC (consensus)"        = function() mosaic_labels(
      CEC(L, distmeasure = DM, normalize = NM, method = ME,
          r = vapply(L, ncol, 1L) %/% 2L, nrclusters = as.list(rep(K, S)),
          linkage = rep("ward", S)), K),
  "intNMF (factor)"        = function() intNMF(L, k = K, nstart = 4)$cluster,
  "Data nuggets"           = function() nugget_cluster(
      do.call(cbind, L), k = K, max_nuggets = 40)$cluster,
  "ABC (voting)"           = function() mosaic_labels(
      M_ABC(L, numsim = 60, NC = K), K)
)

ari <- vapply(methods, function(f)
  tryCatch(round(cluster_agreement(f(), truth)[["ARI"]], 2),
           error = function(e) NA_real_), numeric(1))
sort(ari, decreasing = TRUE)


## ----direct, eval = FALSE-----------------------------------------------------
# ADC(L, distmeasure = "euclidean", normalize = FALSE)
# ADEC(L, distmeasure = DM, normalize = NM, method = ME, nrclusters = 3, t = 10)


## ----simil, eval = FALSE------------------------------------------------------
# WeightedClust(L, type = "data", distmeasure = DM, normalize = NM, method = ME,
#               weight = seq(1, 0, -0.25))
# WonM(L, type = "data", distmeasure = DM, normalize = NM, method = ME,
#      nrclusters = rep(list(2:6), S), linkage = rep("ward", S))
# SNF(L, type = "data", distmeasure = DM, normalize = NM, method = ME, NN = 20, T = 20)
# NEMO(L, k = K, NN = 20)$cluster


## ----graph, eval = FALSE------------------------------------------------------
# HBGF(L, type = "data", distmeasure = DM, normalize = NM, method = ME, nrclusters = 3)
# ClusteringAggregation(L, type = "data", distmeasure = DM, normalize = NM,
#                       method = ME, nrclusters = 3)
# # EnsembleClustering() / EHC() call an external graph-partitioning executable
# # (e.g. METIS); install it and pass its path via 'executable' to use them.


## ----voting, eval = FALSE-----------------------------------------------------
# CEC(L, distmeasure = DM, normalize = NM, method = ME,
#     r = vapply(L, ncol, 1L) %/% 2L, nrclusters = as.list(rep(K, S)),
#     linkage = rep("ward", S))
# ConsensusClustering(L, type = "data", distmeasure = DM, normalize = NM,
#                     method = ME, nrclusters = 3)
# EvidenceAccumulation(L, type = "data", distmeasure = DM, normalize = NM,
#                      method = ME, nrclusters = 3)
# LinkBasedClustering(L, type = "data", distmeasure = DM, normalize = NM,
#                     method = ME, nrclusters = 3, linkBasedMethod = "SRS")
# M_ABC(L, numsim = 200, NC = 3)                 # co-clustering consensus
# M_ABCdist(L, numsim = 200, NC = 3)               # distance-accumulating variant
# M_ABCdeep(L, numsim = 200, NC = 3, ae_epochs = 60)   # autoencoder embedding


## ----hier, eval = FALSE-------------------------------------------------------
# HierarchicalEnsembleClustering(L, type = "data", distmeasure = DM,
#                                normalize = NM, method = ME)


## ----factor, eval = FALSE-----------------------------------------------------
# intNMF(L, k = K, nstart = 8)$cluster
# spectral_clustering(affinity = as.matrix(exp(-as.matrix(D1))), k = K)
# LUCID(G = exposures, Z = L, Y = outcome, K = 3)   # supervised; needs 'LUCIDus'


## ----nuggets------------------------------------------------------------------
dn <- create_data_nuggets(L[[1]], max_nuggets = 25)   # compress objects
dn
w  <- nugget_feature_weights(dn, type = "between")    # feature importance
head(sort(w, decreasing = TRUE))


## ----nuggets2, eval = FALSE---------------------------------------------------
# nugget_cluster(do.call(cbind, L), k = K, max_nuggets = 40)$cluster  # cluster nuggets
# Wkmeans(dn, k = K); Whclust(dn, k = K)                              # weighted k-means / hclust
# M_ABC(L, numsim = 200, NC = 3, weighting = c("nugget", "nugget")) # nugget-weighted


## ----selectk------------------------------------------------------------------
sel <- SelectnrClusters(L, type = "data", distmeasure = DM, normalize = NM,
                        method = ME, nrclusters = 2:6)
sel$Optimal_Nr_of_CLusters


## ----evaluate-----------------------------------------------------------------
fit_a <- ADC(L, distmeasure = "euclidean", normalize = FALSE)
fit_b <- SNF(L, type = "data", distmeasure = DM, normalize = NM, method = ME)
compare_clusterings(fit_a$DistM, fit_b$DistM, NC = 3)$ARI
cluster_agreement(mosaic_labels(fit_b, 3), truth)


## ----dendro, fig.height = 4---------------------------------------------------
fit <- SNF(L, type = "data", distmeasure = DM, normalize = NM, method = ME)
hc  <- stats::as.hclust(fit$Clust)
plot(hc, labels = FALSE, main = "SNF consensus", xlab = "", sub = "")
stats::rect.hclust(hc, k = 3, border = 2:4)


## ----viz, eval = FALSE--------------------------------------------------------
# ComparePlot(list("ADC" = fit_a, "SNF" = fit_b), nrclusters = 3)
# SimilarityHeatmap(fit_b$DistM)
# ProfilePlot(L[[1]], mosaic_labels(fit_b, 3))


## ----charac, eval = FALSE-----------------------------------------------------
# CharacteristicFeatures(L, mosaic_labels(fit_b, 3))
# ReorderToReference(reference = truth, tocompare = mosaic_labels(fit_a, 3))


## ----bio, eval = FALSE--------------------------------------------------------
# dg <- DiffGenes(L, mosaic_labels(fit_b, 3), geneexpr = L[[1]])
# pa <- PathwayAnalysis(L, mosaic_labels(fit_b, 3), geneexpr = L[[1]])
# PlotPathways(pa)


## ----together, eval = FALSE---------------------------------------------------
# L    <- my_sources                                   # list of matrices, objects in rows
# k    <- SelectnrClusters(L, type = "data", distmeasure = rep("euclidean", length(L)),
#                          normalize = rep(FALSE, length(L)),
#                          method = rep(list(NULL), length(L)),
#                          nrclusters = 2:8)$Optimal_Nr_of_CLusters
# fit  <- SNF(L, type = "data", distmeasure = rep("euclidean", length(L)),
#             normalize = rep(FALSE, length(L)), method = rep(list(NULL), length(L)))
# lab  <- mosaic_labels(fit, k)
# CharacteristicFeatures(L, lab)                        # what defines the clusters


## ----session, echo = FALSE----------------------------------------------------
sessionInfo()

