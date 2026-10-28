# MosaiClusteR: run one method per paradigm and compare them
# https://lucp12891.github.io/MosaiClusteR/
#
# No single method is best across all data regimes, so the recommended workflow
# is to run several on the same footing and compare how well each recovers the
# structure. Here each method is cut to K clusters and scored by ARI.

library(MosaiClusteR)

data(mosaic_toy)
L     <- mosaic_toy$List
truth <- mosaic_toy$truth
S     <- length(L)
K     <- 3
DM    <- rep("euclidean", S)
NM    <- rep(FALSE, S)
ME    <- rep(list(NULL), S)

methods <- list(
  "ADC (direct)"           = function() mosaic_labels(
      ADC(L, distmeasure = "euclidean", normalize = FALSE), K),
  "WeightedClust (simil.)" = function() mosaic_labels(
      WeightedClust(L, type = "data", distmeasure = DM, normalize = NM,
                    method = ME, weight = seq(1, 0, -0.25)), K),
  "SNF (similarity)"       = function() mosaic_labels(
      SNF(L, type = "data", distmeasure = DM, normalize = NM, method = ME), K),
  "NEMO (neighbourhood)"   = function() NEMO(L, k = K, NN = 20)$cluster,
  "CEC (consensus)"        = function() mosaic_labels(
      CEC(L, distmeasure = DM, normalize = NM, method = ME,
          r = vapply(L, ncol, 1L) %/% 2L, nrclusters = as.list(rep(K, S)),
          linkage = rep("ward", S)), K),
  "intNMF (factor)"        = function() intNMF(L, k = K, nstart = 4)$cluster,
  "M-ABC (voting)"         = function() mosaic_labels(
      M_ABC(L, numsim = 60, NC = K, seed = 1), K)
)

run_one <- function(name, f) {
  t0  <- proc.time()[["elapsed"]]
  res <- tryCatch(
    data.frame(method = name, status = "ok",
               ARI = round(cluster_agreement(f(), truth)[["ARI"]], 2),
               message = "", stringsAsFactors = FALSE),
    error = function(e)
      data.frame(method = name, status = "failed", ARI = NA_real_,
                 message = conditionMessage(e), stringsAsFactors = FALSE))
  res$seconds <- round(proc.time()[["elapsed"]] - t0, 2)
  res[c("method", "status", "ARI", "seconds", "message")]
}

comparison <- do.call(rbind, Map(run_one, names(methods), methods))
comparison <- comparison[order(-comparison$ARI, na.last = TRUE), ]
rownames(comparison) <- NULL
print(comparison)
