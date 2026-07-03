# =============================================================================
# nemo.R  - NEMO: NEighborhood based Multi-Omics clustering
# Rappoport & Shamir (2019), Bioinformatics. Pure R.
# Per-omic locally-scaled kNN affinities are averaged across the omics in which
# each pair of objects is measured (so partial / missing modalities need no
# imputation), then spectral clustering is applied to the integrated affinity.
# =============================================================================

# locally-scaled, kNN-sparsified affinity for one omic (objects in rows)
.nemo_affinity <- function(X, NN) {
  X  <- as.matrix(X)
  Dm <- as.matrix(stats::dist(X))
  n  <- nrow(Dm)
  NN <- max(1L, min(NN, n - 1L))
  # local scale: mean distance to the NN nearest neighbours
  sigma <- apply(Dm, 1, function(d) mean(sort(d)[2:(NN + 1L)]))
  sigma[sigma <= 0] <- .Machine$double.eps
  W <- exp(-(Dm^2) / outer(sigma, sigma))
  diag(W) <- 0
  # keep only the NN strongest neighbours per row, then symmetrise by max
  for (i in seq_len(n)) {
    thr <- sort(W[i, ], decreasing = TRUE)[NN]
    W[i, W[i, ] < thr] <- 0
  }
  W <- pmax(W, t(W))
  dimnames(W) <- list(rownames(X), rownames(X))
  W
}

#' NEMO: neighborhood-based multi-omics clustering
#'
#' NEMO \insertCite{Wang2014a}{MosaiClusteR} integrates several omics by building
#' a locally-scaled, k-nearest-neighbour affinity network for each one and
#' averaging them, for every pair of objects, over only the omics in which both
#' objects are measured. Partial / missing modalities are therefore handled
#' without imputation. The integrated affinity is clustered spectrally.
#'
#' @param List A list of data matrices over the same objects, \strong{objects in
#'   rows}. Sources may cover different subsets of objects (matched by row name)
#'   for the partial-data setting.
#' @param k Number of clusters (\code{NULL} = eigengap estimate).
#' @param NN Number of neighbours for the affinity graphs and local scaling.
#' @param normalize Logical; row-normalise each affinity to a transition matrix
#'   before integration (NEMO default \code{TRUE}).
#' @return A list of class \code{"NEMO"} with \code{FusedM} (integrated
#'   affinity), \code{DistM} (\code{1 - FusedM}), \code{cluster} (labels) and
#'   \code{Clust} (an \code{agnes} object on \code{DistM} for interoperability).
#' @references \insertCite{Wang2014a}{MosaiClusteR}; Rappoport, N. & Shamir, R.
#'   (2019). NEMO: cancer subtyping by integration of partial multi-omic data.
#'   Bioinformatics 35:3348-3356.
#' @seealso \code{\link{SNF}}, \code{\link{spectral_clustering}}
#' @examples
#' data(mosaic_toy)
#' fit <- NEMO(mosaic_toy$List, k = 3, NN = 15)
#' table(fit$cluster)
#' @export
NEMO <- function(List, k = NULL, NN = 20, normalize = TRUE) {
  if (!is.list(List)) stop("'List' must be a list of data matrices.")
  objs <- Reduce(union, lapply(List, rownames))
  if (is.null(objs))
    stop("NEMO needs row names (objects) on the input matrices.")
  n <- length(objs)

  acc <- matrix(0, n, n, dimnames = list(objs, objs))
  cnt <- matrix(0, n, n, dimnames = list(objs, objs))
  for (X in List) {
    W <- .nemo_affinity(X, NN)
    if (isTRUE(normalize)) {
      rs <- rowSums(W); rs[rs <= 0] <- .Machine$double.eps
      W  <- W / rs
      W  <- (W + t(W)) / 2
    }
    idx <- match(rownames(W), objs)
    acc[idx, idx] <- acc[idx, idx] + W
    cnt[idx, idx] <- cnt[idx, idx] + 1
  }
  cnt[cnt == 0] <- NA
  Fused <- acc / cnt                       # average over omics measuring each pair
  Fused[is.na(Fused)] <- 0
  Fused <- (Fused + t(Fused)) / 2
  diag(Fused) <- 0

  sc <- spectral_clustering(Fused, k = k)

  D <- 1 - Fused / max(Fused[Fused > 0], 1)
  diag(D) <- 0
  Clust <- cluster::agnes(D, diss = TRUE, method = "ward")

  out <- list(FusedM = Fused, DistM = D, cluster = sc$cluster,
              k = sc$k, Clust = Clust)
  attr(out, "method") <- "NEMO"
  class(out) <- c("NEMO", "list")
  out
}
