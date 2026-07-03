# =============================================================================
# M_ABCpp.R
# Multi-source M-ABC -- the canonical, C++-accelerated wrapper (Osang'ir).
# Functions: f.clustABC.MultiSource, M_ABCpp
#
# This retains the ORIGINAL design:
#   * per source we run ABCpp.SingleInMultiple and row-stack the label matrices;
#   * the consensus co-clustering counts are computed by the C++ kernel
#     count_coclustering_cpp when it is available (exactly as in the original,
#     via an exists() guard), otherwise a pure-R double loop is used;
#   * the [0,1] dissimilarity is 1 - Sij ("abc") or sqrt(1 - Sij^2) ("mabc").
#
# In the original this kernel was compiled inline with Rcpp::cppFunction() at
# source() time. As a package, the SAME C++ code lives in
# src/mcpp_sampling.cpp (compiled once at install) and is picked up here by the
# very same exists("count_coclustering_cpp") guard -- so the behaviour and the
# pure-R fallback are unchanged.
#
# Corrections carried over:
#   B6  Sij denominator = co-selection count, not nrow(res)
#   B7  'dissimilarity' argument ("abc" | "mabc") replaces the old distmeth flag
#   B8  alpha is an explicit parameter
#   B9  NC2 / extra args are not forwarded to ABCpp.SingleInMultiple
#
# NEW in 0.1.1: weighting = "nugget" (data nuggets) is passed straight through
#   to ABCpp.SingleInMultiple; nfeat / nugget_type / nugget_args are forwarded.
# =============================================================================

#' Consensus clustering from stacked ABC label matrices
#'
#' Turns a row-stacked matrix of per-iteration cluster labels
#' (\eqn{K \cdot R} rows by \eqn{N} columns) into a co-clustering dissimilarity
#' and applies \code{cluster::agnes}. \eqn{S_{ij}} is the co-clustering
#' probability over iterations selecting both objects (B6). The co-clustering
#' counts use the compiled C++ kernel \code{count_coclustering_cpp}
#' (from \code{src/mcpp_sampling.cpp}) when available, otherwise a pure-R loop.
#'
#' @param res Integer label matrix (\code{0} = object not selected).
#' @param numclust Optional number of clusters to cut the tree into.
#' @param dissimilarity \code{"mabc"} (arc-cosine, default) or \code{"abc"} (linear).
#' @param linkage AGNES agglomeration method.
#' @param alpha Flexible-linkage parameter.
#' @param mds Logical; draw an MDS of the dissimilarities.
#' @return A list with \code{DistM}, \code{Clust} and (if \code{numclust} given)
#'   \code{cut}.
#' @references \insertCite{Amaratunga2008}{MosaiClusteR}
#' @seealso \code{\link{M_ABCpp}}
#' @export
f.clustABC.MultiSource <- function(res,
                                   numclust      = NULL,
                                   dissimilarity = "mabc",
                                   linkage       = "ward",
                                   alpha         = 0.625,
                                   mds           = FALSE) {

  res <- as.matrix(res)
  storage.mode(res) <- "integer"
  N <- ncol(res)

  # co-clustering counts: prefer the C++ kernel if present, else pure-R loop.
  # (Same exists() guard as the original; the kernel now ships in src/.)
  if (exists("count_coclustering_cpp", mode = "function", inherits = TRUE)) {
    counts       <- get("count_coclustering_cpp")(res)
    co_sel       <- counts$co_sel        # times i and j were BOTH selected
    not_co_clust <- counts$not_co_clust  # ... and placed in DIFFERENT clusters
  } else {
    co_sel       <- matrix(0.0, N, N)
    not_co_clust <- matrix(0.0, N, N)
    for (i in seq_len(N)) {
      for (j in i:N) {                    # upper triangle (symmetry)
        both <- (res[, i] != 0L) & (res[, j] != 0L)
        cs   <- sum(both)
        ncc  <- sum(both & (res[, i] != res[, j]))
        co_sel[i, j]       <- cs;  co_sel[j, i]       <- cs
        not_co_clust[i, j] <- ncc; not_co_clust[j, i] <- ncc
      }
    }
  }

  # Sij = co-clustering probability with the correct denominator          [B6]
  Sij <- matrix(0.0, N, N)
  nz  <- co_sel > 0
  Sij[nz] <- (co_sel[nz] - not_co_clust[nz]) / co_sel[nz]

  # Dissimilarity                                                          [B7]
  Dij <- switch(
    dissimilarity,
    abc  = 1.0 - Sij,
    mabc = { M <- 1.0 - Sij^2; M[M < 0] <- 0; sqrt(M) },   # shape-preserving
    stop("'dissimilarity' must be 'abc' or 'mabc'")
  )
  diag(Dij) <- 0.0
  rownames(Dij) <- colnames(Dij) <- colnames(res)

  if (isTRUE(mds)) {
    pts <- stats::cmdscale(stats::as.dist(Dij), k = 2)
    plot(pts, pch = 19, xlab = "Dim 1", ylab = "Dim 2",
         main = paste0("MDS - f.clustABC.MultiSource (", dissimilarity, ")"))
  }

  Clust <- cluster::agnes(Dij, diss = TRUE, method = linkage, par.method = alpha)
  out <- list(DistM = Dij, Clust = Clust)
  if (!is.null(numclust)) {
    nc <- as.integer(numclust)
    if (nc >= 2L && nc < N) out$cut <- stats::cutree(stats::as.hclust(Clust), k = nc)
  }
  attr(out, "method") <- "M-ABC"
  out
}

#' Multi-source Aggregating Bundles of Clusters (M-ABC), C++-accelerated
#'
#' Extends the single-source ABC algorithm of
#' \insertCite{Amaratunga2008}{MosaiClusteR} to \eqn{K \ge 2} sources: runs
#' \code{\link{ABCpp.SingleInMultiple}} per source, row-stacks the per-iteration
#' label matrices, and builds a single consensus dissimilarity via
#' \code{\link{f.clustABC.MultiSource}} (C++ kernel when available).
#'
#' @section Data-nugget weighting:
#' Set \code{weighting = "nugget"} (per source, recycled) to weight features by
#' the data-nugget between-group variance instead of the classic variance; see
#' \code{\link{ABCpp.SingleInMultiple}} and \code{\link{create_data_nuggets}}.
#'
#' @param List A list of \eqn{K} data matrices over the same objects,
#'   \strong{objects (samples) in rows} and features in columns.
#' @param distmeasure,weighting,normalize,bag,numvar,linkage Per-source vectors
#'   (recycled to \eqn{K}) forwarded to \code{\link{ABCpp.SingleInMultiple}}.
#' @param gr Unused; kept for compatibility.
#' @param numsim Iterations per source.
#' @param NC Base-cluster count per source (\code{NULL} = data-driven).
#' @param dissimilarity Consensus dissimilarity, \code{"mabc"} or \code{"abc"}.
#' @param final_linkage,alpha Final clustering controls.
#' @param mds Logical; MDS of the consensus dissimilarities.
#' @param nfeat Feature-subsample size forwarded to
#'   \code{\link{ABCpp.SingleInMultiple}} (\code{NULL} = \eqn{\lfloor\sqrt{G}\rfloor}).
#' @param nugget_type,nugget_args Forwarded to the nugget weighting.
#' @return A list of class \code{"M_ABC"} with \code{DistM} and \code{Clust}.
#' @references \insertAllCited{}
#' @seealso \code{\link{ABCpp.SingleInMultiple}}, \code{\link{f.clustABC.MultiSource}},
#'   \code{\link{create_data_nuggets}}
#' @examples
#' data(mosaic_toy)
#' fit <- M_ABCpp(mosaic_toy$List, numsim = 50, NC = 3,
#'                weighting = c("var", "nugget"))
#' table(stats::cutree(stats::as.hclust(fit$Clust), 3))
#' @export
M_ABCpp <- function(List,
                    distmeasure   = rep("euclidean", length(List)),
                    weighting     = rep("var",       length(List)),
                    normalize     = rep(FALSE,       length(List)),
                    gr            = c(),
                    bag           = rep(TRUE,        length(List)),
                    numsim        = 1000,
                    numvar        = rep(100,         length(List)),
                    linkage       = rep("ward.D2",   length(List)),
                    NC            = NULL,
                    dissimilarity = "mabc",
                    final_linkage = "ward",
                    alpha         = 0.625,
                    mds           = FALSE,
                    nfeat         = NULL,
                    nugget_type   = "between",
                    nugget_args   = list()) {

  K <- length(List)

  # recycle the per-source arguments to length K (as in the original)
  distmeasure <- rep_len(distmeasure, K)
  weighting   <- rep_len(weighting,   K)
  normalize   <- rep_len(normalize,   K)
  bag         <- rep_len(bag,         K)
  numvar      <- rep_len(numvar,      K)
  linkage     <- rep_len(linkage,     K)
  # per-source nugget args: one shared list, or a length-K list of lists
  if (length(nugget_args) && all(vapply(nugget_args, is.list, logical(1))))
    nugget_args <- rep_len(nugget_args, K)
  else
    nugget_args <- rep(list(nugget_args), K)

  Results <- NULL
  for (i in seq_len(K)) {
    message("M_ABCpp: source ", i, " of ", K, " (weighting = ", weighting[i], ") ...")
    norm_i <- if (isFALSE(normalize[[i]])) NULL else normalize[[i]]
    res <- ABCpp.SingleInMultiple(          # the original single-source engine
      data        = List[[i]],
      distmeasure = distmeasure[i],
      weighting   = weighting[i],
      normalize   = norm_i,
      gr          = gr,
      bag         = bag[i],
      numsim      = numsim,
      numvar      = numvar[i],
      linkage     = linkage[i],
      NC          = NC,                     # B9: no NC2 / extra args
      nfeat       = nfeat,
      nugget_type = nugget_type,
      nugget_args = nugget_args[[i]]
    )
    res[is.na(res)] <- 0L
    Results <- rbind(Results, res)
  }

  message("M_ABCpp: building consensus dissimilarity",
          if (exists("count_coclustering_cpp", mode = "function")) " (C++ kernel)" else " (pure R)",
          " ...")
  out <- f.clustABC.MultiSource(
    res           = Results,
    dissimilarity = dissimilarity,
    linkage       = final_linkage,
    alpha         = alpha,                  # B8
    mds           = mds
  )
  class(out) <- c("M_ABC", class(out))
  out
}
