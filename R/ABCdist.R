# =============================================================================
# ABCdist.R
# ABC variant using direct distance-matrix accumulation instead of label
# co-clustering. Steps 1-3 follow Amaratunga et al. (2008); steps 4-5 accumulate
# the mean scaled pairwise distance over co-selected objects.
# Main function: ABCdist.SingleInMultiple
# =============================================================================

#' Single-source ABC with distance accumulation
#'
#' A variant of \code{\link{ABCpp.SingleInMultiple}} that accumulates the mean
#' scaled pairwise distance between co-selected objects across iterations,
#' yielding a \code{[0, 1]} dissimilarity matrix that \code{\link{M_ABCdist}}
#' fuses across sources.
#'
#' @inheritParams ABCpp.SingleInMultiple
#' @return A list with \code{D} (\eqn{N \times N} dissimilarity), \code{NC} and
#'   \code{ng} (per-iteration feature count); if \code{mds = TRUE} also
#'   \code{mds_coords}.
#' @references \insertCite{Amaratunga2008}{MosaiClusteR}
#' @seealso \code{\link{M_ABCdist}}, \code{\link{ABCpp.SingleInMultiple}}
#' @examples
#' data(mosaic_toy)
#' out <- ABCdist.SingleInMultiple(mosaic_toy$List[[1]], numsim = 40)
#' dim(out$D)
#' @export
ABCdist.SingleInMultiple <- function(data,
                                     distmeasure = "euclidean",
                                     weighting   = "var",
                                     normalize   = NULL,
                                     gr          = c(),
                                     bag         = TRUE,
                                     numsim      = 1000,
                                     numvar      = 100,
                                     linkage     = "ward.D",
                                     NC          = NULL,
                                     mds         = FALSE,
                                     nfeat       = NULL,
                                     nugget_type = "between",
                                     nugget_args = list()) {
  data <- as.matrix(data)
  samp_names <- rownames(data)
  X  <- t(data)        # internal: features x objects
  ns <- ncol(X)
  nf <- nrow(X)

  if (is.null(NC)) NC <- max(2L, floor(sqrt(ns)))
  NC <- max(2L, min(as.integer(NC), ns - 1L))
  g  <- .abc_gsize(nf, nfeat)

  zf <- .abc_zf(X, weighting, nugget_type, nugget_args)

  D_sum <- matrix(0.0, ns, ns)
  D_cnt <- matrix(0L,  ns, ns)

  for (i in seq_len(numsim)) {
    if (i %% 100 == 0) message("ABCdist iteration ", i, " / ", numsim)
    ids <- if (bag) sort(unique(sample.int(ns, size = ns, replace = TRUE)))
           else      seq_len(ns)
    if (length(ids) < 2L) next
    idf <- sort(sample.int(nf, size = g, replace = FALSE, prob = zf))
    dat <- X[idf, ids, drop = FALSE]

    d_mat <- as.matrix(Distance(t(dat), distmeasure = distmeasure,
                                normalize = normalize))
    # robustness: Gower distance is NaN when a sampled feature is constant across
    # the selected objects -> fall back to Euclidean for this iteration.
    if (anyNA(d_mat) || any(!is.finite(d_mat))) d_mat <- as.matrix(stats::dist(t(dat)))
    lt   <- which(lower.tri(d_mat), arr.ind = TRUE)
    ii   <- ids[lt[, 1L]]; jj <- ids[lt[, 2L]]; vals <- d_mat[lt]
    D_sum[cbind(ii, jj)] <- D_sum[cbind(ii, jj)] + vals
    D_sum[cbind(jj, ii)] <- D_sum[cbind(jj, ii)] + vals
    D_cnt[cbind(ii, jj)] <- D_cnt[cbind(ii, jj)] + 1L
    D_cnt[cbind(jj, ii)] <- D_cnt[cbind(jj, ii)] + 1L
  }

  D_mean <- ifelse(D_cnt > 0L, D_sum / D_cnt, NA_real_)
  max_d  <- max(D_mean, na.rm = TRUE)
  if (!is.finite(max_d) || max_d <= 0) max_d <- 1
  D_mean[is.na(D_mean)] <- max_d
  D_mat <- D_mean / max_d
  diag(D_mat) <- 0.0
  rownames(D_mat) <- colnames(D_mat) <- samp_names

  out <- list(D = D_mat, NC = NC, ng = g)
  if (isTRUE(mds)) {
    fit <- stats::cmdscale(stats::as.dist(D_mat), k = 2)
    colnames(fit) <- c("Dim1", "Dim2"); out$mds_coords <- fit
  }
  out
}
