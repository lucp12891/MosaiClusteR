# =============================================================================
# intnmf.R  - integrative Non-negative Matrix Factorisation (intNMF)
# Concatenated low-rank integration: a shared basis W (objects x r) and per-omic
# coefficient matrices H_k are found so that X_k ~ W H_k for all sources, via
# multiplicative updates. Objects are clustered from W by consensus over
# restarts. (Chalise & Fridley 2017; review: Zhang et al. 2022.)
# =============================================================================

#' Integrative NMF clustering (intNMF)
#'
#' Joint non-negative matrix factorisation across data sources: a single shared
#' basis matrix \eqn{W} (objects \eqn{\times} rank) and source-specific
#' coefficient matrices \eqn{H_k} minimise \eqn{\sum_k \alpha_k \lVert X_k - W
#' H_k \rVert_F^2} with all factors non-negative. Sources are weighted by
#' \eqn{\alpha_k = \max_j \overline{SS}_j / \overline{SS}_k}. Object clusters are
#' obtained by consensus of the dominant-factor assignment across restarts.
#'
#' @param List A list of data matrices over the same objects, \strong{objects in
#'   rows}. Matrices are shifted to be non-negative internally.
#' @param k Number of clusters (also the factorisation rank).
#' @param nstart Number of random restarts forming the consensus.
#' @param max_iter Maximum multiplicative-update iterations per restart.
#' @param seed Optional RNG seed.
#' @return A list of class \code{"intNMF"} with \code{cluster}, the consensus
#'   \code{DistM}, an \code{agnes} \code{Clust}, the best-fit basis \code{W} and
#'   source weights \code{alpha}.
#' @references Chalise, P. & Fridley, B. L. (2017). Integrative clustering of
#'   multi-level omic data based on non-negative matrix factorization. PLoS ONE
#'   12:e0176278. Reviewed in Zhang et al. (2022), WIREs Comput Stat 14:e1553.
#' @seealso \code{\link{M_ABCpp}}, \code{\link{SNF}}
#' @examples
#' data(mosaic_toy)
#' fit <- intNMF(mosaic_toy$List, k = 3, nstart = 5, seed = 1)
#' table(fit$cluster)
#' @export
intNMF <- function(List, k, nstart = 10, max_iter = 200, seed = NULL) {
  if (!is.list(List)) stop("'List' must be a list of data matrices.")
  if (!is.null(seed)) set.seed(seed)
  Xs <- lapply(List, function(X) {
    X <- as.matrix(X); X <- X - min(X, na.rm = TRUE); X[is.na(X)] <- 0; X })
  n  <- nrow(Xs[[1]])
  r  <- max(2L, as.integer(k))
  eps <- .Machine$double.eps

  ssm   <- vapply(Xs, function(X) mean(X^2), numeric(1))
  alpha <- max(ssm) / pmax(ssm, eps)

  fit_one <- function() {
    W  <- matrix(stats::runif(n * r), n, r)
    Hs <- lapply(Xs, function(X) matrix(stats::runif(r * ncol(X)), r, ncol(X)))
    for (it in seq_len(max_iter)) {
      for (kk in seq_along(Xs)) {
        WtW <- crossprod(W)
        Hs[[kk]] <- Hs[[kk]] * (crossprod(W, Xs[[kk]]) /
                                  pmax(WtW %*% Hs[[kk]], eps))
      }
      num <- matrix(0, n, r); den <- matrix(0, n, r)
      for (kk in seq_along(Xs)) {
        num <- num + alpha[kk] * tcrossprod(Xs[[kk]], Hs[[kk]])
        den <- den + alpha[kk] * W %*% tcrossprod(Hs[[kk]])
      }
      W <- W * (num / pmax(den, eps))
    }
    err <- sum(vapply(seq_along(Xs), function(kk)
      alpha[kk] * sum((Xs[[kk]] - W %*% Hs[[kk]])^2), numeric(1)))
    list(W = W, err = err, lab = max.col(W, ties.method = "first"))
  }

  CO <- matrix(0, n, n); best <- NULL; best_err <- Inf
  for (s in seq_len(nstart)) {
    f <- fit_one()
    CO <- CO + outer(f$lab, f$lab, "==")
    if (f$err < best_err) { best_err <- f$err; best <- f }
  }
  CO <- CO / nstart
  D  <- 1 - CO; diag(D) <- 0
  rownames(D) <- colnames(D) <- rownames(Xs[[1]])
  Clust <- cluster::agnes(D, diss = TRUE, method = "ward")
  cl <- stats::setNames(stats::cutree(stats::as.hclust(Clust), k = r),
                        rownames(Xs[[1]]))

  out <- list(cluster = cl, DistM = D, Clust = Clust,
              W = best$W, alpha = alpha, k = r)
  attr(out, "method") <- "intNMF"
  class(out) <- c("intNMF", "list")
  out
}
