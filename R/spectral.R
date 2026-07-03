# =============================================================================
# spectral.R  - spectral clustering on a similarity/affinity matrix
# Used directly and as the final clustering step of NEMO and similarity fusion.
# =============================================================================

#' Spectral clustering of an affinity matrix
#'
#' Normalised (Ng-Jordan-Weiss) spectral clustering. Given a symmetric,
#' non-negative affinity matrix it forms the normalised Laplacian, embeds the
#' objects in the space of its leading eigenvectors and applies k-means.
#'
#' @param affinity A symmetric \eqn{n \times n} non-negative similarity/affinity
#'   matrix over the objects.
#' @param k Number of clusters. If \code{NULL}, chosen by the largest eigengap
#'   of the normalised Laplacian (searched over \code{2:max_k}).
#' @param max_k Upper bound for the eigengap search when \code{k = NULL}.
#' @param nstart Number of k-means restarts.
#' @return A list with \code{cluster} (named integer labels), \code{k} and the
#'   spectral \code{embedding}.
#' @references Ng, A. Y., Jordan, M. I. & Weiss, Y. (2002). On spectral
#'   clustering: analysis and an algorithm. NIPS 14.
#' @seealso \code{\link{NEMO}}, \code{\link{SNF}}
#' @examples
#' set.seed(1)
#' X <- rbind(matrix(rnorm(20 * 5, 0), 20), matrix(rnorm(20 * 5, 6), 20))
#' A <- exp(-as.matrix(dist(X))^2 / 2)
#' spectral_clustering(A, k = 2)$cluster[1:5]
#' @export
spectral_clustering <- function(affinity, k = NULL, max_k = 10, nstart = 10) {
  W <- as.matrix(affinity)
  W <- (W + t(W)) / 2
  diag(W) <- 0
  n <- nrow(W)
  d <- rowSums(W)
  d[d <= 0] <- .Machine$double.eps
  Dm <- 1 / sqrt(d)
  Lsym <- diag(n) - (Dm * W) * rep(Dm, each = n)   # I - D^-1/2 W D^-1/2
  Lsym <- (Lsym + t(Lsym)) / 2

  ev <- eigen(Lsym, symmetric = TRUE)
  vals <- sort(ev$values)                            # ascending

  if (is.null(k)) {
    upper <- min(max_k, n - 1L)
    gaps  <- diff(vals[seq_len(upper + 1L)])
    k <- which.max(gaps[-1]) + 1L                    # ignore trivial first gap
    k <- max(2L, k)
  }
  k <- max(2L, min(as.integer(k), n - 1L))

  # k smallest eigenvectors of Lsym (largest of the normalised affinity)
  U <- ev$vectors[, order(ev$values)[seq_len(k)], drop = FALSE]
  rn <- sqrt(rowSums(U^2)); rn[rn == 0] <- 1
  U  <- U / rn                                        # row-normalise

  km <- stats::kmeans(U, centers = k, nstart = nstart, iter.max = 50L)
  cl <- stats::setNames(as.integer(km$cluster), rownames(W))
  list(cluster = cl, k = k, embedding = U)
}
