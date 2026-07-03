# =============================================================================
# nuggets.R
# Data nuggets: a representative-sample compression of a data set, and the
# nugget-based feature-weighting scheme used as an alternative to variance
# weighting inside the ABC / M-ABC family.
#
# A data nugget summarises a subset of observations by three quantities
# (Cherasia et al. 2023; Beavers/Cabrera, WCluster):
#   * center : the within-nugget centroid (location)
#   * weight : the number of observations assigned to the nugget (importance)
#   * scale  : the within-nugget variability (spread / internal variance)
#
# In MosaiClusteR data nuggets are used in two ways:
#   1. As a Big-Data compression that any downstream clusterer can run on.
#   2. To derive *feature* weights for M-ABC: a feature that varies a lot
#      *between* nugget centres (weighted by nugget size) carries signal,
#      whereas a feature whose variance is mostly *within* nuggets is noise.
#      The nugget between-group weighted variance therefore replaces the plain
#      feature variance as the basis for the ABC selection probabilities.
# =============================================================================

#' Create data nuggets from a data matrix
#'
#' \code{create_data_nuggets} compresses \code{n} observations into
#' \code{M <= n} representative "nuggets". Each nugget stores a centre
#' (location), a weight (the number of member observations) and a scale
#' (within-nugget variability). Data nuggets preserve the covariance structure
#' of the data - including the periphery - far better than a random sub-sample,
#' which makes them a robust, memory-light surrogate for the full data set when
#' clustering or computing feature weights.
#'
#' @param x A numeric matrix or data frame with observations in the \emph{rows}
#'   and variables in the \emph{columns}.
#' @param max_nuggets Target (maximum) number of nuggets \eqn{M}. Defaults to a
#'   data-driven value (see Details). Ignored when \code{engine = "datanugget"}.
#' @param min_size Minimum number of observations per nugget (native engine
#'   only). Defaults to \code{1}.
#' @param engine Either \code{"native"} (a self-contained k-means based
#'   partition, no extra dependencies) or \code{"datanugget"} (delegates to the
#'   \pkg{datanugget} package via \code{datanugget::create.DN} /
#'   \code{datanugget::refine.DN} when it is installed).
#' @param scale_features Logical; if \code{TRUE} (default) the partition is
#'   computed on standardised features so that no single high-variance variable
#'   dominates the nugget assignment. Centres and scales are always reported on
#'   the original scale.
#' @param seed Optional integer seed for reproducible partitioning.
#' @param ... Additional arguments forwarded to
#'   \code{datanugget::create.DN} when \code{engine = "datanugget"}.
#'
#' @details
#' With the native engine the number of nuggets defaults to
#' \code{min(nrow(x), max(2, floor(nrow(x) / 10)))} when \code{n > 50}, and to
#' \code{n} otherwise (in which case every observation is its own nugget and the
#' nugget-weighted variance reduces to the ordinary variance - a safe fallback
#' for small samples). The partition is obtained with \code{stats::kmeans};
#' missing values are mean-imputed for the partition step only.
#'
#' @return An object of class \code{"data_nugget"}: a list with elements
#'   \item{centers}{\eqn{M \times p} matrix of nugget centres.}
#'   \item{weights}{Length-\eqn{M} integer vector of nugget sizes \eqn{n_i}.}
#'   \item{scales}{\eqn{M \times p} matrix of within-nugget variances.}
#'   \item{membership}{Length-\eqn{n} integer vector mapping each observation to
#'         its nugget.}
#'   \item{x_scale}{Length-\eqn{M} vector of total within-nugget variability
#'         \eqn{\mathrm{tr}(\mathrm{Cov})}.}
#'   \item{engine}{The engine used.}
#'
#' @references
#' \insertRef{Cherasia2023}{MosaiClusteR}
#'
#' \insertRef{WCluster2025}{MosaiClusteR}
#'
#' @seealso \code{\link{nugget_feature_weights}}, \code{\link{M_ABCpp}}
#' @examples
#' set.seed(1)
#' x <- rbind(matrix(rnorm(200 * 5,  0), ncol = 5),
#'            matrix(rnorm(200 * 5,  4), ncol = 5))
#' dn <- create_data_nuggets(x, max_nuggets = 40)
#' dn
#' @export
create_data_nuggets <- function(x,
                                 max_nuggets   = NULL,
                                 min_size      = 1L,
                                 engine        = c("native", "datanugget"),
                                 scale_features = TRUE,
                                 seed          = NULL,
                                 ...) {
  engine <- match.arg(engine)
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  n <- nrow(x)
  p <- ncol(x)
  if (n < 2L) stop("create_data_nuggets: need at least 2 observations (rows).")
  if (!is.null(seed)) set.seed(seed)

  # ---- delegate to the datanugget package when requested & available --------
  if (engine == "datanugget") {
    if (!requireNamespace("datanugget", quietly = TRUE))
      stop("engine = 'datanugget' requires the 'datanugget' package; ",
           "install it or use engine = 'native'.")
    return(.dn_from_datanugget(x, ...))
  }

  # ---- native k-means partition --------------------------------------------
  if (is.null(max_nuggets)) {
    max_nuggets <- if (n > 50L) min(n, max(2L, floor(n / 10))) else n
  }
  M <- max(2L, min(as.integer(max_nuggets), n))

  xp <- x
  if (anyNA(xp)) {
    cm <- colMeans(xp, na.rm = TRUE); cm[is.na(cm)] <- 0
    for (j in seq_len(p)) xp[is.na(xp[, j]), j] <- cm[j]
  }
  xs <- if (isTRUE(scale_features)) scale(xp) else xp
  xs[is.na(xs)] <- 0  # zero-variance columns -> NaN after scale()

  if (M >= n) {
    membership <- seq_len(n)
  } else {
    km <- stats::kmeans(xs, centers = M, iter.max = 30L, nstart = 1L)
    membership <- km$cluster
    # enforce min_size by absorbing tiny nuggets into their nearest neighbour
    if (min_size > 1L) membership <- .absorb_small(xs, membership, min_size)
  }

  ids     <- sort(unique(membership))
  M       <- length(ids)
  centers <- matrix(0, M, p, dimnames = list(NULL, colnames(x)))
  scales  <- matrix(0, M, p, dimnames = list(NULL, colnames(x)))
  weights <- integer(M)
  remap   <- match(membership, ids)

  for (i in seq_len(M)) {
    rows <- which(remap == i)
    weights[i]   <- length(rows)
    sub          <- x[rows, , drop = FALSE]
    centers[i, ] <- colMeans(sub, na.rm = TRUE)
    scales[i, ]  <- if (length(rows) > 1L)
      matrixStats::colVars(sub, na.rm = TRUE) else rep(0, p)
  }
  scales[is.na(scales)] <- 0

  out <- list(
    centers    = centers,
    weights    = weights,
    scales     = scales,
    membership = remap,
    x_scale    = rowSums(scales),
    engine     = "native",
    n_obs      = n
  )
  class(out) <- "data_nugget"
  out
}

# delegate path: build a data_nugget object from the datanugget package --------
.dn_from_datanugget <- function(x, ...) {
  dn <- datanugget::create.DN(x, ...)
  dn <- tryCatch(datanugget::refine.DN(x, dn), error = function(e) dn)
  tab <- dn[["Data Nuggets"]]
  cnm <- grep("^Center", colnames(tab), value = TRUE)
  snm <- grep("^Scale|^Shape", colnames(tab), value = TRUE)
  centers <- as.matrix(tab[, cnm, drop = FALSE])
  weights <- as.numeric(tab[["Weight"]])
  scales  <- if (length(snm)) as.matrix(tab[, snm, drop = FALSE]) else
    matrix(0, nrow(centers), ncol(centers))
  out <- list(centers = centers, weights = weights, scales = scales,
              membership = NULL, x_scale = rowSums(scales),
              engine = "datanugget", n_obs = nrow(x))
  class(out) <- "data_nugget"
  out
}

# absorb nuggets smaller than min_size into the nearest surviving centre -------
.absorb_small <- function(xs, membership, min_size) {
  repeat {
    tab <- table(membership)
    small <- as.integer(names(tab)[tab < min_size])
    if (!length(small)) break
    keep <- as.integer(names(tab)[tab >= min_size])
    if (!length(keep)) break
    cent <- t(vapply(keep, function(k)
      colMeans(xs[membership == k, , drop = FALSE]), numeric(ncol(xs))))
    for (s in small) {
      rows <- which(membership == s)
      d <- as.matrix(stats::dist(rbind(colMeans(xs[rows, , drop = FALSE]), cent)))[1, -1]
      membership[rows] <- keep[which.min(d)]
    }
  }
  membership
}

#' @export
print.data_nugget <- function(x, ...) {
  cat("<data_nugget>\n")
  cat("  engine     :", x$engine, "\n")
  cat("  observations:", x$n_obs, "compressed into", length(x$weights), "nuggets\n")
  cat("  reduction  :", sprintf("%.1f%%", 100 * (1 - length(x$weights) / x$n_obs)), "\n")
  cat("  weight range:", paste(range(x$weights), collapse = " - "), "(obs per nugget)\n")
  invisible(x)
}

#' Feature weights derived from data nuggets
#'
#' Computes a per-feature importance score from a data-nugget
#' (\code{\link{create_data_nuggets}}) object built over the \emph{observations}
#' (samples). This score is the
#' data-nugget analogue of the feature variance used by the ABC / M-ABC family:
#' it rewards features that separate the representative groups while discounting
#' within-nugget noise, giving a robust, Big-Data-friendly weighting.
#'
#' @param nuggets A \code{"data_nugget"} object from
#'   \code{\link{create_data_nuggets}} (built with samples in the rows).
#' @param type One of:
#'   \describe{
#'     \item{\code{"between"}}{(default) size-weighted between-nugget variance of
#'       each feature: \eqn{\sum_i w_i (c_{ij} - \bar c_j)^2 / \sum_i w_i}.}
#'     \item{\code{"ratio"}}{a separability / pseudo-\eqn{F} score,
#'       \eqn{\mathrm{between}_j / (\mathrm{between}_j + \overline{\mathrm{within}}_j)}.}
#'     \item{\code{"inv_scale"}}{inverse mean within-nugget variance, which
#'       down-weights noisy features.}
#'   }
#'
#' @return A named numeric vector of length \code{ncol(centers)} of
#'   non-negative feature scores.
#'
#' @references \insertRef{Cherasia2023}{MosaiClusteR}
#' @seealso \code{\link{create_data_nuggets}}, \code{\link{M_ABCpp}}
#' @examples
#' set.seed(1)
#' x <- cbind(signal = c(rnorm(50, 0), rnorm(50, 5)),  # separates groups
#'            noise  = rnorm(100))                      # pure noise
#' dn <- create_data_nuggets(x, max_nuggets = 20)
#' nugget_feature_weights(dn, type = "between")
#' @export
nugget_feature_weights <- function(nuggets, type = c("between", "ratio", "inv_scale")) {
  type <- match.arg(type)
  if (!inherits(nuggets, "data_nugget"))
    stop("nuggets must be a 'data_nugget' object (see create_data_nuggets).")
  C <- nuggets$centers
  w <- nuggets$weights
  W <- sum(w)
  cbar    <- colSums(C * w) / W
  between <- colSums(w * sweep(C, 2, cbar)^2) / W
  within  <- if (!is.null(nuggets$scales))
    colSums(nuggets$scales * w) / W else rep(0, ncol(C))

  score <- switch(type,
    between   = between,
    ratio     = between / (between + within + .Machine$double.eps),
    inv_scale = 1 / (within + .Machine$double.eps)
  )
  score[!is.finite(score)] <- 0
  stats::setNames(as.numeric(score), colnames(C))
}
