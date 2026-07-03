# =============================================================================
# nugget-clustering.R
# Weighted clustering of (data-nugget) representatives - the clustering side of
# the data-nugget toolkit (WCluster: Dey, Duan, Cabrera & Cheng 2025;
# Cherasia et al. 2023). Weighted k-means minimises the Weighted Within-Cluster
# Sum of Squares (WWCSS); weighted hierarchical clustering uses the weighted
# Ward merge cost SS(A union B) - SS(A) - SS(B).
# =============================================================================

#' Weighted k-means (minimises WWCSS)
#'
#' k-means in which each observation carries a non-negative weight. Cluster
#' centres are weighted means and the objective is the Weighted Within-Cluster
#' Sum of Squares \eqn{\sum_c \sum_{i \in c} w_i \lVert x_i - \mu_c \rVert^2}.
#'
#' @param x Numeric matrix, observations in rows.
#' @param k Number of clusters.
#' @param weights Non-negative weight per observation (default all 1).
#' @param nstart Number of random restarts (best WWCSS kept).
#' @param max_iter Maximum Lloyd iterations per start.
#' @return A list with \code{cluster}, \code{centers} and \code{wwcss}.
#' @references Dey, T., Duan, Y., Cabrera, J. & Cheng, X. (2025). WCluster:
#'   Clustering and PCA with Weights, and Data Nuggets Clustering.
#' @seealso \code{\link{nugget_cluster}}, \code{\link{Whclust}}
#' @examples
#' set.seed(1)
#' x <- rbind(matrix(rnorm(20 * 3, 0), 20), matrix(rnorm(20 * 3, 5), 20))
#' Wkmeans(x, k = 2, weights = runif(40, 1, 5))$cluster[1:5]
#' @export
Wkmeans <- function(x, k, weights = rep(1, nrow(x)), nstart = 10, max_iter = 100) {
  x <- as.matrix(x); n <- nrow(x)
  w <- as.numeric(weights); w[w < 0] <- 0
  k <- max(1L, min(as.integer(k), n))

  wmeans <- function(rows) {
    if (!length(rows)) return(NULL)
    colSums(x[rows, , drop = FALSE] * w[rows]) / sum(w[rows])
  }
  best <- NULL; best_obj <- Inf
  for (s in seq_len(nstart)) {
    centers <- x[sample.int(n, k, prob = w / sum(w)), , drop = FALSE]
    cl <- integer(n)
    for (it in seq_len(max_iter)) {
      D  <- as.matrix(stats::dist(rbind(centers, x)))[seq_len(k), -seq_len(k), drop = FALSE]
      new <- max.col(-t(D), ties.method = "first")
      if (identical(new, cl)) break
      cl <- new
      for (c in seq_len(k)) {
        m <- wmeans(which(cl == c))
        if (!is.null(m)) centers[c, ] <- m
      }
    }
    obj <- sum(vapply(seq_len(k), function(c) {
      rows <- which(cl == c); if (!length(rows)) return(0)
      sum(w[rows] * rowSums(sweep(x[rows, , drop = FALSE], 2, centers[c, ])^2))
    }, numeric(1)))
    if (obj < best_obj) { best_obj <- obj; best <- list(cluster = cl, centers = centers) }
  }
  list(cluster = stats::setNames(best$cluster, rownames(x)),
       centers = best$centers, wwcss = best_obj)
}

# weighted Ward agglomeration -> hclust object (centres in rows, weights w)
.weighted_ward <- function(C, w) {
  C <- as.matrix(C); m <- nrow(C)
  active <- seq_len(m); wt <- w; cen <- C
  id <- -seq_len(m)                       # hclust leaf ids are negative
  merge <- matrix(0L, m - 1L, 2L); height <- numeric(m - 1L)
  ward <- function(a, b)
    (wt[a] * wt[b] / (wt[a] + wt[b])) * sum((cen[a, ] - cen[b, ])^2)
  for (step in seq_len(m - 1L)) {
    best <- Inf; bi <- bj <- NA
    for (ii in seq_along(active)) for (jj in seq_len(ii - 1L)) {
      a <- active[ii]; b <- active[jj]; d <- ward(a, b)
      if (d < best) { best <- d; bi <- a; bj <- b }
    }
    height[step] <- best
    merge[step, ] <- sort(c(id[bi], id[bj]))
    nw  <- wt[bi] + wt[bj]
    cen[bi, ] <- (wt[bi] * cen[bi, ] + wt[bj] * cen[bj, ]) / nw
    wt[bi] <- nw; id[bi] <- step
    active <- setdiff(active, bj)
  }
  ord <- stats::hclust(stats::dist(C), method = "ward.D2")$order  # plausible leaf order
  structure(list(merge = merge, height = cumsum(sort(height)) / seq_along(height),
                 order = ord, labels = rownames(C), method = "weighted.ward",
                 dist.method = "euclidean"), class = "hclust")
}

#' Weighted hierarchical clustering (weighted Ward)
#'
#' Agglomerative clustering of weighted observations using the weighted Ward
#' merge cost. Intended for the small set of data-nugget centres.
#'
#' @param x Numeric matrix, observations (e.g. nugget centres) in rows.
#' @param weights Non-negative weight per observation.
#' @return An \code{\link[stats]{hclust}} object.
#' @references Dey et al. (2025). WCluster.
#' @seealso \code{\link{Wkmeans}}, \code{\link{nugget_cluster}}
#' @examples
#' set.seed(1)
#' x <- rbind(matrix(rnorm(8 * 3, 0), 8), matrix(rnorm(8 * 3, 5), 8))
#' stats::cutree(Whclust(x, weights = rep(2, 16)), k = 2)[1:4]
#' @export
Whclust <- function(x, weights = rep(1, nrow(x))) {
  x <- as.matrix(x)
  if (is.null(rownames(x))) rownames(x) <- as.character(seq_len(nrow(x)))
  .weighted_ward(x, as.numeric(weights))
}

#' Data-nugget clustering
#'
#' Cluster a (potentially large) data set by first compressing its observations
#' into weighted data nuggets (\code{\link{create_data_nuggets}}), clustering the
#' nugget centres with a weighted clusterer, and mapping the nugget labels back
#' to the original observations. This is the clustering analogue of the
#' data-nugget weighting used by \code{\link{M_ABCpp}} and follows
#' \code{DN.Wkmeans} / \code{DN.Whclust} from the WCluster package.
#'
#' @param data A numeric matrix (objects in rows) \emph{or} a precomputed
#'   \code{\link{create_data_nuggets}{data_nugget}} object. For multi-source
#'   data, column-bind the per-source matrices first (the direct/\code{ADC}
#'   strategy) or pass a single fused matrix.
#' @param k Number of clusters.
#' @param method Weighted clusterer for the nugget centres: \code{"kmeans"}
#'   (WWCSS) or \code{"ward"} (weighted Ward).
#' @param max_nuggets,nugget_args Passed to \code{\link{create_data_nuggets}}
#'   when \code{data} is a raw matrix.
#' @return A list of class \code{"nugget_cluster"} with \code{cluster} (labels
#'   for the original objects), \code{nugget_cluster} (labels for the nuggets)
#'   and the \code{nuggets} object.
#' @references Cherasia et al. (2023); Dey et al. (2025), WCluster.
#' @seealso \code{\link{create_data_nuggets}}, \code{\link{Wkmeans}}
#' @examples
#' data(mosaic_toy)
#' X <- do.call(cbind, mosaic_toy$List)        # fuse sources (objects in rows)
#' fit <- nugget_cluster(X, k = 3, max_nuggets = 30)
#' table(fit$cluster)
#' @export
nugget_cluster <- function(data, k, method = c("kmeans", "ward"),
                           max_nuggets = NULL, nugget_args = list()) {
  method <- match.arg(method)
  dn <- if (inherits(data, "data_nugget")) data
        else do.call(create_data_nuggets,
                     c(list(x = as.matrix(data), max_nuggets = max_nuggets), nugget_args))
  if (is.null(dn$membership))
    stop("nugget_cluster needs nugget membership; build nuggets with the ",
         "native engine (create_data_nuggets(..., engine='native')).")

  nlab <- if (method == "kmeans")
    Wkmeans(dn$centers, k = k, weights = dn$weights)$cluster
  else
    stats::cutree(Whclust(dn$centers, weights = dn$weights), k = k)

  obj_lab <- nlab[dn$membership]
  if (inherits(data, "data_nugget")) names(obj_lab) <- NULL
  else names(obj_lab) <- rownames(as.matrix(data))

  out <- list(cluster = obj_lab, nugget_cluster = nlab, nuggets = dn,
              k = k, method = method)
  attr(out, "method") <- paste0("DN.", method)
  class(out) <- c("nugget_cluster", "list")
  out
}
