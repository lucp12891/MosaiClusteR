# =============================================================================
# data-sim.R  - synthetic multi-source data generator + bundled example data
# =============================================================================

#' Simulate a multi-source data set with planted clusters
#'
#' Generates \code{K} feature-by-sample matrices that share a common set of
#' \code{n} samples drawn from \code{g} ground-truth groups. Each source has its
#' own informative features (shifted between groups) plus noise features, so the
#' clusters are only partially visible in any single source - the canonical
#' setting in which multi-source integration helps.
#'
#' @param n Number of samples (objects).
#' @param g Number of ground-truth groups.
#' @param p Per-source feature counts (recycled to length \code{K}).
#' @param informative Per-source number of group-informative features.
#' @param effect Between-group mean shift for informative features.
#' @param K Number of sources.
#' @param seed Optional RNG seed.
#' @return A list with \code{List} (the \code{K} matrices, \strong{objects in
#'   rows}, features in columns) and \code{truth} (the integer group label per
#'   object).
#' @examples
#' sim <- mosaic_sim(n = 30, g = 3, K = 2, seed = 1)
#' lengths(lapply(sim$List, dim))
#' @export
mosaic_sim <- function(n = 60, g = 3, p = c(120, 150), informative = c(20, 25),
                       effect = 2.5, K = 2, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  p           <- rep_len(p, K)
  informative <- rep_len(informative, K)
  truth   <- sort(rep_len(seq_len(g), n))
  s_names <- sprintf("obj%02d", seq_len(n))

  List <- vector("list", K)
  for (k in seq_len(K)) {
    X <- matrix(stats::rnorm(n * p[k]), nrow = n, ncol = p[k])  # objects x features
    inf <- seq_len(informative[k])
    # each informative feature separates a (cyclically chosen) group
    for (j in inf) {
      gj <- ((j - 1) %% g) + 1
      X[truth == gj, j] <- X[truth == gj, j] + effect
    }
    rownames(X) <- s_names
    colnames(X) <- sprintf("src%d_f%03d", k, seq_len(p[k]))
    List[[k]] <- X
  }
  names(List) <- sprintf("source%d", seq_len(K))
  list(List = List, truth = stats::setNames(truth, s_names))
}

#' Toy multi-omics example data
#'
#' A small two-source example produced by \code{\link{mosaic_sim}} with three
#' planted groups, bundled for examples and tests.
#'
#' @format A list with two elements:
#' \describe{
#'   \item{List}{A list of two numeric matrices with the 60 shared objects in
#'     rows and features in columns (120 and 150 features).}
#'   \item{truth}{Integer vector of the 60 ground-truth group labels.}
#' }
#' @examples
#' data(mosaic_toy)
#' str(mosaic_toy, max.level = 2)
"mosaic_toy"
