# =============================================================================
# M_ABCdist.R
# Multi-source M-ABC, distance-accumulation variant. Each source yields a
# [0,1]-scaled mean-distance matrix from ABCdist.SingleInMultiple; these are
# combined as a weighted mean D_multi = sum_k w_k D_k and clustered.
# Function: M_ABCdist
# =============================================================================

#' Multi-source M-ABC with distance accumulation
#'
#' Runs \code{\link{ABCdist.SingleInMultiple}} on each source to obtain a
#' \code{[0, 1]} dissimilarity, combines them as a weighted mean
#' \eqn{D = \sum_k w_k D_k} and applies a final agglomerative clustering.
#'
#' @inheritParams M_ABCpp
#' @param source_weights Optional non-negative source weights (normalised to sum
#'   to one); \code{NULL} = equal weights.
#' @param NC2 Optional number of clusters to cut the final tree into.
#' @return A list of class \code{"M_ABCdist"} with \code{source_D}, \code{DistM},
#'   \code{weights} and \code{Clust}.
#' @references \insertCite{Amaratunga2008}{MosaiClusteR}
#' @seealso \code{\link{M_ABCpp}}, \code{\link{ABCdist.SingleInMultiple}},
#'   \code{\link{M_ABCdist.WC}}
#' @examples
#' data(mosaic_toy)
#' fit <- M_ABCdist(mosaic_toy$List, numsim = 40, weighting = c("var", "nugget"))
#' dim(fit$DistM)
#' @export
M_ABCdist <- function(List,
                      source_weights = NULL,
                      distmeasure    = rep("euclidean", length(List)),
                      weighting      = rep("var",       length(List)),
                      normalize      = rep(list(NULL),  length(List)),
                      gr             = c(),
                      bag            = rep(TRUE,        length(List)),
                      numsim         = 1000L,
                      numvar         = rep(100L,        length(List)),
                      linkage        = rep("ward.D",    length(List)),
                      NC             = NULL,
                      NC2            = NULL,
                      final_linkage  = "ward",
                      alpha          = 0.625,
                      mds            = FALSE,
                      nfeat          = NULL,
                      nugget_type    = "between",
                      nugget_args    = list()) {
  if (!is.list(List)) stop("List must contain at least one data source.")
  K <- length(List)
  .rec <- function(x, nm) rep_len(x, K)
  distmeasure <- .rec(distmeasure); weighting <- .rec(weighting)
  bag <- as.logical(.rec(bag)); numvar <- as.integer(.rec(numvar))
  linkage <- .rec(linkage)
  if (!is.list(normalize)) normalize <- as.list(normalize)
  normalize <- rep_len(normalize, K)
  if (length(nugget_args) && all(vapply(nugget_args, is.list, logical(1))))
    nugget_args <- rep_len(nugget_args, K) else nugget_args <- rep(list(nugget_args), K)

  if (is.null(source_weights)) source_weights <- rep(1.0 / K, K)
  else {
    if (length(source_weights) != K) stop("source_weights must have length ", K, ".")
    source_weights <- source_weights / sum(source_weights)
  }

  List <- lapply(List, as.matrix)
  samp_names <- rownames(List[[1]]); N <- nrow(List[[1]])
  for (k in seq_len(K))
    if (nrow(List[[k]]) != N)
      stop("Source ", k, " has ", nrow(List[[k]]), " rows; expected ", N, ".")

  source_D <- vector("list", K)
  for (k in seq_len(K)) {
    message("M_ABCdist: source ", k, " of ", K,
            " (weight=", round(source_weights[k], 4), ", weighting=", weighting[k], ") ...")
    source_D[[k]] <- ABCdist.SingleInMultiple(
      data = List[[k]], distmeasure = distmeasure[k], weighting = weighting[k],
      normalize = normalize[[k]], gr = gr, bag = bag[k], numsim = as.integer(numsim),
      numvar = numvar[k], linkage = linkage[k], NC = NC, nfeat = nfeat,
      nugget_type = nugget_type, nugget_args = nugget_args[[k]])$D
  }

  D_multi <- Reduce("+", Map("*", source_weights, source_D))
  D_multi <- matrix(as.double(D_multi), N, N); diag(D_multi) <- 0
  rownames(D_multi) <- colnames(D_multi) <- samp_names

  if (isTRUE(mds)) {
    co <- stats::cmdscale(stats::as.dist(D_multi), k = 2L)
    plot(co, pch = 19, col = "darkgreen", xlab = "Dim 1", ylab = "Dim 2",
         main = "MDS - M-ABCdist combined dissimilarity")
  }

  Clust <- cluster::agnes(stats::as.dist(D_multi), diss = TRUE,
                          method = final_linkage, par.method = alpha)
  out <- list(source_D = source_D, DistM = D_multi,
              weights = source_weights, Clust = Clust)
  if (!is.null(NC2)) {
    nc <- as.integer(NC2)
    if (nc >= 2L && nc < N) out$cut <- stats::cutree(stats::as.hclust(Clust), k = nc)
  }
  attr(out, "method") <- "M-ABCdist"
  class(out) <- c("M_ABCdist", "list")
  out
}
