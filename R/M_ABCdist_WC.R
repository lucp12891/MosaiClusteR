# =============================================================================
# M_ABCdist_WC.R
# Multi-source M-ABC + WeightedClust fusion. Per source ABCdist.SingleInMultiple
# yields a [0,1] dissimilarity; WeightedClust() then explores convex
# combinations of the source dissimilarities over a weight grid.
# Function: M_ABCdist.WC
# =============================================================================

#' Multi-source M-ABC fused with WeightedClust
#'
#' Runs \code{\link{ABCdist.SingleInMultiple}} on each source, then feeds the
#' per-source dissimilarities to \code{\link{WeightedClust}} (\code{type =
#' "dist"}) to explore convex combinations over a weight grid and return the
#' focus-weight clustering.
#'
#' @inheritParams M_ABCdist
#' @param StopRange Passed to \code{\link{WeightedClust}}.
#' @param wc_weight Weight grid for \code{\link{WeightedClust}}.
#' @param wc_weightclust Focus weight whose result is returned in \code{Clust}.
#' @param wc_linkage,wc_alpha Final clustering controls for \code{\link{WeightedClust}}.
#' @return A list of class \code{"M_ABCdist_WC"} with \code{source_D}, the full
#'   \code{WC} object, the focus \code{DistM} and \code{Clust}, and the
#'   \code{weight_grid}.
#' @references \insertCite{Amaratunga2008}{MosaiClusteR}
#' @seealso \code{\link{M_ABCdist}}, \code{\link{WeightedClust}}
#' @examples
#' data(mosaic_toy)
#' fit <- M_ABCdist.WC(mosaic_toy$List, numsim = 40,
#'                     wc_weight = seq(1, 0, -0.25), wc_weightclust = 0.5)
#' dim(fit$DistM)
#' @export
M_ABCdist.WC <- function(List,
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
                         StopRange      = FALSE,
                         wc_weight      = seq(1, 0, -0.1),
                         wc_weightclust = 0.5,
                         wc_linkage     = "ward",
                         wc_alpha       = 0.625,
                         mds            = FALSE,
                         nugget_type    = "between",
                         nugget_args    = list()) {
  if (!exists("WeightedClust", mode = "function", inherits = TRUE))
    stop("M_ABCdist.WC() needs WeightedClust(); it ships with MosaiClusteR.")
  K <- length(List)
  if (K < 1L) stop("List must contain at least one data source.")
  distmeasure <- rep_len(distmeasure, K); weighting <- rep_len(weighting, K)
  bag <- as.logical(rep_len(bag, K)); numvar <- as.integer(rep_len(numvar, K))
  linkage <- rep_len(linkage, K)
  if (!is.list(normalize)) normalize <- as.list(normalize)
  normalize <- rep_len(normalize, K)
  if (length(nugget_args) && all(vapply(nugget_args, is.list, logical(1))))
    nugget_args <- rep_len(nugget_args, K) else nugget_args <- rep(list(nugget_args), K)

  List <- lapply(List, as.matrix)
  samp_names <- rownames(List[[1]]); N <- nrow(List[[1]])
  for (k in seq_len(K))
    if (nrow(List[[k]]) != N)
      stop("Source ", k, " has ", nrow(List[[k]]), " rows; expected ", N, ".")

  source_D <- vector("list", K)
  for (k in seq_len(K)) {
    message("M_ABCdist.WC: source ", k, " of ", K, " (weighting=", weighting[k], ") ...")
    source_D[[k]] <- ABCdist.SingleInMultiple(
      data = List[[k]], distmeasure = distmeasure[k], weighting = weighting[k],
      normalize = normalize[[k]], gr = gr, bag = bag[k], numsim = as.integer(numsim),
      numvar = numvar[k], linkage = linkage[k], NC = NC,
      nugget_type = nugget_type, nugget_args = nugget_args[[k]])$D
  }

  WC <- WeightedClust(
    List = source_D, type = "dist", distmeasure = rep(NULL, K),
    normalize = rep(FALSE, K), method = rep(NULL, K), StopRange = StopRange,
    weight = wc_weight, weightclust = wc_weightclust, clust = "agnes",
    linkage = wc_linkage, alpha = wc_alpha)

  DistM <- matrix(as.double(WC$Clust$DistM), N, N,
                  dimnames = list(samp_names, samp_names))
  diag(DistM) <- 0
  Clust <- WC$Clust$Clust

  if (isTRUE(mds)) {
    co <- stats::cmdscale(stats::as.dist(DistM), k = 2L)
    plot(co, pch = 19, col = "darkgreen", xlab = "Dim 1", ylab = "Dim 2",
         main = "MDS - M-ABCdist-WC")
  }

  out <- list(source_D = source_D, WC = WC, DistM = DistM, Clust = Clust,
              weight_grid = names(WC$Results))
  if (!is.null(NC2)) {
    nc <- as.integer(NC2)
    if (nc >= 2L && nc < N) out$cut <- stats::cutree(stats::as.hclust(Clust), k = nc)
  }
  attr(out, "method") <- "M-ABCdist-WC"
  class(out) <- c("M_ABCdist_WC", "list")
  out
}
