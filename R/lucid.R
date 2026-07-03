# =============================================================================
# lucid.R  - wrapper for LUCID (Latent Unknown Clustering Integrating multi-omics
# Data with phenotypic traits). Delegates to the LUCIDus package (Suggests).
# Zhao, Jia, Goodrich & Conti (2024), The R Journal 16(2).
# =============================================================================

#' LUCID: latent unknown clustering integrating omics with an outcome
#'
#' Thin wrapper around \pkg{LUCIDus}' \code{estimate_lucid}. LUCID is a
#' \emph{model-based, supervised} integrative clustering method: it defines a
#' categorical latent variable (the clusters) jointly from upstream exposures
#' \code{G}, multi-omics data \code{Z} and an outcome \code{Y} via an EM
#' algorithm (a quasi-mediation model). Unlike the unsupervised methods in
#' MosaiClusteR it requires an outcome and returns posterior cluster membership.
#'
#' @param G Exposure/predictor matrix (objects in rows).
#' @param Z Omics matrix, or a list of matrices for the parallel/serial models
#'   (objects in rows).
#' @param Y Outcome vector (length = number of objects).
#' @param K Number of latent clusters.
#' @param lucid_model One of \code{"early"}, \code{"parallel"}, \code{"serial"}.
#' @param family Outcome family, \code{"normal"} or \code{"binary"}.
#' @param ... Further arguments forwarded to \code{LUCIDus::estimate_lucid}.
#' @return A list of class \code{"LUCID"} with the fitted \code{model} and the
#'   posterior \code{cluster} assignment per object.
#' @references Zhao, Y., Jia, K., Goodrich, J. & Conti, D. V. (2024). LUCIDus:
#'   an R package for implementing LUCID with phenotypic traits. The R Journal
#'   16(2).
#' @seealso \code{\link{M_ABCpp}}, \code{\link{intNMF}}
#' @examples
#' \dontrun{
#' # requires the 'LUCIDus' package
#' data(mosaic_toy)
#' Z <- mosaic_toy$List[[1]]
#' G <- mosaic_toy$List[[2]][, 1:3]
#' Y <- as.numeric(mosaic_toy$truth)
#' fit <- LUCID(G = G, Z = Z, Y = Y, K = 3)
#' table(fit$cluster)
#' }
#' @export
LUCID <- function(G, Z, Y, K, lucid_model = c("early", "parallel", "serial"),
                  family = c("normal", "binary"), ...) {
  if (!requireNamespace("LUCIDus", quietly = TRUE))
    stop("LUCID() requires the suggested package 'LUCIDus'. ",
         "Install it with install.packages('LUCIDus').")
  lucid_model <- match.arg(lucid_model)
  family      <- match.arg(family)

  fit <- LUCIDus::estimate_lucid(
    lucid_model = lucid_model, G = as.matrix(G), Z = Z, Y = as.matrix(Y),
    K = K, family = family, ...)

  # posterior cluster membership across LUCIDus versions
  post <- fit$inclusion.p
  if (is.null(post)) post <- fit$post.p
  if (is.null(post)) post <- fit$z
  cluster <- if (!is.null(post)) max.col(as.matrix(post), ties.method = "first")
             else fit$pred

  out <- list(model = fit,
              cluster = stats::setNames(as.integer(cluster), rownames(as.matrix(G))),
              K = K)
  attr(out, "method") <- "LUCID"
  class(out) <- c("LUCID", "list")
  out
}
