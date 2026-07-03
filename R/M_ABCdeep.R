# =============================================================================
# M_ABCdeep.R
# Multi-source M-ABC with a deep-learning (autoencoder) Step 4.
# Function: M_ABCdeep
#
# Mirror of M_ABCpp: per source we run ABCdeep.SingleInMultiple and row-stack the
# label matrices, then fuse them into a consensus dissimilarity with the SAME
# f.clustABC.MultiSource (C++ co-clustering kernel when available). Only the
# per-source base clustering differs (autoencoder latent space instead of a raw
# distance -> hclust). The result carries class "M_ABC" so it works with the
# usual helpers (mosaic_labels, ComparePlot, ...).
# =============================================================================

#' Multi-source M-ABC with a deep-learning (autoencoder) Step 4 (M-ABCdeep)
#'
#' The deep-learning counterpart of \code{\link{M_ABCpp}}: runs
#' \code{\link{ABCdeep.SingleInMultiple}} per source (autoencoder embedding +
#' latent clustering for the base partitions), row-stacks the per-iteration label
#' matrices, and builds one consensus dissimilarity via
#' \code{\link{f.clustABC.MultiSource}}. The ABC algorithm of
#' \insertCite{Amaratunga2008}{MosaiClusteR} is preserved end to end; only Step 4
#' is deep. There is no Python/torch dependency and no PCA fallback (see
#' \code{\link{ABCdeep.SingleInMultiple}}).
#'
#' @param List A list of \eqn{K} data matrices over the same objects,
#'   \strong{objects (samples) in rows} and features in columns.
#' @param weighting,normalize,bag,numvar Per-source vectors (recycled to \eqn{K})
#'   forwarded to \code{\link{ABCdeep.SingleInMultiple}}.
#' @param gr Unused; kept for compatibility.
#' @param numsim Iterations per source.
#' @param NC Base-cluster count per source (\code{NULL} = data-driven).
#' @param dissimilarity Consensus dissimilarity, \code{"mabc"} or \code{"abc"}.
#' @param final_linkage,alpha Final clustering controls.
#' @param mds Logical; MDS of the consensus dissimilarities.
#' @param topk_frac Fraction of top-weighted features the autoencoder focuses on
#'   (default 0.2); forwarded to \code{\link{ABCdeep.SingleInMultiple}}.
#' @param latent_dim,hidden_units,ae_epochs,ae_lr,cluster_in_latent,linkage
#'   Autoencoder / latent-clustering controls (scalar or length-\eqn{K}),
#'   forwarded to \code{\link{ABCdeep.SingleInMultiple}}.
#' @param nugget_type,nugget_args Forwarded to the nugget weighting.
#' @param min_samples Warn (per source) when \eqn{N} is below this; the
#'   autoencoder is underdetermined at very small \eqn{N} (default 40).
#' @param seed Optional RNG seed.
#' @return A list of class \code{"M_ABC"} with \code{DistM} and \code{Clust}.
#' @references \insertAllCited{}
#' @seealso \code{\link{ABCdeep.SingleInMultiple}}, \code{\link{M_ABCpp}},
#'   \code{\link{f.clustABC.MultiSource}}
#' @examples
#' data(mosaic_toy)
#' fit <- M_ABCdeep(mosaic_toy$List, numsim = 40, NC = 3, ae_epochs = 40, seed = 1)
#' table(stats::cutree(stats::as.hclust(fit$Clust), 3))
#' @export
M_ABCdeep <- function(List,
                      weighting         = rep("var", length(List)),
                      normalize         = rep(FALSE, length(List)),
                      gr                = c(),
                      bag               = rep(TRUE,  length(List)),
                      numsim            = 500,
                      numvar            = rep(100,   length(List)),
                      NC                = NULL,
                      dissimilarity     = "mabc",
                      final_linkage     = "ward",
                      alpha             = 0.625,
                      mds               = FALSE,
                      topk_frac         = 0.2,
                      latent_dim        = "auto",
                      hidden_units      = 32L,
                      ae_epochs         = 80L,
                      ae_lr             = 1e-3,
                      cluster_in_latent = "ward",
                      linkage           = rep("ward.D2", length(List)),
                      nugget_type       = "between",
                      nugget_args       = list(),
                      min_samples       = 40L,
                      seed              = NULL) {

  K <- length(List)
  if (K < 1L) stop("List must contain at least one data source.")
  if (!is.null(seed)) set.seed(seed)

  # recycle the per-source arguments to length K
  weighting         <- rep_len(weighting,         K)
  normalize         <- rep_len(normalize,         K)
  bag               <- rep_len(bag,               K)
  numvar            <- rep_len(numvar,            K)
  linkage           <- rep_len(linkage,           K)
  hidden_units      <- rep_len(hidden_units,      K)
  ae_epochs         <- rep_len(ae_epochs,         K)
  ae_lr             <- rep_len(ae_lr,             K)
  cluster_in_latent <- rep_len(cluster_in_latent, K)
  latent_dim_list   <- if (length(latent_dim) == 1L) rep(list(latent_dim), K)
                       else as.list(rep_len(latent_dim, K))
  if (length(nugget_args) && all(vapply(nugget_args, is.list, logical(1))))
    nugget_args <- rep_len(nugget_args, K)
  else
    nugget_args <- rep(list(nugget_args), K)

  Results <- NULL
  for (i in seq_len(K)) {
    message("M_ABCdeep: source ", i, " of ", K,
            " (weighting = ", weighting[i], ", latent = ", latent_dim_list[[i]], ") ...")
    norm_i <- if (isFALSE(normalize[[i]])) NULL else normalize[[i]]
    res <- ABCdeep.SingleInMultiple(
      data              = List[[i]],
      weighting         = weighting[i],
      normalize         = norm_i,
      gr                = gr,
      bag               = bag[i],
      numsim            = numsim,
      numvar            = numvar[i],
      NC                = NC,
      topk_frac         = topk_frac,
      latent_dim        = latent_dim_list[[i]],
      hidden_units      = hidden_units[i],
      ae_epochs         = ae_epochs[i],
      ae_lr             = ae_lr[i],
      cluster_in_latent = cluster_in_latent[i],
      linkage           = linkage[i],
      nugget_type       = nugget_type,
      nugget_args       = nugget_args[[i]],
      min_samples       = min_samples
    )
    res[is.na(res)] <- 0L
    Results <- rbind(Results, res)
  }

  message("M_ABCdeep: building consensus dissimilarity",
          if (exists("count_coclustering_cpp", mode = "function")) " (C++ kernel)" else " (pure R)",
          " ...")
  out <- f.clustABC.MultiSource(
    res           = Results,
    dissimilarity = dissimilarity,
    linkage       = final_linkage,
    alpha         = alpha,
    mds           = mds
  )
  attr(out, "method") <- "M-ABCdeep"
  class(out) <- c("M_ABC", class(out))
  out
}
