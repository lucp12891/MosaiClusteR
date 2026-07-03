# =============================================================================
# cluster.R  - single-source base clustering (Layer 2 baseline primitive)
# =============================================================================

#' Single-source base clustering
#'
#' Cluster one data (or distance) matrix with agglomerative hierarchical
#' clustering. This is the Layer-2 baseline primitive of the MoSaIC framework:
#' run it on each tile on its own to see what each modality recovers before
#' integrating. It also underpins \code{\link{HierarchicalEnsembleClustering}}.
#'
#' @param Data A data matrix (objects in rows) or a precomputed distance matrix.
#' @param type \code{"data"} to compute a distance first, or \code{"dist"} if
#'   \code{Data} is already a dissimilarity matrix.
#' @param distmeasure Distance measure (see \code{\link{Distance}}), used when
#'   \code{type = "data"}.
#' @param normalize,method Normalisation controls passed to \code{\link{Distance}}.
#' @param clust Clustering backend; currently \code{"agnes"}.
#' @param linkage Agglomeration method. Default \code{"ward"}.
#' @param alpha Flexible-linkage parameter for \code{cluster::agnes}.
#' @param StopRange Logical; if \code{FALSE} and the dissimilarities fall outside
#'   \code{[0, 1]} they are range-normalised for comparability.
#' @param ... Ignored; accepted for backward compatibility (e.g. \code{gap},
#'   \code{maxK}).
#' @return A list with \code{DistM} (the dissimilarity matrix) and \code{Clust}
#'   (an \code{agnes} object).
#' @seealso \code{\link{Distance}}, \code{\link{M_ABCpp}}
#' @examples
#' data(mosaic_toy)
#' base1 <- Cluster(t(mosaic_toy$List[[1]]), type = "data",
#'                  distmeasure = "euclidean")
#' mosaic_labels(base1, k = 3)[1:6]
#' @export
Cluster <- function(Data, type = c("data", "dist"),
                    distmeasure = "euclidean", normalize = FALSE, method = NULL,
                    clust = "agnes", linkage = "ward", alpha = 0.625,
                    StopRange = TRUE, ...) {
  type <- match.arg(type)
  norm <- if (isFALSE(normalize) || is.null(normalize)) NULL else normalize

  DistM <- if (type == "data")
    as.matrix(Distance(Data, distmeasure = distmeasure, normalize = norm))
  else
    as.matrix(Data)

  if (!StopRange && !(min(DistM) >= 0 && max(DistM) <= 1))
    DistM <- Normalization(DistM, method = "Range")

  Clust <- cluster::agnes(DistM, diss = TRUE, method = linkage, par.method = alpha)
  out <- list(DistM = DistM, Clust = Clust)
  attr(out, "method") <- "Cluster"
  out
}
